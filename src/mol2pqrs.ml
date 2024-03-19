(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

(* convert a .mol2 file _WITH_ partial charges into a .pqrs file;
 * all molecules are converted *)

open Printf

let process_ligand verbose mol =
  let buff = Buffer.create 10_000 in
  let name = Mol2.get_name mol in
  begin
    try
      let vertices, graph = Mol_graph.create mol in
      let dists = Mol_graph.dist_matrix vertices graph in
      let degrees = Mol_graph.degrees graph in
      (if verbose then
         let () = printf "degrees:\n" in
         Mol_graph.print_degrees degrees
      );
      let num_atoms = Mol2.num_atoms mol in
      let rot_bonds = Mol_graph.list_rotatable_bonds mol vertices degrees graph in
      let num_rbonds = L.length rot_bonds in
      bprintf buff "%d:%d:%s\n" num_atoms num_rbonds name;
      let atoms = Mol2.get_atoms mol in
      let all_null_charges = A.for_all Mol2.atom_charge_is_null atoms in
      if all_null_charges then
        Log.error "all null charges: %s" name;
      (* atoms *)
      A.iter (fun atom ->
          bprintf buff "%s\n" (Mol2.pqrs_line_of_atom atom)
        ) atoms;
      (* RotBonds *)
      L.iter (fun bond ->
          let rot_group = Mol_graph.compute_rot_group vertices graph bond in
          Mol_graph.bprint_rot_group buff graph bond rot_group
        ) rot_bonds;
      (* dist matrix *)
      Mol_graph.bprint_dist_matrix buff graph dists
    with Mol_graph.Disconnected_atom -> Log.error "disconnected atom in %s" name
  end;
  Buffer.contents buff

let process_receptor mol =
  let buff = Buffer.create 10_000 in
  let name = Mol2.get_name mol in
  let num_atoms = Mol2.num_atoms mol in
  bprintf buff "%d:%s\n" num_atoms name;
  let atoms = Mol2.get_atoms mol in
  let all_null_charges = A.for_all Mol2.atom_charge_is_null atoms in
  (if all_null_charges then
     let () = Log.fatal "all null charges: %s" name in
     exit 1
  );
  (* atoms *)
  A.iter (fun atom ->
      bprintf buff "%s\n" (Mol2.pqrs_line_of_atom atom)
    ) atoms;
  Buffer.contents buff

(* .xyz output format *)
let process_molecule_xyz counter mol =
  let buff = Buffer.create 10_000 in
  let name = Mol2.get_name mol in
  let num_atoms = Mol2.num_atoms mol in
  bprintf buff "%d\n" num_atoms;
  bprintf buff "%d\t%s\t\n" !counter name;
  incr counter;
  let atoms = Mol2.get_atoms mol in
  (* atoms *)
  A.iter (fun atom ->
      bprintf buff "%s\n" (Mol2.xyz_line_of_atom atom)
    ) atoms;
  Buffer.contents buff

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  Log.(set_prefix_builder short_prefix_builder);
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s\n  \
              [-lig <ligands.mol2>]: input file (incompatible w/ -rec)\n  \
              [-rec <protein.mol2>]: input file (incompatible w/ -lig)\n  \
              [-xyz]: output in .xyz file format (with -lig)\n  \
              [-np <int>]: nprocs (if many ligands and -lig)\n  \
              [-o <filename.pqrs>]: output file\n  \
              [-v]: debug info on stdout\n"
       Sys.argv.(0);
     exit 1);
  let maybe_ligs_fn = CLI.get_string_opt ["-lig"] args in
  let maybe_rec_fn = CLI.get_string_opt ["-rec"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let nprocs = ref (CLI.get_int_def ["-np"] args 1) in
  let xyz_out = CLI.get_set_bool ["-xyz"] args in
  let verbose = CLI.get_set_bool ["-v"] args in
  let counter = ref 0 in
  CLI.finalize (); (* ------------------------------------------------------ *)
  let reader_fun, input_fn =
    match (maybe_ligs_fn, maybe_rec_fn) with
    | Some ligs_fn, None ->
      ((if xyz_out then
          (fun _verb mol -> process_molecule_xyz counter mol)
        else
          process_ligand),
       ligs_fn)
    | None, Some rec_fn ->
      ((fun _verb mol ->
          nprocs := 1; (* protein case: not parallel *)
          process_receptor mol),
       rec_fn)
    | _, _ -> failwith "use either -lig OR -rec" in
  LO.with_infile_outfile input_fn output_fn (fun input output ->
      Parany.run !nprocs ~preserve:true (* keep input order *)
        ~demux:(fun () ->
            try Mol2.read_one input
            with End_of_file -> raise Parany.End_of_input)
        ~work:(function None -> ""
                      | Some mol -> reader_fun verbose mol)
        ~mux:(output_string output)
    )

let () = main ()
