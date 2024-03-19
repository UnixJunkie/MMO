(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

(* Only keep molecules OK for QM and not too flexible *)

open Printf

module IntSet = BatSet.Int

(* C,H,N,O,S,F,Cl *)
let ani2_supported_anums =
  IntSet.of_list [1;  (* H *) 
                  6;  (* C *) 
                  7;  (* N *) 
                  8;  (* O *) 
                  9;  (* F *) 
                  16; (* S *) 
                  17; (* Cl *)]

let atoms_filter anums_set mol =
  A.for_all (fun anum ->
      IntSet.mem anum anums_set
    ) (Mol.get_anums mol)

(* let rot_bonds_filter mol = *)
(*   Mol.num_rbonds mol <= 7 *)

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  Log.(set_prefix_builder short_prefix_builder);
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s\n  \
              [-i <ligands.mol2>]: input file\n  \
              [-o <filtered.mol2>]: output file\n  \
              [-np <int>]: nprocs (if many ligands)\n  \
              [-v]: verbose mode\n"
       Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let nprocs = CLI.get_int_def ["-np"] args 1 in
  let _verbose = CLI.get_set_bool ["-v"] args in
  let rejected = ref 0 in
  CLI.finalize (); (* ------------------------------------------------------ *)
  let mol_ligands = Mol.ligands_of_mol2_file nprocs input_fn in
  let mol2_blocks = Mol2.read_all_mol2_blocks input_fn in
  let n = L.length mol_ligands in
  let m = L.length mol2_blocks in
  (if n <> m then
     let () = Log.fatal "mols: %d <> mol2_blocks: %d" n m in
     exit 1
  );
  LO.with_out_file output_fn (fun out ->
      L.iter2 (fun mol mol2_block ->
          if atoms_filter ani2_supported_anums mol then
            (* && rot_bonds_filter mol then *)
            A.iter (fprintf out "%s\n") mol2_block
          else
            incr rejected
        ) mol_ligands mol2_blocks
    );
  Log.info "rejected %d/%d from %s" !rejected n input_fn

let () = main ()
