(* Copyright (C) 2024, Francois Berenger
   Tsuda laboratory, The University of Tokyo,
   5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.

   Extract ligand-defined binding site
   from lig.mol2 out of rec.mol2 *)

open Printf

module A = BatArray
module CLI = Minicli.CLI
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module V3 = Vector3

module Point_3D = struct
  type t = V3.t
  let dist = V3.dist
end

module BST = Bst.Bisec_tree.Make(Point_3D)

let main () =
  Log.(set_prefix_builder short_prefix_builder);
  Log.color_on ();
  Log.(set_log_level INFO);
  let default_cutoff = 5.0 in (* (Angstrom) around ligand's heavy atoms *)
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s\n  \
              -l <ligand.mol2>: xtal ligand input file\n  \
              -p <protein.mol2>: receptor protein input file\n  \
              [-d <float>]: distance cutoff (default=%.2f)\n  \
              -o <output.pqrs>: ligand-defined binding site output file\n  \
              [-v]: verbose/debug mode\n"
       Sys.argv.(0) default_cutoff;
     exit 1);
  let verbose = CLI.get_set_bool ["-v"] args in
  if verbose then Log.(set_log_level DEBUG);
  let ligand_fn = CLI.get_string ["-l"] args in
  let protein_fn = CLI.get_string ["-p"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let cutoff = CLI.get_float_def ["-d"] args default_cutoff in
  CLI.finalize (); (* ----------------------------------------------------- *)
  LO.with_out_file output_fn (fun output ->
      (* read ligand from MOL2 *)
      let lig = Mol2.read_one_from_file ligand_fn in
      let lig_coords = Mol2.get_coords lig in
      (* index ligand using a BST *)
      let lig_bst = BST.(create 1 Two_bands lig_coords) in
      (* read protein from MOL2 *)
      let prot = Mol2.read_one_from_file protein_fn in
      let prot_atoms = Mol2.get_atoms prot in
      (* output carved out protein atoms to pqrs file *)
      A.iter (fun prot_atom ->
          let xyz = Mol2.get_atom_coord prot_atom in
          let _nearest, dist = BST.nearest_neighbor xyz lig_bst in
          if dist <= cutoff then
            let pqrs = Mol2.pqrs_line_of_atom prot_atom in
            fprintf output "%s\n" pqrs
        ) prot_atoms
    )

let () = main ()
