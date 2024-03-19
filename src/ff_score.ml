(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

(* read in a pqrs file; compute the UFF non-bonded energy for each conformer;
   !!! _WITHOUT_ prior minimization !!! *)

let main () =
  let input_fn = Sys.argv.(1) in
  let ligands = Mol.ligands_of_pqrs_file input_fn in
  L.iter (fun lig ->
      let name = Mol.get_name lig in
      let ene = Mol.ene_intra_UFFNB_brute lig in
      Printf.printf "%s\t%g\n" name ene
    ) ligands

let () = main ()
