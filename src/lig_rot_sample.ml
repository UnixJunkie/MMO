(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
 *
 * rotational SO(3) sampling of a ligand *)

module Quat = Mmo.Quat
module Rot = Mmo.Rot

let main () =
  let argc = A.length Sys.argv in
  (if argc <> 4 then
     let () =
       Printf.eprintf "usage:\n%s num_samples input.mol2 output.mol2\n"
         Sys.argv.(0) in
     exit 1
  );
  let n = int_of_string Sys.argv.(1) in
  let mol2_in_fn = Sys.argv.(2) in (* SINGLE MOLECULE input file *)
  let mol2_out_fn = Sys.argv.(3) in (* output file *)
  let rotations =
    A.map (fun quat ->
        let axis, angle = Quat.to_axis_angle quat in
        Rot.of_axis_angle axis angle
      ) (SO3.sample n) in
  (* (\* read mol2 ligand *\) *)
  (* let _mol2 = LO.with_in_file mol2_in_fn Mol2.read_one in *)
  (* save current center; center ligand *)
  let tmp_pqrs_out_fn = Fn.temp_file "lrs_" ".pqrs" in
  Utls.mol2pqrs 1 mol2_in_fn tmp_pqrs_out_fn;
  let mol2 = Mol2.read_one_from_file mol2_in_fn in
  (* create Mol.read_one_from_pqrs_file ? *)
  let mol = match Mol.ligands_of_pqrs_file tmp_pqrs_out_fn with
    | [x] -> x
    | [] -> failwith ("lig_rot: no ligand in " ^ tmp_pqrs_out_fn)
    | _ -> failwith ("lig_rot: several ligands in " ^ tmp_pqrs_out_fn) in
  Sys.remove tmp_pqrs_out_fn; (* cleanup *)
  let orig_center = Mol.get_center mol in
  (* center, apply rotation then translate back to initial position *)
  let rotated_copies =
    A.map (fun rot ->
        Mol.center_rotate_translate_copy mol rot orig_center
      ) rotations in
  (* write out as mol2 *)
  LO.with_out_file mol2_out_fn (fun out ->
      A.iter (fun mol' ->
          let mol2' = Mol.update_mol2 mol2 mol' in
          Mol2.output_one out mol2'
        ) rotated_copies
    )

let () = main ()
