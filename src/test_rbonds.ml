(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

module Math = Mmo.Math

(* test that all rot bonds can be rotated +/- *)

let () =
  let ligs = Mol.ligands_of_pqrs_file "data/GBA_lig1.pqrs" in
  assert(L.length ligs = 1);
  let lig = L.hd ligs in
  let num_rbonds = Mol.num_rbonds lig in
  let frame = ref 0 in
  let delta = (2. *. Math.pi) /. 20.0 in
  (* turn each in one direction *)
  for i = 0 to num_rbonds - 1 do
    for _j = 1 to 20 do
      Mol.rotate_bond lig i delta;
      Mol.xyz_dump stdout !frame "" lig;
      incr frame
    done
  done;
  (* then each in the other direction *)
  for i = 0 to num_rbonds - 1 do
    for _j = 1 to 20 do
      Mol.rotate_bond lig i (-.delta);
      Mol.xyz_dump stdout !frame "" lig;
      incr frame
    done
  done
