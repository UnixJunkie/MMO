(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

module V3 = Mmo.V3

open Printf

let () =
  let lig = Mol.receptor_of_pqrs_file "data/caffeine.pqrs" in
  (* Tx *)
  let fn = sprintf "/tmp/lig0.bild" in
  let white = Ptable.White_H in
  Mol.to_bild ~verbose:true ~style:white fn lig;
  eprintf "%s\n" fn;
  for i = 1 to 10 do
    Mol.translate_by lig (V3.make 10. 0. 0.);
    let fn = sprintf "/tmp/lig%d.bild" i in
    eprintf "%s\n" fn;
    Mol.to_bild ~verbose:true ~style:white fn lig
  done;
  for i = 11 to 20 do
    Mol.translate_by lig (V3.make 0. 10. 0.);
    let fn = sprintf "/tmp/lig%d.bild" i in
    eprintf "%s\n" fn;
    Mol.to_bild ~verbose:true ~style:white fn lig
  done;
  for i = 21 to 30 do
    Mol.translate_by lig (V3.make 0. 0. 10.);
    let fn = sprintf "/tmp/lig%d.bild" i in
    eprintf "%s\n" fn;
    Mol.to_bild ~verbose:true ~style:white fn lig
  done;
  (* (\* let n = int_of_string Sys.argv.(1) in *\) *)
  (* let n = 1_000_000 in *)
  (* let lig1 = Mol.copy lig  in *)
  (* let lig2 = Mol.copy lig1 in *)
  (* let delta = V3.make 0.001 0.001 0.001 in *)
  (* let dt_in_place, _ = *)
  (*   Utls.time_it (fun () -> *)
  (*       for _i = 1 to n do *)
  (*         Mol.translate_by lig1 delta *)
  (*       done *)
  (*     ) in *)
  (* eprintf "in-place dt: %f\n" dt_in_place; *)
  (* eprintf "in-place center: "; *)
  (* V3.xyz_dump stderr (Mol.geometric_center lig1); *)
  (* eprintf "\n"; *)
  (* let dt_copy, translated_copy = *)
  (*   Utls.time_it (fun () -> *)
  (*       let res = ref lig2 in *)
  (*       for _i = 1 to n do *)
  (*         res := Mol.copy_translate_by !res delta *)
  (*       done; *)
  (*       !res *)
  (*     ) in *)
  (* eprintf "copy dt: %f\n" dt_copy; *)
  (* eprintf "copy center: "; *)
  (* V3.xyz_dump stderr (Mol.geometric_center translated_copy); *)
  (* eprintf "\n" *)
