(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

module Math = Mmo.Math
module Rot = Mmo.Rot
module V3 = Mmo.V3

open Printf

let () =
  let dt = 5.0 in
  let lig = Mol.receptor_of_pqrs_file "data/cartesian_frame.pqrs" in
  Mol.center lig;
  let lig0 = Mol.copy lig in
  (* Rx *)
  let fn = sprintf "/tmp/lig0.bild" in
  let white = Ptable.White_H in
  Mol.to_bild ~verbose:true ~style:white fn lig;
  eprintf "%s\n" fn;
  for i = 1 to 10 do
    Mol.translate_by lig (V3.make dt 0. 0.);
    Mol.inplace_rotate lig (Rot.rx (0.1 *. Math.two_pi));
    let fn = sprintf "/tmp/lig%d.bild" i in
    eprintf "%s\n" fn;
    Mol.to_bild ~verbose:true ~style:white fn lig;
  done;
  (* RMSD between lig0 and lig10 should be ~=0 *)
  let lig10 = Mol.copy lig in
  Mol.center lig10;
  let rmsd0_10 = Mol.rmsd lig0 lig10 in
  eprintf "RMSD lig0 lig10: %f\n" rmsd0_10;
  (* Ry *)
  for i = 11 to 20 do
    Mol.translate_by lig (V3.make 0. dt 0.);
    Mol.inplace_rotate lig (Rot.ry (0.1 *. Math.two_pi));
    let fn = sprintf "/tmp/lig%d.bild" i in
    eprintf "%s\n" fn;
    Mol.to_bild ~verbose:true ~style:white fn lig;
  done;
  (* RMSD between lig10 and lig20 should be ~=0 *)
  let lig20 = Mol.copy lig in
  Mol.center lig20;
  let rmsd10_20 = Mol.rmsd lig10 lig20 in
  eprintf "RMSD lig10 lig20: %f\n" rmsd10_20;
  (* Rz *)
  for i = 21 to 30 do
    Mol.translate_by lig (V3.make 0. 0. dt);
    Mol.inplace_rotate lig (Rot.rz (0.1 *. Math.two_pi));
    let fn = sprintf "/tmp/lig%d.bild" i in
    eprintf "%s\n" fn;
    Mol.to_bild ~verbose:true ~style:white fn lig;
  done;
  (* RMSD between lig20 and lig30 should be ~=0 *)
  let lig30 = Mol.copy lig in
  Mol.center lig30;
  let rmsd20_30 = Mol.rmsd lig20 lig30 in
  eprintf "RMSD lig20 lig30: %f\n" rmsd20_30
