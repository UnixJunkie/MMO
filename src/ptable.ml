(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

(* periodic table *)

module IMap = BatMap.Int
module ISet = BatSet.Int
module SMap = BatMap.String

let all_supported_anums_l = [1; 6; 7; 8; 9; 12; 15; 16; 17; 35; 53]
let all_supported_anums = ISet.of_list all_supported_anums_l
let anum_max = ISet.max_elt all_supported_anums

let symbols = A.make 119 "X" (* X stands for unknown *)
(* only fill in the ones we support *)
let () =
  begin
    symbols.(1)  <- "H" ;
    symbols.(6)  <- "C" ;
    symbols.(7)  <- "N" ;
    symbols.(8)  <- "O" ;
    symbols.(9)  <- "F" ;
    symbols.(12) <- "Mg" ;
    symbols.(15) <- "P" ;
    symbols.(16) <- "S" ;
    symbols.(17) <- "Cl";
    symbols.(35) <- "Br";
    symbols.(53) <- "I" ;
  end

let vdW_radii = A.make 119 nan
(* only fill in the ones we support
   Python code to get those:
---
import rdkit
from rdkit import Chem
ptable = Chem.GetPeriodicTable()
ptable.GetRvdw(6)
--- *)
let () =
  begin
    vdW_radii.(1)  <- 1.2 ;
    vdW_radii.(6)  <- 1.7 ;
    vdW_radii.(7)  <- 1.6 ;
    vdW_radii.(8)  <- 1.55;
    vdW_radii.(9)  <- 1.5 ;
    vdW_radii.(12) <- 2.2 ;
    vdW_radii.(15) <- 1.95;
    vdW_radii.(16) <- 1.8 ;
    vdW_radii.(17) <- 1.8 ;
    vdW_radii.(35) <- 1.9 ;
    vdW_radii.(53) <- 2.1 ;
  end

(* d_ij < a * (r_vdW_i + r_vdW_j) is inspired by Majeux PROTEINS 1999:
 * PROTEINS: Structure, Function, and Genetics 37:88–105 (1999)
 * Exhaustive Docking of Molecular Fragments
 * With Electrostatic Solvation
 * Nicolas Majeux, Marco Scarsi, Joannis Apostolakis,
 * Claus Ehrhardt and Amedeo Caflisch *)
let vdW_clash_alpha = 0.8 (* somewhat arbitrary *)

(* vdW clash parameters squared: given two atomic numbers, return
   x2 = (0.8 * (vdW(i) + vdW(j)))^2
   a faster clash test is: (d_ij)^2 < x2 *)
let vdW_clash_params2: (float array) array = A.make (anum_max + 1) [||]
let () = (* initialize this 2D array *)
  begin
    L.iter (fun anum_i ->
        let vdW_i = vdW_radii.(anum_i) in
        let row = A.make (anum_max + 1) nan in
        vdW_clash_params2.(anum_i) <- row;
        L.iter (fun anum_j ->
            let vdW_j = vdW_radii.(anum_j) in
            let x = vdW_clash_alpha *. (vdW_i +. vdW_j) in
            row.(anum_j) <- x *. x
          ) all_supported_anums_l
      ) all_supported_anums_l
  end

(* min/max radii you can encounter in a molecule *)
let vdW_min = A.min vdW_radii
let vdW_max = A.max vdW_radii

let sym2anu =
  SMap.of_seq
    (BatSeq.of_list
       (* currently supported elements *)
       [("C" ,  6);
        ("H" ,  1);
        ("N" ,  7);
        ("O" ,  8);
        ("P" , 15);
        ("S" , 16);
        ("F" ,  9);
        ("Cl", 17);
        ("Br", 35);
        ("I" , 53);
        ("Mg", 12)])

(* atomic number for the given chemical element *)
let anum_of_symbol (elt: string): int =
  try SMap.find elt sym2anu
  with Not_found ->
    let _ = Log.fatal "Ptable.anum_of_symbol: unsupported element: %s" elt in
    exit 1

let symbol_of_anum (anum: int): string =
  A.unsafe_get symbols anum

let vdW_radius_of_anum (anum: int): float =
  A.unsafe_get vdW_radii anum

let standard_colors_no_H =
  [( 6, "grey"  );
   ( 7, "blue"  );
   ( 8, "red"   );
   (15, "purple");
   (16, "yellow");
   ( 9, "green" );
   (17, "green" );
   (35, "green" );
   (53, "green" )]

type coloring_scheme = No_color
                     | White_H
                     | Pink_H
                     | Cyan_H
                     | Orange_H

(* CPK atom coloring scheme.
   At least for the ligand, this will be helpful *)
let anum2cpk =
  IMap.of_seq
    (BatSeq.of_list
       ((1, "white") :: standard_colors_no_H))

let anum2cpk_pink_H =
  IMap.of_seq
    (BatSeq.of_list
       ((1, "pink") :: standard_colors_no_H))

let anum2cpk_cyan_H =
  IMap.of_seq
    (BatSeq.of_list
       ((1, "cyan") :: standard_colors_no_H))

let anum2cpk_orange_H =
  IMap.of_seq
    (BatSeq.of_list
       ((1, "orange") :: standard_colors_no_H))

let color_of_anum (anum: int): string =
  IMap.find anum anum2cpk

let color_of_anum_pink_H (anum: int): string =
  IMap.find anum anum2cpk_pink_H

let color_of_anum_orange_H (anum: int): string =
  IMap.find anum anum2cpk_orange_H

let color_of_anum_cyan_H (anum: int): string =
  IMap.find anum anum2cpk_cyan_H

exception Unsupported_atom of string

(* atomic number for Sybyl atom types (found in MOL2 files) *)
let anum_of_mol2_type = function
  (* mol2 type; comment *)
  | "H"     (* hydrogen *)
  | "H.t3p" (* hydrogen in Transferable intermolecular Potential (TIP3P) water model *)
  | "H.spc" (* hydrogen in Single Point Charge (SPC) water model *) -> 1
  | "Li"    (* lithium *) -> 3
  | "C.3"   (* carbon sp3 *)
  | "C.2"   (* carbon sp2 *)
  | "C.1"   (* carbon sp *)
  | "C.ar"  (* carbon aromatic *)
  | "C.cat" (* carbocation (C + ) used only in a guadinium group *) -> 6
  | "N.3"   (* nitrogen sp3 *)
  | "N.2"   (* nitrogen sp2 *)
  | "N.1"   (* nitrogen sp *)
  | "N.ar"  (* nitrogen aromatic *)
  | "N.am"  (* nitrogen amide *)
  | "N.pl3" (* nitrogen trigonal planar *)
  | "N.4"   (* nitrogen sp3 positively charged *) -> 7
  | "O.3"   (* oxygen sp3 *)
  | "O.2"   (* oxygen sp2 *)
  | "O.co2" (* oxygen in carboxylate and phosphate groups *)
  | "O.spc" (* oxygen in Single Point Charge (SPC) water model *)
  | "O.t3p" (* oxygen in Transferable Intermolecular Potential (TIP3P) water model *) -> 8
  | "F"     (* fluorine *) -> 9
  | "Na"    (* sodium *) -> 11
  | "Mg"    (* magnesium *) -> 12
  | "Al"    (* aluminum *) -> 13
  | "Si"    (* silicon *) -> 14
  | "P.3"   (* phosphorous sp3 *) -> 15
  | "S.3"   (* sulfur sp3 *)
  | "S.2"   (* sulfur sp2 *)
  | "S.O"   (* sulfoxide sulfur *) | "S.o"
  | "S.O2"  (* sulfone sulfur *) | "S.o2" -> 16
  | "Cl"    (* chlorine *) -> 17
  | "K"     (* potassium *) -> 19
  | "Ca"    (* calcium *) -> 20
  | "Cr.th" (* chromium (tetrahedral) *)
  | "Cr.oh" (* chromium (octahedral) *) -> 24
  | "Mn"    (* manganese *) -> 25
  | "Fe"    (* iron *) -> 26
  | "Co.oh" (* cobalt (octahedral) *) -> 27
  | "Cu"    (* copper *) -> 29
  | "Zn"    (* zinc *) -> 30
  | "Se"    (* selenium *) -> 34
  | "Br"    (* bromine *) -> 35
  | "Mo"    (* molybdenum *) -> 42
  | "Sn"    (* tin *) -> 50
  | "I"     (* iodine *) -> 53
  | other ->
    (* "LP"    (* lone pair *)
       "Du"    (* dummy atom *)
       "Du.C"  (* dummy carbon *)
       "Any"   (* any atom *)
       "Hal"   (* halogen *)
       "Het"   (* heteroatom = N, O, S, P *)
       "Hev"   (* heavy atom (non hydrogen) *) *)
    raise (Unsupported_atom other)
