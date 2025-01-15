(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
 *
 * Useful constants *)

module Math = Mmo.Math

(* water radius *)
let r_H2O = 1.4 (* A *)

let charged_cutoff = 12.0 (* A *)

let charged_cutoff_squared = charged_cutoff *. charged_cutoff

(* https://doi.org/10.1021/ct400065j *)
(* protein electric permitivity *)
let epsilon_prot = 4.0 (* most common value; Scarsi, Majeux and Caflisch PROTEINS 1999 *)

(* water electric permitivity *)
let epsilon_HOH = 78.5 (* at 20 Celsius degrees; Scarsi, Majeux and Caflisch PROTEINS 1999 *)

let room_temp_K = 293.15 (* 20 Celsius in Kelvin *)

(* Boltzmann constant in kcal/(mol*K) *)
let kB = 0.0019872041

let room_temp_beta = 1.0 /. (kB *. room_temp_K)

(* Majeux, Scarsi and Caflisch PROTEINS 2001; left constant from eq. (2) *)
let desolvation = (1./.epsilon_prot -. 1./.epsilon_HOH) /. (8.0 *. Math.pi)

(* cf. https://en.wikipedia.org/wiki/Hartree *)
let hartree_to_kcal_per_mol = 627.5094740631
