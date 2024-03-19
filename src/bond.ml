(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

type t = { left:  int ; (* begin atom index *)
           right: int } (* end atom index *)

(* tip of the rotation axis *)
let right_end b =
  b.right

let create left right =
  { left; right }

let dummy =
  create (-1) (-1)

let to_pair b =
  (b.left, b.right)
