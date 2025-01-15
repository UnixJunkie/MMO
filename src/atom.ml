(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

module V3 = Mmo.V3

type t = { idx: int; pos: V3.t }

let create idx pos =
  { idx; pos }

let dist a1 a2 =
  V3.dist a1.pos a2.pos

let get_idx a =
  a.idx
