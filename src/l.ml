(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

include BatList

let combine3 l1 l2 l3 =
  map2 (fun (x, y) z ->
      (x, y, z)
    ) (combine l1 l2) l3
