(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

let square x =
  x *. x

let pow3 x =
  x *. x *. x

let pow6 x =
  pow3 (x *. x)

let geo_mean x y =
  sqrt (x *. y)

let shift_12A d =
  if d < 12.0 then
    square (1.0 -. (square (d /. 12.0)))
  else 0.0
