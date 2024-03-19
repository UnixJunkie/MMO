(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

type t = { w: float ;
           x: float ;
           y: float ;
           z: float }

let create w x y z =
  { w; x; y; z }

let to_string q =
  Printf.sprintf "{%g;%g;%g;%g}"
    q.w q.x q.y q.z

let of_string s =
  Scanf.sscanf s "{%f;%f;%f;%f}" create

(* https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation *)
let of_axis_angle (v3, theta) =
  Math.(
    let t = cos_sin (theta *. 0.5) in
    V3.{ w = t.cos;
         x = v3.x *. t.sin;
         y = v3.y *. t.sin;
         z = v3.z *. t.sin }
  )

(* https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
 * paragraph: Recovering the axis-angle representation *)
let to_axis_angle ({ w; x; y; z}: t): V3.t * float =
  let mag = sqrt (x*.x +. y*.y +. z*.z) in
  let axis = V3.make (x /. mag) (y /. mag) (z /. mag) in
  let theta = 2.0 *. (atan2 mag w) in
  (axis, theta)
