(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

include Vector3

(* alias *)
let create = make

let cube side =
  make side side side

(* to show a coordinate in the logs *)
let short_str v =
  Printf.sprintf "%g %g %g" v.x v.y v.z

(* for a .xyz trajectory file *)
let xyz_dump out v =
  Printf.fprintf out " %.3f %.3f %.3f\n" v.x v.y v.z

(* [dist2]: squared distance,
   saves one call to sqrt compared to [dist] below *)
let dist2 u v =
  let dx = u.x -. v.x in
  let dy = u.y -. v.y in
  let dz = u.z -. v.z in
  dx *. dx +. dy *. dy +. dz *. dz
[@@inline]

(* the one in library vector3 should be optimized *)
let dist u v =
  sqrt (dist2 u v)
[@@inline]
