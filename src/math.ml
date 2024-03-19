(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

(* cos and sin of an angle *)
type cosin = { cos: float ;
               sin: float }

let cos_sin theta =
  { cos = cos theta ;
    sin = sin theta }

let pi = 4.0 *. (atan 1.0)

let two_pi = 2.0 *. pi

let deg2rad = pi /. 180.0

let rad2deg = 180.0 /. pi

let to_radian x =
  x *. deg2rad

let is_NaN x = match classify_float x with
  | FP_nan -> true
  | _ -> false

(*
(* Blondel and Karplus, JoCC 1996, signed dihedral angle *)
let dihedral_angle i j k l =
  V3.(let f = diff i j in
      let g = diff j k in
      let h = diff l k in
      let a = cross f g in
      let b = cross h g in
      let mag_a = mag a in
      let mag_b = mag b in
      let mag_g = mag g in
      let cos_phi = (dot a b) /. (mag_a *. mag_g) in
      let sin_phi = (dot (cross b a) g) /. (mag_a *. mag_b *. mag_g) in
      atan2 sin_phi cos_phi)
*)

(* wikipedia; second formula using atan2; simpler than Blondel and Karplus *)
let dihedral_angle i j k l =
  let open V3 in
  let u1 = diff j i in
  let u2 = diff k j in
  let u2u3 = cross u2 (diff l k) in
  atan2
    (dot (mult u1 (mag u2)) u2u3)
    (dot (cross u1 u2)      u2u3)

let approx_equal tolerance reference actual =
  abs_float (reference -. actual) <= tolerance

(* avoid division by 0 in case of too small atomic distance *)
let non_zero_dist x =
  if x < 0.01 then
    0.01 (* trick used in D. Horvath's S4MPLE docking program *)
  else
    x
