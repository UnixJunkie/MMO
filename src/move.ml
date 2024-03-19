(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

(* ligand rigid body moves *)

(* ROTATION ---------------------------------------------------------------- *)

(* uniform random angle in [-max_rot..max_rot] *)
let rand_angle rng max_rot =
  (Random.State.float rng (2.0 *. max_rot)) -. max_rot

(* 1) rand. choose Ox, Oy or Oz
 * 2) rotate by theta(radian) on that axis
 * 3) update current rotation *)
let rand_rot rng dr rot =
  (* in [-max_rot..max_rot] *)
  let theta = rand_angle rng dr in
  (* in [0..2] *)
  let axis = Random.State.int rng 3 in
  let rotate_by =
    if      axis = 0 then Rot.rx theta
    else if axis = 1 then Rot.ry theta
    else if axis = 2 then Rot.rz theta
    else assert(false) in
  (* rotations are applied right to left *)
  Rot.mult rotate_by rot

(* completely random initial rotation *)
let rand_rot_full rng =
  rand_rot rng Math.pi (
      rand_rot rng Math.pi (
          rand_rot rng Math.pi (Rot.id ())
        )
    )

(* TRANSLATION ------------------------------------------------------------- *)

(* 1) choose a random 3D direction in ([-1..1], [-1..1], [-1..1])
 * 2) scale it by translation step [dt]
 * 3) add to current pos. *)
let rand_trans rng dt pos =
  (* random direction *)
  let dir =
    V3.make
      ((Random.State.float rng 2.0) -. 1.0)
      ((Random.State.float rng 2.0) -. 1.0)
      ((Random.State.float rng 2.0) -. 1.0) in
  (* new position *)
  V3.add pos (V3.mult dir dt)
