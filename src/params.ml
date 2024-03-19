(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

(* some simulation defaults *)

(* in radian, whole molecule *)
let max_rot = Math.to_radian 15.0

(* in Angstrom, whole molecule *)
let max_trans = 0.15

(* in radian, one rotatable bond only *)
let max_rbond_rot = Math.to_radian 5.0

(* in radian, one rotatable bond only *)
let max_rbond_flip = Math.pi

(* WARNING: grid_step=0.1 comsumes all RAM; 0.25 is too slow to init. *)
let grid_step = 0.5 (* A *)

(* maximal energy value allowed in grid;
   keep high enough for --no-vdw-clash to work *)
let max_E = 100_000.0

(* how many frames before updating some counters *)
let block_size = 100

(* how many frames before we try to flip a rotatable bond
 * once per block <=> max introduced error in acceptance rate=1% *)
let rbond_flip_block = ref block_size
