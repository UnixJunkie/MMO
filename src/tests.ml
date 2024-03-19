(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

let cos_pi_div_4 = cos (0.25 *. Math.pi)

let () =
  Log.(set_log_level INFO);
  (* dihedral angle test
     j-k
     | |
     i l *)
  let i = V3.make 0. 0. 0. in
  let j = V3.make 0. 1. 0. in
  let k = V3.make 1. 1. 0. in

  let l0 = V3.make 1. 2. 0. in
  let l1 = V3.make 1. (1. +. cos_pi_div_4) cos_pi_div_4 in
  let l2 = V3.make 1. 1. 1. in
  let l3 = V3.make 1. (1. -. cos_pi_div_4) cos_pi_div_4 in
  let l4 = V3.make 1. 0. 0. in
  let l5 = V3.make 1. (1. -. cos_pi_div_4) (-.cos_pi_div_4) in
  let l6 = V3.make 1. 1. (-.1.) in
  let l7 = V3.make 1. (1. +. cos_pi_div_4) (-.cos_pi_div_4) in

  let phi0 = Math.dihedral_angle i j k l0 in
  let phi1 = Math.dihedral_angle i j k l1 in
  let phi2 = Math.dihedral_angle i j k l2 in
  let phi3 = Math.dihedral_angle i j k l3 in
  let phi4 = Math.dihedral_angle i j k l4 in
  let phi5 = Math.dihedral_angle i j k l5 in
  let phi6 = Math.dihedral_angle i j k l6 in
  let phi7 = Math.dihedral_angle i j k l7 in

  Log.info "phi0: %.2f" (Math.rad2deg *. phi0);
  Log.info "phi1: %.2f" (Math.rad2deg *. phi1);
  Log.info "phi2: %.2f" (Math.rad2deg *. phi2);
  Log.info "phi3: %.2f" (Math.rad2deg *. phi3);
  Log.info "phi4: %.2f" (Math.rad2deg *. phi4);
  Log.info "phi5: %.2f" (Math.rad2deg *. phi5);
  Log.info "phi6: %.2f" (Math.rad2deg *. phi6);
  Log.info "phi7: %.2f" (Math.rad2deg *. phi7);

  assert (Math.approx_equal 0.001 (1.     *. Math.pi) phi0);
  assert (Math.approx_equal 0.001 (-0.75  *. Math.pi) phi1);
  assert (Math.approx_equal 0.001 (-0.5   *. Math.pi) phi2);
  assert (Math.approx_equal 0.001 (-0.25  *. Math.pi) phi3);
  assert (Math.approx_equal 0.001 (0.0    *. Math.pi) phi4);
  assert (Math.approx_equal 0.001 (0.25   *. Math.pi) phi5);
  assert (Math.approx_equal 0.001 (0.5    *. Math.pi) phi6);
  assert (Math.approx_equal 0.001 (0.75   *. Math.pi) phi7)
