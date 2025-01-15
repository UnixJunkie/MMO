(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
 *
 * Implementation of Super-Fibonacci Spirals:
 * Fast, Low-Discrepancy Sampling of SO(3)
 * Marc Alexa; Proceedings of CVPR, 2022, pp. 8291-8300
 * https://openaccess.thecvf.com/content/CVPR2022/html/\
 * Alexa_Super-Fibonacci_Spirals_Fast_Low-\
 * Discrepancy_Sampling_of_SO3_CVPR_2022_paper.html *)

module Math = Mmo.Math
module Quat = Mmo.Quat
module Rot = Mmo.Rot

(* both constants are from p4 of the paper *)
let phi = sqrt 2.0
let psi = 1.533751168755204288118041

(* core function from Algorithm 1 on p4 of
   doc/Alexa_CVPR_2022_SO3sampling.pdf *)
let super_fibonacci n i =
  let s = (float i) +. 0.5 in
  let t = s /. n in
  let d = Math.two_pi *. s in
  let c_r = sqrt t in
  let c_R = sqrt (1.0 -. t) in
  let alpha = d /. phi in
  let beta  = d /. psi in
  Quat.create
    (c_r *. sin alpha)
    (c_r *. cos alpha)
    (c_R *. sin beta)
    (c_R *. cos beta)

let sample n =
  Array.init n (super_fibonacci (float n))

let rotations (n: int): Rot.t array =
  A.map (fun quat ->
      let axis, angle = Quat.to_axis_angle quat in
      Rot.of_axis_angle axis angle
    ) (sample n)
