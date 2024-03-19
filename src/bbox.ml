(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
 *
 * a bounding box *)

open Printf

type t = { low:  V3.t ; (* lower corner *)
           high: V3.t ; (* high corner *)
           dims: V3.t } (* dimensions *)

let get_dims b =
  b.dims

let get_low b =
  b.low

let get_high b =
  b.high

(* from two vectors *)
let create_2v low high =
  let dims = V3.diff high low in
  assert(V3.(dims.x > 0.0 && dims.y > 0.0 && dims.z > 0.0));
  { low; high; dims }

(* from six floats *)
let create_6f x0 y0 z0 x1 y1 z1 =
  let low  = V3.make x0 y0 z0 in
  let high = V3.make x1 y1 z1 in
  create_2v low high

let dx b =
  b.dims.x

let dy b =
  b.dims.y

let dz b =
  b.dims.z

let is_inside box p =
  V3.(box.low.x  <= p.x && box.low.y  <= p.y && box.low.z  <= p.z &&
      box.high.x >= p.x && box.high.y >= p.y && box.high.z >= p.z)

(* uniform random draw of a 3D point inside box [b] *)
let rand_point_inside rng b =
  V3.make
    (b.low.x +. (Random.State.float rng b.dims.x))
    (b.low.y +. (Random.State.float rng b.dims.y))
    (b.low.z +. (Random.State.float rng b.dims.z))

(* dump to .bild file for chimera *)
let to_bild fn color box =
  Log.info "box in %s" fn;
  LO.with_out_file fn (fun out ->
      fprintf out ".transparency 0.8\n";
      fprintf out ".color %s\n" color;
      V3.(fprintf out ".box %g %g %g %g %g %g\n"
            box.low.x  box.low.y  box.low.z
            box.high.x box.high.y box.high.z)
    )

(* move [p] inside the [box], if needed (periodic boundary condition) *)
let rebox box p =
  (* in the book "computer simulation of liquids" they give
   * a higher performance version, not using if-then-else *)
  if is_inside box p then
    p (* no need *)
  else
    (* we assume the point _just_ escaped the box
       <=> is in a direct neighbor box *)
    let open V3 in
    let p' = V3.diff p box.low in
    let dx = if p'.x < 0.0 then box.dims.x
             else if p'.x > box.dims.x then -.box.dims.x
             else 0.0 in
    let dy = if p'.y < 0.0 then box.dims.y
             else if p'.y > box.dims.y then -.box.dims.y
             else 0.0 in
    let dz = if p'.z < 0.0 then box.dims.z
             else if p'.z > box.dims.z then -.box.dims.z
             else 0.0 in
    V3.add p (V3.make dx dy dz)
