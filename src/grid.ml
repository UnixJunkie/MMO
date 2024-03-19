(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

open Printf

(* lowest corner of this grid is at (0, 0, 0) *)

type t = { step: float ; (* Cubic grid: dx = dy = dz *)
           one_div_step: float ; (* (1 / step) *)
           xs: float array ;
           ys: float array ;
           zs: float array ;
           x_dim: int ;
           y_dim: int ;
           z_dim: int ;
           xy_dim: int }

let get_dims g =
  (g.x_dim,
   g.y_dim,
   g.z_dim)

let num_voxels g =
  g.x_dim * g.y_dim * g.z_dim

let voxel_volume g =
  g.step *. g.step *. g.step

(* bitmask for grid [g] *)
let get_bitmask g =
  Bitv.create (num_voxels g) false

(* the grid must englobe the whole box; a little bit more is OK *)
let num_steps dx length =
  int_of_float (ceil (length /. dx))

let from_box step bbox =
  let one_div_step = 1.0 /. step in
  let x_dim' = num_steps step (Bbox.dx bbox) in
  let y_dim' = num_steps step (Bbox.dy bbox) in
  let z_dim' = num_steps step (Bbox.dz bbox) in
  let x_dim = x_dim' + 1 in
  let y_dim = y_dim' + 1 in
  let z_dim = z_dim' + 1 in
  let xy_dim = x_dim * y_dim in
  let xs = A.of_list (L.frange 0.0 `To (step *. float x_dim') x_dim) in
  let ys = A.of_list (L.frange 0.0 `To (step *. float y_dim') y_dim) in
  let zs = A.of_list (L.frange 0.0 `To (step *. float z_dim') z_dim) in
  { step; one_div_step; x_dim; y_dim; z_dim; xy_dim; xs; ys; zs }

let dummy = { step = nan ;
              one_div_step = nan ;
              x_dim = 0 ;
              y_dim = 0 ;
              z_dim = 0 ;
              xy_dim = 0 ;
              xs = [||] ;
              ys = [||] ;
              zs = [||] }

(* draw a uniform random grid coordinate *)
let rand_coord rng g =
  (Random.State.int rng g.x_dim,
   Random.State.int rng g.y_dim,
   Random.State.int rng g.z_dim)

(* random 3D point at an exact grid coordinate *)
let rand_lattice_point rng grid =
  let ix, iy, iz = rand_coord rng grid in
  V3.make
    ((float ix) *. grid.step)
    ((float iy) *. grid.step)
    ((float iz) *. grid.step)

(* rand. 3D point not necessarily on the latice *)
let rand_point rng g =
  V3.make
    (Random.State.float rng (g.step *. float g.x_dim))
    (Random.State.float rng (g.step *. float g.y_dim))
    (Random.State.float rng (g.step *. float g.z_dim))

(* lower corner coordinate (nearest to origin)
   of the grid voxel containing [p] *)
let coord_of_point g p =
  V3.(int_of_float ((p.x -. g.xs.(0)) /. g.step),
      int_of_float ((p.y -. g.ys.(0)) /. g.step),
      int_of_float ((p.z -. g.zs.(0)) /. g.step))

let str_peek g =
  Printf.sprintf "xs: %g %g %g... %g ys: %g %g %g... %g zs: %g %g %g... %g"
    g.xs.(0) g.xs.(1) g.xs.(2) g.xs.(g.x_dim - 1)
    g.ys.(0) g.ys.(1) g.ys.(2) g.ys.(g.y_dim - 1)
    g.zs.(0) g.zs.(1) g.zs.(2) g.zs.(g.z_dim - 1)

let idx_of_ijk grid i j k =
  i + j * grid.x_dim + k * grid.xy_dim

let ijk_of_idx grid idx =
  let k = idx / grid.xy_dim in
  let j = (idx - k * grid.xy_dim) / grid.x_dim in
  let i = idx - (k * grid.xy_dim + j * grid.x_dim) in
  (i, j, k)

let to_bild fn g =
  LO.with_out_file fn (fun out ->
      let m, n, p = get_dims g in
      (* grid points color *)
      fprintf out ".color orange\n";
      for i = 0 to m - 1 do
        let x = g.xs.(i) in
        for j = 0 to n - 1 do
          let y = g.ys.(j) in
          for k = 0 to p - 1 do
            let z = g.zs.(k) in
            (* only one point every 10 *)
            if (i mod 10) = 0 &&
               (j mod 10) = 0 &&
               (k mod 10) = 0 then
              fprintf out ".sphere %g %g %g 0.3\n" x y z
          done
        done
      done
    )
