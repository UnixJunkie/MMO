(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
 *
 * 3D grid stored as a 1D array, to store pre-calculated components of a FF *)

module V3 = Mmo.V3

open Printf

type t = BA1.ba1f (* grid data *)

(* store in a file *)
let to_ba1_file fn grid arr =
  let dims_fn = fn ^ ".dims" in
  let n = Grid.num_voxels grid in
  assert(n = BA1.dim arr);
  let fd = Unix.(openfile fn [O_RDWR; O_CREAT; O_TRUNC] 0o600) in
  (* create mmapped bigarray *)
  let dst =
    BA.array1_of_genarray
      (Unix.map_file fd BA.Float32 BA.c_layout true [|n|]) in
  (* copy existing bigarray to the (new) one mapped to file *)
  BA1.blit arr dst;
  Unix.close fd;
  (* store dims *)
  Log.info "creating %s" dims_fn;
  LO.lines_to_file dims_fn
    [sprintf "step: %g"  grid.step;
     sprintf "x_dim: %d" grid.x_dim;
     sprintf "y_dim: %d" grid.y_dim;
     sprintf "z_dim: %d" grid.z_dim]

let id x =
  x

(* retrieve grid step and dimensions *)
let parse_dims_file fn =
  LO.with_in_file fn (fun input ->
      let step  = Scanf.sscanf (input_line input) "step: %f"  id in
      let x_dim = Scanf.sscanf (input_line input) "x_dim: %d" id in
      let y_dim = Scanf.sscanf (input_line input) "y_dim: %d" id in
      let z_dim = Scanf.sscanf (input_line input) "z_dim: %d" id in
      (step, x_dim, y_dim, z_dim)
    )

let create grid =
  let n = Grid.num_voxels grid in
  let arr = BA1.create BA.Float32 BA.c_layout n in
  BA1.fill arr 0.0;
  arr

(* restore grid from a bigarray
   >1000 TIMES FASTER THAN UNMARSHAL *)
let of_ba1_file g fn =
  let n = Grid.num_voxels g in
  let fd = Unix.(openfile fn [O_RDONLY] 0o400) in
  (* populate bigarray from file *)
  let ba1 =
    BA.array1_of_genarray
      (Unix.map_file fd BA.Float32 BA.c_layout false [|-1|]) in
  assert(BA1.dim ba1 = n);
  Unix.close fd;
  ba1

(* integer grid coord *)
type ijk = { i: int ;
             j: int ;
             k: int }

(* when you need a value of the type, but will do nothing useful with it *)
let dummy =
  BA1.create BA.Float32 BA.c_layout 0

(* init grid cell w/ already computed 1D index *)
let init_idx arr idx x =
  (* initialization code is not supposed to overwrite grid cells; *)
  (* each grid cell is initialized only once *)
  assert(BA1.get arr idx = 0.0);
  BA1.set arr idx x

(* initialize grid cell *)
let init grid arr i j k x =
  init_idx arr Grid.(i + j * grid.x_dim + k * grid.xy_dim) x

let min_avg_max_sparse arr =
  BA1.(min arr, favg arr, max arr, sparsity arr)

(* --- 3D grid points / voxel ---

  h---g
e/--f/|
| d |/c
a/--b

 * tri-linear interpolation from the grid *)
let trilin grid arr (p: V3.t): float =
  Grid.(
    (* CODE HAND CHECKED SEVERAL TIMES - DO NOT TOUCH *)
    (* (i0, j0, k0): low voxel corner for [p] *)
    let i0 = int_of_float V3.(p.x *. grid.one_div_step) in
    let j0 = int_of_float V3.(p.y *. grid.one_div_step) in
    let k0 = int_of_float V3.(p.z *. grid.one_div_step) in
    (* (i1, j1, k1): high voxel corner for [p] *)
    let i1 = i0 + 1 in
    let j1 = j0 + 1 in
    let k1 = k0 + 1 in
    (* eprintf "%d %d %d %d %d %d\n" i0 j0 k0 i1 j1 k1; *)
    let j0x  = j0 * grid.x_dim  in
    let j1x  = j1 * grid.x_dim  in
    let k0xy = k0 * grid.xy_dim in
    let k1xy = k1 * grid.xy_dim in
    (* grid content at voxel corners *)
    (* try *)
    (* let a = (BA1.get g3d.arr (i0 + j0x + k0xy)) in
       let b = (BA1.get g3d.arr (i1 + j0x + k0xy)) in
       let c = (BA1.get g3d.arr (i1 + j1x + k0xy)) in
       let d = (BA1.get g3d.arr (i0 + j1x + k0xy)) in
       let e = (BA1.get g3d.arr (i0 + j0x + k1xy)) in
       let f = (BA1.get g3d.arr (i1 + j0x + k1xy)) in
       let g = (BA1.get g3d.arr (i1 + j1x + k1xy)) in
       let h = (BA1.get g3d.arr (i0 + j1x + k1xy)) in *)
    (* fp coordinates for low voxel corner *)
    let l = V3.{ x = A.unsafe_get grid.xs i0 ;
                 y = A.unsafe_get grid.ys j0 ;
                 z = A.unsafe_get grid.zs k0 } in
    (* weights; all in [0:1] *)
    let wl = V3.{ x = V3.(p.x -. l.x) *. grid.one_div_step ;
                  y = V3.(p.y -. l.y) *. grid.one_div_step ;
                  z = V3.(p.z -. l.z) *. grid.one_div_step } in
    let wh = V3.{ x = 1.0 -. wl.x ;
                  y = 1.0 -. wl.y ;
                  z = 1.0 -. wl.z } in
    (* volume weights for voxel corners; CONSTR: sum(weights) ~= 1.0 *)
    (* { v_a = (wh.x *. wh.y *. wh.z) ;
         v_b = (wl.x *. wh.y *. wh.z) ;
         v_c = (wl.x *. wl.y *. wh.z) ;
         v_d = (wh.x *. wl.y *. wh.z) ;
         v_e = (wh.x *. wh.y *. wl.z) ;
         v_f = (wl.x *. wh.y *. wl.z) ;
         v_g = (wl.x *. wl.y *. wl.z) ;
         v_h = (wh.x *. wl.y *. wl.z) } *)
    (* final interpolation *)
    ((BA1.unsafe_get arr (i0 + j0x + k0xy)) *. (wh.x *. wh.y *. wh.z) +.
     (BA1.unsafe_get arr (i1 + j0x + k0xy)) *. (wl.x *. wh.y *. wh.z) +.
     (BA1.unsafe_get arr (i1 + j1x + k0xy)) *. (wl.x *. wl.y *. wh.z) +.
     (BA1.unsafe_get arr (i0 + j1x + k0xy)) *. (wh.x *. wl.y *. wh.z) +.
     (BA1.unsafe_get arr (i0 + j0x + k1xy)) *. (wh.x *. wh.y *. wl.z) +.
     (BA1.unsafe_get arr (i1 + j0x + k1xy)) *. (wl.x *. wh.y *. wl.z) +.
     (BA1.unsafe_get arr (i1 + j1x + k1xy)) *. (wl.x *. wl.y *. wl.z) +.
     (BA1.unsafe_get arr (i0 + j1x + k1xy)) *. (wh.x *. wl.y *. wl.z))
    (* with exn -> *)
    (* begin *)
    (*   eprintf "%d %d %d %d %d %d\n" i0 j0 k0 i1 j1 k1; *)
    (*   raise exn *)
    (* end *)
  )

(* fast detection of some vdW atom clashes:
   any of the eight 3D neighbors clashing means vdW clash;
   valid test for an atom centered at [p] *)
let vdW_clash_OR grid bitmask p: bool =
  Grid.(
    (* (i0, j0, k0): low voxel corner for [p] *)
    let i0 = int_of_float V3.(p.x *. grid.one_div_step) in
    let j0 = int_of_float V3.(p.y *. grid.one_div_step) in
    let k0 = int_of_float V3.(p.z *. grid.one_div_step) in
    (* (i1, j1, k1): high voxel corner for [p] *)
    let i1 = i0 + 1 in
    let j1 = j0 + 1 in
    let k1 = k0 + 1 in
    (* eprintf "%d %d %d %d %d %d\n" i0 j0 k0 i1 j1 k1; *)
    let j0x  = j0 * grid.x_dim  in
    let j1x  = j1 * grid.x_dim  in
    let k0xy = k0 * grid.xy_dim in
    let k1xy = k1 * grid.xy_dim in
    (* test the eight 3D neighbors *)
    Bitv.get bitmask (i0 + j0x + k0xy) ||
    Bitv.get bitmask (i1 + j0x + k0xy) ||
    Bitv.get bitmask (i1 + j1x + k0xy) ||
    Bitv.get bitmask (i0 + j1x + k0xy) ||
    Bitv.get bitmask (i0 + j0x + k1xy) ||
    Bitv.get bitmask (i1 + j0x + k1xy) ||
    Bitv.get bitmask (i1 + j1x + k1xy) ||
    Bitv.get bitmask (i0 + j1x + k1xy)
  )

(* stricter version: all eight neighbors must be clashing *)
let vdW_clash_AND grid bitmask p: bool =
  Grid.(
    (* (i0, j0, k0): low voxel corner for [p] *)
    let i0 = int_of_float V3.(p.x *. grid.one_div_step) in
    let j0 = int_of_float V3.(p.y *. grid.one_div_step) in
    let k0 = int_of_float V3.(p.z *. grid.one_div_step) in
    (* (i1, j1, k1): high voxel corner for [p] *)
    let i1 = i0 + 1 in
    let j1 = j0 + 1 in
    let k1 = k0 + 1 in
    (* eprintf "%d %d %d %d %d %d\n" i0 j0 k0 i1 j1 k1; *)
    let j0x  = j0 * grid.x_dim  in
    let j1x  = j1 * grid.x_dim  in
    let k0xy = k0 * grid.xy_dim in
    let k1xy = k1 * grid.xy_dim in
    (* test the eight 3D neighbors *)
    Bitv.get bitmask (i0 + j0x + k0xy) &&
    Bitv.get bitmask (i1 + j0x + k0xy) &&
    Bitv.get bitmask (i1 + j1x + k0xy) &&
    Bitv.get bitmask (i0 + j1x + k0xy) &&
    Bitv.get bitmask (i0 + j0x + k1xy) &&
    Bitv.get bitmask (i1 + j0x + k1xy) &&
    Bitv.get bitmask (i1 + j1x + k1xy) &&
    Bitv.get bitmask (i0 + j1x + k1xy)
  )
