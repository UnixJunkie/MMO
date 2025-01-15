(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
 *
 * Spherical Region of Interest;
 * probably the binding-site center plus 5 to 10A away.
 * CCDC GOLD uses up to 10A away from the binding-site center. *)

module V3 = Mmo.V3

open Printf

type sphere = { c: V3.t;
                out_r: float;
                out_r2: float }

let create c out_r =
  { c; out_r; out_r2 = out_r *. out_r }

(* WARNING: the bild file [fn] is in original PDB coordinates
 *          will need to be translated to sim. box *)
let from_bild fn =
  let ok_lines = LO.filter fn (fun l -> S.starts_with l ".sphere ") in
  match ok_lines with
  | [line] ->
    let x, y, z, out_r =
      Scanf.sscanf line ".sphere %f %f %f %f"
        (fun x y z r -> (x, y, z, r)) in
    create (V3.make x y z) out_r
  | _ ->
    let () = Log.fatal "ROI.from_bild: several sphere lines in: %s" fn in
    exit 1

let get_center s =
  s.c

let get_out_radius s =
  s.out_r

(* to save the one after translation to sim. box *)
let to_bild fn s =
  LO.lines_to_file fn
    [".transparency 0.8";
     ".color blue";
     sprintf ".sphere %s %g" (V3.short_str s.c) s.out_r]

(* should be more visible in chimera compared to the version
   produced by [to_bild] above *)
let dotted_to_bild fn s =
  let dot_surfaced_sphere =
    A.map V3.of_triplet (A.of_list (Genspir.genspir_cartesian 100)) in
  let surfaced_sphere_at_0 =
    let r = get_out_radius s in
    A.map (fun dot -> V3.mult dot r) dot_surfaced_sphere in
  let dots =
    let c = get_center s in
    A.map (V3.add c) surfaced_sphere_at_0 in
  LO.with_out_file fn (fun out ->
      fprintf out ".color yellow\n";
      A.iter (fun dot ->
          fprintf out ".sphere %s 0.1\n" (V3.short_str dot)
        ) dots
    )

let is_inside s p =
  V3.dist2 s.c p < s.out_r2

let is_outside s p =
  V3.dist2 s.c p > s.out_r2

let translate s v =
  { s with c = V3.add s.c v }

(* (x_min, x_max,
 *  y_min, y_max,
 *  z_min, z_max) *)
let get_bounds roi =
  let x, y, z = V3.to_triplet (get_center roi) in
  let radius = get_out_radius roi in
  (x -. radius, x +. radius,
   y -. radius, y +. radius,
   z -. radius, z +. radius)
