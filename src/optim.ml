(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

(* what the optimizer has to deal with for a given ligand *)
(*
type config = { x: float; (* translation *)
                y: float;
                z: float;
                rx: float; (* rotation *)
                ry: float;
                rz: float;
                rbonds: float array (* RotBonds *) }
*)
type config = float array

let create x y z rx ry rz rbonds =
  A.append [|x; y; z; rx; ry; rz|] rbonds

(* For derivative-free local-optimization algorithms:
 * hint for initial step size in each dimension
 * [dt]: translational stepsize in A
 * [dr]: rotational stepsize in radian for any angle (RotBonds included) *)
let step_sizes c dt dr =
  let n = A.length c in
  A.append [|dt; dt; dt; dr; dr; dr|] (A.make (n - 6) dr)

let pi = Math.pi

let get_min_bounds roi tweak_rbonds lig =
  let x_min, _,
      y_min, _,
      z_min, _ = ROI.get_bounds roi in
  let rot_bonds =
    if tweak_rbonds then
      (* flexible ligand *)
      Mol.num_rbonds lig
    else
      (* rigid-ligand *)
      0 in
  A.append
    [|x_min; y_min; z_min|]
    (A.make (3 + rot_bonds) (-.pi))

let get_max_bounds roi tweak_rbonds lig =
  let _, x_max,
      _, y_max,
      _, z_max = ROI.get_bounds roi in
  let rot_bonds =
    if tweak_rbonds then
      (* flexible ligand *)
      Mol.num_rbonds lig
    else
      (* rigid-ligand *)
      0 in
  A.append
    [|x_max; y_max; z_max|]
    (A.make (3 + rot_bonds) pi)

let apply_config centered_lig conf =
  (* rotate bonds *)
  let lig = Mol.copy centered_lig in
  let n = A.length conf in
  (* 0..5 are the rigid body DOFs *)
  for i = 6 to n - 1 do
    Mol.rotate_bond lig (i - 6) conf.(i)
  done;
  (* forbid too elongated ligand;
   * the scoring function will have to return the worst
   * score in case exception is raised here *)
  Mol.check_elongation_exn lig Const.charged_cutoff;
  (* rotate then translate lig *)
  Mol.rotate_then_translate_copy
    lig
    (Rot.r_xyz conf.(3) conf.(4) conf.(5))
    (V3.make   conf.(0) conf.(1) conf.(2))

let rot_trans_of conf =
  (Rot.r_xyz conf.(3) conf.(4) conf.(5),
   V3.make   conf.(0) conf.(1) conf.(2))

let to_string conf =
  (* use '%f' for high precision *)
  Utls.string_of_array conf " " (Printf.sprintf "%f")
