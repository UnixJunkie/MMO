(* optimize weights for the rescoring function *)

module SL = struct
  type t = bool * float (* (label, score) *)
  let create label score =
    (label, score)
  let get_score (_, s) = s
  let get_label (l, _) = l
  (* to do a decreasing sort of a score labels list *)
  let high_score_first_cmp (_, s1) (_, s2) =
    BatFloat.compare s2 s1
end

module ROC = Cpm.MakeROC.Make(SL)

type components = { name: string;
                    inter_elec: float;
                    inter_vdW: float;
                    desolv_prot: float;
                    desolv_lig: float;
                    rbonds: int;
                    heavy_atoms: int }

let parse_detailed_scores_line l =
  Scanf.sscanf l "%s %f %f %f %f %d %d"
    (fun name inter_elec inter_vdW desolv_prot desolv_lig rbonds heavy_atoms ->
       { name;
         inter_elec; inter_vdW;
         desolv_prot; desolv_lig;
         rbonds;
         heavy_atoms }
    )

let load_detailed_scores fn =
  match LO.lines_of_file fn with
  | [] -> assert(false)
  | _header :: lines -> A.of_list (L.map parse_detailed_scores_line lines)

let is_active comps =
  S.starts_with_stdlib ~prefix:"active" comps.name

(* ----- Scoring Functions ----- *)

let inter_desolv_rbonds a b c comps =
  SL.create
    (is_active comps)
    (* invert sign of docking scores *)
    (-.(a *. (comps.inter_elec +. comps.inter_vdW) +.
        b *. (comps.desolv_prot +. comps.desolv_lig) +.
        c *. (float comps.rbonds)))

let inter_elec_inter_vdW_desolv a b c comps =
  SL.create
    (is_active comps)
    (* invert sign of docking scores *)
    (-.(a *. comps.inter_elec +.
        b *. comps.inter_vdW  +.
        c *. (comps.desolv_prot +. comps.desolv_lig)))

let inter_elec_inter_vdW a b comps =
  SL.create
    (is_active comps)
    (* invert sign of docking scores *)
    (-.(a *. comps.inter_elec +.
        b *. comps.inter_vdW))

let inter_desolv a b comps =
  SL.create
    (is_active comps)
    (* invert sign of docking scores *)
    (-.(a *. (comps.inter_elec +. comps.inter_vdW) +.
        b *. (comps.desolv_prot +. comps.desolv_lig)))

let inter_desolvProt_desolveLig a b c comps =
  SL.create
    (is_active comps)
    (* invert sign of docking scores *)
    (-.(a *. (comps.inter_elec +. comps.inter_vdW) +.
        b *. comps.desolv_prot +.
        c *. comps.desolv_lig))

let strain_UB heavy_atoms =
  max 0.0 (0.3 *. float (heavy_atoms - 10))

let interE_desolv_strainUB a b c comps =
  SL.create
    (is_active comps)
    (* invert sign of docking scores *)
    (-.(a *. (comps.inter_elec +. comps.inter_vdW) +.
        b *. (comps.desolv_prot +. comps.desolv_lig) +.
        c *. (strain_UB comps.heavy_atoms)))

(* Ligand Efficiency *)
let full_Einter_LE comps =
  SL.create
    (is_active comps)
    (* invert sign of docking scores *)
    (-.(1.0 *. (comps.inter_elec +. comps.inter_vdW) /. (float comps.heavy_atoms)))

let weight_fun_of_string = function
  | "inter_desolv_rbonds" ->         inter_desolv_rbonds
  | "inter_elec_inter_vdW_desolv" -> inter_elec_inter_vdW_desolv
  | "inter_desolvProt_desolveLig" -> inter_desolvProt_desolveLig
  | "interE_desolv_strainUB" ->      interE_desolv_strainUB
  | "full_Einter_LE" ->              (fun _a _b _c -> full_Einter_LE)
  | "inter_elec_inter_vdW" ->        (fun a b _c -> inter_elec_inter_vdW a b)
  | "inter_desolv" ->                (fun a b _c -> inter_desolv a b)
  | other ->
    failwith ("Optim_weights.weight_fun_of_string: unsupported: " ^ other)

let reweight_then_ROC_AUC w_fun (detailed_scores: components array): float =
  let reweighted = A.map w_fun detailed_scores in
  ROC.auc_a reweighted

let parse_filenames fns =
  A.of_list (S.split_on_char ',' fns)

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let weight_fun = weight_fun_of_string Sys.argv.(1) in
  let filenames = parse_filenames Sys.argv.(2) in
  let num_score_files = A.length filenames in
  Log.info "%d score files" num_score_files;
  let detailed_scores = A.map load_detailed_scores filenames in
  let a = ref 0.0 in
  let b = ref 0.0 in
  let step = 0.05 in
  while !a <= 1.0 do
    while !b <= 1.0 -. !a do
      let c = 1.0 -. (!a +. !b) in
      let roc_AUCs =
        A.map
          (reweight_then_ROC_AUC
             (weight_fun !a !b c)
          ) detailed_scores in
      Printf.printf "%f %f %f %f\n" !a !b c (A.favg roc_AUCs);
      b := !b +. step
    done;
    b := 0.0;
    a := !a +. step
  done

let () = main ()
