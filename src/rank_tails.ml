(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

(* rank-order molecules from a list of scores files;
 * such scores files are produced by exhaustive docking w/ -top option *)

(* parse all scores from left tail of the distribution for given molecule *)
let scores_from_file fn =
  let res =
    A.of_list (LO.map fn (fun line -> Scanf.sscanf line "%f" (fun x -> x))) in
  assert(A.is_sorted res);
  res

(* retrieve molecule name from .scores filename *)
let mol_name_from_scores_fn fn =
  let replaced, mol2_fn = S.replace ~str:fn ~sub:".scores" ~by:".mol2" in
  assert(replaced);
  let mol = Mol2.read_one_from_file mol2_fn in
  Mol2.get_name mol

(* sort score tails from most to least dominating *)
let compare_scores s1 s2 =
  let left, right = Utls.dominate s1 s2 in
  if left > right then -1
  else if left = right then 0
  else 1

let main () =
  let score_files_fn = Sys.argv.(1) in
  let scores_files = LO.lines_of_file score_files_fn in
  let scores = L.map scores_from_file scores_files in
  let names = L.map mol_name_from_scores_fn scores_files in
  let labels = L.map (String.starts_with ~prefix:"active") names in
  let score_label_names = A.of_list (L.combine3 scores labels names) in
  (* we want decreasingly sorted output *)
  A.sort (fun (s1, _l1, _n1) (s2, _l2, _n2) -> compare_scores s2 s1)
    score_label_names;
  (* dump them out *)
  let n = A.length score_label_names in
  A.iteri (fun i (scores, label, name) ->
      let min, avg, max = A.min_avg_max scores in
      Printf.printf "%d %d %g %g %g %s\n"
        (* high rank -> high negative "docking score" for bin/dock_ROC.sh *)
        (i - n) (Utls.int_of_bool label) min avg max name
    ) score_label_names

let () = main ()
