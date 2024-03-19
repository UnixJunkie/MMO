(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

(* compute the number of ROI entry and exit events *)

type state = In (* d < d_min *)
           | Between (* d_min <= d <= d_max *)
           | Out (* d > d_max *)

(* extract distance to ROI center from an energy logfile line *)
let get_dist line =
  float_of_string (S.cut_on_char ' ' 1 line) (* fields start at 0; unlike w/ awk *)

let main () =
  (* the two params could be changed on the CLI; those are default values *)
  let in_dist = 5.0 in
  let out_dist = 10.0 in
  let ene_fn = Sys.argv.(1) in
  let entries = ref 0 in
  let exits = ref 0 in
  let prev_d = ref 1000.0 in
  LO.with_in_file ene_fn (fun input ->
      (* skip comment line *)
      let _header = input_line input in
      try
        while true do
          let line = input_line input in
          let d = get_dist line in
          (if !prev_d > in_dist && d < in_dist then
             incr entries
          );
          (if !prev_d < out_dist && d > out_dist then
             incr exits
          );
          prev_d := d
        done
      with End_of_file -> ()
    );
  Printf.printf "%d %d\n" !entries !exits

let () = main ()
