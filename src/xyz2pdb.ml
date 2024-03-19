(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

(* convert a .xyz file into a .pdb; all molecules are converted *)

open Printf

(* read one molecule from a .xyz file *)
let xyz_read_one counter input =
  let n = int_of_string (input_line input) in
  (* get all the lines, don't parse them yet;
     one comment line + n atom lines *)
  let i = !counter in
  incr counter;
  (i, A.init (n + 1) (fun _i -> input_line input))

let first_time = ref true

let xyz2pdb pdb_out (frame_num, lines) =
  if !first_time && frame_num > 9999 then
    begin
      first_time := false;
      Log.error "xyz2pdb.xyz2pdb: frame_num over 9999"
    end;
  fprintf pdb_out "MODEL\t%d\n" frame_num;
  (* (\* actual PDB format for model number line *\) *)
  (* highest MODEL number supported by PDB format *)  
  (* assert(frame_num <= 9999); *)
  (* fprintf pdb_out "MODEL     %04d\n" frame_num; *)
  let n = A.length lines in
  (* don't allow crazy big ligand *)
  assert(n < 100);
  (* example HETATM line
     |0    5    10   5    20   5    30   5    40   5    50   5    60   5    70   5  |
     ^HETATM   10  N9  CFF A 330       6.240 -34.112 -34.550  1.00147.23           N$ *)
  for i = 1 to n - 1 do
    Scanf.sscanf lines.(i) "%s@ %f %f %f" (fun elt x y z ->
        fprintf pdb_out "HETATM% 5d%4s  LIG A   1    %8.3f%8.3f%8.3f  1.00  0.00%12s\n"
          i elt x y z elt
      )
  done;
  fprintf pdb_out "ENDMDL\n"

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  Log.(set_prefix_builder short_prefix_builder);
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s\n  \
              [-i <filename.xyz>]: input\n  \
              [-o <filename.pdb>]: output\n"
       Sys.argv.(0);
     exit 1);
  let input_fn  = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  CLI.finalize (); (* ------------------------------------------------------ *)
  LO.with_infile_outfile input_fn output_fn (fun input output ->
      let count = ref 0 in
      try
        while true do
          let mol = xyz_read_one count input in
          xyz2pdb output mol
        done
      with End_of_file ->
        Log.info "read: %d" !count
    )

let () = main ()
