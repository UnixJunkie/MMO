(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

(* rename molecules from a .mol2 file; duplicate names will be
 * appended "_%d" (the number of times this name was already seen) *)

open Printf

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  Log.(set_prefix_builder short_prefix_builder);
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s\n  \
              [-i <ligands.mol2>]: input file\n  \
              [-o <renamed.mol2>]: output file\n"
       Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  CLI.finalize (); (* ------------------------------------------------------ *)
  let counter = ref 0 in
  let ht = Ht.create 10007 in
  LO.with_infile_outfile input_fn output_fn (fun input out ->
      try
        while true do
          let mol2_lines = Mol2.get_mol2_block_exn input in
          incr counter;
          assert(mol2_lines.(0) = Mol2.molecule_tag);
          let name' = mol2_lines.(1) in
          let count = Ht.find_default ht name' 0 in
          let name =
            if count = 0 then
              (* first time seen *)
              (Ht.add ht name' 1;
               name')
            else
              (Ht.replace ht name' (count + 1);
               sprintf "%s_%d" name' (count + 1)) in
          fprintf out "%s\n" mol2_lines.(0);
          fprintf out "%s\n" name;
          A.iteri (fun i line ->
              if i >= 2 then
                fprintf out "%s\n" line
            ) mol2_lines
        done
      with End_of_file ->
        Log.info "read %d molecules from %s" !counter input_fn
    )

let () = main ()
