
(* cleanup mol2 file generated by CCDC GOLD *)

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  Log.(set_prefix_builder short_prefix_builder);
  let input_fn  = Sys.argv.(1) in
  let tmp_fn = Fn.temp_file ~temp_dir:"/tmp" "mol2clean_" ".mol2" in
  let output_fn = Sys.argv.(2) in
  (* remove mol2 sections added by CCDC GOLD *)
  LO.with_infile_outfile input_fn tmp_fn
    (fun input output ->
       Mol2.filter_gold_output_mol2 input output);
  (* remove lone pairs added by CCDC GOLD *)
  let cleaned_mols = Mol2.read_all tmp_fn in
  Sys.remove tmp_fn;
  Mol2.output_all output_fn cleaned_mols

let () = main ()
