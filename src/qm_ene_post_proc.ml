(* convert absolute lig_E_intra_QM into relative ligand strain energy *)

let split_name_conf s =
  let l = S.split_on_char '_' s in
  let l' = L.rev l in
  match l' with
  | [] -> assert(false)
  | x :: xs -> (S.concat "_" (L.rev xs), (* mol_name *)
                int_of_string x) (* conf_id *)

let parse_line l =
  Scanf.sscanf l "%s %f" (fun name ene ->
      let name, conf = split_name_conf name in
      (name, conf, ene)
    )

let process_key_value name conf_energies =
  let min_E = L.min (L.map snd conf_energies) in
  (name,
   L.map (fun (conf, ene) -> (conf, ene -. min_E)) conf_energies)

let main () =
  let input_fn = Sys.argv.(1) in
  let ht = Ht.create (LO.length input_fn) in
  LO.iter input_fn (fun line ->
      let name, conf, ene = parse_line line in
      let prev = Ht.find_default ht name [] in
      Ht.replace ht name ((conf, ene) :: prev)
    );
  Ht.iter (fun name conf_energies ->
      let name, conf_strains = process_key_value name conf_energies in
      L.iter (fun (conf, rel_ene) ->
          Printf.printf "%s_%d %g\n" name conf rel_ene
        ) conf_strains
    ) ht

let () = main ()
