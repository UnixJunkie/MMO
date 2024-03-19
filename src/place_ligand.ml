(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
 *
 * place ligand at given x y z a b g [rbond0 [rbond1 ...]]
 * (a, b, g) are the Cartesian angles alpha, beta, gamma
 * if rbonds are not given; only apply rigid body movement *)

let main () =
  let argc = A.length Sys.argv in
  (if argc = 1 then
     let () =
       (*                      0  1          2           3 4 5 6 7 8 9 ... *)
       Printf.eprintf "usage:\n%s input.mol2 output.mol2 x y z a b g [rbonds]\n"
         Sys.argv.(0) in
     exit 1
  );
  let mol2_in_fn  = Sys.argv.(1) in (* SINGLE MOLECULE input file *)
  let mol2_out_fn = Sys.argv.(2) in (* output file *)
  (* 6 rigid body DOFs *)
  let x = float_of_string Sys.argv.(3) in
  let y = float_of_string Sys.argv.(4) in
  let z = float_of_string Sys.argv.(5) in
  let a = float_of_string Sys.argv.(6) in
  let b = float_of_string Sys.argv.(7) in
  let g = float_of_string Sys.argv.(8) in
  (* rbonds *)
  let rbonds_given = argc - 9 in
  Log.info "rbonds given: %d" rbonds_given;  
  (* read in mol2 *)
  let mol2 = match LO.with_in_file mol2_in_fn Mol2.read_one with
    | None -> failwith ("place_ligand: cannot parse " ^ mol2_in_fn)
    | Some x -> x in
  (* create Mol.t *)
  let tmp_pqrs_out_fn = Fn.temp_file "place_lig_" ".pqrs" in
  Utls.mol2pqrs 1 mol2_in_fn tmp_pqrs_out_fn;
  let centered_lig: Mol.t = match Mol.ligands_of_pqrs_file tmp_pqrs_out_fn with
    | [x] -> (Mol.center x; x)
    | [] -> failwith ("place_ligand: no ligand in " ^ tmp_pqrs_out_fn)
    | _ -> failwith ("place_ligand: several ligands in " ^ tmp_pqrs_out_fn) in
  Sys.remove tmp_pqrs_out_fn; (* cleanup *)
  (* check number of RotBonds *)
  let num_rbonds = Mol.num_rbonds centered_lig in
  let rbonds =
    if rbonds_given = 0 then
      (* only rigid body move *)
      A.make num_rbonds 0.0
    else
      A.init rbonds_given (fun i ->
          float_of_string Sys.argv.(9 + i)
        ) in
  (if rbonds_given > 0 && rbonds_given <> num_rbonds then
     let () =
       Log.fatal "place_ligand: %d rbond values but mol has %d"
         rbonds_given num_rbonds in
     exit 1
  );
  (* apply transform *)
  let cfg = Optim.create x y z a b g rbonds in
  let mol = Optim.apply_config centered_lig cfg in
  let mol2' = Mol.update_mol2 mol2 mol in
  (* save as mol2 *)
  Mol2.write_one_to_file mol2_out_fn mol2'

let () = main ()
