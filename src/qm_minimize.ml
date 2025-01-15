(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
 *
 * Minimize ligand conformer using ANI-2 QM approximation *)

module Math = Mmo.Math

open Printf

(* get ready to use some Python bindings *)
let () = Py.initialize ~version:3 ()

(* per rotatable bond *)
type discretization = Degrees of float (* rotational step size *)
                    | Steps of int (* number of steps for one turn *)

let radians_of_discretization disc =
  Math.to_radian (match disc with
      | Degrees d -> d
      | Steps s -> (360.0 /. (float s))
    )

module Defaults = struct
  let nprocs = 1
  let chunk_size = 50
  let degrees = 30
  let steps = 360 / 30
  let max_rbonds = 7
  let confs_cap = 10_000
end

(* the type of things that will be emitted by the parmap demux function *)
type conf_in = { index: int;   (* molecule index in input file *)
                 offset: int;  (* 1st conformer to generate *)
                 length: int } (* number of conformers to generate *)

let rec demux mol_mol2_ligands alpha csize lig_i conf_i =
  let num_ligs = A.length mol_mol2_ligands in
  if !lig_i = num_ligs then
    (* all ligands were processed *)
    raise Parany.End_of_input
  else
    let lig, _mol2 = mol_mol2_ligands.(!lig_i) in
    let max_conf = Mol.count_conformers lig alpha in
    (if !conf_i = 0 then
       Log.info "%s: %d max confs" (Mol.get_name lig) max_conf
    );
    if !conf_i >= max_conf then
      (incr lig_i;
       conf_i := 0;
       demux mol_mol2_ligands alpha csize lig_i conf_i)
    else
      let res = { index = !lig_i;
                  offset = !conf_i;
                  length = min csize (max_conf - !conf_i) } in
      conf_i := !conf_i + csize;
      res

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  Log.(set_prefix_builder short_prefix_builder);
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s\n  \
              [-i <ligands.mol2>]: input file\n  \
              [-o <QM_minimized.mol2>]: output file\n  \
              [-np <int>]: nprocs (default=1)\n  \
              [--cap <int>]: \"cap\" (max number of conformers per molecule; default=%d)\n  \
              [-c <int>]: chunk size (number of conformers scored on one core)\n  \
              [-d <float>: RotBond step in degrees \
              (e.g. %d; incompatible w/ -s)\n  \
              [-r <int>]: max RotBond allowed (default=%d)\n  \
              [-s <int>]: RotBond steps for one turn \
              (e.g. %d; incompatible w/ -d)\n  \
              [--seed <int>] RNG seed\n  \
              [-v]: verbose mode\n"
       Sys.argv.(0)
       Defaults.confs_cap
       Defaults.degrees
       Defaults.max_rbonds
       Defaults.steps;
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let nprocs = CLI.get_int_def ["-np"] args Defaults.nprocs in
  Flags.verbose := CLI.get_set_bool ["-v"] args;
  let maybe_degrees = CLI.get_float_opt ["-d"] args in
  let maybe_rot_steps = CLI.get_int_opt ["-s"] args in
  let _cap = CLI.get_int_def ["--cap"] args Defaults.confs_cap in
  let csize = CLI.get_int_def ["-c"] args Defaults.chunk_size in
  let max_rot_bonds = CLI.get_int_def ["-r"] args Defaults.max_rbonds in
  let _rng = match CLI.get_int_opt ["--seed"] args with
    | None -> RNG.(make (entropy_160b ()))
    | Some seed -> RNG.make [|seed|] in
  CLI.finalize (); (* ------------------------------------------------------ *)
  (if !Flags.verbose then
     Log.(set_log_level DEBUG)
  );
  let rot_step = match (maybe_degrees, maybe_rot_steps) with
    | (None, None) -> failwith "qm_minimize: use either -d OR -s"
    | (Some d, None) -> radians_of_discretization (Degrees d)
    | (None, Some s) -> radians_of_discretization (Steps s)
    | (Some _, Some _) -> failwith "qm_minimize: use either -d OR -s, NOT both" in
  let mol_mol2_ligands =
    let mol_ligands = A.of_list (Mol.ligands_of_mol2_file nprocs input_fn) in
    let mol2_ligands = Mol2.read_all input_fn in
    let all = A.combine mol_ligands mol2_ligands in
    A.filter (fun (mol, _mol2) ->
        let name = Mol.get_name mol in
        let rbonds = Mol.num_rbonds mol in
        if rbonds <= max_rot_bonds then
          true
        else
          let () = Log.warn "%s: ignored (%d RotBonds)" name rbonds in
          false
      ) all in
  let num_ligs = A.length mol_mol2_ligands in
  let rank2name = Ht.create num_ligs in
  A.iteri (fun i (mol, _mol2) ->
      Ht.add rank2name i (Mol.get_name mol)
    ) mol_mol2_ligands;
  let ht = Ht.create num_ligs in
  Log.info "%s: %d molecules" input_fn num_ligs;
  LO.with_out_file output_fn (fun output ->
      let lig_i = ref 0 in
      let conf_i = ref 0 in
      Parany.run nprocs
        ~csize:1 (* csize=1: DO NOT CHANGE; the demux fun handles csize *)
        ~demux:(fun () -> demux mol_mol2_ligands rot_step csize lig_i conf_i)
        ~work:(fun { index; offset; length } ->
            let lig, mol2 = mol_mol2_ligands.(index) in
            let name = Mol.get_name lig in
            let topo_dists_ge3 = Mol.get_topo_dists_ge3 lig in
            (if offset = 0 && Mol.ligand_atoms_clash_fast topo_dists_ge3 lig then
               Log.warn "%s: init conf vdW clash" name
            );
            let anums = Mol.ene_QM_anums lig in
            let num_rbonds = Mol.num_rbonds lig in
            let max_conf = Mol.count_conformers lig rot_step in
            (* generate conformers, score them, return lowest energy one *)
            let i = ref 0 in
            let curr_conf = ref lig in
            let curr_ene = ref (Mol.ene_QM_score anums lig) in
            let steps_per_rot_bond = int_of_float (Math.two_pi /. rot_step) in
            while !i < length && offset + !i < max_conf do
              let conf_id = offset + !i in
              (* DEBUG *)
              Log.debug "%s: conf %d" name conf_id;
              let rbonds = Mol.decode_rbonds_conf num_rbonds steps_per_rot_bond rot_step conf_id in
              let conf = Mol.rotate_bonds lig rbonds in
              (if not (Mol.ligand_atoms_clash_fast topo_dists_ge3 conf) then
                 let ene = Mol.ene_QM_score anums conf in
                 if ene < !curr_ene then
                   (curr_ene := ene;
                    curr_conf := conf)
              );
              incr i
            done;
            (!curr_ene, Mol.update_mol2 mol2 !curr_conf)
          )
        ~mux:(fun (curr_ene, curr_mol2) ->
            let name = Mol2.get_name curr_mol2 in
            let prev_ene, _prev_mol2 =
              try Ht.find ht name
              with Not_found -> (infinity, curr_mol2) in
            if curr_ene < prev_ene then
              let () = Log.info "%s: %.3f" name curr_ene in
              Ht.replace ht name (curr_ene, curr_mol2)
          );
      (* write lowest energy conformer of each molecule to disk
         while preserving molecules input order *)
      for i = 0 to num_ligs - 1 do
        let name = Ht.find rank2name i in
        let ene, mol2 = Ht.find ht name in
        Mol2.output_named_one output (sprintf "%s_E=%f" name ene) mol2
      done
    )

let () = main ()
