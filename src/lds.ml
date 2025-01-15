(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

open Printf

module Math = Mmo.Math
module Rdkit = Rdkit.Rdkit
module Rot = Mmo.Rot
module TopK = Cpm.TopKeeper
module V3 = Mmo.V3

(* get ready to use some Python bindings *)
let () = Py.initialize ~version:3 ()

(* the protein is moved at:
   (margin + bbox.dx/2, margin + bbox.dy/2, margin + bbox.dz/2)
   so that all its coordinates are positive *)
let preprocess_protein prot_bild_fn box_bild_fn margin prot =
  let box = Mol.bounding_box prot in
  let box_dims = Bbox.get_dims box in
  Log.info "tight protein box(dx, dy, dz) = %s" (V3.short_str box_dims);
  let margins = V3.cube (2. *. margin) in
  let sim_box_dims = V3.add box_dims margins in
  let sim_box = Bbox.create_2v V3.origin sim_box_dims in
  let sim_box_center = V3.mult sim_box_dims 0.5 in
  let old_prot_center = Mol.get_center prot in
  Log.(info "%s %s"
         (color green "old prot center:")
         (V3.short_str old_prot_center));
  Mol.translate_to prot sim_box_center;
  Log.(info "%s %s"
         (color green "new prot center:")
         (V3.short_str sim_box_center));
  Mol.to_bild ~verbose:true ~frame:true ~style:Ptable.No_color
    prot_bild_fn prot;
  Log.info "simulation box(dx, dy, dz) = %s"
    (V3.short_str (Bbox.get_dims sim_box));
  Bbox.to_bild box_bild_fn "blue" sim_box;
  sim_box

(* the ligand is placed at the origin;
   so that it is ready for subsequent rotation then translation *)
let preprocess_ligand lig_bild_fn box_bild_fn lig_i lig =
  Mol.translate_to lig V3.origin;
  let box = Mol.bounding_box lig in
  let box_dims = Bbox.get_dims box in
  Log.info "l%d bbox(dx, dy, dz) = %s" lig_i (V3.short_str box_dims);
  Mol.to_bild ~verbose:true ~style:Ptable.White_H lig_bild_fn lig;
  Bbox.to_bild box_bild_fn "red" box

type rotation = Rand_rot
              | No_rot (* frozen orientation *)
              | Fixed_rot of (int * float) (* (axis, angle);
                                              axis in 0,1,2; angle(radian) *)

type translation = Rand_trans
                 | No_trans (* frozen position *)
                 | Fixed_trans of (int * float) (* (axis, trans) axis in 0,1,2;
                                                   trans(Angstrom) *)

let int_float_of_string s =
  Scanf.sscanf s "%d:%f" (fun i x -> (i, x))

(* [temp_K]: temperature in Kelvin *)
let beta temp_K =
  1.0 /. (Const.kB *. temp_K)

let boltzmann_criterion rng beta prev_e curr_e =
  (* always allow energy to decrease *)
  (curr_e <= prev_e) ||
  (* accept increase depending on dice roll and deltaE *)
  ((Random.State.float rng 1.0) < exp (-.(curr_e -. prev_e) *. beta))

let compute_bitmask_brute prot grid =
  let mask = Grid.get_bitmask grid in
  let xdim, ydim, zdim = Grid.get_dims grid in
  let atom_coords = Mol.get_all_atom_coords prot in
  (* sequential version *)
  for i = 0 to xdim - 1 do
    let x = grid.xs.(i) in
    for j = 0 to ydim - 1 do
      let y = grid.ys.(j) in
      for k = 0 to zdim - 1 do
        let z = grid.zs.(k) in
        if Mol.any_atom_nearer atom_coords (V3.create x y z) Const.charged_cutoff then
          let idx = Grid.idx_of_ijk grid i j k in
          Bitv.set mask idx true
      done
    done
  done;
  mask

(* set of grid points which are near enough a protein atom
 * to have non null E_inter values *)
let bitmask_whole_protein nprocs prot_bst grid =
  let mask = Grid.get_bitmask grid in
  let xdim, ydim, zdim = Grid.get_dims grid in
  Log.info "(xdim, ydim, zdim): %d*%d*%d = %#d"
    xdim ydim zdim (Grid.num_voxels grid);
  (* (\* sequential version *\) *)
  (* for i = 0 to xdim - 1 do *)
  (*   let x = grid.xs.(i) in *)
  (*   for j = 0 to ydim - 1 do *)
  (*     let y = grid.ys.(j) in *)
  (*     for k = 0 to zdim - 1 do *)
  (*       let z = grid.zs.(k) in *)
  (*       let query = Mol.Atom.create (-1) (V3.create x y z) in *)
  (*       let _nn, dist = Mol.BST.nearest_neighbor query prot_bst in *)
  (*       if dist < Const.charged_cutoff then *)
  (*         let idx = i + j * grid.x_dim + k * grid.xy_dim in *)
  (*         Bitv.set mask idx true *)
  (*     done *)
  (*   done *)
  (* done; *)
  let i' = ref 0 in
  let demux () =
    if !i' = xdim then
      raise Parany.End_of_input
    else
      let res = !i' in
      incr i';
      res in
  let work i =
    let res = ref [] in
    let x = grid.xs.(i) in
    for j = 0 to ydim - 1 do
      let y = grid.ys.(j) in
      for k = 0 to zdim - 1 do
        let z = grid.zs.(k) in
        let query = Atom.create (-1) (V3.create x y z) in
        let _nn, dist = Mol.BST.nearest_neighbor query prot_bst in
        if dist < Const.charged_cutoff then
          res := (j, k) :: !res (* set bit *)
      done
    done;
    (i, !res) in
  let mux (i, jk_l) =
    L.iter (fun (j, k) ->
        let idx = Grid.idx_of_ijk grid i j k in
        Bitv.set mask idx true
      ) jk_l in
  Parany.run nprocs ~demux ~work ~mux;
  mask

(* set bitmask inside atom's volume to boolean [b] *)
let atom_bitmask_set grid mask xyz radius b =
  let (i, j, k) = Grid.coord_of_point grid xyz in
  (* radius upper bound in number of grid cells *)
  let r_steps = int_of_float (ceil (radius /. grid.step)) in
  let i_min = i - r_steps in
  let j_min = j - r_steps in
  let k_min = k - r_steps in
  let i_max = i + r_steps in
  let j_max = j + r_steps in
  let k_max = k + r_steps in
  let r2 = radius *. radius in
  (* just test bits inside the smallest cube englobing this atom *)
  for i = i_min to i_max do
    let x = grid.xs.(i) in
    for j = j_min to j_max do
      let y = grid.ys.(j) in
      for k = k_min to k_max do
        if V3.(dist2 xyz (make x y grid.zs.(k))) < r2 then
          Bitv.set mask (i + j * grid.x_dim + k * grid.xy_dim) b
      done
    done
  done

(* set of grid points inside the first solvent shell surrounding [mol] *)
let first_solvent_shell grid mol =
  let mask = Grid.get_bitmask grid in
  let atoms = Mol.get_all_atom_coords mol in
  let radii = Mol.get_all_atom_radii mol in
  (* for each atom, SET bits inside r_vdW + r_H2O *)
  A.iter2 (fun xyz radius ->
      atom_bitmask_set grid mask xyz (radius +. Const.r_H2O) true
    ) atoms radii;
  (* for each atom, UNSET bits inside r_vdW *)
  A.iter2 (fun xyz radius ->
      atom_bitmask_set grid mask xyz radius false
    ) atoms radii;
  mask

(* this bitmask can be used for clash detection w/ protein atoms *)
let vdW_volume grid mol =
  let mask = Grid.get_bitmask grid in
  let n = Mol.num_atoms mol in
  (* for each atom, SET bits inside r_vdW *)
  for i = 0 to n - 1 do
    let xyz   = Mol.get_xyz    mol i in
    let r_vdW = Mol.get_radius mol i in
    atom_bitmask_set grid mask xyz r_vdW true
  done;
  mask

type desolvation = { prot: float;
                     lig: float }

(* compute the desolvation penalty of each voxel inside the protein's first
   shell of water (only inside the ROI though);
   implementation of formula (2) p258 of Majeux, Scarsi and Caflisch PROTEINS 2001 *)
let protein_desolv roi grid prot_bst prot_solvent_shell prot =
  let voxel_vol = Grid.voxel_volume grid in
  let n = Grid.num_voxels grid in
  assert(n = Bitv.length prot_solvent_shell);
  let res = A.make n 0.0 in
  (* for each protein solvent shell voxel *)
  let count = ref 0 in
  Bitv.iteri_true (fun idx ->
      (* this voxel's x,y,z *)
      let x_p =
        let i, j, k = Grid.ijk_of_idx grid idx in
        V3.make grid.xs.(i) grid.ys.(j) grid.zs.(k) in
      if ROI.is_inside roi x_p then
        begin
          let neighbors = Mol.(BST.neighbors (Atom.create (-1) x_p)
                                 Const.charged_cutoff prot_bst) in
          L.iter (fun prot_atom ->
              let prot_i = Atom.get_idx prot_atom in
              let x_j = Mol.get_xyz    prot prot_i in
              let q_j = Mol.get_charge prot prot_i in
              let d2 = V3.dist2 x_p x_j in
              let x = q_j /. d2 in
              res.(idx) <- res.(idx) +. (x *. x)
            ) neighbors;
          res.(idx) <- Const.desolvation *. (voxel_vol *. res.(idx));
          incr count
        end
    ) prot_solvent_shell;
  Log.info "Lds.protein_desolv: %d voxels" !count;
  res

(* implementation of formula (2) p258 of Majeux, Scarsi and Caflisch PROTEINS 2001 *)
let desolvation_penalty grid prot_desolv_contribs prot_solvent_shell _prot lig =
  let voxel_vol = Grid.voxel_volume grid in
  let lig_solvent_shell = first_solvent_shell grid lig in
  let desolvated = Bitv.bw_and prot_solvent_shell lig_solvent_shell in
  (* ligand desolvation *)
  let lig_desolv = ref 0.0 in
  let lig_atoms   = Mol.get_all_atom_coords lig in
  let lig_charges = Mol.get_all_charges     lig in
  (* protein desolvation *)
  let prot_desolv = ref 0.0 in
  let count = ref 0 in
  (* for each desolvated voxel *)
  Bitv.iteri_true (fun idx ->
      (* this voxel's x,y,z *)
      let x_p =
        let i, j, k = Grid.ijk_of_idx grid idx in
        V3.make grid.xs.(i) grid.ys.(j) grid.zs.(k) in
      (* ligand contrib *)
      A.iter2 (fun x_j q_j ->
          let d2 = V3.dist2 x_p x_j in
          (* hard cutoff *)
          if d2 < Const.charged_cutoff_squared then
            let x = q_j /. d2 in
            lig_desolv := !lig_desolv +. (x *. x)
        ) lig_atoms lig_charges;
      (* ligand contrib *)
      prot_desolv := !prot_desolv +. (A.unsafe_get prot_desolv_contribs idx);
      incr count
    ) desolvated;
  Log.info "Lds.desolv_contribs: %d voxels" !count;
  { prot = !prot_desolv;
    lig = Const.desolvation *. (voxel_vol *. !lig_desolv) }

let bitmask_ROI_only nprocs roi grid =
  let c = ROI.get_center roi in
  let r2 =
    (* we want ligand's E_inter to be able to become null out of the ROI *)
    let r = ROI.get_out_radius roi +. (Const.charged_cutoff *. 2.0) in
    r *. r in
  let mask = Grid.get_bitmask grid in
  let xdim, ydim, zdim = Grid.get_dims grid in
  Log.info "(xdim, ydim, zdim): %d*%d*%d = %#d"
    xdim ydim zdim (Grid.num_voxels grid);
  let i' = ref 0 in
  let demux () =
    if !i' = xdim then
      raise Parany.End_of_input
    else
      let res = !i' in
      incr i';
      res in
  let work i =
    let res = ref [] in
    let x = grid.xs.(i) in
    for j = 0 to ydim - 1 do
      let y = grid.ys.(j) in
      for k = 0 to zdim - 1 do
        (* grid point inside "augmented" ROI *)
        if V3.dist2 c (V3.create x y grid.zs.(k)) < r2 then
          res := (j, k) :: !res (* set bit *)
      done
    done;
    (i, !res) in
  let mux (i, jk_l) =
    L.iter (fun (j, k) ->
        let idx = Grid.idx_of_ijk grid i j k in
        Bitv.set mask idx true
      ) jk_l in
  Parany.run nprocs ~demux ~work ~mux;
  mask

(* try placing [n] times the ligand in ROI *)
let place_ligand_in_ROI rng prot_bst prot roi centered_ligs n =
  (* We draw uniform random points in the ROI's englobing cube,
     until we have enough points inside the sphere *)
  let box =
    let r = ROI.get_out_radius roi in
    let c = ROI.get_center roi in
    let ray = V3.cube r in
    let low_corner  = V3.diff c ray in
    let high_corner = V3.add  c ray in
    Bbox.create_2v low_corner high_corner in
  let trials = ref 0 in
  let rec draw_one lig_i: (Rot.t * V3.t) =
    let centered_lig = centered_ligs.(lig_i) in
    (* translation draw *)
    let pos = Bbox.rand_point_inside rng box in
    incr trials;
    if !trials = 100_000 then
      (* don't try forever: some binding-sites are unreachable from the solvent
       * region and without the right ligand starting conformer *)
      let () =
        Log.fatal "Lds.place_ligands_in_ROI: ligand %s: start conformer \
                   cannot be placed in ROI after 100k trials"
          (Mol.get_name centered_lig) in
      exit 1
    else if not (ROI.is_inside roi pos) then
      draw_one lig_i (* pos KO *)
    else
      (* rotation draw *)
      let rot = Move.rand_rot_full rng in
      let lig' = Mol.rotate_then_translate_copy centered_lig rot pos in
      if not (Mol.heavy_atom_clash prot_bst prot lig') then
        (rot, pos)
      else
        (* rot & pos KO *)
        draw_one lig_i in
  let res = A.init n draw_one in
  Log.info "Lds.place_ligands_in_ROI: %d w/ %d trials" n !trials;
  res

(* random placing [n] times the ligand in ROI;
   clashes w/ protein are allowed
   (the optimization process has to find a better solution) *)
let drop_ligand_in_ROI rng roi n =
  (* We draw uniform random points in the ROI's englobing cube,
     until we have enough points inside the sphere *)
  let box =
    let r = ROI.get_out_radius roi in
    let c = ROI.get_center roi in
    let ray = V3.cube r in
    let low_corner  = V3.diff c ray in
    let high_corner = V3.add  c ray in
    Bbox.create_2v low_corner high_corner in
  let trials = ref 0 in
  let rec draw_one (): (Rot.t * V3.t) =
    (* translation draw *)
    let pos = Bbox.rand_point_inside rng box in
    incr trials;
    if not (ROI.is_inside roi pos) then
      draw_one () (* OOROI *)
    else
      (* rotation draw *)
      let rot = Move.rand_rot_full rng in
      (rot, pos) in
  let res = A.init n (fun _i -> draw_one ()) in
  Log.info "Lds.drop_ligand_in_ROI: %d w/ %d trials" n !trials;
  res

(* BST indexing of the receptor *)
let index_protein_atoms prot =
  let to_index =
    let coords = Mol.get_all_atom_coords prot in
    A.mapi Atom.create coords in
  Log.info "indexing %d protein atoms" (A.length to_index);
  Mol.BST.create 1 Bst.Bisec_tree.Two_bands to_index

module Defaults = struct
  let acceptance_ratio = 0.5
  let lig_init_pos = V3.origin (* MANDATORY because rotation applied first *)
  let lig_init_rot = Rot.id ()
  let temperature = Const.room_temp_K
end

let mk_temp_dir () =
  let exit_code, temp_dir' = BatUnix.run_and_read "mktemp -d /tmp/lds.XXXXXX" in
  assert(exit_code = Unix.WEXITED 0);
  let temp_dir = BatString.strip temp_dir' in
  Log.info "created out dir: %s" temp_dir;
  temp_dir

let create_dir dir =
  let cmd = "mkdir " ^ dir in
  let ret = Unix.system cmd in
  (if ret <> Unix.WEXITED 0 then
     let () = Log.fatal "command failed: %s; use -f to overwrite" cmd in
     exit 1
  );
  dir

let recreate_dir dir =
  let ret = Unix.system ("rm -rf " ^ dir) in
  assert(ret = Unix.WEXITED 0);
  create_dir dir

let pos_of_string xyz_s =
  Scanf.sscanf xyz_s "%f,%f,%f" V3.make

let write_SAS_to_bild_file fn lig_dotted_SAS =
  Log.info "creating: %s" fn;
  LO.with_out_file fn (fun out ->
      (* solvent-accessible --> blue *)
      fprintf out ".color blue\n";
      A.iter (
        A.iter (fun dot ->
            fprintf out ".sphere %s 0.05\n" (V3.short_str dot)
          )
      ) lig_dotted_SAS
    )

let pre_calculate_FF_component prot_grid bitmask prot ff single_atom_lig =
  (* or the translation below would not be valid *)
  assert(Mol.num_atoms single_atom_lig = 1);
  let dim_i, dim_j, dim_k = Grid.get_dims prot_grid in
  let res = G3D.create prot_grid in
  for i = 0 to dim_i - 1 do
    let x = prot_grid.xs.(i) in
    for j = 0 to dim_j - 1 do
      let y = prot_grid.ys.(j) in
      for k = 0 to dim_k - 1 do
        let idx = Grid.idx_of_ijk prot_grid i j k in
        if Bitv.get bitmask idx then
          begin
            let z = prot_grid.zs.(k) in
            (* single atom ligand translation:
               WARNING: ligand center is out of date after *)
            Mol.set_x_y_z single_atom_lig 0 x y z;
            (* bound grid E_max *)
            G3D.init prot_grid res i j k (min Params.max_E (ff prot single_atom_lig))
          end
      done
    done
  done;
  res

(* factorized version (much more efficient than previous one) *)
let pre_calculate_FF_components_grid prot_grid indexes prot_bst prot single_atom_ligs =
  (* one grid per atom type; created all at once *)
  let res =
    A.init (A.length single_atom_ligs) (fun _i -> G3D.create prot_grid) in
  L.iter (fun idx ->
      let i, j, k = Grid.ijk_of_idx prot_grid idx in
      let xyz =
        V3.make
          prot_grid.xs.(i)
          prot_grid.ys.(j)
          prot_grid.zs.(k) in
      let energies = Mol.ene_inter_UFF_shifted_grid prot_bst prot xyz single_atom_ligs in
      (* bound grid E_max *)
      A.iteri (fun l ene ->
          G3D.init_idx res.(l) idx (min Params.max_E ene)
        ) energies
    ) indexes;
  res

let all_ba1_files_exist lig_fn single_atoms =
  A.for_all (fun single_atom_mol ->
      let name = Mol.get_name single_atom_mol in
      let ba1_fn = sprintf "%s.%s.ba1" lig_fn name in
      if !Flags.no_compress then
        Sys.file_exists ba1_fn
      else
        Sys.file_exists (ba1_fn ^ ".zst")
    ) single_atoms

let pre_compute_FF_components
    nprocs apt lig_fn overwrite_grids grid bitmask prot prot_bst lig_FF_uniq_atoms =
  let res = A.make (A.length lig_FF_uniq_atoms) G3D.dummy in
  let to_do = ref (A.enumerate lig_FF_uniq_atoms) in
  let indexes = Bitv.to_list bitmask in
  let apt_i = ref 0 in
  let overwrite =
    overwrite_grids || not (all_ba1_files_exist lig_fn lig_FF_uniq_atoms) in
  if overwrite then Log.warn "--ogrids or some grid files missing";
  (* compute grids in parallel, but load them from memmapped bigarray files
   * to prevent passing large values (40MB files) between processes
   * (can crash Marshal) *)
  Parany.run nprocs
    ~demux:(fun () ->
        (* some user feedback/progress report in there might be nice *)
        let jobs, rest = A.takedrop !to_do apt.(!apt_i) in
        apt_i := (!apt_i + 1) mod nprocs;
        to_do := rest;
        if A.length jobs = 0 then
          raise Parany.End_of_input
        else
          jobs)
    ~work:(fun enumerated_atoms ->
        let single_atoms = A.map snd enumerated_atoms in
        let grids =
          if overwrite then
            pre_calculate_FF_components_grid grid indexes prot_bst prot single_atoms
          else
            [||] (* this array should be unused *) in
        A.mapi (fun grid_i (i, single_atom_mol) ->
            let name = Mol.get_name single_atom_mol in
            (* this file is not under the output dir, because we want
             * it to persist across runs *)
            let ba1_fn     = sprintf "%s.%s.ba1" lig_fn name in
            let ba1_zst_fn = ba1_fn ^ ".zst" in
            if !Flags.no_compress &&
               Sys.file_exists ba1_fn && (not overwrite_grids) then
              (i, ba1_fn)
            else if Sys.file_exists ba1_zst_fn && (not overwrite_grids) then
              (i, ba1_zst_fn)
            else
              let () = Log.info "creating %s" ba1_fn in
              let g3d = grids.(grid_i) in
              let ba1_dt, () = Utls.time_it (G3D.to_ba1_file ba1_fn grid) g3d in
              Log.info "wrote %s in %.3f (s)" ba1_fn ba1_dt;
              let res =
                if !Flags.no_compress then
                  ba1_fn
                else
                  Utls.zstd_compress_file ba1_fn in
              (if !Flags.verbose then
                 let mini, avg, maxi, sparse = G3D.min_avg_max_sparse g3d in
                 Log.info "(min,avg,max,non0) %s = (%g,%g,%g,%g)"
                   name mini avg maxi sparse
              );
              (i, res)
          ) enumerated_atoms
      )
    ~mux:(
      A.iter (fun (i, zstd_fn) ->
          (* DANGER: uncompressed file is transient;
           * for a given ligands file on one computer, ONLY SINGLE RUN IS OK *)
          let fn =
            if !Flags.no_compress then
              zstd_fn
            else
              Utls.zstd_uncompress_file zstd_fn in
          let ba1_dt, g3d = Utls.time_it (G3D.of_ba1_file grid) fn in
          Log.info "read %s in %.3f (s)" fn ba1_dt;
          if not !Flags.no_compress then
            Sys.remove fn
          ;
          res.(i) <- g3d)
    );
  res

(* accept 'k' or 'M' suffix for the number of steps CLI param *)
let steps_of_string s =
  let n = S.length s in
  let suffix = S.get s (n - 1) in
  let no_suffix = S.left s (n - 1) in
  match suffix with
  | '0' | '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9' ->
    (* no suffix *)
    int_of_string s
  | 'k' ->     1_000 * (int_of_string no_suffix)
  | 'M' -> 1_000_000 * (int_of_string no_suffix)
  | other ->
    let () = Log.fatal "Lds.steps_of_string: unsupported suffix: %c" other in
    exit 1

(* unless --no-interp, those FFs are just used to initialize the grid *)
type ff_style = Brute_global
              | Brute_local  (* local (shifted) version of Brute_global *)
              | Bst_local    (* fast but exact version of Brute_local *)

let ff_style_of_string = function
  | "BrG" -> Brute_global
  | "BrL" -> Brute_local
  | "Bst" -> Bst_local
  | other ->
    let () = Log.fatal "Lds.ff_style_of_string: unsupported: %s" other in
    exit 1

(* adaptive rot. trans. sampling scheme *)
let update_rot_trans_params frame mol_name
    accepts_rejects target_low target_high max_rot max_trans =
  let accept_ratio = SW.get_ratio accepts_rejects in
  (* accept_ratio staying at 0.0 means there is a problem w/ the simulation
   * like a new bug just introduced... *)
  if !Flags.verbose || accept_ratio = 0.0 || accept_ratio = 1.0 then
    Log.warn "f%#d mol:%s AR: %.2f" frame mol_name accept_ratio;
  if accept_ratio <= target_low then
    (* too many collisions: maneuver more delicately/decelerate *)
    (max_trans := 0.95 *. !max_trans;
     max_rot   := 0.95 *. !max_rot)
  else if accept_ratio >= target_high then
    (* accelerate *)
    (max_trans := 1.05 *. !max_trans;
     max_rot   := min Math.pi (1.05 *. !max_rot))

(* adaptive RotBond parameter *)
let update_rbond_rot_param target_low target_high lig i =
  let accept_ratio = SW.get_ratio (Mol.get_rbond_sw lig i) in
  (* accept_ratio staying at 0.0 means there is a problem w/ the simulation
   * like a new bug just introduced... *)
  (* (if accept_ratio = 0.0 || accept_ratio = 1.0 then *)
  (*    Log.warn "Rbond %d AR: %g" i accept_ratio *)
  (* ); *)
  if accept_ratio <= target_low then
    (* too many collisions: maneuver more delicately/decelerate *)
    Mol.set_rbond_dr lig i (0.95 *. Mol.get_rbond_dr lig i)
  else if accept_ratio >= target_high then
    (* accelerate *)
    Mol.set_rbond_dr lig i (min Math.pi (1.05 *. Mol.get_rbond_dr lig i))

let update_rbond_rot_params low high lig =
  let n = Mol.num_rbonds lig in
  for i = 0 to n - 1 do
    update_rbond_rot_param low high lig i
  done

(* to examine if ligand start poses are OK (no clash, inside ROI) *)
let dump_start_ligands out_dir base_lig_fn lig_i centered_ligs poses =
  A.iteri (fun start_i (rot, pos) ->
      let lig = Mol.rotate_then_translate_copy centered_ligs.(start_i) rot pos in
      let fn = sprintf "%s/%s_l%d_s%d.bild"
                 out_dir base_lig_fn lig_i start_i in
      Mol.to_bild ~style:Ptable.White_H fn lig
    ) poses

let reset_run_params rot0 pos0 centered_lig
      max_rot max_trans rot pos best_rot best_pos prev_E best_E prev_lig best_lig
      accepts_rejects =
  max_rot   := Params.max_rot;
  max_trans := Params.max_trans;
  rot := rot0;
  pos := pos0;
  best_rot := (Rot.id ());
  best_pos := V3.origin;
  (* interaction Energy is negative *)
  prev_E := infinity;
  best_E := infinity;
  (* the following _also_ resets the rbonds_params;
     like Mol.reset_rbonds_params does *)
  prev_lig := Mol.rotate_then_translate_copy centered_lig !rot !pos;
  best_lig := !prev_lig;
  SW.reset accepts_rejects

(* hysteresis parameters *)
let target_low  = Defaults.acceptance_ratio -. 0.05
let target_high = Defaults.acceptance_ratio +. 0.05

(* we alternate rot/trans with conformational moves so that
 * rot/trans parameters are decoupled from conformational ones
 * (sampling should be more efficient)
 * #frame mod 2 = 0: whole ligand rigid body move
 * #frame mod 2 = 1: ligand conformational step *)
let is_RIGID_BODY_MOVE = 0
let is_CONFORMER_MOVE = 1

let lig_E_intra_QM xyz_fn =
  let cmd =
    sprintf "./src/QM_energy.py %s /dev/stdout 2>/dev/null | \
             grep -v torchani | tail -1 | cut -f2" xyz_fn in
  Log.info "running: %s" cmd;
  let ret, res = BatUnix.run_and_read cmd in
  if ret <> Unix.WEXITED 0 then
    let () = Log.fatal "Lds.lig_E_intra_QM: command failed: %s" cmd in
    exit 1
  else if res = "nan\n" then
    nan
  else
    try Scanf.sscanf res "%f" (fun x -> x)
    with exn ->
      (Log.fatal "lig_E_intra_QM: could not parse: %s" res;
       raise exn)

(* ligand's (MM level) FF's lig_E_intra is ignored because vdW potential should prevent clashes;
   also, it is super imprecise; cf. many authors *)
let rescore mol_name ene_rescore best_lig out =
  let num_rbonds = Mol.num_rbonds best_lig in
  let num_HAs = Mol.num_heavy_atoms best_lig in
  (* e_inter_exact (not interpolated) might deviate a little from best_score *)
  let dt, e_inter_exact = Utls.time_it ene_rescore best_lig in
  Log.info "Lds.rescore.ene_rescore %s: %fs" mol_name dt;
  fprintf out "%s\t%g\t%d\t%d\n%!" mol_name e_inter_exact num_rbonds num_HAs

type force_field = NULL (* only for lig_E_intra: always 0.0 *)
                 | UFF_NB (* partial UFF implementation;
                             only non-bonded interactions *)
                 | RDKIT_UFF
                 | RDKIT_MMFF
                 | NO_VDW_CLASH
                 | TORCHANI_QM

let requires_rdkit = function
  | NULL | UFF_NB | NO_VDW_CLASH | TORCHANI_QM -> false
  | RDKIT_MMFF | RDKIT_UFF -> true

(* how to get ene_intra for all simulation modes *)
let ene_intra_builder
    intra_FF tweak_rbonds num_rbonds centered_lig maybe_rdkit_conf =
  match intra_FF with
  | NULL -> (fun _lig -> 0.0)
  | UFF_NB ->
    if tweak_rbonds && num_rbonds > 0 then
      Mol.ene_intra_UFFNB_brute
    else (* rigid ligand *)
      let const_ene_intra =
        Mol.ene_intra_UFFNB_brute centered_lig in
      (fun _lig -> const_ene_intra)
  | TORCHANI_QM ->
    let species, init_E = Mol.ene_QM_prepare centered_lig in
    if tweak_rbonds && num_rbonds > 0 then
      (* !!! FBR: WARNING: this assumes that the starting conformer is the lowest E_QM conf. !!! *)
      (fun lig -> (Mol.ene_QM_score species lig) -. init_E)
    else (* rigid ligand *)
      (fun _lig -> init_E)
  | RDKIT_UFF
  | RDKIT_MMFF ->
    let intra_ene = match intra_FF with
      | RDKIT_UFF  -> Mol.ene_intra_UFF_rdkit_exn
      | RDKIT_MMFF -> Mol.ene_intra_MMFF_rdkit_exn
      | _ -> assert(false) in
    let rdkit_conf = Option.get maybe_rdkit_conf in
    if tweak_rbonds && num_rbonds > 0 then
      (fun lig -> intra_ene lig rdkit_conf)
    else (* rigid ligand *)
      let const_ene_intra = intra_ene centered_lig rdkit_conf in
      (fun _lig -> const_ene_intra)
  | NO_VDW_CLASH ->
    let topo_dists_ge3 = Mol.get_topo_dists_ge3 centered_lig in
    if tweak_rbonds && num_rbonds > 0 then
      Mol.ene_intra_no_vdW_clash topo_dists_ge3
    else (* rigid ligand *)
      let const_ene_intra =
        Mol.ene_intra_no_vdW_clash topo_dists_ge3 centered_lig in
      (fun _lig -> const_ene_intra)

let simulate_lig tweak_rbonds intra_FF sequential_run enforce_ROI
    rng nsteps traj_N ene_N prot_fn roi ene_inter ene_rescore
    beta_at_K scores_out maybe_rescore_out out_dir
    base_lig_fn lig_i rot0 pos0 start_i centered_lig maybe_rdkit_conf =
  let ooroi = ref 0 in
  let later_funs = ref [] in
  (* replacement for at_exit for parallel programs *)
  let later f =
    later_funs := f :: !later_funs in
  let later_is_now () =
    L.iter (fun f -> f ()) !later_funs in
  (* initialize run parameters *)
  let max_rot   = ref Params.max_rot   in
  let max_trans = ref Params.max_trans in
  let rot = ref rot0 in
  let pos = ref pos0 in
  let best_rot = ref (Rot.id ()) in
  let best_pos = ref V3.origin in
  let start_conf = Mol.rotate_then_translate_copy centered_lig rot0 pos0 in
  let prev_lig = ref start_conf in
  let best_lig = ref !prev_lig in
  (* in case of flexible ligand, conformational changes accumulate in conf *)
  let conf = ref (Mol.copy centered_lig) in
  let accepts_rejects = SW.create Params.block_size in
  let mol_name =
    sprintf "%s_l%d_s%d"
      (Mol.get_name centered_lig) lig_i start_i in
  let start_lig_fn =
    sprintf "%s/%s_l%d_s%d_start.bild"
      out_dir base_lig_fn lig_i start_i in
  let best_lig_bild_fn =
    sprintf "%s/%s_l%d_s%d_best.bild"
      out_dir base_lig_fn lig_i start_i in
  let stop_lig_fn =
    sprintf "%s/%s_l%d_s%d_stop.bild"
      out_dir base_lig_fn lig_i start_i in
  let out_radius = ROI.get_out_radius roi in
  let roi_center = ROI.get_center roi in
  let num_rbonds = Mol.num_rbonds centered_lig in
  let rbonds_block_size = Params.block_size * num_rbonds in
  let rotate_bond =
    if tweak_rbonds && num_rbonds > 0 then
      let rbf = !Params.rbond_flip_block in
      (fun i rng lig0 ->
         let lig = Mol.copy lig0 in
         (* do not modify significantly the start conformer *)
         let just_rotated =
           if i > 0 && i mod rbf = 0 then
             (* let () = Log.info "rbond flip %d" i in *)
             Mol.flip_rbond rng lig
           else
             (* let () = Log.info "rbond tweak %d" i in *)
             Mol.tweak_rbond rng lig in
         Mol.check_elongation_exn lig Const.charged_cutoff;
         (just_rotated, lig)
      )
    else (* do nothing *)
      (fun _i _rng _lig -> (-1, centered_lig)) in
  let ene_intra =
    ene_intra_builder
      intra_FF tweak_rbonds num_rbonds centered_lig maybe_rdkit_conf in
  let prev_E_intra = ref (ene_intra !prev_lig) in
  let prev_E_inter = ref (ene_inter !prev_lig) in
  let prev_E = ref (!prev_E_inter +. !prev_E_intra) in
  let best_E = ref !prev_E in
  Log.info "f%#d l%d s%d minE: %.2f"
    0 lig_i start_i !prev_E;
  Mol.to_bild ~verbose:true ~style:Ptable.Orange_H ~header:(sprintf "E: %g" !prev_E)
    start_lig_fn !prev_lig;
  let ene_dump =
    let ene_out =
      let ene_fn = sprintf "%s/%s_l%d_s%d_E.txt"
          out_dir base_lig_fn lig_i start_i in
      open_out ene_fn in
    (* file format header comment *)
    fprintf ene_out "#frame dist RMSD E_tot E_inter E_intra AR DR DT ligARs ligDRs\n";
    later (fun () -> close_out ene_out);
    if ene_N > 0 then
      (* detailed prod. run: user want details *)
      (fun step dist e_inter e_intra ar dr dt lig ->
         if step mod ene_N = 0 then
           let rmsd = Mol.heavy_atoms_rmsd start_conf lig in
           let lig_accept_rates = Mol.sprintf_rbonds_sw lig in
           let lig_rot_params   = Mol.sprintf_rbonds_dr lig in
           fprintf ene_out "%d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %f %s %s\n"
             step dist rmsd (e_inter +. e_intra) e_inter e_intra
             ar dr dt lig_accept_rates lig_rot_params
      )
    else
      (* quick prod. run: only log lower E *)
      (fun step dist e_inter e_intra ar dr dt lig ->
         if (e_inter +. e_intra) < !best_E then
           let rmsd = Mol.heavy_atoms_rmsd start_conf lig in
           let lig_accept_rates = Mol.sprintf_rbonds_sw lig in
           let lig_rot_params   = Mol.sprintf_rbonds_dr lig in
           fprintf ene_out "%d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %f %s %s\n"
             step dist rmsd (e_inter +. e_intra) e_inter e_intra
             ar dr dt lig_accept_rates lig_rot_params
      ) in
  let traj_dump =
    if traj_N = 0 then
      (fun _i _mol -> ())
    else
      let traj_out =
        let traj_fn =
          sprintf "%s/%s_l%d_s%d.xyz"
            out_dir base_lig_fn lig_i start_i in
        open_out traj_fn in
      later (fun () -> close_out traj_out);
      if traj_N = 1 then
        (fun i mol -> Mol.xyz_dump traj_out i "" mol)
      else
        (fun i mol ->
           if i mod traj_N = 0 then
             Mol.xyz_dump traj_out i "" mol) in
  (* keep best pose even in case of crash *)
  later (fun () ->
      let best_E_str = sprintf "E: %g" !best_E in
      let prev_E_str = sprintf "E: %g" !prev_E in
      Log.info "l%d s%d %s" lig_i start_i best_E_str;
      fprintf scores_out "%s\t%f\n%!" mol_name !best_E;
      BatOption.may
        (rescore mol_name ene_rescore !best_lig)
        maybe_rescore_out;
      Log.info "l%d s%d best_pos: %s"
        lig_i start_i (V3.short_str !best_pos);
      let a, b, g = Rot.decompose !best_rot in
      Log.info "l%d s%d best_rot: (%g, %g, %g)"
        lig_i start_i a b g;
      Mol.to_bild ~verbose:true ~style:Ptable.Cyan_H
        ~header:best_E_str best_lig_bild_fn !best_lig;
      Mol.to_bild ~verbose:true ~style:Ptable.Pink_H
        ~header:prev_E_str stop_lig_fn !prev_lig;
      Log.(info "########## %s ##########" (color green "CHIMERA COMMAND"));
      (* just show in chimera receptor, lig_start and lig_best;
         lig_stop is saved in the output dir though *)
      eprintf "chimera %s %s %s\n%!" prot_fn start_lig_fn best_lig_bild_fn
    );
  (try
     let rigid_step = ref 0 in
     let conf_step = ref 0 in
     for frame = 0 to nsteps - 1 do
       let xRIGID_MOVE = (frame mod 2 = is_RIGID_BODY_MOVE) in
       let xCONF_MOVE = not xRIGID_MOVE in
       (* tweak/flip one rotatable bond *)
       let just_rotated, conf' =
         if xCONF_MOVE then
           rotate_bond !conf_step rng !conf
         else (-1, !conf) in
       (* rotate then translate ligand *)
       let rot', pos' =
         if xRIGID_MOVE then
           (Move.rand_rot   rng !max_rot   !rot,
            Move.rand_trans rng !max_trans !pos)
         else (!rot, !pos) in
       let lig' = Mol.center_rotate_translate_copy conf' rot' pos' in
       (* it might be nice to output only ligand's frame to traj. file,
        * if some CLI option is present *)
       (* FF eval *)
       (if xCONF_MOVE then
          (* only update E_intra in case of conformational move *)
          prev_E_intra := ene_intra lig'
       );
       prev_E_inter := ene_inter lig';
       let curr_E = !prev_E_inter +. !prev_E_intra in
       (* if hard ROI: outside of ROI -> EMERGENCY STOP *)
       let dist_roi = V3.dist roi_center (Mol.get_center lig') in
       if enforce_ROI then
         if dist_roi > out_radius then
           begin
             (* let rem_steps = nsteps - (frame + 1) in *)
             reset_run_params rot0 pos0 centered_lig max_rot max_trans rot pos
               best_rot best_pos prev_E best_E prev_lig best_lig accepts_rejects;
             incr ooroi
             (* Log.warn "%s: OOROI; %d steps left" mol_name rem_steps *)
           end
       else
         begin
           (* E=0 -> EMERGENCY STOP (there is NO reboxing) *)
           if !prev_E_inter = 0.0 then
             (* ligand went too far away from protein *)
             let rem_steps = nsteps - (frame + 1) in
             reset_run_params rot0 pos0 centered_lig
               max_rot max_trans rot pos best_rot best_pos prev_E best_E
               prev_lig best_lig accepts_rejects;
             Log.error "%s: E=0; %d steps left"
               mol_name rem_steps
           else
             begin
               (* Metropolis criterion *)
               (if (curr_E <= !prev_E) || (* always allow energy to decrease *)
                   (* accept increase depending on dice roll and deltaE *)
                   ((Random.State.float rng 1.0) < exp (-.(curr_E -. !prev_E) *.
                                                          beta_at_K)) then
                  begin
                    (if xRIGID_MOVE then
                       SW.process accepts_rejects true
                     else if just_rotated > -1 then
                       SW.process (Mol.get_rbond_sw conf' just_rotated) true
                    );
                    rot := rot';      (* update pose *)
                    pos := pos';
                    prev_E := curr_E; (* update E *)
                    prev_lig := lig'; (* update lig *)
                    conf := conf';    (* update conformer *)
                  end
                else
                  begin
                    if xRIGID_MOVE then
                      SW.process accepts_rejects false
                    else if just_rotated > -1 then
                      SW.process (Mol.get_rbond_sw conf' just_rotated) false
                  end
               );
               ene_dump frame dist_roi !prev_E_inter !prev_E_intra
                 (SW.get_ratio accepts_rejects) !max_rot !max_trans lig';
               traj_dump frame !prev_lig;
               (* monitor E *)
               (if curr_E < !best_E then
                  begin
                    best_E := curr_E;
                    best_rot := rot';
                    best_pos := pos';
                    best_lig := lig';
                    if sequential_run then
                      Log.info "f%#d l%d s%d minE: %.2f"
                        frame lig_i start_i curr_E
                  end
               );
               (* update dr dt *)
               (if xRIGID_MOVE && !rigid_step > 0 && !rigid_step mod Params.block_size = 0 then
                  update_rot_trans_params frame mol_name accepts_rejects
                    target_low target_high max_rot max_trans
               );
               (* update dr for each RotBond *)
               (if num_rbonds > 0 && tweak_rbonds && xCONF_MOVE &&
                   !conf_step > 0 && !conf_step mod rbonds_block_size = 0 then
                  begin
                    update_rbond_rot_params target_low target_high conf';
                    if Mol.strange_rbond_acceptance_rate conf' then
                      begin (* inspect rbonds rot. params *)
                        Log.warn "f%#d mol:%s rbonds_AR: %s" frame mol_name
                          (Mol.sprintf_rbonds_sw conf');
                        Log.warn "f%#d mol:%s rbonds_dr: %s" frame mol_name
                          (Mol.sprintf_rbonds_dr conf')
                      end
                  end
               );
             end
         end;
       if xRIGID_MOVE then
         incr rigid_step
       else
         incr conf_step
     done
   with Mol.Ene_intra -> Log.error "Mol.Ene_intra %s" mol_name
      | Mol.Too_long -> Log.error "Mol.Too_long %s" mol_name
  );
  Log.info "%s: OOROI: %d" mol_name !ooroi;
  later_is_now ()

(* easy selection of optim. algo. from the CLI *)
let glob_opt_of_int = function
  | 0 -> Nlopt.direct        (* OK; was seen to work for 3RFM:CFF starting from binding-site center *)
  | 1 -> Nlopt.direct_l_rand (* OK *)
  | 2 -> Nlopt.direct_noscal
  | 3 -> Nlopt.direct_l_noscal
  | 4 -> Nlopt.direct_l_rand_noscal
  | 5 -> Nlopt.orig_direct_l (* OK *)
  | 6 -> Nlopt.crs2_lm
  | other -> (Log.fatal "Lds.glob_opt_of_int: unsupported: %d" other;
              exit 1)

let nlopt_minimize tweak_rbonds max_evals roi lig config0 score =
  let ndims = A.length config0 in
  (* hard/stupid stop conditions *)
  let min_dock_score = -1000.0 in (* minimum docking score? *)
  let global = Nlopt.(create (glob_opt_of_int !Flags.glob_opt) ndims) in
  (* Nlopt.stogo -> NLOPT_INVALID_ARGS
   * Nlopt.lbfgs -> NLOPT_FAILURE after 14 steps
   *                (probably because numerical gradient becomes misleading *)
  Nlopt.set_min_objective global score;
  (* parameter bounds *)
  let min_bounds = Optim.get_min_bounds roi tweak_rbonds lig in
  let max_bounds = Optim.get_max_bounds roi tweak_rbonds lig in
  Log.info "min_bounds: %s" (Utls.string_of_array min_bounds " " (sprintf "%.2f"));
  Log.info "cur_values: %s" (Utls.string_of_array config0    " " (sprintf "%.2f"));
  Log.info "max_bounds: %s" (Utls.string_of_array max_bounds " " (sprintf "%.2f"));
  Nlopt.set_lower_bounds global min_bounds;
  Nlopt.set_upper_bounds global max_bounds;
  (* hard/stupid stop conditions *)
  Nlopt.set_stopval global min_dock_score;
  Nlopt.set_maxeval global max_evals;
  let stop_cond, params, min_score = Nlopt.optimize global config0 in
  Log.info "NLopt optimize global: %s" (Nlopt.string_of_result stop_cond);
  (params, min_score)

(* current brute method: score all possible rotations and translations
 * of the ligand; return best one *)
let exhaustive_rigid_ligand_docking
    topk roi prot_grid prot_vdW_mask trans_step rotations centered_lig score =
  (* FBR: add a global flag to allow disabling each optim *)
  let lig_center_clash =
    if Mol.is_ligand_center_vdW_occuppied centered_lig then
      G3D.vdW_clash_AND prot_grid prot_vdW_mask
      (* (fun pos -> *)
      (*    let res = G3D.vdW_clash_AND prot_grid prot_vdW_mask pos in *)
      (*    if res then incr avoid_E_inter; *)
      (*    res) *)
    else
      (fun _pos -> false) in
  let x_min, x_max,
      y_min, y_max,
      z_min, z_max = ROI.get_bounds roi in
  let top_k =
    if topk > 0 then
      TopK.create topk
    else
      TopK.create 1 (* but won't be used *) in
  let top_keep =
    if topk <= 0 then
      ignore
    else
      (fun score -> TopK.add top_k score ()) in
  let g =
    Grid.from_box trans_step
      (Bbox.create_6f
         x_min y_min z_min
         x_max y_max z_max) in
  (* relative (ROI grid) to absolute (binding-site) coordinates *)
  let xs = A.map ((+.) x_min) Grid.(g.xs) in
  let ys = A.map ((+.) y_min) Grid.(g.ys) in
  let zs = A.map ((+.) z_min) Grid.(g.zs) in
  let num_voxels = Grid.num_voxels g in
  Log.info "num_voxels: %#d" num_voxels;
  let num_rotations = A.length rotations in
  let tot = sprintf "%#d" (num_voxels * num_rotations) in
  let best_score = ref infinity in
  let best_pos = ref V3.origin in
  let best_rot_i = ref 0 in
  let rotated_copies =
    A.map (Mol.centered_rotate_copy centered_lig) rotations in
  A.iteri (fun k z ->
      let z_dim = k * g.xy_dim in
      A.iteri (fun j y ->
          let jz_dim = j * g.x_dim + z_dim in
          A.iteri (fun i x ->
              let pos = V3.make x y z in
              (* we are sweaping the whole cube englobing the ROI sphere,
                 but only inside the sphere can E_inter be interpolated *)
              if ROI.is_inside roi pos && not (lig_center_clash pos) then
                (* sweap all rotations *)
                A.iteri (fun rot_i rot_lig ->
                    let lig' = Mol.translate_copy_to rot_lig pos in
                    if not (Mol.protein_ligand_clash prot_grid prot_vdW_mask lig') then
                      let curr_score = score lig' in
                      (* good docking scores are negative numbers *)
                      top_keep (-.curr_score);
                      if curr_score < !best_score then
                        let frame = rot_i + num_rotations * (i + jz_dim) in
                        best_score := curr_score;
                        best_pos := pos;
                        best_rot_i := rot_i;
                        Log.info "f%#d/%s %.2f | %.2f %.2f %.2f %d"
                          frame tot curr_score x y z rot_i
                  ) rotated_copies
            ) xs
        ) ys
    ) zs;
  let top_scores =
    let xs = TopK.high_scores_first top_k in
    (* put back proper score sign *)
    A.of_list (L.map (fun (s, _) -> -.s) xs) in
  (top_scores, !best_score, !best_pos, rotations.(!best_rot_i))

(* numerical approximation of the local derivative for each dimension
 * result is written directly to [res] *)
let approx_local_gradient ene_intra ene_inter centered_lig pos_min_steps params res =
  (* local copy *)
  let copy = A.copy params in
  A.iteri (fun i eps_i ->
      (* before *)
      copy.(i) <- params.(i) -. eps_i;
      let before =
        let lig = Optim.apply_config centered_lig copy in
        ene_intra lig +. ene_inter lig in
      (* after *)
      copy.(i) <- params.(i) +. eps_i;
      let after =
        let lig = Optim.apply_config centered_lig copy in
        ene_intra lig +. ene_inter lig in
      (* reset *)
      copy.(i) <- params.(i);
      (* discrete approximation *)
      res.(i) <- (after -. before) /. (2. *. eps_i)
    ) pos_min_steps

(* docking w/ NLopt *)
let minimize_ligand tweak_rbonds intra_FF max_steps traj_N prot_fn roi
    ene_inter ene_rescore
    scores_out maybe_rescore_out
    out_dir pos_rot_out base_lig_fn lig_i rot0 pos0 start_i
    centered_lig maybe_rdkit_conf =
  let later_funs = ref [] in
  (* replacement for at_exit for parallel programs *)
  let later f =
    later_funs := f :: !later_funs in
  let later_is_now () =
    L.iter (fun f -> f ()) !later_funs in
  (* initialize run parameters *)
  let start_conf = Mol.rotate_then_translate_copy centered_lig rot0 pos0 in
  let prev_lig = ref start_conf in
  let best_lig = ref !prev_lig in
  (* in case of flexible ligand, conformational changes accumulate in conf *)
  let orig_name = Mol.get_name centered_lig in
  let mol_name =
    sprintf "%s_l%d_s%d" orig_name lig_i start_i in
  let start_lig_fn =
    sprintf "%s/%s_l%d_s%d_start.bild"
      out_dir base_lig_fn lig_i start_i in
  let best_lig_bild_fn =
    sprintf "%s/%s_l%d_s%d_best.bild"
      out_dir base_lig_fn lig_i start_i in
  let best_lig_mol2_fn =
    sprintf "%s/%s_l%d_s%d_best.mol2"
      out_dir base_lig_fn lig_i start_i in
  let stop_lig_fn =
    sprintf "%s/%s_l%d_s%d_stop.bild"
      out_dir base_lig_fn lig_i start_i in
  let init_mol2_fn = sprintf "%s/ligs/%d.mol2" out_dir lig_i in
  Log.info "reading %s" init_mol2_fn;
  let init_mol2 = Mol2.read_one_from_file init_mol2_fn in
  let roi_center = ROI.get_center roi in
  let num_rbonds = Mol.num_rbonds centered_lig in
  let ene_intra =
    ene_intra_builder
      intra_FF tweak_rbonds num_rbonds centered_lig maybe_rdkit_conf in
  let prev_E_intra = ref (ene_intra !prev_lig) in
  let prev_E_inter = ref (ene_inter !prev_lig) in
  (* to check if energy decreased during docking *)
  let init_E = !prev_E_inter +. !prev_E_intra in
  let prev_E = ref init_E in
  let best_E = ref !prev_E in
  Log.info "f%#d l%d s%d minE: %.2f" 0 lig_i start_i !prev_E;
  Mol.to_bild
    ~verbose:true ~style:Ptable.Orange_H ~header:(sprintf "E: %g" !prev_E)
    start_lig_fn !prev_lig;
  let ene_dump =
    let ene_out =
      let ene_fn =
        sprintf "%s/%s_l%d_s%d_E.txt"
          out_dir base_lig_fn lig_i start_i in
      open_out ene_fn in
    (* file format header comment *)
    fprintf ene_out "#frame dist RMSD E_tot E_inter E_intra\n";
    later (fun () -> close_out ene_out);
    (* quick prod. run: only log lower E *)
    (fun step dist e_inter e_intra lig ->
       let e_tot = e_inter +. e_intra in
       if e_tot < !best_E then
         let rmsd = Mol.heavy_atoms_rmsd start_conf lig in
         fprintf ene_out "%d %.3f %.3f %.3f %.3f %.3f\n"
           step dist rmsd e_tot e_inter e_intra
    ) in
  let traj_dump =
    if traj_N = 0 then
      (fun _i _mol -> ()) (* OFF *)
    else
      let traj_out =
        let traj_fn =
          sprintf "%s/%s_l%d_s%d.xyz"
            out_dir base_lig_fn lig_i start_i in
        open_out traj_fn in
      later (fun () -> close_out traj_out);
      if traj_N = 1 then
        (fun i mol -> Mol.xyz_dump traj_out i "" mol)
      else
        (fun i mol ->
           if i mod traj_N = 0 then
             Mol.xyz_dump traj_out i "" mol) in
  (* keep best pose even in case of crash *)
  later (fun () ->
      let best_E_str = sprintf "E: %g" !best_E in
      let prev_E_str = sprintf "E: %g" !prev_E in
      Log.info "l%d s%d %s" lig_i start_i best_E_str;
      Mol.to_bild ~verbose:true ~style:Ptable.Cyan_H
        ~header:best_E_str best_lig_bild_fn !best_lig;
      Mol.to_bild ~verbose:true ~style:Ptable.Pink_H
        ~header:prev_E_str stop_lig_fn !prev_lig;
      Log.(info "########## %s ##########" (color green "CHIMERA COMMAND"));
      (* just show in chimera receptor, lig_start and lig_best;
         lig_stop is saved in the output dir though *)
      eprintf "chimera %s %s %s\n%!" prot_fn start_lig_fn best_lig_bild_fn
    );
  (try
     let frame = ref 0 in
     let start_params =
       let x0, y0, z0 = V3.to_triplet pos0 in
       let alpha, beta, gama = Rot.decompose rot0 in
       Optim.create x0 y0 z0 alpha beta gama
         (* rigid-ligand docking handled here *)
         (if tweak_rbonds then A.make num_rbonds 0.0 else [||]) in
     let pos_min_steps = Optim.step_sizes start_params 0.1 (Math.pi /. 180.0) in
     let best_params, best_score =
       nlopt_minimize tweak_rbonds max_steps roi centered_lig start_params
         (fun config grad ->
            let lig' = Optim.apply_config centered_lig config in
            let dist_roi = V3.dist roi_center (Mol.get_center lig') in
            (* FF eval *)
            prev_E_intra := ene_intra lig';
            prev_E_inter := ene_inter lig';
            let curr_E = !prev_E_inter +. !prev_E_intra in
            (if curr_E <= !prev_E then (* always log energy decrease *)
               (prev_E := curr_E; (* update E *)
                prev_lig := lig') (* update lig *)
            );
            ene_dump !frame dist_roi !prev_E_inter !prev_E_intra lig';
            traj_dump !frame !prev_lig;
            (* monitor E *)
            (if curr_E < !best_E then
               (best_E := curr_E;
                best_lig := lig';
                Log.info "f%#d l%d s%d minE: %.2f | %s"
                  !frame lig_i start_i curr_E)
                 (Optim.to_string config)
            );
            incr frame;
            (* compute local gradient *)
            BatOption.may
              (approx_local_gradient ene_intra ene_inter centered_lig pos_min_steps config)
              grad;
            (* local value *)
            curr_E
         ) in
     (if best_score <> !best_E then
        Log.warn "%s: best scores diff: %g %g %g"
          mol_name best_score !best_E (best_score -. !best_E)
     );
     (* best docking score: "nan" <=> molecule could not be minimized *)
     let final_score = if best_score < init_E then best_score else nan in
     fprintf scores_out "%s\t%f\n%!" mol_name final_score;
     let config_str = Optim.to_string best_params in
     fprintf pos_rot_out "%s %s\n%!" orig_name config_str;
     Log.info "%s minE: %f f%#d l%d s%d | %s"
       mol_name final_score !frame lig_i start_i config_str;
     (* rescoring *)
     BatOption.may
       (rescore mol_name ene_rescore !best_lig)
       maybe_rescore_out;
     (* save as mol2 *)
     let mol2 = Mol.update_mol2 init_mol2 !best_lig in
     Log.info "creating %s" best_lig_mol2_fn;
     Mol2.write_one_to_file best_lig_mol2_fn mol2
   with Mol.Ene_intra -> Log.error "Mol.Ene_intra %s" mol_name
      | Mol.Too_long -> Log.error "Mol.Too_long %s" mol_name
  );
  later_is_now ()

let exhaustive_docking topk intra_FF trans_step rotations prot_fn
    roi prot_grid prot_vdW_mask ene_inter
    out_dir scores_out pos_rot_out base_lig_fn lig_i
    centered_lig maybe_rdkit_conf =
  let name = Mol.get_name centered_lig in
  let mol2_out_fn = sprintf "%s/%s_l%d_best.mol2" out_dir base_lig_fn lig_i in
  let top_scores_out_fn = sprintf "%s/%s_l%d_best.scores" out_dir base_lig_fn lig_i in
  let init_mol2_fn = sprintf "%s/ligs/%d.mol2" out_dir lig_i in
  Log.info "reading %s" init_mol2_fn;
  let mol2 = Mol2.read_one_from_file init_mol2_fn in
  (* ligand is considered rigid -> constant E_intra *)
  let const_ene_intra = match intra_FF with
    | NULL -> 0.0
    | TORCHANI_QM ->
      let _species, init_E = Mol.ene_QM_prepare centered_lig in
      init_E
    | NO_VDW_CLASH ->
      Mol.(ene_intra_no_vdW_clash (get_topo_dists_ge3 centered_lig)
             centered_lig)
    | UFF_NB ->
      Mol.ene_intra_UFFNB_brute centered_lig
    | RDKIT_UFF ->
      Mol.ene_intra_UFF_rdkit_exn  centered_lig (Option.get maybe_rdkit_conf)
    | RDKIT_MMFF ->
      Mol.ene_intra_MMFF_rdkit_exn centered_lig (Option.get maybe_rdkit_conf) in
  let score lig =
    const_ene_intra +. (ene_inter lig) in
  let top_scores, min_E, best_pos, best_rot =
    exhaustive_rigid_ligand_docking
      topk roi prot_grid prot_vdW_mask
      trans_step rotations centered_lig score in
  fprintf scores_out "%s\t%f\n%!" name min_E;
  (* create config *)
  let x, y, z = V3.to_triplet best_pos in
  let a, b, g = Rot.decompose best_rot in
  (* log it *)
  let num_rbonds = Mol.num_rbonds centered_lig in
  let rbonds = A.make num_rbonds 0.0 in
  let config = Optim.create x y z a b g rbonds in
  let config_str = Optim.to_string config in
  Log.info "%s minE: %.2f | %s" name min_E config_str;
  fprintf pos_rot_out "%s %s\n%!" name config_str;
  Log.(info "########## %s ##########" (color green "CHIMERA COMMAND"));
  (* just show in chimera receptor and lig_best *)
  eprintf "chimera %s %s\n%!" prot_fn mol2_out_fn;
  (* apply transform *)
  let mol = Optim.apply_config centered_lig config in
  (* save as mol2 *)
  let mol2' = Mol.update_mol2 mol2 mol in
  Mol2.write_one_to_file mol2_out_fn mol2';
  (* dump top scores *)
  LO.with_out_file top_scores_out_fn (fun out ->
      A.iter (fprintf out "%g\n") top_scores
    )

let simulate_isolated_lig ff sequential_run rng nsteps traj_N ene_N
    beta_at_K out_dir base_lig_fn lig_i start_i centered_lig maybe_rdkit_conf =
  let num_rbonds = Mol.num_rbonds centered_lig in
  (if num_rbonds = 0 then
     (Log.fatal "Lds.simulate_isolated_lig: ligand %s is rigid"
        (Mol.get_name centered_lig);
      exit 1)
  );
  let rbonds_block_size = Params.block_size * num_rbonds in
  let later_funs = ref [] in
  (* replacement for at_exit for parallel programs *)
  let later f =
    later_funs := f :: !later_funs in
  let later_is_now () =
    L.iter (fun f -> f ()) !later_funs in
  (* initialize run parameters *)
  let prev_E = ref infinity in (* better E is smaller *)
  let best_E = ref infinity in
  let prev_lig = ref (Mol.copy centered_lig) in
  let best_lig = ref !prev_lig in
  let init_mol2_fn = sprintf "%s/ligs/%d.mol2" out_dir lig_i in
  Log.info "reading %s" init_mol2_fn;
  let init_mol2 = Mol2.read_one_from_file init_mol2_fn in
  let mol_name =
    sprintf "%s_l%d_s%d"
      (Mol.get_name centered_lig) lig_i start_i in
  let start_lig_fn =
    sprintf "%s/%s_l%d_s%d_start.bild"
      out_dir base_lig_fn lig_i start_i in
  let start_lig_mol2_fn =
    sprintf "%s/%s_l%d_s%d_start.mol2"
      out_dir base_lig_fn lig_i start_i in
  let best_lig_bild_fn =
    sprintf "%s/%s_l%d_s%d_best.bild"
      out_dir base_lig_fn lig_i start_i in
  let best_lig_mol2_fn =
    sprintf "%s/%s_l%d_s%d_best.mol2"
      out_dir base_lig_fn lig_i start_i in
  let stop_lig_fn =
    sprintf "%s/%s_l%d_s%d_stop.bild"
      out_dir base_lig_fn lig_i start_i in
  let ene_dump =
    let ene_out =
      let ene_fn =
        sprintf "%s/%s_l%d_s%d_E.txt"
          out_dir base_lig_fn lig_i start_i in
      open_out ene_fn in
    (* file format header comment *)
    fprintf ene_out "#frame E\n";
    later (fun () -> close_out ene_out);
    if ene_N > 0 then
      (* detailed prod. run: user want details about the E landscape *)
      (fun step ene ->
         if step mod ene_N = 0 then
           fprintf ene_out "%d %.3f\n"
             step ene)
    else
      (* quick prod. run: only log lower E *)
      (fun step ene ->
         if ene < !best_E then
           fprintf ene_out "%d %.3f\n"
             step ene) in
  let traj_dump =
    if traj_N = 0 then
      (fun _i _mol -> ())
    else
      let traj_out =
        let traj_fn = sprintf "%s/%s_l%d_s%d.xyz"
            out_dir base_lig_fn lig_i start_i in
        open_out traj_fn in
      later (fun () -> close_out traj_out);
      if traj_N = 1 then
        (fun i mol -> Mol.xyz_dump traj_out i "" mol)
      else
        (fun i mol ->
           if i mod traj_N = 0 then
             Mol.xyz_dump traj_out i "" mol) in
  let rotate_bond =
    let rbf = !Params.rbond_flip_block in
    (fun i rng lig ->
       (* do not modify significantly the start conformer *)
       if i > 0 && i mod rbf = 0 then
         (* let () = Log.info "rbond flip %d" i in *)
         Mol.flip_rbond rng lig
       else
         (* let () = Log.info "rbond tweak %d" i in *)
         Mol.tweak_rbond rng lig) in
  let ene_intra = match ff with
    | NULL -> failwith "Lds.simulate_isolated_lig: FF = NULL"
    | _ -> ene_intra_builder
             ff true num_rbonds centered_lig maybe_rdkit_conf in
  (* keep best pose even in case of crash *)
  later (fun () ->
      let best_E_str = sprintf "E: %g" !best_E in
      let prev_E_str = sprintf "E: %g" !prev_E in
      Log.info "l%d s%d %s" lig_i start_i best_E_str;
      Mol.to_bild ~verbose:true ~style:Ptable.Cyan_H
        ~header:best_E_str best_lig_bild_fn !best_lig;
      (* save as mol2 *)
      let mol2 = Mol.update_mol2 init_mol2 !best_lig in
      Log.info "creating %s" best_lig_mol2_fn;
      Mol2.write_one_to_file best_lig_mol2_fn mol2;
      Mol.to_bild ~verbose:true ~style:Ptable.Pink_H
        ~header:prev_E_str stop_lig_fn !prev_lig;
      Log.(info "########## %s ##########" (color green "CHIMERA COMMAND"));
      (* just show in chimera lig_start and lig_best;
         lig_stop is saved in the output dir though *)
      eprintf "chimera %s %s\n%!" start_lig_fn best_lig_bild_fn
    );
  (try
     for i = 0 to nsteps - 1 do
       let lig' = Mol.copy !prev_lig in
       let just_rotated = rotate_bond i rng lig' in
       (* FF eval *)
       let curr_E = ene_intra lig' in
       (if i = 0 then
          let curr_E_str = sprintf "E: %g" curr_E in
          let mol2 = Mol.update_mol2 init_mol2 lig' in
          Log.info "creating %s" start_lig_mol2_fn;
          Mol2.write_one_to_file start_lig_mol2_fn mol2;
          Mol.to_bild ~verbose:true ~style:Ptable.Orange_H ~header:curr_E_str
            start_lig_fn lig'
       );
       (* Metropolis criterion *)
       (if (curr_E <= !prev_E) || (* always allow energy to decrease *)
           (* accept increase depending on dice roll and deltaE *)
           ((Random.State.float rng 1.0) < exp (-.(curr_E -. !prev_E) *.
                                                  beta_at_K)) then
          begin
            SW.process (Mol.get_rbond_sw lig' just_rotated) true;
            prev_E := curr_E; (* update E *)
            prev_lig := lig' (* update lig *)
          end
        else
          SW.process (Mol.get_rbond_sw lig' just_rotated) false
       );
       ene_dump  i !prev_E;
       traj_dump i !prev_lig;
       (* monitor E *)
       (if curr_E < !best_E then
          begin
            best_E := curr_E;
            best_lig := lig';
            if sequential_run then
              Log.info "f%#d l%d s%d minE: %.2f"
               i lig_i start_i curr_E
          end
       );
       (* update dr for each RotBond *)
       (if i > 0 && i mod rbonds_block_size = 0 then
          begin
            update_rbond_rot_params target_low target_high lig';
            if Mol.strange_rbond_acceptance_rate lig' then
              begin (* inspect rbonds rot. params *)
                Log.warn "mol:%s rbonds_AR: %s" mol_name (Mol.sprintf_rbonds_sw lig');
                Log.warn "mol:%s rbonds_dr: %s" mol_name (Mol.sprintf_rbonds_dr lig')
              end
          end
       );
     done
   with Mol.Ene_intra -> Log.error "Mol.Ene_intra %s" mol_name
  );
  later_is_now ()

let read_count = ref 0

let demux arr () =
  let n = A.length arr in
  if !read_count < n then
    let res = (!read_count, arr.(!read_count)) in
    incr read_count;
    res
  else
    raise Parany.End_of_input

let get_my_rank nprocs =
  if nprocs = 1 then
    0
  else
    Parany.get_rank ()

let fst3 (a, _, _) = a
let snd3 (_, b, _) = b

type bitmask = Whole_protein of Mol.BST.t
             | ROI_only of ROI.sphere

let compute_bitmask nprocs bitmask_fn overwrite_grids roi_or_prot_bst grid =
  let bitmask_dt, bitmask =
    if Sys.file_exists bitmask_fn && (not overwrite_grids) then
      let () = Log.info "loading %s" bitmask_fn in
      Utls.(time_it bitmask_from_file bitmask_fn)
    else
      let () = Log.info "creating %s" bitmask_fn in
      Utls.(time_it (fun x ->
          let res = match roi_or_prot_bst with
            | Whole_protein prot_bst -> bitmask_whole_protein nprocs prot_bst x
            | ROI_only roi -> bitmask_ROI_only nprocs roi x in
          Utls.bitmask_to_file bitmask_fn res;
          res) grid) in
  Log.info (match roi_or_prot_bst with
      | Whole_protein _ -> "whole prot. bitmask dt: %.2f (s)"
      | ROI_only _ -> "enlarged ROI bitmask dt: %.2f (s)"
    ) bitmask_dt;
  bitmask

let bitmask_to_bild color fn g bitmask =
  LO.with_out_file fn (fun out ->
      let m, n, p = Grid.get_dims g in
      (* grid points color *)
      fprintf out ".color %s\n" color;
      for i = 0 to m - 1 do
        let x = g.xs.(i) in
        for j = 0 to n - 1 do
          let y = g.ys.(j) in
          for k = 0 to p - 1 do
            let idx = Grid.idx_of_ijk g i j k in
            if Bitv.get bitmask idx then
              fprintf out ".sphere %g %g %g 0.1\n" x y g.zs.(k)
          done
        done
      done
    )

(* how many cores for each loop; from outer to inner *)
type par_scheme = { ligands: int;
                    starts: int }

let par_scheme_ncores p =
  p.ligands * p.starts

let parse_par_scheme s =
  match L.map int_of_string (S.split_on_char '/' s) with
  | [i; j] ->
    let numcores = Cpu.numcores () in
    let requested = i * j in
    if requested > numcores then
      (Log.fatal "Lds.parse_par_scheme: requested(%d) > available(%d)"
         requested numcores;
       exit 1)
    else
      { ligands = i; starts = j }
  | _ ->
    let () = Log.fatal "Lds.parse_par_scheme: cannot parse: %s" s in
    exit 1

let dump_atom_types fn types =
  Log.info "creating %s" fn;
  LO.with_out_file fn (fun out ->
      A.iteri (fun i (anum, pcharge) ->
          fprintf out "%d (%d, %g)\n" i anum pcharge
        ) types
    )

(* same number of molecules; same names
 * if not, we need to examine the diff *)
let check_mol2_pqrs_files mol2_fn pqrs_fn =
  let mol2_names = Mol2.mol_names_from_file mol2_fn in
  let pqrs_names = Pqrs.mol_names_from_file pqrs_fn in
  Log.info "%s: %d molecules" mol2_fn (SS.cardinal mol2_names);
  Log.info "%s: %d molecules" pqrs_fn (SS.cardinal pqrs_names);
  assert(SS.equal mol2_names pqrs_names)

(* number of atom types to handle per processor for optimal load balancing *)
let load_balance nprocs num_atom_types maybe_apt =
  match maybe_apt with
  | Some x ->
    if nprocs * x >= num_atom_types then
      A.make nprocs x
    else (* the user might want to underload cores to save RAM at runtime
            (in case of too many atom types) *)
      let n = int_of_float (ceil (float num_atom_types /. float x)) in
      Log.warn "Lds.load_balance: saving RAM?";
      A.make n x
  | None ->
    if nprocs = 1 then (* most efficient way for a single processor *)
      let () = Log.warn "nprocs=1; ignoring -apt" in
      [|num_atom_types|]
    else
      (* load all processors equaly, then round-robin until there are
         no more jobs *)
      let per_core = num_atom_types / nprocs in
      let res = A.make nprocs per_core in
      let rem = ref (num_atom_types - (nprocs * per_core)) in
      let i = ref 0 in
      while !rem > 0 do
        res.(!i) <- res.(!i) + 1;
        rem := !rem - 1;
        i := (!i + 1) mod nprocs
      done;
      res

type mode = Ligand_only_MC
          | Protein_ligand_MC
          | Heuristic_docking
          | Exhaustive_rigid_docking

(* return (dx, rotations) for exhaustive rigid ligand docking w/ [steps] in [roi] *)
let exhaustive_params_for_ROI (steps: int) (roi: ROI.sphere): float * Rot.t array =
  (* sample as much the translational and the rotational space;
     steps = num_trans * num_rots w/ num_trans = num_rots *)
  let num_trans_f = ceil (sqrt (float steps)) in
  let num_trans_i = Utls.int_of_float_exn num_trans_f in
  let target_side_dim = num_trans_f ** (1./.3.) in
  let x_min, x_max, _, _, _, _ = ROI.get_bounds roi in
  let dx = (x_max -. x_min) /. target_side_dim in
  let rotations = SO3.rotations num_trans_i in
  Log.info "exhaust. dock. params: dx=%g rots=%d" dx num_trans_i;
  (dx, rotations)

let parse_ext_params s =
  Scanf.sscanf s "%f,%d" (fun f i -> (f, i))

(* return a mol_name to (rot, pos) hashtable *)
let load_pos_rot_from_file fn =
  let res = Ht.create (LO.count fn) in
  LO.iter fn (fun line ->
      Scanf.sscanf line "%s %f %f %f %f %f %f"
        (* Optim.to_string wrote translational parameters
           before rotational ones! *)
        (fun mol_name x y z a b g ->
           if Ht.mem res mol_name then
             failwith
               (sprintf "Lds.load_rot_pos_from_file %s: duplicate name: %s"
                  fn mol_name)
           else
             Ht.add res mol_name (Rot.r_xyz a b g, V3.make x y z)
        )
    );
  res

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  Log.(set_prefix_builder short_prefix_builder);
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s\n  \
              -lig <FILE.mol2>: ligand\n  \
              [--xtal-lig <FILE.mol2>]: crystallographic (reference) ligand\n  \
              -rec <FILE.pqrs>: receptor protein\n  \
              [--par <int/int>]: parallelization scheme (L/S)\n  \
              (default=1/1)\n  \
              [-s <int>]: random seed\n  \
              [-steps <int>[k|M]]: maximum number of frames\n  \
              [-top <int>]: keep top-k scores (for --ext)\n  \
              [-starts <int>]: how many ligand starting positions\n  \
              [-ff]: FF eval. style: {BrG|BrL|Bst*} (default=*)\n  \
              [-T <float>]: temperature in Kelvin (default=%g)\n  \
              [--out-dir <string>]: create output files' directory\n  \
              [--ogrids]: overwrite existing FF grid cache files\n  \
              [--no-interp]: turn OFF FF interpolation (slow)\n  \
              [-f]: overwrite output dir, if any\n  \
              [--traj <int>]: write ligand's traj.xyz output every N frames \
              (default=OFF)\n  \
              [--rigid-ligand]: freeze ligand's rotatable bonds\n  \
              [--no-E-intra]: ignore lig_E_intra (w/ --dock OR --ext)\n  \
              [--intra-QM]: torchani ANI-2 QM approximation\n  \
              [--intra-NB]: lig_E_intra is UFF Non-Bonded interactions only\n  \
              [--intra-UFF]: lig_E_intra is rdkit UFF\n  \
              [--intra-MMFF]: lig_E_intra is rdkit MMFF94\n  \
              [--no-vdW-clash]: lig_E_intra is only vdW clashes disallowed\n  \
              [--no-flip]: disable rbond flip move in MC simulations\n  \
              [--ene <int>]: log E every N frames (default=only log lower)\n  \
              [--hard-ROI]: enforce ROI during simulation (default=false)\n  \
              -roi <FILE.bild>: ROI sphere in original PDB coordinates\n  \
              [--undock]: undocking experiment (no random ligand placement/ROI)\n  \
              [--relax-from <FILE.tsv>]: protein-ligand MC simulation\n  \
              starting from prior docking experiment result\n  \
              [--ligand-only]: ligand-only MC simulation (no protein)\n  \
              [--dock]: protein-ligand docking (ligand minimization in binding-site)\n  \
              [--ext <float,int>]: rigid ligand exhaustive docking params (dx, rots)\n  \
              [-opt [0:6]]: global optim. algo. (for --dock; default=5)\n  \
              [--no-compress]: use much more disk space (but faster simulation startup)\n  \
              [--rand-conf]: randomize all ligand dihedral angles prior to docking\n  \
              [-apt <int>]: atom types per thread (default=auto)\n  \
              [--less-charges]: reduce precision of partial charges (if too many ligands; EXPERIMENTAL)\n  \
              [--rescore]: turn ON final pose rescoring for each ligand\n  \
              detailed scores in <OUT_DIR>/detailed_scores.tsv\n  \
              [-v]: verbose/debug mode\n"
       Sys.argv.(0) Defaults.temperature;
     exit 1);
  let rand_conf = CLI.get_set_bool ["--rand-conf"] args in
  let verbose = CLI.get_set_bool ["-v"] args in
  let no_compress = CLI.get_set_bool ["--no-compress"] args in
  let undock = CLI.get_set_bool ["--undock"] args in
  let rot_pos_fn = CLI.get_string_opt ["--relax-from"] args in
  let relax, name2rot_pos = match rot_pos_fn with
    | None -> (false, Ht.create 0)
    | Some fn -> (true, load_pos_rot_from_file fn) in
  let ligand_only = CLI.get_set_bool ["--ligand-only"] args in
  (* exhaustive docking params *)
  let topk = CLI.get_int_def ["-top"] args 0 in
  let translation_step, rotations = match CLI.get_string_opt ["--ext"] args with
    | None -> (nan, [||])
    | Some s ->
      let dx, num_rots = parse_ext_params s in
      (dx, SO3.rotations num_rots) in
  let num_rotations = A.length rotations in
  let exhaustive_dock = num_rotations > 0 in
  let docking = CLI.get_set_bool ["--dock"] args in
  let enforce_ROI = CLI.get_set_bool ["--hard-ROI"] args in
  let traj_N = CLI.get_int_def ["--traj"] args 0 in
  let tweak_rbonds = CLI.get_reset_bool ["--rigid-ligand"] args in
  let ff = match (CLI.get_set_bool ["--no-E-intra"]   args,
                  CLI.get_set_bool ["--intra-NB"]     args,
                  CLI.get_set_bool ["--intra-UFF"]    args,
                  CLI.get_set_bool ["--intra-MMFF"]   args,
                  CLI.get_set_bool ["--no-vdW-clash"] args,
                  CLI.get_set_bool ["--intra-QM"]     args) with
  | (true,  false, false, false, false, false) -> NULL
  | (false, true , false, false, false, false) -> UFF_NB
  | (false, false, true , false, false, false) -> RDKIT_UFF
  | (false, false, false, true,  false, false) -> RDKIT_MMFF
  | (false, false, false, false, true , false) -> NO_VDW_CLASH
  | (false, false, false, false, false, true ) -> TORCHANI_QM
  | _ -> failwith "Lds.main: which lig_E_intra FF to use?" in
  let no_flip_move = CLI.get_set_bool ["--no-flip"] args in
  let ene_N = CLI.get_int_def ["--ene"] args 0 (* default: only lower E logged *) in
  let no_interp = CLI.get_set_bool ["--no-interp"] args in
  let mol2_lig_fn = CLI.get_string ["-lig"] args in
  let maybe_xtal_lig_mol2 = CLI.get_string_opt ["--xtal-lig"] args in
  let roi_fn = CLI.get_string ["-roi"] args in
  let base_lig_fn = Fn.basename mol2_lig_fn in
  let base_lig_dir = Fn.dirname mol2_lig_fn in
  let rec_fn = CLI.get_string ["-rec"] args in
  let base_rec_fn = Fn.basename rec_fn in
  let para = parse_par_scheme (CLI.get_string_def ["--par"] args "1/1") in
  let maybe_apt = CLI.get_int_opt ["-apt"] args in
  let nsteps = steps_of_string (CLI.get_string ["-steps"] args) in
  let starts = CLI.get_int_def ["-starts"] args 1 in
  let seed_stream = match CLI.get_int_opt ["-s"] args with
    | None -> RNG.(make (entropy_160b ()))
    | Some seed -> RNG.make [|seed|] in
  let temp = CLI.get_float_def ["-T"] args Defaults.temperature in
  let ff_eval_style = ff_style_of_string (CLI.get_string_def ["-ff"] args "Bst") in
  let out_dir_overwrite = CLI.get_set_bool ["-f"] args in
  let overwrite_grids = CLI.get_set_bool ["--ogrids"] args in
  let out_dir = match CLI.get_string_opt ["--out-dir"] args with
    | None -> mk_temp_dir ()
    | Some dir -> (if out_dir_overwrite
                   then recreate_dir
                   else create_dir) dir in
  let rescore = CLI.get_set_bool ["--rescore"] args in
  let maybe_rescore_out =
    if not rescore then
      None
    else
      let rescore_fn = sprintf "%s/detailed_scores.tsv" out_dir in
      let rescore_out = open_out rescore_fn in
      (* file format header *)
      fprintf rescore_out "#name\tE_inter\tRB\tHA\n";
      Some rescore_out in
  Flags.glob_opt := CLI.get_int_def ["-opt"] args 5;
  let less_charges = CLI.get_set_bool ["--less-charges"] args in
  CLI.finalize (); (* ------------------------------------------------------ *)
  let simulation =
    if ligand_only then
      Ligand_only_MC
    else if exhaustive_dock then
      Exhaustive_rigid_docking
    else if docking then
      Heuristic_docking
    else
      Protein_ligand_MC in
  let nprocs = par_scheme_ncores para in
  (if verbose then
     (Log.(set_log_level DEBUG);
      Flags.verbose := true)
  );
  if no_compress then
    Flags.no_compress := true
  ;
  if no_flip_move then Params.rbond_flip_block := max_int;
  (* margin allows E_inter to go down to zero BEFORE ligand escaping
     simulation box. This imposes a constraint on the "longest" ligand
     that can be handled in the simulation box *)
  let margin = Const.charged_cutoff *. 3.0 in
  let beta_at_K = beta temp in
  let lig_pqrs_fn      = sprintf "%s/%s.pqrs"       out_dir base_lig_fn in
  let xtal_lig_pqrs_fn = sprintf "%s/xtal_lig.pqrs" out_dir in
  let xtal_lig_mol2_fn = sprintf "%s/xtal_lig.mol2" out_dir in
  let atom_types_fn = sprintf "%s/atom_types.txt" out_dir in
  Utls.mol2pqrs nprocs mol2_lig_fn lig_pqrs_fn;
  check_mol2_pqrs_files mol2_lig_fn lig_pqrs_fn;
  let ligands = Mol.ligands_of_pqrs_file lig_pqrs_fn in
  let prot = Mol.receptor_of_pqrs_file rec_fn in
  let prot_init_center = Mol.get_center prot in
  (* keep initial ligand indexes *)
  let ligs_dir = sprintf "%s/ligs" out_dir in
  let _ = create_dir ligs_dir in
  let lig_mol2_fns = Mol2.read_all_raw ligs_dir mol2_lig_fn in
  let num_mol2_ligs = A.length lig_mol2_fns in
  let num_ligands = L.length ligands in
  (if num_mol2_ligs <> num_ligands then
     let () =
       Log.fatal "Lds.main: num_mol2_ligs(%d) <> num_ligands(%d)"
         num_mol2_ligs num_ligands in
     exit 1
  );
  let ligands =
    L.mapi (fun i x ->
        if requires_rdkit ff then
          (* each ligand is now isolated in its own mol2 file;
             for perf. reasons *)
          let conf = Rdkit.__init__ ~mol2:lig_mol2_fns.(i) ~i:0 () in
          (* why there is never a call to Rdkit.is_valid_MMFF ?! *)
          if Rdkit.is_valid_UFF conf () then
            Some (i, x, Some conf)
          else
            let () = Log.error "rdkit error on %s" lig_mol2_fns.(i) in
            None
        else
          (* we do not care if the molecule is valid according to RDKit
           * hey, if another chemoinfo toolkit could create a 3D conformer
           * add hydrogens and assign partial charges to it; the molecule
           * is probably not that wrong... *)
          Some (i, x, None)
      ) ligands in
  (* filter out too elongated ligands, compared to simulation box's margin *)
  let ligands =
    L.filter (function
        | None -> false
        | Some (_i, lig, _conf) ->
          try
            Mol.check_elongation_exn lig Const.charged_cutoff;
            true
          with Mol.Too_long ->
            (Log.error "Mol.Too_long %s" (Mol.get_name lig);
             false)
      ) ligands in
  let ligands = L.map Option.get ligands in
  let ligands =
    if less_charges then
      L.map (fun (x, mol, y) ->
          Mol.reduce_partial_charges_precision mol;
          (x, mol, y)
        ) ligands
    else
      ligands in
  (* for FF interpolation from grid points *)
  let atom_types =
    let typ2id = Mol.get_type2id_ht (L.map snd3 ligands) in
    L.iter (fun (_lig_i, lig, _conf) ->
        Mol.assign_FF_types typ2id lig
      ) ligands;
    Mol.types_array typ2id in
  let num_atom_types = A.length atom_types in
  Log.info "lig FF atom types: %d" num_atom_types;
  let apt = load_balance nprocs num_atom_types maybe_apt in
  Log.info "apt: [|%s|]" (Utls.string_of_array apt ";" string_of_int);
  dump_atom_types atom_types_fn atom_types;
  let lig_FF_uniq_atoms = Mol.get_all_vdW_elec_types atom_types in
  Log.info "uniq atoms, charges: %d, %d"
    (A.length (Mol.uniq_anums   atom_types))
    (A.length (Mol.uniq_charges atom_types));
  (* place protein *)
  let prot_bild_fn = sprintf "%s/%s.bild" out_dir base_rec_fn in
  let prot_pdb_fn  = sprintf "%s/%s.pdb"  out_dir base_rec_fn in
  let prot_box_bild_fn = sprintf "%s/%s.box.bild" out_dir base_rec_fn in
  let box = preprocess_protein prot_bild_fn prot_box_bild_fn margin prot in
  (* !!! create BST AFTER preprocess_protein !!! *)
  (* because preprocess_protein moves the protein
     at the center of the simulation box *)
  let prot_bst = index_protein_atoms prot in
  let prot_prod_center = Mol.get_center prot in
  (* how the protein was translated into the simulation box *)
  let init_to_prod = V3.diff prot_prod_center prot_init_center in
  BatOption.may (fun fn ->
      Log.info "translating xtal lig to %s" xtal_lig_mol2_fn;
      Utls.mol2pqrs 1 fn xtal_lig_pqrs_fn;
      let xtal_lig_mol2 = Mol2.read_one_from_file fn in
      match Mol.ligands_of_pqrs_file xtal_lig_pqrs_fn with
      | [xtal_lig] ->
        Mol.translate_by xtal_lig init_to_prod;
        let xtal_lig_mol2' = Mol.update_mol2 xtal_lig_mol2 xtal_lig in
        Mol2.write_one_to_file xtal_lig_mol2_fn xtal_lig_mol2'
      | _ ->
        let () = Log.fatal "no or several molecules in %s" fn in
        exit 1
    ) maybe_xtal_lig_mol2;
  Mol.to_pdb_file prot_pdb_fn prot;
  (if undock then
     (* in an undocking experiment, ligands need to move like the protein
        which was moved into the simulation box *)
     L.iter (fun (_lig_i, lig, _conf) ->
         (* the conf. coordinates are not updated;
          * they will be updated only before FF evaluation, if needed *)
         Mol.translate_by lig init_to_prod
       ) ligands
  );
  let roi = ROI.translate (ROI.from_bild roi_fn) init_to_prod in
  (* dump ROI to disk *)
  let roi_bild_fn        = sprintf "%s/ROI.bild" out_dir in
  let dotted_roi_bild_fn = sprintf "%s/dot_ROI.bild" out_dir in
  ROI.to_bild        roi_bild_fn        roi;
  ROI.dotted_to_bild dotted_roi_bild_fn roi;
  let ff0 = match ff_eval_style with
    | Brute_global -> Mol.ene_inter_UFF_global_brute
    | Brute_local -> Mol.ene_inter_UFF_shifted_brute
    | Bst_local -> Mol.ene_inter_UFF_shifted_bst prot_bst in
  let grid = Grid.from_box Params.grid_step box in
  Log.info "%s" (Grid.str_peek grid);
  let grid_bild_fn = sprintf "%s/%s.grid.bild" out_dir base_rec_fn in
  let prot_solv_bitmask_fn = sprintf "%s/%s.solv_bits.bild" out_dir base_rec_fn in
  let prot_vdW_bitmask_fn  = sprintf "%s/%s.vdW_bits.bild"  out_dir base_rec_fn in
  Log.info "grid file: %s" grid_bild_fn;
  Grid.to_bild grid_bild_fn grid;
  Log.info "computing bitmask...";
  (* pre-compute FF components *)
  let ff1 =
    if no_interp then
      ff0 prot
    else
      (* this file is persistent across runs; unless --ogrids *)
      let bitmask_fn = sprintf "%s.bitmask" rec_fn in
      let bitmask =
        compute_bitmask
          nprocs bitmask_fn overwrite_grids (ROI_only roi) grid in
      let ff_comps =
        let base_lig_fn = sprintf "%s/%s" base_lig_dir base_lig_fn in
        pre_compute_FF_components
          nprocs apt base_lig_fn overwrite_grids grid bitmask prot prot_bst lig_FF_uniq_atoms in
        (* pre_calculate_FF_components nprocs csize grid bitmask prot_bst prot lig_FF_uniq_atoms in *)
      Mol.ene_inter_UFF_interp grid ff_comps in
  Log.info "protein 1st solvent shell bitmask...";
  let dt, prot_solvent_shell = Utls.time_it (first_solvent_shell grid) prot in
  Log.info "took %.3fs" dt;
  let ene_rescore = Mol.ene_inter_UFF_shifted_bst prot_bst prot in
  Log.info "protein vdW volume...";
  let dt, prot_vdW_bitmask =
    Utls.time_it (vdW_volume grid) prot in
  Log.info "took %.3fs" dt;
  (if verbose then
     let () = Log.info "creating %s" prot_solv_bitmask_fn in
     bitmask_to_bild "blue" prot_solv_bitmask_fn grid prot_solvent_shell;
     let () = Log.info "creating %s" prot_vdW_bitmask_fn in
     bitmask_to_bild "grey" prot_vdW_bitmask_fn grid prot_vdW_bitmask;
     let _radii, prot_dotted_SAS = Mol.compute_SAS_dot_surface prot in
     let prot_SAS_bild_fn = sprintf "%s/%s_prot.SAS.bild" out_dir base_lig_fn in
     write_SAS_to_bild_file prot_SAS_bild_fn prot_dotted_SAS
  );
  (* each run must be statistically independent from the others *)
  let rngs =
    let num_ligs = 1 + L.max (L.map fst3 ligands) in
    A.init (starts * num_ligs) (fun _i -> RNG.split seed_stream) in
  let sequential = (nprocs = 1) in
  let scores_out =
    let scores_fn = sprintf "%s/docking_scores.tsv" out_dir in
    open_out scores_fn in
  let pos_rot_out =
    let rot_pos_fn = sprintf "%s/pos_rot.tsv" out_dir in
    open_out rot_pos_fn in
  (* parallel loop over all ligands *)
  Parany.Parmap.pariter para.ligands (fun (lig_i, lig, conf) ->
      let mol_name = Mol.get_name lig in
      let rng0 = rngs.(starts * lig_i) in
      (* place ligand *)
      let lig_bild_fn     = sprintf "%s/%s_l%d.bild"     out_dir base_lig_fn lig_i in
      let lig_box_bild_fn = sprintf "%s/%s_l%d.box.bild" out_dir base_lig_fn lig_i in
      (* needed by undocking *)
      let orig_lig = Mol.copy lig in
      let num_rbonds = Mol.num_rbonds orig_lig in
      (if num_rbonds > 7 || num_rbonds = 0
       then Log.warn "%s: %d rbonds"
       else Log.info "%s: %d rbonds") mol_name num_rbonds;
      (* center ligand *)
      preprocess_ligand lig_bild_fn lig_box_bild_fn lig_i lig;
      let centered_ligs =
        A.init starts
          (if rand_conf then
             (fun _ -> Mol.randomize_conformer lig rng0)
           else
             (fun _ -> Mol.copy lig)) in
      Log.info "lig %d SAS dots..." lig_i;
      let lig_SAS_bild_fn = sprintf "%s/%s_l%d.SAS.bild" out_dir base_lig_fn lig_i in
      (* only compute and store SAS of 1st one *)
      let _radii, lig_dotted_SAS = Mol.compute_SAS_dot_surface centered_ligs.(0) in
      write_SAS_to_bild_file lig_SAS_bild_fn lig_dotted_SAS;
      let start_poses =
        if undock || exhaustive_dock then
          (* even during an undocking simulation, we might want
           * to see several independent trials *)
          A.make starts (Rot.id (), Mol.get_center orig_lig)
        else if docking then
          (* random point inside ROI and fully random rotation *)
          drop_ligand_in_ROI rng0 roi starts
        else if relax then
          (* use best pose found by a prior run of
             rigid-ligand rigid-protein docking no_lig_Eintra *)
          A.make starts (Ht.find_exn (fun x -> x) name2rot_pos mol_name)
        else
          place_ligand_in_ROI rng0 prot_bst prot roi centered_ligs starts in
      dump_start_ligands out_dir base_lig_fn lig_i centered_ligs start_poses;
      (* SIMULATION LOOP ------------------------------------------------------ *)
      Parany.Parmap.array_pariteri para.starts (fun start_i (start_rot, start_pos) ->
          let rng = rngs.(starts * lig_i + start_i) in
          let elapsed, () =
            Utls.time_it
              (match simulation with
               | Ligand_only_MC ->
                 simulate_isolated_lig ff sequential rng nsteps traj_N ene_N
                   beta_at_K out_dir base_lig_fn lig_i start_i centered_ligs.(start_i)
               | Heuristic_docking ->
                 minimize_ligand tweak_rbonds ff
                   nsteps traj_N
                   prot_bild_fn
                   roi ff1 ene_rescore
                   scores_out maybe_rescore_out
                   out_dir pos_rot_out base_lig_fn lig_i
                   start_rot start_pos start_i centered_ligs.(start_i)
               | Exhaustive_rigid_docking ->
                 exhaustive_docking
                   topk
                   ff
                   translation_step
                   rotations
                   prot_bild_fn
                   roi grid prot_vdW_bitmask
                   ff1
                   out_dir scores_out pos_rot_out base_lig_fn lig_i
                   centered_ligs.(start_i)
               | Protein_ligand_MC ->
                 simulate_lig tweak_rbonds ff sequential
                   enforce_ROI rng nsteps traj_N
                   ene_N prot_bild_fn
                   roi ff1 ene_rescore
                   beta_at_K
                   scores_out maybe_rescore_out
                   out_dir base_lig_fn lig_i
                   start_rot start_pos start_i centered_ligs.(start_i)) conf in
          Log.info "%d frames in %.2f (s) @ %.2f (Hz)"
            nsteps elapsed ((float (nsteps * starts)) /. elapsed)
        ) start_poses
    ) ligands;
  close_out scores_out;
  close_out pos_rot_out;
  BatOption.may close_out maybe_rescore_out

let () = main ()
