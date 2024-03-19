(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

module FSet = BatSet.Float
module ISet = BatSet.Int
module Rdkit = Rdkit.Rdkit

open Printf

(* to index all protein atoms *)
module BST = Bst.Bisec_tree.Make(Atom)

type t = { xs: float array; (* atoms X coords *)
           ys: float array; (* atoms Y coords *)
           zs: float array; (* atoms Z coords *)
           mutable center: V3.t;
           q_a: float array; (* partial charges *)
           r_a: float array; (* vdW radii *)
           elt_a: int array; (* atomic numbers *)
           (* dists_a, rbonds_a and rgroups_a will be left uninitialized
              for a protein; only ligands have those *)
           dists_a: int array; (* topological distance matrix *)
           interacting_a: bool array array; (* are two atoms UFF non-bonded
                                               interacting (>= 3 bonds away) *)
           rbonds_a: Bond.t array; (* rotatable bonds *)
           rbonds_dr_a: float array; (* rot. bonds rotational parameters *)
           rbonds_sw_a: SW.t array; (* rot. bonds current acceptance rate *)
           rgroups_a: (int array) array; (* rotatable atom indexes for each RotBond *)
           (* t_a: uninitialized until assign_FF_types is called *)
           t_a: int array; (* FF atom types *)
           name: string } (* /!\ RENAMING MOLECULES IS _FORBIDDEN_ /!\ *)

let create mol_name num_atoms: t =
  let xs     = A.make num_atoms 0.0  in
  let ys     = A.make num_atoms 0.0  in
  let zs     = A.make num_atoms 0.0  in
  let center = V3.origin             in
  let q_a    = A.make num_atoms 0.0  in
  let r_a    = A.make num_atoms 0.0  in
  let elt_a  = A.make num_atoms 0    in
  let dists_a       = [||] in
  let interacting_a = [||] in
  let rbonds_a      = [||] in
  let rbonds_dr_a   = [||] in
  let rbonds_sw_a   = [||] in
  let rgroups_a     = [||] in
  let t_a = A.make num_atoms (-1) in
  { xs; ys; zs; center; q_a; r_a; elt_a; dists_a; interacting_a; rbonds_a;
    rgroups_a; rbonds_dr_a; rbonds_sw_a; t_a; name = mol_name }

(* same molecule, but 3D coordinates in a fresh array (for later modif.) *)
(* !!! ALL MODIFIABLE ARRAYS NEED TO BE COPIED TO FRESH ONES !!! *)
let copy m =
  { m with xs = A.copy m.xs ;
           ys = A.copy m.ys ;
           zs = A.copy m.zs ;
           rbonds_dr_a = A.copy m.rbonds_dr_a ;
           rbonds_sw_a = A.copy m.rbonds_sw_a }

let num_atoms m: int =
  A.length m.r_a

let num_heavy_atoms m: int =
  A.count_matching (fun anum -> anum > 1) m.elt_a

let num_rbonds m: int =
  A.length m.rbonds_a

let get_name m =
  m.name

(* set_* are bound-checked (infrequent, only at molecule's creation time);
   get_* are not (frequent, during FF evaluation) *)

let uget = A.unsafe_get
let uset = A.unsafe_set

let get_radius m i: float =
  uget m.r_a i

let set_radius m i r: unit =
  A.set m.r_a i r

let get_charge m i: float =
  uget m.q_a i

let set_charge m i q: unit =
  A.set m.q_a i q

let get_anum m i: int =
  uget m.elt_a i

let get_ff_type m i: int * float =
  (get_anum m i, get_charge m i)

let get_typ m i: int =
  uget m.t_a i

let get_symbol m i: string =
  Ptable.symbol_of_anum (get_anum m i)

let set_anum m i anum: unit =
  A.set m.elt_a i anum

let get_xs m: float array = m.xs
let get_ys m: float array = m.ys
let get_zs m: float array = m.zs

let get_xyz m i: V3.t =
  V3.make
    (uget m.xs i)
    (uget m.ys i)
    (uget m.zs i)

let get_rbond m i: Bond.t =
  uget m.rbonds_a i

let get_rbond_dr m i: float =
  uget m.rbonds_dr_a i

let set_rbond_dr m i x: unit =
  uset m.rbonds_dr_a i x

let get_rbond_sw m i: SW.t =
  uget m.rbonds_sw_a i

(* check if there is any RotBond with 0.0 or 1.0 acceptance rate *)
let strange_rbond_acceptance_rate m: bool =
  A.exists (fun sw ->
      let x = SW.get_ratio sw in
      (x = 0.0) || (x = 1.0) || (classify_float x = FP_nan)
    ) m.rbonds_sw_a

(* all rbonds_dr as one string *)
let sprintf_rbonds_dr m =
  Utls.string_of_array m.rbonds_dr_a " " (sprintf "%.2f")

(* all acceptance rates as one string *)
let sprintf_rbonds_sw m =
  Utls.string_of_array (A.map SW.get_ratio m.rbonds_sw_a) " " (sprintf "%.2f")

let get_rgroup m i: int array =
  uget m.rgroups_a i

(* distance in bonds between two atoms
   (topological distance on the molecular graph) *)
let get_topo_dist m i j: int =
  uget m.dists_a (i + (j * (num_atoms m)))

(* number of Hs attached to each atom (0 for each H; >=0 for each HA *)
let hydrogen_counts m: int array =
  let n = num_atoms m in
  A.init n (fun i ->
      if get_anum m i = 1 then
        0 (* by definition, an H in a small molecule has no H attached
             (exception: H2), but H2 isn't a ligand so we should not
             encounter it *)
      else
        let res = ref 0 in
        for j = 0 to n - 1 do
          if i <> j && get_anum m j = 1 && get_topo_dist m i j = 1 then
            incr res
        done;
        !res
    )

(* number of HAs attached to each atom *)
let heavy_counts m: int array =
  let n = num_atoms m in
  A.init n (fun i ->
      if get_anum m i = 1 then
        1 (* by definition, an H is attached to a single HA
             (exception: H2), but H2 isn't a ligand so we should not
             encounter it *)
      else
        let res = ref 0 in
        for j = 0 to n - 1 do
          if i <> j && get_anum m j > 1 && get_topo_dist m i j = 1 then
            incr res
        done;
        !res
    )

(* list rotatable bonds which only move hydrogens;
   e.g. -CH3, -NH2, -OH, etc. *)
let non_conformer_changing_rbonds m: Bond.t array =
  let hydrogens = hydrogen_counts m in
  let heavies = heavy_counts m in
  A.filter (fun rbond ->
      let j, k = Bond.to_pair rbond in
      (* single bonds to H are not supposed to be counted as rbond *)
      assert(get_anum m j > 1 &&
             get_anum m k > 1);
      (hydrogens.(j) >= 1 && heavies.(j) = 1) ||
      (hydrogens.(k) >= 1 && heavies.(k) = 1)
    ) m.rbonds_a

(* precompute matrix of Non-Bonded interacting atoms for E_intra *)
let create_interacting_a n dists_a: bool array array =
  A.init n (fun i ->
      A.init n (fun j ->
          dists_a.(i + (j * n)) >= 3
        )
    )

let get_interacting m i: bool array =
  uget m.interacting_a i

(* RMSD without optimal superposition;
   atom numbering must be the same in both molecules;
   WARNING: does _NOT_ take into account molecular symmetries;
   i.e. does NOT return the lowest RMSD;
   would be probably less bad if we were computing the HA
   RMSD only *)
let rmsd (m1: t) (m2: t): float =
  let n = num_atoms m1 in
  assert(n = num_atoms m2);
  let diff_squares =
    A.init n (fun i ->
        V3.mag2
          (V3.diff
             (get_xyz m1 i)
             (get_xyz m2 i)
          )
      ) in
  sqrt (A.favg diff_squares)

let heavy_atoms_rmsd (m1: t) (m2: t): float =
  let num_HA = num_heavy_atoms m1 in
  let sum = ref 0.0 in
  let n = num_atoms m1 in
  for i = 0 to n - 1 do
    if get_anum m1 i > 1 then
      sum := !sum +. V3.(mag2
                           (diff
                              (get_xyz m1 i)
                              (get_xyz m2 i)))
  done;
  sqrt (!sum /. (float num_HA))

let get_all_atom_coords m: V3.t array =
  A.init (num_atoms m) (get_xyz m)

let get_all_atom_radii m: float array =
  m.r_a

let get_all_charges m: float array =
  m.q_a

(* reduce the number of atom types we have to deal with
   when screening _many_ ligands *)
let reduce_partial_charges_precision m: unit =
  let n = num_atoms m in
  for i = 0 to n - 1 do
    m.q_a.(i) <- Utls.reduce_precision m.q_a.(i)
  done

let get_anums m: int array =
  m.elt_a

(* update one 3D coordinate *)
let set_xyz m i v: unit =
  V3.(uset m.xs i v.x;
      uset m.ys i v.y;
      uset m.zs i v.z)

(* like [set_xyz] above, but does not require a V3.t parameter *)
let set_x_y_z m i x y z: unit =
  uset m.xs i x;
  uset m.ys i y;
  uset m.zs i z

(* TO PRE-COMPUTE ALL FF INTERACTION ENERGY VALUES AT GRID POINTS ---------- *)

(* assign an identifier to each (anum, pcharge) pair *)
let get_type2id_ht (mols: t list): (int * float, int) Ht.t =
  (* list of unique (anum, pcharge) pairs *)
  let ht = Ht.create 1013 in
  let j = ref 0 in
  L.iter (fun mol ->
      let n = num_atoms mol in
      for i = 0 to n - 1 do
        let typ = get_ff_type mol i in
        if not (Ht.mem ht typ) then
          let () = Ht.add ht typ !j in
          incr j
      done
    ) mols;
  ht

let get_all_vdW_elec_types (uniq_types: (int * float) array): t array =
  A.mapi (fun i (anum, pcharge) ->
      (* single atom dummy molecule *)
      let mol = create (sprintf "t%d" i) 1 in
      set_x_y_z  mol 0 0. 0. 0.;
      set_charge mol 0 pcharge;
      set_radius mol 0 0.0;
      set_anum   mol 0 anum;
      mol
    ) uniq_types

(* all uniq anums -> single atom virtual molecules w/ only a vdW interaction;
   partial charge is zeroed *)
let uniq_anums (uniq_types: (int * float) array): t array =
  let uniq_anums =
    let anums = (A.map fst uniq_types) in
    ISet.(to_array (of_array anums)) in
  A.map (fun anum ->
      (* single atom dummy molecule w/o partial charge *)
      let mol = create (sprintf "vdW:%d" anum) 1 in
      set_x_y_z  mol 0 0. 0. 0.;
      set_charge mol 0 0.;
      set_radius mol 0 Ptable.vdW_radii.(anum);
      set_anum   mol 0 anum;
      mol
    ) uniq_anums

(* all uniq partial charges -> single atom virtual molecules w/ only
   an electrostatic interaction; vdW parameters are null *)
let uniq_charges (uniq_types: (int * float) array): t array =
  let uniq_charges =
    let charges = (A.map snd uniq_types) in
    FSet.(to_array (of_array charges)) in
  A.map (fun charge ->
      (* single atom dummy molecule w/o vdW parameters *)
      let mol = create (sprintf "q:%f" charge) 1 in
      set_x_y_z  mol 0 0. 0. 0.;
      set_charge mol 0 charge;
      set_radius mol 0 0.;
      set_anum   mol 0 0; (* anum=0: virtual element w/o vdW params *)
      mol
    ) uniq_charges

(* ------------------------------------------------------------------------- *)

let x_min_max m: float * float =
  A.min_max m.xs

let y_min_max m: float * float =
  A.min_max m.ys

let z_min_max m: float * float =
  A.min_max m.zs

(* mean_xyz *)
let get_center m: V3.t =
  m.center

let update_center m =
  m.center <- V3.{ x = A.favg m.xs ;
                   y = A.favg m.ys ;
                   z = A.favg m.zs }

(* smallest atom radius *)
let r_min m: float =
  A.min m.r_a

(* largest atom radius *)
let r_max m: float =
  A.max m.r_a

(* only up to atom coordinates; a ligand's pqrs file will contain
   additional info *)
let pqrs_read_header_coords input =
  let header = input_line input in
  let num_atoms, num_rbonds, mol_name = Pqrs.parse_header header in
  (* Log.info "m:%s atoms:%d" mol_name num_atoms; *)
  let mol = create mol_name num_atoms in
  for i = 0 to num_atoms - 1 do
    let x, y, z, q, r, elt = Pqrs.parse_body (input_line input) in
    set_charge mol i q;
    set_radius mol i r;
    set_anum   mol i (Ptable.anum_of_symbol elt);
    A.set mol.xs i x;
    A.set mol.ys i y;
    A.set mol.zs i z
  done;
  update_center mol;
  (mol, num_rbonds)

let ligand_pqrs_read_one input =
  let mol, num_rbonds = pqrs_read_header_coords input in
  let n = num_atoms mol in
  let rot_bonds = A.make num_rbonds Bond.dummy in
  let rot_groups = A.make num_rbonds [||] in
  for i = 0 to num_rbonds - 1 do
    let axis, in_group = Pqrs.parse_rot_bond n input in
    rot_bonds.(i) <- axis;
    rot_groups.(i) <- Pqrs.rgroup_of_bools axis in_group
  done;
  let dist_mtx = Pqrs.parse_dist_matrix n input in
  let interacting_a = create_interacting_a n dist_mtx in
  { mol with rbonds_a = rot_bonds;
             rbonds_dr_a = A.make num_rbonds Params.max_rbond_rot;
             rbonds_sw_a = A.init num_rbonds (fun _i -> SW.create Params.block_size);
             rgroups_a = rot_groups;
             dists_a = dist_mtx;
             interacting_a = interacting_a }

(* reset all that was learned about those RotBonds flexibility *)
let reset_rbonds_params lig =
  let n = num_rbonds lig in
  for i = 0 to n - 1 do
    lig.rbonds_dr_a.(i) <- Params.max_rbond_rot;
    lig.rbonds_sw_a.(i) <- SW.create Params.block_size
  done

(* read in a protein .pqrs file *)
let receptor_of_pqrs_file (fn: string): t =
  LO.with_in_file fn (fun input ->
      let mol, _num_rbonds = pqrs_read_header_coords input in
      mol
    )

(* read all ligands from a .pqrs file; cf. data/caffeine.pqrs *)
let ligands_of_pqrs_file (fn: string): t list =
  LO.with_in_file fn (fun input ->
      let res, exn = L.unfold_exn (fun () -> ligand_pqrs_read_one input) in
      (if exn <> End_of_file then
         raise exn
      );
      (* force molecule names to be uniq: molecules were initially read from
         a .mol2 file; IT IS NOT ALLOWED TO RENAME THEM *)
      let ht = Ht.create (L.length res) in
      let dups = ref 0 in
      L.iter (fun mol ->
          let name = get_name mol in
          if Ht.mem ht name then
            let _ = incr dups in
            Log.fatal "Mol.ligands_of_pqrs_file: fn: %s duplicate: %s" fn name
          else
            Ht.add ht name ()
        ) res;
      if !dups > 0 then exit 1;
      res
    )

let first_ligand_of_pqrs_file (fn: string): t =
  LO.with_in_file fn ligand_pqrs_read_one

(* read all ligands from a .mol2 file *)
let ligands_of_mol2_file (nprocs: int) (mol2_fn: string): t list =
  let tmp_pqrs_fn= Utls.temp_file "mmcld_" ".pqrs" in
  (* convert to .pqrs file *)
  Utls.mol2pqrs nprocs mol2_fn tmp_pqrs_fn;
  (* read .pqrs file *)
  let ligs = ligands_of_pqrs_file tmp_pqrs_fn in
  Sys.remove tmp_pqrs_fn;
  ligs

(* initialize m.t_a *)
let assign_FF_types type2id m =
  let n = num_atoms m in
  for i = 0 to n - 1 do
    let id = Ht.find type2id (get_ff_type m i) in
    m.t_a.(i) <- id
  done

let types_array ht =
  let n = Ht.length ht in
  let res = A.make n (-1, 0.0) in
  Ht.iter (fun anum_pcharge id ->
      res.(id) <- anum_pcharge
    ) ht;
  res

(* output for UCSF chimera *)
let to_bild
    ?(verbose = false)
    ?(frame = false)
    ?(style = Ptable.White_H)
    ?(header = "")
    (fn: string) (m: t): unit =
  let coloring out =
    Ptable.(match style with
        | No_color -> (fun _anum -> ())
        | White_H -> (* default *)
          (fun anum ->
             fprintf out ".color %s\n"
               (color_of_anum anum))
        | Pink_H -> (* last ligand *)
          (fun anum ->
             fprintf out ".color %s\n"
               (color_of_anum_pink_H anum))
        | Cyan_H -> (* minimized ligand *)
          (fun anum ->
             fprintf out ".color %s\n"
               (color_of_anum_cyan_H anum))
        | Orange_H -> (* start ligand *)
          (fun anum -> fprintf out ".color %s\n"
              (color_of_anum_orange_H anum))
      ) in
  if verbose then
    Log.info "write %s to %s" m.name fn
  ;
  LO.with_out_file fn (fun out ->
      if header <> "" then
        fprintf out ".comment %s\n" header;
      if frame then
        begin
          (* cartesian frame *)
          fprintf out ".color red\n";
          fprintf out ".arrow 0. 0. 0. 1. 0. 0.\n";
          fprintf out ".color green\n";
          fprintf out ".arrow 0. 0. 0. 0. 1. 0.\n";
          fprintf out ".color blue\n";
          fprintf out ".arrow 0. 0. 0. 0. 0. 1.\n";
          fprintf out ".color white\n"
        end;
      let n = num_atoms m in
      for i = 0 to n - 1 do
        let r = get_radius m i in
        if r > 0.0 then
          let p    = get_xyz  m i in
          let anum = get_anum m i in
          coloring out anum;
          V3.(fprintf out ".sphere %g %g %g %g\n" p.x p.y p.z r)
      done
    )

(* append to an already opened .xyz file *)
let xyz_dump out counter str_data m =
  let n = num_atoms m in
  (* num_atoms line *)
  fprintf out "%d\n" n;
  (* comment line: ^frame_number<TAB>mol_name<TAB>str_data\n *)
  fprintf out "%d\t%s\t%s\n" counter (get_name m) str_data;
  (* ^elt x y z$ lines *)
  for i = 0 to n - 1 do
    fprintf out "%s" (get_symbol m i);
    V3.xyz_dump out (get_xyz m i)
  done

let to_xyz_file fn m =
  LO.with_out_file fn (fun out ->
      xyz_dump out 0 "" m
    )

(* move/update coordinates of original Mol2.t object *)
let update_mol2 mol2 m =
  let atoms' =
    A.mapi (fun i a ->
        Mol2.atom_move a
          (uget m.xs i)
          (uget m.ys i)
          (uget m.zs i)
      ) (Mol2.get_atoms mol2) in
  Mol2.move mol2 atoms'

let append_to_open_pdb_file out i m =
  fprintf out "MODEL     %04d\n" i;
  let n = num_atoms m in
  for i = 0 to n - 1 do
    let elt          = get_symbol m i in
    let V3.{x; y; z} = get_xyz    m i in
    fprintf out "ATOM  % 5d%4s  REC A   1    %8.3f%8.3f%8.3f  1.00  0.00%12s\n"
      i elt x y z elt
  done;
  fprintf out "ENDMDL\n"

let to_pdb_file fn m =
  LO.with_out_file fn (fun out ->
      (* model numbers starts at 1 in PDB format *)
      append_to_open_pdb_file out 1 m;
      fprintf out "END\n"
    )

(* sphere completely containing all ligand atom centers.
   Since we already know the ligand's geometric center, this is
   cheaper to compute O(n) compared to the conformer's diameter O(n^2). *)
(* englobing_sphere_radius *)
let radius m: float =
  let n = num_atoms m in
  let center = get_center m in
  let maxi = ref 0.0 in
  for i = 0 to n - 1 do
    maxi := max !maxi (0.01 +. V3.dist center (get_xyz m i))
  done;
  !maxi

(* in case the current conformer's englobing sphere radius
   is greater than the cutoff *)
exception Too_long

let check_elongation_exn mol max_radius =
  if radius mol > max_radius then
    raise Too_long

let translate_by m (V3.{ x; y; z} as v3): unit =
  let n = num_atoms m in
  for i = 0 to n - 1 do
    A.unsafe_set m.xs i (A.unsafe_get m.xs i +. x);
    A.unsafe_set m.ys i (A.unsafe_get m.ys i +. y);
    A.unsafe_set m.zs i (A.unsafe_get m.zs i +. z)
  done;
  m.center <- V3.add m.center v3

(* CONSTRAINT: ligand must already be centered at origin *)
let centered_rotate m r: unit =
  let n = num_atoms m in
  for i = 0 to n - 1 do
    set_xyz m i (Rot.rotate r (get_xyz m i))
  done

(* rotate rotatable bond [i] by [alpha] (in radians) *)
let rotate_bond m i alpha: unit =
  let rbond  = get_rbond  m i in
  let rgroup = get_rgroup m i in
  let center = get_xyz m rbond.right in
  (* rotation that will be applied to all members of the group *)
  let rot =
    let axis =
      V3.(normalize
            (diff center (get_xyz m rbond.left))) in
    Rot.of_axis_angle axis alpha in
  (* rotate rotatable atoms of this group *)
  A.iter (fun i ->
      (* - translate to origin
       * - rotate
       * - translate back
       * - update impacted coords *)
      set_xyz m i (V3.add
                     (Rot.rotate rot (V3.diff (get_xyz m i) center))
                     center)
    ) rgroup;
  (* the geometric center needs updating *)
  update_center m

(* try to locally optimize one dihedral *)
(* NOT ROBUST TO MOLECULE W/O ANY ROTBOND *)
let tweak_rbond rng m: int =
  let i = Random.State.int rng (num_rbonds m) in
  rotate_bond m i (Move.rand_angle rng (get_rbond_dr m i));
  i

(* try to reset one dihedral; ligand conformer search *)
(* THIS MOVE SHOULD BE RARELY USED; if used too often, it will
   mangle with the acceptance ratio statistics *)
(* NOT ROBUST TO MOLECULE W/O ANY ROTBOND *)
let flip_rbond rng m: int =
  let i = Random.State.int rng (num_rbonds m) in
  rotate_bond m i (Move.rand_angle rng Params.max_rbond_flip);
  i

(* prior to docking, it might be interesting to fully
   randomize dihedrals of the starting conformer;
   combined w/ random starts, this might force optimization
   to explore more the search space *)
let randomize_conformer centered_lig rng =
  let lig = copy       centered_lig in
  let n   = num_rbonds centered_lig in
  let angles = A.init n (fun _ -> Move.rand_angle rng Math.pi) in
  A.iteri (rotate_bond lig) angles;
  (* forbid too elongated ligand *)
  check_elongation_exn lig Const.charged_cutoff;
  lig

(* return a rotated copy of ligand [l];
 * CONSTRAINT: [l] must already be centered at the origin *)
let centered_rotate_copy l r: t =
  let lig' = copy l in
  centered_rotate lig' r;
  lig'

let rotate_then_translate_copy lig rot pos =
  let lig' = centered_rotate_copy lig rot in
  translate_by lig' pos;
  lig'

(* WARNING: less efficient than centered_rotate_by *)
let inplace_rotate m r: unit =
  let c = get_center m in
  (* center *)
  translate_by m (V3.neg c);
  (* rotate *)
  centered_rotate m r;
  (* put back in place *)
  translate_by m c

(* same as [inplace_rotate] but using quaternion [q] *)
let inplace_rotate_quat m q =
  inplace_rotate m (Rot.of_quat q)

(* [p] will become the new geometric center *)
let translate_to m p: unit =
  translate_by m (V3.diff p m.center);
  m.center <- p

let translate_copy_to m p: t =
  let lig' = copy m in
  translate_to lig' p;
  lig'

(* translate to the origin so that the molecule is ready to rotate *)
let center m: unit =
  let mean_xyz = get_center m in
  translate_by m (V3.neg mean_xyz);
  m.center <- V3.origin

(* center, then rotate, then translate by delta *)
let center_rotate_translate_copy lig rot pos =
  let lig' = copy lig in
  center lig';
  centered_rotate lig' rot;
  translate_by lig' pos;
  lig'

(* two corners of box *)
let bounding_box m =
  (* this is not the tighest box, but we don't care *)
  let max_rad = r_max m in
  let min_x, max_x = x_min_max m in
  let min_y, max_y = y_min_max m in
  let min_z, max_z = z_min_max m in
  Bbox.create_6f
    (min_x -. max_rad) (min_y -. max_rad) (min_z -. max_rad)
    (max_x +. max_rad) (max_y +. max_rad) (max_z +. max_rad)

(* is the point inside a sphere? *)
let atom_too_near atom_i dot atom_index_center_radii =
  A.exists (fun (j, center, radius) ->
      (atom_i <> j) && (V3.dist2 dot center < radius *. radius)
    ) atom_index_center_radii

(* is any atom [a_i] of [mol] such that
   dist(a_i, query) < threshold
   All molecule atom coordinates must have been previously obtained
   using [get_all_atom_coords mol]. *)
let any_atom_nearer (atom_coords: V3.t array) (query: V3.t) (threshold: float): bool =
  let t2 = threshold *. threshold in
  A.exists (fun xyz ->
      V3.dist2 query xyz < t2
    ) atom_coords

(* WARNING: for MM/GBSA, radii must probably be some Born radii, not vdW *)
let compute_SAS_dot_surface m =
  (* augment all radii by water radius *)
  let aug_radii = A.map ((+.) Const.r_H2O) m.r_a in
  (* reference dot surfaced sphere with N points BUT at origin *)
  let dot_surfaced_sphere =
    (* A.map V3.of_triplet (A.of_list (Genspir.genspir_cartesian 1000)) in *)
    (* reduced to 50 points per atom max; so that we can also run this for the protein *)
    A.map V3.of_triplet (A.of_list (Genspir.genspir_cartesian 50)) in
  let atom_coords = get_all_atom_coords m in
  let unfiltered_dots_per_atom =
    let surfaced_shperes_at_0 =
      A.map (fun r ->
          A.map (fun dot -> V3.mult dot r) dot_surfaced_sphere
        ) aug_radii in
    A.map2 (fun xyz surfaced_sphere_at_0 ->
        A.map (V3.add xyz) surfaced_sphere_at_0
      ) atom_coords surfaced_shperes_at_0 in
  (* filter them *)
  let index_center_radii =
    A.mapi (fun i aug_radius ->
        (i, atom_coords.(i), aug_radius)
      ) aug_radii in
  let dots_per_atom =
    A.mapi (fun i dots ->
        A.filter (fun dot ->
            not (atom_too_near i dot index_center_radii)
          ) dots
      ) unfiltered_dots_per_atom in
  (aug_radii, dots_per_atom)

(* retrieve some precomputed vdW param for
   an atom from the ligand and one from the protein *)
let get_lig_prot mtx lig_j prot_i =
  uget (uget mtx lig_j) prot_i

(* chemical formula of molecule *)
let brute_formula (_m: t): string =
  failwith "Mol.brute_formula: UNIMPLEMENTED"

(* just load anums/species for this ligand *)
let ene_QM_anums lig =
  QM_ene.QM_ene.__init__ ~anums:lig.elt_a ()

(* load anums for this ligand + lig_E_intra_QM *)
let ene_QM_prepare lig =
  let species = QM_ene.QM_ene.__init__ ~anums:lig.elt_a () in
  (* assume init conformer is a reasonably low one *)
  let init_ene_QM = QM_ene.QM_ene.score_conf species ~xs:lig.xs ~ys:lig.ys ~zs:lig.zs () in
  (species, init_ene_QM)

(* score ligand conformer for previously loaded ligand *)
let ene_QM_score species lig =
  QM_ene.QM_ene.score_conf species ~xs:lig.xs ~ys:lig.ys ~zs:lig.zs ()

(* UFF non-bonded energy term between a protein and a ligand
   brute force *)
let ene_inter_UFF_global_brute prot lig =
  let p = num_atoms prot in
  let l = num_atoms lig in
  let sum_elec = ref 0.0 in
  let sum_vdW  = ref 0.0 in
  (* for all protein atoms *)
  for i = 0 to p - 1 do
    let q_i       = get_charge prot i in
    let r_i       = get_xyz    prot i in
    let prot_anum = get_anum   prot i in
    (* for all ligand atoms *)
    for j = 0 to l - 1 do
      let q_j      = get_charge lig j in
      let r_j      = get_xyz    lig j in
      let lig_anum = get_anum   lig j in
      let r_ij = Math.non_zero_dist (V3.dist r_i r_j) in
      let vdw = UFF.vdW_xiDi prot_anum lig_anum in
      let p6 = FF.pow6 (vdw.x_ij /. r_ij) in
      sum_elec := !sum_elec +. ((q_i *. q_j) /. r_ij);
      sum_vdW := !sum_vdW +. (vdw.d_ij *. ((-2.0 *. p6) +. (p6 *. p6)))
    done;
  done;
  (UFF.elec_weight_protein *. !sum_elec) +. !sum_vdW

(* UFF non-bonded energy term between a protein and a ligand; brute force
   "shifted" means it is the "global" version above, but shifted *)
let ene_inter_UFF_shifted_brute prot lig =
  let p = num_atoms prot in
  let l = num_atoms lig  in
  let sum_elec = ref 0.0 in
  let sum_vdW  = ref 0.0 in
  (* for all protein atoms *)
  for i = 0 to p - 1 do
    let q_i       = get_charge prot i in
    let r_i       = get_xyz    prot i in
    let prot_anum = get_anum   prot i in
    (* for all ligand atoms *)
    for j = 0 to l - 1 do
      let r_ij2 = V3.dist2 r_i (get_xyz lig j) in
      (* over 12*12, it's all 0s *)
      if r_ij2 < 144.0 then
        begin
          let q_j      = get_charge lig j in
          let lig_anum = get_anum   lig j in
          let r_ij = Math.non_zero_dist (sqrt r_ij2) in
          let w = FF.shift_12A r_ij in
          let vdw = UFF.vdW_xiDi prot_anum lig_anum in
          let p6 = FF.pow6 (vdw.x_ij /. r_ij) in
          sum_elec := !sum_elec +. w *. ((q_i *. q_j) /. r_ij);
          sum_vdW := !sum_vdW +. w *. (vdw.d_ij *. ((-2.0 *. p6) +. (p6 *. p6)))
        end
    done;
  done;
  (UFF.elec_weight_protein *. !sum_elec) +. !sum_vdW

(* ligand intra-molecular energy; shifted, brute-force
   only non-bonded interactions are taken into account *)
let ene_intra_UFFNB_shifted_brute_DEPRECATED lig =
  let n = num_atoms lig  in
  let sum_elec = ref 0.0 in
  let sum_vdW  = ref 0.0 in
  for i = 0 to n - 2 do
    let q_i         = get_charge lig i in
    let r_i         = get_xyz    lig i in
    let anum_i      = get_anum   lig i in
    let interacting = get_interacting lig i in
    for j = i + 1 to n - 1 do
      (* 1,2 and 1,3 NB-interactions are ignored by UFF;
         starting from 3 bonds away (1,4 interactions),
         they are fully taken into account *)
      if A.unsafe_get interacting j then
        let r_ij = Math.non_zero_dist (V3.dist r_i (get_xyz lig j)) in
        if r_ij < 12.0 then
          let anum_j = get_anum lig j in
          let w = FF.shift_12A r_ij in
          let vdw = UFF.vdW_xiDi anum_i anum_j in
          let p6 = FF.pow6 (vdw.x_ij /. r_ij) in
          sum_elec := !sum_elec +. w *. ((q_i *. (get_charge lig j)) /. r_ij);
          sum_vdW  := !sum_vdW  +. w *. (vdw.d_ij *. ((-2.0 *. p6) +. (p6 *. p6)))
    done
  done;
  (UFF.elec_weight_protein *. !sum_elec) +. !sum_vdW

(* ligand intra-molecular energy; brute-force
   only non-bonded interactions are taken into account *)
let ene_intra_UFFNB_brute lig =
  let n = num_atoms lig  in
  let sum_elec = ref 0.0 in
  let sum_vdW  = ref 0.0 in
  for i = 0 to n - 2 do
    let q_i         = get_charge lig i in
    let r_i         = get_xyz    lig i in
    let anum_i      = get_anum   lig i in
    let interacting = get_interacting lig i in
    for j = i + 1 to n - 1 do
      (* 1,2 and 1,3 NB-interactions are ignored by UFF;
         starting from 3 bonds away (1,4 interactions),
         they are fully taken into account *)
      if A.unsafe_get interacting j then
        let r_ij = Math.non_zero_dist (V3.dist r_i (get_xyz lig j)) in
        let anum_j = get_anum lig j in
        let vdw = UFF.vdW_xiDi anum_i anum_j in
        let p6 = FF.pow6 (vdw.x_ij /. r_ij) in
        sum_elec := !sum_elec +. (q_i *. (get_charge lig j)) /. r_ij;
        sum_vdW  := !sum_vdW  +. vdw.d_ij *. ((-2.0 *. p6) +. (p6 *. p6))
    done
  done;
  (UFF.elec_weight_protein *. !sum_elec) +. !sum_vdW

exception Ene_intra (* to catch any python layer error *)

(* let rdkit do the full UFF evaluation of current conformer *)
let ene_intra_UFF_rdkit_exn lig conf =
  try
    Rdkit.update_coords conf ~xs:lig.xs ~ys:lig.ys ~zs:lig.zs ();
    Rdkit.score_UFF conf ()
  with _exn -> raise Ene_intra

(* let rdkit do the full MMFF evaluation of current conformer *)
let ene_intra_MMFF_rdkit_exn lig conf =
  try
    Rdkit.update_coords conf ~xs:lig.xs ~ys:lig.ys ~zs:lig.zs ();
    Rdkit.score_MMFF conf ()
  with _exn -> raise Ene_intra

type ele_vdW = { mutable sum_elec: float;
                 mutable sum_vdW: float }

(* UFF non-bonded energy term between a protein and a ligand
   shifted + bst-indexed protein.
   electrostatic and vdW parts are kept separated for the need
   of rescoring functions *)
let ene_inter_UFF_shifted_bst_components prot_bst prot lig =
  let l = num_atoms lig in
  let res = { sum_elec = 0.0; sum_vdW = 0.0 } in
  (* for all ligand atoms *)
  for j = 0 to l - 1 do
    let q_j      = get_charge lig j in
    let r_j      = get_xyz    lig j in
    let lig_anum = get_anum   lig j in
    (* get all "neighbor" protein atoms *)
    let query = Atom.create (-1) r_j in
    (* fixed 12A cutoff distance *)
    let indexes = BST.neighbors query 12.0 prot_bst in
    L.iter (fun prot_atom ->
        (* for all neighbor protein atoms *)
        let i = Atom.get_idx prot_atom in
        let q_i       = get_charge prot i in
        let r_i       = get_xyz    prot i in
        let prot_anum = get_anum   prot i in
        let r_ij = Math.non_zero_dist (V3.dist r_i r_j) in
        (* shift function weight *)
        let w = FF.shift_12A r_ij in
        let vdw = UFF.vdW_xiDi lig_anum prot_anum in
        let p6 = FF.pow6 (vdw.x_ij /. r_ij) in
        res.sum_elec <- res.sum_elec +. w *. ((q_i *. q_j) /. r_ij);
        res.sum_vdW <- res.sum_vdW +. w *. (vdw.d_ij *. ((-2.0 *. p6) +. (p6 *. p6)))
      ) indexes
  done;
  res.sum_elec <- res.sum_elec *. UFF.elec_weight_protein;
  res

let ene_inter_UFF_shifted_bst prot_bst prot lig =
  let res = ene_inter_UFF_shifted_bst_components prot_bst prot lig in
  res.sum_elec +. res.sum_vdW

(* specific version for single-atom ligands; all at same coordinate
   /!\ FOR GRID INIT ONLY /!\ *)
let ene_inter_UFF_shifted_grid prot_bst prot r_j ligs =
  let n = A.length ligs in
  let sum_elec = A.make n 0.0 in
  let sum_vdW  = A.make n 0.0 in
  (* all "neighbor" protein atoms at fixed cutoff distance *)
  let indexes = BST.neighbors (Atom.create (-1) r_j) 12.0 prot_bst in
  L.iter (fun prot_atom ->
      (* for each neighbor protein atoms *)
      let i = Atom.get_idx prot_atom in
      let q_i       = get_charge prot i in
      let r_i       = get_xyz    prot i in
      let prot_anum = get_anum   prot i in
      let r_ij = Math.non_zero_dist (V3.dist r_i r_j) in
      (* shift function weight *)
      let w = FF.shift_12A r_ij in
      (* for each single_atom ligands *)
      A.iteri (fun l lig ->
          let vdw = UFF.vdW_xiDi (get_anum lig 0) prot_anum in
          let p6 = FF.pow6 (vdw.x_ij /. r_ij) in
          uset sum_elec l ((uget sum_elec l) +. w *. ((q_i *. (get_charge lig 0)) /. r_ij));
          uset sum_vdW  l ((uget sum_vdW  l) +. w *. (vdw.d_ij *. ((-2.0 *. p6) +. (p6 *. p6))))
        ) ligs
    ) indexes;
  A.map2 (fun sum_elec sum_vdW ->
      UFF.elec_weight_protein *. sum_elec +. sum_vdW
    ) sum_elec sum_vdW

(* maximum factorization of the work for grid initialization;
   For each (and only once per) exact grid coordinate:
   get not too far protein atom indexes,
   distance of each to r_j and associated weight *)
let indexes_rij_w prot_bst prot r_j =
  (* all "neighbor" protein atoms at fixed cutoff distance *)
  let indexes = BST.neighbors (Atom.create (-1) r_j) 12.0 prot_bst in
  let n = L.length indexes in
  let res = A.make n (0, 0., 0.) in
  (* for each neighbor protein atom *)
  L.iteri (fun res_i prot_atom ->
      let prot_i = Atom.get_idx prot_atom in
      let r_ij = Math.non_zero_dist (V3.dist (get_xyz prot prot_i) r_j) in
      (* shift function weight *)
      let w = FF.shift_12A r_ij in
      res.(res_i) <- (prot_i, r_ij, w)
    ) indexes;
  res

(* UFF non-bonded energy term between a protein and a ligand
   interpolated from grid points *)
let ene_inter_UFF_interp grid ff_comps lig =
  (* for each ligand atom, get the contribution from the corresp. FF comp. *)
  let l = num_atoms lig in
  let res = ref 0.0 in
  (* for all ligand atoms *)
  for j = 0 to l - 1 do
    res := !res +. G3D.trilin grid (uget ff_comps (get_typ lig j)) (get_xyz lig j)
  done;
  !res

exception Atom_clash

(* do not use: was inlined
let atom_clash_DEPRECATED xyz_i anum_i xyz_j anum_j =
   V3.dist2 xyz_i xyz_j < Ptable.vdW_clash_params2.(anum_i).(anum_j)
*)

(* a clash is two atoms which are at least three bonds away
   on the molecular graph and whose 3D cartesian distance is less than
   a * (vdW_i + vdW_j) *)
let ligand_atoms_clash verbose m =
  let log =
    if verbose then
      (fun i j -> Log.warn "%s: vdW clash: %s:%d Vs %s%d (d_topo=%d)"
          (get_name m)
          (get_symbol m i) i
          (get_symbol m j) j
          (get_topo_dist m i j))
    else
      (fun _i _j -> ()) in
  try
    let n = num_atoms m in
    for i = 0 to n - 2 do
      let xyz_i  = get_xyz  m i in
      let anum_i = get_anum m i in
      for j = i + 1 to n - 1 do
        (* this list of atom pairs with get_topo_dist >= 3
           could be calculated once and for all, if we have
           to analyze many conformers of the same ligand *)
        if get_topo_dist m i j >= 3 &&
           V3.dist2 xyz_i (get_xyz m j) <
           Ptable.vdW_clash_params2.(anum_i).(get_anum m j) then
          let () = log i j in
          raise Atom_clash
      done
    done;
    false
  with Atom_clash -> true

(* array of atom pair indexes with [get_topo_dist m i j >= 3];
   each array element is of the form:
   (i, [j0, j1, ..., jn]) *)
let get_topo_dists_ge3 m: (int * int array) array =
  let n = num_atoms m in
  A.init (n - 1) (fun i ->
      (* i in [0 .. n-2] *)
      let js = ref [] in
      for j = i + 1 to n - 1 do
        if get_topo_dist m i j >= 3 then
          js := j :: !js
      done;
      (i, A.of_list (L.rev !js))
    )

(* high performance version of [ligand_atoms_clash] *)
let ligand_atoms_clash_fast topo_dists_ge3 m =
  try
    A.iter (fun (i, js) ->
        let xyz_i = get_xyz m i in
        let clash_params2 = Ptable.vdW_clash_params2.(get_anum m i) in
        A.iter (fun j ->
            if V3.dist2 xyz_i (get_xyz m j) < clash_params2.(get_anum m j) then
              raise Atom_clash
          ) js
      ) topo_dists_ge3;
    false
  with Atom_clash -> true

(* for conformer generation via dihedral angles only *)
let rotate_bonds centered_lig rbonds =
  let lig = copy centered_lig in
  let n = A.length rbonds in
  assert(n = num_rbonds centered_lig);
  A.iteri (fun i x ->
      rotate_bond lig i x
    ) rbonds;
  lig

let decode_rbonds_conf
    (rbonds: int) (steps_per_rot_bond: int) (rot_step: float) (conf_id: int)
  : float array =
  let powers = (* from high to low *)
    A.init rbonds (fun i -> BatInt.pow steps_per_rot_bond (rbonds - (i + 1))) in
  let rest = ref conf_id in
  A.map (fun p ->
      let integ = !rest / p in
      rest := !rest - (integ * p);
      (float integ) *. rot_step
    ) powers

(* upper bound on the number of conformers if using dihedral
   rotational step [alpha] in radians; there might be less, if
   vdW clashing ones are removed *)
let count_conformers m alpha =
  let rbonds = num_rbonds m in
  (* int_of_float <=> floor, because we will never go over 2*pi *)
  let steps_per_rot_bond = int_of_float (Math.two_pi /. alpha) in
  BatInt.pow steps_per_rot_bond rbonds

(* generate all conformers, using [alpha] (in radians) rotational step;
   vdW clashing conformers are discarded;
   if |conformers| > [cap]; only a random partition of size [cap]
   will be returned.
   This can trigger OOM if greedy algorithm. *)
let sample_conformers rng cap m alpha: t array =
  let name = get_name m in
  let topo_dists_ge3 = get_topo_dists_ge3 m in
  (if ligand_atoms_clash_fast topo_dists_ge3 m then
     Log.warn "%s: init conf vdW clash" name
  );
  let steps_per_rot_bond = int_of_float (Math.two_pi /. alpha) in
  let max_num_confs = count_conformers m alpha in
  (if max_num_confs > 10_000 then
     Log.warn "Mol.sample_conformers: %s: %d confs" name max_num_confs
  );
  (* each configuration of dihedrals is encoded as a number
     in base [per_rot_bond] *)
  let confs_to_generate =
    if max_num_confs <= cap then
      A.init max_num_confs (fun i -> i) (* all conf_ids *)
    else
      (* random partition of size [cap] *)
      let to_shuffle = A.of_list (L.range 1 `To (max_num_confs - 1)) in
      A.shuffle ~state:rng to_shuffle;
      A.sub (* always keep conf_id=0 (the initial conformer) *)
        (A.append [|0|] to_shuffle) 0 cap in
  (* non vdW clashing confs are Some, clashing ones are Nones *)
  let rbonds = num_rbonds m in
  let maybes =
    A.map (fun conf_id ->
        let rbonds = decode_rbonds_conf rbonds steps_per_rot_bond alpha conf_id in
        let maybe_clashing = rotate_bonds m rbonds in
        if ligand_atoms_clash_fast topo_dists_ge3 maybe_clashing then
          None
        else
          Some maybe_clashing
      ) confs_to_generate in
  (* return the array of remaining conformers *)
  A.map Option.get (A.filter Option.is_some maybes)

let ene_intra_no_vdW_clash topo_dists_ge3 lig =
  if ligand_atoms_clash_fast topo_dists_ge3 lig then
    Params.max_E
  else
    0.0

(* check if any heavy atom from the ligand is too near from
 * any protein heavy atom *)
let heavy_atom_clash prot_bst prot lig =
  try
    A.iteri (fun i anum ->
        (* only check for heavy atoms *)
        (if anum > 1 then
           let lig_xyz  = get_xyz  lig i in
           let lig_anum = get_anum lig i in
           let query = Atom.create (-1) lig_xyz in
           let neighbors = BST.neighbors query (2.0 *. Ptable.vdW_max) prot_bst in
           L.iter (fun prot_atom ->
               (* for all neighbor protein atoms *)
               let j = Atom.get_idx prot_atom in
               if get_anum prot j > 1 &&
                  V3.dist2 lig_xyz (get_xyz prot j) <
                  Ptable.vdW_clash_params2.(lig_anum).(get_anum prot j) then
                 raise Atom_clash
             ) neighbors
        )
      ) lig.elt_a;
    false
  with Atom_clash -> true

(* is any of [lig]'s atoms vdW clashing with the protein?
   WARNING: this is not an exact test; some clashes are not
   detected; but if this returns true, you can be sure there is a clash *)
let protein_ligand_clash prot_grid prot_vdW_mask lig =
  try
    let n = num_atoms lig in
    for i = 0 to n - 1 do
      if G3D.vdW_clash_OR prot_grid prot_vdW_mask (get_xyz lig i) then
        raise Atom_clash
    done;
    false
  with Atom_clash -> true

(* is [lig]'s geometric center vdW occuppied by one of its atoms;
   this would constrain where we can translate (then try all rotations)
   this ligand in a protein's ligand-binding site
   (only to non protein vdW occuppied grid cells) *)
let is_ligand_center_vdW_occuppied lig =
  try
    let c = get_center lig in
    let n = num_atoms  lig in
    for i = 0 to n - 1 do
      if V3.dist c (get_xyz lig i) < get_radius lig i then
        raise Atom_clash
    done;
    false
  with Atom_clash -> true

let heavy_atom_indexes m =
  let res = ref [] in
  let n = num_atoms m in
  for i = 0 to n - 1 do
    if get_anum m i > 1 then (* HA *)
      res := i :: !res
  done;
  A.of_list (L.rev !res)

exception Found

(* there is at least one protein atom within [dist] from a ligand heavy atom *)
let lig_nearby_prot prot_bst dist lig =
  let n = num_atoms lig in
  try
    for i = 0 to n - 1 do
      if get_anum lig i > 1 then
        (* Heavy Atom *)
        let xyz_i = get_xyz lig i in
        let query = Atom.create (-1) xyz_i in
        let _nearest, nearest_dist = BST.nearest_neighbor query prot_bst in
        (* Log.info "Mol.lig_nearby_prot: nearest: %.2f" dist; *)
        if nearest_dist < dist then
          raise Found
    done;
    false
  with Found -> true

(* upper bound of ligand strain upon docking (in kcal/mol);
   Jain, JMC 2023 - https://doi.org/10.1021/acs.jmedchem.2c01744 *)
let ligand_strain_UB mol =
  max 0.0 (0.3 *. float ((num_heavy_atoms mol) - 10))
