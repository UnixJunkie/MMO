(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

(* MOL2 parsing *)

open Printf

let molecule_tag = "@<TRIPOS>MOLECULE"
let atoms_tag = "@<TRIPOS>ATOM"
(* sometimes, there is a block between atoms and bonds;
 * this block describes atom attributes *)
let united_atoms_tag = "@<TRIPOS>UNITY_ATOM_ATTR"
let bonds_tag = "@<TRIPOS>BOND"
(* end_bonds_tag might be CCDC GOLD output specific *)
let end_bonds_tag = "@<TRIPOS>SUBSTRUCTURE"

type gold_mol2_reader_state = Waiting_header
                            | Reading_molecule

(* keep only some mol2 "sections" from a .mol2 file created by CCDC GOLD
   of course they use a non-standard file format... *)
let filter_gold_output_mol2 input output =
  let state = ref Waiting_header in
  try
    while true do
      let line = input_line input in
      match !state with
      | Waiting_header ->
        if line = molecule_tag then
          (state := Reading_molecule;
           fprintf output "%s\n" line)
      | Reading_molecule ->
        if line = end_bonds_tag then
          (* no need for this tag in the output file *)
          state := Waiting_header
        else
          fprintf output "%s\n" line
    done
  with End_of_file -> ()

type atom = { x: float;
              y: float;
              z: float;
              q: float;
              r: float;
              symb: string; (* chemical elt. symbol *)
              (* we don't record the id (an index starting from 1) *)
              name: string;
              typ: string }

let is_lone_pair a =
  a.typ = "LP"

let get_atom_coord a =
  V3.make a.x a.y a.z

(* update atom's coordinates *)
let atom_move a x' y' z' =
  { a with x = x'; y = y'; z = z' }

let atom_charge_is_null atom =
  atom.q = 0.0

let pqrs_line_of_atom a =
  (* pqrs files ignore: atom_id, atom_name and atom_type *)
  sprintf "%g %g %g %g %g %s" a.x a.y a.z a.q a.r a.symb

let xyz_line_of_atom a =
  sprintf "%s %g %g %g" a.symb a.x a.y a.z

type bond = { src: int;
              dst: int;
              bo: float;
              typ: string }

let bond_create src dst bo typ =
  { src; dst; bo; typ }

let get_bond_src b =
  b.src

let get_bond_dst b =
  b.dst

let get_bond_order b =
  b.bo

let dummy_bond =
  bond_create (-1) (-1) nan

(* one molecule read from a .mol2 file *)
type t = { name: string;
           atoms: atom array;
           bonds: bond array }

let remove_lone_pairs m =
  let lone_pairs = A.map is_lone_pair m.atoms in
  let num_LPs = A.count_matching (fun x -> x) lone_pairs in
  (if num_LPs > 0 then
     Log.warn "Mol2.remove_lone_pairs: %d LPs in %s" num_LPs m.name
  );
  (* filter out bonds connecting a LP *)
  let bonds' =
    A.filter (fun bond ->
        not (lone_pairs.(bond.src) || lone_pairs.(bond.dst))
      ) m.bonds in
  (* filter out LP "atoms" *)
  let atoms' =
    A.filter (fun a -> not (is_lone_pair a)) m.atoms in
  { m with atoms = atoms'; bonds = bonds' }

let move m atoms' =
  { m with atoms = atoms' }

let get_name m =
  m.name

let num_atoms m =
  A.length m.atoms

let num_bonds m =
  A.length m.bonds

let get_atoms m =
  m.atoms

(* return all atom coordinates as an array *)
let get_coords m =
  A.map get_atom_coord (get_atoms m)

let get_bonds m =
  m.bonds

let create name atoms bonds =
  { name; atoms; bonds }

(* bond order from mol2 type *)
let bond_order = function
  | "1" -> 1.
  | "2" -> 2.
  | "3" -> 3.
  | "ar" -> 1.5
  | "am" -> 1. (* same as rdkit *)
  | "du" -> failwith "Mol2.bond_order: dummy bond"
  | "un" -> failwith "Mol2.bond_order: unknown bond"
  | "nc" -> failwith "Mol2.bond_order: not connected bond"
  | other -> failwith ("Mol2.bond_order: " ^ other)

(* mol2 specification:
   'atom_id atom_name x y z atom_type \
    [subst_id [subst_name [charge [status_bit]]]]' *)
let parse_atom_line line =
  try
    Scanf.sscanf line " %d %s %f %f %f %s %d %s %f"
      (fun _id name x y z typ _resnum _resname q ->
         try
           let anum = Ptable.anum_of_mol2_type typ in
           let r = Ptable.vdW_radius_of_anum anum in
           let symb = Ptable.symbol_of_anum anum in
           { x; y; z; q; r; symb; name; typ }
         with (Ptable.Unsupported_atom other) as exn ->
           (* CCDC GOLD introduces lone pairs in its mol2 output files
            * we filter them out later and correct the number of atoms and bonds
            * at that time *)
           if other = "LP" then
             let r = nan in
             let symb = "LP" in
             { x; y; z; q; r; symb; name; typ }
           else
             raise exn
      )
  with exn ->
    (Log.error "Mol2.parse_atom_line: cannot parse: %s" line;
     raise exn)

let write_atom_line (out: out_channel) (id: int) (a: atom): unit =
  fprintf out "%7d %-8s%10.4f%10.4f%10.4f %-8s  1 <0>     %10.4f\n"
    id a.name a.x a.y a.z a.typ a.q

let write_bond_line (out: out_channel) (id: int) (b: bond): unit =
  (* bond indexes start from 1 *)
  fprintf out "%6d%5d%5d %s\n" id (b.src + 1) (b.dst + 1) b.typ

let parse_bond_line line =
  try Scanf.sscanf line " %d %d %d %s"
        (fun _id src' dst' typ ->
           (* mol2 indexes start at 1; computer science ones at 0 *)
           let src = src' - 1 in
           let dst = dst' - 1 in
           let bo = bond_order typ in
           { src; dst; bo; typ }
        )
  with exn ->
    (Log.error "Mol2.parse_bond_line: cannot parse: %s" line;
     raise exn)

let parse_atoms_bonds_line line =
  try Scanf.sscanf line " %d %d %s"
        (fun num_atoms num_bonds _whatever -> (num_atoms, num_bonds))
  with exn ->
    (Log.fatal "Mol2.parse_atoms_bonds_line: cannot parse: %s" line;
     raise exn)

let write_atoms_bonds_line out m =
  fprintf out "%5d%6d%6d%6d%6d\n" (num_atoms m) (num_bonds m) 0 0 0

(* a single molecule in MOL2 format *)
let get_mol2_block_exn input: string array =
  let res = ref [] in
  begin
    try
      let line = ref (input_line input) in
      (if !line = molecule_tag then
         (* first time this file is being read *)
         res := !line :: !res
       else
         res := !line :: molecule_tag :: !res
      );
      line := input_line input;
      while !line <> molecule_tag do
        res := !line :: !res;
        line := input_line input
      done
    with End_of_file -> ()
  end;
  let arr = A.of_list (L.rev !res) in
  if A.length arr = 0 then
    raise End_of_file
  else
    arr

let read_all_mol2_blocks (fn: string): (string array) list =
  LO.with_in_file fn (fun input ->
      let res, exn = L.unfold_exn (fun () -> get_mol2_block_exn input) in
      if exn <> End_of_file then raise exn;
      res
    )

(* parse one molecule from an already opened input file *)
let parse_one_mol2_block mol2_block: t option =
  (* Log.error "###############################################"; *)
  (* A.iter (Log.error "%s") mol2_block; *)
  (* Log.error "///////////////////////////////////////////////"; *)
  if mol2_block.(0) <> molecule_tag then
    (Log.error "Mol2.parse_one_mol2_block: doesn't start w/ molecule_tag";
     None)
  else
    try
      (* next line is molecule's name *)
      let name = S.strip mol2_block.(1) in
      (* Log.warn "name: %s |mol2_block|=%d" name (A.length mol2_block); *)
      let atoms_bonds_abc = mol2_block.(2) in
      let n_atoms, n_bonds = parse_atoms_bonds_line atoms_bonds_abc in
      (* Log.warn "%d %d" n_atoms n_bonds; *)
      let i = ref 3 in
      (* look for atoms start *)
      while mol2_block.(!i) <> atoms_tag do
        incr i
      done;
      (* skip atoms_tag line *)
      incr i;
      let atoms =
        try
          A.init n_atoms (fun _ ->
              let res = parse_atom_line mol2_block.(!i) in
              incr i;
              res)
        with _exn -> [||] in
      let line = ref mol2_block.(!i) in
      (if !line = united_atoms_tag then
         begin (* skip united atoms lines *)
           Log.warn "Mol2.parse_one_mol2_block: %s in '%s'" united_atoms_tag name;
           while mol2_block.(!i) <> bonds_tag do
             incr i
           done;
           line := mol2_block.(!i);
         end
      );
      if !line <> bonds_tag then
        (Log.error "Mol2.parse_one_mol2_block: skipping no bonds_tag molecule '%s': %s"
           name !line;
         None)
      else
        let bonds =
          try
            incr i; (* skip bonds header *)
            A.init n_bonds (fun _ ->
                let res = parse_bond_line mol2_block.(!i) in
                incr i;
                res)
          with _exn -> [||] in
        let res = create name atoms bonds in
        if num_atoms res = n_atoms && num_bonds res = n_bonds then
          Some (remove_lone_pairs res)
        else
          (Log.error
             "Mol2.parse_one_mol2_block: skipping erroneous molecule '%s'"
             name;
           None)
    with Invalid_argument _maybe_index_out_of_bounds ->
      (Log.error "Mol2.parse_one_mol2_block: skipping erroneous molecule";
       None)

let read_one input: t option =
  parse_one_mol2_block (get_mol2_block_exn input)

let read_all fn: t array =
  LO.with_in_file fn (fun input ->
      let maybes, exn = L.unfold_exn (fun () -> read_one input) in
      if exn <> End_of_file then raise exn;
      let expected = L.length maybes in
      let certains = L.count_matching Option.is_some maybes in
      if certains = expected then
        A.of_list (L.map Option.get maybes)
      else
        (Log.fatal
           "Mol2.read_all: expected %d but only %d readable molecules in %s"
           expected certains fn;
         exit 1)
    )

(* output one mol2 block *)
let output_named_one out name mol =
  fprintf out "%s\n" molecule_tag;
  fprintf out "%s\n" name;
  write_atoms_bonds_line out mol;
  (* standard header lines before atoms block *)
  fprintf out "SMALL\nUSER_CHARGES\n\n";
  fprintf out "%s\n" atoms_tag;
  A.iteri (fun i a ->
      (* atom indexes start from 1 *)
      write_atom_line out (i + 1) a
    ) mol.atoms;
  fprintf out "%s\n" bonds_tag;
  A.iteri (fun i b ->
      (* bond indexes start from 1 *)
      write_bond_line out (i + 1) b
    ) mol.bonds

let output_one out mol =
  output_named_one out mol.name mol

let output_all fn mols =
  LO.with_out_file fn (fun out ->
      A.iter (output_one out) mols
    )

let write_one_to_file fn mol =
  LO.with_out_file fn (fun out -> output_one out mol)

(* in case there are several molecules in that file, all
 * but the first one are ignored... *)
let read_one_from_file fn =
  LO.with_in_file fn (fun input ->
      match read_one input with
      | Some x -> x
      | None -> (Log.fatal "Mol2.read_one_from_file: could not read %s" fn;
                 exit 1))

(* separate each molecule to its own mol2 file;
   return the array of such filenames *)
let read_all_raw out_dir mol2_fn =
  let all_mols =
    LO.with_in_file mol2_fn (fun input ->
        let res, exn = L.unfold_exn (fun () -> get_mol2_block_exn input) in
        if exn <> End_of_file then raise exn;
        A.of_list res
      ) in
  let n = A.length all_mols in
  Log.info "read %d molecules from %s" n mol2_fn;
  A.mapi (fun i mol2_lines ->
      let fn = sprintf "%s/%d.mol2" out_dir i in
      Utls.lines_to_file_a fn mol2_lines;
      fn
    ) all_mols

(* extract all molecule names from given .mol2 file *)
let mol_names_from_file fn =
  let res = ref SS.empty in
  LO.with_in_file fn (fun input ->
      try
        while true do
          match read_one input with
          | None -> ()
          | Some mol ->
            let name = get_name mol in
            res := SS.add name !res
        done
      with End_of_file -> ()
    );
  !res
