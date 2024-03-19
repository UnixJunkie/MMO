(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

(* extract all molecule names from given .pqrs file *)
let mol_names_from_file fn =
  LO.fold fn (fun acc line ->
      if S.contains line ':' then
        let name =
          Scanf.sscanf line "%d:%d:%s"
            (fun _num_atoms _num_bonds name -> name) in
        SS.add name acc
      else
        acc
    ) SS.empty

let parse_header line: int * int * string =
  try (* ligand header w/ rot-bonds *)
    Scanf.sscanf line "%d:%d:%s" (fun num_atoms num_bonds name ->
        (num_atoms, num_bonds, name)
      )
  with _exn ->
    begin
      try (* receptor header (rot-bonds are not mentioned) *)
        Scanf.sscanf line "%d:%s" (fun n name -> (n, -1, name))
      with exn ->
        let () = Log.error "Mol.parse_pqrs_header: cannot parse: %s" line in
        raise exn
    end

let parse_body line: float * float * float * float * float * string =
  try Scanf.sscanf line "%f %f %f %f %f %s"
        (fun x y z q r elt -> (x, y, z, q, r, elt))
  with exn ->
    let () =
      Log.error "Molecule.parse_pqrs_body: cannot parse: %s" line in
    raise exn

(* parse one rot-bond line
 * e.g. '5-9= 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0' *)
let parse_rot_bond num_atoms input =
  let line = input_line input in
  try Scanf.sscanf line "%d-%d= %s@\n" (fun beg_atom end_atom rot_group ->
      let bool_strings = S.split_on_char ' ' rot_group in
      let rot_group = A.of_list (L.map Utls.bool_of_string bool_strings) in
      assert (A.length rot_group = num_atoms);
      (* convention: the rotation axis is from fixed atom to movable one;
       * to be precise, BOTH ATOMS ON THE AXIS NEED NOT ROTATE, only atoms after
       * the axis need *)
      assert (rot_group.(beg_atom) <> rot_group.(end_atom));
      let axis =
        if rot_group.(end_atom) then
          Bond.create beg_atom end_atom
        else
          Bond.create end_atom beg_atom in
      (axis, rot_group)
    )
  with exn ->
    let () = Log.fatal "Mol.parse_rot_bond: cannot parse: %s" line in
    raise exn

let parse_dist_matrix num_atoms input =
  let res = A.make (num_atoms * num_atoms) 0 in
  for i = 0 to num_atoms - 1 do
    let line = input_line input in
    let dists = S.split_on_char ' ' line in
    assert(L.length dists = num_atoms);
    L.iteri (fun j dist_str ->
        let dist = int_of_string dist_str in
        let k = i + (j * num_atoms) in
        res.(k) <- dist;
        (* check diagonal *)
        if i = j then
          assert(dist = 0)
      ) dists
  done;
  res

(* extract indexes of atoms in this rot. group *)
let rgroup_of_bools bond bools =
  let res = ref [] in
  A.iteri (fun i in_group ->
      (* exclude axis end from the group: it will never rotate *)
      if in_group && Bond.right_end bond <> i then
        res := i :: !res
    ) bools;
  A.of_list (L.rev !res)
