(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

(* compute the molecular graph plus some operations
 * - find rings
 * - compute topological distance matrix
 * - find RotBonds *)

module G = Graph.Pack.Graph
module PathCheck = G.PathCheck
module IntSet = BatSet.Int

open Printf

(* Integer-Weighted Edge *)
module IWE = struct
  type edge = G.E.t
  type t = int
  let weight e = G.E.label e
  let zero = 0
  let add = (+)
  let sub = (-)
  let compare = compare
end

(* implements all_pairs_shortest_paths *)
module Johnson = Graph.Path.Johnson(G)(IWE)

type t = G.t

let create (mol: Mol2.t): (G.V.t array * t) =
  let num_atoms = Mol2.num_atoms mol in
  (* we need to keep vertices because ocamlgraph is inefficient at
     retrieving them... *)
  let vertices = A.init num_atoms G.V.create in
  let res = G.create ~size:num_atoms () in
  A.iter (G.add_vertex res) vertices;
  A.iter (fun bond ->
      let edge =
        G.E.create
          vertices.(Mol2.get_bond_src bond)
          1 (* at distance one bond *)
          vertices.(Mol2.get_bond_dst bond) in
      G.add_edge_e res edge
    ) (Mol2.get_bonds mol);
  (vertices, res)

(* an atom not connected to the rest of the molecule in a .mol2 file
   will crash dist_matrix below *)
exception Disconnected_atom

let dist_matrix vertices mol_graph =
  let n = G.nb_vertex mol_graph in
  let res = A.make (n * n) 0 in
  let all_dists = Johnson.all_pairs_shortest_paths mol_graph in
  (* the diagonal is already all 0s so we don't update it *)
  for i = 0 to n - 2 do
    let v1 = vertices.(i) in
    for j = i + 1 to n - 1 do
      let v2 = vertices.(j) in
      let d =
        try Johnson.HVV.find all_dists (v1, v2)
        with Not_found -> raise Disconnected_atom in
      (* matrix is diagonal symetric *)
      res.(i + (j * n)) <- d;
      res.(j + (i * n)) <- d
    done
  done;
  res

let print_dist_matrix out mol_graph mtx =
  let n = G.nb_vertex mol_graph in
  for i = 0 to n - 1 do
    for j = 0 to n - 1 do
      let d = mtx.(i + (j * n)) in
      fprintf out (if j = 0 then "%d" else " %d") d
    done;
    fprintf out "\n"
  done

let bprint_dist_matrix buff mol_graph mtx =
  let n = G.nb_vertex mol_graph in
  for i = 0 to n - 1 do
    for j = 0 to n - 1 do
      let d = mtx.(i + (j * n)) in
      bprintf buff (if j = 0 then "%d" else " %d") d
    done;
    Buffer.add_char buff '\n'
  done

(* compute the degree of each node/vertex/atom *)
let degrees g =
  let n = G.nb_vertex g in
  let res = A.make n 0 in
  G.iter_edges (fun v1 v2 ->
      let i = G.V.label v1 in
      let j = G.V.label v2 in
      res.(i) <- res.(i) + 1;
      res.(j) <- res.(j) + 1
    ) g;
  res

let print_degrees degs =
  A.iteri (fun i d ->
      printf (if i = 0 then "%d" else " %d") d
    ) degs;
  printf "\n"

let is_ring_bond vertices g b: bool =
  let v1 = vertices.(Mol2.get_bond_src b) in
  let v2 = vertices.(Mol2.get_bond_dst b) in
  let edge = G.find_edge g v1 v2 in
  G.remove_edge_e g edge;
  let res = PathCheck.check_path (PathCheck.create g) v1 v2 in
  (* restore the original molecular graph *)
  G.add_edge_e g edge;
  res

let smallest_set s1 s2 =
  if IntSet.cardinal s1 <= IntSet.cardinal s2 then
    s1
  else
    s2

(* a rotatable bond is a single bond, not in any ring, between two
 * atoms where none of them is a terminal atom (terminal means degree=1) *)
let list_rotatable_bonds mol vertices degrees g: Mol2.bond list =
  let res = ref [] in
  A.iter (fun b ->
      (* tests in increasing order of complexity *)
      if Mol2.(b.bo) = 1.0 &&
         degrees.(Mol2.(b.src)) > 1 &&
         degrees.(Mol2.(b.dst)) > 1 &&
         not (is_ring_bond vertices g b) then
        res := b :: !res
    ) (Mol2.get_bonds mol);
  L.rev !res

let smallest_group groups =
  let ht = Ht.create 2 in
  A.iter (fun x ->
      let prev_count = Ht.find_default ht x 0 in
      Ht.replace ht x (prev_count + 1)
    ) groups;
  let key_values = Ht.bindings ht in
  match key_values with
  | [(g1, g1_card); (g2, g2_card)] ->
    if g1_card <= g2_card then
      g1
    else
      g2
  | _ -> failwith "Mol_graph.smallest_group: not two groups"

(* compute the set of connected atoms which are impacted by a rotation
 * of RotBond [b]; !!! [b] MUST BE A ROTATABLE BOND !!! *)
let compute_rot_group vertices g b: IntSet.t =
  (* rm RotBond *)
  let v1 = vertices.(Mol2.get_bond_src b) in
  let v2 = vertices.(Mol2.get_bond_dst b) in
  let edge = G.find_edge g v1 v2 in
  G.remove_edge_e g edge;
  (* initially, each atom in a different group *)
  let n = G.nb_vertex g in
  let groups = A.init n (fun i -> i) in
  (* loop until groups are stable *)
  let changed = ref true in
  while !changed do
    changed := false;
    G.iter_edges (fun v1 v2 ->
        let i = G.V.label v1 in
        let j = G.V.label v2 in
        let g_i = groups.(i) in
        let g_j = groups.(j) in
        if g_i < g_j then
          (groups.(j) <- g_i;
           changed := true)
        else if g_j < g_i then
          (groups.(i) <- g_j;
           changed := true)
      ) g
  done;
  (* restore initial molecular graph *)
  G.add_edge_e g edge;
  (* movable = smallest_group *)
  let smallest_group_id = smallest_group groups in
  A.fold_lefti (fun acc i x ->
      if x = smallest_group_id then
        IntSet.add i acc
      else
        acc
    ) IntSet.empty groups

let print_rot_group out g b rot_group =
  let i = Mol2.get_bond_src b in
  let j = Mol2.get_bond_dst b in
  fprintf out "%d-%d=" i j;
  let n = G.nb_vertex g in
  for i = 0 to n - 1 do
    if IntSet.mem i rot_group then
      fprintf out " 1"
    else
      fprintf out " 0"
  done;
  fprintf out "\n"

let bprint_rot_group buff g b rot_group =
  let i = Mol2.get_bond_src b in
  let j = Mol2.get_bond_dst b in
  bprintf buff "%d-%d=" i j;
  let n = G.nb_vertex g in
  for i = 0 to n - 1 do
    if IntSet.mem i rot_group then
      bprintf buff " 1"
    else
      bprintf buff " 0"
  done;
  bprintf buff "\n"
