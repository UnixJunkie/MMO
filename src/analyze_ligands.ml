
(* count how many ligands have their geometric center HA occupied;
   to determine if we will exploit this later to
   skip all rotation searches at a given grid point while
   exhaustive docking *)

let pqrs_in_fn = Sys.argv.(1) in
let ligands = Mol.ligands_of_pqrs_file pqrs_in_fn in
let occuppied = L.filter Mol.is_ligand_center_vdW_occuppied ligands in
let occ = L.length occuppied in
let tot = L.length ligands in
Printf.printf "%d/%d\n" occ tot
