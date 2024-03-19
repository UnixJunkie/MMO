
(* list the names of ligands which only have non-conformer-changing RotBonds *)

let main () =
  let pqrs_in_fn = Sys.argv.(1) in
  let ligands = Mol.ligands_of_pqrs_file pqrs_in_fn in
  L.iter (fun lig ->
      if Mol.num_rbonds lig = A.length (Mol.non_conformer_changing_rbonds lig) then
        Printf.printf "%s\n" (Mol.get_name lig)
    ) ligands

let () = main ()
