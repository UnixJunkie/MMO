(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
 *
 * UFF, a Full Periodic Table Force Field for Molecular Mechanics and
 * Molecular Dynamics Simulations *)

module IMap = BatMap.Int

let anum_xi_Di = (* currently supported elements *)
  [|( 0, (0.0  , 0.0  )); (* virtual chemical element w/o vdW interactions *)
    ( 1, (2.886, 0.044));
    ( 6, (3.851, 0.105));
    ( 7, (3.660, 0.069));
    ( 8, (3.500, 0.060));
    ( 9, (3.364, 0.050));
    (12, (3.021, 0.111));
    (15, (4.147, 0.305));
    (16, (4.035, 0.274));
    (17, (3.947, 0.227));
    (35, (4.189, 0.251));
    (53, (4.500, 0.339))|]

(* 332.0637: constant from the UFF paper so that ene is in kcal/mol *)
let elec_weight_protein = 332.0637 /. Const.epsilon_prot

(* len[0..118] = 119
   there is no chemical element w/ anum=0; but it is useful
   to not have to shift all atomic numbers by 1; plus we
   use it so store a special virtual chemical elements
   which has no vdW interactions *)
let anums = 119

type xij_dij = { x_ij: float ;
                 d_ij: float }

let xidi_params = A.make (anums * anums) { x_ij = nan; d_ij = nan }

let vdW_xiDi anu1 anu2 =
  A.get xidi_params ((anu1 * anums) + anu2)

(* initialize [xidi_params] *)
let () =
  A.iter (fun     (anu1, (xi1, di1)) ->
      A.iter (fun (anu2, (xi2, di2)) ->
          let i = (anu1 * anums) + anu2 in
          let x_ij = FF.geo_mean xi1 xi2 in
          let d_ij = FF.geo_mean di1 di2 in
          xidi_params.(i) <- { x_ij; d_ij }
        ) anum_xi_Di
    ) anum_xi_Di
