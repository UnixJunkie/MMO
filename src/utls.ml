(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

open Printf

let time_it (f: 'a -> 'b) (x: 'a): float * 'b =
  let start = Unix.gettimeofday () in
  let res = f x in
  (Unix.gettimeofday() -. start, res)

let bitmask_to_file fn mask =
  LO.with_out_file fn (fun out ->
      fprintf out "%s\n" (Bitv.M.to_string mask)
    )

let bitmask_from_file fn =
  LO.with_in_file fn (fun input ->
      Bitv.M.of_string (input_line input)
    )

let zstd_compress_file fn =
  let cmd = sprintf "zstd --rm -qf %s" fn in (* remove quiet force *)
  let dt, ret = time_it Unix.system cmd in
  (if ret <> Unix.WEXITED 0 then
     let () = Log.fatal "Utls.zstd_compress_file: command failed: %s" cmd in
     exit 1
  );
  Log.info "compressed %s in %.3f s" fn dt;
  (fn ^ ".zst")

let zst_sfx = ".zst"
let zst_sfx_len = S.length zst_sfx

let zstd_uncompress_file fn =
  assert(S.ends_with fn zst_sfx);
  let cmd = sprintf "zstd -dqfk %s" fn in (* decompress quiet force keep *)
  let dt, ret = time_it Unix.system cmd in
  (if ret <> Unix.WEXITED 0 then
     let () = Log.fatal "Utls.zstd_uncompress_file: command failed: %s" cmd in
     exit 1
  );
  Log.info "uncompressed %s in %.3f s" fn dt;
  (* uncompressed filename *)
  S.rchop ~n:zst_sfx_len fn

(* the stdlib one requires "true"/"false" *)
let bool_of_string = function
  | "1" -> true
  | "0" -> false
  | other -> failwith ("Mol.bool_of_string: " ^ other)

let int_of_bool = function
  | true -> 1
  | false -> 0

let string_of_array a sep to_str =
  let res = ref [] in
  let n = A.length a in
  for i = n - 1 downto 0 do
    res := to_str (A.unsafe_get a i) :: !res
  done;
  S.concat sep !res

let lines_to_file_a fn a =
  LO.with_out_file fn (fun out ->
      A.iter (fprintf out "%s\n") a
    )

(* ligand pqrs file is produced at runtime *)
let mol2pqrs nprocs in_mol2_fn out_pqrs_fn =
  let cmd = sprintf "mol2pqrs -np %d -lig %s -o %s" nprocs in_mol2_fn out_pqrs_fn in
  Log.info "running: %s" cmd;
  let ret = Unix.system cmd in
  if ret <> Unix.WEXITED 0 then
    let () = Log.fatal "Utls.mol2pqrs: command failed: %s" cmd in
    exit 1

(* int_of_float, but enforces that [x] has no fractional part *)
let int_of_float_exn (x: float): int =
  let frac, integ = modf x in
  assert(frac = 0.0);
  int_of_float integ

(* does [s1] dominate [s2]? [s1] and [s2] MUST BE incr. sorted;
 * return (s1_victories, s2_victories);
 * dominating means having lower docking scores; related to Mann-Whithey U
 * WARNING: brute-force O(N^2) *)
let dominate_BRUTE_FORCE (s1: float array) (s2: float array): int * int =
  let n = A.length s1 in
  assert(n = A.length s2);
  let left = ref 0 in
  let right = ref 0 in
  for i = 0 to n - 1 do
    let x = A.unsafe_get s1 i in
    for j = 0 to n - 1 do
      let y = A.unsafe_get s2 j in
      if x < y then
        incr left
      else if y < x then
        incr right
    done
  done;
  (!left, !right)

(* faster version: O(N*log(N)) *)
let dominate s1 s2 =
  let n = A.length s1 in
  let m = A.length s2 in
  assert(n = m);
  let left = ref 0 in
  let right = ref 0 in
  A.iter (fun (x: float) -> match A.bsearch BatFloat.ord s2 x with
      | `All_bigger -> left  := !left  + n
      | `All_lower  -> right := !right + n
      | `At i -> (left  := !left  + (n - (i + 1));
                  right := !right + i)
      | `Just_after i -> (left  := !left  + (n - (i + 1));
                          right := !right + (i + 1))
      | `Empty -> assert(false)
    ) s1;
  (!left, !right)

let temp_file prfx sfx =
  Filename.temp_file ~temp_dir:(Filename.get_temp_dir_name ()) prfx sfx

(* keep only two digits after the decimal-point character;
   rounds away from zero *)
let reduce_precision x =
  float (int_of_float (x *. 100.0 +.
                       (if x >= 0.0 then 0.5 else -.0.5)
                      )) /. 100.0
