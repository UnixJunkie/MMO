(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

let verbose = ref false

(* for computers with lots of disk space *)
let no_compress = ref false

(* Nlopt global optimization algorithm choice *)
let glob_opt = ref 0
