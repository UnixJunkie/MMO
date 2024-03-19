(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

include Dolog.Log

(* ANSI terminal colors for UNIX *)
let black   = "\027[30m"
let red     = "\027[31m"
let green   = "\027[32m"
let yellow  = "\027[33m"
let blue    = "\027[34m"
let magenta = "\027[35m"
let cyan    = "\027[36m"
let white   = "\027[37m"
let default = "\027[39m"
let reset   = "\027[0m"

let color c s =
  Printf.sprintf "%s%s%s" c s reset
