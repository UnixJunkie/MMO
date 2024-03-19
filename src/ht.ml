(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

include BatHashtbl

let find_exn string_of_key ht k =
  try find ht k
  with Not_found ->
    let () = Log.fatal "Ht.find_exn: Not_found: %s" (string_of_key k) in
    raise Not_found
