(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

(* x in [0..255] *)
let rand_byte rng =
  Random.State.int rng 256

let split rng_state =
#if OCAML_VERSION >= (5, 0, 0)
      Random.State.split rng_state
#else
  (* this is not a very good idea, but allows to compile
   * w/ ocaml < 5.0.0; the two RNGs might be correlated *)
  Random.State.make (A.init 20 (fun _i -> rand_byte rng_state))
#endif

(* 160 bits of entropy read from /dev/urandom *)
let entropy_160b (): int array =
  LO.with_in_file "/dev/urandom" (fun input ->
      A.init 20 (fun _i -> input_byte input)
    )

let make = Random.State.make
