(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

(* float bigarray used as a shared memory (shm) *)

module BA  = BatBigarray
module BA1 = BA.Array1

type ba1f = (float, BA.float64_elt, BA.c_layout) BA.Array1.t

type t = { tmp_fn: string;
           tmp_fd: Unix.file_descr;
           ba1: ba1f }

let create (n: int): t =
  let tmp_fn = Filename.temp_file "MMCLD_shm_" "" in
  Log.info "Shm: %s" tmp_fn;
  let tmp_fd = Unix.(openfile tmp_fn [O_RDWR] 0o600) in
  let ba1 =
    BA.array1_of_genarray
      (Unix.map_file tmp_fd BA.Float64 BA.c_layout true [|n|]) in
  { tmp_fn; tmp_fd; ba1 }

(* copy to a new float array *)
let to_fresh_array (x: t): float array =
  BA1.to_array x.ba1

let destroy (x: t): unit =
  Unix.close x.tmp_fd;
  Unix.unlink x.tmp_fn
