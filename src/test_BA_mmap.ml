(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

module BA  = Bigarray
module BA1 = BA.Array1

open Printf

(* partial init of BA by child process *)
let child_fun ba start stop =
  for i = start to stop - 1 do
    BA1.set ba i (float i)
  done

let () =
  let tmp_fn = Filename.temp_file "" "" in
  let fd = Unix.(openfile tmp_fn [O_RDWR] 0o600) in
  let n = 10 in
  let ba1 = BA.array1_of_genarray (Unix.map_file fd BA.Float64 BA.c_layout true [|n|]) in
  let pid0 = Unix.fork () in
  (match pid0 with
   | 0 -> (* 1st child *)
     (child_fun ba1 0 (n/2);
      exit 0)
   | _ -> (* father *) ()
  );
  let pid1 = Unix.fork () in
  (match pid1 with
   | 0 -> (* 2nd child *)
     (child_fun ba1 (n/2) n;
      exit 0)
   | _ -> (* father *) ()
  );
  let _ = Unix.wait () in
  let _ = Unix.wait () in
  printf "father can see:\n";
  for i = 0 to n - 1 do
    printf "%f\n" (BA1.get ba1 i)
  done
