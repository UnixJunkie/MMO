(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

include BatArray

(* divide array into [nparts]
   (the last part(s) might have less elements than previous ones) *)
let separate (nparts: int) (arr: 'a array): 'a array array =
  if nparts = 1 then
    [|arr|]
  else
    let n = length arr in
    assert(n > 0 && nparts > 1);
    (* max length of an output array *)
    let len = int_of_float (ceil ((float n) /. (float nparts))) in
    init nparts (fun i ->
        let start = i * len in
        let stop = Stdlib.min (start + len) n in
        let delta = stop - start in
        if start >= n || start + delta > n then
          [||]
        else
          sub arr start delta
      )

let concat_a (arrs: 'a array array): 'a array =
  let n = fold (fun acc arr -> acc + length arr) 0 arrs in
  let res = make n (arrs.(0)).(0) in
  let i = ref 0 in
  iter (fun arr ->
      iter (fun x ->
          res.(!i) <- x;
          incr i
        ) arr
    ) arrs;
  res

(* <=> BatList.takedrop *)
let takedrop arr m =
  let n = length arr in
  if m < n then
    (left arr m, sub arr m (n - m))
  else
    (arr, [||])

let enumerate arr =
  mapi (fun i x -> (i, x)) arr

let is_sorted a =
  let n = length a in
  if n <= 1 then
    true
  else
    let i = ref 1 in
    while !i < n && unsafe_get a (!i - 1) <= unsafe_get a !i do
      incr i
    done;
    (!i = n)

let min_avg_max a =
  (min a, favg a, max a)
