(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

include BA.Array1

type ba1f = (float, BA.float32_elt, BA.c_layout) BA.Array1.t

let iter f a =
  let n = dim a in
  for i = 0 to n - 1 do
    f (unsafe_get a i) (* this one is safe *)
  done

let min a =
  let res = ref infinity in
  iter (fun x ->
      res := BatFloat.min x !res
    ) a;
  !res

let max a =
  let res = ref neg_infinity in
  iter (fun x ->
      res := BatFloat.max x !res
    ) a;
  !res

let sum a =
  let res = ref 0.0 in
  iter (fun x ->
      res := !res +. x
    ) a;
  !res

let favg a =
  (sum a) /. (float (dim a))

let sparsity a =
  let non_zeroes = ref 0 in
  iter (fun x ->
      if x <> 0.0 then
        incr non_zeroes
    ) a;
  (float !non_zeroes) /. (float (dim a))
