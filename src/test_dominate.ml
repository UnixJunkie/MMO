
(* tests for Utls.dominate and dominate_brute *)

let get_rands rng =
  let res = A.init 1000 (fun _i -> Random.State.float rng 1.0) in
  A.sort BatFloat.compare res;
  res

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let rng = Random.State.make_self_init () in
  let s1 = get_rands rng in
  let s2 = get_rands rng in
  let left, right = Utls.dominate [|0.;1.5;2.5;3.5|] [|0.;1.;2.;3.|] in
  assert(left = 6 && right = 9);
  let dt0, (l0, r0) = Utls.time_it (Utls.dominate_BRUTE_FORCE s1) s2 in
  Log.info "dt=%f; l0, r0: %d, %d" dt0 l0 r0;
  assert(l0 + r0 = 1_000_000);
  let dt1, (l1, r1) = Utls.time_it (Utls.dominate s1) s2 in
  Log.info "dt=%f; l1, r1: %d, %d" dt1 l1 r1;
  assert(l1 + r1 = 1_000_000);
  assert(l0 = l1);
  assert(r0 = r1)

let () = main ()
