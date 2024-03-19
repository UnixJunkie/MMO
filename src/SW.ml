(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
 *
 * a Sliding Window *)

type t = { size: int; (* window length *)
           accepts: int ref;
           rejects: int ref;
           (* an accept event is true; a reject is false *)
           events: bool Queue.t }

let create n =
  { size = n;
    accepts = ref 0;
    rejects = ref 0;
    events = Queue.create () }

(* process an event *)
let process sw evt =
  Queue.push evt sw.events;
  (if evt then
     incr sw.accepts
   else
     incr sw.rejects
  );
  (* maintain window size during cruise *)
  if Queue.length sw.events > sw.size then
    let too_old = Queue.pop sw.events in
    if too_old then
      decr sw.accepts
    else
      decr sw.rejects

let get_ratio sw =
  float !(sw.accepts) /. float (!(sw.accepts) + !(sw.rejects))

(* if we need to reset statistics *)
let reset sw =
  sw.accepts := 0;
  sw.rejects := 0;
  Queue.clear sw.events
