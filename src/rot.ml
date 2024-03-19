(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

type t = { a: float ; b: float ; c: float ;
           d: float ; e: float ; f: float ;
           g: float ; h: float ; i: float }

let create a b c d e f g h i =
  {a; b; c;
   d; e; f;
   g; h; i}

let id (): t =
  { a=1.0 ; b=0.0 ; c=0.0 ;
    d=0.0 ; e=1.0 ; f=0.0 ;
    g=0.0 ; h=0.0 ; i=1.0 }

(* rotation matrix for rotation around X axis by theta(radians)
 * cf. http://motion.pratt.duke.edu/RoboticSystems/3DRotations.html
 * Wolfram alpha calls this a "Direction cosine matrix" *)
let rx theta: t =
  Math.(
    let t = cos_sin theta in
    { a=1.0; b=0.0      ; c=0.0   ;
      d=0.0; e=t.cos    ; f=t.sin ;
      g=0.0; h=(-.t.sin); i=t.cos }
  )

(* rotation matrix for rotation around Y axis by theta(radians) *)
let ry theta: t =
    Math.(
      let t = cos_sin theta in
      {a=t.cos; b=0.0; c=(-.t.sin) ;
       d=0.0  ; e=1.0; f=0.0       ;
       g=t.sin; h=0.0; i=t.cos     }
    )

(* rotation matrix for rotation around Z axis by theta(radians) *)
let rz theta: t =
  Math.(
    let t = cos_sin theta in
    { a=t.cos    ; b=t.sin; c=0.0 ;
      d=(-.t.sin); e=t.cos; f=0.0 ;
      g=0.0      ; h=0.0  ; i=1.0 }
  )

(* composition of three successive rotations by
 * [a] around Ox, then [b] around Oy, then [g] around Oz
 * R = Rz * Ry * Rx (applied from right to left)
 * formula from https://stackoverflow.com/questions/39128589/decomposing-rotation-matrix-x-y-z-cartesian-angles *)
let r_xyz alpha' beta' gamma' =
  Math.(
    let alpha = cos_sin alpha' in
    let beta  = cos_sin beta'  in
    let gamma = cos_sin gamma' in
    { a= beta.cos *. gamma.cos;
      b= gamma.cos *. alpha.sin *. beta.sin -. alpha.cos *. gamma.sin;
      c= alpha.sin *. gamma.sin +. alpha.cos *. gamma.cos *. beta.sin;
      d= beta.cos *. gamma.sin;
      e= alpha.cos *. gamma.cos +. alpha.sin *. beta.sin *. gamma.sin;
      f= alpha.cos *. beta.sin *. gamma.sin -. gamma.cos *. alpha.sin;
      g= -.beta.sin;
      h= beta.cos *. alpha.sin;
      i= alpha.cos *. beta.cos }
  )

(* return (alpha, beta, gamma) from rotation matrix [rot];
 * i.e. decomposition into Cartesian angles
 * formula also from https://stackoverflow.com/questions/39128589/decomposing-rotation-matrix-x-y-z-cartesian-angles *)
let decompose rot =
  let beta = atan2 (-.rot.g) (sqrt (rot.a *. rot.a +. rot.d *. rot.d)) in
  let c_beta = cos beta in
  (atan2 (rot.h /. c_beta) (rot.i /. c_beta), beta,
   atan2 (rot.d /. c_beta) (rot.a /. c_beta))

let mult (r1: t) (r2: t): t =
  (* r1: [|r1.a; r1.b; r1.c; *)
  (*       r1.d; r1.e; r1.f; *)
  (*       r1.g; r1.h; r1.i|] *)
  (* r2: [|r2.a; r2.b; r2.c; *)
  (*       r2.d; r2.e; r2.f; *)
  (*       r2.g; r2.h; r2.i|] *)
  { a= r1.a *. r2.a +. r1.b *. r2.d +. r1.c *. r2.g ;
    b= r1.a *. r2.b +. r1.b *. r2.e +. r1.c *. r2.h ;
    c= r1.a *. r2.c +. r1.b *. r2.f +. r1.c *. r2.i ;

    d= r1.d *. r2.a +. r1.e *. r2.d +. r1.f *. r2.g ;
    e= r1.d *. r2.b +. r1.e *. r2.e +. r1.f *. r2.h ;
    f= r1.d *. r2.c +. r1.e *. r2.f +. r1.f *. r2.i ;

    g= r1.g *. r2.a +. r1.h *. r2.d +. r1.i *. r2.g ;
    h= r1.g *. r2.b +. r1.h *. r2.e +. r1.i *. r2.h ;
    i= r1.g *. r2.c +. r1.h *. r2.f +. r1.i *. r2.i }

(* apply rotation [t] to 3D vector [v] *)
let rotate (r: t) (v: V3.t): V3.t =
  V3.{ x = r.a *. v.x +. r.b *. v.y +. r.c *. v.z ;
       y = r.d *. v.x +. r.e *. v.y +. r.f *. v.z ;
       z = r.g *. v.x +. r.h *. v.y +. r.i *. v.z }

(* cf. doc/rotation-matrix-of-quaternion.svg
 * extracted from wikipedia UNVERIFIED *)
let of_quat ({w; x; y; z}: Quat.t): t =
  let p2 = V3.make (x *. x) (y *. y) (z *. z) in
  {a = 1. -. 2. *. (p2.y +. p2.z);
   b =       2. *. (x *. y -. z *. w);
   c =       2. *. (x *. z +. y *. w);
   d =       2. *. (x *. y +. z *. w);
   e = 1. -. 2. *. (p2.x +. p2.z);
   f =       2. *. (y *. z -. x *. w);
   g =       2. *. (x *. z -. y *. w);
   h =       2. *. (y *. z +. x *. w);
   i = 1. -. 2. *. (p2.x +. p2.y)}

(* let axis_angle_of_quat {w = qr; x = qi; y = qj; z = qk} = *)
(*   failwith "not implemented yet" *)

(* rotation matrix from axis-angle rotation
 * p18 of doc/Quaternions_Orientation.pdf *)
let of_axis_angle ({ x; y; z }: V3.t) (theta: float): t =
(*  Math.(
    let t = cos_sin theta in
    let p2 = V3.make (x *. x) (y *. y) (z *. z) in
    let oneMct = 1.0 -. t.cos in
    { a= p2.x +. t.cos *. (1.0 -. p2.x);
      b= x *. y *. oneMct +. z *.t.sin;
      c= x *. z *. oneMct -. y *. t.sin;
      d= x *. y *. oneMct -. z *. t.sin;
      e= p2.y +. t.cos *. (1.0 -. p2.y);
      f= y *. z *. oneMct +. x *. t.sin;
      g= x *. z *. oneMct +. y *. t.sin;
      h= y *. z *. oneMct -. x *. t.sin;
      i= p2.z +. t.cos *. (1.0 -. p2.z) }
    ) *)
  let t = Math.cos_sin theta in
  let oneMct = 1.0 -. t.cos in
  { a = t.cos +. x*.x *. oneMct;
    b = x *. y *. oneMct -. z *. t.sin;
    c = x *. z *. oneMct +. y *. t.sin;
    d = x *. y *. oneMct +. z *. t.sin;
    e = t.cos +. y*.y *. oneMct;
    f = y *. z *. oneMct -. x *. t.sin;
    g = x *. z *. oneMct -. y *. t.sin;
    h = y *. z *. oneMct +. x *. t.sin;
    i = t.cos +. z*.z *. oneMct }

(* the trace is the sum of the matrix diagonal *)
let trace (r: t): float =
  r.a +. r.e +. r.i

(* http://motion.pratt.duke.edu/RoboticSystems/3DRotations.html
 * r11 r12 r13   a b c
 * r21 r22 r23 = d e f
 * r31 r32 r33   g h i *)
let to_axis_angle (r: t): V3.t * float =
  let theta = acos (((trace r) -. 1.0) *. 0.5) in
  (* WARNING: singularity if theta=0 *)
  let axis =
    let divisor = 1.0 /. (2.0 *. (sin theta)) in
    V3.make
      ((r.h -. r.f) *. divisor)
      ((r.c -. r.g) *. divisor)
      ((r.d -. r.b) *. divisor) in
  (axis, theta)

(* quaternion conversion *)
let to_quat r =
  Quat.of_axis_angle (to_axis_angle r)

(* as a quaternion: only 4 floats to log *)
let to_string r =
  Quat.to_string (to_quat r)

let of_string s =
  of_quat (Quat.of_string s)
