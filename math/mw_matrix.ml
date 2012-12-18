exception RecoverUnderflow

let iden n = (Array.init n (fun i ->
                Array.init n (fun j -> if i=j then 1. else 0.)))

class triQRDecomposition (diag:float list) (sub:float list) =
    object (self)
        val _diag = Array.of_list diag
        val  _sub = Array.of_list (0. :: sub)
        val _eigenVectors = iden (List.length diag)

        val mutable _iter = 0

        method getEigenValues =
            Array.to_list _diag

        method getEigenVectors =
            _eigenVectors

        initializer 
            let print_list (a: float array) =
                let pp x = 
                    print_float x;
                    print_string(", ");
                in
                Array.iter pp a 
            in

            let n = List.length diag in

            let offDiagIsZero k e =
                abs_float e.(k) < epsilon_float
            in

            (* FIXME: Should use Array.mapi *)
            for k = n-1 downto 1 do
                while not (offDiagIsZero k _sub) do
                    let l = ref (pred k) in

                    while (!l > 0 && not (offDiagIsZero !l _sub)) do
                        l := pred !l
                    done;
                    _iter <- succ _iter;

                    let q = ref _diag.(!l) in

                    (* calculated eigenvalue of 2x2 sub matrix of
                     * [ d_[k-1] e_[k] ]
                     * [  e_[k]  d_[k] ]
                     * which is closer to d_[k+1].
                     * FLOATING_POINT_EXCEPTION
                     *)

                    let t1 = ( sqrt (0.25 *. (_diag.(k) *. _diag.(k) +. _diag.(k-1) *. _diag.(k-1)) -. 0.5 *. _diag.(k-1) *. _diag.(k) +. _sub.(k) *. _sub.(k)) ) in
                    let t2 = 0.5 *. (_diag.(k) +. _diag.(k-1)) in

                    let lambda = (
                        if ( abs_float (t2 +. t1 -. _diag.(k)) < abs_float (t2 -. t1 -. _diag.(k))) then t2 +. t1
                        else t2 -. t1) in

                    q := !q -.  (if (k == (n - 1)) then 1.25 else 1.0) *. lambda;


                    (* the QR transformation *)
                    let sine = ref 1.0 in
                    let cosine = ref 1.0 in
                    let u = ref 0.0 in

                    try
                        for i = !l + 1 to k do
                            let h = ref (!cosine *. _sub.(i)) in
                            let p = ref (!sine *. _sub.(i)) in

                            _sub.(i-1) <- sqrt (!p *. !p +. !q *. !q);
                            if (_sub.(i-1) != 0.0) then begin
                                sine := !p /. _sub.(i-1);
                                cosine := !q /. _sub.(i-1);

                                let g = ref (_diag.(i-1) -. !u) in
                                let t = ref ((_diag.(i) -. !g) *. !sine +. 2. *. !cosine *. !h) in

                                u := !sine *. !t;
                                _diag.(i-1) <- !g +. !u;
                                q := !cosine *. !t -. !h;

                                for j = 0 to n - 1 do (* FIXME! *)
                                    let tmp = _eigenVectors.(j).(i-1) in
                                    _eigenVectors.(j).(i-1) <- !sine *. _eigenVectors.(j).(i) +. !cosine *. tmp;
                                    _eigenVectors.(j).(i) <- !cosine *. _eigenVectors.(j).(i) -. !sine *. tmp;
                                done;

                            end else begin
                                (* recover from underflow *)
                                _diag.(i-1) <- _diag.(i-1) -. !u;
                                _sub.(!l) <- 0.0;
                                raise RecoverUnderflow
                                (* recoverUnderflow = true; *)
                            end;
                        done;

                        _diag.(k) <- _diag.(k) -. !u;
                        _sub.(k) <- !q;
                        _sub.(!l) <- 0.0;

                    with RecoverUnderflow ->
                        print_string("!!!!!!")
                done
            done

            (* sort (eigenvalues, eigenvectors)

            std::vector<std::pair<Real, std::vector<Real> > > temp(n);
            std::vector<Real> eigenVector(ev_.rows());

            for (i=0; i<n; i++) {
                if (ev_.rows() > 0)
                    std::copy(ev_.column_begin(i),
                              ev_.column_end(i), eigenVector.begin());
                temp[i] = std::make_pair(d_[i], eigenVector);
            }
            std::sort(temp.begin(), temp.end(),
                      std::greater<std::pair<Real, std::vector<Real> > >());
            // first element is positive
            for (i=0; i<n; i++) {
                d_[i] = temp[i].first;
                Real sign = 1.0;
                if (ev_.rows() > 0 && temp[i].second[0]<0.0)
                    sign = -1.0;
                for (Size j=0; j<ev_.rows(); ++j) {
                    ev_[j][i] = sign * temp[i].second[j];
                }
            }
            *)
    end;;




                    (*
                Real q = d_[l];
                if (strategy != NoShift) {
                    const Real t1 = std::sqrt(
                      0.25*(d_[k]*d_[k] + d_[k-1]*d_[k-1])
                      - 0.5*d_[k-1]*d_[k] + e[k]*e[k]);
                    const Real t2 = 0.5*(d_[k]+d_[k-1]);

                    const Real lambda =
                        (std::fabs(t2+t1 - d_[k]) < std::fabs(t2-t1 - d_[k]))?
                        t2+t1 : t2-t1;

                    if (strategy == CloseEigenValue) {
                        q-=lambda;
                    } else {
                        q-=((k==n-1)? 1.25 : 1.0)*lambda;
                    }
                }
                *)







(* FOR FUTURE USE, MAYBE

module Linear = struct
let len = Array.length and aini = Array.init

let print m mat = let i = len mat and j = len mat.(0) in
print_string m; print_newline ();
  for i = 0 to i-1 do for j = 0 to j-1 do print_float mat.(i).(j); print_char ' ' done;
  print_newline () done

(* OCaml code to invert a matrix *)
exception Singular
let inv mat =
let l = len mat in
let am = Array.mapi (fun m row -> (Array.append row
      (aini l (fun n -> if m=n then 1. else 0.)))) mat in
for i = 0 to l-1 do (let im = ref 0 and mv = ref (abs_float am.(i).(i)) in
   for j = i+1 to l-1 do (let ae = abs_float am.(j).(i) in
       if (!mv < ae) then (mv := ae; im := j)) done;
   if !mv = 0. then raise Singular;
   if !im > i then (for n = i to (2*l - 1) do
      (let s = am.(i).(n) in am.(i).(n) <- am.(!im).(n); am.(!im).(n) <- s) done);
   let r = 1. /. am.(i).(i) in
   for j = i to 2*l - 1 do (am.(i).(j) <- r *. am.(i).(j)) done;
   for k = i+1 to l-1 do (let f = am.(k).(i) in
      for j = i+1 to 2*l - 1 do (am.(k).(j) <- am.(k).(j) -. f *. am.(i).(j))
      done); done) done;
for i = 0 to l-1 do (for j = i+1 to l-1 do (let p = am.(i).(j) in
      for k = i+1 to 2*l - 1 do
         (am.(i).(k) <- am.(i).(k) -. am.(j).(k) *. p) done) done) done;
Array.map (fun row -> Array.sub row l l) am

(* Matrix Multiply *)
exception LopSided
let mul a b =
   let l1 = len a and l2 = len a.(0)
   and l3 = len b and l4 = len b.(0) in
   if l2 <> l3 then raise LopSided;
   aini l1 (fun i -> aini l4 (fun k ->
     let c = ref 0. in for j = 0 to l2-1 do
        c := !c +. a.(i).(j) *. b.(j).(k) done; !c))

let determ a = (let s = len a in let w = Array.make (s * s) 0. in
  for i = 0 to s-1 do (for j = 0 to s-1 do (w.(s*i + j) <- a.(i).(j)) done) done;
  let rec det b k = (if k = 1 then w.(b)
    else (let rec piv j = if j=k then 0.
     else (if w.(b + s*j) = 0. then piv (j+1)
       else (if j <> 0 then (let q = b + j*s in
         for n = 0 to k-1 do (let t = w.(q + n) in
             w.(q + n) <- w.(b+n); w.(b+n) <- -. t) done);
         (let t = 1. /. w.(b) in for i = 1 to k-1 do (
            let q = b + i*s in let a = t *. w.(q) in
               for r=1 to k-1 do (w.(q+r) <- w.(q+r) -. a *. w.(b+r)) done) done);
          w.(b) *. (det (b+s+1) (k-1)))) in piv 0)) in det 0 s)

exception Ragged
let trans a = let m = len a and n = len a.(0) in
  for i = 1 to m-1 do if len a.(i) <> n then raise Ragged done;
   aini n (fun i -> aini m (fun j -> a.(j).(i)))

let ip a b = let n = len a and s = ref 0. in
   if n <> len b then raise Ragged;
   for i = 0 to n-1 do s := !s +. a.(i) *. b.(i) done; !s

let mtv a b = let l1 = len a.(0) and l2 = len b in
   if l1 <> l2 then raise LopSided;
   aini (len a) (fun j -> ip a.(j) b)

let vneg a = aini (len a) (fun j -> -. a.(j))

let sm sc a = aini (len a) (fun j -> sc *. a.(j))

let vadd a b = let l = (len a) in assert (l = (len b)); aini l (fun j -> a.(j) +. b.(j))

let iden n = (Array.init n (fun i -> Array.init n (fun j -> if i=j then 1. else 0.)))

let gs dp bas = let n = len bas and m = len bas.(0) in
  let nb = Array.make n [||] in for i = 0 to n-1 do nb.(i) <- Array.copy bas.(i);
  for k = 0 to i-1 do let ip = dp nb.(i) nb.(k) in for j = 0 to m-1 do
      nb.(i).(j) <- nb.(i).(j) -. ip *. nb.(k).(j) done done;
    nb.(i) <- sm (1. /. (sqrt (dp nb.(i) nb.(i)))) nb.(i) done; nb

let eigen m = let n = len m in let m2 = aini n (fun j -> Array.copy m.(j)) in
(let id = iden n in (let mOd = ref 0. and ip = ref (-1, -1) in
  while mOd := 0.; ip := (-1, -1);
  for i=0 to n-1 do for j=i+1 to n-1 do let q = (abs_float m2.(i).(j)) in
     if q > !mOd then (mOd := q; ip := (i, j)) done done;
     !mOd > 0.00000000001 do let (i, j) = !ip in
     let th = (0.5 *. (atan (2. *. (m2.(i).(j) /. (m2.(i).(i) -. m2.(j).(j)))))) in
     (let c = cos th and s = sin th in let twst m = (let tmp = vadd (sm c m.(i)) (sm s m.(j)) in
       m.(j) <- vadd (sm (-. s) m.(i)) (sm c m.(j)); m.(i) <- tmp) in
    for k=0 to n-1 do (let tmp = c *. m2.(k).(i) +. s *. m2.(k).(j) in
         m2.(k).(j) <- (-. s) *. m2.(k).(i) +. c *.  m2.(k).(j);  m2.(k).(i) <- tmp) done;
    (twst m2; twst id)) done;
     aini n (fun j -> m2.(j).(j), id.(j))))

end;;
*)



(*

    TqrEigenDecomposition::TqrEigenDecomposition(const Array& diag,
                                                 const Array& sub,
                                                 EigenVectorCalculation calc,
                                                 ShiftStrategy strategy)
    : iter_(0), d_(diag),
      ev_((calc == WithEigenVector)? d_.size() :
          (calc == WithoutEigenVector)? 0 : 1, d_.size(), 0)
    {
        Size n = diag.size();

        QL_REQUIRE(n == sub.size()+1, "Wrong dimensions");

        Array e(n, 0.0);
        std::copy(sub.begin(),sub.end(),e.begin()+1);
        Size i;
        for (i=0; i < ev_.rows(); ++i) {
            ev_[i][i] = 1.0;
        }

        for (Size k=n-1; k >=1; --k) {
            while (!offDiagIsZero(k, e)) {
                Size l = k;
                while (--l > 0 && !offDiagIsZero(l,e));
                iter_++;

                Real q = d_[l];
                if (strategy != NoShift) {
                    // calculated eigenvalue of 2x2 sub matrix of
                    // [ d_[k-1] e_[k] ]
                    // [  e_[k]  d_[k] ]
                    // which is closer to d_[k+1].
                    // FLOATING_POINT_EXCEPTION
                    const Real t1 = std::sqrt(
                                          0.25*(d_[k]*d_[k] + d_[k-1]*d_[k-1])
                                          - 0.5*d_[k-1]*d_[k] + e[k]*e[k]);
                    const Real t2 = 0.5*(d_[k]+d_[k-1]);

                    const Real lambda =
                        (std::fabs(t2+t1 - d_[k]) < std::fabs(t2-t1 - d_[k]))?
                        t2+t1 : t2-t1;

                    if (strategy == CloseEigenValue) {
                        q-=lambda;
                    } else {
                        q-=((k==n-1)? 1.25 : 1.0)*lambda;
                    }
                }

                // the QR transformation
                Real sine = 1.0;
                Real cosine = 1.0;
                Real u = 0.0;

                bool recoverUnderflow = false;
                for (Size i=l+1; i <= k && !recoverUnderflow; ++i) {
                    const Real h = cosine*e[i];
                    const Real p = sine*e[i];

                    e[i-1] = std::sqrt(p*p+q*q);
                    if (e[i-1] != 0.0) {
                        sine = p/e[i-1];
                        cosine = q/e[i-1];

                        const Real g = d_[i-1]-u;
                        const Real t = (d_[i]-g)*sine+2*cosine*h;

                        u = sine*t;
                        d_[i-1] = g + u;
                        q = cosine*t - h;

                        for (Size j=0; j < ev_.rows(); ++j) {
                            const Real tmp = ev_[j][i-1];
                            ev_[j][i-1] = sine*ev_[j][i] + cosine*tmp;
                            ev_[j][i] = cosine*ev_[j][i] - sine*tmp;
                        }
                    } else {
                        // recover from underflow
                        d_[i-1] -= u;
                        e[l] = 0.0;
                        recoverUnderflow = true;
                    }
                }

                if (!recoverUnderflow) {
                    d_[k] -= u;
                    e[k] = q;
                    e[l] = 0.0;
                }
            }
        }

        // sort (eigenvalues, eigenvectors),
        // code taken from symmetricSchureDecomposition.cpp
        std::vector<std::pair<Real, std::vector<Real> > > temp(n);
        std::vector<Real> eigenVector(ev_.rows());
        for (i=0; i<n; i++) {
            if (ev_.rows() > 0)
                std::copy(ev_.column_begin(i),
                          ev_.column_end(i), eigenVector.begin());
            temp[i] = std::make_pair(d_[i], eigenVector);
        }
        std::sort(temp.begin(), temp.end(),
                  std::greater<std::pair<Real, std::vector<Real> > >());
        // first element is positive
        for (i=0; i<n; i++) {
            d_[i] = temp[i].first;
            Real sign = 1.0;
            if (ev_.rows() > 0 && temp[i].second[0]<0.0)
                sign = -1.0;
            for (Size j=0; j<ev_.rows(); ++j) {
                ev_[j][i] = sign * temp[i].second[j];
            }
        }
    }

    // see NR for abort assumption as it is
    // not part of the original Wilkinson algorithm
*)



