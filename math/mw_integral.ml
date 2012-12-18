open Mw_matrix
open Printf

class gaussianLaguerrePolynomial = 
    object (self)

        method value n x =
            let rec cal n x =
                match n with
                    0 -> 1.
                |   1 -> x -. (self#alpha 0)
                |   _ -> (x -. (self#alpha n)) *. (cal (n-1) x) -. (self#beta (n-1)) *. (cal (n-2) x )
            in
            cal n x

        method alpha n = 
            float (2 * n + 1)

        method beta n = 
            float (n * n)

        method mu_0 =
            1.

        method w x =
            exp (-. x)

    end;;

let sumProduct (a: float list) (b:float list) =
    List.fold_left ( +. ) 0. (List.map2 ( *. )  a b)
;;

let print_list (a: float list) =
    let pp x = 
        print_float x;
        print_string(" ");
        print_newline()
    in
    List.iter pp a 
;;

(* Returns [i, j) *)
let rec range i j = if i > (j-1) then [] else i :: (range (i+1) j);;

class gaussianLaguerreIntegration n = 
    object (self)

        val _poly = new gaussianLaguerrePolynomial
        val mutable _weight = ([] : float list)
        val mutable _x = ([] : float list)

        method integrate f =
            (* Not necessary to provide range and default is [0, +oo) *) 
            sumProduct _weight (List.map f _x)


        initializer
            let diag = List.map _poly#alpha (range 0 n) in
            let sub = List.map (fun x -> sqrt (_poly#beta x)) (range 1 n) in

            let tqr = new triQRDecomposition diag sub in

            _x <- tqr#getEigenValues;

            let eigenVectors = tqr#getEigenVectors in
            let mu_0 = _poly#mu_0 in
            (* print_list _x; *)

            let first_ev = Array.to_list tqr#getEigenVectors.(0) in

            let t = List.map (fun x -> x *. x *. mu_0) first_ev in
            _weight <- List.map2 ( /. ) t (List.map _poly#w _x);

            (*
            print_newline();
            print_list _weight
            *)

    end;;
