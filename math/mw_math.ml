open Utils
open Random


let pi = 4.0 *. atan 1.0;; (* FIXME! *)

let pow x n =
    exp ((float_of_int n) *. log(x) )
;;

let close_enough x y = 
    let diff = abs_float (x -. y) in
    let tolerance = 42. *. epsilon_float in 
    (diff <= tolerance *. (abs_float x)) || (diff <= tolerance *. (abs_float y))
;;

let standard_normal_density x =
    1. /. sqrt(2. *. pi) *. exp( -. x *. x /. 2.)
;;

class cumulativeNormalDistribution =
    object (self)
        method value x =

            (* 
             * Unknown, probably bad normal CDF approx. from
             * http://eternallearningq.wordpress.com/2012/04/01/cumulative-normal-distribution-in-qkdb/
             *
             * TODO: Look at quantlib's implementation for the Error Function.
             *)

            let a1 = 0.31938153 and 
                a2 = -0.356563782 and 
                a3 = 1.781477937 and 
                a4 = -1.821255978 and 
                a5 = 1.330274429 in

            let l  = abs_float(x) in
            let k  = 1.0 /. (1.0 +. 0.2316419 *. l) in
            let w  = ref (1.0-.1.0/.sqrt(2.0*.pi)*.exp(-.l*.l/.2.0)*.(a1*.k+.a2*.k*.k+.a3*. (pow k 3)+.a4*.(pow k 4)+.a5*.(pow k 5))) in
            if (x < 0.0) then  w := 1.0 -. !w ;
            !w

    end;;

class inverseCumulativeNormal = 
    object (self)

        (* Coefficients for the rational approximation *)
        val a1_ = -3.969683028665376e+01
        val a2_ =  2.209460984245205e+02
        val a3_ = -2.759285104469687e+02
        val a4_ =  1.383577518672690e+02
        val a5_ = -3.066479806614716e+01
        val a6_ =  2.506628277459239e+00

        val b1_ = -5.447609879822406e+01
        val b2_ =  1.615858368580409e+02
        val b3_ = -1.556989798598866e+02
        val b4_ =  6.680131188771972e+01
        val b5_ = -1.328068155288572e+01

        val c1_ = -7.784894002430293e-03
        val c2_ = -3.223964580411365e-01
        val c3_ = -2.400758277161838e+00
        val c4_ = -2.549732539343734e+00
        val c5_ =  4.374664141464968e+00
        val c6_ =  2.938163982698783e+00

        val d1_ =  7.784695709041462e-03
        val d2_ =  3.224671290700398e-01
        val d3_ =  2.445134137142996e+00
        val d4_ =  3.754408661907416e+00

        (* Define Break Points *)
        val x_low_ = 0.02425
        val x_high_= 1.0 -. 0.02425 


        method oneSample =
            Random.init  (int_of_float(Sys.time() *. 100000.));
            let rnd = Random.float 1.0 in
            self#standard_value rnd


        method tail_value x =
            let z = ref 0. in
            if x <= 0.0 || x >= 1.0 then 
                (* try to recover if due to numerical error *)
                if close_enough x 1.0 then
                    z := float_of_int max_int (* largest value available *)
                else if abs_float(x) < epsilon_float then
                    z := float_of_int min_int (* largest negative value available *)
                else 
                    raise (MW_FAIL ("inverseCumulativeNormal(" ^ (string_of_float x) ^ ") undefined: must be 0 < x < 1")) 
                    (*
                     *
                    print_newline()
                     * *)
            else begin
                if x < x_low_ then begin
                    (* Rational approximation for the lower region 0<x<u_low *)
                    z := sqrt( -2.0 *. log(x) );
                    z := (((((c1_*. !z +.c2_)*. !z +.c3_)*. !z +.c4_)*. !z +.c5_)*. !z +.c6_) /. ((((d1_*. !z +.d2_)*. !z +.d3_)*. !z +.d4_)*. !z +.1.0);

                end else begin
                    (* Rational approximation for the upper region u_high<x<1 *)
                    z := sqrt( -2.0 *. log(1.0 -. x) );

                    z := -.(((((c1_ *. !z +. c2_) *. !z +. c3_) *. !z +. c4_) *. !z +.
                    c5_) *. !z +. c6_) /. ((((d1_ *. !z +. d2_) *. !z +. d3_) *. !z +. d4_) *. !z +. 1.0);
                end
            end;

            !z


        method standard_value x =
            let z = ref 0. in
            if x < x_low_ || x_high_ < x then
                z := self#tail_value x
            else begin
                    (* Rational approximation for central region. *)
                    z := x -. 0.5;
                    let r = !z *. !z in
                    z := (((((a1_*.r+.a2_)*.r+.a3_)*.r+.a4_)*.r+.a5_)*.r+.a6_)*.
                    !z /. (((((b1_*.r+.b2_)*.r+.b3_)*.r+.b4_)*.r+.b5_)*.r+.1.0);
            end;

            (*
             * The relative error of the approximation has absolute value less
             * than 1.15e-9.  One iteration of Halley's rational method (third
             * order) gives full machine precision.
             * error (f_(z) - x) divided by the cumulative's derivative 
             *)

            let f_ = (new cumulativeNormalDistribution)#value in (* wow! OCAML *)
            let r = ((f_ !z) -. x) *. (sqrt 2.) *. (sqrt pi) *. (exp (0.5 *. !z *. !z)) in
            (* Halley's method *)
            z := !z -. r /. (1. +. 0.5 *. !z *. r);
            
            !z
    end;;



