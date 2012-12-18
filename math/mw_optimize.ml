open Gsl
open Root.Bracket

class brent =
    object(self)

        method solve f ?(max_iter=5000) ?(accuracy=1e-12) (initial:float) min max =

            let solver = Root.Bracket.make Root.Bracket.BRENT f min max in
            let r = ref 0. in

            let rec proc i = function
            | true -> !r
            | _ when i >= max_iter -> nan
            | _ ->
                iterate solver; 
                r := root solver;

                let (x_lo, x_hi) = interval solver in
                let status = Root.test_interval ~lo:x_lo ~up:x_hi ~epsabs:0. ~epsrel:accuracy in
                proc (succ i) status
            in  
            proc 1 false;


    end;;
