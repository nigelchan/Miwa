open Utils
open Mw_math
open Payoff
open Process


class pathsGenerator (process:stochasticProcess) expiry timeSteps = 
    object (self)
        val _process = process
        val _timeSteps = timeSteps
        val _timeSpacing = expiry /. (float timeSteps) (* FIXME: Should be replaced by _timeGrid later *)

        method generateOnePath =

            let x0 = _process#x0 in
            let dim = _process#factors in

            let simulate n = 
                let rec simulate_once n (accum: multiPath) =
                    if List.length accum == n then accum
                    else begin
                        let z = Array.init dim (fun n -> (new inverseCumulativeNormal)#oneSample) in
                        let s_t_dt = _process#evolve (List.hd accum) z _timeSpacing in
                        simulate_once n (s_t_dt :: accum);
                    end;
                in 
                simulate_once n [x0]
            in 

            simulate (timeSteps + 1)
    end;;



            
class pathPricer (payoff:europeanPayoff) discountFactor =

    object (self)
        val _discountFactor = discountFactor
        val _payoff = payoff

        method getValue (a_path:multiPath) =
            _discountFactor *. _payoff#price(a_path)

    end;;

