open Utils
open Gsl
open Fun

open Mw_math
open Mw_optimize

open Engine
open Process
open Payoff 
open Model

class europeanOption ?(engine: engine option)
                     (payoff:europeanPayoff)
                     (expiry:float) =
    object (self)

        val mutable _engine = engine
        val _payoff = payoff

        method getPrice =
            match _engine with
                None -> raise (MW_FAIL "No engine is provided")
            |   Some e -> e#calculate

        method setEngine e =
            _engine <- Some e

    end;;


let blackFormula ?(discount=1.) (optionType:optionType) (f:float) (k:float) (stdev:float)  =

    (* Note: stdev is vol * sqrt(t) *)

    let d1 = log( f /. k ) /. stdev  +. 0.5 *. stdev  in
    let d2 = d1 -. stdev in

    let cdf = new cumulativeNormalDistribution in


    let optionInd = match optionType with
                    |   Call -> 1.
                    |   Put -> -1. in

    let nd1 = cdf#value (d1 *. optionInd) in
    let nd2 = cdf#value (d2 *. optionInd) in

     discount *. optionInd *. (f *. nd1 -. k *. nd2)

;;

(* An instrument to which hestonCalibratedModel calibrates *)
class hestonCalibrateInstrument (impVol:float)
                                (s0:float)
                                (strike:float)
                                (rfRate:float)
                                (expiry:float) =

    let db_payoff = new europeanPayoff Call strike in 

    object (self)

        val _impVol = impVol
        val _s0 = s0
        val _strike = strike
        val _rfRate = rfRate
        val _expiry = expiry

        val mutable _marketValue = 0.

        val discount = (1. /. ((1. +. rfRate) ** expiry) ) (* Should be replaced by forward rate object and its methods *)
        val instrument = new europeanOption db_payoff expiry


        method blackPrice vol = 
            let stdev = vol *. (sqrt expiry) in
            blackFormula Call _s0 (_strike *. discount) stdev

        method calculateImpVol modelPrice minVol maxVol = 
            (*
            printf "Solving for implied Vol...";
            print_newline();
            *)

            let f x = (self#blackPrice x) -. modelPrice in 
            let solver = new brent in
            solver#solve f ~max_iter:5000 ~accuracy:1e-12 _impVol minVol maxVol

        method modelValue =
            instrument#getPrice

        method marketValue =
            _marketValue

        method setEngine e =
            let ee = Oo.copy e in

            ee#setupArguments {optionType = Call; strike = _strike; expiry = _expiry};
            instrument#setEngine ee


        method calibrationError =
            (*
            printf "Calculating calibrationError...";
            print_newline();
            *)

            let lowerPrice = self#blackPrice 0.001 in
            let upperPrice = self#blackPrice 10. in
            let modelPrice = ref self#modelValue in

            (*
            print_float !modelPrice;
            print_newline();
            *)

            let calibratedImpVol = ref 0. in
            if (!modelPrice <= lowerPrice) then
                calibratedImpVol := 0.001
            else
                if (!modelPrice >= upperPrice) then
                    calibratedImpVol := 10.0
                else
                    calibratedImpVol := self#calculateImpVol !modelPrice 0.001 10.;

            (*
            printf "calibratedImpVol = %f" !calibratedImpVol;
            print_newline();
            *)

            !calibratedImpVol -. _impVol

    end;;

class optimizationConstraint =
    object (self)


        method test (x: float array) =
            let flag = ref true in
            if x.(0) <= 0. then flag := false;
            if x.(1) <= 0. then flag := false;
            if x.(2) <= 0. then flag := false;
            if x.(3) > 1. || x.(3) < -1. then flag := false;
            if x.(4) <= 0. then flag := false;

            !flag
    end;;

class costFunction (model:hestonProcess model) (instruments: hestonCalibrateInstrument array) =
    object(self) 

        val _model = model
        val _instruments = instruments 

        method gradient (params: float array) =
            let p = Array.length params in
            let n = Array.length _instruments in
            let j = Matrix.create n p in

            let dx = 1e-4 in (* FIXME *)

            let f = ref (self#individualCosts params) in

            for l=0 to pred p do
                let xdx = Vector.of_array params in
                xdx.{l} <- xdx.{l} +. dx;

                let f_forward = ref (self#individualCosts (Vector.to_array xdx)) in

                for k=0 to pred n do
                    j.{k, l} <- (!f_forward.(k) -. !f.(k)) /. dx
                done
            done;
            j

        (* returns sum of all cost values *) 
        method totalCost (params: float array) =
            _model#setParams params;

            let sumCost (accum:float) (instru:hestonCalibrateInstrument) =
                accum +. instru#calibrationError;
            in
            Array.fold_left sumCost 0. _instruments 

        (* return an array of cost value *)
        method individualCosts (params: float array) = 
            (*
            print_array params;
            *)
            _model#setParams params;
            Array.map (fun (instru:hestonCalibrateInstrument) -> instru#calibrationError) _instruments

        method totalCostAndGradient (params: float array) =
            ()

    end;;

class optimizationProblem (costfcn: costFunction) (optConstraint: optimizationConstraint) (initialValues: float array) =
    object (self)

        val _costfcn = costfcn
        val _optConstraint = optConstraint
        val _initialValues = initialValues

        method optConstraint = _optConstraint
        method costfcn = _costfcn
        method initialValues = _initialValues
    end;;



class levenbergMarquardt (problem: optimizationProblem) =
    object(self)

        val _problem = problem
        val _initialCostValues = (problem#costfcn)#individualCosts problem#initialValues
        val _initialGradient = (problem#costfcn)#gradient problem#initialValues

        (* returns the GSL's multi_fun_fdf 
         *
         *  type multi_fun_fdf = { 
         *      multi_f   : x:vector -> f:vector -> unit ;
         *      multi_df  : x:vector -> j:Matrix.matrix -> unit ;
         *      multi_fdf : x:vector -> f:vector -> j:Matrix.matrix -> unit ;
         *    }
         *
         *)

        method getFunctions =
            let opt_f ~x ~f =
                let n = Vector.length f in
                let x_array = Vector.to_array x in
                (* assert(Array.length y = n); *)

                let diff = Vector.copy x in
                Vector.sub diff (Vector.of_array problem#initialValues);

                if ((_problem#optConstraint)#test x_array) && (Array.fold_left ( +. )  0. (Vector.to_array diff)) >= 1e-4  then
                    let indCostsArray = ref ((_problem#costfcn)#individualCosts x_array) in
                    for i=0 to pred n do
                        f.{i} <- !indCostsArray.(i) (* FIXME: write a wrapper that take an Array and SET it in Bigarray *)
                    done
                else
                    for i=0 to pred n do
                        f.{i} <- _initialCostValues.(i) (* FIXME *)
                    done
            in

            let opt_df ~x ~j =

                let (n,p) = Matrix.dims j in 
                let x_array = Vector.to_array x in

                assert(Vector.length x = p);

                let diff = Vector.copy x in
                Vector.sub diff (Vector.of_array problem#initialValues);
                let diff_value = (Array.fold_left ( +. )  0. (Vector.to_array diff)) in

                if ((_problem#optConstraint)#test x_array) && diff_value >= 1e-4 then
                    let grad = ref ((_problem#costfcn)#gradient x_array) in
                    for l=0 to pred p do
                        for k=0 to pred n do
                            j.{k, l} <- !grad.{k, l}
                        done
                    done
                else
                    for l=0 to pred p do
                        for k=0 to pred n do
                            j.{k, l} <- _initialGradient.{k, l}
                        done
                    done
            in

            let opt_fdf ~x ~f ~j =
                opt_f ~x ~f;
                opt_df ~x ~j
            in

            {
              multi_f   = opt_f;
              multi_df  = opt_df;
              multi_fdf = opt_fdf;
            }

        method minimize =
            self#run _initialCostValues _problem#initialValues


        method run (y: float array) xinit =

            let print_state n p =
                let x = Vector.create p in
                let f = Vector.create n in
                fun iter s ->
                    Multifit_nlin.get_state s ~x ~f ();
                    Printf.printf "iter: %3u x = %15.8f % 15.8f % 15.8f % 15.8f % 15.8f"
                    iter x.{0} x.{1} x.{2} x.{3} x.{4};
                    print_newline()
            in

            let n = Array.length y in
            let p = Array.length xinit in
            let solver = Multifit_nlin.make Multifit_nlin.LMSDER ~n ~p (self#getFunctions) (Vector.of_array xinit) in  

            let rec proc iter = 
                Multifit_nlin.iterate solver; 
                let print_state = print_state n p in
                print_state iter solver; 

                let status = Multifit_nlin.test_delta solver ~epsabs:1e-4 ~epsrel:1e-4 in
                let maxiter = 5000 in

                match status with
                |   true -> Printf.printf "\nstatus = converged\n"
                |   false when iter >= maxiter -> Printf.printf "\nstatus = too many iterations\n"
                |   false -> proc (succ iter)
            in  
            Printf.printf "%-14s" " ";
            Printf.printf "%15s %15s %15s %15s %15s" "Theta" "Kappa" "Sigma" "Rho" "v0";
            print_newline();
            proc 1 ; 

    end;;

class hestonCalibratedModel (process: hestonProcess) =
    object (self)

        inherit ['a] model process


        method calibrate (calibrateInstruments: hestonCalibrateInstrument array) =
            printf "Start calibration:";
            print_newline();

            let optConstraint = new optimizationConstraint in (* FIXME *)

            let f = new costFunction (self :> hestonProcess model) calibrateInstruments in
            let op = new optimizationProblem f optConstraint _params in

            let lm = new levenbergMarquardt op in

            lm#minimize 

        method value =
            ()

        method updateProcessArguments (params: float array) =
            _process#updateArguments params


        initializer 
            _params <- Array.make 5 0.;
            _params.(0) <- _process#theta;
            _params.(1) <- _process#kappa;
            _params.(2) <- _process#sigma;
            _params.(3) <- _process#rho;
            _params.(4) <- _process#v0


    end;;
