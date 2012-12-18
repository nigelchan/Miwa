open Complex

open Utils
open Mw_math
open Mw_integral
open Process
open Engine

open Model;;

let integrator = new gaussianLaguerreIntegration 50;; (* FIXME *)

class analyticHestonEuropeanEngine ?(arguments: vanillaArguments option) (model:hestonProcess model) =
    object (self)
        inherit engine ?arguments model


        method calculate =

            let _s0     = (model#process)#s0 in
            let _rfRate = (model#process)#rfRate in
            let _v0     = (model#process)#v0  in
            let _rho    = (model#process)#rho  in
            let _kappa  = (model#process)#kappa  in
            let _theta  = (model#process)#theta  in
            let _sigma  = (model#process)#sigma  in

            (*
            printf "_v0    = %f\n" _v0   ;
            printf "_rho   = %f\n" _rho     ;
            printf "_kappa = %f\n" _kappa   ;
            printf "_theta = %f\n" _theta   ;
            printf "_sigma = %f\n" _sigma   ;
            *)


            let (expiry, strike, optionType) = match _arguments with
            |   None -> raise (MW_FAIL "No argument is provided")
            |   Some a -> a.expiry, a.strike, a.optionType in
            

            let mulc c x =
                {re = c *. x.re ; im = c *. x.im}
            in

            let discount = (1. /. ((1. +. _rfRate) ** expiry) ) in (* Should be replaced by forward rate object and its methods *)
            
            let big_p j phi = 

                let logStock = Pervasives.log _s0 in
                let logStrike = Pervasives.log strike in
                let ratio = discount in  (* should be discount / dividendFac *)
                let dd = logStock -. Pervasives.log ratio in
                let sigmaSQ = _sigma *. _sigma in
                let rsigma = _rho *. _sigma in
                let t0 = _kappa -. (if (j==0) then (_rho *. _sigma) else 0.) in

                let rpsig = rsigma *. phi in

                let t1 = {re = t0; im = -. rpsig} in

                let d = sqrt (sub (mul t1 t1) (mul (mulc (sigmaSQ *. phi) one)
                                                   {re = -.phi; im = if (j==0) then 1. else -1.}) ) in

                let ex = exp (neg (mulc expiry d)) in

                let ret = ref 0. in

                if _sigma > epsilon_float then begin

                    let p = div (sub t1 d) (add t1 d) in
                    let g = log (div (sub one (mul p ex))
                                     (sub one p)) in

                    let big_C = mulc ((_kappa *. _theta) /. sigmaSQ)
                                     (sub (mulc expiry (sub t1 d))
                                          (mulc 2.0 g)) in

                    let big_D = div (mul (sub t1 d) (sub one ex))
                                    (mulc sigmaSQ
                                          (sub one (mul ex p))) in



                    let iux =  {re = 0.; im = phi *. (dd -. logStrike)} in

                    ret := (exp (add (add (mulc _v0 big_D) big_C) iux )).im /. phi;
                              
                end else begin
                    let td = mul (div (mulc phi one) (mulc 2.0 t1))
                                 {re = -.phi; im = if (j==0) then 1. else -1.} in

                    let p = mul td (div (mulc sigmaSQ one) (add t1 d)) in
                    let g = mul p (sub one ex) in

                    let big_D = div (mul td (sub one ex))
                                    (sub one (mul p ex)) in

                    let big_C = mulc (_kappa *. _theta) 
                                     (sub (mulc expiry td) (mulc (2.0 /. sigmaSQ) g)) in

                    let iux =  {re = 0.; im = phi *. (dd -. logStrike)} in

                    ret := (exp ( add (add (mulc _v0 big_D) big_C) iux )).im /. phi

                end;
                !ret
            in



            let big_p_0 = (integrator#integrate (big_p 0)) /. pi in
            let big_p_1 = (integrator#integrate (big_p 1)) /. pi in

            match optionType with
                Call -> _s0 *. (big_p_0 +. 0.5) -. strike *. discount *. (big_p_1 +. 0.5)
            |   Put -> _s0 *. (big_p_0 -. 0.5) -. strike *. discount *. (big_p_1 -. 0.5)
    end;;

(*
let s = 36. in 
let sigma = 0.2 in
let r = 0.06 in 
let arg = {optionType = Put; strike = 40.; expiry = 1.} in

let hp = new hestonProcess s r (sigma *. sigma) 0. 1. (sigma *. sigma) 0.001 in 

let eng = new analyticHestonEuropeanEngine hp arg in

print_float(eng#calculate);
*)
