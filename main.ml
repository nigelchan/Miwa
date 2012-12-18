open Utils

open Mw_math
open Mw_stat
open Mw_optimize

open Model
open Process
open Payoff

open CalibratedModel
open AnalyticHeston
open Engine
open McEngine;;


let strikeArray = [| 3400.;3600.;3800.;4000.;4200.;4400.;4500.;4600.;4800.;5000.;5200.;5400.;5600. |] in
let expirydayArray = [| 13.; 41.; 75.; 165.; 256.; 345.; 524.; 703. |] in
let expiryArray = (Array.map (fun e -> e /. 365. ) expirydayArray) in
let volSurface =
     [| [| 0.6625;0.4875;0.4204;0.3667;0.3431;0.3267;0.3121;0.3121 |];
        [| 0.6007;0.4543;0.3967;0.3511;0.3279;0.3154;0.2984;0.2921 |];
        [| 0.5084;0.4221;0.3718;0.3327;0.3155;0.3027;0.2919;0.2889 |];
        [| 0.4541;0.3869;0.3492;0.3149;0.2963;0.2926;0.2819;0.2800 |];
        [| 0.4060;0.3607;0.3330;0.2999;0.2887;0.2811;0.2751;0.2775 |];
        [| 0.3726;0.3396;0.3108;0.2781;0.2788;0.2722;0.2661;0.2686 |];
        [| 0.3550;0.3277;0.3012;0.2781;0.2781;0.2661;0.2661;0.2681 |];
        [| 0.3428;0.3209;0.2958;0.2740;0.2688;0.2627;0.2580;0.2620 |];
        [| 0.3302;0.3062;0.2799;0.2631;0.2573;0.2533;0.2504;0.2544 |];
        [| 0.3343;0.2959;0.2705;0.2540;0.2504;0.2464;0.2448;0.2462 |];
        [| 0.3460;0.2845;0.2624;0.2463;0.2425;0.2385;0.2373;0.2422 |];
        [| 0.3857;0.2860;0.2578;0.2399;0.2357;0.2327;0.2312;0.2351 |];
        [| 0.3976;0.2860;0.2607;0.2356;0.2297;0.2268;0.2241;0.2320 |] |] in

let s0 = 4468.17 in
let rfRate = 0.036 in

(******************************
 *                            *
 * Phase I: Model Calibration *
 *                            *
 ******************************)

printf "[DATA] S: %0.2f" s0;
print_newline();

printf "[DATA] Implied Volatility Surface:";
print_newline();

printf "%s" "K\\ T"; 
Array.iter (fun x -> printf "%15.4f" x) expiryArray;
print_newline();

printf "============================================================================================================================"; 
print_newline();

Array.iteri (fun i volCurve ->
    printf "%4.0f" strikeArray.(i);
    Array.iter (fun vol -> printf "%15.4f" vol) volCurve;
    print_newline();
) volSurface;
print_newline();




(* Generate instruments Array *)
let instrumentsArray = ref [||] in 

Array.iteri (fun k volCurve ->
    Array.iteri (fun t impVol ->
        let temp = new  hestonCalibrateInstrument impVol s0 strikeArray.(k) rfRate expiryArray.(t) in
        instrumentsArray := Array.append !instrumentsArray [|temp|]
    ) volCurve 
) volSurface;

(* Initial Values *)
let v0 = 0.1 in
let kappa = 1.0 in
let theta = 0.1 in
let sigma = 0.5 in
let rho = -0.5 in


let process = new hestonProcess s0 rfRate v0 rho kappa theta sigma in 
let model = new hestonCalibratedModel process in
let hs_analytic_e = new analyticHestonEuropeanEngine (model :> hestonProcess model) in 

Array.iter (fun (instr:hestonCalibrateInstrument) -> instr#setEngine (hs_analytic_e :> engine)) !instrumentsArray;

model#calibrate !instrumentsArray;




(****************************
 *                          *
 * Phase II: Option Pricing *
 *                          *
 ****************************)

print_newline();
print_newline();
printf "Option pricing (K = 4200; T = 0.4583)...";
print_newline();
print_newline();


let cal_strike = 4200. in 
let cal_expiry = 0.4583 in 

let db_payoff = new europeanPayoff Call cal_strike in 
let arg = {optionType = Call; strike = cal_strike; expiry = cal_expiry} in

let db_hsi = new europeanOption db_payoff cal_expiry in


(*********************
 *                   *
 * Crude Monte Carlo *
 *                   *
 *********************)

printf "Using Crude Monte Carlo Engine...";
print_newline();

let bs = new blackScholesMertonProcess s0 rfRate 0.3 in 
let bs_model = new blackScholeMertonModel bs in

let bs_e = new mcEuropeanEngine ~arguments:arg bs_model db_payoff cal_expiry 1 10000 in (* FIXME *)

db_hsi#setEngine bs_e;

printf "Call price = %f" db_hsi#getPrice;
print_newline();
print_newline();


(************************
 *                      *
 * Heston Semi Analytic *
 *                      *
 ************************)

printf "Using Heston Semi Analytic Engine...";
print_newline();

db_hsi#setEngine hs_analytic_e;
hs_analytic_e#setupArguments arg;

printf "Call price = %f" db_hsi#getPrice;
print_newline();
print_newline();


(**********************
 *                    *
 * Heston Monte Carlo *
 *                    *
 **********************)

printf "Using Heston Monte Carlo Engine...";
print_newline();

let hs_montecarlo_e = new mcEuropeanEngine ~arguments:arg (model:> hestonProcess model) db_payoff cal_expiry 50 1000 in (* FIXME *)
db_hsi#setEngine hs_montecarlo_e;

printf "Call price = %f" db_hsi#getPrice;
print_newline();
print_newline();





(*
(* let hp = new hestonProcess s0 r (sigma *. sigma) 1. 1. (sigma *. sigma) 0.001 in  
let hs_analytic_e = new analyticHestonEuropeanEngine hp arg in
*)
let pg = new pathsGenerator bs in
let pp = new pathPricer test_payoff (1. /. (pow (1. +. r) 1) ) in
let a_path = pg#generateOnePath in
List.iter (fun x -> print_float x; print_newline()) (List.rev a_path);
print_float (pp#getValue a_path);;



let a_path = pg#generateOnePath in
let accumulator = new statistics in
List.map accumulator#push a_path;
print_float accumulator#average;;

print_float (bs#evolve 23000. 0.0001 rnd);;


print_float(icdf#standard_value 0.90);;
print_float(icdf#standard_value 0.95);;
print_float(icdf#standard_value 0.975);;
*)


                    (*
temp = dmax1(epsfcn,MACHEP);
eps = std::sqrt(temp);
ij = 0; 
for( j=0; j<n; j++ )
    {    
    temp = x[j];
    h = eps * std::fabs(temp);
    if(h == zero)
        h = eps; 
    x[j] = temp + h; 
    fcn(m,n,x,wa,iflag);
    if( *iflag < 0) 
        return;
    x[j] = temp;
    for( i=0; i<m; i++ )
        {    
        fjac[ij] = (wa[i] - fvec[i])/h;
        ij += 1;    /* fjac[i+m*j] */
        }    
    } 


let a = x.{0} in
let lambda = x.{1} in

(* Jacobian matrix J(i,j) = dfi / dxj, *)
(* where fi = (Yi - yi)/sigma[i],      *)
(*       Yi = A * exp(-lambda * i) + b  *)
(* and the xj are the parameters (A,lambda,b) *)

let e = exp (-. lambda *. float i) in
let s = sigma.(i) in
j.{i, 0} <- e /. s ;
j.{i, 1} <- float (-i) *. a *. e /. s ;
j.{i, 2} <- 1. /. s
done
*)

(*
 Array xt(n);
        std::copy(x, x+n, xt.begin());
        // constraint handling needs some improvement in the future:
        // starting point should not be close to a constraint violation
            const Array& tmp = currentProblem_->values(xt);
            std::copy(tmp.begin(), tmp.end(), fvec);
        } else {
            std::copy(initCostValues_.begin(), initCostValues_.end(), fvec);
        }   

                let a = x.{0} in
                let lambda = x.{1} in
                let b = x.{2} in

                

                    let yi = a *. exp (-. lambda *. (float i)) +. b in

                    f.{i} <- (yi -. y.(i)) /. sigma.(i)
                done
        (*
        (new levenbergMarquardt)#minimize (data ()) [| 1.0;  0.;  0.; |]

        *)


            (*
             *
let solv (y, sigma) xinit = 
  assert(Array.length sigma = n) ;
  let print_state = print_state n p in
  Printf.printf "\nsolver: %s\n" (Multifit_nlin.name s) ;
  print_state 0 s ; 

  let pos = Vector.create 3 in
  Multifit_nlin.position s pos ;
  let covar = Matrix.create p p in
  Multifit_nlin.covar s ~epsrel:0. covar ;
  Printf.printf
    "A      = %.5f +/- %.5f\n" pos.{0} (sqrt covar.{0, 0}) ;
  Printf.printf
    "lambda = %.5f +/- %.5f\n" pos.{1} (sqrt covar.{1, 1}) ;
  Printf.printf
    "b      = %.5f +/- %.5f\n" pos.{2} (sqrt covar.{2, 2}) 

let _ = 
             *
             *
             *
             * *)
*)
