(* Stochastic Process Abstract Class *)
(* Can be multi-dimensional. Therefore, everything has to be an Array *)
class virtual stochasticProcess (x0:float array) (rfRate:float) =
    object (self)

        val _x0 = x0  (* Inital values -- an Array *)
        val _rfRate = rfRate

        method x0 = _x0 
        method factors = Array.length _x0
        method rfRate = _rfRate

        (* move forward from s_t to s_{t + dt} *) 
        method virtual evolve : float array -> float array -> float -> float array

    end;;

class blackScholesMertonProcess (s0:float) (rfRate:float) (vol:float) =
    object (self)
        inherit stochasticProcess [|s0|] rfRate

        val _vol = vol

        method drift dt =
            ( _rfRate -. 1. /. 2. *. _vol *. _vol ) *. dt

        method diffusion dt dW =
            _vol *. (sqrt dt) *. dW

        method evolve st dW dt =
            [|st.(0) *. ( exp ((self#drift dt) +. (self#diffusion dt dW.(0))) )|]

        method s0 = _x0.(0)
    end;;

(* dS = \mu S dt + \sqrt{v} S dZ_1
 * dv = -\kappa (v - \theta) dt + \sigma \sqrt{v} dZ_2 *)
class hestonProcess (s0:float) (rfRate:float)
                    (v0:float) (rho:float)
                    (kappa:float) (theta:float) (sigma:float) =

    object (self)

        inherit stochasticProcess [|s0; v0|] rfRate

        val mutable _rho    = rho 
        val mutable _kappa  = kappa 
        val mutable _theta  = theta 
        val mutable _sigma  = sigma 

        (*
        method drift dt = 
        method diffusion dt dW = 
        *)

        method updateArguments (params: float array) = 
            _theta <- params.(0);
            _kappa <- params.(1);
            _sigma <- params.(2);
            _rho <- params.(3);
            _x0.(1) <- params.(4)


        method evolve (xt:float array) (dW:float array) (dt:float) =

            let vol = if (xt.(1) > 0.) then (sqrt xt.(1)) else 0. in (* Truncation at 0 *)
            let mu = _rfRate -. 0.5 *. vol *. vol in
            let dWstar = _rho *. dW.(0) +. (sqrt (1. -. _rho *.  _rho)) *. dW.(1) in

            let nu = _kappa *. (_theta -. vol *. vol) in
            
            let ret = Array.create 2 0. in

            ret.(0) <- xt.(0) *. exp (mu *. dt +. vol *. dW.(0)  *. (sqrt dt));
            ret.(1) <- xt.(1) +. nu *. dt +. _sigma *. vol *. (sqrt dt) *. dWstar;

            ret
            
        method rho    = _rho 
        method kappa  = _kappa 
        method theta  = _theta 
        method sigma  = _sigma 
        method s0 = _x0.(0)
        method v0 = _x0.(1)

    end;;

