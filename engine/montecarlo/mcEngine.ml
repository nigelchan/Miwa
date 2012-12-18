open Utils
open Mw_math
open Mw_stat
open Engine
open Payoff
open Process
open McTools

open Model

class virtual mcEngine ?(arguments:vanillaArguments option) 
                       (model:'a model)
                       (payoff:europeanPayoff)
                       (expiry:float)
                       (timeSteps:int)
                       (nSamples:int) =

    object (self)
        inherit engine ?arguments model


        val _pathsGenerator = new pathsGenerator (model#process :> stochasticProcess) expiry timeSteps
        val _pathPricer = new pathPricer payoff (1. /. ((1. +. (model#process)#rfRate) ** expiry) ) (* Should be replaced by forward rate object and its methods *)
        val _expiry = expiry 
        val _timeSteps = timeSteps
        val _nSamples = nSamples

        val mutable _mean = 0.
        val mutable _errorEst = 0.

        val mutable _accumulator = new statistics

        method virtual calculate: float
    end;;

class mcEuropeanEngine ?(arguments: vanillaArguments option) 
                       (model:'a model)
                       (payoff:europeanPayoff)
                       (expiry:float)
                       (timeSteps:int)
                       (nSamples:int) =
    object (self)
        inherit mcEngine ?arguments model payoff expiry timeSteps nSamples

        method calculate (* nSamples *) =
            let paths = replicate _nSamples (fun () -> _pathsGenerator#generateOnePath) in 
            let price = List.map _pathPricer#getValue paths in
            _accumulator#pushList price;
            _accumulator#average

        (*
        method getMean = 
        method getErrorEstimate =
        *)
    end;;
