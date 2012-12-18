open Process

class virtual ['a] model (process: 'a) =
    object (self)

        val _process = process

        val mutable _params = ([||] : float array)

        method process =
            _process

        method setParams (params: float array) =
            _params <- params;
            self#updateProcessArguments params

        method virtual updateProcessArguments: float array -> unit

    end;;


class blackScholeMertonModel (process: blackScholesMertonProcess) =
    object (self)
        inherit ['a] model process
        method updateProcessArguments (params:float array) =
            ()
    end;;
