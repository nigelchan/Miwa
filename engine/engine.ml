open Utils
open Model

class virtual engine ?(arguments: vanillaArguments option) (model:'a model) = (* Should be general model *)
    object (self)

        val _model = model
        val _process = model#process 
        val mutable _arguments = arguments

        method virtual calculate: float 
        method setupArguments (arguments:vanillaArguments) =
            _arguments <- Some arguments
    end;;
