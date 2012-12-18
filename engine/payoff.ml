open Utils

class europeanPayoff optionType strike =
    object (self)

    val _strike = strike
    val _optionType = optionType

    method price path =
        (* only requires first element of the path element S_t because this is a single asset instrument *)

        match _optionType with
            Call -> max ((List.hd path).(0) -. _strike) 0.0 
        |   Put  -> max (_strike -. (List.hd path).(0)) 0.0

    end;;
