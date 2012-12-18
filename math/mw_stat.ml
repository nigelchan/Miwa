class statistics =
    object (self)
        val mutable _container = ([] : float list)
        method pushElement x =
            _container <- x :: _container

        method pushList x =
            _container <- List.append x _container;

        method average =
            (List.fold_left ( +. ) 0. _container) /. (float_of_int (List.length _container))
    end;;
