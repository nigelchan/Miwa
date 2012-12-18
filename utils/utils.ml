open Printf


(* EXCEPTION *)
exception MW_FAIL of string


(* TYPES *)
type optionType = Call | Put;;
type path = float list;;
type multiPath = float array list;;

type vanillaArguments = {optionType: optionType; strike: float; expiry: float};;  (* expiry should be (last) date instead *)



(* FUNCTIONS *) 
let printf = Printf.printf;; 
let print_array a = Array.map (fun aelt -> printf "%f" aelt; print_newline()) a;;

let rec replicate n f = 
    let rec do_once n accum =
        (*
        if (List.length accum) mod 100 == 0 then begin
            printf "%u" (List.length accum);
            print_newline();
        end;
        *)
        if List.length accum == n then accum
        else do_once n (f() :: accum)
    in 
    do_once n []
;;
