#!/bin/bash

config=$1

[[ ${2:-0} -eq 1 || $# -eq 0 ]] && { make savehist.exe || exit 1 ; }
[[ ${3:-0} -eq 1 || $# -eq 0 ]] && { make drawhist.exe || exit 1 ; }

[[ $config != *.conf ]] && {
    echo "usage: ./run.sh [conf] ([savehist] [drawhist] [varlist])"
    exit 1
}

t_variables='_forest,,_pt0p1,_eta5'
[[ $# -ge 4 ]] && t_variables=$4
IFS=',' ; v_variables=($t_variables) ; unset IFS;

[[ ${2:-0} -eq 1 ]] && {
    ./savehist.exe $config
    for v in ${v_variables[@]} ; do
        echo $v
        ./savehist.exe $config $v &
    done
    wait
}

[[ ${3:-0} -eq 1 ]] && {
    echo $t_variables
    ./drawhist.exe $config $t_variables
}
