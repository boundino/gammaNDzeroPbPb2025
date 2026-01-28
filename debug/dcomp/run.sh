#!/bin/bash

VARS=(
    # Dmass
    # Dalpha
    # Dchi2cl
    # Ddtheta
    # Dtrk1Pt
    # Dtrk2Pt
    dls
)

CONFIGS=(
    configs/23_gammaN_Cut23PAS.conf
    configs/25_gammaN_Cut23PAS.conf
    configs/25noccf_gammaN_Cut23PAS.conf
    configs/23_Ngamma_Cut23PAS.conf
    configs/25_Ngamma_Cut23PAS.conf
    configs/25noccf_Ngamma_Cut23PAS.conf
)

[[ $# -eq 0 || ${1:-0} -eq 1 ]] && { make savehist || exit 1 ; }
[[ $# -eq 0 || ${2:-0} -eq 1 ]] && { make drawhist || exit 1 ; }

for var in ${VARS[@]} ; do
    [[ ${1:-0} -eq 1 ]] && {
        for conf in ${CONFIGS[@]} ; do
            echo $conf $var
            ./savehist $conf $var &
        done
        wait
    }

    [[ ${2:-0} -eq 1 ]] && {
        ./drawhist configs/23_gammaN_Cut23PAS.conf,configs/25noccf_gammaN_Cut23PAS.conf,configs/25_gammaN_Cut23PAS.conf $var
        ./drawhist configs/23_Ngamma_Cut23PAS.conf,configs/25noccf_Ngamma_Cut23PAS.conf,configs/25_Ngamma_Cut23PAS.conf $var
    }
done
