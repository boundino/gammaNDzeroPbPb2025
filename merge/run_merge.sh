#!/bin/bash

# skim=0 ; label=''
skim=1 ; label='_Dsize-gt-0_hltZDCOr_12ePD'
ntotal=-1
[[ $ntotal -gt 0 ]] && label=${label}'_nf-'$ntotal
inputdirs=(
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260123_PbPbUPC_HIForward*_Drej-pasor.root"
    
    #####################
    # Example input patterns
    # The automatic output name is always the last part of the input 
    #####################
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260123_PbPbUPC_HIForward*_Drej-pasor.root" # Asterisk
    # /eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260123_2023PbPbUPC_Jan2024ReReco_20260201Forest_HIForward0_Drej-pasor # Input directory
    # /export/d00/scratch/jwang/L1PbPb2021/crab_L1_20200322_HIZeroBiasReducedFormat_HIRun2018A_v1_2/ # Input directiory ending with /
    # /eos/cms/store/group/phys_heavyions/wangj/Forest2023/Lcpks/Pythia8_LcToKsPr_prompt_gN-PhotonA_Pthat2-pt0p9_PbPb_5362GeV/crab_HiForest_251030_LcToKsPr_prompt_gN-PhotonA_Pthat2_v0/251031_005700/0000,/eos/cms/store/group/phys_heavyions/wangj/Forest2023/Lcpks/Pythia8_LcToKsPr_prompt_gN-PhotonA_Pthat2-pt0p9_PbPb_5362GeV/crab_HiForest_251030_LcToKsPr_prompt_gN-PhotonA_Pthat2_v0.root # Forced output name
)

########################################
## >>> do not change lines below >>>  ##
########################################

tmp=$(date +%y%m%d%H%M%S)
cp merge.cc merge_${tmp}.cc
make merge_${tmp}.exe || { rm merge_${tmp}.cc ; exit 1 ; }
mkdir -p filelists

for i in "${inputdirs[@]}" ; do

    echo
    echo -e "\e[32m$i\e[0m"

    input="${i%,*}"
    input=${input%/}

    output="${i#*,}"
    [[ "$i" != *,* ]] && {
        inputplain=${input//'*'/''}
        IFS='/'; subdir=($inputplain); unset IFS;
        request=${subdir[${#subdir[@]}-1]} ## last path
        output=${inputplain%${request}*}$request
    }
    inputtag=${output##*/} ; inputtag=${inputtag%.root} ; # inputtag independent of label

    [[ $output != *.root ]] && { output=${output}.root ; }
    output=${output/.root/${label}.root}
    echo $output

    [[ "$input" != *.root ]] && { input=${input}/*.root ; }
    filelist=filelists/files_${inputtag}.txt
    echo $filelist
    rm -f $filelist
    ls $input -d > $filelist 
    [[ -s $filelist ]] || { rm $filelist ; echo -e "\e[31mwarning: no valid input files, skip.\e[0m" ; continue ; }

    [[ ${1:-0} -eq 1 ]] && {
        willrun=1
        [[ -f $output ]] && {
            echo "error: output $output exits. "
            rewrite=
            while [[ $rewrite != 'y' && $rewrite != 'n' ]] ; do
                echo "remove output? (y/n):"            
                read rewrite
                if [[ $rewrite == 'y' ]] ; then { rm $output ; echo "$output removed" ; } ;
                elif [[ $rewrite == 'n' ]] ; then { echo "please change output file name" ; willrun=0 ; } ;
                fi
            done
        }

        [[ $willrun -eq 0 ]] && continue
        ./merge_${tmp}.exe $output $filelist $skim $ntotal;
    }
done

echo
rm -v merge_${tmp}.*
