#!/bin/bash

INPUTS=(
    # Event variables
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260212-hfle_2023PbPbUPC_Jan2024ReReco_20260212Forest_HIForward0_Dpt-2_Trig-2.root,/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260212-hfle_2023PbPbUPC_Jan2024ReReco_20260212Forest_HIForward1_Dpt-2_Trig-2.root,/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260212-hfle_2023PbPbUPC_Jan2024ReReco_20260212Forest_HIForward2_Dpt-2_Trig-2.root;2023 (Jan24 Reco);d23-rJan24;"
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260426-yrefmva_2023PbPbUPC_Feb2025ReReco_20260521Forest_HIForward0_Dpt-2_Trig-2.root,/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260426-yrefmva_2023PbPbUPC_Feb2025ReReco_20260521Forest_HIForward1_Dpt-2_Trig-2.root,/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260426-yrefmva_2023PbPbUPC_Feb2025ReReco_20260521Forest_HIForward2_Dpt-2_Trig-2.root;2023 (Feb25 Reco);d23-rFeb25;"
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260426-yrefmva_PbPbUPC_HIForward0_Dpt-2.root;2025;d25-rp;"

    # D mesons
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260212-hfle_2023PbPbUPC_Jan2024ReReco_20260212Forest_HIForward_Dpt-2_Trig-2_Dsize_xbr.root;2023 (Jan24 Reco);d23-rJan24"
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260426-yrefmva_2023PbPbUPC_Feb2025ReReco_20260521Forest_HIForward_Dpt-2_Trig-2_Dsize.root;2023 (Feb25 Reco);d23-rFeb25"
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260426-yrefmva_PbPbUPC_HIForward_Dpt-2_Dsize_12ePD.root;2025;d25-rp"

    # Empty BX
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260212-hfle_HiForest_260218_HIEmptyBX_HIRun2023A_PromptReco_v2.root;2023 EmptyBX;d23-rp"
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260212-hfle_HiForest_260218_HIEmptyBX_HIRun2025A_PromptReco_v1.root;2025 EmptyBX;d25-rp"
)

VARS=(
    # nTrackInAcceptanceHP
    # ZDCsumPlus
    # ZDCsumMinus
    # HFEMaxPlusforest
    # HFEMaxMinusforest
    # HFEMaxPlusforest-zoom
    # HFEMaxMinusforest-zoom
    nVtx
    
    #
    # Dmass
    # Dalpha
    # Ddls
    # Dalpha-Low
    # Dtrk1Pt
    # Dtrk2Pt
    # Dchi2cl
    # Dtrk1ptrel
    # Dtrk2ptrel
    # Dtrk1nhit
    # Dtrk2nhit
)

CUTEVTS=(
    # "isL1ZDCOr && ZDCgammaN && HFEMaxPlus_forest < 16 && cscTightHalo2015Filter && selectedVtxFilter;ZDC Xn0n (#gammaN), gap, |v_{z}| < 15, cscHaloTight%%clusComp for 2023;gammaN;2025"
    # "isL1ZDCOr && ZDCgammaN && gapgammaN && selectedBkgFilter && selectedVtxFilter;ZDC Xn0n (#gammaN), gap, |v_{z}| < 15, cscHaloTight%%clusComp for 2023;gammaN;2023"
    # "isL1ZDCOr && ZDCgammaN && gapgammaN && cscTightHalo2015Filter && selectedVtxFilter;ZDC Xn0n (#gammaN), gap, |v_{z}| < 15, cscHaloTight;gammaN-noccf;2023"

    "isL1ZDCOr && ZDCgammaN && selectedBkgFilter && selectedVtxFilter;ZDC Xn0n (#gammaN), |v_{z}| < 15, cscHaloTight%%clusComp for 2023;gammaN-nogap;2023"
    "isL1ZDCOr && ZDCgammaN && cscTightHalo2015Filter && selectedVtxFilter;ZDC Xn0n (#gammaN), |v_{z}| < 15, cscHaloTight%%clusComp for 2023;gammaN-nogap;2025"

    # "isL1ZDCOr && ZDCgammaN && gapgammaN && selectedBkgFilter && selectedVtxFilter && nVtx <= 3;ZDC Xn0n (#gammaN), gap, |v_{z}| < 15, cscHaloTight%%clusComp, nVtx <= 3 for 2023;gammaN-nvtx3;2023"

    # "isL1ZDCOr;HLT_ZDCOr && PV filter && ZDC0nOr (1500 GeV);nocut;"
    # "isNotBptxOR;HLT_HIL1NotBptxOR;isNotBptxOR;"
)
##

CUTDS=(
    "1;;" # for event variables 
    "Dpt>2 && Dpt<5 && fabs(Dy)<2;;-Dpre"
    "Dpt>2 && Dpt<5 && fabs(Dy)<2 && Dtrk1PtErr/Dtrk1Pt<0.1 && Dtrk2PtErr/Dtrk2Pt<0.1 && DpassCut23PAS && (Dtrk1PixelHit+Dtrk1StripHit)>=11 && (Dtrk2PixelHit+Dtrk2StripHit)>=11;DpassCut23PAS;-D23pas"
)

make savehist.exe calchists.exe drawhists.exe || exit 1

for var in "${VARS[@]}" ; do
    for cutdstr in "${CUTDS[@]}" ; do
        IFS=';' ; cutdtags=($cutdstr) ; unset IFS ; cutd=${cutdtags[0]} ; cutd_tex=${cutdtags[1]} ; cutd_tag=${cutdtags[2]} ; 
        [[ ($var == D* && ${cutd} != 1) || ($var != D* && ${cutd} == 1) ]] || { continue ; }

        # event cut
        for cutevtstr in "${CUTEVTS[@]}" ; do
            IFS=';' ; cutevttags=($cutevtstr) ; unset IFS ; cutevt=${cutevttags[0]} ; cutevt_tex=${cutevttags[1]} ; cutevt_tag=${cutevttags[2]} ; cutevt_year=${cutevttags[3]}

            cut_tag=$cutevt_tag$cutd_tag
            cutstr=$cutevt" && "$cutd";"$cutevt_tex"%%"$cutd_tex";"$cut_tag
            
            echo -e "\033[33m"$var" \033[33;2m("$cut_tag")\033[0m \033[2m"$cutstr"\033[0m"

            compare_list=''
            tag_list=''
            for inputstr in "${INPUTS[@]}" ; do
                IFS=';' ; inputtags=($inputstr) ; unset IFS ; input_tag=${inputtags[2]} ; 
                [[ x$input_tag == x ]] && { echo "warning: missed input_tag. skip." ; continue ; }
                [[ ($var == D* && ${inputtags[0]} == *Dsize*) || ($var != D* && ${inputtags[0]} != *Dsize*) ]] || { continue ; }
                [[ ($cutevt_year == 2025 && $input_tag == d23*) || ($cutevt_year == 2023 && $input_tag == d25*) ]] && { continue ; } #

                echo -e "    \033[33m"$var" \033[33;2m("$input_tag")\033[0m"

                itag="rootfiles/"$var"/"$cut_tag"-"$input_tag"_savehist" # 
                echo "    "$itag

                [[ ${1:-0} -eq 1 ]] && {
                    ./savehist.exe "$inputstr" "$cutstr" "$var" $itag
                }

                [[ ${2:-0} -eq 1 ]] && {
                    ./calchists.exe $itag".root"
                }

                itag=${itag/_savehist/_calchist}
                compare_list=$compare_list","$itag".root"
                tag_list=$tag_list"_"$cut_tag"-"$input_tag
                itag=''
            done
            compare_list=${compare_list#,}
            tag_list=${tag_list#_}

            echo $compare_list
            echo $tag_list
            [[ ${3:-0} -eq 1 ]] && {
                ./drawhists.exe "$compare_list" "$tag_list" 0
            }
        done

    done
done
