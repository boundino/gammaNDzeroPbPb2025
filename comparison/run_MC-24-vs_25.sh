#!/bin/bash

INPUTS=(
    # D mesons
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2024PbPb/Dzero_260714-gen_HiForest_260328_prompt_GNucleusToD0-PhotonBeamA_Bin-Pthat0_Kpi_trkpt0p1_Drej-genmatched_Dpt-2_Dsize.root;2024 MC (#gammaN);mc24-BeamA"
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPbMC/Dzero_260714-gen_HiForest_260904_prompt_GNucleusToD0-BeamA_SoftQCD_KPi_2025_trkpt0p1_Drej-genmatched_Dpt-2_Dsize.root;2025 MC (#gammaN);mc25-BeamA"

    # Gen D mesons & Event variables
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2024PbPb/Dzero_260714-gen_HiForest_260328_prompt_GNucleusToD0-PhotonBeamA_Bin-Pthat0_Kpi_trkpt0p1_Drej-genmatched_Dpt-2.root;2024 MC (#gammaN);mc24-BeamA"
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPbMC/Dzero_260714-gen_HiForest_260904_prompt_GNucleusToD0-BeamA_SoftQCD_KPi_2025_trkpt0p1_Drej-genmatched_Dpt-2.root;2025 MC (#gammaN);mc25-BeamA"
)

VARS=(
    # Dmva_BDT
    Dmass
    Dalpha
    Ddls
    Ddls_2D
    Dtrk1Pt
    Dtrk2Pt
    Dchi2cl
    Dtrk1ptrel
    Dtrk2ptrel
    Dpt
    Dy

    # ZDCsumPlus
    # ZDCsumMinus
    # nTrackInAcceptanceHP
    # HFEMaxPlusforest
    # HFEMaxMinusforest
    # HFEMaxPlusforest-zoom
    # HFEMaxMinusforest-zoom
    # nVtx

    # Gpt
    # Gy
    # Gtk1pt
    # Gtk2pt
    # Gtk1eta
    # Gtk2eta
)

CUTEVTS=(
    "1;No event selections;noevtcut"
)
##

CUTDS=(
    "1;;" # for event variables 
    "Dgen==23333;Signal D#scale[0.6]{#lower[-0.7]{0}} #rightarrow K#pi;-Dgen"
    # "Dgen==23333 && TMath::Abs(Dtrk1PtErr/Dtrk1Pt)<0.1 && TMath::Abs(Dtrk2PtErr/Dtrk2Pt)<0.1 && TMath::Abs(Dtrk1Eta) < 2.4 && TMath::Abs(Dtrk2Eta) < 2.4 && Dtrk1Pt > 0.5 && Dtrk2Pt > 0.5 && Dchi2cl > 0.05 && (DsvpvDistance/DsvpvDisErr) > 1. && DsvpvDisErr>1.e-8 && DsvpvDisErr_2D>1.e-8;Signal D#scale[0.6]{#lower[-0.7]{0}} #rightarrow K#pi && Precuts;-Dgenprecut"
    "GisSignalCalc;Gen signal;-Gsignal"
    # "Dpt>0;;-Dnocut"
    # "Dtrk1PtErr/Dtrk1Pt<0.1 && Dtrk2PtErr/Dtrk2Pt<0.1 && DpassCut23PAS && (Dtrk1PixelHit+Dtrk1StripHit)>=11 && (Dtrk2PixelHit+Dtrk2StripHit)>=11;DpassCut23PAS;-D23pas"
)

make savehist.exe calchists.exe drawhists.exe || exit 1

for cutdstr in "${CUTDS[@]}" ; do
    IFS=';' ; cutdtags=($cutdstr) ; unset IFS ; cutd=${cutdtags[0]} ; cutd_tex=${cutdtags[1]} ; cutd_tag=${cutdtags[2]} ; 

    # event cut
    for cutevtstr in "${CUTEVTS[@]}" ; do
        IFS=';' ; cutevttags=($cutevtstr) ; unset IFS ; cutevt=${cutevttags[0]} ; cutevt_tex=${cutevttags[1]} ; cutevt_tag=${cutevttags[2]} ; cutevt_year=${cutevttags[3]}

        cut_tag=$cutevt_tag$cutd_tag
        cutstr=$cutevt" && "$cutd";"$cutevt_tex"%%"$cutd_tex";"$cut_tag
        
        for var in "${VARS[@]}" ; do
            [[ ($var == D* && ${cutd} == D*) || ($var == G* && ${cutd} == G*) || (($var != D* && $var != G*) && ${cutd} == 1) ]] || { continue ; }

            echo -e "\033[33m"$var" \033[33;2m("$cut_tag")\033[0m \033[2m"$cutstr"\033[0m"

            compare_list=''
            tag_list=''
            for inputstr in "${INPUTS[@]}" ; do
                IFS=';' ; inputtags=($inputstr) ; unset IFS ; input_tag=${inputtags[2]} ; 
                [[ x$input_tag == x ]] && { echo "warning: missed input_tag. skip." ; continue ; }

                echo -e "    \033[33m"$var" \033[33;2m("$input_tag")\033[0m"

                itag="rootfiles/"$cut_tag"/"$var"_"$input_tag"_savehist" # 
                echo "    "$itag

                [[ ${1:-0} -eq 1 ]] && {
                    ./savehist.exe "$inputstr" "$cutstr" "$var" $itag &
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
        wait

    done
done
