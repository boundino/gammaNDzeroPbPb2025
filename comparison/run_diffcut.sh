#!/bin/bash

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
    # # Dalpha-Low
    # Dtrk1Pt
    # Dtrk2Pt
    # Dchi2cl
    # Dtrk1ptrel
    # Dtrk2ptrel
    # Dtrk1nhit
    # Dtrk2nhit
)

CUTDS=(
    "1;;"
    "Dpt>2 && Dpt<5 && fabs(Dy)<2;;-Dpre"
    "Dpt>2 && Dpt<5 && fabs(Dy)<2 && Dtrk1PtErr/Dtrk1Pt<0.1 && Dtrk2PtErr/Dtrk2Pt<0.1 && DpassCut23PAS && (Dtrk1PixelHit+Dtrk1StripHit)>=11 && (Dtrk2PixelHit+Dtrk2StripHit)>=11;DpassCut23PAS;-D23pas"
)

make drawhists.exe || exit 1

for var in "${VARS[@]}" ; do
    for cutdstr in "${CUTDS[@]}" ; do
        IFS=';' ; cutdtags=($cutdstr) ; unset IFS ; cutd=${cutdtags[0]} ; cutd_tex=${cutdtags[1]} ; cutd_tag=${cutdtags[2]} ; 
        [[ ($var == D* && ${cutd} != 1) || ($var != D* && ${cutd} == 1) ]] || { continue ; }

        draw_list=(
            # gammaN${cutd_tag}-d23-rJan24,gammaN${cutd_tag}-d23-rFeb25,gammaN${cutd_tag}-d25-rp";"0
            # gammaN${cutd_tag}-d23-rJan24,gammaN${cutd_tag}-d23-rFeb25";"0
            # gammaN${cutd_tag}-d23-rJan24,gammaN${cutd_tag}-d25-rp";"0

            gammaN-nogap${cutd_tag}-d23-rFeb25,gammaN-nogap${cutd_tag}-d25-rp";"0 # for HFEmax
            # gammaN-noccf${cutd_tag}-d23-rJan24,gammaN-noccf${cutd_tag}-d23-rFeb25,gammaN${cutd_tag}-d25-rp";"0
            # gammaN${cutd_tag}-d23-rJan24,gammaN${cutd_tag}-d23-rFeb25,gammaN${cutd_tag}-d25-rp";"0
        )
        for items in "${draw_list[@]}" ; do
            IFS=';' ; draw_opts=($items) ; unset IFS ; do_save_png=${draw_opts[1]}

            compare_list=
            tag_list=
            IFS=',' ; draw_tags=(${draw_opts[0]}) ; unset IFS ;
            for itag in "${draw_tags[@]}" ; do
                compare_list=$compare_list",rootfiles/"$var"/"$itag"_calchist.root"
                tag_list=$tag_list"_"$itag
            done
            compare_list=${compare_list#,}
            tag_list=${tag_list#_}
            echo $compare_list
            echo $tag_list
            [[ ${1:-0} -eq 1 ]] && {
                ./drawhists.exe "$compare_list" "$tag_list" $do_save_png
            }
        done
    done
done
