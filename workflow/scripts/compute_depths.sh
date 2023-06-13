#!/usr/bin/env bash
set -euxo pipefail

mkdir -p depths

populations=(ada_grc ana_tur arm_arm car_aut car_svn_hrv cau_tur cec_grc cyp_cyp ibe_esp_eus ibe_esp_north ibe_esp_south ibe_esp_west lig_ita mac_bgr mac_mkd_grc mac_rou_mda mel_che mel_dnk mel_imn mel_irl mel_rus rut_mlt)
for population in "${populations[@]}"; do
    (samtools merge -u - results/map/filt/"${population}".*.{1..16}.cram \
    | samtools depth -a - \
    | cut -f 3 \
    | sort -n \
    | uniq -c \
    | awk '{print $2"\t"$1}' \
    > depths/"$population".depths.tsv) &
done; wait

