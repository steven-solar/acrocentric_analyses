#!/bin/bash

awk '$1 == "mGorGor1" && $2 ~ /hsa9/ && $3 ~ /hsa18/' /data/Phillippy2/projects/acro_comparisons/acrocentric_analyses/PHR_analysis/define_PHRs/all.tsv | awk '{bp+=$4; w_pid+=$4*$5;} END {print bp, w_pid/bp}'
awk '$1 == "mGorGor1" && $2 ~ /hsa9/ && $3 ~ /hsa18/' /data/Phillippy2/projects/acro_comparisons/acrocentric_analyses/PHR_analysis/define_PHRs/exclude_censat_telo_subterm.tsv | awk '{bp+=$4; w_pid+=$4*$5;} END {print bp, w_pid/bp}'

awk '$1 == "mGorGor1" && (($2 ~ /hsa2a/ && ($3 ~ /hsa13/ || $3 ~ /hsa15/)) || ($2 ~ /hsa13/ && $3 ~ /hsa15/))' all.tsv | awk '{bp+=$4; w_pid+=$4*$5;} END {print bp, w_pid/bp}'
awk '$1 == "mGorGor1" && (($2 ~ /hsa2a/ && ($3 ~ /hsa13/ || $3 ~ /hsa15/)) || ($2 ~ /hsa13/ && $3 ~ /hsa15/))' exclude_censat_telo_subterm.tsv | awk '{bp+=$4; w_pid+=$4*$5;} END {print bp, w_pid/bp}'