#!/bin/bash

mkdir -p out

for species in mSymSyn1; do #mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 HG002
	echo $species
	for status in exclude_censat_telo_subterm; do
		echo $status
		cat $species/*.*.$status.bed \
			| awk '$5 >= 990 {print}' \
			| awk '{split($1, a, "_"); split($4, b, "_"); if (a[1] != b[1]) {print}}' \
			| bedtools sort -i - -g /data/Phillippy2/projects/acro_comparisons/refs/$species/ref.fa.fai \
			| bedtools merge -i - \
			> $species.PHRs.$status.bed
	done
	echo unfiltered
	cat $(ls $species/*.*.bed | grep -v exclude) \
		| awk '$5 >= 990 {print}' \
		| awk '{split($1, a, "_"); split($4, b, "_"); if (a[1] != b[1]) {print}}' \
		| bedtools sort -i - -g /data/Phillippy2/projects/acro_comparisons/refs/$species/ref.fa.fai \
		| bedtools merge -i - \
		> $species.PHRs.unfiltered.bed
done

for species in mSymSyn1; do #mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 HG002
	for status in unfiltered exclude_censat exclude_censat_telo_subterm; do
		python phr_relationship_graph.py $status.tsv $species out/$species.$status.svg
	done
done

for species in mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 HG002; do
	awk '{sum+=$3-$2} END {print sum}' $species.PHRs.unfiltered.bed
done

for species in mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 HG002; do
	awk '{sum+=$3-$2} END {print sum}' $species.PHRs.exclude_censat_telo_subterm.bed
done

for species in mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 HG002; do
	bedtools subtract -a $species.PHRs.unfiltered.bed -b /data/Phillippy2/projects/acro_comparisons/acrocentric_analyses/genome_features/rdna_arrays_final/$species.rdna_arrays.bed \
		| bedtools intersect -a - -b /data/Phillippy2/projects/acro_comparisons/acrocentric_analyses/genome_features/censat_final/${species}_merged.bed \
		| awk '{sum+=$3-$2} END {print sum}'
done
