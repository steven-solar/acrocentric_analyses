#!/bin/bash

ml load R

plot() {
	paf=$1
	out=$2
	method=$3
	strand_bed=$4
	unstranded_bed=$5
	echo "Rscript SVbyEye.R $paf $out $method $strand_bed $unstranded_bed"
	Rscript SVbyEye.R $paf $out $method $strand_bed $unstranded_bed
}

# for paf in ../out/individual_pafs/*.minimap.paf; do
# 	bn=$(basename $paf)
# 	name="${bn%.paf}"
# 	hsa=$(basename $paf| cut -d. -f2)
# 	mkdir -p plots/single/$hsa
# 	while [ `jobs -p | wc -l` -ge 32 ]
# 		do
# 			sleep 5
# 		done
# 	plot $paf plots/single/$hsa/$name.png plotMiro &
# done

# for paf in ../out/all_vs_all/*.wfmash_ava.paf; do
# 	bn=$(basename $paf)
# 	name="${bn%.paf}"
# 	hsa=$(basename $paf| cut -d. -f1)
# 	method=$(basename $paf| cut -d. -f2)
# 	if [[ $hsa == "hsa2a" || $hsa == "hsa2b" ]]; then
# 		continue
# 	fi
# 	mkdir -p plots/ava/$method
# 	while [ `jobs -p | wc -l` -ge 32 ]
# 		do
# 			sleep 5
# 		done
# 	# awk '{print }
# 	plot $paf plots/ava/$method/$hsa.png plotAVA &
# done

samples=(chm13 mPanTro3 mPanPan1 mGorGor1 mPonAbe1 mPonPyg2)

# grep -v track /data/Phillippy2/projects/acro_comparisons/ape_autosomes_figure/panelA/tracks/sat/chm13.bed | awk -v OFS='\t' -v sample=chm13 -v i=0 '{print i"."sample"."$0}' > sats.bed

# for i in {1..5}; do
# 	sample=${samples[$i]}
# 	short_id=${sample:1:-1}
# 	grep -v track /data/Phillippy2/projects/acro_comparisons/ape_autosomes_figure/panelA/tracks/sat/${short_id}_rm.cenSatv1.2.bed | awk -v OFS='\t' -v sample=$sample -v i=$i '{print i"#"sample"#"$0}' >> sats.bed
# done

# for paf in ../out/merged_pafs/*.wfmash_50kb_merged.paf; do
# 	bn=$(basename $paf)
# 	name="${bn%.paf}"
# 	hsa=$(basename $paf| cut -d. -f2)
# 	method=$(basename $paf| cut -d. -f3)
# 	# echo Rscript SVbyEye.R $paf plots/merged/$method/$hsa.$method.allsats.pdf $method $bed
# 	mkdir -p plots/merged/$method
# 	while [ `jobs -p | wc -l` -ge 32 ]
# 		do
# 			sleep 5
# 		done
# 	# plot $paf plots/merged/final/$hsa.pdf plotAVA sat.final.stranded.tsv sat.final.unstranded.tsv &> $hsa.out &
# 	plot $paf plots/merged/final/$hsa.pdf plotAVA sat.tsv &> $hsa.out &
# done

for paf in ../out/merged_pafs/*.wfmash_50kb_merged.filt.paf; do
	bn=$(basename $paf)
	name="${bn%.paf}"
	hsa=$(basename $paf| cut -d. -f2)
	method=$(basename $paf| cut -d. -f3)
	mkdir -p plots/merged/$method
	while [ `jobs -p | wc -l` -ge 32 ]
		do
			sleep 5
		done
	# plot $paf plots/merged/final/$hsa.pdf plotAVA sat.final.stranded.tsv sat.final.unstranded.tsv &> $hsa.out &
	plot $paf plots/merged/final/$hsa.filt.pdf plotAVA sat.tsv &> $hsa.out &
done

wait `jobs -p`
