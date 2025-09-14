#!/bin/bash

ml load R

plot() {
	paf=$1
	out=$2
	method=$3
	bed=$4
	echo "Rscript SVbyEye.R $paf $out $method $bed"
	Rscript SVbyEye.R $paf $out $method $bed
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

# for paf in ../out/individual_pafs/*/*.wfmash_50kb.paf; do
# 	bn=$(basename $paf)
# 	name="${bn%.paf}"
# 	hsa=$(basename $paf| cut -d. -f2)
# 	mkdir -p plots/single/$hsa
# 	echo $name
# 	while [ `jobs -p | wc -l` -ge 16 ]
# 		do
# 			sleep 5
# 		done
# 	# awk '{print }
# 	plot $paf plots/single/$hsa/$name.pdf plotMiro sat.tsv &> logs/$name.out &
# done

pafs=(/data/Phillippy2/projects/acro_comparisons/synteny_plots/out/individual_pafs/hsa2a/chr2.hsa2a.mPanPan1.mPonAbe1.wfmash_50kb.paf
	/data/Phillippy2/projects/acro_comparisons/synteny_plots/out/individual_pafs/hsa2b/chr2.hsa2b.mGorGor1.mPonAbe1.wfmash_50kb.paf
	/data/Phillippy2/projects/acro_comparisons/synteny_plots/out/individual_pafs/hsa14/chr14.hsa14.mGorGor1.mPonAbe1.wfmash_50kb.paf
	/data/Phillippy2/projects/acro_comparisons/synteny_plots/out/individual_pafs/hsa15/chr15.hsa15.mPanTro3.mPonAbe1.wfmash_50kb.paf
	/data/Phillippy2/projects/acro_comparisons/synteny_plots/out/individual_pafs/hsa18/chr18.hsa18.chm13.mPonAbe1.wfmash_50kb.paf)

for paf in ${pafs[@]}; do
	bn=$(basename $paf)
	name="${bn%.paf}"
	hsa=$(basename $paf| cut -d. -f2)
	mkdir -p plots/single/$hsa
	echo $name
	while [ `jobs -p | wc -l` -ge 16 ]
		do
			sleep 5
		done
	plot $paf plots/single/$hsa/$name.pdf plotMiro sat.tsv &> logs/$name.out &
done

wait `jobs -p`

# plot /data/Phillippy2/projects/acro_comparisons/synteny_plots/out/individual_pafs/hsa18/chr18.hsa18.chm13.mPanPan1.wfmash_50kb.paf plots/single/hsa18/chr18.hsa18.chm13.mPanPan1.wfmash_50kb.png plotMiro sats.4col.bed
