#!/bin/bash

# chrom_num() {
# 	chrom=$1
# 	species=$2
# 	if [[ $species == "mSymSyn1" || $species == "HG002" ]]; then
# 		chrom_name=$(echo $chrom | cut -d'_' -f1 | sed 's/chr//')
# 	else
# 		chrom_name=$(echo $chrom | cut -d'_' -f3 | sed 's/hsa//')
# 	fi
# 	echo $chrom_name
# }

# get_pair_relationship() {
# 	mashmap_bed=$1
# 	species=$2
# 	ref_chrom=$(basename $mashmap_bed | cut -d. -f1)
# 	tgt_chrom=$(basename $mashmap_bed | cut -d. -f2)
# 	ref_chrom_num=$(chrom_num $ref_chrom $species)
# 	tgt_chrom_num=$(chrom_num $tgt_chrom $species)
# 	if [[ $ref_chrom_num == $tgt_chrom_num ]]; then
# 		return
# 	fi
# 	awk -v OFS='\t' -v ref_chrom=$ref_chrom -v tgt_chrom=$tgt_chrom '{ 
# 		if ($5 >= 990) {
# 			tot_bp+=($3-$2); 
# 			tot_pid+=($5 * ($3-$2));
# 		} 
# 		} END {
# 			if (tot_bp > 0) {
# 				print ref_chrom, tgt_chrom, tot_bp, tot_pid/tot_bp/10;
# 			} else {
# 				print ref_chrom, tgt_chrom, 0, 0;
# 			}
# 		}' $mashmap_bed >> $species/$species.with_censat.tsv
# 	echo $species $ref_chrom $tgt_chrom done
# }

# get_pair_relationship_v2() {
# 	mashmap_bed=$1
# 	species=$2
# 	out=$3
# 	ref_chrom=$(basename $mashmap_bed | cut -d. -f1)
# 	tgt_chrom=$(basename $mashmap_bed | cut -d. -f2)
# 	awk -v OFS='\t' -v ref_chrom=$ref_chrom -v tgt_chrom=$tgt_chrom '{ 
# 		if ($5 >= 990) {
# 			tot_bp+=($3-$2); 
# 			tot_pid+=($5 * ($3-$2));
# 		} 
# 		} END {
# 			if (tot_bp > 0) {
# 				print ref_chrom, tgt_chrom, tot_bp, tot_pid/tot_bp/10;
# 			} else {
# 				print ref_chrom, tgt_chrom, 0, 0;
# 			}
# 		}' $mashmap_bed >> $out
# 	echo $species $ref_chrom $tgt_chrom done
# }

# phr_graph() {
# 	species=$1
# 	python phr_relationship_graph_v2.py all.with_censat.tsv $species $species/$species.relationships.with_censat.png > $species/$species.with_censat.norm.tsv
# }

phr_graph_v2() {
	species=$1
	python phr_relationship_graph_v2.py all.tsv $species out/$species.relationships.all.svg # > $species/$species.all.norm.tsv
	python phr_relationship_graph_v2.py exclude_censat.tsv $species out/$species.relationships.exclude_censat.svg # > $species/$species.exclude_censat.norm.tsv
	python phr_relationship_graph_v2.py exclude_censat_telo_subterm.tsv $species out/$species.relationships.exclude_censat_telo_subterm.svg # > $species/$species.exclude_censat_telo_subterm.norm.tsv
}

# mkdir -p out
# rm -f all_species.all.tsv all_species.exclude_censat.tsv all_species.exclude_censat_telo_subterm.tsv
# rm -f out/*.png

# for species in mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 mSymSyn1 chm13 HG002; do
# 	mkdir -p $species
# 	rm -f $species/$species.all.tsv $species/$species.exclude_censat.tsv $species/$species.exclude_censat_telo_subterm.tsv
# 	ref=../../refs/$species/ref.fa.gz
# 	echo $species

# 	seqtk telo -d 10000 $ref | cut -f1-3 > $species/$species.telo.bed
# 	python extend_subterms.py $species/$species.telo.bed $ref.fai 200000 > $species/$species.telo.ext.bed
# 	awk '$4 ~/subTerm/' ../../genome_features/censat/$species.censat.bed | bedtools merge -d 500000 -i - > $species/$species.subterm_exclude.bed
# 	python extend_subterms.py $species/$species.subterm_exclude.bed $ref.fai 500000 > $species/$species.subterm_exclude.ext.bed
# 	cut -f1-3 ../../genome_features/censat/$species.censat.bed > $species/$species.censat.3col.bed
# 	bedtools sort -g $ref.fai -i $species/$species.censat.3col.bed | bedtools merge -i - > $species/$species.exclude_censat.bed
# 	cat $species/$species.exclude_censat.bed $species/$species.telo.ext.bed $species/$species.subterm_exclude.ext.bed | bedtools sort -g $ref.fai -i - | bedtools merge -i - > $species/$species.exclude_censat_telo_subterm.bed

# 	mashmap_dir=../mashmap_heatmaps/$species
# 	while read -r ref_chrom; do
# 		while read -r mashmap_bed; do
# 			while [ `jobs -p | wc -l` -gt 31 ]; do
# 				sleep 1
# 			done

# 			tgt_chrom=$(basename $mashmap_bed | cut -d. -f2)
# 			ref_chrom_num=$(chrom_num $ref_chrom $species)
# 			tgt_chrom_num=$(chrom_num $tgt_chrom $species)
# 			if [[ $ref_chrom_num == $tgt_chrom_num ]]; then
# 				continue
# 			fi

# 			bedtools subtract -a $mashmap_bed -b $species/$species.exclude_censat.bed > $mashmap_dir/$ref_chrom.$tgt_chrom.exclude_censat.bed
# 			bedtools subtract -a $mashmap_bed -b $species/$species.exclude_censat_telo_subterm.bed > $mashmap_dir/$ref_chrom.$tgt_chrom.exclude_censat_telo_subterm.bed

# 			get_pair_relationship_v2 $mashmap_bed $species $species/$species.all.tsv &
# 			get_pair_relationship_v2 $mashmap_dir/$ref_chrom.$tgt_chrom.exclude_censat.bed $species $species/$species.exclude_censat.tsv &
# 			get_pair_relationship_v2 $mashmap_dir/$ref_chrom.$tgt_chrom.exclude_censat_telo_subterm.bed $species $species/$species.exclude_censat_telo_subterm.tsv &
# 		done < <(ls $mashmap_dir/$ref_chrom.*.bed | grep -E -v "random|chrM|MT|EBV|utig|morph")
# 	done < <(cut -f1 $ref.fai | grep -E -v "random|chrM|MT|EBV|utig|morph")

# 	wait `jobs -p` 

# 	awk -v OFS='\t' -v species=$species '{print species, $0}' $species/$species.all.tsv >> all_species.all.tsv
# 	awk -v OFS='\t' -v species=$species '{print species, $0}' $species/$species.exclude_censat.tsv >> all_species.exclude_censat.tsv
# 	awk -v OFS='\t' -v species=$species '{print species, $0}' $species/$species.exclude_censat_telo_subterm.tsv >> all_species.exclude_censat_telo_subterm.tsv
# done

for species in mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 mSymSyn1 chm13 HG002; do
	echo $species
	phr_graph_v2 $species &
done

wait `jobs -p`
