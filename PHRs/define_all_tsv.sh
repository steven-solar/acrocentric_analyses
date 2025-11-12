#!/bin/bash


for species in mSymSyn1; do #mPanTro3 mPanPan1 mGorGor1 mPonAbe1 mPonPyg2 HG002
	echo $species
	mkdir -p $species
	ref=/data/Phillippy2/projects/acro_comparisons/refs/$species/ref.fa
	# seqtk telo -d 10000 $ref | cut -f1-3 > $species.telo.bed
	# python extend_subterms.py $species.telo.bed $ref.fai 200000 > $species.telo.ext.bed
	awk '$4 ~/subTerm/' /data/Phillippy2/projects/acro_comparisons/acrocentric_analyses/genome_features/censat_final/${species}_merged.bed | bedtools merge -d 500000 -i - > $species.subterm_exclude.bed
	python extend_subterms.py $species.subterm_exclude.bed $ref.fai 500000 > $species.subterm_exclude.ext.bed

	telo=$species.telo.ext.bed
	subterm=$species.subterm_exclude.ext.bed

	# cat /data/Phillippy2/projects/acro_comparisons/acrocentric_analyses/genome_features/censat_final/${species}_merged.bed \
	# 	| cut -f1-3 \
	# 	| bedtools sort -g $ref.fai -i - \
	# 	| bedtools merge -i - \
	# 	> $species.censat.bed	
	cat $telo $subterm $species.censat.bed | bedtools sort -g $ref.fai -i - | bedtools merge -i - > $species.censat_telo_subterm.bed

	echo "done with making exclusion beds"

	while read -r bed; do
		echo $bed
		filename=$(basename "$bed")
		chrom1=$(echo "$filename" | cut -d'.' -f1)
		chrom2=$(echo "$filename" | cut -d'.' -f2)
		ln -sf $bed $species/$chrom1.$chrom2.bed
		bedtools subtract -a $bed -b $species.censat.bed > $species/$chrom1.$chrom2.exclude_censat.bed
		bedtools subtract -a $bed -b $species.censat_telo_subterm.bed > $species/$chrom1.$chrom2.exclude_censat_telo_subterm.bed
	done < <(ls /data/Phillippy2/projects/acro_comparisons/acrocentric_analyses/PHR_analysis/mashmap_heatmaps/$species/*.*.bed | grep -v -E "EBV|MT|chrM" | grep -v exclude)

	echo "done w filtering and linking beds"

	rm -f $species.PHRs.exclude_censat_telo_subterm.tsv
	rm -f $species.PHRs.exclude_censat.tsv
	rm -f $species.PHRs.unfiltered.tsv

	for bed in $species/*.*.*bed; do
		echo $bed

		if [[ $bed == *"exclude_censat_telo_subterm"* ]]; then
			status="exclude_censat_telo_subterm"
		elif [[ $bed == *"exclude_censat"* ]]; then
			status="exclude_censat"
		else
			status="unfiltered"
		fi
		echo $status
		filename=$(basename "$bed")
		base=${filename%%.exclude_censat_telo_subterm.bed}
		base=${base%%.exclude_censat.bed}
		base=${base%%.bed}
		chrom1=${base%%.*}
		chrom2=${base#*.}

		python get_bp_pid_from_bed.py $bed | grep -E -v "morph|chrM|MT|EBV|random" >> $species.PHRs.$status.tsv
	done
done

for species in mSymSyn1; do #HG002 mPanTro3 mPanPan1 mGorGor1 mPonAbe1 mPonPyg2 
	for status in unfiltered exclude_censat exclude_censat_telo_subterm; do
		python squish_dups.py $species.PHRs.$status.tsv > $species.PHRs.$status.squished.tsv
	done
done

for status in unfiltered exclude_censat exclude_censat_telo_subterm; do
	rm -f $status.tsv
	cat *.PHRs.$status.squished.tsv > $status.tsv
done
