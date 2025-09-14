#!/bin/bash

ml load mashmap wfmash

submit() {
	log=$1
	fa_1=$2
	fa_2=$3
	s=$4
	out=$5
	jid=$(sbatch --time=01:00:00 --mem=16g --cpus-per-task=4 --partition=norm,quick -o $log wfmash.sh $fa_1 $fa_2 $s $out)
	echo $jid $log wfmash.sh $fa_1 $fa_2 $s $out
}
mkdir -p logs out refs out/individual_pafs

samples=(chm13 mPanTro3 mPanPan1 mGorGor1 mPonAbe1 mPonPyg2)
ancient_acros_human=(chr2 chr2 chr9 chr13 chr14 chr15 chr18 chr21 chr22 chrY)
ancient_acros_primates=(hsa2a hsa2b hsa9 hsa13 hsa14 hsa15 hsa18 hsa21 hsa22 hsaY)

for i in {0..9}; do
	acro_primate=${ancient_acros_primates[$i]}
	acro_human=${ancient_acros_human[$i]}
	echo "Processing $acro_primate $acro_human"
	rm -f refs/$acro_primate.fa.gz
	mkdir -p out/individual_pafs/$acro_primate
	for a in {0..5}; do
		species=${samples[$a]}
		mkdir -p refs/$species
		echo $species
		if [[ $species == "chm13" ]]; then
			ref=refs/$species/$acro_human.fa.gz
			gunzip < refs/$species/$acro_human.fa.gz | sed "s/^>/\>${a}#${species}#/g" | bgzip -@8 > refs/$species/$acro_human.rename.fa.gz
			samtools faidx refs/$species/$acro_human.rename.fa.gz
			rename_ref=refs/$species/$acro_human.rename.fa.gz
		else
			ref=refs/$species/*$acro_primate.fa.gz
			chrom=$(basename $ref | cut -d. -f1)
			if [[ $species == "mPanPan1" || $species == "mPanTro3" ]]; then
				if [[ $acro_human == "chr9" ]]; then
					gunzip < $ref | sed "s/^>/\>${a}#${species}#/g" | sed "s/hsa9/hsa9#rev/g" | seqtk seq -r - | bgzip -@8 > refs/$species/$chrom.rename.fa.gz
					samtools faidx refs/$species/$chrom.rename.fa.gz
				elif [[ $acro_human == "chrY" ]]; then
					gunzip < $ref | sed "s/^>/\>${a}#${species}#/g" | sed "s/hsaY/hsaY#rev/g" | seqtk seq -r - | bgzip -@8 > refs/$species/$chrom.rename.fa.gz
					samtools faidx refs/$species/$chrom.rename.fa.gz
				else
					gunzip < $ref | sed "s/^>/\>${a}#${species}#/g" | bgzip -@8 > refs/$species/$chrom.rename.fa.gz
					samtools faidx refs/$species/$chrom.rename.fa.gz
				fi
			else
				gunzip < $ref | sed "s/^>/\>${a}#${species}#/g" | bgzip -@8 > refs/$species/$chrom.rename.fa.gz
				samtools faidx refs/$species/$chrom.rename.fa.gz
			fi
			rename_ref=refs/$species/$chrom.rename.fa.gz
		fi
		echo merge
		cat $rename_ref >> refs/$acro_primate.fa.gz
		samtools faidx refs/$acro_primate.fa.gz
	done
	for a in {0..5}; do
		for b in {0..5}; do
			species_1=${samples[$a]}
			species_2=${samples[$b]}
			if [[ $species_1 == $species_2 ]]; then
				continue
			fi
			echo "Comparing $species_1 and $species_2 $acro_human"
			if [[ $species_1 == "chm13" ]]; then
				fa_1=refs/$species_1/$acro_human.rename.fa.gz
			else
				fa_1=$(ls refs/$species_1/*$acro_primate.rename.fa.gz)
			fi
			if [[ $species_2 == "chm13" ]]; then
				fa_2=refs/$species_2/$acro_human.rename.fa.gz
			else
				fa_2=$(ls refs/$species_2/*$acro_primate.rename.fa.gz)
			fi
			if [ ! -f $fa_1.fai ]; then samtools faidx $fa_1; fi
			if [ ! -f $fa_2.fai ]; then samtools faidx $fa_2; fi
			while [ `jobs -p | wc -l` -ge 32 ]; do 
				sleep 5
			done
			if [[ $acro_human == "chr9" || $acro_human == "chrY" ]]; then
				if [[ $species_1 == "mPanPan1" || $species_1 == "mPanTro3" || $species_2 == "mPanPan1" || $species_2 == "mPanTro3" ]]; then
					submit logs/$acro_human.$acro_primate.$species_1.$species_2.wfmash_50kb.out $fa_1 $fa_2 50000 out/individual_pafs/$acro_primate/$acro_human.$acro_primate.$species_1.$species_2.wfmash_50kb.paf &
				fi
			fi
		done
	done
done

wait `jobs -p`
