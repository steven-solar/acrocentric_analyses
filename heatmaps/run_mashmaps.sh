#!/bin/bash

if [ -z ${SLURM_CPUS_PER_TASK+x} ]; then
    cpus=32
else
    cpus=$SLURM_CPUS_PER_TASK
fi

MAX_JOBS=$cpus
USER_ID=solarsj

mkdir -p logs
ml load mashmap ucsc

run_job() {
	sample=$1
	fa_1=$2
	fa_2=$3
	chr_1=$(basename $fa_1 | cut -d. -f1)
	chr_2=$(basename $fa_2 | cut -d. -f1)
	jid=$(sbatch --time=2-00:00:00 --partition=norm --mem=121g --cpus-per-task=32 -o logs/$sample/$sample.$chr_1.$chr_2.out mashmap.sh $sample $fa_1 $fa_2 $sample)
	echo "jobid: $jid sample: $sample chr_1: $chr_1 chr_2: $chr_2"
}

run_to_bed() {
	sample=$1
	chr_1=$2
	chr_2=$3
	python get_bed.py $sample/$chr_1.$chr_2.mashmap > $sample/$chr_1.$chr_2.bed
}

make_bigbed() {
	sample=$1
	chrom=$2
	cat $sample/*.$chrom.bed | grep -E -v "random|chrM|morph|utig|MT|EBV" > $sample/$chrom.bed
	bedToBigBed -tab -type=bed9 $sample/$chrom.bed $ref_dir/$sample/ref.fa.fai seqIdy/$sample/$chrom.bb
}

ref_dir=$1

while read -r sample; do 
	echo "---$sample---"
	mkdir -p $sample $sample logs/$sample
	for fa_1 in $ref_dir/$sample/*.fa; do
		for fa_2 in $ref_dir/$sample/*.fa; do
			if [ "$fa_1" != "$fa_2" ]; then
				if [[ $fa_1 != *"random"* && $fa_2 != *"random"* && $fa_1 != *"MT"* && $fa_2 != *"MT"* && $fa_1 != *"morph"* && $fa_2 != *"morph"* ]]; then
					chr_1=$(basename $fa_1 | cut -d. -f1)
					chr_2=$(basename $fa_2 | cut -d. -f1)
					echo "$chr_1 to $chr_2"
					while [ `jobs -p | wc -l` -ge ${MAX_JOBS} ]
					do
						sleep 5
					done
					while [ $(squeue -u $USER_ID | grep mashmap | wc -l) -ge 250 ]
					do
						sleep 300
					done
					run_job $sample $fa_1 $fa_2 &
				fi
			fi
		done
	done
done < samples.txt

for sample in chm13 HG002 mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 mSymSyn1; do
	mkdir -p seqIdy/$sample
	echo -----$sample-----
	while read -r chrom1; do
		while read -r chrom2; do
			if [[ $chrom1 == $chrom2 ]]; then
				continue
			fi
			echo $chrom1 to $chrom2
			while [ `jobs -p | wc -l` -ge ${MAX_JOBS} ]; do
				sleep 2
			done
			run_to_bed $sample $chrom1 $chrom2 &
		done < <(cut -f1 $ref_dir/$sample/ref.fa.fai | grep -E -v "morph|utig|random|chrM|chrMT|EBV|MT")
	done < <(cut -f1 $ref_dir/$sample/ref.fa.fai | grep -E -v "morph|utig|random|chrM|chrMT|EBV|MT")
done

wait `jobs -p`

for sample in chm13 HG002 mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 mSymSyn1; do
	echo -----$sample-----
	while read chrom; do
		make_bigbed $sample $chrom &
	done < <(cut -f1 $ref_dir/$sample/ref.fa.fai | grep -E -v "morph|utig|random|chrM|chrMT|EBV|MT")
done

wait `jobs -p`
