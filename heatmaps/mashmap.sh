#!/bin/bash

if [ -z ${SLURM_CPUS_PER_TASK+x} ]; then
    cpus=32
else
    cpus=$SLURM_CPUS_PER_TASK
fi

sample=$1
fa_1=$2
fa_2=$3
out=$4
chr_1=$(basename $fa_1 | cut -d. -f1)
chr_2=$(basename $fa_2 | cut -d. -f1)
echo "$chr_1 to $chr_2"

mashmap --threads $cpus -r $fa_2 -q $fa_1 -t $cpus -M --pi 80 -s 10000 -o $out/$chr_1.$chr_2.mashmap
python get_bed.py $out/$chr_1.$chr_2.mashmap > $out/$chr_1.$chr_2.bed
