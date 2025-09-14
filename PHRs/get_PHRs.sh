#!/bin/bash

if [ -z ${SLURM_CPUS_PER_TASK+x} ]; then
    cpus=32
else
    cpus=$SLURM_CPUS_PER_TASK
fi

MAX_JOBS=$cpus

pwd=$(pwd)
species=$1
awk -v OFS='\t' '$5 >= 990 {print}' | bedtools merge -i - -c 4,5,6 -o unique,mean,mode