#!/bin/bash

fa_1=$1
fa_2=$2
s=$3
out=$4

echo "wfmash -t 16 -s $s -p 85 $fa_2 $fa_1 > $out"

wfmash \
	-t 16 \
	-s $s \
	-p 85 \
	$fa_2 $fa_1 \
	> $out

query=$(head -n1 $out | cut -f1)
q_len=$(head -n1 $out | cut -f2)
ref=$(head -n1 $out | cut -f6)
r_len=$(head -n1 $out | cut -f7)

echo "$query	$q_len	0	1	+	$ref	$r_len	0	1	1	1	1	gi:f:1	bi:f:1	md:f:1	cg:z:1=" >> $out
echo "$query	$q_len	$((q_len-1))	$q_len	+	$ref	$r_len	$((r_len-1))	$r_len	1	1	1	gi:f:1	bi:f:1	md:f:1	cg:z:1=" >> $out
