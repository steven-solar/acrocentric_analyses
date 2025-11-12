#!/bin/bash

for species in mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 mSymSyn1 HG002; do
	for status in unfiltered exclude_censat exclude_censat_telo_subterm; do
		# echo -e "$species\t$status\t$(grep $species $status.tsv | awk '{sum+=$4} END {print sum}')"
		python phr_relationship_graph.py $status.tsv $species out/$species.$status.png
	done
done