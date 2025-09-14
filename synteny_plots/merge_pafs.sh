#!/bin/bash

samples=(chm13 mPanTro3 mPanPan1 mGorGor1 mPonAbe1 mPonPyg2)
ancient_acros_human=(chr2 chr2 chr9 chr13 chr14 chr15 chr18 chr21 chr22 chrY)
ancient_acros_primates=(hsa2a hsa2b hsa9 hsa13 hsa14 hsa15 hsa18 hsa21 hsa22 hsaY)

for i in {0..9}; do
	acro_primate=${ancient_acros_primates[$i]}
	acro_human=${ancient_acros_human[$i]}
	echo "Processing $acro_primate $acro_human"
	cat out/individual_pafs/$acro_primate/added.paf | sed 's/#/\./g' > out/merged_pafs/$acro_human.$acro_primate.wfmash_50kb_merged.paf
	for a in {0..5}; do
		for b in {0..5}; do
			species_1=${samples[$a]}
			species_2=${samples[$b]}
			if [[ $species_1 == $species_2 ]]; then
				continue
			fi
			# awk -v OFS='\t' -v ref_species=$species_1 -v query_species=$species_2 '{print query_species"#"$1,$2,$3,$4,$5,ref_species"#"$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22}' out/individual_pafs/$acro_human.$acro_primate.$species_2.$species_1.minimap.paf >> out/merged_pafs/$acro_human.$acro_primate.minimap.paf
			# cat out/individual_pafs/$acro_primate/$acro_human.$acro_primate.$species_1.$species_2.minimap.paf | sed 's/#/\./g' >> out/merged_pafs/$acro_human.$acro_primate.minimap_merged.paf
			# cat out/individual_pafs/$acro_primate/$acro_human.$acro_primate.$species_1.$species_2.wfmash_10kb.paf | sed 's/#/\./g' >> out/merged_pafs/$acro_human.$acro_primate.wfmash_10kb_merged.paf
			cat out/individual_pafs/$acro_primate/$acro_human.$acro_primate.$species_1.$species_2.wfmash_50kb.paf | sed 's/#/\./g' >> out/merged_pafs/$acro_human.$acro_primate.wfmash_50kb_merged.paf
			# cat out/individual_pafs/$acro_primate/$acro_human.$acro_primate.$species_1.$species_2.wfmash_25kb.paf | sed 's/#/\./g' >> out/merged_pafs/$acro_human.$acro_primate.wfmash_25kb_merged.paf
			# cat out/individual_pafs/$acro_human.$acro_primate.$species_1.$species_2.wfmash.paf >> out/merged_pafs/$acro_human.$acro_primate.wfmash.paf
		done
	done
done
