#!/bin/bash

samples=(chm13 mPanTro3 mPanPan1 mGorGor1 mPonAbe1 mPonPyg2 mSymSyn1)

grep -v track /data/Phillippy2/projects/acro_comparisons/synteny_plots/SVbyEye/test/extdata/censat/chm13v2.0.cenSatv2.1.bed | awk -v OFS='\t' -v sample=chm13 -v i=0 '{print i"."sample"."$0}' > sats.bed
grep -v track /data/Phillippy2/projects/acro_comparisons/synteny_plots/SVbyEye/test/extdata/censat/chm13v2.0.SatelliteStrandv2.1.bed | awk -v OFS='\t' -v sample=chm13 -v i=0 '{print i"."sample"."$0}' > sats.strand.bed
awk -v OFS='\t' -v sample=chm13 -v i=0 '{print i"."sample"."$0}' /data/Phillippy2/projects/acro_comparisons/ape_autosomes_figure/chr_p_q_lens/cens/chm13.cens.bed > cens.bed
awk -v OFS='\t' -v sample=chm13 -v i=0 '{print i"."sample"."$0}' $CHM13.fai > all.fai


for i in {1..6}; do
	sample=${samples[$i]}
	echo $sample
	censat=/data/Phillippy2/projects/acro_comparisons/synteny_plots/SVbyEye/test/extdata/censat/$sample.cenSatv2.0.bed
	censat_strand=/data/Phillippy2/projects/acro_comparisons/synteny_plots/SVbyEye/test/extdata/censat/$sample.cenSatStrandv2.0.bed
	cens=/data/Phillippy2/projects/acro_comparisons/ape_autosomes_figure/chr_p_q_lens/cens/$sample.cens.bed
	fai=/data/Phillippy2/projects/acro_comparisons/refs/$sample/ref.fa.fai
	if [[ $sample == "mPanTro3" || $sample == "mPanPan1" ]]; then
		grep -v track $censat | grep -E "hsa9|hsaY" | python invert_bed.py $fai | awk -v OFS='\t' -v sample=$sample -v i=$i '{print i"."sample"."$1".rev",$2,$3,$4,$5,$6,$7,$8,$9}' >> sats.bed
		grep -v track $censat | grep -E -v "hsa9|hsaY" | awk -v OFS='\t' -v sample=$sample -v i=$i '{print i"."sample"."$0}' >> sats.bed
		grep -v track $censat_strand | grep -E "hsa9|hsaY" | python invert_bed.py $fai | awk -v OFS='\t' -v sample=$sample -v i=$i '{print i"."sample"."$1".rev",$2,$3,$4,$5,$6,$7,$8,$9}' >> sats.strand.bed
		grep -v track $censat_strand | grep -E -v "hsa9|hsaY" | awk -v OFS='\t' -v sample=$sample -v i=$i '{print i"."sample"."$0}' >> sats.strand.bed
		grep -E "hsa9|hsaY" $cens | python invert_bed.py $fai | awk -v OFS='\t' -v sample=$sample -v i=$i '{print i"."sample"."$1".rev",$2,$3}' >> cens.bed
		grep -E -v "hsa9|hsaY" $cens | awk -v OFS='\t' -v sample=$sample -v i=$i '{print i"."sample"."$0}' >> cens.bed
		cat $fai | grep -E "hsa9|hsaY" | awk -v OFS='\t' -v sample=$sample -v i=$i '{print i"."sample"."$1".rev",$2,$3}' >> all.fai
		cat $fai | grep -E -v "hsa9|hsaY" | awk -v OFS='\t' -v sample=$sample -v i=$i '{print i"."sample"."$0}' >> all.fai
	else
		grep -v track $censat | awk -v OFS='\t' -v sample=$sample -v i=$i '{print i"."sample"."$0}' >> sats.bed
		grep -v track $censat_strand | awk -v OFS='\t' -v sample=$sample -v i=$i '{print i"."sample"."$0}' >> sats.strand.bed
		awk -v OFS='\t' -v sample=$sample -v i=$i '{print i"."sample"."$0}' $cens >> cens.bed
		awk -v OFS='\t' -v sample=$sample -v i=$i '{print i"."sample"."$0}' $fai >> all.fai
	fi
done

echo done

echo -e "seqnames\tstart\tend\tname\tID\tstrand\tcolor" > sat.tsv
cat sats.bed | grep -v -w ct | python filter_bed.py | awk -v OFS='\t' '{print $1,$2,$3,$4,$1,$6,$9}' >> sat.tsv

echo -e "seqnames\tstart\tend\tname\tID\tstrand\tcolor" > sat.strand.tsv
cat sats.strand.bed | grep -v -w ct | python filter_bed.py | awk -v OFS='\t' '{print $1,$2,$3,$4,$1,$6,$9}' >> sat.strand.tsv

grep -i active_hor sats.bed > cens.bed
bedtools intersect -a sats.strand.bed -b cens.bed | grep -i AS_strand | grep -w -v -E "ct|gap" | python filter_bed.py > sats.all.bed
bedtools subtract -A -a sats.bed -b cens.bed | grep -i -v active_hor | grep -w -v -E "ct|gap" | python filter_bed.py >> sats.all.bed
bedtools sort -i sats.all.bed -g all.fai > sats.final.bed
echo -e "seqnames\tstart\tend\tname\tID\tstrand\tcolor" > sat.final.tsv
awk -v OFS='\t' '{print $1,$2,$3,$4,$1,$6,$9}' sats.final.bed >> sat.final.tsv

echo -e "seqnames\tstart\tend\tname\tID\tstrand\tcolor" > sat.final.stranded.tsv
awk -v OFS='\t' '$6 == "+" || $6 == "-" {print $1,$2,$3,$4,$1,$6,$9}' sats.final.bed >> sat.final.stranded.tsv
echo -e "seqnames\tstart\tend\tname\tID\tstrand\tcolor" > sat.final.unstranded.tsv
awk -v OFS='\t' '$6 == "." {print $1,$2,$3,$4,$1,$6,$9}' sats.final.bed >> sat.final.unstranded.tsv


grep chm13 sat.final.tsv | sed 's/0\.chm13\.//g' | awk -v OFS='\t' 'NR>1 {print $1,$2,$3,$4,'0',$6,$2,$3,$7}' > chm13.final.bed
grep mPonAbe1 sat.final.tsv | sed 's/0\.chm13\.//g' | awk -v OFS='\t' 'NR>1 {print $1,$2,$3,$4,'0',$6,$2,$3,$7}' > mPonAbe1.final.bed
