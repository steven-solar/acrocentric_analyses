#!/bin/bash

minimap2 -c -x asm20 --eqx --MD -t 2 45S.fa $species.fa > $species.to.45S.paf
rustybam invert $species.to.45S.paf > 45S.to.$species.paf
awk -v OFS='\t' '{if($11>=10000 && $10/$11 >= 0.95) print $6,$8,$9, $10/$11}' 45S.to.$species.paf > 45S.to.$species.filt.bed

minimap2 -c -x asm20 --eqx --MD -t 2 DJ_palindrome.fa $species.fa > $species.to.DJ_palindrome.paf
rustybam invert $species.to.DJ_palindrome.paf > DJ_palindrome.to.$species.paf
awk -v OFS='\t' '{if($11>=100000 && $10/$11 >= 0.9) print $6,$8,$9, $10/$11}' DJ_palindrome.to.$species.paf > DJ_palindrome.to.$species.filt.bed
