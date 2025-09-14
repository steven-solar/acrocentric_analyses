#!/bin/bash

echo "all"
python summarize_phrs.py all.tsv
echo ""

echo "exclude_censat"
python summarize_phrs.py exclude_censat.tsv
echo ""

echo "exclude_censat_telo_subterm"
python summarize_phrs.py exclude_censat_telo_subterm.tsv
