#!/usr/bin/env python3
import sys
from collections import defaultdict
import re

# Usage: python combine_haplotypes_avg.py input.tsv > output.tsv

def simplify(chrom):
    """
    Simplify chromosome names:
    - If it has _hsaXX, return chrXX (keep hsa number only)
    - Else strip haplotype suffixes (_MATERNAL, _PATERNAL, _mat, _pat, _hap1, etc.)
    """
    m = re.search(r'_hsa(\d+\w*|X|Y)$', chrom)
    if m:
        chrom=m.group(1)
        chrom = re.sub(r'chr', '', chrom)
        return f"{m.group(1)}"
    # otherwise, strip haplotype suffixes
    chrom = re.sub(r'_(MATERNAL|PATERNAL|hap1|hap2)$', '', chrom)
    chrom = re.sub(r'^(chr|hsa)_', '', chrom)
    return chrom

groups = defaultdict(lambda: [0, 0, 0])  # (sample, (chr1, chr2)) -> [count, total_bp, weighted_pid_sum]

with open(sys.argv[1]) as f:
	for line in f:
		if not line.strip():
			continue
		sample, chr1, chr2, bp, pid = line.strip().split()
		bp = float(bp)
		pid = float(pid)

		s1 = simplify(chr1)
		s2 = simplify(chr2)
		if s1 == s2:
			continue
		key = tuple(sorted([s1, s2]))

		rec = groups[(sample, key)]
		rec[0] += 1
		rec[1] += bp
		rec[2] += bp * pid  # for weighted PID

# print(groups)
# print average bp and weighted PID
for (sample, (c1, c2)), (n, total_bp, weighted_pid_sum) in groups.items():
	# print(sample, c1, c2, n, total_bp, weighted_pid_sum)
	if total_bp > 0:
		avg_bp = total_bp / n
		avg_pid = weighted_pid_sum / total_bp
	else:
		avg_bp = 0
		avg_pid = 0
	print(f"{sample}\t{c1}\t{c2}\t{avg_bp:.0f}\t{avg_pid:.4f}")
