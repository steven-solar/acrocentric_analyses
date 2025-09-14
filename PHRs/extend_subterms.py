import sys

subterms = open(sys.argv[1])
fai = open(sys.argv[2])
added_len = int(sys.argv[3])

chrom_lens = dict()
for line in fai:
	line_split = line.strip().split('\t')
	chrom, length = line_split[0], int(line_split[1])
	chrom_lens[chrom] = length

for line in subterms:
	line_split = line.strip().split('\t')
	chrom, start, end = line_split[0], int(line_split[1]), int(line_split[2])
	if start != 0:
		start=max(0, start-added_len)
	if end != chrom_lens[chrom]:
		end=min(chrom_lens[chrom], end+added_len)
	print(chrom, start, end, sep='\t')
