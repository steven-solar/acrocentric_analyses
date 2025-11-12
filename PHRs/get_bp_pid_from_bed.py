import sys

bed_file = sys.argv[1]
bed_file_split = bed_file.split('/')
species=bed_file_split[-2]
base_name=bed_file_split[-1]
base_name_split=base_name.split('.')
chrom1=base_name_split[0]
chrom2=base_name_split[1]
if base_name_split[2] == 'bed':
	status='unfiltered'
elif base_name_split[2] == 'exclude_censat':
	status='exclude_censat'
elif base_name_split[2] == 'exclude_censat_telo_subterm':
	status='exclude_censat_telo_subterm'
 
tot_bp=0
tot_pid=0
for line in open(bed_file):
	line_split=line.split('\t')
	if float(line_split[4]) >= 990:
		pid=float(line_split[4])/10
		bp=int(line_split[2])-int(line_split[1])
		tot_bp+=bp
		tot_pid+=(pid*bp)

if tot_bp == 0:
	print('\t'.join([species, chrom1, chrom2, str(tot_bp), str(tot_pid)]))
else:
	tot_pid /= tot_bp
	print('\t'.join([species, chrom1, chrom2, str(tot_bp), str(tot_pid)]))