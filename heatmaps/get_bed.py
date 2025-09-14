import sys
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def percent_identity_to_rgb(percent_identity, min_percent_identity=80):
	if percent_identity == 0:
		return "255,255,255"
	if percent_identity <= 1:
		percent_identity=percent_identity*100
	if percent_identity < min_percent_identity:
		return "255,255,220"
	# Define the colormap
	cmap = plt.get_cmap('YlOrRd')

	# Define the normalization range
	norm = mcolors.Normalize(vmin=min_percent_identity, vmax=100)

	# Get the RGBA color corresponding to the percent identity
	rgba = cmap(norm(percent_identity))

	# Convert RGBA to RGB
	rgb = tuple(int(255 * x) for x in rgba[:3])
	return str(rgb[0]) + "," + str(rgb[1]) + "," + str(rgb[2])

mashmap_fp = sys.argv[1]
mashmap_f = open(mashmap_fp)

for line in mashmap_f:
	line_split = line.strip().split('\t')
	chrom, start, end = line_split[0], int(line_split[2]), int(line_split[3])
	idy = float(line_split[12].split(':')[2])
	if idy<=1:
		idy*=100
	q_chr, q_start, q_end = line_split[5], int(line_split[7]), int(line_split[8])
	score_str = str(round(idy,1))
	name = f"{q_chr}:{q_start}-{q_end}_{score_str}%"
	strand="."
	thick_start = start
	thick_end = end
	print('\t'.join([chrom, str(start), str(end), name, str(int(round(10*idy,0))), strand, str(thick_start), str(thick_end), percent_identity_to_rgb(idy)]))
