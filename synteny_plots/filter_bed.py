import pandas as pd
import sys

# def merge_bed_intervals(bed_df):
# 	"""Merge adjacent BED intervals with the same strand and satellite type."""
# 	merged = []
# 	current = bed_df.iloc[0]

# 	for _, row in bed_df.iterrows():
# 		if (
# 			row['chrom'] == current['chrom'] and
# 			row['strand'] == current['strand'] and
# 			row['name'] == current['name'] and
# 			row['start'] <= current['end'] + 1
# 		):
# 			# Merge intervals
# 			current['end'] = max(current['end'], row['end'])
# 		else:
# 			# Append current interval to merged list
# 			merged.append(current)
# 			current = row
	
# 	# Add the last interval
# 	merged.append(current)
# 	return pd.DataFrame(merged)

# def filter_intervals(bed_df):
# 	filtered = []

# 	for i, row in bed_df.iterrows():
# 		if i == 0:
# 			# Always keep the first row; check proximity for the second row
# 			filtered.append(row)
# 		else:
# 			# Check proximity to the previous row
# 			prev_row = bed_df.iloc[i - 1]
# 			if (
# 				row['chrom'] == prev_row['chrom'] and
# 				abs(row['start'] - prev_row['end']) <= 50000
# 			):
# 				# Keep both rows if within 50kb
# 				if prev_row not in filtered:
# 					filtered.append(prev_row)
# 				filtered.append(row)
# 			elif (row['end'] - row['start']) >= 50000:
# 				filtered.append(row)

# 	return pd.DataFrame(filtered)

# bed = pd.read_csv(sys.stdin, sep="\t", header=None)

# print(bed.head())
# print(bed.shape)
# print(bed.shape[1])

# if bed.shape[1] < 3 or bed.shape[1] > 9:
# 	sys.stderr.write("Error: BED file must have between 3 and 9 columns.\n")
# 	sys.exit(1)

# # Assign standard column names for known BED structure
# columns = ["chrom", "start", "end", "name", "score", "strand", "extra1", "extra2", "extra3"][:bed.shape[1]]
# bed.columns = columns
# bed = bed.sort_values(by=['chrom', 'start', 'end', 'strand']).reset_index(drop=True)
# merged_bed = merge_bed_intervals(bed)

# filtered_bed = filter_intervals(merged_bed)

# # Output to stdout
# bed.to_csv(sys.stdout, sep="\t", index=False, header=False)
import sys

def parse_bed_line(line):
	"""Parse a single BED line into its components."""
	fields = line.strip().split("\t")
	chrom = fields[0]
	start = int(fields[1])
	end = int(fields[2])
	name = fields[3]
	score = fields[4]
	strand = fields[5]
	thick_start = fields[6]
	thick_end = fields[7]
	rgb = fields[8]
	return chrom, start, end, name, score, strand, thick_start, thick_end, rgb

def merge_bed_entries(entries):
	"""Merge adjacent or overlapping BED entries with the same col4 and col5 values."""
	merged = []
	current_entry = entries[0]

	for i in range(1, len(entries)):
		chrom, start, end, name, score, strand, thick_start, thick_end, rgb = entries[i]
		curr_chrom, curr_start, curr_end, curr_name, curr_score, curr_strand, curr_ts, curr_te, curr_rgb = current_entry

		# Check if the current entry can be merged with the next, check strand as only the cens will have strand info bc of the way the bed file is formatted
		if (
			chrom == curr_chrom and name == curr_name and strand == curr_strand and
			start <= curr_end + 10000  # Overlapping or within 10kb
		):
			current_entry = (curr_chrom, curr_start, max(curr_end, end), curr_name, curr_score, curr_strand, curr_ts, curr_te, curr_rgb)
		else:
			merged.append(current_entry)
			current_entry = (chrom, start, end, name, score, strand, thick_start, thick_end, rgb)

	# Add the last entry
	merged.append(current_entry)
	return merged

def filter_bed_entries(entries):
	"""Filter entries based on size and proximity conditions."""
	large_entries = [
		(chrom, start, end)
		for chrom, start, end, name, score, strand, thick_start, thick_end, rgb in entries
		if end - start >= 10000
	]
	filtered = []

	for entry in entries:
		chrom, start, end, name, score, strand, thick_start, thick_end, rgb = entry
		size = end - start

		# If the entry is at least 50kb, always keep it
		if size >= 50000:
			filtered.append(entry)
			continue

		# Otherwise, check proximity to neighbors
		has_large_neighbor = any(
			large_chrom == chrom and (abs(start - large_end) <= 10000 or abs(large_start - end) <= 10000)
			for large_chrom, large_start, large_end in large_entries
		)
			
		if has_large_neighbor:
			filtered.append(entry)

		# If no neighboring entry is within 10kb, discard this entry
	return filtered

def rgb_to_hex(rgb):
	"""Convert an RGB string to a hexadecimal string."""
	r, g, b = map(int, rgb.split(","))
	return "{:02x}{:02x}{:02x}".format(r, g, b)

def apply_color(entries):
	recolored_entries = []
	for entry in entries:
		chrom, start, end, name, score, strand, thick_start, thick_end, rgb = entry
		lower_name = name.lower()
		if 'ct' in lower_name: #first so it gets overriden by anything else
			rgb = '255,255,255'
		if 'gap' in lower_name:
			rgb = '224,224,224'
		if 'hsat1b' in lower_name or 'hsat2' in lower_name or 'gsat' in lower_name or 'censat' in lower_name:
			rgb = '0,204,204'
		if 'dhor' in lower_name or 'mixedalpha' in lower_name or 'hor' in lower_name or 'mon' in lower_name:
			rgb = '255,146,0'
		if 'active_hor' in lower_name or 'as_strand' in lower_name: #after so it overwrites color for active_hor
			rgb = '153,0,0'
		if 'hsat1a' in lower_name:
			rgb = '0,222,96'
		if 'hsat3' in lower_name:
			rgb = '51,81,137'
		if 'bsat' in lower_name:
			rgb = '250,153,255'
		if 'rdna' in lower_name:
			rgb = '0,0,0'
		if 'sst1' in lower_name:
			rgb = '172,51,199'
		if 'acro1' in lower_name:
			rgb = '102,47,144'
		rgb = rgb_to_hex(rgb)
		recolored_entries.append((chrom, start, end, name, score, strand, thick_start, thick_end, rgb))
	return recolored_entries
		
def restrand_entries(entries):
	restranded_entries = []
	for entry in entries:
		chrom, start, end, name, score, strand, thick_start, thick_end, rgb = entry
		if strand == '.':
			strand = '*'
		entry = (chrom, start, end, name, score, strand, thick_start, thick_end, rgb)
		restranded_entries.append(entry)
	return restranded_entries
	
def main():
	"""Main function to process BED data from stdin."""
	bed_entries = []

	# Read from stdin
	for line in sys.stdin:
		if line.strip():
			bed_entries.append(parse_bed_line(line))

	# Sort entries by chromosome and start position
	bed_entries.sort(key=lambda x: (x[0], x[1]))

	# Merge entries
	merged_entries = merge_bed_entries(bed_entries)

	# Filter entries
	filtered_entries = filter_bed_entries(merged_entries)

	recolored_entries = apply_color(filtered_entries)

	# restranded_entries = restrand_entries(recolored_entries)

	# Output the results
	for entry in recolored_entries:
		print("\t".join(map(str, entry)))

if __name__ == "__main__":
	main()
