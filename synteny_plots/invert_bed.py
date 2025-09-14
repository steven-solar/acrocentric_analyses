import pandas as pd
import sys

fai_file = sys.argv[1]

chrom_lengths = pd.read_csv(fai_file, sep="\t", header=None, usecols=[0, 1], names=["chrom", "length"])
chrom_lengths_dict = dict(zip(chrom_lengths["chrom"], chrom_lengths["length"]))

bed = pd.read_csv(sys.stdin, sep="\t", header=None)

if bed.shape[1] < 3 or bed.shape[1] > 9:
	sys.stderr.write("Error: BED file must have between 3 and 9 columns.\n")
	sys.exit(1)

# Assign standard column names for known BED structure
columns = ["chrom", "start", "end", "name", "score", "strand", "extra1", "extra2", "extra3"][:bed.shape[1]]
bed.columns = columns

# Filter rows to ensure chromosomes exist in the lengths file
bed = bed[bed["chrom"].isin(chrom_lengths_dict.keys())]

# Invert the start and end positions based on chromosome lengths
bed["new_start"] = bed.apply(lambda row: chrom_lengths_dict[row["chrom"]] - row["end"], axis=1)
bed["new_end"] = bed.apply(lambda row: chrom_lengths_dict[row["chrom"]] - row["start"], axis=1)

# Adjust the strand if present
if "strand" in bed.columns:
	bed["strand"] = bed["strand"].map({"+": "-", "-": "+", ".": "."})

# Update the DataFrame with the new positions
bed["start"] = bed["new_start"]
bed["end"] = bed["new_end"]
bed = bed.drop(columns=["new_start", "new_end"])

# Output to stdout
bed.to_csv(sys.stdout, sep="\t", index=False, header=False)
