import pandas as pd

# Load the file
df = pd.read_csv("all_species.exclude_censat_telo_subterm.tsv", sep="\t", header=None)

# Keep only the first row of each unique (col1, col2, col3) combination
df_filtered = df.drop_duplicates(subset=[0, 1, 2], keep="first")

# Save the output
df_filtered.to_csv("exclude_censat_telo_subterm.tsv", sep="\t", index=False, header=False)
