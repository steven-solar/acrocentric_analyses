import sys
import pandas as pd

# Load TSV (no header)
filename = sys.argv[1]
df = pd.read_csv(filename, sep='\t', header=None)
df.columns = ['species', 'ref_chrom', 'tgt_chrom', 'bp', 'pid']

# Convert to appropriate types
df['bp'] = pd.to_numeric(df['bp'])
df['pid'] = pd.to_numeric(df['pid'])

# Filter out rows with 0 shared_seq and 0 pid
filtered = df[~(df['bp'] == 0)]
filtered = df[~(df['species'] == 'chm13')]


# Group by species
grouped = filtered.groupby('species')

# Compute average shared sequence and weighted PID
summary = grouped.apply(lambda g: pd.Series({
    'avg_bp': g['bp'].mean(),
    'weighted_pid': (g['bp'] * g['pid']).sum() / g['bp'].sum()
})).reset_index()

summary['avg_bp'] = summary['avg_bp'].round(0).astype(int)
summary['weighted_pid'] = summary['weighted_pid'].map('{:.2f}'.format)

# Get row with max shared_seq per species
max_rows = filtered.loc[filtered.groupby('species')['bp'].idxmax()].reset_index(drop=True)
max_rows['bp'] = max_rows['bp'].round(0).astype(int)
max_rows['pid'] = max_rows['pid'].map('{:.2f}'.format)

# Disable scientific notation for printing
pd.set_option('display.float_format', '{:.2f}'.format)

print("Summary per species:\n", summary.to_string(index=False))
print("\nMax bp row per species:\n", max_rows.to_string(index=False))

total_bp = df["bp"].sum()
avg_bp = total_bp/6
avg_mb = avg_bp / 1000000

weighted_pid = (df["bp"] * df["pid"]).sum() / total_bp

# Output with formatting
print()
print(f"Avg shared Mbp: {avg_mb:,.1f} Mb")
# print("Average shared bp per row:", round(avg_bp / len(df)))
print("Weighted average PID: {:.2f}".format(weighted_pid))