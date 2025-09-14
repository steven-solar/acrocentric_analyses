import pandas as pd

df1 = pd.read_csv('sats.bed', sep='\t', header=None)
df2 = pd.read_csv('sats.strand.bed', sep='\t', header=None)

merged_df = pd.merge(df1, df2, on=[0, 1, 2])
print(merged_df.head())
merged_df = merged_df[[0,1,2,'3_x','4_x','5_y','6_x','7_x','8_x']]
print(merged_df.head())
merged_df.to_csv('merged_sats.bed', sep='\t', header=False, index=False)
