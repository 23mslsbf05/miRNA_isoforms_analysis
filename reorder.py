import pandas as pd

# Load metadata.txt (TSV file) and get the desired column order from its first column
metadata = pd.read_csv('metadata.txt', sep='\t', header=None)
desired_order = metadata.iloc[:, 0].tolist()

# Load the merged.csv file
merged = pd.read_csv('merged_isomir_counts.csv')

# Ensure 'mir' is first
desired_order = [col for col in desired_order if col != 'mir']  # remove if already present
ordered_columns = ['mir'] + [col for col in desired_order if col in merged.columns]

# Append remaining columns not in metadata.txt or 'mir'
remaining_columns = [col for col in merged.columns if col not in ordered_columns]
final_columns = ordered_columns + remaining_columns

# Reorder and save
merged[final_columns].to_csv('reordered.csv', index=False)

