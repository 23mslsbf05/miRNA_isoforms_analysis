import pandas as pd

# Step 1: Read the input file
df = pd.read_csv('sig_DESeq2_results.csv')

# Step 2: Modify 'gene' column to remove everything before and including the first '_'
df['gene'] = df['gene'].str.replace(r'^.*?_', '', regex=True)

# Step 3: Keep only rows where 'gene' appears at least twice
gene_counts = df['gene'].value_counts()
df = df[df['gene'].isin(gene_counts[gene_counts >= 2].index)]

# Step 4: Keep only rows where 'gene' has different 'threshold' values
# This is done by grouping and checking for more than one unique threshold per gene
filtered_df = df.groupby('gene').filter(lambda x: x['threshold'].nunique() > 1)

# Step 5: Write to output file
filtered_df.to_csv('de_isomir.csv', index=False)
