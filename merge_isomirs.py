#!/usr/bin/env python3

import os
import glob
import pandas as pd

# Step 1: Find all isomiR count files
files = glob.glob("*_isomir_counts.csv")

merged_df = None

for file in files:
    sample_name = file.replace("_isomir_counts.csv", "")
    df = pd.read_csv(file)

    # Step 2: Handle duplicate 'mir' values
    if df['mir'].duplicated().any():
        df['mir'] = df.apply(
            lambda row: f"{row['mir']}_{row['precursor']}" if df['mir'].duplicated(keep=False)[row.name] else row['mir'],
            axis=1
        )

    # Step 3: Rename freq column to sample name
    df = df[['mir', 'freq']].rename(columns={'freq': sample_name})

    # Step 4: Merge using INNER JOIN to keep only common 'mir'
    if merged_df is None:
        merged_df = df
    else:
        merged_df = pd.merge(merged_df, df, on='mir', how='inner')

# Step 5: Save merged result
merged_df.to_csv("merged_isomir_counts.csv", index=False)
print("Merged output saved to merged_isomir_counts.csv")

