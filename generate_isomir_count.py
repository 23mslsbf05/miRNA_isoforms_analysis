import pandas as pd
import argparse

# -------------------- Function to remove suffix --------------------
def remove_suffix(seq, add):
    if str(seq).endswith(str(add)):
        return seq[:-len(add)]
    return seq

# -------------------- Function to check lowercase or zero --------------------
def is_lower_or_zero(val):
    val = str(val)
    return val == '0' or val.islower()

# -------------------- Function to move 'mir' column to first --------------------
def move_mir_column_first(df):
    cols = list(df.columns)
    if 'mir' in cols:
        cols.insert(0, cols.pop(cols.index('mir')))
        return df[cols]
    else:
        raise ValueError("'mir' column not found in DataFrame.")

# -------------------- Argument Parsing --------------------
parser = argparse.ArgumentParser(description='Process isomiR data from miraligner output.')
parser.add_argument('input_file', help='Path to the input .mirna file')
parser.add_argument('output_file', help='Path for the output CSV file')
args = parser.parse_args()

input_file = args.input_file
output_file = args.output_file

# -------------------- Step 1: Load and preprocess input --------------------
df = pd.read_csv(input_file, sep='\t')
df['seq'] = df.apply(lambda row: remove_suffix(row['seq'], row['add']), axis=1)
df = df.sort_values(by=['mir', 'freq'], ascending=[True, False])
df['group_flag'] = df['t5'].apply(is_lower_or_zero)

group_df = df[df['group_flag']]
grouped = group_df.groupby(['mir', 'precursor'], as_index=False).agg({
    'seq': 'first',
    'add': 'first',
    't5': 'first',
    'freq': 'sum'
})
final_df = pd.concat([grouped, df[~df['group_flag']]], ignore_index=True)
final_df = final_df.sort_values(by=['mir', 'freq'], ascending=[True, False])

# -------------------- Step 2: Combine by 'precursor' and 't5' --------------------
grouped_df = final_df.groupby(['precursor', 't5'], as_index=False).agg({
    'freq': 'sum',
    'seq': 'first',
    'add': 'first',
    'mir': 'first'
})

# -------------------- Step 3: Modify 'mir' column with 't5' prefix unless t5 == '0' --------------------
def prefix_mir_with_t5(t5, mir):
    return f"{t5}_{mir}" if str(t5) != '0' else str(mir)

grouped_df['mir'] = grouped_df.apply(lambda row: prefix_mir_with_t5(row['t5'], row['mir']), axis=1)

# -------------------- Step 4: Move 'mir' column to first --------------------
grouped_df = move_mir_column_first(grouped_df)

# -------------------- Final Output --------------------
grouped_df.to_csv(output_file, index=False)
print(f"Final isomiR count saved to '{output_file}'")

