import pandas as pd
from Bio import SeqIO

# Load DESeq2 results and filter for 'Up' and 'Down'
df = pd.read_csv("sig_DESeq2_results.csv")

# Filter only rows where threshold is 'Up' or 'Down'
df_filtered = df[df["threshold"].isin(["Up", "Down"])]

# Create sets of gene names for Up and Down
up_genes = set(df_filtered[df_filtered["threshold"] == "Up"]["gene"])
down_genes = set(df_filtered[df_filtered["threshold"] == "Down"]["gene"])

# Parse the fasta file
fasta_file = "utr3_sequences.fa"
up_seqs = []
down_seqs = []

for record in SeqIO.parse(fasta_file, "fasta"):
    for gene in up_genes:
        if gene in record.id:
            up_seqs.append(record)
            break
    for gene in down_genes:
        if gene in record.id:
            down_seqs.append(record)
            break

# Write to output FASTA files
SeqIO.write(up_seqs, "utr3_up.fa", "fasta")
SeqIO.write(down_seqs, "utr3_down.fa", "fasta")
