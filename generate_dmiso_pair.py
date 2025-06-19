import os
import glob
import csv
from Bio import SeqIO

# Input files
query_fasta = "query.fa"
utr3_fasta = "utr3_sequences.fa"
input_files = glob.glob("*.csv")

# Load query.fa sequences into a dictionary
query_seqs = {record.id: str(record.seq) for record in SeqIO.parse(query_fasta, "fasta")}

# Load utr3_sequences.fa sequences into a dictionary
utr3_seqs = {record.id: str(record.seq) for record in SeqIO.parse(utr3_fasta, "fasta")}

# Open output file
with open("dmiso_pair.tsv", "w") as out_file:
    writer = csv.writer(out_file, delimiter='\t')
    # Header
    writer.writerow(["File", "ID1", "ID2", "QuerySeq", "UTR3Seq"])
    
    # Process each CSV file
    for file in input_files:
        base_name = os.path.splitext(os.path.basename(file))[0]
        with open(file, "r") as f:
            for line in f:
                parts = line.strip().split(",")
                if len(parts) < 3:
                    continue
                id1 = parts[0]
                id2 = parts[2]

                query_seq = query_seqs.get(id1, "")
                utr3_seq = utr3_seqs.get(id2, "")

                if query_seq and utr3_seq:
                    writer.writerow([base_name, id1, id2, query_seq, utr3_seq])
