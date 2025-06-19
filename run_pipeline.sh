#!/bin/bash

# Check if GNU parallel is installed
if ! command -v parallel &> /dev/null; then
    echo "GNU parallel is not installed. Install it and try again."
    exit 1
fi

# Check if at least one input FASTQ file is given
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <input_fastq1> [input_fastq2 ...]"
    exit 1
fi

# Function to process a single file
process_sample() {
    fq_file="$1"
    sample_name=$(basename "$fq_file" .fastq)

    echo "Processing $sample_name..."

    # Step 1: Trimming
    trim_galore -j 1 -e 0.2 -q 10 --small_rna --stringency 6 --fastqc "$fq_file"

    trimmed_file="${sample_name}_trimmed.fq"

    # Step 2: Collapse
    seqcluster collapse -f "$trimmed_file" -o .

    collapsed_file="${sample_name}_trimmed_trimmed.fastq"

    # Step 3: Convert to FASTA
    fasta_file="${sample_name}.fa"
    seqtk seq -A "$collapsed_file" > "$fasta_file"

    # Step 4: miraligner
    miraligner_output="miraligner_${sample_name}"
    miraligner -minl 16 -sub 0 -s rno -i "$fasta_file" \
        -db . \
        -pre . \
        -o "$miraligner_output"

    # Step 5: Generate isomiR count
    python generate_isomir_count.py "${miraligner_output}.mirna" "${sample_name}_isomir_counts.csv"

    echo "Finished processing $sample_name."
}

export -f process_sample

# Run the pipeline in parallel on 4 jobs
parallel -j 4 process_sample ::: "$@"

# Wait for all parallel jobs to finish before merging
wait

# Final merge step using Python
python merge_isomirs.py

