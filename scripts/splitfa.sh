#!/bin/bash

# Count the number of sequences in the FASTA file
num_sequences=$(grep -c "^>" $1)

# Calculate the number of sequences per file
sequences_per_file=$((num_sequences / $2))
remainder=$((num_sequences % $2))

# Initialize counters
seq_count=0
file_count=1

# Read the FASTA file and split it into 100 files
awk -v sequences_per_file="$sequences_per_file" -v remainder="$remainder" -v file_count="$file_count" -v seq_count="$seq_count" '
BEGIN { filename = sprintf("output_%d.fasta", file_count) }
{
    if ($0 ~ /^>/) {
        seq_count++
        # Start a new file if the limit for the current file is reached
        if ((seq_count > sequences_per_file && file_count <= $2 - remainder) || 
            (seq_count > sequences_per_file + 1 && file_count > $2 - remainder)) {
            file_count++
            seq_count = 1
            filename = sprintf("output_%d.fasta", file_count)
        }
    }
    print > filename
}' $1
