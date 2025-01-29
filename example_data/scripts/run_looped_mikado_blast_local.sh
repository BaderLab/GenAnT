#!/bin/bash


i=$1

n=`basename $i .fasta`


wd=$outDir/transcript_selection/blast

cd $wd


ASSEMBLY=$outDir/assembly/assembly.softmasked.fa

UNIDB=$dataDir/uniprot_sprot

mkdir -p $wd/$n

cd $UNIDB

# Paths to files
QUERY_FILE=$wd/$i   # Replace with your FASTA file
SPROT_DB=$UNIDB"/uniprot_sprot_nodup.fasta"             # Replace with the location of Swiss-Prot BLAST DB
OUTPUT_DIR=$wd/$n

blastx -max_target_seqs 5 -num_threads 20 -query "$QUERY_FILE" \
 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop" \
 -db "$SPROT_DB" -evalue 0.000001 -out $OUTPUT_DIR/blast_results.tsv
