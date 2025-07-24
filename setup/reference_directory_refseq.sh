#!/bin/bash

# examples are commented

wd=$1 # /.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/test/GenomeAnnotationTutorial/data/references/mmus_GRC39
GAT=$2 # /.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline
fa=$3 # GCF_000001635.27_GRCm39_genomic.fna
gff=$4 # GCF_000001635.27_GRCm39_genomic.gff

cd $wd


# clean GFF

# clean GFF1 for moving gene ID to gene symbol

gffread $gff --keep-genes -o $prefix.gffread.gff

sed '/^>/ s/ .*//' $fa | sed 's/[ryswkmbdhv]/N/gi' > $prefix.clean.fa

# Index FASTA file
samtools faidx $prefix.cleanhead.fa

echo "make protein faa for orthofinder"

gffread -y $prefix".protein.faa" -g $fa $prefix".gffread.gff"

# filters for gff's with pseudogenes containing premature stop codons. 
sed '/^>/!s/[.*]//g' $prefix".protein.faa" > $prefix".nostop.protein.faa" 


echo "make bed12 for TOGA"

$GAT/external/kent/gff3ToGenePred $prefix".gffread.gff" $prefix".genePred"

$GAT/external/kent/genePredToBed $prefix".genePred" $prefix".bed"

echo "Make gene tables for TOGA, Orthofinder, and gene symbol mapping"

Rscript --vanilla $GAT/setup/ncbi_gene_key.R -g $gff

awk 'FNR==NR { if (NR > 1) keep[$2]; next } $4 in keep' $prefix".isoforms.tsv" $prefix".bed" > $prefix".toga.bed"

