#!/bin/bash

# examples are commented

wd=$1 # /.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/test/GenomeAnnotationTutorial/data/references/mmus_GRC39_embl
GAT=$2 # /.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/test
fa=$3 # Mus_musculus.GRCm39.dna_sm.primary_assembly.fa
gff=$4 # Mus_musculus.GRCm39.113.gff3

prefix=`basename $gff .gff3`

# clean GFF

# clean GFF1 for moving gene ID to gene symbol

gffread $gff --keep-genes -o $prefix.gffread.gff

sed '/^>/ s/ .*//' $fa > $prefix.cleanname.fa

# Index FASTA file
samtools faidx $prefix.cleanname.fa

cut -f1 $prefix.gffread.gff | grep -v "^#" | sort | uniq > $prefix.gff.chroms.txt

samtools faidx $prefix.cleanname.fa $(cat ${prefix}.gff.chroms.txt) > $prefix.clean.fa

rm $prefix.cleanname.fa

# Index FASTA file
samtools faidx $prefix.clean.fa

echo "Clean gff file for reference"

# make protein faa for orthofinder

gffread -y $prefix".protein.faa" -g $fa $prefix".gffread.gff"

# filters for gffs with pseudogenes containing premature stop codons. 
sed  '/^>/!s/[.*]//g'  $prefix".protein.faa" > $prefix".nostop.protein.faa" 

# make bed12 for TOGA

$GAT/external/kent/gff3ToGenePred $prefix".gffread.gff" $prefix".genePred"

$GAT/external/kent/genePredToBed $prefix".genePred" $prefix".bed"

# Make gene tables for TOGA, Orthofinder, and gene symbol mapping

Rscript --vanilla $GAT/setup/ensembl_gene_key.R -g $gff

awk 'FNR==NR { if (NR > 1) keep[$2]; next } $4 in keep' $prefix".isoforms.tsv" $prefix".bed" > $prefix".toga.bed"

