#!/bin/bash

# examples are commented

wd=$1 # /.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/test/GenomeAnnotationTutorial/data/references/mmus_GRC39
GAT=$2 # /.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline
fa=$3 # GCF_000001635.27_GRCm39_genomic.fna
gff=$4 # GCF_000001635.27_GRCm39_genomic.gff

cd $wd

echo "Clean gff file for reference"

prefix=`basename $gff .gff`

# clean GFF

gffread $gff --keep-genes -o $prefix".gffread.gff"

echo "make protein faa for orthofinder"

gffread -y $prefix".protein.faa" -g $fa $prefix".gffread.gff"

# filters for gff's with pseudogenes containing premature stop codons. 
sed '/^>/!s/[.*]//g' $prefix".protein.faa" > $prefix".nostop.protein.faa" 


echo "make bed12 for TOGA"

$GAT/external/kent/gff3ToGenePred $prefix".gffread.gff" $prefix".genePred"

$GAT/external/kent/genePredToBed $prefix".genePred" $prefix".bed"

awk 'FNR==NR { if (NR > 1) keep[$2]; next } $4 in keep' $prefix".isoforms.tsv" $prefix".bed" > $prefix".toga.bed"


echo "Make gene tables for TOGA, Orthofinder, and gene symbol mapping"

Rscript --vanilla $GAT/setup/ncbi_gene_key.R -g $gff
