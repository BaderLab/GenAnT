#!/bin/bash

# examples are commented

wd=$1 # /.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/test/GenomeAnnotationTutorial/data/references/mmus_GRC39_embl
GAT=$2 # /.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/test
fa=$3 # Mus_musculus.GRCm39.dna_sm.primary_assembly.fa
gff=$4 # Mus_musculus.GRCm39.113.gff3

echo "Clean gff file for reference"

prefix=`basename $gff .gff3`

# clean GFF

# clean GFF1 for moving gene ID to gene symbol

# awk 'BEGIN {OFS="\t"} {$9 = gensub(/transcript:|gene:/, "", "g", $9); print}' input.gff > output.gff

# awk 'BEGIN {OFS="\t"} {$9 = gensub(/transcript:|gene:/, "", "g", $9); print}' $gff > $prefix.clean.gff3

awk -F'\t' 'BEGIN {OFS="\t"} {
  if (NF >= 9) {
    for (i = 1; i <= 8; i++) printf "%s\t", $i;
    sub(/transcript:/, "", $9);
    sub(/gene:/, "", $9);
    gsub(/transcript:/, "", $9);
    gsub(/gene:/, "", $9);
    print $9
  } else {
    print
  }
}' $gff > $prefix.clean.gff3

gffread $prefix.clean.gff3 --keep-genes -o $prefix.gffread.gff


# removes additional fasta headers that get confused with different tools
sed '/^>/ s/ .*//' $fa > $prefix.cleanhead.fa

# replaces any sequence that is not a clear amino acid with N (missing)
awk '/^>/ {print; next} {gsub(/[^TCAGNtcagn]/, "N"); print}' $prefix.cleanhead.fa > $prefix.cleanseq.fa


# Index FASTA file
samtools faidx $prefix.cleanseq.fa
# find chromosomes in GFF and fasta
cut -f1 $prefix.gffread.gff | grep -v "^#" | sort | uniq > $prefix.gff.chroms.txt

# include sequences that contain a gene only
samtools faidx $prefix.cleanseq.fa $(cat ${prefix}.gff.chroms.txt) > $prefix.clean.fa

# tidy extra fasta file

rm $prefix.cleanhead.fa
rm $prefix.cleanseq.fa



# make protein faa for orthofinder

gffread -y $prefix".protein.faa" -g $fa $prefix".gffread.gff"

# filters for gff's with pseudogenes containing premature stop codons. 
sed  '/^>/!s/[.*]//g'  $prefix".protein.faa" > $prefix".nostop.protein.faa" 

# make bed12 for TOGA

$GAT/external/kent/gff3ToGenePred $prefix".gffread.gff" $prefix".genePred"

$GAT/external/kent/genePredToBed $prefix".genePred" $prefix".bed"

# Make gene tables for TOGA, Orthofinder, and gene symbol mapping

Rscript --vanilla $GAT/setup/ensembl_gene_key.R -g $prefix.clean.gff3

awk 'FNR==NR { if (NR > 1) keep[$2]; next } $4 in keep' $prefix".clean.isoforms.tsv" $prefix".bed" > $prefix".clean.bed"
