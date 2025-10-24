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

# clean GFF1 for moving gene ID to gene symbol

awk '$7 == "+" || $7 == "-" ' $prefix.gff > $prefix.strand.gff

gffread $prefix.strand.gff -o $prefix.gffread.gff

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

prefix=`basename $gff .gff`

# clean GFF

echo "make protein faa for orthofinder"

gffread -y $prefix".protein.faa" -g $fa $prefix".gffread.gff"

# filters for gff's with pseudogenes containing premature stop codons. 
# Secondary filtering for unrecognized protein
sed -i '/^>/!s/[.*]//g' $prefix".protein.faa" 
sed -i '/^>/! s/\./X/g' $prefix".protein.faa" 


echo "make bed12 for TOGA"

$GAT/external/kent/gff3ToGenePred $prefix".gffread.gff" $prefix".genePred"

$GAT/external/kent/genePredToBed $prefix".genePred" $prefix".bed"

echo "Make gene tables for TOGA, Orthofinder, and gene symbol mapping"f

Rscript --vanilla $GAT/setup/mikado_gene_key.R -g $gff


awk 'FNR==NR { if (NR > 1) keep[$2]; next } $4 in keep' $prefix".isoforms.tsv" $prefix".bed" > $prefix".toga.bed"

# further cleaning of toga for certain GFF files with length 0 CDS annotations (found in GCF_048164855.1_mCanLup2.hap1)


awk '$8 <= $7 { next } { print }' $prefix.toga.bed > $prefix.toga.clean.bed

cut -f4 $prefix.isoforms.tsv $prefix.toga.clean.bed | sort -u > toga.transcripts.txt

awk 'FNR==NR { keep[$1]; next }
     FNR==1 { print; next }
     ($2 in keep)' toga.transcripts.txt $prefix.isoforms.tsv > $prefix.isoforms.clean.tsv

rm toga.transcripts.txt
