#!/bin/bash

# outDir=$1 # /scratch/example_isoseq
# tutorialDir=$2 # /.mounts/labs/simpsonlab/users/dsokolowski/projects/GenAnT
# orthofinderTab=$3 # $tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.table.txt 
# refTogaIsoform=$4 # $tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.isoforms.toga.tsv
# refTogaIsoform2=$5 # "none"


cd $outDir

mkdir -p gene_symbol ; cd gene_symbol

# Copy annotation gff over
cp $outDir/ncRNA_analysis/full_annotation.gff ./

# Copy toga output over
cp $outDir/transcript_selection/toga.gffread.gff ./

# Copy liftoff output over

cp $outDir/transcript_selection/liftoff.gffread.gff ./


#
## preprocess liftoff, TOGA, and mikado
#

# mikado
grep -P "\texon\t" full_annotation.gff > mikado.exon.gff

# liftoff
grep -P "\texon\t" liftoff.gffread.gff > liftoff.exon.gff

# toga.r1
grep -P "\texon\t" toga.r1.gffread.gff > toga.r1.exon.gff

# toga.r2

grep -P "\texon\t" toga.r2.gffread.gff > toga.r2.exon.gff

# liftoff
bedtools intersect -a mikado.exon.gff -b liftoff.exon.gff -wo > mikado.liftoff.exon.txt

# toga.r1
bedtools intersect -a mikado.exon.gff -b toga.r1.exon.gff -wo > mikado.toga.r1.exon.txt

# toga.r2
bedtools intersect -a mikado.exon.gff -b toga.r2.exon.gff -wo > mikado.toga.r2.exon.txt


# Liftoff

cut -f1-9 mikado.liftoff.exon.txt > liftoff_overlap.mikadoInfo.gff
cut -f10-18 mikado.liftoff.exon.txt > liftoff_overlap.liftoffInfo.gff

# TOGA (r1)

cut -f1-9 mikado.toga.r1.exon.txt > toga_overlap.r1.mikadoInfo.gff
cut -f10-18 mikado.toga.r1.exon.txt > toga_overlap.r1.togaInfo.gff

# TOGA (r2)
cut -f1-9 mikado.toga.r2.exon.txt > toga_overlap.r2.mikadoInfo.gff
cut -f10-18 mikado.toga.r2.exon.txt > toga_overlap.r2.togaInfo.gff#
## preprocess liftoff, TOGA, and mikado
#

# mikado
grep -P "\texon\t" full_annotation.gff > mikado.exon.gff

# liftoff
grep -P "\texon\t" liftoff.gffread.gff > liftoff.exon.gff

# toga.r1
grep -P "\texon\t" toga.r1.gffread.gff > toga.r1.exon.gff

# toga.r2

grep -P "\texon\t" toga.r2.gffread.gff > toga.r2.exon.gff

# liftoff
bedtools intersect -a mikado.exon.gff -b liftoff.exon.gff -wo > mikado.liftoff.exon.txt

# toga.r1
bedtools intersect -a mikado.exon.gff -b toga.r1.exon.gff -wo > mikado.toga.r1.exon.txt

# toga.r2
bedtools intersect -a mikado.exon.gff -b toga.r2.exon.gff -wo > mikado.toga.r2.exon.txt


# Liftoff

cut -f1-9 mikado.liftoff.exon.txt > liftoff_overlap.mikadoInfo.gff
cut -f10-18 mikado.liftoff.exon.txt > liftoff_overlap.liftoffInfo.gff

# TOGA (r1)

cut -f1-9 mikado.toga.r1.exon.txt > toga_overlap.r1.mikadoInfo.gff
cut -f10-18 mikado.toga.r1.exon.txt > toga_overlap.r1.togaInfo.gff

# TOGA (r2)
cut -f1-9 mikado.toga.r2.exon.txt > toga_overlap.r2.mikadoInfo.gff
cut -f10-18 mikado.toga.r2.exon.txt > toga_overlap.r2.togaInfo.gff

# Copy orthofinder output over

cp $outDir/orthofinder/orthofinder_protein.tsv ./
cp $orthofinderTab ./reference.table.txt

cp $refTogaIsoform ./reference.toga.r1.table.txt
cp $refTogaIsoform2 ./reference.toga.r2.table.txt

# Start the gene symbol table and add in results from LiftOff and TOGA

Rscript --vanilla $tutorialDir/scripts/MakeGeneSymbolTableLiftOffTOGA.R

# Add orthofinder genes

Rscript --vanilla $tutorialDir/scripts/AddOrthoFinderGenes.R

# Generate final gff

Rscript --vanilla $tutorialDir/scripts/FormatFinalGFF.R


mv full_annotation.geneSymbols.gff $outDir
