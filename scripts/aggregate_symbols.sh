#!/bin/bash

outDir=$1
orthofinderTab=$2
snakeDir=$3

cd $outDir

mkdir -p gene_symbol ; cd gene_symbol

# Copy annotation gff over
cp $outDir/ncRNA_analysis/full_annotation.gff ./

# Copy toga output over
cp $outDir/transcript_selection/toga.gffread.gff ./

# Copy liftoff output over

cp $outDir/transcript_selection/liftoff.gffread.gff ./

# preprocess liftoff, TOGA, and mikado

grep -P "\tmRNA\t|\tlncRNA\t" full_annotation.gff > mikado.mRNA.lncRNA.gff
grep -P "\tmRNA\t|\tlnc_RNA\t|\ttranscript\t" liftoff.gffread.gff > liftoff.mRNA.lncRNA.gff
grep -P "\tmRNA\t|\tlnc_RNA\t|\ttranscript\t" toga.gffread.gff > toga.mRNA.lncRNA.gff

bedtools intersect -a mikado.mRNA.lncRNA.gff -b liftoff.mRNA.lncRNA.gff -wo > mikado.liftoff.mRNA.lncRNA.txt
bedtools intersect -a mikado.mRNA.lncRNA.gff -b toga.mRNA.lncRNA.gff -wo > mikado.toga.mRNA.lncRNA.txt

# Liftoff

cut -f1-9 mikado.liftoff.mRNA.lncRNA.txt > liftoff_overlap.mikadoInfo.gff
cut -f10-18 mikado.liftoff.mRNA.lncRNA.txt > liftoff_overlap.liftoffInfo.gff

# TOGA

cut -f1-9 mikado.toga.mRNA.lncRNA.txt > toga_overlap.mikadoInfo.gff
cut -f10-18 mikado.toga.mRNA.lncRNA.txt > toga_overlap.togaInfo.gff

# Copy orthofinder output over

cp $outDir/orthofinder/orthofinder_protein.tsv ./
cp $orthofinderTab ./reference.table.txt

# Start the gene symbol table and add in results from LiftOff and TOGA

Rscript --vanilla $snakeDir/scripts/MakeGeneSymbolTableLiftOffTOGA.R

# Add orthofinder genes

Rscript --vanilla $snakeDir/scripts/AddOrthoFinderGenes.R

# Generate final gff

Rscript --vanilla $snakeDir/scripts/FormatFinalGFF.R


mv full_annotation.geneSymbols.gff $outDir