#!/bin/bash

outDir=$1
orthofinderTab=$2
snakeDir=$3


cd $outDir

mkdir -p gene_symbol ; cd gene_symbol

# Copy annotation gff over
cp $outDir/ncRNA_analysis/full_annotation.gff ./

# Copy toga output over
cp $outDir/transcript_selection/toga.r1.gffread.gff ./

# Copy liftoff output over

cp $outDir/transcript_selection/liftoff.gffread.gff ./

# [ -f source_file ] && cp source_file destination/

if [ -f $outDir/transcript_selection/toga.r2.gffread.gff ]; then
    cp $outDir/transcript_selection/toga.r2.gffread.gff ./

else
	touch ./toga.r2.gffread.gff

fi

#
## preprocess liftoff, TOGA, and mikado
#

# mikado
grep -P "\tmRNA\t|\tlncRNA\t" full_annotation.gff > mikado.mRNA.lncRNA.gff

# liftoff
grep -P "\tmRNA\t|\tlnc_RNA\t|\ttranscript\t" liftoff.gffread.gff > liftoff.mRNA.lncRNA.gff

# toga.r1
grep -P "\tmRNA\t|\tlnc_RNA\t|\ttranscript\t" toga.r1.gffread.gff > toga.r1.mRNA.lncRNA.gff

# toga.r2

grep -P "\tmRNA\t|\tlnc_RNA\t|\ttranscript\t" toga.r2.gffread.gff > toga.r2.mRNA.lncRNA.gff

# liftoff
bedtools intersect -a mikado.mRNA.lncRNA.gff -b liftoff.mRNA.lncRNA.gff -wo > mikado.liftoff.mRNA.lncRNA.txt

# toga.r1
bedtools intersect -a mikado.mRNA.lncRNA.gff -b toga.r1.mRNA.lncRNA.gff -wo > mikado.toga.r1.mRNA.lncRNA.txt

# toga.r2
bedtools intersect -a mikado.mRNA.lncRNA.gff -b toga.r2.mRNA.lncRNA.gff -wo > mikado.toga.r2.mRNA.lncRNA.txt


# Liftoff

cut -f1-9 mikado.liftoff.mRNA.lncRNA.txt > liftoff_overlap.mikadoInfo.gff
cut -f10-18 mikado.liftoff.mRNA.lncRNA.txt > liftoff_overlap.liftoffInfo.gff

# TOGA (r1)

cut -f1-9 mikado.toga.r1.mRNA.lncRNA.txt > toga_overlap.r1.mikadoInfo.gff
cut -f10-18 mikado.toga.r1.mRNA.lncRNA.txt > toga_overlap.r1.togaInfo.gff

# TOGA (r2)
cut -f1-9 mikado.toga.r2.mRNA.lncRNA.txt > toga_overlap.r2.mikadoInfo.gff
cut -f10-18 mikado.toga.r2.mRNA.lncRNA.txt > toga_overlap.r2.togaInfo.gff

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
