#!/bin/bash

# outDir=$1 # /scratch/example_isoseq

cd $outDir/transcript_selection

ALIGNMENTS=$outDir/RNAseq_alignment

ASSEMBLY=$outDir/assembly/assembly.softmasked.fa

samtools merge $outDir/transcript_selection/merged.bam $ALIGNMENTS/*merged.bam

portcullis full -t 20 $ASSEMBLY $outDir/transcript_selection/merged.bam -o $outDir/transcript_selection
