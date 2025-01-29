#!/bin/bash

cd $outDir/transcript_selection

ALIGNMENTS=$outDir/RNAseq_alignment

ASSEMBLY=$outDir/assembly/assembly.softmasked.fa

samtools merge $outDir/transcript_selection/merged.bam $ALIGNMENTS/*merged.bam

portcullis full -t 20 $ASSEMBLY $outDir/transcript_selection/merged.bam -o $outDir/transcript_selection
