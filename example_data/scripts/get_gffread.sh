#!/bin/bash

kent_bin=$externalDir/kent
export PATH="$kent_bin:$PATH"

cd $outDir

mkdir -p transcript_selection ; cd transcript_selection

PREFIX=$target
togaDir=$outDir/toga_out
liftoff=$outDir/liftoff
stringtieDir=$outDir/stringtie_mergeonly
brakerDir=$outDir/braker
ASMPATH=$outDir/assembly/assembly.softmasked.fa

gffread $brakerDir/braker.gtf --keep-genes -o braker.gffread.gff
gffread -y braker.faa -g $ASMPATH braker.gffread.gff

gffread $liftOffDir/liftoff.gff --keep-genes -o liftoff.gffread.gff
gffread -y liftoff.faa -g $ASMPATH liftoff.gffread.gff

mv $stringtieDir/stringtie.merged.gtf stringtie_merged.gffread.gff
gffread -y stringtie_merged.faa -g $ASMPATH stringtie_merged.gffread.gff

bedToGenePred $togaDir/query_annotation.bed $togaDir/query_annotation.genePred 

genePredToGtf file $togaDir/query_annotation.genePred $togaDir/query_annotation.gtf

gffread $togaDir/query_annotation.gtf --keep-genes -o toga.gffread.gff
gffread -y toga.faa -g $ASMPATH toga.gffread.gff
