#!/bin/bash

outDir=$1
externalDir=$2
target=$3
refToga=$4
refTogaBed=$5
refTogaIsoform=$6
round=$7

cd $outDir

mkdir -p toga_out$round

cd $outDir/toga_out$round

# export PATH="$externalDir/TOGA:$PATH"

kent_bin=$externalDir/kent
# export PATH="$kent_bin:$PATH"

refPrefix=$refToga

prefix=$target

supplyDir=$externalDir/TOGA/supply

ln -s $externalDir/TOGA/CESAR2.0 ./

# theTogaBin=/.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/toga_take3

$externalDir/TOGA/toga.py $outDir/cactus_aln$round/target-to-ref.chain $refTogaBed $outDir/cactus_aln$round/$refPrefix.2bit $outDir/cactus_aln$round/$prefix.2bit --project_name $prefix"_toga" --nextflow_config_dir $externalDir/TOGA/nextflow_config_files/ --cesar_binary $externalDir/TOGA/CESAR2.0/cesar --isoforms $refTogaIsoform

togaDir=$outDir/toga_out$round

$kent_bin/bedToGenePred $togaDir/$prefix"_toga"/query_annotation.bed $togaDir/$prefix"_toga"/query_annotation.genePred

$kent_bin/genePredToGtf file $togaDir/$prefix"_toga"/query_annotation.genePred $togaDir/$prefix"_toga"/query_annotation.gtf

gffread $togaDir/$prefix"_toga"/query_annotation.gtf --keep-genes -o $outDir/transcript_selection/toga.$round.gffread.gff
gffread -y $outDir/transcript_selection/toga.$round.faa -g $outDir/assembly/assembly.softmasked.fa $outDir/transcript_selection/toga.$round.gffread.gff
