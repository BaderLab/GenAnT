#!/bin/bash
#$ -l h_vmem=24G,h_rt=60:00:00,h_stack=32M
#$ -pe smp 8

cd $outDir

mkdir -p toga_out

cd $outDir/toga_out

source $sourceDir annotation_tutorial

kent_bin=$externalDir/kent

refPrefix=$refToga

prefix=$target

supplyDir=$externalDir/TOGA/supply

ln -s $externalDir/TOGA/CESAR2.0 ./

$externalDir/TOGA/toga.py $outDir/cactus_aln/target-to-ref.chain $refTogaBed $outDir/cactus_aln/$refPrefix.2bit $outDir/cactus_aln/$prefix.2bit --project_name $prefix"_toga" --nextflow_config_dir $externalDir/TOGA/nextflow_config_files/ --cesar_binary $externalDir/TOGA/CESAR2.0/cesar --isoforms $refTogaIsoform

# $externalDir/TOGA/toga.py $outDir/cactus_aln/target-to-ref.chain $refTogaBed $outDir/cactus_aln/$refPrefix.2bit $outDir/cactus_aln/$prefix.2bit --project_name $prefix"_toga" --nextflow_config_dir $externalDir/TOGA/nextflow_config_files/ --cesar_binary $externalDir/TOGA/CESAR2.0/cesar --isoforms $refTogaIsoform

togaDir=$outDir/toga_out

$kent_bin/bedToGenePred $togaDir/$prefix"_toga"/query_annotation.bed $togaDir/$prefix"_toga"/query_annotation.genePred

$kent_bin/genePredToGtf file $togaDir/$prefix"_toga"/query_annotation.genePred $togaDir/$prefix"_toga"/query_annotation.gtf

gffread $togaDir/$prefix"_toga"/query_annotation.gtf --keep-genes -o $outDir/transcript_selection/toga.gffread.gff
gffread -y $outDir/transcript_selection/toga.faa -g $outDir/assembly/assembly.softmasked.fa $outDir/transcript_selection/toga.gffread.gff
