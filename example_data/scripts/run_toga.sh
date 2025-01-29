#!/bin/bash
#$ -l h_vmem=24G,h_rt=60:00:00,h_stack=32M
#$ -pe smp 8

cd $outDir

mkdir -p toga_out

cd $outDir/toga_out

source $sourceDir annotation_tutorial

export PATH="$externalDir/TOGA:$PATH"

refPrefix=$refToga

prefix=$target

supplyDir=$externalDir/TOGA/supply

ln -s $externalDir/TOGA/CESAR2.0 ./

toga.py $outDir/cactus_aln/target-to-ref.chain $refTogaBed $outDir/cactus_aln/$refPrefix.2bit $outDir/cactus_aln/$prefix.2bit --project_name $prefix"_toga" --nextflow_config_dir $externalDir/TOGA/nextflow_config_files/ --cesar_binary $externalDir/TOGA/CESAR2.0/cesar --isoforms $refTogaIsoform
