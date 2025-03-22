#!/bin/bash
#$ -l h_vmem=24G,h_rt=20:00:00,h_stack=32M
#$ -pe smp 8

outDir=$1
target=$2
refToga=$3
refTogaFa=$4

mkdir -p $outDir/cactus_aln

cd $outDir/cactus_aln

line1="($target:1.0,$refToga:1.0);"
line2=$target$'\t'$outDir/assembly/assembly.softmasked.fa
line3=$refToga$'\t'$refTogaFa

echo $line1$'\n'$line2$'\n'$line3 > cactus_in.txt

