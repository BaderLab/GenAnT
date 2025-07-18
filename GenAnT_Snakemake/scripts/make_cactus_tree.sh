#!/bin/bash

outDir=$1
target=$2
refToga=$3
refTogaFa=$4
round=$5

mkdir -p $outDir/cactus_aln$round

cd $outDir/cactus_aln$round

line1="($target:1.0,$refToga:1.0);"
line2=$target$'\t'$outDir/assembly/assembly.softmasked.fa
line3=$refToga$'\t'$refTogaFa

echo $line1$'\n'$line2$'\n'$line3 > cactus_in.txt

