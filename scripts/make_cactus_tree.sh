#!/bin/bash

# outDir=$1 # /scratch/example_isoseq
# refTogaFa=$2 $tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.fna # directory in /data -- adding species can be done with scripts in /utils
# refToga=$3 # mouse
# target=$4 # "example"

mkdir -p $outDir/cactus_aln

cd $outDir/cactus_aln

line1="($target:1.0,$refToga:1.0);"
line2=$target$'\t'$outDir/assembly/assembly.softmasked.fa
line3=$refToga$'\t'$refTogaFa

echo $line1$'\n'$line2$'\n'$line3 > cactus_in.txt

