#!/bin/bash

# outDir=$1 # /scratch/example_isoseq
# refLiftOffGff=$2 # $tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.gffread.gff # directory in /data -- adding species can be done with scripts in /utils
# refLiftOffFa=$3 # $tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.fna # directory in /data -- adding species can be done with scripts in /utils
cd $outDir

mkdir -p liftoff

cd liftoff

liftoff -g $refLiftOffGff $outDir/assembly/assembly.softmasked.fa $refLiftOffFa -o ./liftoff.gff -u ./unmapped.liftoff.txt -copies

gffread ./liftoff.gff --keep-genes -o $outDir/transcript_selection/liftoff.gffread.gff
                                                                                                         
