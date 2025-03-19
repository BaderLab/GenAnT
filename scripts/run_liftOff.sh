#!/bin/bash

cd $outDir

mkdir -p liftoff

cd liftoff

liftoff -g $refLiftOffGff $outDir/assembly/assembly.softmasked.fa $refLiftOffFa -o ./liftoff.gff -u ./unmapped.liftoff.txt -copies

gffread ./liftoff.gff --keep-genes -o $outDir/transcript_selection/liftoff.gffread.gff
                                                                                                         
