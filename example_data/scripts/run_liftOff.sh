#!/bin/bash
#$ -l h_vmem=24G,h_rt=200:00:00,h_stack=32M
#$ -pe smp 16

cd $outDir

mkdir -p liftoff

cd liftoff

liftoff -g $refLiftOffGff $outDir/assembly/assembly.softmasked.fa $refLiftOffFa -o ./liftoff.gff -u ./unmapped.liftoff.txt -copies
~                                                                                                                                          
