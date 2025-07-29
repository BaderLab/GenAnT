#!/bin/bash

# outDir=$1 # /scratch/example_isoseq
# mirmachineClade=$2 # "Mammalia"
# species=$3 # "heterocephalus_glaber"


cd $outDir

mkdir -p mirmachine ; cd mirmachine

ASSEMBLY=$outDir/assembly/assembly.softmasked.fa

samtools faidx $ASSEMBLY

MirMachine.py -n $mirmachineClade -s $species --genome $ASSEMBLY -m deutero --cpu 20
