#!/bin/bash

cd $outDir

mkdir -p mirmachine ; cd mirmachine

ASSEMBLY=$outDir/assembly/assembly.softmasked.fa

samtools faidx $ASSEMBLY

MirMachine.py -n $mirmachineClade -s $species --genome $ASSEMBLY -m deutero --cpu 20
