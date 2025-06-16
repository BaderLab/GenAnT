#!/bin/bash

rfamDir=$dataDir/Rfam

cd $outDir

mkdir -p blast_seed

cd blast_seed

blastn -db $rfamDir/Rfam_nodup -query $outDir/assembly/assembly.softmasked.fa -evalue 1e-2 -outfmt 6 -out assembly.Rfam.blastn -num_threads 24

blastn -db $rfamDir/Rfam_nodup -query $outDir/assembly/assembly.softmasked.fa -evalue 1e-6 -outfmt 6 -out assembly.Rfame6.blastn -num_threads 24


awk -F "\t" '{print $1 "\t" $7 "\t" $8 "\t" $2}' assembly.Rfam.blastn > assembly.rfam.bed
