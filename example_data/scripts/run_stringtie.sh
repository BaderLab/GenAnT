#!/bin/bash
#$ -l h_vmem=24G,h_rt=20:00:00,h_stack=32M
#$ -pe smp 8

export PATH="$externalDir/stringtie:$PATH"

cd $outDir
cd RNAseq_alignment

for i in *.bam ; do

b=`basename $i .bam`

stringtie $i -l $b -o $b".gtf" -p 8 --conservative

done


