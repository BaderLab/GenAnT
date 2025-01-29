#!/bin/bash
#$ -l h_vmem=24G,h_rt=20:00:00,h_stack=32M
#$ -pe smp 8

export PATH="$externalDir/stringtie:$PATH"

cd $outDir
cd RNAseq_alignment

stringtie --merge *gtf -o stringtie.merged.gtf 



