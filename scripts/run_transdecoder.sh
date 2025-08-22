#!/bin/bash

mkdir -p $outDir/transcript_selection/transdecoder

cd $outDir/transcript_selection/transdecoder

TransDecoder.LongOrfs -t $outDir/transcript_selection/mikado_prepared.fasta -m 100 -O "transdecoder"
TransDecoder.Predict -t $outDir/transcript_selection/mikado_prepared.fasta --retain_long_orfs_length 100 -O "transdecoder"
