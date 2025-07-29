#!/bin/bash

# outDir=$1 # /scratch/example_isoseq
# tutorialDir=$2 # /.mounts/labs/simpsonlab/users/dsokolowski/projects/GenAnT

externalDir=$tutorialDir/external

cd $outDir

mkdir -p interproscan ; cd interproscan

$externalDir/my_interproscan/interproscan-5.69-101.0/interproscan.sh -i $outDir/transcript_selection/mikado_lenient.faa -f tsv -dp
