#!/bin/bash

# outDir=$1 # /scratch/example_isoseq
# tutorialDir=$2 # /.mounts/labs/simpsonlab/users/dsokolowski/projects/GenAnT


cd $outDir/transcript_selection

mkdir -p blast

cd blast

ln -s $outDir/transcript_selection/mikado_prepared.fasta ./

bash $tutorialDir/scripts/splitfa.sh mikado_prepared.fasta 100

unlink mikado_prepared.fasta


# when local
for i in *fasta ; do bash $tutorialDir/scripts/run_looped_mikado_blast_local.sh $i ; done

# when subbed

# for i in *fasta ; do qsub -N $i -P simpsonlab -cwd -V -v I=$i $outDir/scripts/run_looped_mikado_blast.sh ; done

