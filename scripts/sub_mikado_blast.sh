#!/bin/bash
#$ -l h_vmem=14G,h_rt=10:00:00,h_stack=32M

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

