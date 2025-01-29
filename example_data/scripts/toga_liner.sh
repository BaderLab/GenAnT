#!/bin/bash
#$ -l h_vmem=24G,h_rt=200:00:00,h_stack=32M
#$ -pe smp 16

cd $outDir/scripts


# bash make_cactus_tree.sh 
# bash cactus_align_and_chain_sif.sh 
bash run_toga.sh
