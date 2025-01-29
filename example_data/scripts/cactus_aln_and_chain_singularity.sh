#!/bin/bash
#$ -l h_vmem=24G,h_rt=200:00:00,h_stack=32M
#$ -pe smp 16


cactus_venv=$externalDir/cactus-bin-v2.8.4/venv-cactus-v2.8.4/bin/activate
cactus_bin=$externalDir/cactus-bin-v2.8.4/bin
kent_bin=$externalDir/kent

source $cactus_venv

export PATH="$cactus_bin:$PATH"
export PATH="$kent_bin:$PATH" 

cd $outDir/cactus_aln

cactus jobstore cactus_in.txt target_ref.hal --binariesMode local --realTimeLogging True --batchSystem single_machine --workDir $outDir/cactus_aln

j=$target # target
k=$refToga # reference

halStats --genomes target_ref.hal
# Anc0 MolAen mouse

hal2fasta target_ref.hal $k | faToTwoBit stdin $k.2bit
hal2fasta target_ref.hal $j | faToTwoBit stdin $j.2bit

halStats --bedSequences $k target_ref.hal > $k.bed
halStats --bedSequences $j target_ref.hal > $j.bed


halLiftover --outPSL target_ref.hal $j \
      $j.bed $k /dev/stdout | \
      pslPosTarget stdin $k"-to-"$j".psl"

axtChain -psl -linearGap=loose $k"-to-"$j".psl" $k.2bit $j.2bit "target-to-ref.chain"

