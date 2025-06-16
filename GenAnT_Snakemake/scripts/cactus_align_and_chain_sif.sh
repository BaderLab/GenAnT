#!/bin/bash

# {params.outDir} {params.externalDir} {params.TogaDir} {params.target} {params.refToga}
outDir=$1
externalDir=$2
TogaDir=$3
target=$4
refToga=$5
round=$6

# singularity build cactus.v2.9.3.sif docker://quay.io/comparative-genomics-toolkit/cactus:v2.9.3

# singularity exec cactus.v2.9.3.sif hal2fasta halLiftover

CACTUS_SIF=$externalDir/singularity_images/cactus.v2.9.3.sif

kent_bin=$externalDir/kent

# export PATH="$kent_bin:$PATH" 

cd $outDir/cactus_aln$round

wd=$outDir/cactus_aln$round
targetDir=$outDir/assembly
singularity exec --bind ${wd},${TogaDir},${targetDir} ${CACTUS_SIF} cactus \
jobstore $wd/cactus_in.txt $wd/target_ref.hal --binariesMode local \
--realTimeLogging True --batchSystem single_machine --workDir $wd

j=$target # target
k=$refToga # reference

# Anc0 MolAen mouse

singularity exec --bind ${wd} ${CACTUS_SIF} hal2fasta $wd/target_ref.hal $k > $wd/$k.fa
$kent_bin/faToTwoBit $k.fa $k.2bit
singularity exec --bind ${wd} ${CACTUS_SIF} hal2fasta $wd/target_ref.hal $j > $wd/$j.fa
$kent_bin/faToTwoBit $j.fa $j.2bit

singularity exec --bind ${wd} ${CACTUS_SIF} halStats --bedSequences $k $wd/target_ref.hal > $wd/$k.bed
singularity exec --bind ${wd} ${CACTUS_SIF} halStats --bedSequences $j $wd/target_ref.hal > $wd/$j.bed


singularity exec --bind ${wd} ${CACTUS_SIF} halLiftover --outPSL target_ref.hal $j \
      $j.bed $k /dev/stdout | \
      $kent_bin/pslPosTarget stdin $k"-to-"$j".psl"

$kent_bin/axtChain -psl -linearGap=loose $k"-to-"$j".psl" $k.2bit $j.2bit "target-to-ref.chain"

