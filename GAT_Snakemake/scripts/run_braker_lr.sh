#!/bin/bash

outDir=$1
externalDir=$2
target=$3
condaDir=$4
dataDir=$5

mkdir -p $outDir/braker_lr

cd $outDir/braker_lr

species_suffix=$RANDOM

prefix=$target

bams=`ls $outDir/ISOseq_alignment/*.merged.bam -m` # get bam files separated by csv
bams2=$(echo $bams | sed 's/ //g')

assembly=$outDir/assembly/assembly.softmasked.fa # /mHetGlaV3.soft.fa
protDir=$dataDir/braker_protein # /Vertebrata.fa
configPath=$externalDir/Augustus/config

BRAKER_SIF=$externalDir/singularity_images/braker3_lr.sif 

SINGULARITY_CACHEDIR=$outDir/braker_lr/cache
SINGULARITY_TMPDIR=$outDir/braker_lr/tmp

mkdir -p $SINGULARITY_CACHEDIR
mkdir -p $SINGULARITY_TMPDIR

bamDir=$outDir/ISOseq_alignment

wd=$outDir/braker_lr

singularity exec --bind ${bamDir},${wd},${PWD},${assembly},${protDir},${configPath} ${BRAKER_SIF} braker.pl --AUGUSTUS_CONFIG_PATH=$configPath --genome=$assembly --prot_seq=$protDir/Vertebrata.fa --bam=$bams2 --workingdir=${wd} --species=$prefix$species_suffix --threads 16  &> $wd/brakerlog.log

gffread $wd/braker.gtf --keep-genes -o $outDir/transcript_selection/braker.lr.gffread.gff
