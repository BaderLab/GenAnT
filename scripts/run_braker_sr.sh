#!/bin/bash

mkdir -p $outDir/braker_sr

cd $outDir/braker_sr

species_suffix=$RANDOM

prefix=$target

bams=`ls $outDir/RNAseq_alignment/*.merged.bam -m` # get bam files separated by csv
bams2=$(echo $bams | sed 's/ //g')

assembly=$outDir/assembly/assembly.softmasked.fa # /mHetGlaV3.soft.fa
protDir=$dataDir/braker_protein # /Vertebrata.fa
configPath=$externalDir/Augustus/config

BRAKER_SIF=$externalDir/singularity_images/braker3.sif 

SINGULARITY_CACHEDIR=$outDir/braker_sr/cachedir
SINGULARITY_TMPDIR=$outDir/braker_sr/tmpdir

mkdir -p $SINGULARITY_CACHEDIR
mkdir -p $SINGULARITY_TMPDIR

wd=$outDir/braker_sr

bamDir=$outDir/RNAseq_alignment

singularity exec --bind ${bamDir},${wd},${PWD},${assembly},${protDir},${configPath} ${BRAKER_SIF} braker.pl --AUGUSTUS_CONFIG_PATH=$configPath --genome=$assembly --prot_seq=$protDir/$brakerOdbFaa --bam=$bams2 --workingdir=${wd} --species=$prefix$species_suffix --threads 16  &> $wd/brakerlog.log

gffread $wd/braker.gtf --keep-genes -o $outDir/transcript_selection/braker.sr.gffread.gff
