#!/bin/bash

outDir=$1
externalDir=$2
target=$3
dataDir=$4

mkdir -p $outDir/braker_noRNA

cd $outDir/braker_noRNA

species_suffix=$RANDOM

prefix=$target

assembly=$outDir/assembly/assembly.softmasked.fa # /mHetGlaV3.soft.fa
protDir=$dataDir/braker_protein # /Vertebrata.fa
configPath=$condaDir/config

BRAKER_SIF=$externalDir/singularity_images/braker3.sif 

SINGULARITY_CACHEDIR=$outDir/braker_noRNA/cache
SINGULARITY_TMPDIR=$outDir/braker_noRNA/tmp

mkdir -p $SINGULARITY_CACHEDIR
mkdir -p $SINGULARITY_TMPDIR

wd=$outDir/braker_noRNA

singularity exec --bind ${bamDir},${wd},${PWD},${assembly},${protDir},${configPath} ${BRAKER_SIF} braker.pl --AUGUSTUS_CONFIG_PATH=$configPath --genome=$assembly --prot_seq=$protDir/Vertebrata.fa --workingdir=${wd} --species=$prefix$species_suffix --threads 16  &> $wd/brakerNoRNAseq.log

gffread $wd/braker.gtf --keep-genes -o $outDir/transcript_selection/braker.noRNA.gffread.gff


# singularity exec --bind ${bamDir},${wd},${PWD},${assembly},${protDir},${configPath} ${BRAKER_SIF} fix_gtf_ids.py --gtf $wd/braker.gtf --out braker.fix.gtf
