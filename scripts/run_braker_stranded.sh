#!/bin/bash

# outDir=$1 # /scratch/example_isoseq
# tutorialDir=$2 # /.mounts/labs/simpsonlab/users/dsokolowski/projects/GenAnT
# target=$3 # "example"
# brakerOdbFaa=$4 # "Vertebrata.fa"


externalDir=$tutorialDir/external
dataDir=$tutorialDir/data

mkdir -p $outDir/braker_sr 
bash make_RNAseq_stranded.sh $outDir 16 # "16" represents 16 threads for samtools since it's the same number of threads as in braker, but this is flexible.

cd $outDir/braker_sr

species_suffix=$RANDOM

prefix=$target

# bams=`ls $outDir/RNAseq_alignment/*.merged.bam -m` # get bam files separated by csv
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

singularity exec --bind ${bamDir},${wd},${PWD},${assembly},${protDir},${configPath} ${BRAKER_SIF} braker.pl --AUGUSTUS_CONFIG_PATH=$configPath --genome=$assembly --prot_seq=$protDir/$brakerOdbFaa --bam=plus.bam,minus.bam --stranded=+,- --workingdir=${wd} --species=$prefix$species_suffix --threads 16  &> $wd/brakerlog.log

gffread $wd/braker.gtf --keep-genes -o $outDir/transcript_selection/braker.sr.gffread.gff
