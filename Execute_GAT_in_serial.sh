#!/bin/bash

# Reccomended: 64G mem, 16 cores, 72h runtime
# With parralelization (e.g., the snakemake workflow we're finalizing), this example runs in under 20h.

# Install Script

module load singularity # if you do not have a singularity module on your cluster, do what else you need to have singularity in your environment

# parameters to find coding directories needed for tutorial in general (fix per cluster/machine)

##
### Values you need to change
##
export sourceDir=/.mounts/labs/simpsonlab/users/dsokolowski/miniconda3/bin/activate
export tutorialDir=/.mounts/labs/simpsonlab/users/dsokolowski/projects/GenomeAnnotationTutorial
export condaDir=/.mounts/labs/simpsonlab/users/dsokolowski/miniconda3/envs/annotation_tutorial

##
###
## 

##
### These values assume you have downloaded and processed mmus_GRC39 into data/references/mmus_GRC39
### The step-by-step guide to install this is in https://github.com/BaderLab/GenomeAnnotationTutorial/blob/main/setup/Preprocess%20Reference%20Species.md
### Or can be run with:
# mkdir -p mmus_GRC39 ; cd mmus_GRC39
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz
# for i in *.gz ; do gunzip $i ; echo $i ; done
# bash preprocess_reference_from_refseq.sh \
# $tutorialDir/data/references/human_T2T_NCBI \
# $tutorialDir \
# GCF_000001635.27_GRCm39_genomic.fna \
# GCF_000001635.27_GRCm39_genomic.gff
##


export refToga=mouse # reference species name for TOGA
export TogaDir=$tutorialDir/data/references/mmus_GRC39 # Needed to bind toga sif
export refTogaFa=$tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.fna # Files generated from setup/reference_directory_refseq.sh for grc399. You can change to other reference directories or ensembl with scripts in setup.
export refTogaBed=$tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.toga.bed # Files generated from setup/reference_directory_refseq.sh for grc399. You can change to other reference directories or ensembl with scripts in setup.
export refTogaIsoform=$tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.isoforms.toga.tsv # Files generated from setup/reference_directory_refseq.sh for grc399. You can change to other reference directories or ensembl with scripts in setup.

export refLiftOff=mouse # reference species name for liftOff. In this tutorial we use the same reference genome for TOGA and LiftOff
export refLiftOffFa=$tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.fna # directory in /data -- adding species can be done with scripts in /utils
export refLiftOffGff=$tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.gffread.gff # directory in /data -- adding species can be done with scripts in /utils

export orthofinderFA=$tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.nostop.protein.faa # directory name in ~/data
export orthofinderTab=$tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.table.txt # directory name in ~/data

##
###
##

##
### Values you do not need to change to run the tutorial. Many of these values will need to be updated for your own work though
##

export externalDir=$tutorialDir/external
export dataDir=$tutorialDir/data
export scriptsDir=$tutorialDir

source $sourceDir annotation_tutorial

# Parameters describing your assembly + annotation
export outDir=$tutorialDir/example_nmr
export target="example"
export species="heterocephalus_glaber"
export assemblyFile=$tutorialDir/data/example_data/example_data/NMRchr28.fa
export rnaseqDir=$tutorialDir/data/example_data/example_data/RNAseq_alignment
export isoseqDir=$tutorialDir/data/example_data/example_data/ISOseq_alignment
export customGFF="none" # if not none, then it should be the path to your gff file
export customRef="FALSE" # is the custom gff a reference to be upgraded. Switch to "TRUE" if it is a reference.
export liftoffRef="FALSE" # is the gff file from liftOff a reference assembly. Switch to "TRUE" if it is a reference. Note, if you have a custom reference gff then this will be over-written.
export MaskedAssemblyFile="none" # change to your softmasked .fa if that exists (e.g., /path-to-mastedasm/assembly.softmasked.fasta)
export MaskedAssemblyAnnotation="none" # change to your repeat annotation bed file if that exists (e.g., /path-to-mastedasm/assembly.softmasked.filteredRepeats.bed)
# Parameters describing your reference assemblies + annotation

##
### Executing tutorial
##

# Step 0 - directory structure
mkdir -p $outDir ; cd $outDir # go into directory

mkdir -p assembly
cp $assemblyFile assembly/assembly.fa

cd $outDir

mkdir -p RNAseq_alignment
cp $rnaseqDir/*.merged.bam RNAseq_alignment

mkdir -p ISOseq_alignment
cp $isoseqDir/*.merged.bam ISOseq_alignment

# Step 1: Repeat annotation and masking

bash $scriptsDir/scripts/run_earl_grey.sh

# Step 2: Generating protein-coding gene models

mkdir -p transcript_selection

# liftoff, toga, and stringtie, can be run in parallel


bash $scriptsDir/scripts/process_custom_gff.sh
bash $scriptsDir/scripts/run_liftOff.sh

bash $scriptsDir/scripts/run_stringtie_flexible.sh

bash $scriptsDir/scripts/sub_braker_flexible.sh # because we have sr and lr RNAseq, this script submits the following two scripts
# First script ran from sub_braker_flexible with sr RNA-seq and ISO-seq $scriptsDir/scripts/run_braker_sr.sh
# Second script ran from sub_braker_flexible with sr RNA-seq and ISO-seq $scriptsDir/scripts/run_braker_lr.sh

# Scripts for toga
bash $scriptsDir/scripts/make_cactus_tree.sh 
bash $scriptsDir/scripts/cactus_align_and_chain_sif.sh 
bash $scriptsDir/scripts/run_toga.sh

# Step 3: Combining and filtering gene models

bash $scriptsDir/scripts/12_mikado_configure_and_prepare.sh

# next three scripts can run in parallel
bash $scriptsDir/scripts/sub_mikado_blast.sh # normally splits fa into 100 jobs
bash $scriptsDir/scripts/get_junctions.sh
bash $scriptsDir/scripts/run_transdecoder.sh
 
# once previous three scripts are finished.
bash $scriptsDir/scripts/45_mikado_serialize_pick.sh

# Step 4: Annotating non-coding RNA genes

bash $scriptsDir/scripts/ncRNA_analysis.sh

bash $scriptsDir/scripts/run_mirmachine.sh

bash $scriptsDir/scripts/ncRNA_postprocess.sh

# Step 5: Sequence-similarity-based transfer of gene symbols

bash $scriptsDir/scripts/run_orthofinder.sh

bash $scriptsDir/scripts/gene_symbol_tables.sh


# Step 6: additional annotations

bash $scriptsDir/scripts/run_interproscan.sh

