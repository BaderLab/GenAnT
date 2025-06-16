#!/bin/bash
#$ -l h_vmem=24G,h_rt=20:00:00,h_stack=32M
#$ -pe smp 1

# Install Script

cd /.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/example_isoseq

module load singularity

# parameters to find coding directories needed for tutorial in general (fix per cluster/machine)

export sourceDir=/.mounts/labs/simpsonlab/users/dsokolowski/miniconda3/bin/activate
export tutorialDir=/.mounts/labs/simpsonlab/users/dsokolowski/projects/GenomeAnnotationTutorial
export externalDir=$tutorialDir/external
# export externalDir=/.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/test/GenomeAnnotationTutorial/external

# tool-specific parameters to delineate clade.
export dfamDB="rodentia" # change to a different repeatmasker class if desired, e.g., "arthropoda" for insects
export brakerOdbFaa="Vertebrata.fa" # change to a different odb12 protein set. e.g.,  "Arthropoda.fa" # https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/Arthropoda.fa.gz to get arthopods
export mikadoScore="mammalian.yaml" # Use a different scoring file for a different clade "HISTORIC/dmelanogaster_scoring.yaml"
export mirmachineClade="Mammalia" # "mirmachine-clade" # "Glires" 

export dataDir=$tutorialDir/data
export scriptsDir=$tutorialDir
export condaDir=/.mounts/labs/simpsonlab/users/dsokolowski/miniconda3/envs/annotation_tutorial

source $sourceDir annotation_tutorial

# Parameters describing your assembly + annotation
export outDir=$tutorialDir/example_isoseq
export target="example"
export species="heterocephalus_glaber"
export assemblyFile=$tutorialDir/data/example_data/example_data/NMRchr28.fa
export rnaseqDir=$tutorialDir/data/example_data/example_data/RNAseq_alignment
export isoseqDir=$tutorialDir/data/example_data/example_data/ISOseq_alignment
export customGFF="none" # if not none, then it should be the path to your gff file
export customRef="FALSE" # is the custom gff a reference to be upgraded. Switch to "TRUE" if it is a reference.
export liftoffRef="FALSE" # is the gff file from liftOff a reference assembly. Switch to "TRUE" if it is a reference. Note, if you have a custom reference gff then this will be over-written.
export MaskedAssemblyFile="none"
export MaskedAssemblyAnnotation="none"
# Parameters describing your reference assemblies + annotation

export refToga=mouse
export TogaDir=$tutorialDir/data/references/mmus_GRC39 # Possibly needed to bind toga sif
export refTogaFa=$tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.fna # directory in /data -- adding species can be done with scripts in /utils
export refTogaBed=$tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.toga.bed
export refTogaIsoform=$tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.isoforms.toga.tsv


export refLiftOff=mouse
export refLiftOffFa=$tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.fna # directory in /data -- adding species can be done with scripts in /utils
export refLiftOffGff=$tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.gffread.gff # directory in /data -- adding species can be done with scripts in /utils

export orthofinderFA=$tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.nostop.protein.faa # directory name in ~/data
export orthofinderTab=$tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.table.txt # directory name in ~/data

# Step 0 - directory structure
mkdir -p $outDir ; cd $outDir # go into directory

mkdir -p assembly
# cp $assemblyFile assembly/assembly.fa

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

# Repeat the LiftOff scripts as many times as you need, updating the variables below. We reccomend limiting the number of reference assemblies (we limit to 3) to avoid mikado linking together weak annotations into a long superlocus.

	## export refLiftOff=mouse
	## export refLiftOffFa=$tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.fna # directory in /data -- adding species can be done with scripts in /utils
	## export refLiftOffGff=$tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.gffread.gff # directory in /data -- adding species can be done with scripts in /utils
# Note// We have not written adding gene symbols for a second and third round of liftOff into GenAnT/1.0.0. You can manually add these with MakeGeneSymbolTableLiftOffTOGA.R. We will add this functionality into a later update.

# bash $scriptsDir/scripts/run_stringtie_flexible.sh

# bash $scriptsDir/scripts/sub_braker_flexible.sh # because we have sr and lr RNAseq, this script submits the following two scripts
bash $scriptsDir/scripts/run_braker_sr.sh
bash $scriptsDir/scripts/run_braker_lr.sh



# Scripts for toga
bash $scriptsDir/scripts/make_cactus_tree.sh 
bash $scriptsDir/scripts/cactus_align_and_chain_sif.sh 
bash $scriptsDir/scripts/run_toga.sh

##
 # Repeat the TOGA scripts as many times as you need, updating the variables below. We reccomend limiting the number of reference assemblies (we limit to 2) to avoid mikado linking together weak annotations into a long superlocus.
 	# Write the output as toga.r2.gffread.gff for the these gene models to automatically be added to our gene symbol table. Otherwise add the gene symbols with MakeGeneSymbolTableLiftOffTOGA.R like in LiftOff.

 # export refToga=mouse
 # export TogaDir=$tutorialDir/data/references/mmus_GRC39 # Possibly needed to bind toga sif
 # export refTogaFa=$tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.fna # directory in /data -- adding species can be done with scripts in /utils
 # export refTogaBed=$tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.toga.bed
 # export refTogaIsoform=$tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.isoforms.toga.tsv

##

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

# We currently commented out line 42 "cp $refTogaIsoform2 ./reference.toga.r2.table.txt". Uncomment this line if you did two rounds of TOGA
bash $scriptsDir/scripts/aggregate_symbols.sh


# Step 6: additional annotations

bash $scriptsDir/scripts/run_interproscan.sh

