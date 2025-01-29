#!/bin/bash
#$ -l h_vmem=24G,h_rt=200:00:00,h_stack=32M
#$ -pe smp 16

# Install Script

cd /.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/example_run

module load singularity


# parameters to find coding directories needed for tutorial in general (fix per cluster/machine)

export sourceDir=/.mounts/labs/simpsonlab/users/dsokolowski/miniconda3/bin/activate
export externalDir=/.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/external
export dataDir=/.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/data
export condaDir=/.mounts/labs/simpsonlab/users/dsokolowski/miniconda3/envs/annotation_tutorial

source $sourceDir annotation_tutorial


# Parameters describing your assembly + annotation
export outDir=/.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/example_run
export target="example"
export assemblyFile=/.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/data/example_data/NMRchr28.fa
export rnaseqDir=/.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/data/example_data/RNAseq
export repeatBed=/.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/data/example_data/NMRchr28.filteredRepeats.bed # only required when inputting a masked assembly, otherwise this is made in step 1

# Parameters describing your reference assemblies + annotation

export refToga=mm10
export TogaDir=/.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/data/references/mm10
export refTogaFa=/.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/data/references/mm10/mm10.fa # directory in /data -- adding species can be done with scripts in /utils
export refTogaBed=/.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/data/references/mm10/mm10.v25.for_toga.bed
export refTogaIsoform=/.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/data/references/mm10/mm10.v25.for_toga.isoforms.tsv

export refLiftOff=mm10
export refLiftOffFa=/.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/data/references/mm10/mm10.fa # directory in /data -- adding species can be done with scripts in /utils
export refLiftOffGff=/.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/data/references/mm10/mm10.gff # directory in /data -- adding species can be done with scripts in /utils

export orthofinderRef=mmus # directory name in ~/data

# Step 0 - directory structure
mkdir -p $outDir ; cd $outDir # go into directory

mkdir -p assembly
cp $assemblyFile assembly/assembly.fa

cd $outDir

mkdir -p RNAseq_alignment
 cp $rnaseqDir/*.merged.bam RNAseq_alignment

export PATH="$condaDir/share/RepeatMasker/":$PATH
export PATH="$condaDir/bin/":$PATH

# earlGrey -g $outDir/assembly/assembly.fa -s rattus_rattus -o . -t 50 -r rodentia -d yes

# Step 1: Repeat annotation and masking

	# This example assumes earl grey is run to save >100 of hours of runtime
	# the earl grey command is `earlGrey -g $assemblyFile -s heterocephalus_glaber -o . -t 50 -r rodentia -d yes
	# This step would normally copy the masked fa and the filtered repeats to "assembly"

cp $assemblyFile assembly/assembly.softmasked.fa
cp $repeatBed assembly/assembly.filteredRepeats.bed

# Step 2: Generating protein-coding gene models

# liftoff, toga, and stringtie, can be run in parallel
bash $outDir/scripts/run_liftOff.sh
bash $outDir/scripts/toga_liner.sh
bash $outDir/scripts/run_braker.sh

bash $outDir/scripts/run_stringtie.sh
bash $outDir/scripts/run_stringtie_merge.sh #after run_stringtie is complete

# Step 3: Combining and filtering gene models

bash $outDir/scripts/get_gffread.sh
bash $outDir/scripts/12_mikado_configure_and_prepare.sh

# next three scripts can run in parallel
bash $outDir/scripts/sub_mikado_blast.sh # normally splits fa into 100 jobs. 
bash $outDir/scripts/run_porticullis.sh 
bash $outDir/scripts/run_transdecoder.sh

# once previous three scripts are finished.
bash $outDir/scripts/45_mikado_serialize_pick.sh

# Step 4: Annotating non-coding RNA genes



# Step 5: Sequence-similarity-based transfer of gene symbols

bash $outDir/scripts/run_orthofinder.sh
