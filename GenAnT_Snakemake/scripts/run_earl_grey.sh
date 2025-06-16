#!/bin/bash

outDir=$1
externalDir=$2
species=$3
MaskedAssemblyFile=$4
MaskedAssemblyAnnotation=$5
dfamDB=$6

cd $outDir

if [[ $MaskedAssemblyFile = "none" ]] ; then

	mkdir -p earl_grey ; cd earl_grey

	SINGULARITY_CACHEDIR=$outDir/earl_grey/cache
	SINGULARITY_TMPDIR=$outDir/earl_grey/tmp
	EARLGREY_SIF=$externalDir/singularity_images/earlgrey.sif


	mkdir -p $SINGULARITY_CACHEDIR
	mkdir -p $SINGULARITY_TMPDIR

	singularity exec --bind ${outDir} ${EARLGREY_SIF} earlGrey -g $outDir/assembly/assembly.fa -s $species -o . -t 50 -r $dfamDB/$brakerOdbFaa -d yes

	cd $outDir

	cp earl_grey/$species"_EarlGrey"/$species"_summaryFiles/"$species.softmasked.fasta $outDir/assembly/assembly.softmasked.fa
	cp earl_grey/$species"_EarlGrey"/$species"_summaryFiles/"$species.filteredRepeats.bed $outDir/assembly/assembly.filteredRepeats.bed

else

	cp $MaskedAssemblyFile $outDir/assembly/assembly.softmasked.fa
	cp $MaskedAssemblyAnnotation $outDir/assembly/assembly.filteredRepeats.bed
fi

