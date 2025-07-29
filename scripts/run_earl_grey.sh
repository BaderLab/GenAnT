#!/bin/bash

# outDir=$1 # /scratch/example_isoseq
# MaskedAssemblyFile=$2 # "none"
# MaskedAssemblyAnnotation=$3 # "none"
# species=$4 # "heterocephalus_glaber"
# dfamDB=$5 # "rodentia"

cd $outDir

if [[ $MaskedAssemblyFile = "none" ]] ; then

	mkdir -p earl_grey ; cd earl_grey

	earlGrey -g $outDir/assembly/assembly.fa -s $species -o . -t 16 -r $dfamDB -d yes

	cd $outDir

	cp earl_grey/$species"_EarlGrey"/$species"_summaryFiles/"$species.softmasked.fasta $outDir/assembly/assembly.softmasked.fa
	cp earl_grey/$species"_EarlGrey"/$species"_summaryFiles/"$species.filteredRepeats.bed $outDir/assembly/assembly.filteredRepeats.bed

else

	cp $MaskedAssemblyFile $outDir/assembly/assembly.softmasked.fa
	cp $MaskedAssemblyAnnotation $outDir/assembly/assembly.filteredRepeats.bed
fi

