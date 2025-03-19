#!/bin/bash

cd $outDir

if [[ $MaskedAssemblyFile = "none" ]] ; then

	mkdir -p earl_grey ; cd earl_grey

	earlGrey -g $outDir/assembly/assembly.fa -s $species -o . -t 50 -r rodentia -d yes

	cd $outDir

	cp earl_grey/$species"_EarlGrey"/$species"_summaryFiles/"$species.softmasked.fasta $outDir/assembly/assembly.softmasked.fa
	cp earl_grey/$species"_EarlGrey"/$species"_summaryFiles/"$species.filteredRepeats.bed $outDir/assembly/assembly.filteredRepeats.bed

else

	cp $MaskedAssemblyFile $outDir/assembly/assembly.softmasked.fa
	cp $MaskedAssemblyAnnotation $outDir/assembly/assembly.filteredRepeats.bed
fi

