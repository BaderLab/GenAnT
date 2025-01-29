#!/bin/bash

cd $outDir

mkdir -p earl_grey ; cd earl_grey

export PATH="$condaDir/share/RepeatMasker/":$PATH
export PATH="$condaDir/bin/":$PATH

earlGrey -g $outDir/assembly/assembly.fa -s $species -o . -t 50 -r rodentia -d yes

cp earl_grey/$species"_EarlGrey"/$species"_summaryFiles/"$species.softmasked.fasta $outDir/assembly/assembly.softmasked.fa
cp earl_grey/$species"_EarlGrey"/$species"_summaryFiles/"$species.filteredRepeats.bed $outDir/assembly/assembly.filteredRepeats.bed
