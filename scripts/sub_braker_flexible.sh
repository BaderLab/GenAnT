#!/bin/bash

# No RNA-seq data
# No RNA-seq data

cd $outDir

mkdir -p $outDir/stringtie_out 

if [[ $(ls -A $outDir/RNAseq_alignment | wc -l) -eq 0 && $(ls -A $outDir/ISOseq_alignment | wc -l) -eq 0 ]] ; then

	echo "We did not find short read or long read RNA-seq data (aligned bam files)"

	bash $tutorialDir/scripts/run_braker_noRNAseq.sh

	echo "" > $outDir/transcript_selection/braker.sr.gffread.gff # dummy variable for later
	echo "" > $outDir/transcript_selection/braker.lr.gffread.gff # dummy variable for later

fi

if [[ $(ls -A $outDir/RNAseq_alignment | wc -l) -gt 0 && $(ls -A $outDir/ISOseq_alignment | wc -l) -gt 0 ]] ; then

	echo "We detected short read and long read RNAseq data"

	bash $tutorialDir/scripts/run_braker_sr.sh # each script is normally 100+ hours and therefore should be qsubbed.
	bash $tutorialDir/scripts/run_braker_lr.sh

	echo "" > $outDir/transcript_selection/braker.noRNA.gffread.gff # dummy variable for later

fi

if [[ $(ls -A $outDir/RNAseq_alignment | wc -l) -gt 0 && $(ls -A $outDir/ISOseq_alignment | wc -l) -eq 0  ]] ; then

	echo "Only detected short read RNAseq data." 

	bash $tutorialDir/scripts/run_braker_sr.sh

	echo "" > $outDir/transcript_selection/braker.noRNA.gffread.gff # dummy variable for later
	echo "" > $outDir/transcript_selection/braker.lr.gffread.gff # dummy variable for later

fi

if [[ $(ls -A $outDir/RNAseq_alignment | wc -l) -eq 0 && $(ls -A $outDir/ISOseq_alignment | wc -l) -gt 0  ]] ; then

	echo "Only detected ISOseq data." 

	bash $tutorialDir/scripts/run_braker_lr.sh

	echo "" > $outDir/transcript_selection/braker.noRNA.gffread.gff # dummy variable for later
	echo "" > $outDir/transcript_selection/braker.sr.gffread.gff # dummy variable for later

fi
