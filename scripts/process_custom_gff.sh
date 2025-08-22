#!/bin/bash

$outDir

if [[ "$customGFF" != "none" ]] ; then 
	echo "We detected a custom gff file"

	gffread $customGFF --keep-genes -o $outDir/transcript_selection/custom.gffread.gff
else

	echo "We did not detect a custom gff"
	echo "" > $outDir/transcript_selection/custom.gffread.gff

fi
