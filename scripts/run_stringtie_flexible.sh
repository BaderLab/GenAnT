#!/bin/bash

mkdir -p $outDir/stringtie_out

# No RNA-seq data
if [[ $(ls -A $outDir/RNAseq_alignment | wc -l) -eq 0 && $(ls -A $outDir/ISOseq_alignment | wc -l) -eq 0 ]] ; then

	echo "We did not find short read or long read RNA-seq data (aligned bam files)"
	echo "" > $outDir/transcript_selection/stringtie.gffread.gff # dummy variable for later

fi

if [[ $(ls -A $outDir/RNAseq_alignment | wc -l) -gt 0 && $(ls -A $outDir/ISOseq_alignment | wc -l) -gt 0 ]] ; then

	echo "We detected short read and long read RNAseq data, will perform stringite short, long, and mixed depending on data availability for each tissue"

	cd $outDir

	Rscript --vanilla $tutorialDir/scripts/organize_stringtie_runs.R # find exact matchng bam files in the ISOseq_alignment and RNAseq_alignment directories

	cd $outDir/stringtie_out

	if [[ $(wc -l < mixed.txt) -gt 0 ]] ; then

		while IFS= read -r line; do
			echo "$line"
			i=$line
			b=`basename $i .bam`
			$externalDir/stringtie/stringtie --mix $outDir/RNAseq_alignment/$i $outDir/ISOseq_alignment/$i -l $b -o $outDir/stringtie_out/$i".mix.gtf" -p 8 --conservative

		done < mixed.txt

	fi


	if [[ $(wc -l < isoseq_only.txt) -gt 0 ]] ; then

		while IFS= read -r line; do
   		 	echo "$line"
   		 	i=$line
   			b=`basename $i .bam`

   			$externalDir/stringtie/stringtie $outDir/ISOseq_alignment/$i -l $b -L -o $outDir/stringtie_out/$i".lr.gtf" -p 8 --conservative

		done < isoseq_only.txt

	fi

		if [[ $(wc -l < rnaseq_only.txt) -gt 0 ]] ; then  

		while IFS= read -r line; do # run stringtie for the tissues with short read RNA-seq

			echo "$line"
			i=$line

			b=`basename $i .bam`

			$externalDir/stringtie/stringtie $outDir/RNAseq_alignment/$i -l $b -o $outDir/stringtie_out/$i".sr.gtf" -p 8 --conservative

		done < rnaseq_only.txt

	fi

	cd $outDir

	$externalDir/stringtie/stringtie --merge -o $outDir/stringtie_out/stringtie.merged.gtf $outDir/stringtie_out/*gtf # merge results

	gffread $outDir/stringtie_out/stringtie.merged.gtf --keep-genes -o $outDir/transcript_selection/stringtie.gffread.gff # make compatible with downstreat steps




fi

if [[ $(ls -A $outDir/RNAseq_alignment | wc -l) -gt 0 && $(ls -A $outDir/ISOseq_alignment | wc -l) -eq 0  ]] ; then

	echo "Only detected short read RNAseq data." 

	cd $outDir/RNAseq_alignment
	b=`basename $i .bam`

	for i in *.bam ; do $externalDir/stringtie/stringtie $i -l $b -o $outDir/stringtie_out/$i".gtf" -p 8 --conservative ; done

	$externalDir/stringtie/stringtie --merge -o $outDir/stringtie_out/stringtie.merged.gtf $outDir/stringtie_out/*gtf

	gffread $outDir/stringtie_out/stringtie.merged.gtf --keep-genes -o $outDir/transcript_selection/stringtie.gffread.gff



fi

if [[ $(ls -A $outDir/RNAseq_alignment | wc -l) -eq 0 && $(ls -A $outDir/ISOseq_alignment | wc -l) -gt 0  ]] ; then

	echo "Only detected long read RNAseq data (ISOseq)."
	cd $outDir/RNAseq_alignment
	b=`basename $i .bam`

	for i in *.bam ; do $externalDir/stringtie/stringtie $i -l $b -o $outDir/stringtie_out/$i".gtf" -p 8 --conservative ; done

	$externalDir/stringtie/stringtie --merge -o $outDir/stringtie_out/stringtie.merged.gtf $outDir/stringtie_out/*gtf

	gffread $outDir/stringtie_out/stringtie.merged.gtf --keep-genes -o $outDir/transcript_selection/stringtie.gffread.gff

fi
