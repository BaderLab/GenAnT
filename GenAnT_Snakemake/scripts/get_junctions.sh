#!/bin/bash

outDir=$1

cd $outDir/transcript_selection


if [[ $(ls -A $outDir/RNAseq_alignment | wc -l) -eq 0 && $(ls -A $outDir/ISOseq_alignment | wc -l) -eq 0 ]] ; then

	echo "We did not find short read or long read RNA-seq data (aligned bam files)"
	echo "" > junctions.final.bed
fi

if [[ $(ls -A $outDir/RNAseq_alignment | wc -l) -gt 0 && $(ls -A $outDir/ISOseq_alignment | wc -l) -gt 0 ]] ; then

	echo "We detected short read and long read RNAseq data"


	# Run Portcullis (shor read RNAseq data)

	RNA_aln=$outDir/RNAseq_alignment
	ASSEMBLY=$outDir/assembly/assembly.softmasked.fa

	samtools merge $outDir/transcript_selection/merged.sr.bam $RNA_aln/*merged.bam
	samtools index $outDir/transcript_selection/merged.sr.bam

	portcullis full -t 20 $ASSEMBLY $outDir/transcript_selection/merged.sr.bam -o $outDir/transcript_selection


	ISO_aln=$outDir/ISOseq_alignment
	samtools merge $outDir/transcript_selection/merged.lr.bam $ISO_aln/*merged.bam
	samtools index $outDir/transcript_selection/merged.lr.bam

	regtools junctions extract -s XS -o $outDir/transcript_selection/lr_junctions.bed $outDir/transcript_selection/merged.lr.bam

	cat $outDir/transcript_selection/3-filt/portcullis_filtered.pass.junctions.bed $outDir/transcript_selection/lr_junctions.bed > junctions.final.us.bed

	bedtools sort -i junctions.final.us.bed > junctions.final.bed

fi

if [[ $(ls -A $outDir/RNAseq_alignment | wc -l) -gt 0 && $(ls -A $outDir/ISOseq_alignment | wc -l) -eq 0  ]] ; then

	echo "Only detected short read RNAseq data." 


	# Run Portcullis (shor read RNAseq data)

	RNA_aln=$outDir/RNAseq_alignment
	ASSEMBLY=$outDir/assembly/assembly.softmasked.fa

	samtools merge $outDir/transcript_selection/merged.sr.bam $RNA_aln/*merged.bam

	portcullis full -t 20 $ASSEMBLY $outDir/transcript_selection/merged.sr.bam -o $outDir/transcript_selection

	mv $outDir/transcript_selection/3-filt/portcullis_filtered.pass.junctions.bed $outDir/transcript_selection/junctions.final.bed

fi

if [[ $(ls -A $outDir/RNAseq_alignment | wc -l) -eq 0 && $(ls -A $outDir/ISOseq_alignment | wc -l) -gt 0  ]] ; then

	echo "Only detected ISOseq data." 

	ISO_aln=$outDir/ISOseq_alignment
	samtools merge $outDir/transcript_selection/merged.lr.bam $ISO_aln/*merged.bam

	regtools junctions extract -s XS -o $outDir/transcript_selection/junctions.final.bed $outDir/transcript_selection/merged.lr.bam
fi

