#!/bin/bash


outDir=$1
threads=$2

# This code is a slightly more generalized version of https://www.biostars.org/p/92935/ with the purpose for splitting stranded RNA-seq files (https://github.com/Gaius-Augustus/BRAKER/issues/7) for braker3

cd $outDir/RNAseq_alignment
mkdir -p $outDir/braker_sr

braker_out=$outDir/braker_sr

for i in *.bam ; do samtools index $i ; done


# Get the bam file from the command line
# DATA=$1
# prefix=`$i .sorted.bam`

samtools merge -f $braker_out/merged.bam *.bam
samtools index $braker_out/merged.bam

DATA=$braker_out/merged.bam

# Forward strand.
#
# 1. alignments of the second in pair if they map to the forward strand
# 2. alignments of the first in pair if they map to the reverse  strand
#
samtools view -@ $threads -b -f 128 -F 16 $DATA > $braker_out/fwd1.bam
samtools index $braker_out/fwd1.bam

samtools view -@ $threads -b -f 80 $DATA > $braker_out/fwd2.bam
samtools index $braker_out/fwd2.bam

#
# Combine alignments that originate on the forward strand.
#
samtools merge -@ $threads -f $braker_out/fwd.bam $braker_out/fwd1.bam $braker_out/fwd2.bam
samtools index $braker_out/fwd.bam

# Reverse strand
#
# 1. alignments of the second in pair if they map to the reverse strand
# 2. alignments of the first in pair if they map to the forward strand
#
samtools view -@ $threads -b -f 144 $DATA > $braker_out/rev1.bam
samtools index $braker_out/rev1.bam

samtools view -@ $threads -b -f 64 -F 16 $DATA > $braker_out/rev2.bam
samtools index $braker_out/rev2.bam

#
# Combine alignments that originate on the reverse strand.
#
samtools merge -@ $threads -f $braker_out/rev.bam $braker_out/rev1.bam $braker_out/rev2.bam
samtools index $braker_out/rev.bam

rm $braker_out/rev1.bam
rm $braker_out/rev2.bam
rm $braker_out/fwd1.bam
rm $braker_out/fwd2.bam
rm $braker_out/merged.bam
