#!/bin/bash
#$ -l h_vmem=14G,h_rt=10:00:00,h_stack=32M
#$ -pe smp 16

cd $outDir

mkdir -p ncRNA_analysis ; cd ncRNA_analysis

rfamDir=$dataDir/Rfam


blastn -db $rfamDir/Rfam_nodup \
-query $outDir/assembly/assembly.fa \
-evalue 1e-2 -max_hsps 6 -outfmt 6 \
-num_threads 16 \
-out assembly.rfam.blastn

awk -F "\t" '{print $1 "\t" $7 "\t" $8 "\t" $2}' assembly.rfam.blastn > assembly.rfam.bed

grep -P "\tncRNA\t" $outDir/transcript_selection/mikado_lenient.gff > mikado_annotation.ncRNA.gff

gffread mikado_annotation.ncRNA.gff --bed | cut -f1-4 > mikado_annotation.ncRNA.bed

grep "gbkey=ncRNA" $outDir/transcript_selection/liftoff.gffread.gff | grep -P "RNA\t" > liftoff_annotation.ncRNA.gff

grep "gbkey=ncRNA" $outDir/transcript_selection/custom.gffread.gff | grep -P "RNA\t" > custom_annotation.ncRNA.gff

bedtools subtract -A -a $outDir/transcript_selection/stringtie.gffread.gff -b  $outDir/transcript_selection/mikado_lenient.gff | grep -P "\ttranscript\t" > stringtie_annotation.ncRNA.gtf

gffread stringtie_annotation.ncRNA.gtf --bed | cut -f1-4 > stringtie_annotation.ncRNA.bed

grep -i "RNA" $outDir/assembly/assembly.filteredRepeats.bed | cut -f1-4 > earlGrey_annotation.ncRNA.bed


cat assembly.rfam.bed \
 mikado_annotation.ncRNA.bed \
 liftoff_annotation.ncRNA.bed \
 stringtie_annotation.ncRNA.bed \
 earlGrey_annotation.ncRNA.bed > assembly_ncRNA_seed.bed

bedtools sort -i assembly_ncRNA_seed.bed > assembly_ncRNA_seed.s.bed

bedtools merge -i assembly_ncRNA_seed.s.bed > assembly_ncRNA_seed.m.bed

# awk -F'\t' '!seen[$1 FS $2 FS $3]++' assembly_ncRNA_seed.s.bed | sponge assembly_ncRNA_seed.s.bed

bedtools getfasta -fi $outDir/assembly/assembly.fa -bed assembly_ncRNA_seed.m.bed -fo assembly_ncRNA_seed.m.fasta

cmscan --cpu 16 -Z 1 --cut_ga --rfam --nohmmonly --tblout assembly.tblout \
-o assembly.cmscan --verbose --fmt 2 \
--clanin $rfamDir/Rfam.clanin $rfamDir/Rfam.cm assembly_ncRNA_seed.m.fasta

