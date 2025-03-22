#!/bin/bash

orthofinderData=$dataDir/orthofinder_ref/$orthofinderRef

cd $outDir

mkdir -p orthofinder

cd orthofinder

ln -s $outDir/transcript_selection/mikado_lenient.gff ./

echo "GFF3 file: mikado_lenient.gff"
echo "Name of gene list out: mikado_lenient.txt"

grep -P "\tgene\t" mikado_lenient.gff > mikado_lenient.txt
sed -i 's/.*ID=\(.*\);Name.*/\1/' mikado_lenient.txt

ln -s $orthofinderData/*txt ./

mkdir -p protein_seqs

cd protein_seqs 

ln -s $outDir/transcript_selection/mikado_lenient.faa ./
ln -s $orthofinderData/*protein.faa ./

cd ../

orthofinder -t 20 -a 20 -o $target"_orthofinder" -f protein_seqs
