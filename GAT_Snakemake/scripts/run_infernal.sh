#!/bin/bash

rfamDir=$dataDir/Rfam

cd $outDir

mkdir -p blast_seed

cd blast_seed

bedtools getfasta -fi $outDir/assembly/assembly.softmasked.fa -bed assembly.rfam.bed -fo ncRNA_full_seed.fa

cmscan --cpu 32 -Z 1 --cut_ga --rfam --nohmmonly --tblout assembly.tblout -o assembly.cmscan --verbose --fmt 2 --clanin $rfamDir/Rfam.clanin $rfamDir/Rfam.cm ncRNA_full_seed.fa

perl $externalDir/infernal-tblout2gff.pl --cmscan --fmt2 --desc assembly.tblout > assembly.infernal.gff3

