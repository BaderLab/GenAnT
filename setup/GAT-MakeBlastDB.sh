#!/bin/bash

inDir=$1 # path to "GenomeAnnotationTutorial"

cd $inDir/data/uniprot_sprot

seqkit rmdup -s < uniprot_sprot.fasta > uniprot_sprot_nodup.fasta

makeblastdb -in uniprot_sprot_nodup.fasta -dbtype prot -out uniprot_sprot_nodup.fasta -title "UniPro Sprot  database without duplicated sequences" -parse_seqids 

cd ../Rfam

seqkit rmdup -s < Rfam.fa > Rfam_nodup.fa

makeblastdb -in Rfam_nodup.fa -dbtype nucl -out Rfam_nodup -title "Rfam database without duplicated sequences" -parse_seqids

