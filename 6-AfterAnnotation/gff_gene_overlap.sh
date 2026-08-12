#!/bin/bash

##
### Helper script to intersect gene symbols between a reference gff and target gff 
### The target gff (typically the mikado gff) and a reference gff (e.g., custom GFF),
### and the output bed file. This will map up the mikado ID to the best matching gene ID in the reference GFF
### get_gene_overlap.sh mikado.gff custom.gff custom_mikado.tsv 
##

target=$1
ref=$2
out=$3

# 1. Extract gene features only, convert to BED (0-based, half-open)
awk 'BEGIN{OFS="\t"} $3=="gene" {
    match($9, /ID=([^;]+)/, id); match($9, /gene_name=([^;]+)/, name)
    gname = (name[1] != "") ? name[1] : id[1]
    print $1, $4-1, $5, gname, ".", $7
}' ${ref} | sort -k1,1 -k2,2n > genes1.bed

awk 'BEGIN{OFS="\t"} $3=="gene" {
    match($9, /ID=([^;]+)/, id); match($9, /gene_name=([^;]+)/, name)
    gname = (name[1] != "") ? name[1] : id[1]
    print $1, $4-1, $5, gname, ".", $7
}' ${target} | sort -k1,1 -k2,2n > genes2.bed

# 2. Intersect, keeping full info from both + overlap length
bedtools intersect -a genes1.bed -b genes2.bed -wao > overlap_matches.tsv

sort -k10,10 -k13,13nr overlap_matches.tsv | awk '!seen[$10]++' > ${out}.tsv
