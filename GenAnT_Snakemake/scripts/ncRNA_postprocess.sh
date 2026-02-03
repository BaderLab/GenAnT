#!/bin/bash

outDir=$1
dataDir=$2
externalDir=$3
snakeDir=$4

cd $outDir

mkdir -p ncRNA_analysis ; cd ncRNA_analysis

rfamDir=$dataDir/Rfam

# convert to gff format
perl $externalDir/infernal-tblout2gff.pl --cmscan --fmt2 assembly.tblout > infernal.gff

# copy Rfam family information
cp $dataDir/Rfam/family.txt ./

# reformat the GFF so that the feature types are more recognizable (i.e., add bioid)
Rscript --vanilla $snakeDir/scripts/rfamConversion.R

# isolate the lncRNA features 
grep -P '\tlncRNA\t' infernal.types.gff > infernal.types.lncRNA.gff

bedtools intersect -v -a infernal.types.lncRNA.gff -b $outDir/transcript_selection/mikado_lenient.gff -s > infernal.lncRNA.notInMikado.gff

gffread -F infernal.lncRNA.notInMikado.gff -o infernal.lncRNA.New.gff

bedtools intersect -a infernal.types.lncRNA.gff -b $outDir/transcript_selection/mikado_lenient.gff -wo -s > infernal.lncRNA.InMikado.gff

cp  $outDir/transcript_selection/mikado_lenient.gff ./

# Extend existing gene models that directly overlap with cmscan ncRNA annoated conserved gene segments



echo "Infernal found lncRNAs that have some overlap with existing mikado genes."
Rscript --vanilla $snakeDir/scripts/ExtendlncRNAs.R


# Add in the lncRNAs with no overlap with any existing gene models

cat mikado.infernal.lncRNALabeled.polished.gff infernal.lncRNA.New.ExonID.gff > mikado.infernal.gff

grep -P "\texon\t" mikado.infernal.gff > mikado.lncLabeled.exons.gff

grep -v -P "\tlncRNA\t" infernal.types.gff > infernal.types.noLncRNA.gff



cp $outDir/mirmachine/mirna.filtered.gff ./mirmachine.gff

sed 's/gene_id/ID/g' mirmachine.gff > mirmachine.id.gff
sed -i 's/sequence_with_30nt.*/gbkey=ncRNA/g' mirmachine.id.gff


bedtools subtract -A -a infernal.types.noLncRNA.gff -b mirmachine.id.gff -s > infernal.noLnc.noMir.gff

cat mirmachine.id.gff infernal.noLnc.noMir.gff > short_ncRNAs.gff

bedtools subtract -A -a short_ncRNAs.gff -b mikado.lncLabeled.exons.gff -s > short_ncRNAs.noOverlap.gff

sed -i 's/E-value/evalue/g' short_ncRNAs.noOverlap.gff

grep -v '^#'  short_ncRNAs.gff > short_ncRNAs.nohead.gff
awk 'BEGIN{OFS="\t"} { # Make consistent and generic names for ncRNA 
    split($9,a,";");
    id=a[1]"."NR;
    $9=id";Name="a[1]";"a[2]";"a[3];
    print
}' short_ncRNAs.nohead.gff | \
awk 'BEGIN{OFS="\t"} {
    if($0 ~ /^#/){ next }          # skip comments
    attr=$9; # Remove tabs in 9th column
    for(i=10;i<=NF;i++){ attr=attr" "$i }
    $9=attr;
    NF=9;
    print
}' > short_ncRNA.unclean.gff

gffread -F short_ncRNA.unclean.gff -o short_ncRNAs.polished.gff 

sed -i 's/RFamID=/rfamid=/g; s/E-value=/evalue=/g; s/Name=ID=/Name=/g; s/ID=ID=/ID=/g' short_ncRNAs.polished.gff



cat mikado.infernal.gff short_ncRNAs.polished.gff > full_annotation.unsorted.gff



sed -i 's/RFamID=/rfamid=/g;  s/CDS_infernal_product/cds_infernal_product/g' full_annotation.unsorted.gff

bedtools sort -i full_annotation.unsorted.gff > full_annotation1.gff

mv full_annotation1.gff full_annotation.gff


