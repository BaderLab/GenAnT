#!/bin/bash

cd $outDir

mkdir -p ncRNA_analysis ; cd ncRNA_analysis

rfamDir=$dataDir/Rfam

# convert to gff format
perl $externalDir/infernal-tblout2gff.pl --cmscan --fmt2 assembly.tblout > infernal.gff

# copy Rfam family information
cp $dataDir/Rfam/family.txt ./

# reformat the GFF so that the feature types are more recognizable (i.e., add bioid)
Rscript --vanilla $scriptsDir/scripts/rfamConversion.R

# isolate the lncRNA features 
grep -P '\tlncRNA\t' infernal.types.gff > infernal.types.lncRNA.gff

bedtools intersect -v -a infernal.types.lncRNA.gff -b $outDir/transcript_selection/mikado_lenient.gff > infernal.lncRNA.notInMikado.gff

# If any gene models have populated infernal.lncRNA.notInMikado.gff, append them to mikado.gff.

cat $outDir/transcript_selection/mikado_lenient.gff infernal.lncRNA.notInMikado.gff > mikado.infernal.gff

# label anything in mikado.infernal.gff with information from Infernal.

bedtools intersect -a mikado.infernal.gff -b infernal.types.lncRNA.gff -wo > mikado.infernal.lncRNALabeled.txt

cut -f1-9 mikado.infernal.lncRNALabeled.txt > mikado.infernal.lncRNALabeled.mikadoInfo.gff
cut -f10-18 mikado.infernal.lncRNALabeled.txt > mikado.infernal.lncRNALabeled.infernalInfo.gff


Rscript --vanilla $scriptsDir/scripts/RenamelncRNAs.R

# Integrate these newly formatted lncRNAs with the rest of Mikado's gene models by first subtracting the lncRNA features from the original Mikado features. 

bedtools subtract -A -a mikado.infernal.gff -b mikado.infernal.lncRNALabeled.polished.gff > mikado.noLnc.gff

cat mikado.noLnc.gff mikado.infernal.lncRNALabeled.polished.gff > mikado.lncLabeled.gff

grep -P "\texon\t" mikado.lncLabeled.gff > mikado.lncLabeled.exons.gff

grep -v -P "\tlncRNA\t" infernal.types.gff > infernal.types.noLncRNA.gff



cp $outDir/mirmachine/results/predictions/filtered_gff/$species.PRE.gff ./mirmachine.gff

sed 's/gene_id/ID/g' mirmachine.gff > mirmachine.id.gff
sed -i 's/sequence_with_30nt.*/gbkey=ncRNA/g' mirmachine.id.gff


bedtools subtract -A -a infernal.types.noLncRNA.gff -b mirmachine.id.gff > infernal.noLnc.noMir.gff

cat mirmachine.id.gff infernal.noLnc.noMir.gff > short_ncRNAs.gff

bedtools subtract -A -a short_ncRNAs.gff -b mikado.lncLabeled.exons.gff > short_ncRNAs.noOverlap.gff

sed -i 's/E-value/evalue/g' short_ncRNAs.noOverlap.gff

Rscript --vanilla $scriptsDir/scripts/AddFeaturesNcRNAs.R

cat mikado.lncLabeled.gff short_ncRNAs.polished.gff > full_annotation.unsorted.gff

bedtools sort -i full_annotation.unsorted.gff > full_annotation.gff

