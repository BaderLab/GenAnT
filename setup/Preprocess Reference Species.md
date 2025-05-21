---
title: "InstallAndDownload"
author: "Dustin Sokolowski and Zoe Clarke"
date: "2025-03-10"
output: html_document
---

# Building a reference directory

We use a reference genome three times in our tutorial. We use it to find gene models with LiftOff, to find gene models with TOGA, and to annotate gene symbols with OrthoFinder. The workflow below describes how to preprocess a reference genome from NCBI (RefSeq) or from Ensembl so that the assembly is compatable with each tool in our tutorial. These files may also be directly included in the config for our pipeline. We use the mouse GRC39 assembly for this example.

## General fasta and gff format:
While below build reference directories for NCBI and Ensembl, any reference assembly and annotation will do. For each internal tool to work, please follow the guidelines below:
1. Do not have spaces and special characters (e.g., "|") in the geneID, transcriptID or chromosome name. TOGA and braker will both fail. TOGA fails outright and Braker "cleans" the genome assembly to converting spaces to underscores, meaning your RNAseq data wont match your genome assembly anymore.
2. For liftoff and TOGA, your gff file will be pre-processed with gffread. Make sure that the gene ID and transcript ID within the gffread processed gff file matches the gene ID and transcript ID that you want. It can sometimes get confused. If it does not match, simply load the gffread.gff file into R with rtracklayer and switch the column out for the ID's that you'd like.

## From Refseq

First, travel to the NCBI FTP containing your species and annotation of interest. For the mouse genome, you would travel to the following URL:

https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39

In `~/data/references` make the directory where you will build the reference species:

```
mkdir -p mmus_GRC39 ; cd mmus_GRC39
```

Download the softmasked assembly and annotation (GFF) for that species. These are used for LiftOff, TOGA, and OrthoFinder.

Download assembly:

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz
```

Download GFF:

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz
```

Unzip each file

```
for i in *.gz ; do gunzip $i ; echo $i ; done
```

#### Running Preprocessing Script

While we explain each step below, you can also perform this full analysis with our premade script:
`preprocess_reference_from_refseq.sh`. We also expect the `annotation_tutorial` environment to be active.

This script has four positional arguments:

- working directory (i.e., path to human_T2T_NCBI)
- /scripts directory (i.e., path to the scripts in this tutorial)
- genomic FASTA
- genomic GFF

With our example, the script would be executed as shown below.

Here `~` is the path to the cloned GenomeAnnotationTutorial. In our case it would be:
`/.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/test`

The directory: `~GenomeAnnotationTutorial/data/references/mmus_GRC39` has two files:
`GCF_000001635.27_GRCm39_genomic.fna` and `GCF_000001635.27_GRCm39_genomic.gff`

```
  bash preprocess_reference_from_refseq.sh \
  ~/data/references/human_T2T_NCBI \
  ~ \
  GCF_000001635.27_GRCm39_genomic.fna \
  GCF_000001635.27_GRCm39_genomic.gff
```

#### Inside Preprocessing Script 

The script itself is here, and each line is explained:

Preprocess the downloaded gff file with `GFFRead`. GFF files allow for some flexibility in format, even in NCBI and Ensembl (e.g., custom annotations in model organisms), and `GFFRead` ensures their compatibility with LiftOff and TOGA.

```
gffread GCF_000001635.27_GRCm39_genomic.gff --keep-genes -o GCF_000001635.27_GRCm39_genomic.gffread.gff
```

Generate an amino acid FASTA from CDS regions (for OrthoFinder)

```
gffread -y GCF_000001635.27_GRCm39_genomic.protein.faa -g GCF_000001635.27_GRCm39_genomic.fna GCF_000001635.27_GRCm39_genomic.gffread.gff
```

Filter for pseudogenes containing premature stop codons. This is not always needed but the `.nostop.protein.faa` file will be more generalizable.

```
sed '/^>/!s/[.*]//g' GCF_000001635.27_GRCm39_genomic.protein.faa > GCF_000001635.27_GRCm39_genomic.nostop.protein.faa 
```

Covert GFF file into a BED12 file compatible with TOGA.

```
~/external/kent/gff3ToGenePred GCF_009914755.1_T2T-CHM13v2.0_genomic.gffread.gff GCF_009914755.1_T2T-CHM13v2.0_genomic.genePred

~/external/kent/genePredToBed GCF_009914755.1_T2T-CHM13v2.0_genomic.genePred GCF_009914755.1_T2T-CHM13v2.0_genomic.bed
```

Generate gene-keys.

* An "isoforms.tsv" file for TOGA. This file is a key between gene ID and each transcript ID.
* A "table.txt" file for OrthoFinder. This file is a key between the predicted protein FASTA (FAA) file heading and gene symbol.
* A "genekey.txt" file. This is a key between gene ID, transcript ID, and gene symbol.

We do this using R.

```
library(rtracklayer)

# Read in the GFF file using rtracklayer
gff <- readGFF("GCF_009914755.1_T2T-CHM13v2.0_genomic.gff")

gff_transcript <- gff_transcript <- gff[gff$type == "transcript",]

# Get gene ID and transcript ID keys
# if this is a character already it doesn't change anything, but Rtracklayer sometimes loads "Parent" in as a list.

gff_transcript$gene_id <- unlist(gff_transcript$Parent)
gff_transcript$trans_id <- unlist(gff_transcript$ID)

key <- gff_transcript[,c("gene_id","trans_id","gene")]
colnames(key) <- c("geneID","transcriptID","geneSymbol")

write.table(key[,c("geneSymbol","transcriptID")],file="GCF_009914755.1_T2T-CHM13v2.0_genomic.table.txt",
quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

write.table(key[,c("geneID","transcriptID")],file="GCF_009914755.1_T2T-CHM13v2.0_genomic.isoforms.tsv",
quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

write.table(key,file="GCF_009914755.1_T2T-CHM13v2.0_genomic.genekkey.tsv",
quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
```

The other is a key between gene name and protein ID.

The final step (and it's only needed for specific GFF files) is to remove transcripts that are not compatible with kentutils. In our practice this is specific custom annotations in model organisms (e.g., grc39 -mouse-, dm6 -drosophila melanogaster- genomes.

```
awk 'FNR==NR { if (NR > 1) keep[$2]; next } $4 in keep' GCF_009914755.1_T2T-CHM13v2.0_genomic.isoforms.tsv GCF_009914755.1_T2T-CHM13v2.0_genomic.bed > GCF_009914755.1_T2T-CHM13v2.0_genomic.toga.bed

```


## From ENSEMBL

First travel to the ENSEMBL FTP containing your species and annotation of interest.

https://useast.ensembl.org/info/data/ftp/index.html

In `~/data/references` make the directory where you will build the reference species:

```
mkdir -p mmus_GRC39_embl ; cd mmus_GRC39_embl
```

Download the softmasked assembly and annotation (GFF) for that species. These files are used for LiftOff, TOGA, and OrthoFinder.

Download assembly. Don't forget to use the assembly with the"sm" option, as we need a softmasked assembly for TOGA.

```
wget https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz
```

Download GFF:

```
wget https://ftp.ensembl.org/pub/release-113/gff3/mus_musculus/Mus_musculus.GRCm39.113.gff3.gz
```

Unzip each file

```
for i in *.gz ; do gunzip $i ; echo $i ; done
```

#### Running Preprocessing Script

The EMBL pre-processing script is identical to the RefSeq script - the exception that it uses a different R script to account for Ensembl-specific GFF file headers.

While we explain each step below, you can also perform this full analysis with our premade script:
`reference_directory_ensembl.sh`. 
This script has six positional arguments:

- working directory (i.e., path to human_T2T_NCBI)
- /scripts directory (i.e., path to the scripts in this tutorial)
- genomic FASTA
- genomic GFF

reference_directory_ensembl.sh 

With our example, the script would be executed as:

```
  bash reference_directory_ensembl.sh \
  ~/data/references/mmus_GRC39_embl \
  ~ \
  Mus_musculus.GRCm39.dna_sm.primary_assembly.fa \
  Mus_musculus.GRCm39.113.gff3
```
  
This script and workflow is nearly identical to ncbi. The few differences are.

* Ensembl has the `gff3` suffix instead of `gff`
* Ensembl uses `Name` to denote gene symbol, while ncbi uses `gene`
* Ensembl tags gene and transcript with `gene:` and `transcript:` while Refseq uses `gene-` and `rna-`.
* Ensembl includes transcript number in gene symbol, NCBI does not.
