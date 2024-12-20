# i. Installation

### i.a Conda evnironment
  This tutorial requires a series of conda packages. All of which are found on 
  bioconda or conda-forge. These conda packages are listed below. Installing
  this combined environment can be performed using the following yml file:
  
```
  conda env create -f annotation_tutorial.yml

```
  
conda packages: earlgrey, liftoff, infernal, blast, busco, transdecoder, diamond,
bedtools, mirmachine, nextflow, orthofinder



### i.b Singularity images

  This tutorial uses two singularity images: one for `mikado` and one for 
  `braker`. The below code should sucessfully download `braker3.sif` and 
  `mikado.2.3.5rc2.sif`.
  
We expect that these binaries are all installed within a directory called 
  `external`.
  
```
# Obtain braker3.sif file
singularity build braker3.sif docker://teambraker/braker3:latest

# Obtain mikado.2.3.5rc2.sif file in mikado_singularity/cache/mikado.2.3.5rc2.sif

VERSION=2.3.5rc2
singularity exec docker://gemygk/mikado:v${VERSION} mikado -h

```

### i.c Binaries

  This tutorial uses four package binaries to complete the tutorial. Namely,
  it uses binaries for `TOGA`,  `kent` utulities, `stringtie`, and 
  `interproscan`.
  
  We expect that these binaries are all installed within subdirectories
  of a directory called  `external`. (see `i.e` data and package structure)
  
  To successfully run TOGA, we need two binaries. One is the `CACTUS` aligner
  and the other is the `TOGA` binary itself. Because we are only looking at 
  two species pairwise alignments, `CACTUS` only needs to be run in "local" 
  mode. As a result, the basic installation of CACTUS without additional 
  formatting for different schedulers (PBS, SLURM, SGE etc.). As such, we
  recommend installing and configuring CACTUS separately from this tutorial 
  if you want to use CACTUS for 3+ species multiple sequence alignments, as 
  the basic installation with no changes may not play nicely with your 
  scheduler.
  
  
```

# Assuming you are in ~/external, where "~" is the path to your annotation pipeline

# Download and install CACTUS aligner. 

wget https://github.com/ComparativeGenomicsToolkit/cactus/releases/download/v2.8.4/cactus-bin-v2.8.4.tar.gz
tar -xzf cactus-bin-v2.8.4.tar.gz
cd cactus-bin-v2.8.4

# download and install TOGA
git clone https://github.com/hillerlab/TOGA.git
cd TOGA
python3 -m pip install -r requirements.txt --user
./configure.sh
./run_test.sh micro



```

The `kent` binaries are universally valuable packages used for converting
genomic file types to one another. Accordingly, a number of these binaries are
required in this tutorial.

```

# Assuming you are in ~/external, where "~" is the path to your annotation pipeline

mkdir -p kent

cd kent

# Binaries for generating chain files for TOGA

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chainSwap

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/axtChain

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslPosTarget

# Binaries for converting TOGA outputs, and making reference annotations for TOGA

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToGenePred
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/genePredToBed
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/gtfToGenePred
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gff3ToGenePred

# some clusters may require an extra step to ensure that each binary is executable
for i in `ls` ; do chmod +x $i ; done


```

The `stringtie` binary is required for individuals with RNA-seq data. 
Many different installation sources can be found at their website
(https://ccb.jhu.edu/software/stringtie/), however we found the easiet approach
was to clone their github repository.

```
  git clone https://github.com/gpertea/stringtie
  cd stringtie
  make release
```

Finally, interproscan is used for annotating functional domains for proteins

```

mkdir my_interproscan
cd my_interproscan
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.69-101.0/interproscan-5.69-101.0-64-bit.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.69-101.0/interproscan-5.69-101.0-64-bit.tar.gz.md5

# Recommended checksum to confirm the download was successful:
md5sum -c interproscan-5.69-101.0-64-bit.tar.gz.md5

# Must return *interproscan-5.69-101.0-64-bit.tar.gz: OK*
# If not - try downloading the file again as it may be a corrupted copy

```
  
### i.d data packages

Genome annotation heavily relies on comparing aligned sequences between a reference
dataset (e.g., another species assembly and annotation, publicly avalable
protein sequences etc.) and your target species. As such, a considerable amount
of data needs to be downloaded for this pipeline to work. We store data that
we expect to be useful in most genome annotation pipelines. We provide 
instructions for downloading the equivalent datasets for different clades
in  `ii. utility and preparation`. 

We expect your working directory is `~/data`, where ~ is the path where 
you store the data and packages for genome annotation.

#### uniprot_sprot
The uniprot_sprot database provides experimentally validated protein sequences
to cross-reference candidate transcripts against. This datbase is used in 
Step 3: transcript selection, and Step 5: gene symbol calling. 

```
mkdir -p uniprot_sprot

cd uniprot_sprot

wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

gunzip uniprot_sprot.fasta.gz

```

#### Rfam database
The Rfam database is required for finding ncRNA genes with known secondary 
structures and for identifying them into their RNA family. 

```

mkdir -p Rfam

cd Rfam

wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/Rfam.fa
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz


```

#### Reference genomes and annotations.


We assume that reference genome and annotations are in a directory called 
`Ref`, found in `~/data`

```
mkdir -p Ref # directory containing reference genomes.

```

TOGA allows for homologous gene transfer between a a quality reference genome
and annotation to a target genome within the same order (e.g., Rodentia). For
example, mouse genes successfuly transfer to the naked mole-rat, despite 
an estimated 70 million years diverged. 

TOGA requires a soft masked assembly fasta file, a genome annotation file where
transcripts are extracted and converted to bed12 format, and a key matching
gene IDs to transcript IDs.

LiftOff allows for a more comprehensive (e.g., non-coding genes) homologous 
gene transfer between animals in the same genus (e.g., Marmotta). LiftOff
relies on minimap alignments and therefore does not expect gaps, resulting 
in more stringent requirements for seequence similarity.

LiftOff requires an assembly fasta file and a genome annotation file in the
traditional gtf format.

The authors of TOGA have pre-processed reference annotations for mouse (mm10)
and humans (hg38). The code below pulls the data needed for mouse and human 
annotations. Processing of data required for other species to be the reference
genome (or an updated mouse and human reference) 
can be found in `ii. utility and preparation`.

Toga inputs - assuming you already sucessfully installed TOGA.

#### mouse mm10

```

mkdir -p Ref/mm10
cd Ref/mm10

# you can probably link here as well but I have found that linked directories
# do not always play nicely with bound directories in singularity images in my 
# testing. TOGA is not run in a singularity but I figured it's worth keeping this open.
cp ~/external/TOGA/supply/mm10.*toga* ./

# Getting the assembly from UCSC, which is where TOGA took the mm10 assembly from.
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip mm10.fa.gz

# Get assembly gtf file for liftOff
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.knownGene.gtf.gz

gunzip mm10.knownGene.gtf.gz

```


#### Human hg38

```

mkdir -p Ref/hg38
cd Ref/hg38


cp ~/external/TOGA/supply/hg38.*toga* ./

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.fa.gz
gunzip hg38.fa.gz

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz
gunzip hg38.knownGene.gtf.gz

```

#### Orthofinder Reference

The information for a reference species for orthofinder, our primary tool
dedicated to finding gene symbols, is a protein fasta, translated cds fasta,
and a key of the gene symbol and the protein ID for a species. The protein fasta
and translated cds fasta can be obtained from RefSeq. For example, the the files 
required for mouse (GRCm39) would be in their RefSeq ftp
`https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/`

We are assuming you are in ~/data

```

mkdir -p mmus_GRC39

cd mmus_GRC39

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_protein.faa.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_translated_cds.faa.gz

for i in *.gz ; do gunzip $i ; done

```

There is a script in `utils` that generates the symbol - protein ID table required for orthofinder.

Using the same mouse example as above

```

cd mmus_GRC39

bash ~/utils/refseq_gene_protin_table.sh GCF_000001635.27_GRCm39_translated_cds.faa GCF_000001635_table.txt

```

#### BUSCO Databases

We use BUSCO to evaluate the total number of universal single-copy orthologs (BUSCO)
in a genome assembly, and to identify how many of those BUSCO genes are caputured
in each annotation. BUSCOdbs are specific to clades. For example, the buscodb
for new world rodents (`glires_odb10`), may find more BUSCO genes 
than for mammalia overall `mammalia_odb10`. The code below stores
`mammalia` in default as that applies to the species that this tutorial is 
built for.

Any may be installed from `https://busco-data.ezlab.org/v4/data/lineages/`.

```

mkdir -p buscodb

cd buscodb


wget https://busco-data.ezlab.org/v4/data/lineages/mammalia_odb10.2020-09-10.tar.gz

tar -zxvf *tar.gz -C ./

```

### i.e data and package structure.

Our tutorial expects a specific data structure for two directores: `./data` 
and `./external`. 

#### Data directory structure. 

`./data` stores the reference databases, reference genomes, and example files
to run this tutorial. Below lists the expected directories and their components

* `example_data`: This directory contains the files required to run this tutorial.
To save on run time and memory, The assembly that we use is already softmaksed by Earl Grey
(see Step 1). We also use RNA-seq data from two tissues from Bens et al. 2018 (PMID: 30068345).
Scripts to aid in processing these RNA-seq data from fastq files are found in `utils`. We use the 
`hisat` tool for RNA-seq alignment and processing. Other aligners (e.g., `STAR`) may be used but
typically have additional parameters to be compatible with the transcriptome assembler, StringTie.
Again, to save on runtime and memory, we are only annotate chromosmome 28 of the NMR, and focusing 
on two tissues where reads have already been filtered to only align to chromosome 28.
  * NMRchr28.fa. The soft masked genome assembly for NMR chromosome 28 (assembled using PacBio Hifi data).
  * thyroid.chr28.merged.bam. RNA-seq data from Bens et al. 2018 (PMID: 30068345) aligned to chromosome 28. Thyroid tissue.
    * skin.chr28.merged.bam. RNA-seq data from Bens et al. 2018 (PMID: 30068345) aligned to chromosome 28. Skin tissue.
    * NMRchr28.filteredRepeats.bed. Bed file providing repeat annotations for NMRchr28.fa. This file is an output of `Step 1`, but is copied here as it is one line of code and saves > 48h of runtime.
    
    
* `Rfam`: the downloaded Rfam fasta and cm files, compressed and converted into
a blastdb. Scripts to generate these files are in `utils`. 
* `uniprot_sprot`: the downloaded protein fasta file from the uniprot_sprot database
and converted into a blastdb. Scripts to generate these files are in `utils`. 

* `references`: This directory contains a series of directories for each assembly.Note, we use a wildcard between prefix (mm10) and filetype (e.g., .fa), so these directories should not store additional files (e.g., multiple assemblies, multiple annotations etc.).
We have mm10 and hg38 (the assemblies with pre-processed data from `TOGA`) already stored.
  * Each directory requires a number of files. We will use `mm10` as an example:
    * mm10*.fa: A softmasked genome assembly of the reference species.
    * mm10*.fa.fai: The fasta index for the same reference assembly (generated with `samtools faidx mm10.fa`).
    * mm10*.gtf: The reference genome annotation. This anotation requires the `gene` type, so pulling the gtf from RefSeq will not work.
    * mm10*.for_toga.bed: The bed12 required for a genome annotation with TOGA. 
    This annotation can be created for a custom annotation with `utils`.
    * mm10*.for_toga.isoforms.tsv: a key showing a key between the geneID 
    and transcriptID. 
    This annotation can be created for a custom annotation with `utils`.
* `orthofinder_ref`: This directory contains the reference transcripts for orthofinder,
our primary method for identifying gene symbols. Note that we use wildcards so do not include extra files in each directory.
 * Each directory requires a number of files. We will use `mmus` as an example:
    * mmus.*_protein.faa: the protein fasta file from the reference annotation.
    * mmus.*_translated_cds.faa: the transated cds fasta file from the reference annotation.
    * mmus*_table.txt: a key between gene symbols and protein ID. This file is
    generated using a custom `awk` script in `utils`.
* `buscodb`: a directory that contains a busco database that has been downloaded and unzipped by the user.
The database can be updated to the most approprirate linear for your species, however we have the `mammalia_odb10`
as a default. We expect these folders to be downloaded and formatted in the manner described in `BUSCO Databases`.



#### External directory structure. 

Our pipeline integrates pre-existing tools into a continuous workflow with a small number of custom
scripts to allow for the compatibility between tools and data formats. For this tutorial to run properly,
we expect a singularity image for braker and mikado, as well as binaries for stringtie, toga, cactus, and kent-utils. Details on how to download all of these packages are found in `i.b Singularity images`

`./external` expects the following data directories.

* `braker_singularity`: A directory containing `braker3.sif`. Installation can include other test scripts etc. but the `.sif` file is the only mandatory item in this directory.
* `mikado_singularity`: Installation of the mikado singularity as described in `i.b` puts the .sif image in a cache subdirectory. As a result, we expect `mikado_singularity` to contain `cache/mikado.2.3.5rc2.sif`.
* `cactus-bin-v2.8.4`: A directory with a built and configured binary for the CACTUS aligner. Following instructions in `i.c Binaries` generates this directory.
* `kent`: A directory with executible binaries for the kent utilities. Following instructions in `i.c Binaries` generates this directory.
* `TOGA`: A directory with executible binaries for the TOGA annotation tool. Following instructions in `i.c Binaries` generates this directory.
* `stringtie`: A directory with executible binaries for the stringtie transcript assembler. Following instructions in `i.c Binaries` generates this directory.s

# ii. Utilities

In this tutorial, we expect some pre-processing of publicly available databases (e.g., the uniprot-sprot protein sequence database into a blastdb). We also expect RNA-seq data aligned as sorted bam files (and merged within each tissue). While not required, we also expect that users may use a different reference species for `TOGA` than human or mouse. Lastly, a small amount of custom-preprocessing is required to generate the required input tables for orthofinder, and a small amount of pre-processing may be required to generate transcriptID-genesymbol keys to annotate gene symbols using LiftOff and TOGA. The `utils` directory contains scripts to aid with pre-processing of these data before running the tutorial.

### Required before running the tutorial for the first time.

* `rfam_preprocess.sh`: This script takes the `Rfam.fa` file downloaded in `i`, removes redundant sequences, generates a `blastdb`, and processes covariance models of RNA models to be compatible with infernal. This script should be run after the `Rfam.fa` file is downloaded the first time, and does not need to be rerun unless the fasta file is re-downloaded (e.g., due to an Rfam update).

* `uniprot_preprocess.sh`: This script takes the `uniprot_sprot_nodup.fasta` downloaded in `i`, removes redundant sequences, generates a `blastdb`. Like the `Rfam.fa` database,  This script should be run after the `uniprot_sprot.fasta` file is downloaded the first time, and does not need to be rerun unless the fasta file is re-downloaded (e.g., due to an uniprot update).

### RNA sequencing alignment information.

These scripts are not required but help users make `bam` files that are compatible with this tutorial. If you choose to use your own bam files, we expect the following:
* bam files are merged for each tissue.
* bam files are compatible with the `stringtie` transcript assembler tool. For example, if you choose to use the `STAR` aligner, it needs the `--outSAMstrandField intronMotif` flag. If you choose to use the scripts in `utils`, then reads with be aligned with `hisat` using minimal filtering. The following scripts are run. 

Assuming that data is not aligned, the below scripts serve as a guideline.
* `prep_assembly_for_annotation.sh`. To align RNA-seq data with `hisat2`, an assembly index needs to be detectable. This is a very simple script where `$assembly` is the path the the assembly fasta file and `prefix` is a short-hand name for the assembly you will use in your annotation.

```
samtools faidx $assembly
hisat2-build $assembly $prefix
```

* `RNAseq_alignment.sh`. This script aligns and post-processes paried end RNA-seq data with hisat2.
For a paired end RNA-seq file, we align reads with hisat2, filter the reads with `samtools view`, and then sort the reads with `samtools sort`. The script is submitted for each fastq file. The important elements of the script is below.

```

hisat2 -x $assemblyDir/$prefix -1 $rnaDir/$b"1.fastq.gz" -2 $rnaDir/$b"2.fastq.gz" -S $b".sam"

samtools view -bSq 20 -f 0x2 $b".sam" > $b".bam"

samtools sort -o $b".sorted.bam" $b".sam"

samtools index $b".sorted.bam"

```

* `RNAseq_merge.sh`. This script is a companion for `RNAseq_alignment` and `prep_assembly_for_annotation`. This script uses `samtools merge` to combine bam files originating from the same tissue. This script uses wildcards to combine sorted bam files.

```

samtools merge -o $tissue.merged.bam *$tissue*.sorted.bam


```

### Prep annotation for TOGA

The below script uses a combination of the binaries in the `kent` utilities and a custom R script to generate the `bed12` annotation and `isoform.tsv` files required for `TOGA`. This is only used if you plan to use an assembly that is not `mm10` or `hg38` for TOGA annotations.

The `make_TOGA_ref.sh` script requires prefix for the species gtf file `gtfref`. We expect the gtf file to be a `$assembly.genomic.gtf` downloaded from NCBI. The R script uses the `rtracklayer` package to load in the same gtf file and make a key between the `gene ID` and `transcript ID`, which is then used in `TOGA` annotations.

The most important elements of the script are.

```

gtfToGenePred $gtfref $prefix.genePred -ignoreGroupsWithoutExons
genePredToBed $prefix.genePred $prefix.bed

Rscript --vanilla utils/NCBI_gtf_to_TOGA_input.R

```

### Prep table for Orthofinder

The below script is required to make a key between protein and transcript tables for Orthofinder. It is only required if you are identifying orthologous gene symbols and the reference species is not human or mouse.

* `refseq_gene_protin_table.sh` Is an `awk` script that extracts and maps the transcript and protein IDs in an NCBI annotated assembly to one another. It expects two arguments. The first being the Translated CDS file from RefSeq and the second being the  Name of output TSV file.

# iii. config.yaml

We have set up this tutorial so that every script pulls from the config file built into the tutorial. Therefore, if the config file is correctly set up, then the only requirement is to run the bash scripts stored in the tutorial. Generating this file should be the most time-consuming aspect of the tutorial.

Most of the config file revolves around setting the directories where your reference and target files are, as well as the file prefixes. Some tools use prebuilt models that are tuned to different clades (e.g., Earl Grey, Mirmachine, Busco). The config file sets these prebuilt models as `mammalia`, but changing these parameters may help improve annotations (e.g., Glires models improve naked mole-rat miRNA models). 

### The YAML file

Below is the default example of the yaml file:

```
### CHECK configPath=/.mounts/labs/simpsonlab/users/dsokolowski/miniconda3/envs/genome_eval/config ###


# Configuration file for Annotation pipeline
#
required_directories:
  # Directories required for annotation pipeline
  outDir: /.mounts/labs/simpsonlab/users/dsokolowski/projects/T2T_assembly/annotations/assembly
  # Directory where directory structure will be generated from.
  condaDir: /.mounts/labs/simpsonlab/users/dsokolowski/miniconda3/envs/vertebrate_annotation_tutorial
  # Directory where conda environment is saved to
  sourceDir: /.mounts/labs/simpsonlab/users/dsokolowski/miniconda3/bin/activate
  # Directory for activating conda environment
  externalDir: /.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/external
  # Path to directory to programs that were not downloaded with conda environment
  dataDir: /.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/data/
  # Path to directory with data files used in annotation pipeline
  togaRef: mm10
  # Reference directry must contain assembly fasta file (of the same name, e.g., mm10.fa), A qc.bed file and an qc.isoforms.txt file denoting transcripts required for TOGA (see /utils to generate annotation files for species not in ~/data/references).
  # If you do cannot use an additional reference assembly using LiftOff. Then this species must also include the "genomic.gtf" file for this assembly.
  # We link assemblies and annotations with *fa, *bed, and *txt, so having multiple fasta and gtf files in the assembly may cause issues.
optional_directories:
  EarlGrey_complete: "no"
  # Path to $species"_summaryFiles". We copy masked assembly and filtered annotation.bed files for later analysis
  rnaseqDir: NULL
  # Path to directory containing bam files from paired-end RNA-seq data. Scripts to generate bam files from fastq files using hisat2 are found in /utils.
  # We expect the suffix of each file to be "*sorted.bam"
  isoseqDir: NULL
  # Path to directory containing bam files from paired-end RNA-seq data. Scripts to generate bam files from fastq files using XXX are found in /utils.
  liftoffRef: NULL
  # Path to directory with two files: the assembly fasta (e.g., ref_squirrel.fa) and the genomic annotation (e.g., reg_squirrel.gtf). We link assemblies and annotations with *fa and *gtf, so having multiple fasta and gtf files in the assembly may cause issues.
  singularity_cache: NULL
  # Path to directory for singularity-cache. If NULL, we will build the temporary directory in the same working directories where the singularity is being used (e.g., in transcript_selection).
  singularity_tmpdir: NULL
  # Path to directory for singularity-tmp. If NULL, we will build the temporary directory in the same working directories where the singularity is being used (e.g., in transcript_selection).
required_files_options:
  # Names of required (and optional) input files
  prefix: "T2T_primary"
  # Name of species prefix, this will be used in output files
  assemblyfile: /.mounts/labs/simpsonlab/users/dsokolowski/projects/T2T_assembly/Pup3-NMR.curated.reordered.mt.fa
  # assembly to be annotated (with path). This will be copied into "assembly.fa" file or "assembly.softmasked.fa" if masked=True
filenames_optional:
  # Optional files for higher quality annotations
  Liftoff_annotation:
options:
  EarlGrey:
    masked: False
    # If masked = True, then earl-grey is skipped and 
    # Parameters specifically for Earl Grey
    species: "heterocephalus_glaber"
    # Species name for earl grey
    RepeatDb: "Rodentia"
  buscoDb: mammalia_odb10
    # If user wants to use a different BUSCOdb for genome evaluation (e.g., glires_odb10 for new world rodents), download it into ~/external/busco_db/ and change the name of the directory in `buscoDb: mammalia_odb10`
  MirMachine:
    node: mammal
    # If user wants to use a more specific family (Glires) to capture more miRNA then switch family. Use `MirMachine.py --print-all-nodes` to check available families,
  max_threads: 50
    # Maximum number of threads available for pipeline


```

### clade-specific parameters.

* mirmachine. Mirmachine detects miRNAs in an assembly and has prebuilt models based on clade. `MirMachine.py --print-all-nodes` lets you see which clades are compatible with MirMachine. The most compatible clade should be put in the `MirMachine: node:` section of the config.yaml.
* RepeatDb. Specifying the species or clade can improve repeat annotations. Many, many options are available. For the NMR and Marmot we used Rodentia
* species. The name of your species of interest with a txid. The full species name in lowercase also works (e.g., heterocephalus_glaber). This parameter does not influence the output of earl grey but just directory structure.
* buscoDb. Changed to whatever buscoDb is a better fit for your species. No need to change if you use mammalia.
  
  
# Tutorial flow. 

Assuming code is correctly installed, the `/data` and `/external` directories are correctly strucutred and the `config.yaml` is correctly ordered, the rest of this genome annotation pipeline occurs by executing bash scripts.

The steps below describe each script in more detail but this script sumamrizes what lines of code needs to be run and when.

```
# Step 1: Reeat annotation.

qsub run_earl_grey.sh # note, takes ~150 hours for a mouse genome, may be easier to run separately.

# Step 2: Generating protein-coding gene models. Each annotation tool in step 2 can be run independently.

## LiftOff

qsub run_liftOff.sh

## TOGA

qsub run_toga_wrapper.sh # ~72 hours for mouse --> NMR transfer

## Braker

qsub run_braker.sh  # > 72 hours for NMR.

  # If you do not have RNA-seq data, run_braker.sh will call run_braker_noRNA.sh, which implements braker protein (braker2)
  
## StringTie 

  # If you have RNAseq data
for i in *.merged.bam ; do qsub run_stringtie_only.sh -v I=$i ; echo $i ; done

  # once the stringtie step is complete

qsub stringtie_mergeonly.sh

  # Flagging no RNAseq data in the config.yaml file will skip this step and change Step 3 into the `noRNAseq` path.

## Step 3: Combining and filtering gene models

  # Combining and filtering gene models involves setting up and implementing the mikado transcript selector.

get_gffread.sh # extract gtf files from each annotation.

  # Flagging noRNAseq data will run get_gffread_noRNAseq.sh. This excludes processing the stringtie outputs.

qsub 12_mikado_configure_and_prepare.sh # prepare mikado yaml file and prepare combined gtf before filtering.

# once the previous line is complete, run:

qsub run_portcullis.sh # This step is skipped without RNAseq, and junction detection is skipped.

qsub run_transdecoder.sh

qsub sub_mikado_blast.sh # splits the transcript fasta file into 100 files and submits run_looped_mikado_blast.sh
  # this submits 100 jobs as default, which runs the pipeline in ~1.5 hours assuming all 100 jobs can run simulaneously. 
  
# Once portcullis transdecoder and mikado_blast are all finished.

qsub 45_mikado_serialize_pick.sh

## Step 4: Annotating non-coding RNA genes

qsub blast_seed.sh # blast assembly against the Rfam database

qsub ncRNA_combine_infernal.sh # once seeding is done build ncRNA covariance models.
  # NOTE TO SELF: MAKE THE FIRST FEW LINES FLEXIBLE FOR REPEATS

Rscript --vanilla infernal_biotype_map.R # Map biotype to the final gtf file still must write

qsub run_mirmachine.sh  

Rscript --vanilla reintegrate_ncRNA.R # All gene models are present, make the final gtf file.
  
## Step 5: Sequence-similarity-based transfer of gene symbols

# gene symbols are anottated using TOGA and LiftOff while gene models are transferred. The final step is orthofinder.

qsub run_orthofinder.sh

Rscript --vanilla gene_symbol_table.R

Rscript --vanilla populate_final_gtf.R
  
  
```

# Step 1. Repeat Annotation

Repeat annotation involves identifying the simple repeats, satelites, and transposable elements in a genome assembly. This process has been streamlined by the tool Earl Grey and now uses a single line of code.

```
earlGrey -g $assemblyDir/$assembly -s $species -o . -t 50 -r $RepeatDb -d yes
```

Assuming you are working in this tutorial and are using the `config.yaml` file, all you need to do is submit the `run_earl_grey.sh` job.
* Note. Earl Grey can take a long time. For example, for the NMR, it took 150 hours for a 45% repetitve genome at 32 cores and over 200 Gb. of RAM.

Assuming you are working in this tutorial and are inputting a masked genome (e.g., you ran Earl Grey already, you masked using another tool, then we expect two outputs).

* A softMasked fasta file (like the file used in this tutorial)
* A bed file providing the co-ordinates for repeats in columns 1-3, and the repeat ID in the fourth column. These annotations can help in non-coding RNA gene finding later in the tutorial.

In this instance, `run_earl_grey.sh` simply copies your masked assembly and bed file into the `/assembly` directory that the tutorial reads future files out of.

Lastly, while beyond the scope of this tutorial. `Earl Grey` outputs tonnes of useful suff, particularly in the `speices_EarlGrey/species_summaryFiles/` directory. We reccomend checking out these files to get a comprehensive breakdown to the repeat strucutre of your genome assembly (with publication quality figures).

# Step 2: Generating protein-coding gene models

In this tutorial, we integrate gene models predicted by four tools, `liftOff`, `TOGA`, `StringTie`, and `Braker`.
To generate the final set of protein coding gene models in your species. If you have a genome annotation from a close relative, `liftOff` will likely do a good job of transferring many non-coding gene models as well. Similarly `StringTie` generates gene models using RNA-seq data, meaning that it will simultaneously also identify candidate lncRNA gene models that are evaluated in Step 3.

If you do not have a closely-related species to perform homologous gene transfer with `liftOff`, then `liftOff` will be performed using the same reference genome as TOGA. Most of these models will be filtered out in `Step 3`, but we have found that the few models that it adds are worth keeping overall.

If you do not have any RNA-seq data, then StringTie is skipped. Braker is also run in Braker-Protein (BrakerP, Braker2) mode rather than Braker3 (using RNA-seq and Protein evidence).

### 2a) liftOff

LiftOff is very easy to implement and transferring gene models between closely related rodents takes less that 12 hours on a desktop. The main line of code is below. The `-copies` parameter helps transfer duplicated gene models between assemblies (e.g., the target species has a tandem duplication while the reference does not). LiftOff also runs without a masked fasta file.

```
liftoff -g $refGTF $assembly $refFa -o ./$prefix.liftoff.gtf -u ./$prefix.unmapped.liftoff.txt -copies

```

Using our tutorial + config file, `liftOff` is run with `run_liftOff.sh`.

### 2b) TOGA

TOGA is run in three stages listed below. TOGA also expects both the reference and target assembly fasta's to be soft masked. 

#### 2.b.1) prepare files
Format the reference annotation if neccesary (see ``make_TOGA_ref.sh` in `/utils/`)

The pairwise genome alignment and TOGA scripts can be run in serial with `run_toga_wrapper.sh`

#### 2.b.2) pairwise genome alignment
Align the reference to target species using the `CACTUS` alginer and generate a `chain` file (i.e., a pariwse alignment file that allows for simultaneous gaps). Because we are only generating an alignment between two species, our script runs `CACTUS` in single_machine (i.e., local, CPU...) mode, rather than with any flow control. We found this much easier than having each user configure cactus for their own cluster, and because of the runtime of Braker (which occurs simultanously to TOGA), barely impacts the total runtime of the tutorial. This held true for large genome such as the grey whale. 

Using the tutorial + config, this alignment is performed with:

* `make_cactus_tree.sh`. A simple text processing script that generates an input text file for CACTUS from the species names, directories, and fasta files stored in `config.yaml`

* `cactus_aln_and_chain.sh`. Generate an alignment .hal file with the cactus aligner and then convert the hal file to a chain file with the instructions provided in https://github.com/ComparativeGenomicsToolkit.
* The main line to run the cactus aligner is:

```
cactus jobstore cactus_in.txt target_ref.hal --binariesMode local --realTimeLogging True --batchSystem single_machine --workDir $outdir/cactus_aln
```
#### 2.b.3) transfer annotations
Once steps 1 and 2 are run, running TOGA is completed with the `toga.py` script. An example using mm`0 as a reference is below:

```
toga.py $outDir/cactus_aln/target-to-ref.chain $supplyDir/mm10.wgEncodeGencodeCompVM25.bed $outDir/cactus_aln/$refPrefix.2bit $outDir/cactus_aln/$prefix.2bit --project_name $prefix"_toga" --nextflow_config_dir $externalDir/TOGA/nextflow_config_files/ --cesar_binary $externalDir/TOGA/CESAR2.0/cesar --isoforms $supplyDir/mm10.wgEncodeGencodeCompVM25.isoforms.txt
```

In our pipeline, this is generated with `run_toga.sh`.

Note. Like `Earl Grey`, TOGA has some helpful standalone utility. Three files that we always find useful are listed below, but all of the tables in the output hold some utility.

* loss_summ_data.tsv: prediction of transcript activity and intactness in your species. This file helped find genes involved in the evolution of blindness and sperm motility in naked mole-rats.
* inact_mut_data.txt: a detailed list of what elements of each gene sucessfully and unsuccessfullly transfered between species. An example of a gene from mouse --> NMR is below:
```
GENE: ENSMUST00000075303
# ENSMUST00000075303	204	1	0	Deleted exon	-	masked	DEL_1
# ENSMUST00000075303	204	3	123	ATG	ATG->ATG	masked	ATG_2
# ENSMUST00000075303	204	3	134	ATG	ATG->ATG	masked	ATG_3
# ENSMUST00000075303	204	4	157	ATG	ATA->ATG	masked	ATG_4
# ENSMUST00000075303	204	4	169	ATG	ATG->ATG	masked	ATG_5
# ENSMUST00000075303	204	6	236	ATG	ATG->ATG	masked	ATG_6
# ENSMUST00000075303	204	6	241	ATG	ATG->ATG	masked	ATG_7
# ENSMUST00000075303	204	6	246	ATG	ATG->ATG	masked	ATG_8
# ENSMUST00000075303	204	7	265	ATG	ATA->ATG	masked	ATG_9
# ENSMUST00000075303	204	7	276	ATG	ATG->ATG	masked	ATG_10
# ENSMUST00000075303	204	8	315	ATG	ATG->ATG	masked	ATG_11
# ENSMUST00000075303	204	8	330	ATG	ATG->ATG	masked	ATG_12
# ENSMUST00000075303	204	8	409	ATG	GTG->ATG	masked	ATG_13
# ENSMUST00000075303	204	INTACT_PERC_IGNORE_M 0.8833333333333333
# ENSMUST00000075303	204	INTACT_PERC_INTACT_M 0.8833333333333333
# ENSMUST00000075303	204	INTACT_CODONS_PROP 0.8833333333333333
# ENSMUST00000075303	204	OUT_OF_CHAIN_PROP 0.0
# ENSMUST00000075303	204	MIDDLE_IS_INTACT TRUE
# ENSMUST00000075303	204	MIDDLE_IS_PRESENT TRUE
```
* `orthology_classification.tsv`: helps dileate orthologs and paralogs in your annotation.

### 2C) StringTie

StringTie is run in three stages listed below, with no masking required. This stage is skipped if no RNA-seq data is required.

#### 2.C.1) prepare files
This pipeline expects Merged bam files for each tissue, with scripts to generate these data in `/utils`

#### 2.C.2) Run stringtie

For each merged bamfile, generate a de novo annotation using stringtie.

`for i in *.merged.bam ; do qsub run_stringtie_only.sh -v I=$i ; echo $i ; done`

#### 2.C.3) Merge stringtie models.


After `StringTie` is run, we merge annotations from each tissue: `RNAseq_merge.sh`

StringTie merge takes the the longest and most complete isoforms from these tissues, which filters some short tissue-specific isoforms. Our testing found that applying StringTie merge rather than naively concatenating gtf generated gene models that were closer to those generated by `NCBI` and `ENSEMBL`.


### 2D) Braker

StringTie is run in three stages listed below, with no masking required. This stage is skipped if no RNA-seq data is required.

Braker uses a combination of ab initio gene prediction with predicted gene models from publicly available protein data and RNA-seq data. If RNA-seq data is available, this tutorial uses Braker3 to integrate all datatypes to generate high quality gene models. Otherwise, we use braker2 to integrate ab initio gene prediction with protein data. The computational challenge is getting braker installed, but the singularity issue fixes most of these issues. 
Given the installation and proper config file, braker sucessfully runs with: `run_braker.sh`

# Step 3: Combining and filtering gene models

We combine the gene models from Step 2 using the mikado transcript selector. Briefly, mikado prepares a fasta and gtf file of potential transcripts. It then collects feature information for each transcript by analyzing the GFF files (e.g., number of exons, presenece of UTRs, start/stop codon, etc.) and by interpreting the output of different bioinformatic tools. In this tutorial we use `transdecoder` to evaluate open reading frames, `portcullis` to evaluate CDS regions and splice junctions (RNA-seq required), and `blastp` to evaluate the similarity between a predicted CDS and existing protein sequences in the uniprot swiss-prot database. Mikado then uses this list of features for each potential transcript to pick the most likely transcripts for each loci, thereby establishing gene models for your assemby.

#### gtf processing: 
Mikado expects each genome annotation to be in a gff file with relatitely consistent formatting. The `gffread` tool processes these data in one line of code. Here is an example of braker below.

`gffread $brakerDir/braker.gtf --keep-genes -o braker.gffread.gff`

We use `get_gffread.sh` to post-process gene models from each tool.

#### mikado configure and prepare

We then programatically make the text file that mikado converts into their input yaml filewith a custom R script called `make_mikado_input.R`. We then run `mikado configure` and `mikado prepare` with the transcript scoring parameters trained on mammals `mammalian.yaml`. 

Postprocessig, mikado configure, and mikado prepare are all completed with: `12_mikado_configure_and_prepare.sh `

#### generate transcript features.

Three tools, `portcullis`, `transdecoder`, and `blastp` are used to evaluate junctions and CDS', open reading frames, and protein homology respectfully. `portcullis` relies on RNA-seq data and is therefore skipped if users lack RNA-seq data.

Portcullis is run with default parameters and is executed in this pipeline with: `run_portcullis.sh`
`portcullis full -t 20 $ASSEMBLY merged.bam -o $outDir/transcript_selection`

Transdecoder is run with default parameters and is executed in this pipeline with: `run_transdecoder.sh`
```
TransDecoder.LongOrfs -t $outDir/transcript_selection/mikado_prepared.fasta -m 100 -O "transdecoder"
TransDecoder.Predict -t $outDir/transcript_selection/mikado_prepared.fasta --retain_long_orfs_length 100 -O "transdecoder"
```
The `-m` parameter influences what ORF length should be considered an mRNA gene vs. lncRNA gene with a tiny ORF. we used the default `-m 100` from existing literature but changing this number drastically impacts the number of transcripts that are given the ncRNA or mRNA I.D. from transcripts derived from `StringTie`.

Blastp run with default parameters and is executed in this pipeline with: `sub_mikado_blast.sh`.
This script splits the mikado transcript fasta file into 100 components, and submits each sub-fasta as it's own blastp run, thereby increasing runtime 100 fold (assuming all jobs can be run simultaneously). The script run on each fasta is `run_looped_mikado_blast.sh`, and runs blastx against the uniprot database.

```
blastx -max_target_seqs 5 -num_threads 20 -query "$QUERY_FILE" \
 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop" \
 -db "$SPROT_DB" -evalue 0.000001 -out $OUTPUT_DIR/blast_results.tsv

```

#### serialize and picks

The tools to compute transcript phenotypes have been completed. There information are then parsed and organized into a coherent table and mikado database using `mikado serialize`, before these data are used to select the highest quality transcript for each loci using `mikado pick`.

In our pipeline, `mikado serialize` and `mikado pick` are run with `45_mikado_serialize_pick.sh`.

The key lines of this script are:

```
singularity exec --bind ${outDir},${wd},${ASSEMBLYDIR},${UNIDB} ${MIKADO_SIF} mikado serialise \
 -p 20 --start-method spawn \
 --orfs $outDir"/transcript_selection/transdecoder/transdecoder/mikado_prepared.fasta.transdecoder.bed" \
 --transcripts $outDir"/transcript_selection/mikado_prepared.fasta" \
 --tsv $outDir"/transcript_selection/blast/blast_results.tsv" \
 --json-conf $outDir"/transcript_selection/config.yaml" \
 --genome_fai $outDir"/assembly/assembly.softmasked.fa.fai" \
 -od $outDir"/transcript_selection" \
 --log $outDir"/transcript_selection/mikado_serialise.log" \
 --blast-targets $UNIDB"/uniprot_sprot_nodup.fasta" \
 --max-target-seqs 5 \
 --junctions $outDir"/transcript_selection/3-filt/portcullis_filtered.pass.junctions.bed"
 
 singularity exec --bind ${outDir},${wd},${ASSEMBLYDIR},${UNIDB} ${MIKADO_SIF} mikado pick \
 --json-conf $outDir"/transcript_selection/config.yaml" \
 -db $outDir"/transcript_selection/mikado.db" \
 --mode lenient \
$outDir"/transcript_selection/mikado_prepared.gtf" \
 --log $outDir"/transcript_selection/mikado_pick.log"  \
 --loci-out $outDir"/transcript_selection/mikado_lenient.gff" \
 --fasta $outDir"/assembly/assembly.softmasked.fa" \
 --no-purge


```

# Step 4: Annotating non-coding RNA genes

At this stage, we have annotated all transcripts with the potential to generate an mRNA, and have parsed these transcripts into lprotein coding genes and lncRNA genes. As such, we have yet to annotate short RNA genes (miRNA, tRNA, snRNA etc.) that could not be a transcript.

The Rfam database contains covariance models for thousands of RNA gene families. By analyzing alignments and predicted secondary strucutres against covariance models generated in the Rfam database, one could, in principle, use the entire assembly as a query space for ncRNA genes. For a mammalian genome however, this would likely take months. Instead, we first "seed" the genome assembly, which is described as identifying regons of the assembly that contain any evidence of an ncRNA gene being present there. These sequences are then tested against the Rfam database. Skipping seeding rarerly, if ever, identifies extra ncRNA genes and takes thousands of times longer.

#### 4.A) Seeding:

Seeding is pedominantly performed with a blastn analysis of the entire assembly against the Rfam database.
This is peroformed with `blast_seed.sh`.

The main line of code is: 

```
blastn -db $rfamDir/Rfam_nodup -query $outDir/assembly/$assembly -evalue 1e-6 -max_hsps 6 -max_target_seqs 6 -outfmt 6 -out assembly.Rfam.blastn -num_threads 24

awk -F "\t" '{print $1 "\t" $7 "\t" $8 "\t" $2}' assembly.Rfam.blastn > assembly.rfam.bed
```

Some ncRNA seeding is also completed as a biproduct of Steps 1-3. Other seeds arise from:
  * TE-derived tRNA and sRNA genes. This is why we need `NMRchr28.filteredRepeats.bed` in the tutorial.
  * ncRNA genes models predicted in mikado.
  * ncRNA genes transferred from a high quality reference genome with LiftOff.

We combine the bed files generated from these three sources with the Rfam blast hits found in `blast_seed.sh`. This merged bed file is the seeding used for covariance matrix evaluation. 

This integration of bed files is the first half of the `ncRNA_combine_infernal.sh` script.

#### 4.B) Evaluation:

With seeding completed, we use the Infernal tool to evaluate whether a sequence is derived from a ncRNA gene from a known RNA family. This analysis is perfromed with `ncRNA_combine_infernal.sh`

The command in this script that performs these evaluations is:

```
cmscan --cpu 32 -Z 1 --cut_ga --rfam --nohmmonly --tblout assembly.tblout -o assembly.cmscan --verbose --fmt 2 --clanin $rfamDir/Rfam.clanin $rfamDir/Rfam.cm ncRNA_full_seed.fa

```
  
Lastly, we use a key between RNA family and biotype from Rfam to assign the biotype for each RNA gene. (e.g., map the S6 family to the snRNA biotype). 

```
Rscript --vanilla infernal_biotype_map.R # Map biotype to the final gtf file still must write
                                         # Would this be added at the final stage to make the final gtf?
                                          
```

MirMachine contains clade-specific miRNA covariance models that, in our testing, has the potential to identify 2x more bonadife miRNA genes than using Infernal against Rfam alone. Your specific precomputed clade can be found with ``MirMachine.py --print-all-nodes` (see `/utils`) and included in the config.yaml file.

Exeucute MirMachine on our tutorial with `qsub run_mirmachine.sh`. The main line of code is:

```
MirMachine.py -n Mammalia -s $species --genome $ASSEMBLY -m deutero --cpu 20
```

##### GFF integration

At this stage all mRNA, lncRNA, and short ncRNA models have been identified, but they are spread across three gtf files (mikado_lenient.gff, infernal.gtf, and miRNA.gtf). We have written an R script that ensures that gene IDs and symbols from infernal.gtf and miRNA.gtf.

```
Rscript --vanilla reintegrate_ncRNA.R # All gene models are present, make the final gtf file.
```


# Step 5: Sequence-similarity-based transfer of gene symbols

Gene symbols are anottated using TOGA and LiftOff while gene models are transferred. The final step is orthofinder.

#### Orthofinder

Orthofinder is built to infer ortholgos and paralogs, estimate gene duplication events, and infer conservation across many species. In this tutorial, we limit Orthofinder to two species to allow for flexibility and computational efficiency.

Unless you are calling orthologs between mouse and human, orthofinder requires some preprocessing of the reference species (see Orthofinder Reference). Mouse and human are pre-downloaded for this tutorial. 

In this tutorial we execute `qsub run_orthofinder.sh`

The main line of code executed here is:

```
orthofinder -t 20 -a 20 -o mikado_grc39 -f protein_seqs
```

Like other pipelines used in this tutorial (e.g., Earl Grey, TOGA), Orthofinder returns a number of useful files that are not all used in this specific pipeline. For example, `Comparative_Genomics_Statistics/Statistics_Overall.tsv` gives a summary of the genes that found orthologues. As well as `Comparative_Genomics_Statistics/OrthologuesStats_one-to-one.tsv`, which should show the number of high confidence gene symbols that will come from this analysis.

The table used to identify protein coding gene symbols is in: `orthofinder/mikado_ref/Results_*/Orthologues/Orthologues_*_protein/*_protein__v__mikado_lenient.tsv` where the wildcard is the reference species.

#### Compile remaining gene symbols


Additional gene symbols are pulled from the output of LiftOff and of TOGA. We extract the mRNA types from the output of `LiftOff` and `TOGA` and use `bedtools intersect` with our gene models to identify which `LiftOff`- and `TOGA`- derived gene symbols match our integrated gene models. Then, we use the custom R script, `findGeneSymbols.R` to extract gene symbols from the output of orthofinder. Lastly, we pull these data files together and combine them into a single table, where rows are gene IDs and columns are the predicted gene symbol from Orthofinder, TOGA, LiftOff, and Infernal (for ncRNA) using `gene_symbol_table.R`.

The post-processing scripts described above are run with: `qsub combine_gene_symbols.sh`

#### Add gene symbols:

Lastly, we use a custom gtf script to populate the annotation gff with the gene symbols.

Assuming you have a annotations from a close relative (i.e., the full use of LiftOff), then gene symbols are prioritized as such.

  * If multiple tools agree, then that symbol is chosen. 
  * Otherwise: LiftOff > TOGA > Orthofinder
  * Otherwise (no full use of LiftOff): TOGA > OrthoFinder > LiftOff

`Rscript --vanilla populate_final_gtf.R`

## Horray! You have annotated a mammalian genome

## Step X

There are a few extra steps that may help annotate assemblies in certain cases. These steps are currently under development.

### lncRNA gene symbol identification.

Many lncRNA genes are well characterized and have well studied gene symbols. Performing a recipricol `blastn` analysis of your lncRNA genes against a `blastdb` of these well-characterized lncRNA genes (e.g., `XIST`) could help add some high quality gene symbols. These symbols will likely be transferred with `LiftOff` assuming you can use an annotation of a related species.

### TCR, IGG, and MHC gene annotation

These gene families exist in clusters that can be easily identified after completing our annotation, because a number of genes from each cluster will be identified and provided a gene symbol (assuming the genome assembly resolves these regions). Therefore, by observing these annotations on a genome browser, once could pull out the gene family clusters.

Extracting these clusters, converting them into a translated fasta format (via an ORF finder), and BLASTing them against known domains for these `TCR, IGG, and MHC` may help find some other gene models. This is only true for the V and sometimes C domains of these genes. J domains require single-cell pacbio sequencing, as they are <40nt sequences found on the same mRNA as the V gene. 
