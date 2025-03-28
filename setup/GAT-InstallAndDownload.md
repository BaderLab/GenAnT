---
title: "InstallAndDownload"
author: "Dustin Sokolowski and Zoe Clarke"
date: "2025-03-10"
output: html_document
---

Our mammalian genome annotation tutorial includes (1) a step-by-step tutorial, (2) plug-and-chug scripts, and (3) a Snakemake workflow that facilitates end-to-end genome annotation using a config file. The following installation guide prepares your machine for any one of these workflows. A complete and smooth installation process can be performed by following all of the steps under "Plug-and-Chug", where you set up a Conda environment and then run two scripts that download and install the various tools and create a BLAST database. "Step-by-step" is a more tutorial-like installation process that guides you in detail through the installation process. Choose one to run.

Both the plug-and-chug and Snakemake workflows expect specific directory structure and for specific programs to be installed. The step-by-step tutorial is more flexible.

## Plug-and-Chug

We provide scripts to streamline the installation process.

```
# Clone the tutorial repo.

git clone https://github.com/BaderLab/GenomeAnnotationTutorial.git

cd GenomeAnnotationTutorial

# Create the conda environment from the repo

conda env create -f annotation_tutorial.yml

conda activate annotation_tutorial

# make sure that SINGULARITY_CACHEDIR and SINGULARITY_TMPDIR are not attached to root.
# export SINGULARITY_CACHEDIR=/path-to/singularitycache
# export SINGULARITY_TMPDIR=/path-to/singularitytemp

# run install and download script. This requires the internet.

bash GAT-InstallAndDownload.sh /path-to-GenomeAnnotationTutorial

# build blastDBs

bash GAT-MakeBlastDB.sh /path-to-GenomeAnnotationTutorial
```

It may be easier to edit `GAT-MakeBlastDB.sh` to submit it as a job on your cluster.

Lastly, if your cluster is not `slurm`, you will need to update the nextflow config files described in `##### TOGA config` (or in the TOGA repo itself - https://github.com/hillerlab/toganextflow_config_files)

## Step-by-step

### 1. Installation

#### 1a - Directory structure.

Clone the tutorial repo.

`git clone https://github.com/BaderLab/GenomeAnnotationTutorial.git`

Create a `/data` and a `/external` directory within the `GenomeAnnotationTutorial` directory. `/data` will hold public dataests
and reference genomes, while `/external` will hold binaries and Singularity images
used in the tutorial and workflow.

```
cd GenomeAnnotationTutorial
mkdir -p data
mkdir -p external
```

We will eventually be installing multiple Singularity images. To prepare for this, we need to make sure Singularity will run without error in our environment. Having `$SINGULARITY_CACHEDIR` and `SINGULARITY_TMPDIR` will throw a failure to install images error for many non-admin accounts as their roots typically allow for ~10 or ~100Gb or space. To combat this, we can create new cache and temporary directories for Singularity somewhere in a space owned by the user, and reassign the variables.

```
mkdir -p external/singularity_images
mkdir -p external/singularity_images/singularitycache
mkdir -p external/singularity_images/singularitytemp

export SINGULARITY_CACHEDIR=/path-to/singularitycache
export SINGULARITY_TMPDIR=/path-to/singularitytemp
```

#### 1b - Conda environment

This tutorial requires a series of conda packages. All of which are found on bioconda or conda-forge. These conda packages are listed below. Installing this combined environment can be performed using the following YML file:
  
```
conda env create -f annotation_tutorial.yml
```
  
The YML file is on the GitHub, but can also be copied from here:

```
name: annotation_tutorial  # Name of the environment
channels:             # Channels to search for packages
  - conda-forge
  - bioconda
  - defaults
dependencies:         # List of packages and their versions
  - earlgrey=5.1.0  # check if that's still on bioconda -- blast comes from here otherwise there's an issue
  - liftoff=1.6.3
  - infernal=1.1.2
  - busco=5.7.1
  - transdecoder=5.7.1
  - portcullis=1.2.4
  - minimap2=2.28
  - diamond=2.1.10
  - bedtools=2.31.1
  - mirmachine=0.2.13
  - orthofinder=3.0.1b1
  - hisat2=2.2.1 
  - augustus=3.5.0
  - gffread=0.12.7
  - nextflow=24.10.4
  - regtools=1.0.0
  - seqkit=2.9.0
```
  
#### 1c - Singularity images

This tutorial uses two Singularity images: one for `mikado`, two for `braker` (the short read RNAseq and ISOseq options), and one for the `cactus` aligner.
  

```
cd external/singularity_images
```


Install the relevant Singularity images.

```
singularity build braker3.sif docker://teambraker/braker3:latest

singularity build braker3_lr.sif docker://teambraker/braker3:isoseq

singularity build cactus.v2.9.3.sif docker://quay.io/comparative-genomics-toolkit/cactus:v2.9.3

singularity build mikado_gat.sif docker://risserlin/mikado:ubuntu22_mikado2.3.2
```

Optionally clear the Singularity and temporary cache to save space.

```
rm -r singularitycache/*
rm -r singularitytemp/*
```

#### 1d - Binaries

This tutorial uses four package binaries to complete the tutorial. Namely, it uses binaries for `TOGA`,  `kent` utulities, `stringtie`, and `interproscan`.
  
We expect that these binaries are all installed within subdirectories of a directory called  `external`. (see `i.e` data and package structure)
  
##### TOGA

The documentation of `TOGA` requires that you also install the binary for the Cactus aligner. We elected to use their Singularity image to promote the stability of this tutorial across version upgrades. As such, the `cactus` aligner is already installed at this point.

```
# Navigate back to "external"
cd ..

# download and install TOGA
git clone https://github.com/hillerlab/TOGA.git
cd TOGA
python3 -m pip install -r requirements.txt --user
./configure.sh
```

TOGA should be tested at this point as well. TOGA requires NextFlow, so the test will require a Conda environment with NextFlow (preferable our `annotation_tutorial` environment) Conda environment to be active

```
conda activate annotation_tutorial
./run_test.sh micro
```

##### TOGA config

TOGA requires a cluster specific config file to run nextflow. These files are in `TOGA/nextflow_config_files` and are built as default for `slurm`.

As such, if you are using a `slurm` cluster, you probably do not need to alter any of these files and you can move to the next installation step.

The first script is `call_cesar_config_template.nf `

```
// SLURM config file for CESAR jobs
// since CESAR have various memory requirements, this
// is just a template, TOGA will fill this itself
// depending on the CESAR job bucket
process.executor = 'slurm'
process.time = '24h'  // mostly 8h is enough, just for robustness
process.memory = "${_MEMORY_}G"  // to be replaced
process.cpus = 1  // CESAR utilizes a single core only
executor.queueSize = 1000  // nextflow default is 100 - too few
```

For SGE, we use the code below. I found that my SGE preferred the `process {}` and
`executor {}` syntax, but both should be acceptable.

```
process {
    executor = "sge"
    penv = "smp"
    memory = '10G'
    cpus = '1'
    time = '10h'
    clusterOptions = { "-V -l h_vmem=10G -V -P simpsonlab -l h_stack=32M -l h_rt=10:00:00" }
}

executor {
    name = "sge"
    queueSize = 1000
    queueStatInterval = "10s"
}
```

The next file that may need to be edited is `extract_chain_features_config.nf`, which looks like

```
// SLURM config for chain features extraction jobs
// relatively lightweighted jobs
process.executor = 'sge'
process.memory = '10G'
process.time = '1h'
process.cpus = 1
```

On SGE we use:

```
process {
    executor = "sge"
    penv = "smp"
    memory = '10G'
    cpus = '1'
    time = '10h'
    clusterOptions = { "-V -l h_vmem=10G -V -P simpsonlab -l h_stack=32M -l h_rt=10:00:00" }
}

executor {
    name = "sge"
    queueSize = 1000
    queueStatInterval = "10s"
}
```

The last file is `extract_chain_features_queue.nf`

```
// SLURM config for chain features extraction jobs
// relatively lightweighted jobs
process.executor = 'sge'
process.memory = '10G'
process.time = '1h'
process.cpus = 1
```

On SGE we use:

```
process {
    executor = "sge"
    penv = "smp"
    memory = '10G'
    cpus = '1'
    time = '10h'
    clusterOptions = { "-V -l h_vmem=10G -V -P simpsonlab -l h_stack=32M -l h_rt=10:00:00" }
}

executor {
    name = "sge"
    queueSize = 1000
    queueStatInterval = "10s"
}
```

##### Kent Utils

The `kent` binaries are universally valuable packages used for converting
genomic file types to one another. Accordingly, a number of these binaries are
required in this tutorial.

If these binaries don't work (e.g., genePredToGtf throws an error), it could be that you have a very new version of ubuntu If so, redownload the binaries after removing `.v369` from the link for the most updated version. (e.g., wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chainSwap)

```
# Navigate back to "external"
cd ..

mkdir -p kent

cd kent

# Binaries for generating chain files for TOGA

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/chainSwap

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/axtChain

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/faToTwoBit

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/pslPosTarget

# Binaries for converting TOGA outputs, and making reference annotations for TOGA

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/genePredToGtf
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bedToGenePred
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/genePredToBed
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/gtfToGenePred
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/gff3ToGenePred

# some clusters may require an extra step to ensure that each binary is executable
for i in `ls` ; do chmod +x $i ; done

# Navigate out of "kent"
cd ..
```

##### StringTie 

The `StringTie` binary is required for individuals with RNA-seq data. Many different installation sources can be found at their website (https://ccb.jhu.edu/software/stringtie/), however we found the easiet approach was to clone their github repository.

```
git clone https://github.com/gpertea/stringtie
cd stringtie
make release
cd ..
```

##### Interproscan

Finally, interproscan is used for annotating functional domains for proteins

```
mkdir my_interproscan
cd my_interproscan
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.69-101.0/interproscan-5.69-101.0-64-bit.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.69-101.0/interproscan-5.69-101.0-64-bit.tar.gz.md5

# Recommended checksum to confirm the download was successful:
md5sum -c interproscan-5.69-101.0-64-bit.tar.gz.md5
tar -xvzf interproscan-5.69-101.0-64-bit.tar.gz

cd interproscan-5.69-101.0-64-bit

python3 setup.py -f interproscan.properties

cd ../..

# Must return *interproscan-5.69-101.0-64-bit.tar.gz: OK*
# If not - try downloading the file again as it may be a corrupted copy
```

Lastly, we use a postprocessing perl script for infernal to make a GFF file from the ncRNAs. This script is downloaded directly to `/external`.

```
# Download perl script to convert output of infernal to tblout2gff
wget https://raw.githubusercontent.com/nawrockie/jiffy-infernal-hmmer-scripts/master/infernal-tblout2gff.pl
chmod +x infernal-tblout2gff.pl
```

### 2. Data packages

Some publicly available datasets (e.g., UniProt, Rfam) are required to complete the end to end tutorial. These datasets are downloaded into `GenomeAnnotationTutorial/data`.

```
cd ../data
```

#### uniprot_sprot

The uniprot_sprot database provides experimentally validated protein sequences
to cross-reference candidate transcripts against. This datbase is used in 
`Step 3: transcript selection`. We also remove duplicate fasta headings and generate a blastdb, which is required for the blast step in transcript selection. 

```
mkdir -p uniprot_sprot

cd uniprot_sprot

wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

gunzip uniprot_sprot.fasta.gz

seqkit rmdup -s < uniprot_sprot.fasta > uniprot_sprot_nodup.fasta

makeblastdb -in uniprot_sprot_nodup.fasta -dbtype prot -out uniprot_sprot_nodup.fasta -title "UniPro Sprot  database without duplicated sequences" -parse_seqids

cd ..
```

#### Rfam database

The Rfam database is required for finding ncRNA genes with known secondary 
structures and for identifying them into their RNA family. We also make a blastdb for this database for seeding.

```
mkdir -p Rfam

cd Rfam

wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/Rfam.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz
gunzip Rfam.fa.gz

# Install Clanin
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin

# Compress Rfam covariance models -- required for cmscan
cmpress Rfam.cm

# Add Rfam family file for annotations
https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/family.txt.gz
wget family.txt

seqkit rmdup -s < Rfam.fa > Rfam_nodup.fa

makeblastdb -in Rfam_nodup.fa -dbtype nucl -out Rfam_nodup -title "Rfam database without duplicated sequences" -parse_seqids

cd ..
```

#### Vertrbrate OrthoDB (Odb) protein stequences

BRAKER uses the Odb database to help predict protein locations in an assembly.
These sequences therefore must be downloaded in `/data`

```
mkdir -p braker_protein
cd braker_protein

wget https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/Vertebrata.fa.gz

gunzip Vertebrata.fa.gz

cd ..
```

#### Example data

We include a chromosome from our naked mole-rat assembly (plus RNA-seq and ISO-seq data) to try this tutorial.

These data can be pulled from the zenodo repository associated with this tutorial from ~/data

```
wget https://zenodo.org/records/14962941/files/example_data.tar.gz
tar -xvzf example_data.tar.gz
```

The last step in the setup is building a reference genome for a species. We provide a markdown to build this reference in "PreprocessReferenceSpecies.md". Please move to this markdown to continue building the workflow.




