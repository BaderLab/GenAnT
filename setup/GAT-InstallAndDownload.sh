#!/bin/bash

# clone repo. We assume this is done.
# git clone https://github.com/BaderLab/GenomeAnnotationTutorial.git
# module load singularity
# conda activate annotation_tutorial

# We assume that your cachedir is somewhere that can hold a lot of data (i.e., not in a root directory)

# mkdir -p singularitycache
# mkdir -p singularitytemp

# export SINGULARITY_CACHEDIR=/path-to/singularitycache
# export SINGULARITY_TMPDIR=/path-to/singularitytemp
#  

# Make data and external directory

cd $1

cd GenomeAnnotationTutorial

mkdir -p data
mkdir -p external

# build external directories
cd external

# Singularity images

mkdir -p singularity_images

cd singularity_images



echo "install singularities"

echo "braker3"
singularity build braker3.sif docker://teambraker/braker3:latest

echo "braker3 - long read"
singularity build braker3_lr.sif docker://teambraker/braker3:isoseq

echo "braker3 - cactus"

singularity build cactus.v2.9.3.sif docker://quay.io/comparative-genomics-toolkit/cactus:v2.9.3

echo "mikado_gat"

singularity build mikado_gat.sif docker://risserlin/mikado:ubuntu22_mikado2.3.2

cd ../../

echo "Removing cached and temporary space"
rm -r singularitycache
rm -r singularitytemp

cd external

echo "building TOGA binary"

git clone https://github.com/hillerlab/TOGA.git
cd TOGA
python3 -m pip install -r requirements.txt --user
./configure.sh

cd ../


echo "building kent utils directory"

mkdir -p kent

cd kent

# Binaries for generating chain files for TOGA

  # if these binaries don't work (e.g., genePredToGtf throws an error), it could be that you have a very new version of ubuntu
  # If so, redownload the binaries after removing `.v369` from the link for the most updated version.
  # e.g., wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chainSwap

  
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

cd ../

echo "building stringtie"

  git clone https://github.com/gpertea/stringtie
  cd stringtie
  make release

cd ../

echo "building interproscan"

mkdir my_interproscan
cd my_interproscan
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.69-101.0/interproscan-5.69-101.0-64-bit.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.69-101.0/interproscan-5.69-101.0-64-bit.tar.gz.md5

# Recommended checksum to confirm the download was successful:
md5sum -c interproscan-5.69-101.0-64-bit.tar.gz.md5
tar -xvzf interproscan-5.69-101.0-64-bit.tar.gz

cd interproscan-5.69-101.0

python3 setup.py -f interproscan.properties

cd ../

cd ../../external

# infernal table output to gff file postprocessing perl script.

# Download perl script to convert output of infernal to tblout2gff
wget https://raw.githubusercontent.com/nawrockie/jiffy-infernal-hmmer-scripts/master/infernal-tblout2gff.pl
chmod +x infernal-tblout2gff.pl


##
### Moving onto data
##

cd ../data


echo "Downloading uniprot_sprot database. We convert this into a blastdb separately"

mkdir -p uniprot_sprot

cd uniprot_sprot

wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

gunzip uniprot_sprot.fasta.gz

# seqkit rmdup -s < uniprot_sprot.fasta > uniprot_sprot_nodup.fasta

# makeblastdb -in uniprot_sprot_nodup.fasta -dbtype prot -out uniprot_sprot_nodup.fasta -title "UniPro Sprot  database without duplicated sequences" -parse_seqids 


echo "Downloading rfam database. We convert this into a blastdb separately"

cd ../

mkdir -p Rfam

cd Rfam

wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/Rfam.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.fa
gunzip Rfam.cm.gz

# Install Clanin
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin

# Compress Rfam covariance models -- required for cmscan
cmpress Rfam.cm

cd ../

# seqkit rmdup -s < Rfam.fa > Rfam_nodup.fa

# makeblastdb -in Rfam_nodup.fa -dbtype nucl -out Rfam_nodup -title "Rfam database without duplicated sequences" -parse_seqids

mkdir -p braker_protein

cd braker_protein

wget https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/Vertebrata.fa.gz

gunzip Vertebrata.fa.gz


