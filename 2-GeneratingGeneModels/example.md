# Example: Generating gene models for chromosome 28 of the Naked Mole Rat

## Installations

### LiftOff

Create a conda environment and check that liftoff is working

```
conda create -n liftoff_env liftoff=1.6.3
conda activate liftoff_env
liftoff -h
conda deactivate liftoff
```

Alternatively, there is a Docker available for this version of LiftOff, although we will assume you are using Conda. Here is how you can pull the Docker container:

```
docker pull staphb/liftoff
docker run -v "$(pwd)":/tmp staphb/liftoff liftoff -h
```

### TOGA

Navigate to the directory "external"

```
cd ../external
```

Create a Conda environment for TOGA with Python v3.11 (recommended by the authors) and NextFlow

```
conda create -n toga_env_conda python=3.11 nextflow
conda activate toga_env_conda
```

Create a virtual environment with python3 for installing of TOGA's requirements

```
python3 -m venv toga_env_pip
```

Install the CACTUS aligner and unzip

```
wget https://github.com/ComparativeGenomicsToolkit/cactus/releases/download/v2.8.4/cactus-bin-v2.8.4.tar.gz
tar -xzf cactus-bin-v2.8.4.tar.gz
```

Enter the CACTUS directory, make some modifications to virtual environment and reactivate

```
cd cactus-bin-v2.8.4
#
printf "export PATH=$(pwd)/bin:\$PATH\nexport PYTHONPATH=$(pwd)/lib:\$PYTHONPATH\nexport LD_LIBRARY_PATH=$(pwd)/lib:\$LD_LIBRARY_PATH\n" >> ../toga_env_pip/bin/activate
#
source ../toga_env_pip/bin/activate
```

Install various tools with pip

```
python3 -m pip install -U setuptools pip wheel
python3 -m pip install -U .
python3 -m pip install -U -r ./toil-requirement.txt
```

Exit directory and then test CACTUS

```
cd ..
cactus --help
```

Download the binary for TOGA and install requirements

```
git clone https://github.com/hillerlab/TOGA.git
cd TOGA
python3 -m pip install -r requirements.txt --user
./configure.sh
```

Test that TOGA is working and exit directory

```
./run_test.sh micro
cd ..
```

Install tools from the Comparative Genomics Toolkit (aka "Kent"). These are the binaries for generating chain files for TOGA:

```
mkdir -p kent
cd kent
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chainSwap
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/axtChain
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslPosTarget
```

Download the binaries for converting TOGA outputs and making reference annotations for TOGA:

```
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToGenePred
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/genePredToBed
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/gtfToGenePred
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gff3ToGenePred
```

Ensure each binary is executable

```
for i in `ls` ; do chmod +x $i ; done
```

Deactivate virtual environments until using TOGA and get back to `external`

```
deactivate # Deactivates Python3 environment
conda deactivate # Deactivates Conda environment
cd ..
```

### HISAT2 + StringTie

Create a Conda environment for HISAT2. Install SAMTools in this environment, as well

```
conda create -n hisat2_env hisat2=2.1.0
conda activate hisat2_env
hisat2 -h
deactivate
```

Alternatively, pull a Docker container for HISAT2 and test that it's working

```
docker pull nanozoo/hisat2
docker run -v "$(pwd)":/tmp nanozoo/hisat2 hisat2 -h
```

Clone StringTie's Github

```
 git clone https://github.com/gpertea/stringtie
```

Navigate into their repository and run `make release`, then go back to `external` and test StringTie

```
cd stringtie
make release
cd ..
stringtie --help
```

### BRAKER3

Pull BRAKER3 SIF file, list the files in `external` to see that it is there, and test that it's working

```
singularity build braker3.sif docker://teambraker/braker3:latest
ls
singularity exec braker3.sif braker.pl --help
```

Test BRAKER3 using their provided tests and open the log files to make sure they ran properly

```
singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test1.sh .
singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test2.sh .
singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test3.sh .
export BRAKER_SIF=./braker3.sif
bash test1.sh
bash test2.sh
bash test3.sh
```

After looking at the log files (and any other outputs), move them to a directory to clean things up

```
mkdir braker3_tests
mv test* braker3_tests
```

## Download a reference genome

We're going to use the mouse genome as a reference for the gene liftover methods, LiftOff and TOGA. Download the mouse genome and its annotation from RefSeq.

```
cd ../example_data
mkdir mouse_reference
cd mouse_reference
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz  
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz
gunzip *.gz
```

Return to this directory

```
cd ../../2-GeneratingGeneModels
```

## Run LiftOff

Make a directory to store the results and navigate to this directory

```
mkdir liftoff_example_results
cd liftoff_example_results
```

Activate Conda environment

```
conda activate liftoff_env
```

Run LiftOff; this takes about one hour with one thread

```
nohup liftoff \
 ../../example_data/NMRchr28.fa \
 ../../example_data/mouse_reference/GCF_000001635.27_GRCm39_genomic.fna \
 -g ../../example_data/mouse_reference/GCF_000001635.27_GRCm39_genomic.gff \
 -o nmr_chr28_liftoff.gff \
 -copies -flank 0.5 >& nohup.liftoff.out
```

You can count the number of genes and other features like this:

```
grep -c -P "\tgene\t" nmr_chr28_liftoff.gff
```

Get out of the LiftOff directory and deactivate the environment

```
cd ..
conda deactivate
```

## Run TOGA

First, look at the CACTUS configuration file provided in the `example_data` folder. This points towards the NMR chr28 which was soft-masked by Earl Grey, and the mouse reference genome downloaded from RefSeq. Note that these are relative file paths from this current working directory.

```
less ../example_data/two_species_cactus_config.txt
```

Align the two genome FASTA files with CACTUS, storing files in the directory `tmp1` (do not create this directory in advance or CACTUS will fail); this takes about 7 hours.

```
mkdir cactus_example_results
nohup cactus ./tmp1 \
 ../example_data/two_species_cactus_config.txt \
 ./cactus_example_results/target_ref.hal \
 --binariesMode local >& cactus_example_results/nohup.cactus.out
```

The output is `./cactus_example_results/target_ref.hal`. Convert this to a HAL file using tools from the Comparative Genomics Toolkit. First, set your path to the CACTUS and Kent bins to make these tools easier to access.

```
cactusbin=$(realpath ..)/external/cactus-bin-v2.8.4/bin
kentbin=$(realpath ..)/external/kent
```

Now make the 2bit file for the reference (mouse); takes less than a minute

```
$cactusbin/hal2fasta cactus_example_results/target_ref.hal mouse | $kentbin/faToTwoBit stdin cactus_example_results/reference.2bit
```

And the 2bit file for the target (naked mole rat chr28)

```
$cactusbin/hal2fasta cactus_example_results/target_ref.hal NMR_chr28 | $kentbin/faToTwoBit stdin cactus_example_results/target.2bit
```

Convert the HAL file to a BED file for both reference and target

```
$cactusbin/halStats --bedSequences mouse cactus_example_results/target_ref.hal > cactus_example_results/reference.bed
$cactusbin/halStats --bedSequences NMR_chr28 cactus_example_results/target_ref.hal > cactus_example_results/target.bed
```

Turn the HAL and BED files into a PSL file; takes less than a minute

```
$cactusbin/halLiftover --outPSL cactus_example_results/target_ref.hal NMR_chr28 \
 cactus_example_results/target.bed mouse /dev/stdout | \
 $kentbin/pslPosTarget stdin cactus_example_results/reference-to-target.psl
```

Finally, convert the PSL file to a chain file which can be used for TOGA; takes less than a minute

```
$kentbin/axtChain -psl -linearGap=loose \
 cactus_example_results/reference-to-target.psl \
 cactus_example_results/reference.2bit \
 cactus_example_results/target.2bit \
 cactus_example_results/reference-to-target.chain
```

There is one last step before running TOGA, which is to convert the reference (mouse) GFF file to a BED file with only protein-coding genes

```
gffread ../example_data/mouse_reference/GCF_000001635.27_GRCm39_genomic.gff \
 --keep-genes -C --no-pseudo --bed -o \
 ../example_data/mouse_reference/GCF_000001635.27_GRCm39_genomic.coding.bed
```

Only keep the first 12 columns to perfectly follow the BED12 format

```
cut -f1-12 ../example_data/mouse_reference/GCF_000001635.27_GRCm39_genomic.coding.bed > ../example_data/mouse_reference/GCF_000001635.27_GRCm39_genomic.coding.12.bed
```

Now use the BED file output by GFFRead to create isoforms.txt by isolating gene and transcript IDs

```
cut -f 13 GCF_000001635.27_GRCm39_genomic.coding.bed | sed 's/.*geneID=//;s/;gene_name.*//' > isoforms.txt
paste isoforms.txt <(cut -f 4 GCF_000001635.27_GRCm39_genomic.coding.bed) | sponge isoforms.txt
```

Create isoforms file for TOGA (these formatting steps specifically work with a RefSeq GFF file)

```
grep -P "\tmRNA\t" GCF_000001635.27_GRCm39_genomic.gff > ../example_data/mouse_reference/isoforms.txt
sed -i 's/.*Parent=gene-//;s/;gbkey.*//' ../example_data/mouse_reference/isoforms.txt
sed -i 's/;Dbxref[^;]*;Name=/\t/' ../example_data/mouse_reference/isoforms.txt
```

Once all of this file conversion has finished, it's time to run TOGA. First, create a variable pointing to the TOGA and CESAR installations and a new output directory

```
togabin=$(realpath ..)/external/TOGA
cesarbin=$(realpath ..)/external/TOGA/CESAR2.0/cesar
mkdir toga_example_results
```

Run TOGA

```
nohup $togabin/toga.py \
 cactus_example_results/reference-to-target.chain \
 ../example_data/mouse_reference/GCF_000001635.27_GRCm39_genomic.coding.12.bed \
 cactus_example_results/reference.2bit cactus_example_results/target.2bit \
 --project_dir toga_example_results \
 --project_name mouse_to_NMR_chr28 \
 --isoforms ../example_data/mouse_reference/isoforms.txt \
 --cesar_binary $cesarbin >& toga_example_results/nohup.toga.out
```

## Run HISAT2 + StringTie

We have provided BAM files with RNA-seq data from the naked mole-rat, tissues skin and thyroid. Since BAM files mean that the reads are already aligned to the genome, we will not need to run HISAT2 here. You can find a template of how to run HISAT2 in the README for this section. Instead we will just make sure that the BAM files are sorted, and then go right on to running StringTie.



