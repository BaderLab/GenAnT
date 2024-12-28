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
source toga_env_pip/bin/activate
```

Download the binary for TOGA and install requirements

```
cd ../external
git clone https://github.com/hillerlab/TOGA.git
cd TOGA
python3 -m pip install -r requirements.txt --user
./configure.sh
```

Test that TOGA is working

```
./run_test.sh micro
```

Install the CACTUS aligner back in `external` and unzip

```
cd ..
wget https://github.com/ComparativeGenomicsToolkit/cactus/releases/download/v2.8.4/cactus-bin-v2.8.4.tar.gz
tar -xzf cactus-bin-v2.8.4.tar.gz
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
```

### HISAT2 + StringTie

Pull a Docker container for HISAT2 and test that it's working

```
docker pull nanozoo/hisat2
docker run -v "$(pwd)":/tmp nanozoo/hisat2 hisat2 -h
```

Alternatively, create a Conda environment

```
conda create -n hisat2_env hisat2=2.1.0
conda activate hisat2_env
hisat2 -h
deactivate
```

Clone StringTie's Github

```
 git clone https://github.com/gpertea/stringtie
```

Navigate into their repository and run `make release`, then go back to `external`

```
cd stringtie
make release
cd ..
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

We're going to use the mouse genome as a reference for the gene liftover methods, LiftOff and TOGA. Download the mouse genome and its annotation from RefSeq

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

Run LiftOff (started at 6:23pm Dec 28)

```
nohup liftoff \
 ../../example_data/NMRchr28.fa \
 ../../example_data/mouse_reference/GCF_000001635.27_GRCm39_genomic.fna \
 -g ../../example_data/mouse_reference/GCF_000001635.27_GRCm39_genomic.gff \
 -o nmr_chr28_liftoff.gff \
 -copies -flank 0.5 >& nohup.liftoff.out
```





