# Example: Generating gene models for chromosome 28 of the Naked Mole Rat

## Installations

### LiftOff

Pull a Docker container for LiftOff and test that the container is working

```
docker pull staphb/liftoff
docker run -v "$(pwd)":/tmp staphb/liftoff liftoff -h
```

Alternatively, create a conda environment and check that liftoff is working (this version of LiftOff is specified to match the Docker container)

```
conda create -n liftoff_env liftoff=1.6.3
conda activate liftoff_env
liftoff -h
conda deactivate liftoff
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

## Download a reference genome
