# Example data

Once everything is installed as outlined in the README document, you should be able to run the tutorial on some example data. The example can be run following method 2 (submitting a script) or method 3 (running the Snakemake pipeline). However, before running the example, you will need to download the example data from Zenodo. We include a chromosome from our naked mole-rat assembly (plus RNA-seq and ISO-seq data) to try this tutorial. We assume that the following `wget` command will be performed in /path-to-GAT/GenAnT (see /setup/GAT-InstallAndDownload.md for details). 

```
wget https://zenodo.org/records/14962941/files/example_data.tar.gz
tar -xvzf example_data.tar.gz
```

You will also need the reference genome for the naked mole-rat. In `/path-to-GenomeAnnotationTutorial/GenAnT/data/references`:

```
mkdir -p HetGlaV2_female ; cd HetGlaV2_female

wget https://ftp.ensembl.org/pub/release-115/fasta/heterocephalus_glaber_female/dna/Heterocephalus_glaber_female.Naked_mole-rat_maternal.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-115/gff3/heterocephalus_glaber_female/Heterocephalus_glaber_female.Naked_mole-rat_maternal.115.gff3.gz
for i in *.gz ; do gunzip $i ; echo $i ; done
  bash ../../../setup/reference_directory_ensembl.sh \
  . \ # path to the reference genome directory
  ~/GenAnT \ #  path to GenomeAnnotationTutorial ( `GenAnT` included)
  Heterocephalus_glaber_female.Naked_mole-rat_maternal.dna_sm.toplevel.fa \
  Heterocephalus_glaber_female.Naked_mole-rat_maternal.115.gff3
```

## Method 2 (script submission)

Assuming everything is properly set up, running the tutorial without flow control involves submitting the execute script with two positional arguments:
`bash Execute_GAT_in_serial.sh path-to-GAT path-to-Conda`
For us this looks like

```
bash Execute_GAT_in_serial.sh \
/.mounts/labs/simpsonlab/users/dsokolowski/projects/GenAnT \
/.mounts/labs/simpsonlab/users/dsokolowski/miniconda3
```
Even annotating one chromosome is relatively resource intensive, so we reccomend submitting this as a job (Recommended: 64G mem, 16 cores, 72h runtime)

Lastly, this script has a `module load singularity`. If you access singularity differently (e.g., apptainer, conda environment etc.) then replace that line with what you need to have singularity accessible.

## Method 3 (Snakemake)

Running the example data using the snakemake pipeline takes slightly more work:

There is a config file (`config_example.yaml`) located in the `GenAnT_Snakemake` directory which will need to be modified to run the example data. To run the tutorial with the example data, you need to change: "/path-to-conda/miniconda3/" to the path to the miniconda directory where the `annotation_tutorial` environment lives (e.g., /.mounts/labs/simpsonlab/users/dsokolowski/miniconda3/) and change "path-to-GenAnT/" to the path where you cloned "GenAnT" (e.g., /.mounts/labs/simpsonlab/users/dsokolowski/projects/GenAnT). This will have to be done for each line in `config_example.yaml` that uses the GenAnT file path. Let's walk through this line by line.

Here is what the example configuration file looked like on my computer after I changed all the variables:

```
sourceDir: "/opt/miniconda3/bin/activate"
externalDir: "/home/baderlab/zclarke/GenAnT/external
dataDir: "/home/baderlab/zclarke/GenAnT/data"

# Parameters describing your assembly + annotation
outDir: "/scratch8/badertmp/GenAnT_Example_stranded"
target: "example"
species: "heterocephalus_glaber"

assemblyFile: "/home/baderlab/zclarke/GenAnT/data/example_data/NMRchr28.fa"
MaskedAssemblyFile: "none"
MaskedAssemblyAnnotation: "none"
rnaseqDir: "/home/baderlab/zclarke/GenAnT/data/example_data/RNAseq_alignment"
isoseqDir: "/home/baderlab/zclarke/GenAnT/data/example_data/ISOseq_alignment"
customRef: "FALSE"
liftoffRef: "FALSE"

# Parameters describing your reference assemblies + annotation

refToga: "mouse"
TogaDir: "/home/baderlab/zclarke/GenAnT/data/references/mmus_GRC39"
refTogaFa: "/home/baderlab/zclarke/GenAnT/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.fna"
refTogaBed: "/home/baderlab/zclarke/GenAnT/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.toga.bed"
refTogaIsoform: "/home/baderlab/zclarke/GenAnT/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.isoforms.tsv"

refToga2: "NMR"
TogaDir2: "/home/baderlab/zclarke/GenAnT/data/references/HetGlaV2_female"
refTogaFa2: "/home/baderlab/zclarke/GenAnT/data/references/HetGlaV2_female/Heterocephalus_glaber_female.Naked_mole-rat_maternal.115.clean.fa"
refTogaBed2: "/home/baderlab/zclarke/GenAnT/data/references/HetGlaV2_female/Heterocephalus_glaber_female.Naked_mole-rat_maternal.115.clean.bed"
refTogaIsoform2: "/home/baderlab/zclarke/GenAnT/data/references/HetGlaV2_female/Heterocephalus_glaber_female.Naked_mole-rat_maternal.115.clean.isoforms.tsv"

refLiftOff: "NMR"
refLiftOffFa: "/home/baderlab/zclarke/GenAnT/data/references/HetGlaV2_female/Heterocephalus_glaber_female.Naked_mole-rat_maternal.115.clean.fa"
refLiftOffGff: "/home/baderlab/zclarke/GenAnT/data/references/HetGlaV2_female/Heterocephalus_glaber_female.Naked_mole-rat_maternal.115.gffread.gff"

orthofinderFA: "/home/baderlab/zclarke/GenAnT/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.nostop.protein.faa"
orthofinderTab: "/home/baderlab/zclarke/GenAnT/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.table.txt"

# Tool specific parameters
dfamDB: "rodentia"
brakerOdbFaa: "Vertebrata.fa"
mikadoScore: "mammalian.yaml"
mirmachineClade: "Mammalia"
stranded: "TRUE"
rnaseqCov: 5
isoseqCov: 1
stringtieMerge: "FALSE"

# Adding GFF files that you have already computed
customGFF: "none"
customLiftoff: "none"
customTOGA: "none"
customTOGA2: "none"
customBraker: "none"
customStringtie: "none"
```

Otherwise, follow the steps in "3. A Snakemake pipeline".

