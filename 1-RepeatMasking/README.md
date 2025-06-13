# 1. Repeat masking

The first step in genome annotation is to identify and mask repetitive regions. These can make up over half of a mammalian genome, and can cause trouble when generating genome annotations. For instance, repeat regions may interfere with sequence alignment, as they create an intractable number of alignment matches; repeats may also contain open reading frames (ORFs), and annotation software may mistake these ORFs as genes (therefore increasing the false positive rate when gene models are generated). Therefore it is important that a genome sequence is masked so that annotation software doesn't attempt to place gene models in these regions. Genomes can be soft-masked (repeat regions turned from uppercase letters to lowercase letters in the FASTA file) or hard-masked (repeat regions converted into strings of capital Ns). Soft-masking is generally recommended, and more leniant as it allows for gene-models initiated in non-repeat regions to extend into repeat regions.

Before repeat masking, it's best to check if your genome came with repeats already masked. You can do this by peaking into your genome file using `less genome.fa`, and check if you can see any Ns or lowercase letters.

Note. the scripts in scripts/ expect that variables listed are in your path (e.g., with export). The snakemake pipeline has the equivalent scripts but with positional arguments.

#### Earl Grey

Repeat masking can be done with [Earl Grey](https://github.com/TobyBaril/EarlGrey). Earl Grey integrates multiple common repeat masking tools such as RepeatMasker, which maps repetitive elements from a database, and RepeatModeler, which identifies repeats de novo.  It also uses multiple tools such as cd-hit-est, LTR_finder, rcMergeRepeats, and custom scripts to identify, annotate, filter, and aggregate repeat regions genome wide. Earl Grey is a command-line tool that can be run with a single line of code in a Unix environment, and produces figure-quality summaries of a genome’s transposable element landscape in conjunction with repeats annotated in general feature format (GFF) which are required for downstream analysis. Earl Grey relies on databases of repeat elements, such as DFam, that are used to identify repeats in your genome. The user can specify what clade of species they are working with, which indicates which repeat database Earl Grey should use.

Input to Earl Grey is the FASTA file from your species, the name of your species, and the output directory. It is also helpful to specify the search term used for RepeatMasker with `-r`, which indicates which set of repeats to look for (e.g. “eukarya”). `-d` is a flag that indicates whether or not you would like soft-masking. If you put "yes", Earl Grey will output a soft-masked genome that you can use directly for subsequent annotation steps. The output is stored in multiple folders, with the most important information located in `summaryFiles`. 

```
earlGrey \
 -g your_genome.fasta \
 -s your_species_name \
 -o ./output_directory \
 -r repeat_clade \
 -d yes \
 -t number_of_threads
```
In our pipeline using the example data (a small chromosome from the naked molerat, this command looks as follows:

```
outDir=/path-to-output-directory/
species="heterocephalus_glaber"
scripts/run_earl_grey.sh
```

The main part of the script looks as follows

```
outDir=/path-to-output-directory/
species="heterocephalus_glaber"

earlGrey \
-g $outDir/assembly/assembly.fa \
 -s $species \
-o . \
-t 50 \
-r rodentia \
-d yes
```

Earl Grey is a wrapper of many repeat annotation tools, resulting in many directories. Below is the output directory of Earl Grey when run on the example:

![earlgrey_output](https://github.com/user-attachments/assets/614231d2-5b03-485a-a73b-b9b3183849fb)

As expected, everything needed for downstream genome annotation is in the “heterocephalus_glaber_summaryFiles” directory, shown below:

![earlgrey_summaryFiles](https://github.com/user-attachments/assets/ef2e997d-5ce8-4035-bbd0-6eaf3194cc9b)

The files needed for the rest of the pipeline are “heterocephalus_glaber.filteredRepeats.bed”, which contains the curated set of repeats, and “heterocephalus_glaber.softmasked.fasta”, which is the fasta file used in Braker and TOGA annotations. “heterocephalus_glaber.filteredRepeats.bed” is used in non-coding RNA seeding, as scRNA, srpRNA, and “SINE/tRNA-RTE”  are called as part of repeat elements in earl grey. The other files are useful for summarizing and analyzing repeats.

EarlGrey 5.1.0 (the version at the time of GenAnT’s release) recomputes TE divergence after finalizing TE libraries. This step is very time-consuming, and we have experienced instances on a slow node where this step alone can run for 200 hours. You do not need to wait for these divergences to be recalculated if you don’t want to. If you navigate to the directory below, you can extract the “heterocephalus_glaber.filteredRepeats.bed” and then mask your assembly with `bedtools maskfasta`. 

```
ls earl_grey/heterocephalus_glaber_EarlGrey/heterocephalus_glaber_mergedRepeats/looseMerge/
```

Below is the `species_mergedRepeats/looseMerge` directory in EarlGrey where repeat annotations are stored:

![earlgrey_annotationDirectory](https://github.com/user-attachments/assets/09285585-7d7e-436f-9013-f47c5794ebb2)

#### Earl Grey: installing/running/troubleshooting

- We have had success running Earl Grey on a desktop and high performance compute cluster
- Earl Grey can easily be installed using conda (e.g. `conda create -n earlgrey -c conda-forge -c bioconda earlgrey=5.1.0`); an error causing Earl Grey to crash mid-run required an update to Numpy (`pip install numpy --upgrade`)
- We have moved to running the earlgrey singularity container to simplify installation and running. `scripts/run_earl_grey.sh` executes earl grey through the singularity image.
- Earl Grey takes multiple days to run (be prepared for up to a week)
- If Earl Grey runs out of time, you can rerun it with the same command and directories and it picks up (more-or-less) where it left off.
- Earl Grey does not like spaces in any directory names
- In our testing `forksys:  Program terminated by a signal 9.` usually means an out of RAM issue
