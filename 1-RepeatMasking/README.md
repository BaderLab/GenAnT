# 1. Repeat masking

The first step in genome annotation is to identify and mask repetitive regions. These can make up over half of a mammalian genome, and can cause trouble when generating genome annotations. For instance, repeat regions may interfere with sequence alignment, as they create an intractable number of alignment matches; repeats may also contain open reading frames (ORFs), and annotation software may mistake these ORFs as genes (therefore increasing the false positive rate when gene models are generated). Therefore it is important that a genome sequence is masked so that annotation software doesn't attempt to place gene models in these regions. Genomes can be soft-masked (repeat regions turned from uppercase letters to lowercase letters in the FASTA file) or hard-masked (repeat regions converted into strings of capital Ns). Soft-masking is generally recommended, and more leniant as it allows for gene-models initiated in non-repeat regions to extend into repeat regions.

Before repeat masking, it's best to check if your genome came with repeats already masked. You can do this by peaking into your genome file using `less genome.fa`, and check if you can see any Ns or lowercase letters.

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

	earlGrey \
-g $outDir/assembly/assembly.fa \
 -s $species \
-o . \
-t 50 \
-r rodentia \
-d yes


```

#### Earl Grey: installing/running/troubleshooting

- We have had success running Earl Grey on a desktop and high performance compute cluster
- Earl Grey can easily be installed using conda (e.g. `conda create -n earlgrey -c conda-forge -c bioconda earlgrey=5.1.0`); an error causing Earl Grey to crash mid-run required an update to Numpy (`pip install numpy --upgrade`)
- Earl Grey takes multiple days to run (be prepared for up to a week)
- Earl Grey does not like spaces in any directory names
