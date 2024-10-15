# How to annotate a mammalian genome

This is a tutorial on how to annotate a newly sequenced mammalian genome. This takes the user from their FASTA sequence to a high-quality GFF file annotated with gene symbols. We break the process into four main steps:

1. Repeat masking
2. Generating gene models
3. Combining and filtering gene models
4. Predicting gene function by annotating gene symbols

We recommend tools and best practices for each step, providing code to help the user execute each task. It is up to the user to install the tools that we recommend for this pipeline, however we make note of challenging installation processes that we have encountered with certain tools and how we overcame these challenges.

Generally, genome annotation does not have a comparable ground truth, so we use different sources of evidence to annotate the most likely gene models. These gene models are considered hypotheses for where the genes are located on the genome, but false positives and false negatives will always exist. This pipeline uses existing tools and quality-checking software to try to minimise both of these error rates to create a high-quality annotation.

### Notes on computational requirements

We expect the user to be familiar with installing and running command-line tools, as genome annotation relies on such tools. Additional familiarity with R may be helpful for some more advanced tasks. Many tools can be run on a desktop, but some are very computationally intensive and require the use of a high-performance compute cluster (e.g. Compute Canada). The speed of many tools will improve if they have access to multiple threads and can therefore run tasks in parallel. To check how many threads you can specify when running tools, check the documentation of your compute cluster or run `nproc` on your local desktop.

Because genome annotation relies on a number of different command-line tools, we recommend creating or using unique environments for each tool on your machine, if possible. One tool may rely on one version of a piece of software, while another tool may rely on a more recent or older version; if tools share the same environment, such conflicts may prevent each tool from running properly. Virtual environments can be used or created with tools like Docker, Conda, or Python.

Docker containers exist for certain tools, and typically mean that a tool is packaged with all of its requirements and is ready for you to use; you can check for them by running `docker search name_of_tool`. If a Docker container is listed, you can typically use it by running `docker run -v "$(pwd)":/tmp name_of_container command`. `docker run` means that you are running the container, `-v` indicates where you are mounting the volume, which essentially gives Docker a temporary space on your machine to store data (here just given the placeholder `/tmp`), `name_of_container` is replaced by the name of whatever container you want to try, and `command` indicates the command of the tool you want to use. 

Conda environments can be created using the command `conda create -n name_of_environment required_package_1 required_package_2 ...` where additional packages separated by a space can replace the elipses. Conda considers all the packages required by the user, and creates an isolated environment where all package versions should be compatable. If the user only needs one particular package with all of its dependencies, Conda will automatically find and install all dependencies into the environment when only the one needed is specified. Environments can be activated by running `conda activate name_of_environment`.

## Repeat masking

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
#### Earl Grey: installing/running/troubleshooting

- We have had success running Earl Grey on a desktop and high performance compute cluster
- Earl Grey can easily be installed using conda (e.g. `conda create -n earlgrey -c conda-forge -c bioconda earlgrey=4.4.0`); an error causing Earl Grey to crash mid-run required an update to Numpy (`pip install numpy --upgrade`)
- Earl Grey takes multiple days to run (be prepared for up to a week)
- Earl Grey does not like spaces in any directory names

## Generating gene models

Gene models are hypotheses about the locations of genes and their corresponding features (e.g. mRNA, exons, introns) on the genome. These hypotheses are supported by a variety of evidence, including RNA-alignment information, the presence of ORFs, protein sequence conservation, gene structure and order along the sequence. As gene models are hypotheses to the locations and structures of real genes, annotations may contain both false positives (e.g. a random ORF-like sequence) and false negatives (e.g. a real gene that was missed in the annotation process). To reduce these errors, it is important to use high quality evidence for the existence of genes, and good annotation tools that perform well.

We thoroughly describe two complementary approaches to generate high-quality gene models: homology-based annotation and RNA-and-protein-alignment-based annotation (Figure 1). Broadly, homology-based annotations assume that thousands of gene models will be shared between a reference species (e.g. mouse) and a target species (e.g. woodchuck), at the level of DNA sequence similarity and gene structure (e.g., number of exons). In contrast, transcript assembly through functional annotation assumes that the location of uniquely mapping paired end RNA-seq data represents an expressed region of the genome and a candidate for a gene model.

It is important to note that resulting annotations vary depending on (1) the quality of the genome sequence being annotated, (2) the quality of the evidence provided to inform the annotation (e.g. RNA-seq, homology), and (3) the quality of the bioinformatic tool applied. Therefore, it is best to perform many different annotations, test each one for quality (described in box 1), and choose the best results that can be used in step three, which involves integrating gene models into a single, complete set. Although tools may change over time, homology- and alignment-based annotations should both be generated, tested, and combined.

### Homology-based annotation

Homology-based genome annotation is the derivation of gene models in your species from homologous gene models found in other species. Most gene structures and sequences are conserved across related species, making homologous alignments from a reference species with high-quality gene structures an accurate and computationally efficient method to annotate your species. Many tools are capable of performing homology-based annotation, but we have had the most luck with [LiftOff](https://github.com/agshumate/Liftoff) and [the Tool to infer Orthologs from Genome Alignments (TOGA)](https://github.com/hillerlab/TOGA).

#### Finding annotations for liftover

Before performing homology-based annotation, one needs to decide which genome to pull information from. The best homology-based annotations according to various quality metrics (e.g. BUSCO, GffCompare) occur when the reference genome is high-quality, with a high-quality annotation from a closely related species. To find such a genome, we recommend searching RefSeq or ENSEMBL for a few of the most closely related annotated species to yours, and trying liftover on those species. 

High-quality genomes tend to have smaller contig or scaffold numbers (i.e. the genome sequence is divided into larger chunks), ideally close to the number of chromosomes found in the species, smaller L50s (the smallest number of contigs that make up half of the genome sequence), and larger N50s (the smallest contig length that half of the genome sequence is contained in, of the largest contigs). High-quality annotations can be assumed if the genome is annotated by RefSeq or ENSEMBL. RefSeq annotations are also evaluated for quality using BUSCO, where a curated set of single-copy orthologs is compared to the gene models identified in the annotation; a high quality annotation is indicated by a single-copy ortholog detection rate close to 100%, and missing or fragment orthologs close to 0%. We recommend the user searching these databases for a few of the most closely-related species, comparing these genome statistics, and selecting the assembly and annotation (or multiple) with the most favorable statistics.

#### LiftOff

LiftOff is a gene liftover tool that aligns gene sequences from the reference genome to the target genome using a single line of Unix code, making it quick and easy to use. It uses minimap2 to align the genes to the genome with high accuracy and with relatively low computational resources, so the tool can be run on a desktop computer. LiftOff is a command-line tool that takes a FASTA file and GFF/GTF file from a reference species, and the FASTA file from the target species, and creates a GFF/GTF output file for your species based on the reference annotations. It also provides the user with a list of unmapped genes ("unmapped_features.txt"), which may indicate alignment challenges. The "-copies" flag indicates that LiftOff will look for additional gene copies in the new genome. Because LiftOff is so quick and easy to use, the user can easily use LiftOff to generate annotations from multiple reference species and compare the resulting annotation quality.

```
liftoff \
 -g annotation_of_related_species.gff \
 your_genome.fasta \
 genome_of_related_species.fasta \
 -o output_annotation.gff \
 -u unmapped_features.txt \
 -copies \
 -p number_of_threads
```

#### LiftOff: installing/running/troubleshooting

- Some Docker containers exist and work well for LiftOff (e.g. `docker run staphb/liftoff liftoff`)
- LiftOff also works well with a Conda environment (e.g. `conda create -n liftoff -c bioconda liftoff python=3`)

#### The Tool to infer Orthologs from Genome Alignments (TOGA)

TOGA accurately annotated genes across vertebrates with higher rates of divergence. TOGA relies on a chain file connecting the reference and target species, which is a file that indicates which sections of the reference genome align to which sections of the target genome. This allows TOGA to use synteny to inform its annotation liftover, which improves accuracy considering groups of genes are often conserved across species. Chain files are developed by post-processing whole genome alignments between two species, typically using executable binary scripts developed to improve the compatibility between genomic data types and [the UCSC genome browser](https://github.com/ucscGenomeBrowser/kent). We recommend using the [CACTUS alignment tool](https://github.com/ComparativeGenomicsToolkit/cactus) when generating the initial alignments between two distantly related species. Preparing the data for TOGA and running TOGA is a multistep process:

1. Align genomes with Cactus

   Cactus takes a text configuration file as input, which is a two-species phylogenetic tree of the reference and target species. A template of such a file is as follows, replacing “target” and “ref” with the species names and files:
   ```
   (target:1.0,ref:1.0);
   target       /path-to-target/target.soft.fa
   mouse      /path-to-reference/reference.soft.fa
   ```
   Since the FASTA files of each species are listed in the config file, these do not need to be specified as addition input to CACTUS. Note that each FASTA file is expected to be soft-masked. Cactus outputs a file ending in `.hal`, which stores information about the alignment. CACTUS requires you to specify a temporary directory where Cactus stores large quantities of files while it's running. This temporary directory will change depending on what system you are using to run Cactus. On a local desktop, a temporary directory may simply by `/tmp`, whereas a high performance compute cluster may have a designated temporary directory to use, such as `$SCRATCH/tmp`. Cactus can then be run as follows:

   ```
   cactus $SCRATCH/tmp \
   two_species_cactus_config.txt \
    target_ref.hal \
    --binariesMode local
   ```

2. Convert HAL file to chain file

   This is a multistep process also described by the ComparativeGenomicsToolkit [here](https://github.com/ComparativeGenomicsToolkit/hal/blob/chaining-doc/doc/chaining-mapping.md). The first step involves converting both the reference and target FASTA files to a compressed 2bit format. This can be done using additional tools that are accessible in the [Cactus Github repository](https://github.com/ComparativeGenomicsToolkit/cactus) in this directory: `/path-to-cactus/external/cactus-bin-v2.2.3/bin`. We can set this as a variable to make the tools easier to access.

   `cactusbin=/path-to-cactus/external/cactus-bin-v2.2.3/bin`

   Each FASTA file can be converted to 2bit with the following two commands below. Each `hal2fasta` command requires the HAL file output by CACTUS as input as well as the reference or target FASTA file. The output is directly piped into `faToTwoBit` ("stdin" indicates that `faToTwoBit` takes the piped input) which outputs a compressed [2bit file](https://genome.ucsc.edu/FAQ/FAQformat.html#format7).

   ```
   $cactusbin/hal2fasta target_ref.hal name_of_reference | faToTwoBit stdin reference.2bit
   $cactusbin/hal2fasta target_ref.hal name_of_target | faToTwoBit stdin target.2bit
   ```

   Convert the alignments stored in the HAL file to a BED file using the `halStats` command.

   ```
   $cactusbin/halStats --bedSequences name_of_reference target_ref.hal > reference.bed
   $cactusbin/halStats --bedSequences name_of_target target_ref.hal > target.bed
   ```

   Next, create pairwise alignments, which are stored in a resulting PSL file. This can be done using `halLiftover`. The `--outPSL` flag indicates that the output will be a PSL file; the command takes the HAL file, the name of the target species, the target BED file (created in the previous step), and the name of the reference species as input. The output is specified as `/dev/stdout` which means the output will be printed to the screen. This output is piped into the `pslPosTarget` command which forces the alignments to the positive strand, and outputs the results into a PSL file.

   ```
   $cactusbin/halLiftover --outPSL target_ref.hal name_of_target \
      target.bed name_of_reference /dev/stdout | \
      $cactusbin/pslPosTarget stdin reference-to-target.psl
   ```
   
   Finally, the PSL file can be converted to a chain file using the `axtchain` command from [ucscGenomeBrowser](https://github.com/ucscGenomeBrowser/kent) (sometimes referred to as KentUtils). This bins the alignments at various depths, generalizing the alignment so that instead of storing alignments at specific base pairs, they are stored as blocks of homologous regions. `axtChain` takes the flag `-psl` to indicate PSL input, the recommended parameter setting `-linearGaps=loose`, the PSL file, and both 2bit files as input. The output is the chain alignment file that can now be used for TOGA.

   ```
   axtChain -psl -linearGap=loose reference-to-target.psl reference.2bit target.2bit reference-to-target.chain
   ```
   
3. Perform homology-based annotation with TOGA

   Now that the input files have been prepared and processed, TOGA can be run with one line of UNIX code. The inputs to TOGA are the chain file created in the previous step, the 2bit files for both the reference and the target also created in the previous step, and transcript annotations from the reference species in [BED12](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format. GTF files can be converted to BED12 files using tools available from [ucscGenomeBrowser](https://github.com/ucscGenomeBrowser/kent). This is a two step process: First, convert the GTF file to a genePred file by performing `gtfToGenePred annotation.gtf annotation.genePred`. Then, convert the genePred file to a BED12 file with `genePredToBed annotation.genePred annotation.bed`.

   Isoform data from the reference species is highly recommended when running TOGA. These data are provided in a two-column TSV file with a header. The left column is the gene ID and the right column is the transcript ID; a single gene can be associated with multiple transcripts. This can be created directly from the BED or GFF annotation file of the reference species (we have provided a script XXX that can create this). 
   
   ```
   toga.py \
   reference-to-target.chain \
   reference_annotation.bed \
   reference.2bit target.2bit \
   --project_name ref_to_target \
   --isoforms isoforms.tsv
   ```

   The output of TOGA...

#### TOGA and associated tools: installing/running/troubleshooting

- TOGA takes a lot of work to install and run, but its results are worth it! Benchmarks have shown that TOGA performs very well compared to other methods and its preservation of gene symbols can be quite useful.
- TOGA may not work if run on a desktop; most genomes are large enough that they require a high performance compute cluster
- Cactus and TOGA can easily be installed with Conda

   
   

