# 2. Generating gene models

Gene models are hypotheses about the locations of genes and their corresponding features (e.g. mRNA, exons, introns) on the genome. These hypotheses are supported by a variety of evidence, including RNA-alignment information, the presence of ORFs, protein sequence conservation, gene structure and order along the sequence. As gene models are hypotheses to the locations and structures of real genes, annotations may contain both false positives (e.g. a random ORF-like sequence) and false negatives (e.g. a real gene that was missed in the annotation process). To reduce these errors, it is important to use high quality evidence for the existence of genes, and good annotation tools that perform well.

We thoroughly describe two complementary approaches to generate high-quality gene models: homology-based annotation and RNA-and-protein-alignment-based annotation (Figure 1). Broadly, homology-based annotations assume that thousands of gene models will be shared between a reference species (e.g. mouse) and a target species (e.g. woodchuck), at the level of DNA sequence similarity and gene structure (e.g., number of exons). In contrast, transcript assembly through functional annotation assumes that the location of uniquely mapping paired end RNA-seq data represents an expressed region of the genome and a candidate for a gene model.

It is important to note that resulting annotations vary depending on (1) the quality of the genome sequence being annotated, (2) the quality of the evidence provided to inform the annotation (e.g. RNA-seq, homology), and (3) the quality of the bioinformatic tool applied. Therefore, it is best to perform many different annotations, test each one for quality (described in box 1), and choose the best results that can be used in step three, which involves integrating gene models into a single, complete set. Although tools may change over time, homology- and alignment-based annotations should both be generated, tested, and combined.

Note. the scripts in scripts/ expect that variables listed are in your path (e.g., with export). The snakemake pipeline has the equivalent scripts but with positional arguments.

### Homology-based annotation

Homology-based genome annotation is the derivation of gene models in your species from homologous gene models found in other species. Most gene structures and sequences are conserved across related species, making homologous alignments from a reference species with high-quality gene structures an accurate and computationally efficient method to annotate your species. Many tools are capable of performing homology-based annotation, but we have had the most luck with [LiftOff](https://github.com/agshumate/Liftoff) and [the Tool to infer Orthologs from Genome Alignments (TOGA)](https://github.com/hillerlab/TOGA).

#### Finding annotations for liftover

Before performing homology-based annotation, one needs to decide which genome to pull information from. The best homology-based annotations according to various quality metrics (e.g. BUSCO, GffCompare) occur when the reference genome is high-quality, with a high-quality annotation from a closely related species. To find such a genome, we recommend searching RefSeq or ENSEMBL for a few of the most closely related annotated species to yours, and trying liftover on those species. 

High-quality genomes tend to have smaller contig or scaffold numbers (i.e. the genome sequence is divided into larger chunks), ideally close to the number of chromosomes found in the species, smaller L50s (the smallest number of contigs that make up half of the genome sequence), and larger N50s (the smallest contig length that half of the genome sequence is contained in, of the largest contigs). High-quality annotations can be assumed if the genome is annotated by RefSeq or ENSEMBL. RefSeq annotations are also evaluated for quality using BUSCO, where a curated set of single-copy orthologs is compared to the gene models identified in the annotation; a high quality annotation is indicated by a single-copy ortholog detection rate close to 100%, and missing or fragment orthologs close to 0%. We recommend the user searching these databases for a few of the most closely-related species, comparing these genome statistics, and selecting the assembly and annotation (or multiple) with the most favorable statistics.

Our scripts and instructions to download and prepare a reference genome are found in /setup.

for a RefSeq genome, the instructions are here:
```
reference_directory_refseq.md
```
and can be executed with it's paired script `reference_directory_refseq.sh`

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

Using the example data in our tutorial (`run_liftoff.sh`), the liftoff command looks like this:

Building a reference directory for this tutorial is executed with `/setup/reference_directory_refseq.sh`

```
tutorialDir=/path-to-GAT/GenomeAnnotationTutorial

refLiftOffGff=$tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.gffread.gff

refLiftOffGff=$tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.fna

outDir=/path-to-output-dir/

gffutils-cli create $refLiftOffGff -o referencegff_db # builds reference db in output directory so you can use the same reference on multiple assemblies simultaneously.

liftoff -db referencegff_db $outDir/assembly/assembly.softmasked.fa $refLiftOffFa -o ./liftoff.gff -u ./unmapped.liftoff.txt -copies

```

#### LiftOff: installing/running/troubleshooting

- Some Docker containers exist and work well for LiftOff (e.g., `docker run staphb/liftoff liftoff`)
- LiftOff also works well with a Conda environment (e.g. `conda create -n liftoff -c bioconda liftoff python=3`)

#### The Tool to infer Orthologs from Genome Alignments (TOGA)

TOGA accurately annotated genes across vertebrates with higher rates of divergence. TOGA relies on a chain file connecting the reference and target species, which is a file that indicates which sections of the reference genome align to which sections of the target genome. This allows TOGA to use synteny to inform its annotation liftover, which improves accuracy considering groups of genes are often conserved across species. Chain files are developed by post-processing whole genome alignments between two species, typically using executable binary scripts developed to improve the compatibility between genomic data types and [the UCSC genome browser](https://github.com/ucscGenomeBrowser/kent). We recommend using the [CACTUS alignment tool](https://github.com/ComparativeGenomicsToolkit/cactus) when generating the initial alignments between two distantly related species. Preparing the data for TOGA and running TOGA is a multistep process:

1. Align genomes with Cactus

   Cactus takes a species tree file as input, which is a two-species phylogenetic tree of the reference and target species. A template of such a file is as follows, replacing “target” and “ref” with the species names and files:
   ```
   (target:1.0,ref:1.0);
   target       /path-to-target/target.soft.fa
   ref      /path-to-reference/reference.soft.fa
   ```

We use the `scripts/make_cactus_tree.sh` script, which inputs file names and paths and creates the species tree for you. It assumes your softmasked.fasta is copied to `$outDir/assembly/assembly.softmasked.fa`

```
tutorialDir=/path-to-GAT/GenomeAnnotationTutorial
outDir=/path-to-output-dir/
target="NMR"
refToga="mouse"
refTogaFa=$tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.fna

scripts/make_cactus_tree.sh
```

   Since the FASTA files of each species are listed in the config file, these do not need to be specified as additional input to CACTUS. Note that each FASTA file is expected to be soft-masked. CACTUS outputs a file ending in `.hal`, which stores information about the alignment. CACTUS requires you to specify a temporary directory where CACTUS stores large quantities of files while it's running. This temporary directory will change depending on what system you are using to run CACTUS. On a local desktop, a temporary directory may simply be `/tmp`, whereas a high-performance compute cluster may have a designated temporary directory to use, such as `$SCRATCH/tmp`. CACTUS can then be run as follows:

   ```
    cactus \
    $SCRATCH/tmp \
    cactus_in.txt \
    target_ref.hal \
    --workDir ./ \
    --binariesMode local \
    --maxCores 32 \
    --maxMemory 64G \
    --realTimeLogging \
    --batchSystem single_machine # can remove if you have cactus working on your system for other projects, but is only required to switch for alignments of >2 species.
   ```

1. Convert HAL file to chain file

   This is a multistep process also described by the ComparativeGenomicsToolkit [here](https://github.com/ComparativeGenomicsToolkit/hal/blob/chaining-doc/doc/chaining-mapping.md). The first step involves converting both the reference and target FASTA files to a compressed 2bit format. This can be done using additional tools that are accessible in the [CACTUS Github repository](https://github.com/ComparativeGenomicsToolkit/cactus) in this directory: `/path-to-cactus/external/cactus-bin-v2.2.3/bin`. We can set this as a variable to make the tools easier to access.

   `cactusbin=/path-to-cactus/external/cactus-bin-v2.2.3/bin`

   We will also be using additional tools from the Comparitive Genomics Toolkit. In our installation guide, we have these installed in a directory called `kent`. It will also be helpful to set a variable pointing to this location.

   `kentbin=/path-to-kent/`

   Each FASTA file can be converted to 2bit with the following two commands below. Each `hal2fasta` command requires the HAL file output by CACTUS as input as well as the reference or target FASTA file. The output is directly piped into `faToTwoBit` ("stdin" indicates that `faToTwoBit` takes the piped input) which outputs a compressed [2bit file](https://genome.ucsc.edu/FAQ/FAQformat.html#format7).

   ```
   $cactusbin/hal2fasta target_ref.hal name_of_reference | $kentbin/faToTwoBit stdin reference.2bit
   $cactusbin/hal2fasta target_ref.hal name_of_target | $kentbin/faToTwoBit stdin target.2bit
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
   $kentbin/axtChain -psl -linearGap=loose reference-to-target.psl reference.2bit target.2bit reference-to-target.chain
   ```

  Performing the CACTUS alignment and generating the chain file (2. and 3.) can be performed with the `scripts/cactus_align_and_chain_sif.sh`. We switched to a singularity image as we found it was easier to use across systems. This script is executed as:
  ```
  outDir=/path-to-output-dir/
  tutorialDir=/path-to-GAT/GenomeAnnotationTutorial
  externalDir=$tutorialDir/external
  TogaDir=$tutorialDir/GenomeAnnotationTutorial/data/references/mmus_GRC39 # needed to bind directory to singularity
  target="NMR"
  refToga="mouse"
  refTogaFa=/path-to-GenomeAnnotationTutorial/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.fna


  scripts/cactus_align_and_chain_sif.sh

  ```
4. Perform homology-based annotation with TOGA

   Now that the input files have been prepared and processed, TOGA can be run with one line of UNIX code. The inputs to TOGA are the chain file created in the previous step, the 2bit files for both the reference and the target also created in the previous step, and transcript annotations from the reference species in [BED12](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format. GTF files can be converted to BED12 files using tools available from [ucscGenomeBrowser](https://github.com/ucscGenomeBrowser/kent). This is a two step process: First, convert the GTF file to a genePred file by performing `gtfToGenePred annotation.gtf annotation.genePred`. Then, convert the genePred file to a BED12 file with `genePredToBed annotation.genePred annotation.bed`.

   Isoform data from the reference species is highly recommended when running TOGA. These data are provided in a two-column TSV file with a header. The left column is the gene ID and the right column is the transcript ID; a single gene can be associated with multiple transcripts. This can be created directly from the BED or GFF annotation file of the reference species. The toga BED and isoforms files are generated using the "reference_directory*.sh" scripts in `setup` that best match your reference genome.
   
   ```
   toga.py \
   reference-to-target.chain \
   reference_annotation.bed \
   reference.2bit target.2bit \
   --project_name ref_to_target \
   --isoforms isoforms.tsv
   ```

   The output of TOGA contains many informative files that describe the performance of the annotation, with the main output being query_annotations.bed. This is a BED12 formatted file that contains the annotations with their predicted orthologs. In order to continue with this file, it will have to be converted to GTF and then gff format. This can be done using tools from the [comparative genomics toolkit](https://github.com/ComparativeGenomicsToolkit). First, convert the BED file to a Gene Prediction or "GenePred" file, a table file format often used with the UCSC Genome Browser. This serves as a temporary format that can then be converted to a GTF file using another tool in the toolkit. Lastly, we use gffread to convert this output into a consistently formatted gff3 file for `3-CombineAndFilter` The code is as follows:

   ```
   bedToGenePred toga_output_directory/query_annotation.bed toga_output_directory/query_annotation.genePred

   genePredToGtf file toga_output_directory/query_annotation.genePred toga_output_directory/query_annotation.gtf
   ```

   This gives you the output GTF from TOGA in the TOGA output directory. To make sure that the GTF file is formatted nicely in a way that will definitely work with other tools, you may wish to use GFFRead to quickly clean up the file, like so:

   ```
   gffread query_annotation.gtf --keep-genes -o toga_query_annotation.gffread.gff
   ```

   This gives you a GFF file that we have called toga_query_annotation.gffread.gff to help you keep track that this is your GFF file from TOGA.

   Using scripts within the annotation tutorial, TOGA is run with

 ```
  outDir=/path-to-output-dir/
  tutorialDir=/path-to-GAT/GenomeAnnotationTutorial
  externalDir=$tutorialDir/external
  scriptsDir=$tutorialDir/scripts
  target="NMR"
  refToga="mouse"
  refTogaBed=$tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.toga.bed 
  refTogaIsoform=$tutorialDir/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.isoforms.toga.tsv 

  scripts/run_toga.sh

   ```

We expect that the most typical instance for modern genome annotation is having a combination of RNA-seq and ISO-seq data, where there is not perfect overlap in tissues. A script wrapping combinations of RNA-seq and ISO-seq data is run with: `run_stringtie_flexible.sh`

```
scriptsDir=$tutorialDir/scripts
outDir=/path-to-output-dir/
scripts/run_stringtie_flexible.sh
```


#### TOGA and associated tools: installing/running/troubleshooting

- TOGA may not work if run on a desktop; most genomes are large enough that they require a high performance compute cluster
- Cactus and TOGA can easily be installed with Conda

### Transcript assembly using RNA-sequencing data

RNA- and/or protein-sequence alignment data can be used to inform gene models. Alignment-based methods work by aligning RNA or protein sequences to the genome to determine the location of transcribed and/or protein-coding genes. The specific tools used to perform alignment-based annotation depend on the sequencing data available to the user. If the user has access to high-quality RNA-seq data (e.g. 2 x 100bp paired-end sequencing ideally from as many tissues as possible), [HISAT2](https://daehwankimlab.github.io/hisat2/) and [StringTie2](https://github.com/skovaka/stringtie2) can be used to create an annotation or "transcript assembly" directly from these data. Otherwise, [BRAKER3](https://github.com/Gaius-Augustus/BRAKER) can be used with shorter, lower-quality, or no RNA-seq reads and a database of protein sequences to create an annotation.

#### Finding publicly available RNA-seq data

If you can't generate your own RNA-seq data, there may be publicly available data for your species you can use. One place where you can search for RNA-seq data is the [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra). In the SRA search bar, type `"species name" AND "rna seq"` to find RNA-seq for your species (e.g. `"mus musculus" AND "rna seq"` if you were annotating the mouse genome). If you see an RNA-seq dataset that fits your criteria and that you would like to download, find the experiment accession number that is listed two lines below the link to the data (often begins with "SRX"). SRA Toolkit (downloadable at https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) can use that accession number to download the raw FASTQ files for that dataset.

```
prefetch accession_number
# Creates folder with .sra file inside; folder and SRA file have the same name
fasterq-dump --split-files file_name/file_name.sra
```

The output are the FASTQ files that were submitted for that experiment. For paired-end RNA-seq data, this often means that the output will be two FASTQ files, representing each end of the paired-end reads. Such paired files should be processed together downstream, for instance, when being aligned to the genome. These FASTQ files can then be used for HISAT2 + StringTie2 or BRAKER3.

#### HISAT2 + StringTie2

To generate gene models directly from RNA-seq alignment, the RNA-seq reads first need to be aligned to your genome sequence. This can be done with the genome aligner HISAT2. In order to align the reads, an index needs to be generated for the genome you are annotating with HISAT2. The input is your genome FASTA sequence and the number of threads you would like to use to parallelize the process. HISAT2 outputs multiple files that either end in `.ht2` or `.ht2l` depending on the size of your genome.

```
hisat2-build \
 -p number_of_threads \
 your_genome.fasta \
 base_name_of_genome_index
```

Now you can run HISAT2 to align the RNA-seq reads to the genome you are annotating. The input required is the genome index created by `hisat2-build` and the input RNA-seq files in FASTQ format (or FASTA format if specified with `-f`). The `--dta` flag reports alignments tailored for transcript assemblers (as in this case). The `-1` and `-2` indicate the mates for paired-end RNA-seq. HISAT2 outputs a SAM alignment file.

```
hisat2-align \
 -p number_of_threads \
 --dta \
 -x base_name_of_genome_index \
 -1 first_mate.fastq \
 -2 second_mate.fastq \
 -S name_of_sam_alignment.sam
```

Now that the reads are aligned, StringTie2 can be used to generate gene models (this language is used for simplicity - StringTie2 finds "transcripts" rather than "genes" since it is based on aligning RNA-seq to the genome, but the concept is similar to gene model generation). However, StringTie2 takes BAM files as input, which are a binary version of SAM files. Therefore, you must first convert your SAM file(s) to BAM file(s) with [SAMtools](http://www.htslib.org/). The only input is the SAM file created with HISAT2. `-S` specifies SAM input, `-h` includes the header in the SAM output, `-u` indicates uncompressed BAM output (saves time when converting SAM to BAM).

```
samtools view \
 -@ number_of_threads \
 -Shu \
 -o name_of_sam_alignment.bam \
 name_of_bam_alignment.sam
```

StringTie2 also requires the BAM files to be sorted by reference position. This can also be done with SAMtools. 

```
samtools sort \
 -@ number_of_threads \
 -o name_of_sorted_bam_alignment.bam \
 name_of_bam_alignment.bam
```

Now the sorted BAM file can be used to predict gene models with StringTie2. The output is a GTF file.

```
stringtie \
 name_of_sorted_bam_alignment.bam \
 -o stringtie_output.gtf \
 -p number_of_threads
```


ISO-seq annotations are also performed with stringtie, and parameters are not drastically changed.

We use `minimap2` and `samtools` to perform and filter ISO-seq alignments. We assume your ISO-seq files are in fastq format (generated from something like https://ccs.how/, https://isoseq.how/getting-started.html).


```
outDir=/path-to-output-dir/
i=isoseq.tissue.fastq.gz

b=`basename $i .fastq.gz`

minimap2 -ax splice:hq -uf $outDir/assembly/assembly.fa $i > $b".sam" 

samtools view -bSq 1 $b".sam" > $b".bam"

samtools sort -o $b".sorted.bam" $b".bam"
```

Running `--stringtie` with ISO-seq data instead of RNA-seq requires in `-L` parameter

```
stringtie \
 -L \
 name_of_sorted_bam_alignment.bam \
 -o stringtie_output.gtf \
 -p number_of_threads
```

If you have RNA-seq and ISO-seq for the same tissue, transcripts can be generated by adding the `--mix` parameter.

```
stringtie \
   --mix tissue.RNAseq.bam tissue.ISOseq.bam \
   -l output \
   -o stringtie_out/"tissue.mix.gtf" \
   -p 8 --conservative

```

Modern genome annotations usually have a combination of RNA-seq and ISO-seq data, and the tissues will not perfectly overlap betweend datatypes (e.g, you use public RNA-seq or can have more tissues due to cost constraints).

Annotations are greatly improved with RNA-seq from multiple tissues, meaning that you'll probably have multiple StringTie outputs. These outputs are combined with stringtie --merge.

```
stringtie --merge -o $outDir/stringtie_out/stringtie.merged.gtf $outDir/stringtie_out/*gtf
```

Our script assumes that the RNAseq data lives in $outDir/RNAseq_alignment and the ISOseq data lives in $outDir/ISOseq_alignment, and that experiments of the same tissue have the same bam file name. For example your `kidney` experiments would live in:
* RNA-seq: `$outDir/RNAseq_alignment/kidney.bam`
* ISO-seq: `$outDir/ISOseq_alignment/kidney.bam`
`run_stringtie_flexible.sh` also deals with merging and directory structure.

```
outDir=/path-to-output-dir/
tutorialDir=/path-to-GAT/GenomeAnnotationTutorial
externalDir=$tutorialDir/external
scriptsDir=$tutorialDir/scripts
scripts/run_stringtie_flexible.sh
```



The only features in the output GTF file are transcripts and exons, with no prediction of coding sequences (typically indicated by "CDS" in the third column of the GTF file). Because of this, the output cannot be easily converted into a protein sequence and tested with BUSCO. Although this is not ideal, testing the quality and completeness of a genome annotation with BUSCO is not necessary if it will be combined with additional annotation sets and filtered using [Mikado](https://mikado.readthedocs.io/en/stable/) (explained later). If you used poor-quality or very short RNA-seq data, however (not recommended), there is a risk of generating short, fragmented, monoexonic transcripts. You can check to see if your annotation has many short, monoexonic transcripts using a summary statistics calculator provided by Mikado, `mikado util stats annotation.gff output_summary.tsv`, where `annotation.gff` is replaced by whatever annotation you want the summary statistics for, and `output_summary.tsv` is whatever you name the output summary statistics file. You can compare your summary statistics to those of another mammalian genome annotated by RefSeq or Ensembl. If the statistics are similar, this indicates an annotation that is likely of higher quality. However, if you notice that the average number of exons per transcript is very low and the number of monoexonic transcripts is very high in the genome you are annotating, this indicates that many of the gene models may be short or fragmented, and should potentially be excluded from the final annotation set or run through a pass of very stringent filtering with Mikado (tool explained later, noisy RNA-seq e.g. https://mikado.readthedocs.io/en/stable/Tutorial/Adapting/#case-study-2-noisy-rna-seq-data).



#### BRAKER3

If you don't have access to RNA-seq data or your RNA-seq reads are short (single-end short reads or paired-end reads shorter than 2x100bp), [BRAKER3](https://github.com/Gaius-Augustus/BRAKER) can be used to generate an annotation. BRAKER3 integrates RNA-seq alignment information with protein data and *ab initio* gene prediction. *Ab initio* gene predictors are mathematical models that are fed existing gene models to train their algorithms (i.e., the algorithms learn which aspects of genome structure are associated with different gene model features), so that they can then discover new gene models in genome sequences. The RNA sequences come from the species being annotated, whereas the protein sequences are typically from an online database of homologous sequences, like [OrthoDB](https://www.orthodb.org/). Internally, BRAKER3 uses HISAT2 to align the short RNA-seq reads to the genome, StringTie2 to create candidate gene models from these alignments, and [ProtHint](https://github.com/gatech-genemark/ProtHint) to predict CDS regions using these protein alignments. These data are then used as “hints” i.e., estimations of CDS region and intron placements) when generating *ab initio* gene models with [GeneMark-ETP](https://github.com/gatech-genemark/GeneMark-ETP) and [Augustus](https://github.com/Gaius-Augustus/Augustus). Finally, BRAKER3 can also identify tRNAs, snoRNAs, and UTRs.

The input to BRAKER3 is the soft-masked genome you wish to annotate (`your_genome.fasta`), the RNA sequences you wish to align, and the protein sequence database (`orthodb.fa`). If you generated the previous annotation using HISAT2 + StringTie, you would have already aligned RNA-seq to the genome, which would otherwise be done internally by BRAKER3. Therefore, the sorted BAM file(s) that you used as input to StringTie2 can also be used as input for BRAKER3 (e.g. `rna1.bam,rna2.bam` with more comma-separated files listed if you have additional BAM files). `name_of_your_species` is whatever name you want to call you species to distinguish the output; `number_of_cores` is the same as the number of threads used for other tools (they are slightly different ways to describe essentially the same thing); and `--gff3` indicates that the desired output is a GFF3 file.

In the context of our tutorial, the `orthodb.fa` used is the `Vertebrata.fa` file from orthoDb. (see https://github.com/BaderLab/GenomeAnnotationTutorial/blob/main/setup/GAT-InstallAndDownload.md for details)

```
singularity exec braker3.sif braker.pl \
 --threads=number_of_threads \
 --species=name_of_your_species \
 --genome=your_genome.fasta \
 --bam=rna1.bam,rna2.bam \
 --prot_seq=~/GenomeAnnotationTutorial/data/braker_protein/Vertebrata.fa \ # use different odb fasta if neccessary.
 --gff3
```

If you do not have RNA-seq data and wish to run BRAKER3 in "protein mode", simply remove the `--bam` flag specifying the RNA-seq data.

Running Braker3 with ISO-seq data uses the same syntax as with traditional RNA-seq data, but it uses a different singularity image. Namely, you use `braker3_lr.sif` from  `docker://teambraker/braker3:isoseq` instead of `braker3.sif` from `docker://teambraker/braker3:latest`

As such, braker3 with ISO-seq is executed with:

```
singularity exec braker3_lr.sif braker.pl \
 --threads=number_of_threads \
 --species=name_of_your_species \
 --genome=your_genome.fasta \
 --bam=iso1.bam,iso2.bam \
 --prot_seq=~/GenomeAnnotationTutorial/data/braker_protein/Vertebrata.fa \ # use different odb fasta if neccessary.
 --gff3
```
In step 3, braker.gff and braker_lr.gff are separate gff files input into `mikado`. In our testing, merging braker and braker_lr with mikado provides a nearly identical GFF file to if you integrate the same files with TSEBRA.

Running braker through our tutorial is performed with the `scripts/run_braker*.sh` scripts.

#### BRAKER3: installing/running/troubleshooting

- In our experience, BRAKER is most easily installed and implemented using the Singularity container that the BRAKER authors maintain: `singularity build braker3.sif docker://teambraker/braker3:latest`
- If installing Braker through other methods (e.g. a conda environment) then the `singularity exec braker3.sif` in the command is unnecessary
- You likely need to specify where the `augustus config directory` is. Installing the tutorial with `setup` would have this directory in `~/GenomeAnnotationTutorial/external/Augustus/config`. 
- We have found that the GFF file output by BRAKER3 has some formatting issues that can be fixed by running GFFRead, e.g. `gffread braker.gtf --keep-genes -o braker.gffread.gff`
