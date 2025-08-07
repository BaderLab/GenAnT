# How to evaluate annotation quality

Genome annotation quality across species is constantly improving, and no genome annotation is perfect. For each annotation generated, it is useful to perform a quality assessment to determine how well certain annotation tools worked for your data.

#### BUSCO

The most common way to assess the completeness and quality of the annotation, is to use [BUSCO](https://busco.ezlab.org/) to compare the gene models found in your genome to a curated set of single-copy orthologs for all domains of life stored in the OrthoDB database. To do this, a GFF file and FASTA file can be translated into a FASTA file of protein sequences (e.g. by using GFFRead). The inputs to GFFRead are the genome FASTA file and annotation GFF file of the annotation that you wish to translate into protein sequences, and `protein_sequences.faa` is whatever you decide to name the output.

```
gffread -y protein_sequences.faa -g genome.fasta annotation.gff
```

BUSCO can be run on this FASTA file in protein mode, which functionally scans the protein sequence file for thousands of conserved genes. BUSCO requires the user to pick a pre-existing lineage sequence database to be used for its comparison (e.g. “Glires”); one should specify an available lineage closest to your species of interest. You can see what lineages are available by running:

```
busco --list-datasets
```

BUSCO returns statistics indicating if the expected protein sequences are found, fragmented, or missing. BUSCO scores can be compared across annotations as a judge of quality, with higher BUSCO scores indicating higher quality genomes. BUSCO can be run as follows, requiring the protein sequences created from the GFF file as input, as well as the lineage (e.g. `glires` or the more general `mammalia`), `-m prot` indicating that we are looking at protein sequences, and `-c` indicating the number of threads to speed up the process.

```
busco -i protein_sequences.faa \
 -l lineage \
 -m prot \
 -c number_of_threads
```

BUSCO outputs a directory of results, the most important of which are the percentage of single-copy orthologs that were captured; this statistic is also output to the screen. It's important to note, that the maximum BUSCO score (i.e. the maximum percentage of single-copy orthologs) an annotation can have is equal to the percentage that have been captured in the genome sequence that you are annotating. For example, if your genome sequence only has a BUSCO score of 98% (meaning that the sequence itself is missing BUSCO genes), then the BUSCO score based on the annotation has to be less than 98%. Here is how you can find the genome sequence BUSCO score:

```
busco -i genome.fasta \
 -l lineage \
 -m genome \
 -c nuumber_of_threads
```

#### Mikado util stats

It is also helpful to analyze feature characteristics of a particular annotation (e.g. average exons per transcript, number of monoexonic transcripts, gene lengths) as outliers of these statistics may indicate that there are inaccuracies. For instance, if an RNA-sequencing-alignment-based annotation has a large number of monoexonic transcripts compared to a homology-based annotation, this suggests that the former annotation may be fragmented into artificially small transcripts. Mikado comes with [a command that performs such a summary](https://mikado.readthedocs.io/en/stable/Usage/Utilities/#stats), and outputs a TSV file of summary statistics. The input is a GFF file (no FASTA file required).

```
mikado util stats annotation.gff annotation_summary.tsv`
```

If you see that an annotation has a bunch of short, monoexonic transcripts, there are a couple things you can do. First, you can stringently filter the annotation (e.g. using Mikado or another filtering tool) to remove gene models that are short and/or only contain a single exon. We describe how to use Mikado later in this tutorial, but additional guidance for this specific scenario can also be found [here](https://mikado.readthedocs.io/en/stable/Tutorial/Adapting/#case-study-2-noisy-rna-seq-data). However, such an annotation may have been generated as a result of using very short RNA-seq reads to assemble transcripts. In this case, the quality of the gene models would likely be improved if longer RNA-seq reads were available to rerun the annotation. If you decide to go ahead with these short, fragmented gene models, it may create certain challenges. For example, predicting orthologs would be very challenging; paralogs may be falsely predicted, or a predicted gene symbol may not be identified at all when the fragment really is part of a larger gene.

#### GffCompare

If working with different annotation sets where one is higher quality than another, it may be useful to compare these annotations directly. Different GFF files mapping to the same genome assembly can be compared with [GffCompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml). Briefly, GffCompare inputs a “query” GFF to a “reference” GFF, and outputs a parseable text file (“.stats”) describing how well the base pairs, exons, introns, and transcripts match each other. It can also be valuable if the researcher has a set of experimentally validated or manually curated gene models for their species, as they can compare these validated gene models to the gene models generated in their genome-wide annotations.

```
gffcompare -r reference_annotation.gff query_annotation.gff -o result_prefix
```

#### Synteny

Genomes contain collinear regions called syntenic blocks that are conserved across large evolutionary time spans. In the context of a reference and target species for genome annotation (e.g. mouse and woodchuck), these syntenic blocks typically contain a large number of genes in both species, and the orientation of these genes are often the same. This phenomena is used in some gene annotation pipelines (e.g. TOGA), however it can also be used to manually identify missing annotations or misassembly by comparing genome browser snapshots between the reference and target species. Similarly, while manual in nature, aligning functional genomics data to your assembly and viewing these data with your genome annotation on a genome browser is a crucial sanity check in determining whether genome annotations were correctly performed.

#### Functional genomics alignment

Aligning functional genomics data to your assembly and viewing these data with your genome annotation on a genome browser is a crucial sanity check in determining whether genome annotations were correctly performed. Such data can be the RNA-seq reads that have been aligned to the genome with HISAT2 earlier in the pipeline. The genome FASTA file, BAM alignment files, and the genome annotation GFF file can all be loaded into a genome browser at the same time, and junctions revealed by the BAM files will often align well with exons indicated by the annotation. If clear junctions exist but appropriate gene models are missing, this may indicate a mistake in the annotation. Common genome browsers include...
