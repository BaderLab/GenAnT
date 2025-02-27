# 5. Assign gene symbols to protein-coding genes and some lncRNA genes

At this point in the tutorial, we now have a GFF file that contains the gene-models of protein-coding and non-coding RNA genes, and predicted gene symbols and products for the latter. This leaves the final step in the tutorial, which is assigning gene symbols and predicting the function of protein-coding genes and some lncRNA genes. We can use sequence similarity between our target species and closely-related species to identify predicted orthologs of the gene models from the target species. Assigning gene symbols is challenging because most genes in a mammalian genome originated from another gene (e.g. tandem duplication, gene fusion, translocation), meaning that many genes have at least one paralogous gene with high sequence similarity in exons. Three tools can be used to predict gene symbols in the target species: LiftOff, TOGA, and [OrthoFinder](https://github.com/davidemms/OrthoFinder).

### LiftOff and TOGA

You may have already used LiftOff and TOGA earlier in the tutorial, making it fairly straightforward to use their results. Both LiftOff and TOGA annotate the target species’ gene structures and assign reference gene symbols to the target with a high rate of agreement with the gene symbols found in Ensembl annotations. Therefore, these tools can be used to predict gene symbols for the final, integrated annotation. We transfer gene models by overlapping transcripts derived from TOGA and LiftOff to the final gene models, before transfering the gene symbol to the Mikado-filtered gene identifier (ID).

First, we will want to isolate the protein-coding genes and lncRNA from `full_annotation.gff`, which we created in step 4. We can do this with a simple `grep` statement, with the `-P` flag meaning that we are activating Perl to recognize the tab symbols on either side of "mRNA" and "lncRNA" so that we only isolate features where mRNA or lncRNA make up column 3. The `|` acts as an "or" statement, separating the patterns we wish to extract. This can be output to a file called `mikado.mRNA.lncRNA.gff`.

```
grep -P "\tmRNA\t|\tlncRNA\t" full_annotation.gff > mikado.mRNA.lncRNA.gff
```

Similarly, we will want to pull out the transcript features found by both LiftOff and TOGA. Although we referred to them by different names earlier, for simplicity we will refer to the outputs of LiftOff and TOGA from step 2 as `liftoff.gff` and `toga.gff`.

First, extract mRNA and lncRNA from `liftoff.gff` and `toga.gff` using `grep`. The mRNA features may have either "mRNA" or "transcript" in column 3 depending on if a GFF3 or GTF file was used as input to LiftOff.

```
grep -P "\tmRNA\t|\tlnc_RNA\t|\ttranscript\t" liftoff.gff > liftoff.mRNA.lncRNA.gff
grep -P "\tmRNA\t|\tlnc_RNA\t|\ttranscript\t" toga.gff > toga.mRNA.lncRNA.gff
```

Now that we have isolated the transcripts from Mikado, LiftOff, and TOGA, we can use BEDTools to determine which gene models overlap with each other. We will determine (1) which LiftOff gene models overlap with which Mikado gene models to assign LiftOff gene symbols to the Mikado gene models, and (2) which TOGA gene models overlap with which Mikado gene models to assign gene symbols from TOGA. We cab deternube these overlaps using `bedtools intersect` with `-a` pointing at the Mikado GFF file, and `-b` pointing at the LiftOff and TOGA gene models. `-wo` specifies that the output will contain both A and B entries with the number of base pairs overlap between the two features, with only A features with overlap being reported.

```
bedtools intersect -a mikado.mRNA.lncRNA.gff -b liftoff.mRNA.lncRNA.gff -wo > mikado.liftoff.mRNA.lncRNA.txt
bedtools intersect -a mikado.mRNA.lncRNA.gff -b toga.mRNA.lncRNA.gff -wo > mikado.toga.mRNA.lncRNA.txt
```

We don't need all of the information spit out by BEDTools, so we can make the following adjustments. First, let's extract the first nine columns and then the 10th to 18th columns from the BEDTools outputs to make two separate GFF files.

LiftOff:

```
cut -f1-9 mikado.liftoff.mRNA.lncRNA.txt > liftoff_overlap.mikadoInfo.gff
cut -f10-18 mikado.liftoff.mRNA.lncRNA.txt > liftoff_overlap.liftoffInfo.gff
```

TOGA:

```
cut -f1-9 mikado.toga.mRNA.lncRNA.txt > toga_overlap.mikadoInfo.gff
cut -f10-18 mikado.toga.mRNA.lncRNA.txt > toga_overlap.togaInfo.gff
```

These can now be read into an R notebook called `MakeGeneSymbolTableLiftOffTOGA.rmd`, which processes these files to output a table of which LiftOff, TOGA, and non-coding RNA gene symbols from step 4 align to the different gene IDs. The file is a tab-delimited (TSV) file called `gene_symbols.tsv`.

### OrthoFinder

OrthoFinder, a tool that maps sequence-similarity relationships between proteins across two or more species based on their sequences, can also be used to identify predicted orthologs. OrthoFinder builds gene trees, considers gene duplication events, is considered to be one of the most accurate orthologs inference methods, and was used for gene naming in the DNA zoo annotation project. OrthoFinder outputs lists of protein-protein sequence-similarity relationships that can be used to infer gene-gene relationships.

OrthoFinder is especially helpful for gene symbol labeling if you wish to label the annotated target genome with gene symbols from a species that you DIDN'T use for homology-based annotation (since you won't have the outputs of LiftOff or TOGA). To run OrthoFinder, you need a FASTA file of protein sequences from both your reference and target species. You can even use multiple reference species to analyse more orthologous relationships, but we'll just assume two species.

Get protein sequences from target species with GFFRead:

```
gffread -y target_proteins.faa -g target_genome.softmasked.fasta -S full_annotation.gff
```

Let's say you wish to compare the protein sequences of your target species to that of human. Download human FASTA and GFF and unzip.

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
#
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz
#
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz
#
gunzip GCF_000001405.40_GRCh38.p14_genomic.gff.gz
```

Use these files to get protein sequences with GFFRead:

```
gffread -y reference_proteins.faa -g GCF_000001405.40_GRCh38.p14_genomic.fna -S GCF_000001405.40_GRCh38.p14_genomic.gff
```

Put protein sequences in their own directory called `protein_seqs`

```
mkdir protein_seqs
mv target_proteins.faa protein_seqs
mv reference_proteins.faa protein_seqs
```

Run OrthoFinder on protein sequences. `-t` and `-a` both specify the number of threads for different processes (sequence search and analysis respectively), `-o` specifies the name of the output directory, `-f` points to the folder of protein sequences.

```
orthofinder -t number_of_threads -a number_of_threads -o orthofinder -f protein_seqs
```

The OrthoFinder results we are interested in are pairwise protein-protein relationships between the reference and target species. They can be found in `./orthofinder/Results_date_of_run/Orthologues/Orthologues_target_proteins/target_proteins__v__reference_proteins.tsv`. These can be read into an R notebook called `AddOrthoFinderGenes.Rmd`.



