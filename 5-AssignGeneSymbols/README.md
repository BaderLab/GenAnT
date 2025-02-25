# 5. Assign gene symbols to protein-coding genes

At this point in the tutorial, we now have a GFF file that contains the gene-models of protein-coding and non-coding RNA genes, and predicted gene symbols and products for the latter. This leaves the final step in the tutorial, which is assigning gene symbols and predicting the function of protein-coding genes. We can use sequence similarity between our target species and closely-related species to identify predicted orthologs of the gene models from the target species. Assigning gene symbols is challenging because most genes in a mammalian genome originated from another gene (e.g. tandem duplication, gene fusion, translocation), meaning that many genes have at least one paralogous gene with high sequence similarity in exons. Three tools can be used to predict gene symbols in the target species: LiftOff, TOGA, and [OrthoFinder](https://github.com/davidemms/OrthoFinder).

### LiftOff and TOGA

You may have already used LiftOff and TOGA earlier in the tutorial, making it fairly straightforward to use their results. Both LiftOff and TOGA annotate the target species’ gene structures and assign reference gene symbols to the target with a high rate of agreement with the gene symbols found in Ensembl annotations. Therefore, these tools can be used to predict gene symbols for the final, integrated annotation. We transfer gene models by overlapping transcripts derived from TOGA and LiftOff to the final gene models, before transfering the gene symbol to the Mikado-filtered gene identifier (ID).

First, we will want to isolate the protein-coding genes from `full_annotation.gff`, which we created in step 4. At this point, it only matters that we isolate the mRNA transcript features, which we can do with a simple `grep` statement, with the `-P` flag meaning that we are activating Perl to recognize the tab symbols on either side of "mRNA" so that we only isolate features where mRNA makes up column 2. This can be output to a file called `mikado.mRNA.gff`.

```
grep -P "\tmRNA\t" full_annotation.gff > mikado.mRNA.gff
```

Similarly, we will want to pull out the transcript features found by both LiftOff and TOGA. Although we referred to them by different names earlier, for simplicity we will refer to the outputs of LiftOff and TOGA from step 2 as `liftoff.gff` and `toga.gff`.

First, extract mRNA and lncRNA from `liftoff.gff` using `grep`. If you used a GFF3 file for LiftOff (i.e. rather than a GTF file), then the features we wish to extract will also have "mRNA" in column 2. We can also extract "lnc_RNA" features, just in case a coding region ("CDS" feature) was wrongly assigned to one of Mikado's gene models, but it's actually a lncRNA that also got captured by LiftOff. The `|` acts as an "or" statement, separating the two patterns we wish to extract and `-P` works as described above.

```
grep -P "\tmRNA\t|\tlnc_RNA\t" liftoff.gff > liftoff.mRNA.lncRNA.gff
```

Now we will extract the desired features from `toga.gff`. TOGA only retains protein-coding genes, and we don't have a feature type in column 2. However, all `toga.gff` contains are gene and exon features, with only one transcript-per-exon assumed. TOGA does indicate if a feature is a gene by having "ID" and "geneID" attributes in column 9 as metadata. We can therefore isolate any lines that contain "geneID="

```
grep "geneID=" toga.gff > toga.mRNA.gff
```

Now that we have isolated the transcripts from Mikado, LiftOff, and TOGA, we can use BEDTools to determine which gene models overlap with each other. We will determine (1) which LiftOff gene models overlap with which Mikado gene models to assign LiftOff gene symbols to the Mikado gene models, and (2) which TOGA gene models overlap with which Mikado gene models to assign gene symbols from TOGA. We cab deternube these overlaps using `bedtools intersect` with `-a` pointing at the Mikado GFF file, and `-b` pointing at the LiftOff and TOGA gene models. `-wo` specifies that the output will contain both A and B entries with the number of base pairs overlap between the two features, with only A features with overlap being reported.

First LiftOff:

```
bedtools intersect -a mikado.mRNA.gff -b liftoff.mRNA.lncRNA.gff -wo > mikado.liftoff.mRNA.lncRNA.txt
```

Now TOGA:

```
bedtools intersect -a mikado.mRNA.gff -b toga.mRNA.gff -wo > mikado.toga.mRNA.txt
```

### OrthoFinder

OrthoFinder, a tool that maps sequence-similarity relationships between proteins across two or more species based on their sequences, can also be used to identify predicted orthologs. OrthoFinder builds gene trees, considers gene duplication events, is considered to be one of the most accurate orthologs inference methods49, and was used for gene naming in the DNA zoo annotation project. OrthoFinder outputs lists of protein-protein sequence-similarity relationships that can be used to infer gene-gene relationships.
