# 5. Assign gene symbols to protein-coding genes

At this point in the tutorial, we now have a GFF file that contains the gene-models of protein-coding and non-coding RNA genes, and predicted gene symbols and products for the latter. This leaves the final step in the tutorial, which is assigning gene symbols and predicting the function of protein-coding genes. We can use sequence similarity between our target species and closely-related species to identify predicted orthologs of the gene models from the target species. Assigning gene symbols is challenging because most genes in a mammalian genome originated from another gene (e.g. tandem duplication, gene fusion, translocation), meaning that many genes have at least one paralogous gene with high sequence similarity in exons. Three tools can be used to predict gene symbols in the target species: LiftOff, TOGA, and [OrthoFinder](https://github.com/davidemms/OrthoFinder).

### LiftOff and TOGA

You may have already used LiftOff and TOGA earlier in the tutorial, making it fairly straightforward to use their results. Both LiftOff and TOGA annotate the target species’ gene structures and assign reference gene symbols to the target with a high rate of agreement with the gene symbols found in Ensembl annotations. Therefore, these tools can be used to predict gene symbols for the final, integrated annotation. We transfer gene models by overlapping transcripts derived from TOGA and LiftOff to the final gene models, before transfering the gene symbol to the Mikado-filtered gene identifier (ID).

### OrthoFinder

OrthoFinder, a tool that maps sequence-similarity relationships between proteins across two or more species based on their sequences, can also be used to identify predicted orthologs. OrthoFinder builds gene trees, considers gene duplication events, is considered to be one of the most accurate orthologs inference methods49, and was used for gene naming in the DNA zoo annotation project. OrthoFinder outputs lists of protein-protein sequence-similarity relationships that can be used to infer gene-gene relationships.
