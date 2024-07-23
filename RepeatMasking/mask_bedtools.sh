# You can use the output of Earl Grey to mask your FASTA file
# The `-soft` flag indicates soft masking, and the inputs are indicated with `-fi` for the input FASTA file,  and `-bed` for the GFF or BED file of transposable elements
# The output is your masked FASTA file, with the name indicated by `-fo`.

bedtools maskfasta \
 -soft \
 -fi your_genome.fasta \
 -bed ranges_to_mask.gff \
 -fo your_genome_softmasked.fasta
