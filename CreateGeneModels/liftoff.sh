# Liftoff lifts a genome annotation from a closely-related species to your species
# Input is the GFF or GTF and FASTA file from a closely-related species, as well as the FASTA file of your species
# Output is a GFF annotation file for your species, a list of unmapped genes, and some intermediate files used by the program
# Note that `-copies` looks for gene copies in your species.

liftoff \
 -g annotation_of_related_species.gff \
 your_genome.fasta \
 genome_of_related_species.fasta \
 -o output_annotation.gff \
 -u unmapped_features.txt \
 -copies\
 -p number_of_threads
