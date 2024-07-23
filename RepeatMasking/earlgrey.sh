# Repeat masking can be done using the tools Earl Grey and BedTools, used sequentially
# Input to Earl Grey is the FASTA file from your species, the name of your species, and the output directory
# It is also helpful to specify the search term used for RepeatMasker with `-r`, which indicates which set of repeats to look for (e.g. “arthropoda”) 
# These options are found in the RepeatMasker documentation
# The output is stored in multiple folders, with the most important information located in `summaryFiles`

earlGrey \
 -g your_genome.fasta \
 -s your_species_name \
 -o ./output_directory \
 -r repeat_clade \
 -d yes \
 -t number_of_threads
 
# One of the outputs of Earl Grey is a GFF3 and BED file of transposable element annotations
# If you have access to such a file, you can use it to mask a genome with `bedtools maskfasta`, which quickly creates a masked FASTA file from the GFF3 or BED annotations
