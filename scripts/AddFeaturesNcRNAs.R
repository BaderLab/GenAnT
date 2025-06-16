library(rtracklayer)
library(dplyr)

##
###   This script takes the combined output of short non-coding RNAs from MirMachine and Infernal, and adds gene and exon features to the RNA features already described.
##

gff <- BiocGenerics::as.data.frame( rtracklayer::readGFF("short_ncRNAs.noOverlap.gff"))



# First, remove ".PRE" from the end of the MirMachine IDs


gff$ID <- gsub("\\.PRE$", "", gff$ID)


# Now we want to make sure that each ID is unique, while also preserving the predicted gene symbol. To do this, we'll create a new column called "predicted_gene_symbol" which is a copy of the "ID" column, and then make the IDs unique by adding "-copy#" to however many copies of that particular ID there are.


gff$predicted_gene_symbol <- gff$ID
gff$ID <- make.unique(gff$ID, sep = "-copy")

# Finally, we'll add our gene and exon features and stick these into the appropriate rows. We'll do this by looping through the different rows of the GFF, copying the current row that gene ID is located in, and making appropriate modifications.

# For gene rows:
# - The "type" attribute is changed to "ncRNA_gene"
# - A "gene_biotype" column is added becomes whatever the transcript type is (e.g. "microRNA")
# - The "gbkey" attribute is switched to NA

# For transcript rows:
# - The "ID" gets ".1" added to it
# - A "Parent" column points to the gene ID

# For exon rows:
# - The "ID" gets ".1.exon1" added to it
# - A "Parent" column points to the transcript ID
# - The "gbkey" attribute is switched to NA
# - The "type" attribute is changed to "exon"

# Initialize features that do not yet exist
gff$Parent <- NA
gff$gene_biotype <- NA
# Loop through gene IDs
for(id in gff$ID) {
  # Determine row gene ID is in
  i <- which(gff$ID == id)
  # Extract that row
  row <- gff[i,]
  # Mapy a copy to create gene feature and assign attributes
  gene <- row
  gene$type <- "ncRNA_gene"
  gene$gene_biotype <- row$type
  gene$gbkey <- NA
  # Make a copy to create transcript feature and assign atrributes
  transcript <- row
  transcript$ID <- paste(transcript$ID, ".1", sep = "")
  transcript$Parent <- gene$ID
  # Make a copy to create exon feature and assign attributes
  exon <- row
  exon$ID <- paste(transcript$ID, ".exon1", sep = "")
  exon$Parent <- transcript$ID
  exon$gbkey <- NA
  exon$type <- "exon"
  # Add back to GFF
  if(i == 1) {
    gff <- rbind(gene,
                 transcript,
                 exon,
                 gff[(i+1):nrow(gff),])
  } else {
    gff <- rbind(gff[1:(i-1),],
                 gene,
                 transcript,
                 exon,
                 gff[(i+1):nrow(gff),])
  }
}
head(gff)


# Reorder so that the "Parent" column is second in the metadata


gff <- gff %>%
  dplyr::relocate(Parent, .after = ID)


# Remove any rows of NAs

gff <- gff[!is.na(gff$start),]


# The GFF file can now be saved


rtracklayer::export.gff3(gff, "short_ncRNAs.polished.gff", format = "gff3")

