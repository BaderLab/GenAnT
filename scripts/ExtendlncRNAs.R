# This notebook reads in Mikado gene models labeled by Infernal, and assigns Infernal gene names to all gene models.

## Read in file


library(rtracklayer)
library(dplyr)


# Read in GFF file of Mikado features that overlap with Infernal lncRNAs


mikado <- BiocGenerics::as.data.frame(rtracklayer::readGFF("mikado_lenient.gff"))
# Add metadata columns that will be filled with lncRNA and later short ncRNA
mikado$predicted_gene_symbol <- NA
mikado$infernal_product <- NA
mikado$gbkey <- NA
mikado$gene_biotype <- NA

head(mikado)


# Now read in Infernal lncRNA information 

infernal <- try(BiocGenerics::as.data.frame(read.delim("infernal.lncRNA.InMikado.gff",header=FALSE,as.is=TRUE)))

if(class(infernal)[[1]] != "try-error") {
  
  infernal <- infernal[infernal$V12 != "superlocus",]
  infernal <- infernal[infernal$V7 == infernal$V16,] # redundant, ensure strand
  infernal <- infernal[infernal$V1 == infernal$V10,] # redundant, ensure chr
  # Get Infernal IDs with a strand-specific overlap to the lncRNA
  
  
  infernalID <- gsub("ID=","",unlist(lapply(strsplit(infernal$V18,";"),function(x) return(x[1]))))
  
  # Get other relevant infernal information to add to mikado gff
  thesplit <- as.data.frame(do.call("rbind",strsplit(infernal$V9,";")))
  colnames(thesplit) <- c("ID","evalue","rfamid","product","gbkey")
  for(i in colnames(thesplit)) {
    thesplit[,i] <- gsub(paste0(i,"="),"",thesplit[,i])
  }
  thesplit$product <- gsub("\t", " ", thesplit$product)
  
  for(i in 1:length(infernalID)) {
    if(infernal$V4[i] < infernal$V13[i])  { # extend to the left 
      mikado[mikado$ID == infernalID[i],"start"] <- infernal$V4[i]
      message(paste0(infernalID[i], " is extended from Infernal lncRNA at start"))
    }
    if(infernal$V5[i] > infernal$V15[i])  { # extend to the right 
      mikado[mikado$ID == infernalID[i],"end"] <- infernal$V5[i]
      message(paste0(infernalID[i], " is extended from Infernal lncRNA at end"))
      
    }
    if(infernal$V12[i] %in% c("gene","mRNA","ncRNA_gene","ncRNA")) {
      
      # Add in correct metadata
      mikado$predicted_gene_symbol[i] <- thesplit$ID[i]
      mikado$infernal_product[i] <- thesplit$product[i]
      if(mikado$type[i] %in% c("gene","lncRNA_gene")) {
        mikado$gene_biotype[i] <- "lncRNA"
      } else {
        mikado$gbkey[i] <- "ncRNA"
      }
    }
  }
}
                                          

rtracklayer::export.gff3(mikado, "mikado.infernal.lncRNALabeled.polished.gff", format = "gff3")

###
####
###


mikado <- BiocGenerics::as.data.frame(rtracklayer::readGFF("infernal.lncRNA.notInMikado.gff"))


gene_name <- mikado[]

mikado[] <- lapply(mikado, function(x) if (is.factor(x)) as.character(x) else x)
mikado_test <- mikado[1:5,]
mikado_test$seqid <- "scaffold_4"
mikado <- as.data.frame(rbind(mikado,mikado_test))

for(i in unique(scaffold_id)) {
  mikado_id <- mikado[scaffold_id == i,]
}

if(nrow(mikado) < 3) {
  file.create("infernal.lncRNA.New.ExonID.gff")
}

mikado$predicted_gene_symbol <- unlist(lapply(strsplit(mikado$ID,"_"), function(x) return(x[1])))


mikado$infernal_product <- mikado$product
mikado$gene_biotype <- NA

head(mikado[mikado$source == "cmscan",])


# The first thing we want to do is make sure that all separate "cmscan" lncRNAs have a unique gene and transcript ID. We will assume that any of these lncRNAs that are on the same contig are actually multiple exons derived from the same lncRNA gene, and that lncRNAs on separate contigs are different genes. Therefore, we will loop through any lncRNAs that share the same ID, and determine if they are on the same or different contigs. If they are on the same contig, they will share the same ID; if they are on different contigs, they will get a unique ID.


# Isolate unique Infernal gene symbols with no Mikado overlap
genes <- unique(mikado$predicted_gene_symbol[mikado$source == "cmscan"])
genes <- genes[!(is.na(genes))]
if(length(genes) < 1) {
  file.create("infernal.lncRNA.New.ExonID.gff")
}
# Loop through genes
for(gene in genes) {
  # Gather the different contigs that the gene is on
  contigs <- unique(mikado$seqid[mikado$predicted_gene_symbol == gene])
  # If there is more than one contig, loop through the different contigs
  contigs <- contigs[!is.na(contigs)]
  if(length(contigs) > 1) {
    # Start a counter at 0
    i <- 0
    for(contig in contigs) {
      # Only add modifications at the second contig onward
      if(i > 0){
        # Create new unique ID name
        new_id <- paste(gene, "-copy", i, sep = "")
        # Replace the current ID for this gene and contig with new_id
        mikado$ID[mikado$ID == gene & mikado$seqid == contig] <- rep(new_id, length(mikado$ID[mikado$ID == gene & mikado$seqid == contig]))
      }
      # Increase counter assuming we are working with the next lncRNA copy
      i <- i + 1
    }
  }
}
# Look at all of the unique gene symbols to see if any copies were made
# (They will have "copy#" after the gene symbol)
unique(mikado$predicted_gene_symbol[mikado$source == "cmscan"])
mikado$type <- "gene"
mikado$ID <- make.unique(mikado$ID,sep="-copy")






# Now we need to create lncRNA and gene IDs for these features, especially combining the multi-exonic features. We can treat everything that currently exists as an exon, and create transcript and gene features above these features. Note that the following code assumes that every feature listed above has "lncRNA" in the "type" column, and "cmscan" in the source column (which should be true - this is just a sanity check).


# Add gbkey=NA, gene_biotype=NA, Parent=NA to all cmscan rows just in case there was no overlap with Mikado
# so that rbind will work
mikado$gbkey[mikado$source == "cmscan"] <- rep(NA, length(mikado$source == "cmscan"))
mikado$gene_biotype[mikado$source == "cmscan"] <- rep(NA, length(mikado$source == "cmscan"))
mikado$Parent[mikado$source == "cmscan"] <- rep(NA, length(mikado$source == "cmscan"))
# Isolate unique Infernal gene symbols with no Mikado overlap
genes <- unique(mikado$predicted_gene_symbol[mikado$source == "cmscan"])
genes <- genes[!is.na(genes)]
# Loop through genes
for(gene in genes) {
  # Gather the different contigs that the gene is on
  contigs <- unique(mikado$seqid[mikado$predicted_gene_symbol == gene])
  contigs <- contigs[!is.na(contigs)]
  # Loop through contigs
  for(contig in contigs){
    # Get indexing in terms of where these occur in the dataframe to reduce code
    i <- mikado$predicted_gene_symbol == gene & mikado$seqid == contig
    # First, duplicate one of the "exon" rows - this will become our gene feature
    gene_feat <- mikado[i,][1,]
    # Change type to lncRNA_gene
    gene_feat$type <- "gene"
    # Add gene_biotype=lncRNA
    gene_feat$gene_biotype <- "lncRNA"
    # Make the start value the lowest of all "exon" start values
    gene_feat$start <- min(mikado$start[i])
    # Make the end value the highest of all "exon" end values
    gene_feat$end <- max(mikado$end[i])
    # Create a transcript feature by copying the gene feature
    transcript <- gene_feat
    # Make the gbkey ncRNA and set gene_biotype=NA
    transcript$gbkey <- "ncRNA"
    transcript$gene_biotype <- NA
    # Add .1 to the transcript ID
    transcript$ID <- paste(transcript$ID, ".1", sep = "")
    # Add Parent pointing to gene
    transcript$Parent <- gene_feat$ID
    # Change the type back to lncRNA
    transcript$type <- "lncRNA"
    # Now change the lncRNA designation to exon for the specified gene symbol and contig in the loop
    mikado$type[i] <- rep("exon", length(mikado$type[i]))
    # Add a "Parent ID" pointing to the new transcript ID that we created
    mikado$Parent[i] <- rep(transcript$ID, length(mikado$type[i]))
    # Now change ID to end with .exon1, .exon2, etc.
    for(exon in 1:sum(i)) {
      # Create new ID with paste
      mikado$ID[which(i)[exon]] <- paste(mikado$ID[which(i)[exon]],".1.exon",exon, sep = "")
    }
    # Add gene and transcript features to mikado dataframe, inserting it just before the "exons"
    # First line needs an exception to prevent duplicating an exon
    if(which(i)[1] - 1 == 0) {
      mikado <- rbind(gene_feat,
                      transcript,
                      mikado[which(i)[1]:nrow(mikado), ])
    } else {
      mikado <- rbind(mikado[1:(which(i)[1] - 1), ],
                      gene_feat,
                      transcript,
                      mikado[which(i)[1]:nrow(mikado), ])
    }
  }
}

rtracklayer::export.gff3(mikado, "infernal.lncRNA.New.ExonID.gff", format = "gff3")

