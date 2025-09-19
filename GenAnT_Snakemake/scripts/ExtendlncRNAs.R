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

infernal <- BiocGenerics::as.data.frame(read.delim("infernal.lncRNA.InMikado.gff",header=FALSE,as.is=TRUE))
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

rtracklayer::export.gff3(mikado, "mikado.infernal.lncRNALabeled.polished.gff", format = "gff3")

