library(optparse)
library(rtracklayer)

option_list <- list(
  make_option(c("-g", "--gff"), type = "character", default = "FALSE",
              help = "The gff file for the referencea assembly",
              metavar = "character")
)

# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

gffname <- opt[[1]]

prefix <- gsub(".gff","",gffname)

gff <- readGFF(gffname)

gff_transcript <- gff_transcript <- gff[gff$type == "transcript" | gff$type == "mRNA" ,]

# Get gene ID and transcript ID keys
# if this is a character already it doesn't change anything, but Rtracklayer sometimes loads "Parent" in as a list.

gff_transcript$gene_id <- unlist(gff_transcript$Parent)
gff_transcript$trans_id <- unlist(gff_transcript$ID)

key <- gff_transcript[,c("gene_id","trans_id","gene")]
colnames(key) <- c("geneID","transcriptID","geneSymbol")

# Remove deleting the `gene-` and `rna-` tag as it makes an incompatibility with toga isoforms

# key$geneID <- gsub("gene-","",key$geneID)
# key$transcriptID <- gsub("rna-","",key$transcriptID)

write.table(key[,c("geneSymbol","transcriptID")],file=paste0(prefix,".table.txt"),
            quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

write.table(key[,c("geneID","transcriptID")],file=paste0(prefix,".isoforms.tsv"),
            quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

write.table(key,file=paste0(prefix,".genekkey.tsv"),
            quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
