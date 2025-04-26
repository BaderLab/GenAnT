library(optparse)
library(rtracklayer)

option_list <- list(
  make_option(c("-g", "--gff"), type = "character", default = "FALSE",
              help = "The gff file for the reference assembly",
              metavar = "character")
)

# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

gffname <- opt[[1]]

prefix <- gsub(".gff3","",gffname)

gff <- readGFF(gffname)

gff_transcript <- gff_transcript <- gff[gff$type == "transcript",]

key <- gff_transcript[,c("gene_id","transcript_id","gene_name")]
colnames(key) <- c("geneID","transcriptID","geneSymbol")


key$geneSymbol <- sub("-(?!.*-).*", "", key$geneSymbol,perl = TRUE)

write.table(key[,c("geneSymbol","transcriptID")],file=paste0(prefix,".table.txt"),
            quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

write.table(key[,c("geneID","transcriptID")],file=paste0(prefix,".isoforms.tsv"),
            quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

write.table(key,file=paste0(prefix,".genekkey.tsv"),
            quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
