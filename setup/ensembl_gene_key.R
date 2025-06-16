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



gff_transcript <- gff_transcript <- gff[gff$type == "transcript" | gff$type == "mRNA",]

# gff_transcript <- apply(gff_transcript,2,unlist)

key <- gff_transcript[,c("Parent","ID","Name")]

key <- apply(key,2,unlist)

colnames(key) <- c("geneID","transcriptID","geneSymbol")

key <- as.data.frame(key,stringAsFactors=FALSE)

key$geneSymbol <- sub("-(?!.*-).*", "", key$geneSymbol,perl = TRUE)

write.table(key[,c("geneSymbol","transcriptID")],file=paste0(prefix,".table.txt"),
            quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

write.table(key[,c("geneID","transcriptID")],file=paste0(prefix,".isoforms.tsv"),
            quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

write.table(key,file=paste0(prefix,".genekkey.tsv"),
            quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
