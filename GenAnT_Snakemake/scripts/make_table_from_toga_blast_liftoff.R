setwd("/Users/dsokolowski/Desktop/NMR_multi_assembly/annotation_zoe/t2t_symbol/")

library(rtracklayer)
library(pbapply)
mikado <- readGFF("mikado_lenient.gff")

toga <- read.table("mikado.toga.mRNA.txt",header=FALSE,as.is=TRUE,sep="\t")
liftoff <- read.table("mikado.liftoff.mRNA.txt",header=FALSE,as.is=TRUE,sep="\t")

gff_input <- list(toga = toga, liftoff = liftoff)

symbol_output <- list()
for(n in names(gff_input)) {
  
  gff_olap <- gff_input[[n]]
  getID <- function(y) {
    z <- y[grep("ID=",y)]
    out <- gsub("ID=","",z)
    return(out)
  }
  mikado_id <- unlist(lapply(strsplit(gff_olap$V9,";"),getID))
  mmus_ID <- unlist(lapply(strsplit(gff_olap$V18,";"),getID))
  symbol_olap <- as.data.frame(cbind(mikado_id,mmus_ID,gff_olap$V19),stringsAsFactors=FALSE)
  colnames(symbol_olap) <- c("ID","Ref","overlap")
  symbol_olap$overlap <- as.numeric(symbol_olap$overlap)
  
  unique <- unique(symbol_olap$ID)
  symbol_unique <- c()
  symbol_assign <- do.call("rbind",lapply(as.list(unique), function(x) {
    ol_in <- symbol_olap[symbol_olap$ID == x,]
    ol_out <- ol_in[which.max(ol_in$overlap)[1],]
    return(ol_out)
  }))
  symbol_assign$Ref <- unlist(lapply(strsplit(symbol_assign$Ref,"\\."), function(x) return(x[1])))
  symbol_output[[n]] <- symbol_assign
  print(n)
}

toga_fin <- symbol_output$toga
liftoff_fin <- symbol_output$liftoff

recip_blast <- read.table("reciprocal_blast_combined.txt", header=FALSE, as.is=TRUE,sep="\t")

colnames(recip_blast) <- c("parent","symbol")

mikado_mRNA <- mikado[mikado$type == "mRNA",]
gene_symbol <- matrix("",nrow = nrow(mikado_mRNA), ncol=4)
rownames(gene_symbol) <- mikado_mRNA$ID

colnames(gene_symbol) <- c("parent","blast","toga","liftoff")
gene_symbol <- as.data.frame(gene_symbol)
gene_symbol$parent <- unlist(mikado_mRNA$Parent)


rownames(recip_blast) <- recip_blast$parent
gene_symbol[rownames(recip_blast),"blast"] <- recip_blast$symbol

rownames(toga_fin) <- toga_fin$ID
gene_symbol[rownames(toga_fin),"toga"] <- toga_fin$Ref

rownames(liftoff_fin) <- liftoff_fin$ID
gene_symbol[rownames(liftoff_fin),"liftoff"] <- liftoff_fin$Ref

# mmus to gene

mmusID <- read.table("mmus_gene_ID_key.txt",header=T,as.is=T,sep="\t")
mmusID <- mmusID[!duplicated(mmusID$Transcript.stable.ID),]

# HetGla to gene

59172+12948

liftoff_named <- readGFF("T2T_lift_named.gtf")
liftoff_transcript <- liftoff_named[liftoff_named$type == "transcript",]

liftID <- as.data.frame(cbind(liftoff_transcript$transcript_id,liftoff_transcript$gene_name))
colnames(liftID) <- c("Transcript.stable.ID","Gene.name")

write.table(liftID,file = "HetGla_liftOff_id.txt", quote=F,row.names=F,col.names=T,sep="\t")
write.table(mmusID, file = "mmus_id.txt",quote=F,row.names=F,col.names=T,sep="\t")
table(duplicated(liftoff_transcript$transcript_id))
gene_symbol$toga <- unlist(pblapply(gene_symbol$toga,function(x) {
  index <- which(x == mmusID$Transcript.stable.ID)
  if(length(index) > 0) {
    x <- mmusID[index,"Gene.name"]
  }
  return(x)

}))

gene_symbol$liftoff <- unlist(pblapply(gene_symbol$liftoff,function(x) {
  index <- which(x == liftID$Transcript.stable.ID)
  if(length(index) > 0) {
    x <- liftID[index,"Gene.name"]
  }
  return(x)
  
}))

gene_symbol$toga <- toupper(gene_symbol$toga)

gene_symbol$consensusNov11 <- gene_symbol$liftoff
gene_symbol$consensusNov11[is.na(gene_symbol$consensusNov11)] <- ""
gene_symbol$consensusNov11[gene_symbol$consensusNov11 == ""] <- gene_symbol$toga[gene_symbol$consensusNov11 == ""]
gene_symbol$consensusNov11[gene_symbol$consensusNov11 == ""] <- gene_symbol$blast[gene_symbol$consensusNov11 == ""]

mikado$Parent <- unlist(mikado$Parent)



mikado$ID2 <- gsub("mikado.","",mikado$ID)
mikado$ID2 <- unlist(lapply(strsplit(mikado$ID2,"\\."),function(x) return(x[1])))

gene_symbol$parent1 <- gsub("mikado.","",gene_symbol$parent)

for(i in unique(gene_symbol$parent1)) {
  gp1 <- which(i == mikado$ID2)
  mikado[gp1,"gene_name"] <- gene_symbol[gene_symbol$parent1 == i,"consensusNov11"][1]
  if((which(i==unique(gene_symbol$parent1)) %% 1000) == 0) print(i)
}

write.table(gene_symbol, file = "mikado_gene_symbol_table.txt", quote = FALSE,
            row.names=TRUE, col.names = TRUE,sep="\t")

rtracklayer::export.gff3(mikado,con = "mikado_lenient_symbol.gff")
mikado <- readGFF("mikado_lenient_symbol.gff")
mik_name = unique(mikado$gene_name)[!unique(mikado$gene_name) %in% unique(liftoff_named$gene_name)]
library(gprofiler2)
syms <- gost(mik_name,custom_bg = unique(mikado$gene_name),evcodes = T)
res1 <- syms$result[syms$result$source %in% c("GO:BP","GO:MF","GO:CC"),]

