fwdBlast <- read.table("forward_blast.tsv",header=FALSE,as.is=TRUE,sep="\t")

revBlast <- read.table("reverse_blast.tsv",header=FALSE,as.is=TRUE,sep="\t")

inter <- intersect(fwdBlast$V1,revBlast$V2)

reciprocalBlast <- fwdBlast[fwdBlast$V1 %in% inter,]


write.table(reciprocalBlast, file = "reciprocal_blast.tsv",quote=FALSE,row.names=FALSE,col.names = FALSE,sep="\t")

reciprocal_blast <- read.table("reciprocal_blast.tsv",header=F,as.is=T,sep="\t")

# parent gene

parent_key <- c()
for(i in unique(reciprocal_blast$V1)) {
  g1 <- reciprocal_blast[reciprocal_blast$V1 == i,]
  if(nrow(g1) > 1) {
    g1 <- g1[which.min(g1$V11),]
  }
  gene <- lapply(strsplit(g1$V2,"\\|"), function(x) return(x[3]))[[1]]
  gene <- lapply(strsplit(gene,"_"),function(x)return(x[1]))[[1]]
  parent_key <- rbind(parent_key,c(g1$V1,gene))
}

colnames(parent_key) <- c("parent","symbol")

write.table(parent_key, file = "reciprocal_blast_parent_symbol.txt", quote=FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


