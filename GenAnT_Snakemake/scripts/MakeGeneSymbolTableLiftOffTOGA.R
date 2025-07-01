library(rtracklayer)
library(dplyr)

mikado <- as.data.frame(readGFF("full_annotation.gff"))

mikado <- mikado[mikado$type %in% c("exon"),]


# load overlapping genes
liftoff_mikadoInfo <- as.data.frame(readGFF("liftoff_overlap.mikadoInfo.gff"))

# gene models from exon models
liftoff_mikadoInfo$mikadoGene <- unlist(liftoff_mikadoInfo$Parent)
liftoff_mikadoInfo$mikadoGene <- sub("\\.[^\\.]*$", "", liftoff_mikadoInfo$mikadoGene)

# get transcripts that overlap
liftoff_liftoffInfo <- as.data.frame(readGFF("liftoff_overlap.liftoffInfo.gff"))

# Lmao ensembl
liftoff_liftoffInfo$liftOffGene <- gsub("transcript:","",unlist(liftoff_liftoffInfo$Parent))

# pair matching exons into matrix
liftOffGene <- as.data.frame(cbind(liftoff_mikadoInfo$mikadoGene,liftoff_liftoffInfo$liftOffGene),stringsAsFactors=FALSE)

# make them one element
combo <- paste0(liftOffGene$V1,":",liftOffGene$V2)
liftOffGene$combo <- combo

# how many exons overlap per gene pair from mikado and liftoff
db1 <- as.data.frame(table(combo),stringsAsFactors=FALSE)
db_genes <- as.data.frame(do.call("rbind",strsplit(db1$combo,":")), stringsAsFactors=FALSE)
db_genes$combo <- db1$combo
db_genes$Frequency <- db1$Freq

# extract the mikado-liftoff pair with the highest number of exons overlapping
db_genes <- db_genes[order(db_genes$Frequency,decreasing = TRUE),]
db_genes_noDup <- db_genes[!duplicated(db_genes$V1),]

liftoff_df <- as.data.frame(db_genes_noDup[,1:2])
colanmes(liftoff_df) <- colnames("ID", "gene")

#liftoff_df <- data.frame(mikado_id = liftoff_mikadoInfo$ID, liftoff_gene = liftoff_liftoffInfo$gene)

liftoff_df$mikado_id <- gsub("\\.[^.]*$", "", liftoff_df$mikado_id)

liftoff_df <- dplyr::distinct(liftoff_df)

liftoff_df <- liftoff_df[order(liftoff_df$mikado_id,
                               liftoff_df$liftoff_gene),]

liftoff_df <- plyr::ddply(liftoff_df,
                          "mikado_id",
                          summarize,
                          liftoff_gene = paste(liftoff_gene, collapse = ";"))

liftoff_df$liftoff_gene <- make.unique(liftoff_df$liftoff_gene, sep = "-copy")

liftoff_df$liftoff_gene <- gsub("gene-","",liftoff_df$liftoff_gene)

mikado_df <- dplyr::left_join(mikado_df,
                              liftoff_df,
                              by = "mikado_id")

toga_dfs <- list()
for(i in c("r1","r2")) {

toga_mikadoInfo <- as.data.frame(readGFF(paste0("toga_overlap.",i,".mikadoInfo.gff")))

toga_mikadoInfo$mikadoGene <- unlist(toga_mikadoInfo$Parent)
toga_mikadoInfo$mikadoGene <- sub("\\.[^\\.]*$", "", toga_mikadoInfo$mikadoGene)

toga_togaInfo <- as.data.frame(readGFF(paste0("toga_overlap.",i,".togaInfo.gff")))
toga_togaInfo$togaGene <- gsub("transcript:","",unlist(toga_togaInfo$Parent))

if(nrow(toga_togaInfo) < 5) next

togaGene <- as.data.frame(cbind(toga_mikadoInfo$mikadoGene,toga_togaInfo$togaGene),stringsAsFactors=FALSE)

combo <- paste0(togaGene$V1,":",togaGene$V2)
togaGene$combo <- combo

toga_db1 <- as.data.frame(table(combo),stringsAsFactors=FALSE)
toga_db_genes <- as.data.frame(do.call("rbind",strsplit(toga_db1$combo,":")), stringsAsFactors=FALSE)
toga_db_genes$combo <- toga_db1$combo
toga_db_genes$Frequency <- toga_db1$Freq

toga_db_genes <- toga_db_genes[order(toga_db_genes$Frequency,decreasing = TRUE),]

toga_db_genes_noDup <- toga_db_genes[!duplicated(toga_db_genes$V1),]

toga_mikadoInfo <- as.data.frame(toga_mikadoInfo[,1:2])

colnames(toga_mikadoInfo) <- c("ID", "gene")
##
###
##

##
###
##

toga_df <- data.frame(mikado_id = toga_mikadoInfo$ID, toga_gene = toga_togaInfo$gene)

toga_df$toga_gene <- gsub("\\.[^.]*$", "", toga_df$toga_gene)

toga_df$mikado_id <- gsub("\\.[^.]*$", "", toga_df$mikado_id)



##


toga_df <- dplyr::distinct(toga_df)


toga_df <- toga_df[order(toga_df$mikado_id,
                         toga_df$toga_gene),]

toga_df <- plyr::ddply(toga_df,
                          "mikado_id",
                          summarize,
                          toga_gene = paste(toga_gene, collapse = ";"))

toga_df[[paste0("toga",i,"_gene")]] <- make.unique(toga_df$toga_gene, sep = "-copy")

toga_df[[paste0("toga",i,"_gene")]] <- gsub("gene-","",toga_df$toga_gene)

toga_dfs[[i]] <- toga_df

}

##

for(i in 1:length(toga_dfs)) {
 toga_df <- toga_dfs[[i]]
 mikado_df <- dplyr::left_join(mikado_df,
                               toga_df,
                               by = "mikado_id")
}



write.table(mikado_df, file = "gene_symbols.tsv",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

mikado_df_noCopies <- mikado_df
for(i in grep("_gene",colnames(mikado_df_noCopies))) {
  mikado_df_noCopies[,i] <- gsub("-copy.*", "", mikado_df_noCopies[,i] )
}

write.table(mikado_df_noCopies, file = "gene_symbols_noCopies.tsv",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
