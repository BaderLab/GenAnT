library(rtracklayer)
library(dplyr)

mikado <- as.data.frame(readGFF("full_annotation.gff"))

mikado <- mikado[mikado$type %in% c("gene", "lncRNA_gene"),]

mikado_df <- data.frame(mikado_id = mikado$ID, ncRNA_gene = mikado$predicted_gene_symbol)

liftoff_mikadoInfo <- as.data.frame(readGFF("liftoff_overlap.mikadoInfo.gff"))
liftoff_liftoffInfo <- as.data.frame(readGFF("liftoff_overlap.liftoffInfo.gff"))

liftoff_liftoffInfo$gene <- unlist(liftoff_liftoffInfo$Parent)

liftoff_df <- data.frame(mikado_id = liftoff_mikadoInfo$ID, liftoff_gene = liftoff_liftoffInfo$gene)

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
toga_togaInfo <- as.data.frame(readGFF(paste0("toga_overlap.",i,".togaInfo.gff")))

if(nrow(toga_togaInfo) < 5) next

toga_df <- data.frame(mikado_id = toga_mikadoInfo$ID, toga_gene = toga_togaInfo$ID)

toga_df$toga_gene <- gsub("\\.[^.]*$", "", toga_df$toga_gene)

toga_df$mikado_id <- gsub("\\.[^.]*$", "", toga_df$mikado_id)

#
# reference.toga.r1.table.txt
reference <- read.delim(paste0("reference.toga.",i,".table.txt"), header = TRUE, sep = "\t")

colnames(reference) <- c("geneSymbol", "transcriptID")
rownames(reference) <- reference$transcriptID

toga_ids <- toga_df$toga_gene

toga_gene <- unlist(lapply(toga_ids,function(x) reference[x,"geneSymbol"]))

toga_df$toga_gene <- toga_gene


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
