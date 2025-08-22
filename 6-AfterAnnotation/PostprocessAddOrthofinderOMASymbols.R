library(plyr)
library(tidyverse)

# Parse arguments

ortho_out <- paste0("orthofinder_protein.tsv") 
mikado_df <- read.table("gene_symbols_noCopies.tsv", sep = "\t", header = TRUE)


orthofinder_pairs <- read.table(ortho_out,
                                header = TRUE,
                                sep = "\t")
colnames(orthofinder_pairs)[2:3] <- c("mikado_id","OrthoOMA_gene")
head(orthofinder_pairs)

# Remove Orthogroup column, as it is not needed

orthofinder_pairs$Orthogroup <- NULL

# In each column, you will see many protein or RNA feature names grouped together by commas in each row (if you followed the tutorial up to this point, you are likely looking at human RNA feature names that are formatted like "rna-id" in the right column). We need to expand by those values to get one protein name in each column.

orthofinder_pairs <- tidyr::separate_longer_delim(orthofinder_pairs, 
                                                  "mikado_id", 
                                                  delim = ", ")
orthofinder_pairs <- tidyr::separate_longer_delim(orthofinder_pairs, 
                                                  "OrthoOMA_gene", 
                                                  delim = ", ")

# Remove decimals from "mikado_id" column as they represent splice variants and come from the same gene and remove duplicate rows

orthofinder_pairs$mikado_id <- gsub("\\.[0-9]+$","",
                                    orthofinder_pairs$mikado_id)
orthofinder_pairs$OrthoOMA_gene <- gsub("\\.[0-9]+$","",
                                           orthofinder_pairs$OrthoOMA_gene)
orthofinder_pairs <- dplyr::distinct(orthofinder_pairs)

## Gene name conversion

reference <- read.delim("reference.table.txt", header = TRUE, sep = "\t")

colnames(reference) <- c("gene", "OrthoOMA_gene")
reference$OrthoOMA_gene <- gsub("\\.[0-9]+$","",
                                   reference$OrthoOMA_gene)

head(reference)

# Use `left_join` from `dplyr` to add the "reference" data frame to the "orthofinder_pairs" data frame

orthofinder_pairs <- dplyr::left_join(orthofinder_pairs, reference, by = "OrthoOMA_gene")
head(orthofinder_pairs)

# We can now replace the "OrthoOMA_gene" column with "gene", and then remove the "gene" column. After this, we will also remove redundant rows since different transcript names may share the same gene name.

orthofinder_pairs$OrthoOMA_gene <- orthofinder_pairs$gene
orthofinder_pairs$gene <- NULL
orthofinder_pairs <- dplyr::distinct(orthofinder_pairs)
head(orthofinder_pairs)

# Same as with the LiftOff and TOGA results, we will combine multiple gene symbols that represent the same Mikado ID with semi-colons, and add copy numbers to identical gene symbols that represent different Mikado IDs.

orthofinder_pairs <- orthofinder_pairs[order(orthofinder_pairs$mikado_id,
                                             orthofinder_pairs$OrthoOMA_gene),]

# Now collapse by Mikado ID and add the semi-colons if necessary:

orthofinder_pairs <- plyr::ddply(orthofinder_pairs,
                                 "mikado_id",
                                 summarize,
                                 OrthoOMA_gene = paste(OrthoOMA_gene, collapse = ";"))
head(orthofinder_pairs)

orthofinder_pairs$OrthoOMA_gene <- make.unique(orthofinder_pairs$OrthoOMA_gene, sep = "-copy")

## Add to gene symbol table



# Use `left_join` to add the new gene symbols to both tables, removing the copies from the "no copies" table

mikado_df <- dplyr::left_join(mikado_df,
                              orthofinder_pairs,
                              by = "mikado_id")

write.table(mikado_df, file = "postprocess_gene_symbols_OrthoOMA.tsv",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
