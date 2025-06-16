# This R notebook reads in `full_annotation.gff` and the table of gene symbols, `gene_symbols_noCopies.tsv`, and creates a clean GFF file with gene symbols. Many metadata columns are removed that may have been useful in intermediate steps but that are likely unnecessary going forward.


library(rtracklayer)
library(dplyr)


# First, read in the GFF file


gff <- as.data.frame(readGFF("full_annotation.gff"))


# Let's start by only keeping the most important columns. You can add any additional attributes you wish to keep, but much of this information was used to show evidence for the support of each gene model which you likely won't need directly in the GFF (and can instead find this information in some of the intermediate files, like tables output by Mikado). We'll definitely keep columns 1-to-8, in addition to the attributes listed below. The `which` function pulls the column numbers from the GFF file if the column names occur in the vector below. Note that we are choosing not to keep "predicted_gene_symbol" which has the predicted ncRNA genes from Infernal and MirMachine because we're just going to replace those with the same ncRNA gene stored in the ncRNA column of the gene symbol table.


gff <- gff[,c(1:8, which(colnames(gff) %in%
                           c("ID", "Name", "Parent", "alias", "primary",
                             "evalue","RFamID","product","gbkey",
                             "gene_biotype","infernal_product")))]


# Mikado also creates a feature called "superlocus" which groups together features based on their location in the genome. "Superlocus" is not a usual feature and may cause issues for you in the future, so you may which to remove superlocus rows from the GFF file. To do this, we can look for any rows where the "type" column matches the word "superlocus", and use the `!` symbol to indicate that we don't want to include rows with that feature.


gff <- gff[!gff$type == "superlocus",]



# Now let's read in `gene_symbols_noCopies.tsv`, which we will use to add gene symbols to the GFF file. We want the "no copies" file, because we can choose exactly how we want to add the copy numbers when considering all of the different gene symbol sources.


symbols <- read.table("gene_symbols_noCopies.tsv", sep = "\t", header = TRUE)


# The goal here is to populate the "Name" slot, as this is what will be recognized by other tools that interpret GFF files as the name of that particular gene. Right now, "Name" is almost always populated with a unique "Mikado" ID, and we want to replace it with a gene symbol.

# One strategy to do this is a "hierarchical" strategy, where you rank the different columns in the gene symbol table from best to worst, and you prioritise the "best" gene symbols to first fill in the "Name" slot. If there is an NA for a particular Mikado ID, you can then use the second-best column, etc. until all the columns have been accessed.

# The first step to this will be to create...


symbols$hierarchical_gene <- rep(NA, nrow(symbols))
for (j in 1:nrow(symbols)) {
  # First: if there is an OrthoFinder gene in table, use it
  if (!is.na(symbols$orthofinder_gene[j])) {
    symbols$hierarchical_gene[j] <- as.character(symbols$orthofinder_gene[j])
  }
  # If no OrthoFinder gene, use TOGA gene
  else if (!is.na(symbols$togar1_gene[j])) {
    symbols$hierarchical_gene[j] <- as.character(symbols$togar1_gene[j])
  }
  else if (!is.null(symbols$togar2_gene[j])) { # make sure that the column exists.
    if(!is.na(symbols$togar2_gene[j]))  symbols$hierarchical_gene[j] <- as.character(symbols$togar2_gene[j])
  }
  # If no TOGA gene, use LiftOff gene
  else if (!is.na(symbols$liftoff_gene[j])) {
    symbols$hierarchical_gene[j] <- as.character(symbols$liftoff_gene[j])
  }
  # If no LiftOff gene, use ncRNA gene
  else if (!is.na(symbols$ncRNA_gene[j])) {
    symbols$hierarchical_gene[j] <- as.character(symbols$ncRNA_gene[j])
  }
}


# An alternative to this hierarchical approach is to just pick your favourite column and use the gene names from there.

# Either way, let's assume that you now want to use the `hierarchical_gene` column to populate your gene names in the GFF file. Importantly, this gene name should be unique or else problems may arise in downstream analysis (since many downstream tools expect unique gene symbols, as exist in human or mouse). We can add a column to `symbols` which we'll call "unique_gene", and make the genes unique by adding `-copy#` to repeating gene symbols. We can do this using the `make.unique` function in R. This function sometimes adds unique identifiers to `NA` values, so we'll first grab indices for only non-NA values.


#
##  The below loop simplifies any orthofinder results in (e.g., hierarchical_gene) using their overlap with the gene models from TOGA and liftOff
#

ortho_copies <- grep(";",symbols$hierarchical_gene)

symbols_simplify <- symbols[ortho_copies,]

liftoff_simple <- strsplit(symbols_simplify$liftoff_gene,";")
toga_simple <- strsplit(symbols_simplify$togar1_gene,";")
if("togar2_gene" %in% colnames(symbols)) {
  togar2_simple <- strsplit(symbols_simplify$togar1_gene,";")
}
hier_simple <- strsplit(symbols_simplify$hierarchical_gene,";")

simple <- c()

for(i in 1:length(liftoff_simple)) {
  # Get all potential gene names from each method
  hier <- hier_simple[[i]]
  lift <- liftoff_simple[[i]]
  tog <- toga_simple[[i]]
  genes <- list(lift=lift,tog=tog)
  if("togar2_gene" %in% colnames(symbols)) {
    tog2 <- togar2_simple[[i]]
    genes[["tog2"]] <- tog2
  }
  

  # get number of times each gene shows up
  all_elements <- unique(unlist(genes))
  
  # Count in how many vectors each element appears
  element_counts <- sapply(all_elements, function(x) {
    sum(sapply(genes, function(v) x %in% v))
  })
  
  
  if(!any(hier %in% all_elements)) {
    # If the Orthofinder output never shows up then tough luck
    simple[i] <- symbols_simplify$hierarchical_gene[i]
  } else {
    # if it does, then we reduce the orthofinder genes symbols to those that have the greatest overlap with the other methods.
    for(i in sort(unique(element_counts),decreasing = FALSE)) {
     inter <- intersect(hier,names(element_counts)[element_counts == i])
     if(length(inter) > 0) {
       out1 <- hier[hier %in% tog]
       simple[i] <- paste(out1, collapse = ";")
     }
    }
  }
}
  

symbols$hierarchical_gene[ortho_copies] <- simple

symbols$unique_gene <- symbols$hierarchical_gene
non_na_indices <- !is.na(symbols$hierarchical_gene)
symbols$unique_gene[non_na_indices] <- make.unique(symbols$hierarchical_gene[non_na_indices], sep = "-copy")



# If there are any NAs in the "unique" column, replace these values with the gene ID in the "mikado_id" column.


symbols$unique_gene[is.na(symbols$unique_gene)] <- symbols$mikado_id[is.na(symbols$unique_gene)]

# Once we have these unique values, we can assign them as names to the GFF file. For simplicity, we're going to assign them to the "Name" slot for both RNA and gene features. In e.g. a RefSeq GFF file, this would be similar to the "gene" slot which is the same for mRNA and gene features, whereas the actual name of the feature changes depending on its "Type". Right now, our focus is on simplicity and interpretability which is why we're going for the first option, which should be effective for most downstream purposes.

# You may also wish to store the other gene symbols in the GFF file under labels like "TOGA_gene", etc. To accomplish this all in one step, we can start by combining the gene symbols table, `symbols`, with the GFF dataframe. We will want to add the gene symbols wherever "mikado_id" matches the "ID" (in the case of gene features) or "Parent" (in the case of RNA features) columns.

# To do this, let's create a column in the GFF also called "mikado_id". Whenever a feature is a gene, we can add the value from the "ID" column, and whenever a feature is an RNA (mRNA or lncRNA) we can add the value from the "Parent" column.


gff$mikado_id <- NA
gff$mikado_id[gff$type %in% c("gene","lncRNA_gene")] <- as.character(unlist(gff$ID[gff$type %in% c("gene","lncRNA_gene")]))
gff$mikado_id[gff$type %in% c("mRNA","lncRNA")] <- as.character(unlist(gff$Parent[gff$type %in% c("mRNA","lncRNA")]))


# Now let's use `dplyr`'s `left_join` function to add the gene symbols data frame to the GFF file.

gff <- dplyr::left_join(gff, symbols,
                        by = "mikado_id")


# We no longer need the "mikado_id" column, so we can set that to `NULL` to get rid of it, and then we can also change the "Name" column to the symbols that are stored in the "unique_gene" column, and then nullify that column, as well. Features don't have to have a name, so short ncRNA genes and RNA features will use their ID as a name.


gff$mikado_id <- NULL
gff$Name <- gff$unique_gene
gff$unique_gene <- NULL


# Before we export the GFF, we need to do some finicky formatting things to make sure the GFF looks right. All columns need to be character vectors (or else blanks will be exported with `=character(0)` after the attribute name in the metadata). This can be done by using `lapply` on the GFF, applying the `as.character` function, and will turn all blanks into `character(0)`, which can then be replaced by `NA`.


gff[] <- lapply(gff, as.character) 
gff[gff == "character(0)"] <- NA


# We can now export the GFF file and the updated gene symbol table.


rtracklayer::export.gff3(gff, "full_annotation.geneSymbols.gff", format = "gff3")

write.table(symbols, file = "gene_symbols_full.tsv",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
