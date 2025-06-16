library(optparse)
library(rtracklayer)

option_list <- list(
  make_option(c("-g", "--gff"), type = "character", default = "FALSE",
              help = "The gff file for the reference assembly after initial cleaning",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

gffname <- opt[[1]]

# $prefix.full.gffread.gff

prefix <- gsub(".full.gffread.gff","",gffname)

# Read in GFF

gff <- as.data.frame(rtracklayer::readGFF(gffname))

# Isolate the ID, Name, and type columns. This will create a key that we will use later to keep track of which names match with which IDs

key <- gff[,c("ID","Name","type")]

# Only keep gene and transcript features. We can do this by using `grepl` to return a vector of TRUE or FALSE values and using this to isolate the rows in the key. The following `grepl` statement will match any feature that contains "gene", "transcript", or "RNA" (so will match "pseudogene", "ncRNA", etc.).

key <- key[grepl("gene|transcript|RNA", key$type),]

# Since genes and pseudogenes share the same gene name but we want gene IDs to be unique, we can add "pseudo" after the gene name if something is a pseudogene.

key$Name[key$type == "pseudogene"] <- gsub("$","-pseudo",
                                           key$Name[key$type == "pseudogene"])

# Remove anything that has NA in the "Name" column

key <- key[!is.na(key$Name),]

# If a feature is not a gene, add "rna-" in front of the name

key$Name[grepl("transcript|RNA", key$type)] <- gsub("^","rna-",
                                                    key$Name[grepl("transcript|RNA", key$type)])

# For any features that still don't have a unique gene name, add a copy number to the end

key$Name <- make.unique(key$Name, sep = "-copy")

# Now we want to replace the IDs and Parent features with the new unique names. Whenever there is an exact match between something in the ID or Parent column in the GFF with the ID in the key, the new name will replace the old ID or Parent feature.

# First, create a named vector for fast lookup
id_to_name <- setNames(key$Name, key$ID)
# Replace values in gff$ID if they match any key$ID
gff$ID <- ifelse(gff$ID %in% key$ID, id_to_name[gff$ID], gff$ID)
# Unlist Parent
gff$Parent <- sapply(gff$Parent, function(x) {
  if (length(x) == 0) NA_character_ else as.character(x)
})
# Replace values in gff$Parent if they match any key$ID
gff$Parent <- ifelse(gff$Parent %in% key$ID, id_to_name[gff$Parent], gff$Parent)

rtracklayer::export.gff3(gff,
                         paste0(prefix,".gffread.F.nameIDs.gff3"),
                         format = "gff3")
