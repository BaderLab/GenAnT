library(optparse)

option_list <- list(
  make_option(c("-c", "--customRef"), type = "character", default = "FALSE",
              help = "Whether the custom gff file should be considered the reference assembly",
              metavar = "character"),
  make_option(c("-l", "--liftOffRef"), type = "character", default = "FALSE",
              help = "Whether the liftOff gff file should be considered the reference assembly",
              metavar = "character")
)

# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

custom_ref <- as.logical(toupper(opt[[1]]))
liftoff_ref <- as.logical(toupper(opt[[2]]))

gff <- list.files(pattern="*gff")
dir <- getwd()

rows <- c()
for(i in gff) {
  n = gsub(".gffread.gff","",i)
  

  therow <- c(paste0(dir,"/",i),n, "True","0","False","False")
  
  if(grepl("custom",tolower(n)[1])) {
    if(custom_ref) {
      therow <- c(paste0(dir,"/",i),n, "True","0","True","False")
    }
  }
  
  if(grepl("liftoff",tolower(n)[1])) {
    if(liftoff_ref & !custom_ref) {
      therow <- c(paste0(dir,"/",i),n, "True","0","True","False")
    }
  }
  
  rows <- rbind(rows,therow)
}

write.table(rows, file = "mikado_input_sheet.txt",quote=F,row.names=F,
            col.names=F,sep="\t")
