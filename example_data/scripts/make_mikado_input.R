gff <- list.files(pattern="*gff")
dir <- getwd()

rows <- c()
for(i in gff) {
  n = gsub(".gffread.gff","",i)
  therow <- c(paste0(dir,"/",i),n, "True","0","False","False")
  rows <- rbind(rows,therow)
}

write.table(rows, file = "mikado_input_sheet.txt",quote=F,row.names=F,
            col.names=F,sep="\t")
