library(rtracklayer)

options(scipen=99999999)

Repeats <- read.table(list.files(pattern = "filteredRepeats.bed"),
                      header=F,as.is=T,sep="\t")


Repeat2 <- Repeats[grep("RNA",Repeats$V4),]

write.table(Repeat2[,1:3], file = "ncRNA.fromRepeat.bed",quote=F,
            row.names = FALSE,col.names=FALSE,sep="\t")

seeding <- read.table("assembly.rfam.bed",header=F,as.is=T,sep="\t")

write.table(seeding[,1:3],file="assembly.rfam.s.bed",quote=F,row.names=F,col.names=F,sep="\t")

mikado_lenient <- readGFF("mikado_lenient.gff")

mikado_ncRNA <- mikado_lenient[grep("ncRNA",mikado_lenient$type),]


write.table(mikado_ncRNA[,c(1,4,5)], file = "mikado_ncRNA.bed",quote=F,
            row.names=F,col.names=F,sep="\t")

liftOff_asm <- readGFF(list.files(pattern = "liftoff.gtf"))
biotype <- colnames(liftOff_asm)[grep("biotype",colnames(liftOff_asm))][1]
lifted_RNA <- liftOff_asm[grep("RNA",liftOff_asm[,biotype]),]


write.table(lifted_RNA[,c(1,4,5)], file = "liftOff_ncRNA.bed",quote=F,
            row.names=F,col.names=F,sep="\t")


