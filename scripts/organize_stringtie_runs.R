RNA_aln <- list.files(pattern = "*.bam",path = "RNAseq_alignment")
iso_aln <- list.files(pattern = "*.bam",path = "ISOseq_alignment")

inter <- intersect(RNA_aln,iso_aln)


RNAonly <- RNA_aln[!(RNA_aln %in% inter)]
ISOonly <- iso_aln[!(iso_aln %in% inter)]

write.table(inter, file = "stringtie_out/mixed.txt",quote=F,row.names = F,col.names=F,sep="\t")

write.table(RNAonly, file = "stringtie_out/rnaseq_only.txt",quote=F,row.names = F,col.names=F,sep="\t")

write.table(ISOonly, file = "stringtie_out/isoseq_only.txt",quote=F,row.names = F,col.names=F,sep="\t")
