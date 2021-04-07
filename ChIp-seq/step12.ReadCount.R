#' @description  从bam文件中提取一致性peak的read count矩阵 {\BioninfoR<FastCountBam>}
#' @author shijian

library(BioinforR)
merge.bed <- read.table(file="/pub6/Temp/sj/GSE77737/Chip-seq/mergeBed/merge.bed",sep="\t",header=F,fill=T,quote=NULL,stringsAsFactors=F)
merge.bed[,4] <- rep("*",length(merge.bed[,1]))
merge.bed[,5] <- paste0("peak",1:length(merge.bed[,1]))
rownames(merge.bed) <- merge.bed[,5]
bamFilepaths <- list.files("/pub6/Temp/sj/GSE77737/Chip-seq/Bowtie2",pattern = ".sort.bam$",full.names = T,recursive = T)
bamNames <- list.files("/pub6/Temp/sj/GSE77737/Chip-seq/Bowtie2",pattern = ".sort.bam$",full.names = F,recursive = T)
library(stringr)
bamNames <- as.character(str_split(bamNames,"\\.",simplify=T)[,1])
mergeBedReadCount <- FastCountBam(bamFilepaths=bamFilepaths,
                                 ranges = merge.bed,
                                 bamNames = bamNames,
                                 mode="Union",
                                 singleEnd=T,
                                 ignore.strand=T)
save(mergeBedReadCount,file="/pub6/Temp/sj/GSE77737/Chip-seq/mergeBed/mergeBedReadCount.rda")
