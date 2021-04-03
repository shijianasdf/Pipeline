#' @description SRA ENA 大表读入整理,个别情况下根据实际情况自定义字段名字
#' @author  shi jian
#' SRA大表是以逗号分隔的文件，而ENA大表是标准的数据框，我们用R对两个文件进行读入和整理
#' 提取整合成我们需要的字段信息
#' 该程序不可能无限可重复，需要观察提取自定义字段，所以无法封装成函数，但是一些代码有借鉴意义

#读入ENA大表
ENA <- read.table("C:/Users/SJ/Downloads/filereport_read_run_PRJNA311329_tsv.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
library(stringr)
tt <- str_split(ENA$sample_title," ",simplify = T)
tt <- cbind.data.frame(tt[,1],paste(tt[,2],tt[,3],tt[,4]))
tt[,1] <- str_split(tt[,1],"_",simplify = T)[,1]
colnames(tt) <- c("cellLine","chromState")
tt$cellLine

#读入SRA大表
SRAl <- readLines(con = "C:/Users/SJ/Downloads/SraRunTable.txt")
library(stringr)
#逐行读入数据，每行字符串以逗号分割，但是分号内的逗号不分割
colns <- str_split(SRAl[1],pattern=",(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)",simplify = T)[1,]
sra <- str_split(SRAl,pattern=",(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)",simplify = T)[-1,]
colnames(sra) <- colns 
sra <- as.data.frame(sra)
sra$Cell_Line
head(sra)
SRA <- sra
head(SRA)

#以SRRid为变量，融合两个大表
identical(SRA[match(ENA$run_accession,SRA$Run),]$Run,ENA$run_accession)
SRA <- SRA[match(ENA$run_accession,SRA$Run),]
ENA <- cbind.data.frame(ENA,SRA$Cell_Line,SRA$Tissue,tt)
SampleInfo <- ENA
head(SampleInfo)
colnames(SampleInfo)[3:4] <- c("secondary_sample_accession","sample_accession")
save(SampleInfo,file="D:/Rsources/Project/StudySingleCell/data/SampleInfo.rda")
write.table(SampleInfo,file="D:/Rsources/Project/StudySingleCell/data/SampleInfo.txt",sep="\t",quote=F,row.names=F,col.names=T)

