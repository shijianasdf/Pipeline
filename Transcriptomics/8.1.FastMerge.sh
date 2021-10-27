#' @description merge htseq sample count to matrix
#' @author shi jian
#' @param countsDir htseq count 输出目录
#' @return 合并后的表达矩阵
#' 

mergeCounts <- function(countsDir){
  library(dplyr)
  library(biomaRt)
  library(curl)
  countsfile <- list.files(countsDir,full.names = T,recursive = T)
  counts <- list.files(countsDir,recursive = T)
  patients <- str_split(counts,".",simplify = T)[,dim(str_split(dirname(bamfiles),".",simplify = T))[2]]
  options(stringsAsFactors = FALSE)
  expr.matrix <- read.table(countsfile[1],sep = "\t",col.names = c("gene_id",patients[1]))
  if(length(patients) >= 2){
    for(i in 2:length(patients)){
      tempCount <- read.table(countsfile[i],sep = "\t",col.names = c("gene_id",patients[i]))
      expr.matrix <- full_join(expr.matrix,tempCount,by="gene_id")
    }
    expr.matrix<- expr.matrix[-1:-5,]
    ENSEMBL <- gsub("\\.\\d*", "", expr.matrix$gene_id)
    row.names(expr.matrix) <- ENSEMBL
    #convert 转换
    #。。。。。。。
  }
  return(expr.matrix)
}
