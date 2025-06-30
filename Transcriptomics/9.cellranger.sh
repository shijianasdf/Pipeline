#' cellranger gex pipeline.sh
#' @author shijian
#' @param outDir 输出目录，在这个目录霞下会继续生成每个样本的具体目录
#' @param outfilepath 输出命令行文件
#' @param inputDir fastq文件夹地址
#' @param pattern 文件格式 默认提取.fastq.gz格式文件
#' @param creatBam 是否创建bam文件
#' @param ref 参考基因组地址
#' cd ~/project/scMutDB/exp_matrix/SRP239174
#' mkdir SRR10809462
#' cd SRR10809462
#' cellranger count --id=mRNA  \
#' --fastqs=/home/shijian/project/scMutDB/fastq/SRP239174 \
#' --sample=SRR10809462  \
#' --create-bam=true  \
#' --transcriptome=/home/shijian/refData/refdata-gex-GRCh38-2024-A

cellRanger <- function(outDir,outfilepath,inputDir,
                       pattern=".fastq.gz$",createBam=c("true","false")[1],
                       ref="/home/shijian/refData/refdata-gex-GRCh38-2024-A"){
  #创建输出目录
  # if(!file.exists(outDir)){
  #   dir.create(outdir)
  # }
  
  library(stringr)
  
  if(!is.null(inputDir)){
    # outDir <- "D:"
    # inputDir <- "G:/肾透明细胞癌 尹志豪/"
    # pattern <- ".fastq.gz$" 
    t.files <- list.files(inputDir, pattern = pattern, full.names = TRUE, recursive = TRUE)
    sampleNames <- list.files(inputDir, pattern = pattern, full.names = F, recursive = TRUE)
    sampleNames <- gsub(".*/", "", sampleNames)
    sampleNames <- unique(stringr::str_split(sampleNames,"_",simplify = T)[,1])
    
    #创建每个样本具体的输出目录
    if(!file.exists(outDir)){
      stop("please input outDir")
    }
    
    outDir <- file.path(outDir,sampleNames)
    prefix_command <- "cellranger count --id=mRNA"
    commands <- c()
    for(i in 1:length(outDir)){
      if(!file.exists(outDir[i])){
        dir.create(outDir[i])
      }
      t.command <- paste("cd",outDir[i],sep=" ")
      commands <- c(commands,t.command)
      cellranger.command <- paste0(prefix_command, " --fastqs=", inputDir, " --sample=", 
                                   sampleNames[i], " --create-bam=",createBam, " --transcriptome=", ref)
      commands <- c(commands,cellranger.command)
    }
    
  }else{
    stop("please input inputDir!")
  }
  
  writeLines(commands,con=outfilepath)
  return(commands)

}

cellRanger(outDir="/data/shijian/project/scMutDB/exp_matrix/SRP239174",
           outfilepath="/data/shijian/project/NGSCommand/cellranger_gex.sh",
           inputDir="/data/shijian/project/scMutDB/fastq/SRP239174",
           pattern=".fastq.gz$",createBam=c("true","false")[1],
           ref="/home/shijian/refData/refdata-gex-GRCh38-2024-A")
#nohup bash /data/shijian/project/NGSCommand/cellranger_gex.sh > cellranger.log 2>&1 &   



