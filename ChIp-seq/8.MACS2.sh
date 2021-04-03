#' @description  MACS2 call peaks
#' @author shi jian
#' 
#' macs2 callpeak -t SRR1042593.sorted.bam -c SRR1042594.sorted.bam  -f BAM -B -g hs -n Xu_MUT_rep1 
#' parameters -t -c -f:c("AUTO","BAM","SAM","BED",...) -g:c("hs","mm","ce","dm"),-B
#'            --outdir:输出结果目录  -n:项目名字 -q -p --broad:宽peak --broad-cutoff:
#' 
#' @param  bamDir bowtie比对，samtools sort之后的sort.bam目录地址
#' @param  SampleInfo SRR样本注释大表地址，"D:/Rsources/Project/StudySingleCell/data/SampleInfo.rda"
#'         需要3列：SRSid sample_accession，表观修饰类型(chromState)，以及细胞系类型(cellLine)
#' @param  outfilepath .sh输出路径
#' @param  outDir MACS2工作目录
#' @param  fileType 文件类型
#' @param  isBroad 是否识别宽peaks
#' @param  broad.cutoff 识别宽peak的阈值
#' @param  genome 基因组大小
FastMACS2 <- function(bamDir,
                  SampleInfo,
                  outfilepath,
                  outDir,
                  fileType=c("AUTO","BAM","SAM","BED")[2],
                  isBroad=FALSE, 
                  broad.cutoff = 0.1,
                  genome=c("hs","mm","ce","dm")[1]){
  #macs2 callpeak -t SRR1042593.sorted.bam -c SRR1042594.sorted.bam  -f BAM -B -g hs -n Xu_MUT_rep1 
  #创建MACS2工作目录
  if(!file.exits(outDir)){
    dir.create(outDir)
  }
  #导入SRR注释大表,需要SRRid，染色质信号种类（H3K27AC等），细胞系
  SampleInfo <- get(load(SampleInfo))
  # SampleInfo$sample_accession
  # SampleInfo$cellLine
  # SampleInfo$chromState
  
  bamfiles <- list.files(bamDir,pattern = ".sort.bam$",full.names = T)
  sampleNames <- list.files(bamDir,pattern = ".sort.bam$",full.names = F)
  sampleNames <- gsub(".sort.bam$","",sampleNames)
  table(SampleInfo$chromState)
  
  macs2CMD_pre <- "macs2 callpeak"
  
  # 根据注释文件的"sample_accession","cellLine","chromState"以及待分析样本的sampleNames
  # 找到MACS 所需的实验组和对照组样本，并做成表格形式
  exc_SampleInfo <- SampleInfo$sample_accession[SampleInfo$sample_accession %in% sampleNames,]
  control <- exc_SampleInfo[grep("input",exc_SampleInfo$chromState),c("sample_accession","cellLine","chromState")]
  colnames(control) <- c("contrl","cellLine","control.chromState")
  treat <- exc_SampleInfo[!grepl("input",exc_SampleInfo$chromState),c("sample_accession","cellLine","chromState")]
  colnames(treat) <- c("treat","cellLine","treat.chromState")
  temp <- merge(treat,control,by="cellLine",all.x=T)
  temp$treatpath <- file.path(bamDir,paste0(temp$treat,".sort.bam"))
  temp$controlpath <- file.path(bamDir,paste0(temp$control,".sort.bam"))
  commands <- c()
  for(i in 1:dim(temp)[1]){
    if(!is.na(temp$control[i])){
      macs2CMD_pre <- paste(macs2CMD_pre,"-t",temp$treatpath[i],"-c",temp$controlpath[i]),"-f",fileType,"-B -g",genome,"-n",paste(temp$cellLine,temp$treat.chromState[i],sep="."))
    }else{
      macs2CMD_pre <- paste(macs2CMD_pre,"-t",temp$treatpath[i],"-f",fileType,"-B -g",genome,"-n",paste(temp$cellLine,temp$treat.chromState[i],sep="."))
    }
    if(isBroad == T){
      macs2CMD_pre <- paste(macs2CMD_pre,"--broad","--broad-cutoff",broad.cutoff)
    }
    macs2CMD <-  paste(macs2CMD_pre,"--outdir",temp$treat[i])
    print(macs2CMD)
    commands <- c(commands,macs2CMD)
    system(macs2CMD)
  }
  writeLines(commands,con=outfilepath)
  return(commands)
}
