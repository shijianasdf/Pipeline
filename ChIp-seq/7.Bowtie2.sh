#' @description  bowtie2 align
#' @author shi jian
#' 构建索引 bowtie2-build mm10.fa mm10
#' 单末端比对 bowtie2 -p 10 -x genome_index -U input.fq -S output.sam
#' 双末端比对 bowtie2 -p 6  -x mm10 -1 example_1.fastq -2 example_2.fastq -S example.sam
#' sam转bam samtools sort example.sam > example.bam
#' bam进行sort samtools sort my.bam -o my.sort.bam
#' 对sort后的bam进行索引 samtools index my.sort.bam
#'
#' @param FastqDir fastq所在文件夹
#' @param index 参考基因组索引 "/pub6/Temp/sj/GSE77737/GRCh38_noalt_as/GRCh38_noalt_as",
#'               注意我的索引文件放在/pub6/Temp/sj/GSE77737/GRCh38_noalt_as目录下，
#'               但是需要输入的参数为/pub6/Temp/sj/GSE77737/GRCh38_noalt_as/GRCh38_noalt_as。最后一个GRCh38_noalt_as指的是共用文件名
#'               下载地址: https://benlangmead.github.io/aws-indexes/bowtie,也可以自己生成
#' @param outDir 比对输出目录
#' @param pairEND 单双末端
#' @param is.trimGalore 是否经过trimGalore处理
#' @param isBAm 是否生成BAM文件
#' @param threads 进程数

#'bowtie2(FastqDir="/pub6/Temp/sj/GSE77737/DNase-seq/trimGalore_result",
#'        index="/pub6/Temp/sj/GSE77737/GRCh38_noalt_as/GRCh38_noalt_as",
#'        outDir="/pub6/Temp/sj/GSE77737/DNase-seq/Bowtie2",
#'        outfilepath = "/pub6/Temp/sj/GSE77737/NGScommands/bowtie2.sh",
#'        pairEND = F,
#'        is.trimGalore = T,
#'        isBAM=TRUE,
#'        threads=24)
Fastbowtie2 <- function(FastqDir,
                    index,
                    outDir,
                    outfilepath,
                    pairEND = T,
                    is.trimGalore = T,
                    isBAM=TRUE,
                    extraParameters = NULL,
                    threads=24
){
  #生成输出目录
  if(!file.exists(outDir)){
    dir.create(outDir)
  }
  
  commands <- c()
  if(pairEND ==T){
    if(is.trimGalore == T){
      fastqFiles1 <- list.files(FastqDir,pattern="_1_val_1.fq.gz$",full.names=T,recursive = T)
      fastqFiles2 <- list.files(FastqDir,pattern="_2_val_2.fq.gz$",full.names=T,recursive = T)
      sampleNames <- list.files(FastqDir, pattern = "_1_val_1.fq.gz$", full.names = F)
      sampleNames <- gsub("_1_val_1.fq.gz$", "", sampleNames)
    }else{
      fastqFiles1 <- list.files(FastqDir,pattern="_1.fastq.gz$",full.names=T,recursive = T)
      fastqFiles2 <- list.files(FastqDir,pattern="_2.fastq.gz$",full.names=T,recursive = T)
      sampleNames <- list.files(FastqDir, pattern = "_1.fastq.gz$", full.names = F)
      sampleNames <- gsub("_1.fastq.gz$", "", sampleNames)
    }
    #可以加bowtie2的参数
    runAlign_prefix <- paste("bowtie2 -p",threads,"-x",index)
    
    if( !is.null(extraParameters) ){
      runAlign_prefix <- paste(runAlign_prefix,extraParameters)
    }
    
    for(i in 1:length(fastqFiles1)){
      runAlign <- paste(runAlign_prefix,"-1",fastqFiles1[i],"-2",fastqFiles2[i],"-S", file.path(outDir,paste0(sampleNames[i],".sam")))
      print(runAlign)
      commands <- c(commands,runAlign)
      system(runAlign)
    }
    
  }else{
    if(is.trimGalore == T){
      fastqFiles <- list.files(FastqDir,pattern="_trimmed.fq.gz$",full.names=T,recursive = T)
      sampleNames <- list.files(FastqDir,pattern="_trimmed.fq.gz$",full.names=F,recursive = T)
      sampleNames <- gsub("_trimmed.fq.gz$","",sampleNames)
    }else{
      fastqFiles <- list.files(FastqDir,pattern=".fastq.gz$",full.names=T,recursive = T)
      sampleNames <- list.files(FastqDir,pattern=".fastq.gz$",full.names=F,recursive = T)
      sampleNames <- gsub(".fastq.gz$","",sampleNames)
    }
    
    #可以加bowtie2的参数
    runAlign_prefix <- paste("bowtie2 -p",threads,"-x",index)
    
    if( !is.null(extraParameters) ){
      runAlign_prefix <- paste(runAlign_prefix,extraParameters)
    }
    
    for(i in 1:length(fastqFiles)){
      runAlign <- paste(runAlign_prefix,"-U",fastqFiles[i],"-S", file.path(outDir,paste0(sampleNames[i],".sam")))
      print(runAlign)
      commands <- c(commands,runAlign)
      system(runAlign)
    }
    
  }
  
  if(isBAM == T){
    samFiles <- list.files(outDir,pattern = ".sam$",full.names = T,recursive = T)
    a <- sapply(samFiles,sam2bam)
	a <- as.vector(a)
    commands <- c(commands,a)
  }
  writeLines(commands,con=outfilepath)
  return(commands)
}


#' @description  sam to bam and build index for bam
#' @param samFile .sam文件路径+文件名
#' @author shi jian
#' samtools view -S -b -o /pub6/Temp/sj/GSE77737/DNase-seq/Bowtie2/SRS1282565.bam /pub6/Temp/sj/GSE77737/DNase-seq/Bowtie2/SRS1282565.sam
#' samtools sort /pub6/Temp/sj/GSE77737/DNase-seq/Bowtie2/SRS1282565.bam -o /pub6/Temp/sj/GSE77737/DNase-seq/Bowtie2/SRS1282565.sort.bam
#' samtools index /pub6/Temp/sj/GSE77737/DNase-seq/Bowtie2/SRS1282565.sort.bam
sam2bam <- function(samFile){
  print("converting sam to bam ...")
  # samtools view -S -b -o my.bam my.sam
  bamFile <- gsub("\\.sam$", "\\.bam", samFile)
  convCmd <- paste("samtools view -S -b -o", bamFile, samFile)
  system(convCmd)
  # To create index
  print("sorting first ...")
  # samtools sort my.bam -o my.sort.bed
  sortBamFile <- file.path(dirname(bamFile),paste0(unlist(strsplit(basename(bamFile),"\\."))[1],".sort.bam") )
  sortCmd <- paste("samtools sort", bamFile, "-o", sortBamFile)
  system(sortCmd)
  print("then creating index ...")
  # samtools index my.sorted.bam
  createCmd <- paste("samtools index", sortBamFile)
  system(createCmd)
  commands <- c(convCmd,sortCmd,createCmd)
  return(commands)
}
