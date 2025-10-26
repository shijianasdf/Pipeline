#'作为一款比对软件，建index肯定是必不可少的一步
#'--genomeDir 索引文件夹
#'--genomeFastaFiles 参考基因组fasta序列
#'--sjdbGTFfile 基因注释信息gtf文件
#' @USAGE STAR --runThreadN 20 \
#'--runMode genomeGenerate \
#'--genomeDir /Users/shijian/mydata/bulkRNA/0.reference/star.index/ \
#'--genomeFastaFiles /Users/shijian/mydata/bulkRNA/0.reference/GRCh38.p13.genome.fa \
#'--sjdbGTFfile /Users/shijian/mydata/bulkRNA/0.reference/gencode.v38.chr_patch_hapl_scaff.annotation.gtf \
#'--sjdbOverhang 100

#' @description star pipeline
#' @author shijian
#' @param sampleInfo 样本注释信息 主要有三列(样本名,fastq1地址，fastq2地址) Patient1 /Users/shijian/mydata/bulkRNA/1.fastq/P4_1.fq.gz /Users/shijian/mydata/bulkRNA/1.fastq/P4_2.fq.gz
#' @param sample 样本名信息列
#' @param fastq1 fastq1地址信息列
#' @param fastq2 fastq2地址信息列
#' @param index star index目录
#' @param outfilepath 输出命令行.sh地址
#' @param outDir 输出结果目录
#' @return 一个列表 包含样本注释信息和命令行；同时生成.sh文件用于批处理；生成star比对结果文件夹
#' @usage STAR --runThreadN 10 \ 线程数
#'  --genomeDir /Users/shijian/mydata/bulkRNA/0.reference/star.index/ \ 参考基因组索引目录
#'  --twopassMode Basic \
#'  --readFilesIn /Users/shijian/mydata/bulkRNA/1.fastq/P1_1.fq.gz /Users/shijian/mydata/bulkRNA/1.fastq/P1_2.fq.gz \ 双末端fastq
#'  --readFilesCommand zcat \ 如果文件是fq.gz 需要指定zcat 注意如果是macs需要设置为gzcat
#'  --outSAMtype BAM SortedByCoordinate \ 返回sort后的bam文件
#'  --outFileNamePrefix /Users/shijian/mydata/bulkRNA/3.bam/P132/P132 指定输出目录及前缀
runStar <- function(sampleInfo,
                    sample,
                    fastq1,
                    fastq2,
                    index,
                    outfilepath,
                    outDir){
  library(dplyr)
  #列名重命名
  sampleInfo <- sampleInfo[,c(sample,fastq1,fastq2)] 
  colnames(sampleInfo) <- c("sample","fastq1","fastq2")
  commands <- c()
  command <- paste0("STAR --runThreadN 10 --genomeDir ", index ," --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts --readFilesIn ")
  for(i in 1:length(sampleInfo$sample)){
    dir <- file.path(outDir,sampleInfo$sample[i])
    if(!file.exists(dir))
      dir.create(dir,recursive = T)
    tempcommand <- paste0(command,sampleInfo$fastq1[i]," ",sampleInfo$fastq2[i],
                          " --readFilesCommand gzcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ",
                          paste0(dir,"/",sampleInfo$sample[i]))
    commands <- c(commands,tempcommand)
  }
  writeLines(commands,con = outfilepath)  
  return(list(commands=commands,sampleInfo=sampleInfo))
}
starRes <- runStar(sampleInfo=sampleInfo,
                    sample="sample",
                    fastq1="fastq1",
                    fastq2="fastq2",
                    index="/Users/shijian/mydata/bulkRNA/0.reference/star.index/",
                    outfilepath="/Users/shijian/mydata/bulkRNA/star.sh",
                    outDir="/Users/shijian/mydata/bulkRNA/3.bam")
#bash /Users/shijian/mydata/bulkRNA/star.sh > log.txt &


#' 升级了，对上面那版有改进，更简单
runStar <- function(fastqDir,
                    pattern = "fq.gz",
                    index,
                    outfilepath,
                    outDir){
  library(dplyr)
  fastqfiles <- list.files(fastqDir,full.names = T,recursive = T,pattern = pattern)
  fastqnames <- list.files(fastqDir,recursive = T,pattern = pattern)
  fastq1 <- fastqfiles[seq(1,length(fastqfiles),by=2)] #奇数
  fastq2 <- fastqfiles[seq(2,length(fastqfiles),by=2)] #偶数
  sample <- stringr::str_split(fastqnames[seq(1,length(fastqfiles),by=2)],"_",simplify = T)[,1]
  #生成sampleInfo
  sampleInfo <- data.frame(sample=sample,fastq1=fastq1,fastq2=fastq2)
  #列名重命名
  #sampleInfo <- sampleInfo[,c(sample,fastq1,fastq2)] 
  #colnames(sampleInfo) <- c("sample","fastq1","fastq2")
  commands <- c()
  command <- paste0("STAR --runThreadN 10 --genomeDir ", index ," --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts --readFilesIn ")
  for(i in 1:length(sampleInfo$sample)){
    dir <- file.path(outDir,sampleInfo$sample[i])
    if(!file.exists(dir))
      dir.create(dir,recursive = T)
    # tempcommand <- paste0(command,sampleInfo$fastq1[i]," ",sampleInfo$fastq2[i],
    #                       " --readFilesCommand gzcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ",  
    #                       paste0(dir,"/",sampleInfo$sample[i]))  macos 是gzcat
    tempcommand <- paste0(command,sampleInfo$fastq1[i]," ",sampleInfo$fastq2[i],
                          " --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ",  
                          paste0(dir,"/",sampleInfo$sample[i])) #linux是zcat
    
    commands <- c(commands,tempcommand)
  }
  writeLines(commands,con = outfilepath)  
  return(list(commands=commands,sampleInfo=sampleInfo))
}
starRes <- runStar(fastqDir = "/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/3.trimGalore/trimGalore_result",
                   index="/data/shijian/refData/mouse_reference/mouse_index/star.index",
                   outfilepath="/data/shijian/project/NGSCommand/star.sh",
                   outDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/5.star_bam")

#'star比对后没有bam文件的索引,需要利用samtools先建立索引,因为后续htseq-count会利用索引
#'循环执行samtools index SRR11050949Aligned.sortedByCoord.out.bam >
runSamtools <- function(bamDir,
                        outfilepath){
  bamfiles <- list.files(bamDir,full.names = T,recursive = T,pattern = ".bam$")
  bamnames <- list.files(bamDir,recursive = T,pattern = ".bam$")
  out.dir <- dirname(bamfiles)
  command <- "samtools index "
  commands <- c()
  for(i in 1:length(bamfiles)){
    tempcommand <- paste0(command,bamfiles[i])
    commands <- c(commands,tempcommand)
  }
  writeLines(commands,con = outfilepath) 
  return(commands)
}
runSamtools("/Users/shijian/mydata/bulkRNA/3.bam",
            outfilepath="/Users/shijian/mydata/bulkRNA/samtools.sh")
#bash /Users/shijian/mydata/bulkRNA/samtools.sh > log_samtools.txt &
