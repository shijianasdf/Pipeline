#' @description bam2bigwig Measure read count in genomic bins{\deeptools<bamCoverage>} 
#' @author shi jian
#' @usage  
#' bamCoverage -bs 25 -p max --ignoreDuplicates --normalizeUsingRPKM -b /pub5/xiaoyun/BigData/CellLines/MDA-MB-231/RNA-seq/GSE38790_GSM949589/2.TOPHAT/accepted_hits.bam -o /pub5/xiaoyun/BigData/CellLines/MDA-MB-231/RNA-seq/GSE38790_GSM949589/4.GenomeBrowser/accepted_hits.bam.bw
#' @param bamDir sort.bam文件地址
#' @param SampleInfo SRR文件注释信息
#' @param outfilepath .sh输出地址
#' @param outDir bw输出文件夹
#' @param bin 开窗大小
#' @param type 标化方式
#' @param threads 进程数
#' @param pattern 文件匹配模式
#' @param extraParameters 额外参数设置
#' @export
Fastbam2bw <- function(bamDir,
                   SampleInfo,
                   outfilepath,
                   outDir,
                   bin,
                   type=c("RPKM","CPM","BPM","RPGC","None"),
                   threads="max",
				   extraParameters=NULL,
                   pattern=".sort.bam$"){
  
    bamFiles <- list.files(bamDir,pattern = pattern,full.names = T,recursive = T) 
    SampleNames <- list.files(bamDir,pattern = pattern,full.names = F,recursive = T) 
    library(stringr)
    SampleNames <- str_split(sampleNames,"\\.",simplify = T)[,1]
    SampleInfo <- get(load(SampleInfo))
    
    if(!file.exists(outDir)){
      dir.create(outDir)
    }
    
    commands <- c()
    for(i in 1:length(bamFiles)){
      pos <- which(SampleInfo$sample_accession == SampleNames[i])
      bw.name <- paste(SampleNames[i],str_replace_all(SampleInfo$CellLine[pos]," ","."),str_replace_all(SampleInfo$chromState[pos]," ","_") ,sep="_")
      runCMD <- paste("bamCoverage -bs",bin,"-p",threads,"--normalizeUsing",type,"-b",bamFiles[i],"-o",file.path(outDir,paste0(bw.name,".bw")))
      if(!is.null(extraParameters)){
	    runCMD <-paste(runCMD,extraParameters) 
	  }
	  commands <- c(commands,runCMD)
    }

    writeLines(commands,con = outfilepath)
    return(commands)
 
}
Fastbam2bw(bamDir="/pub6/Temp/sj/GSE77737/DNase-seq/Bowtie2",
                   SampleInfo="/pub6/Temp/sj/GSE77737/SampleInfo.rda",
                   outfilepath="/pub6/Temp/sj/GSE77737/NGScommands/DNasebam2bw1.sh",
                   outDir="/pub6/Temp/sj/GSE77737/DNase-seq/GenomeBrowser",
                   bin=25,
                   type=c("RPKM","CPM","BPM","RPGC","None")[1],
                   threads="max",
				   extraParameters=NULL,
                   pattern=".sort.bam$")
./DNasebam2bw.sh
