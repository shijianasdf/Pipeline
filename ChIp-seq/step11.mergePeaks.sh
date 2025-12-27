#' @description merge peaks {\bedtools<merge>} 每个样本的peak进行合并，得到一致性peak
#' @author  shijian
#' @usage 
#' sortBed -i A.bed > sort.bed
#' bedtools merge -i sort.bed > merge.sort.bed
#' @param  bedDir bed文件所在目录
#' @param  outDir 输出目录
#' @param  outname 输出merge后bed的文件名
#' @param  outfilepath .sh输出路径
#' @param  pattern 文件匹配模式

#' @export
mergePeaks <- function(bedDir,
                       outDir,
                       outname,
                       outfilepath,
                       pattern=".narrowPeak$"){
  #创建工作目录
  if(!file.exists(outDir)){
    dir.create(outDir)  
  }
  
  commands <- c()
  #将所有peak复制到工作目录
  bedFiles <- list.files(bedDir,pattern = pattern,full.names = T,recursive = T) 
  for(i in 1:length(bedFiles)){
    cp.cmd <- paste("cp",bedFiles[i],outDir)
    commands <- c(commands,cp.cmd)
    #system( cp.cmd )
  }
  
  #通过管道操作cat *.bed成一个bed文件并进行sort和merge
  runCMD <- paste0( "cat ",outDir,"/*.narrowPeak |", " sort -k1,1 -k2,2n | bedtools merge > ",file.path(outDir,outname) )
  commands <- c(commands,runCMD)
  print(runCMD)
  #system(runCMD)
  writeLines(commands,con=outfilepath)
  return(commands)
}
mergePeaks(bedDir="/pub6/Temp/sj/GSE77737/Chip-seq/MACS2",
           outDir="/pub6/Temp/sj/GSE77737/Chip-seq/mergeBed",
           outname="merge.bed",
           outfilepath = "/pub6/Temp/sj/GSE77737/NGScommands/Chip_mergePeaks.sh",
           pattern=".narrowPeak$")
./Chip_mergePeaks.sh 2>> mergePeaks.log
