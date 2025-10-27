#RSeQC 判断RNA-seq建库类型
#工作原理
#该脚本通过统计那些唯一映射到某个基因的某个外显子上的读段，检查这些读段是如何映射的。它特别关注那些映射到有链注释的基因（即明确知道该基因在正链还是反链上）的读段。
#它会计算两类读段的比例：
#“1++，1--，2+-，2-+”：读段的方向与基因所在链的方向一致（对应 strand=yes）。
#“1+-，1-+，2++，2--”：读段的方向与基因所在链的方向相反（对应 strand=reverse）。
infer_experiment.py -r <bed_file> -i <bam_file>
#-r 或 --refgene：参考基因模型的BED文件。这是一个关键文件，需要提前准备。
#-i 或 --input-file：你的 RNA-seq 比对后的 BAM 文件。
# 假设你有一个gtf to bed文件 <bedops>
gtf2bed < /data/shijian/refData/mouse_reference/gencode.vM37.chr_patch_hapl_scaff.annotation.gtf > /data/shijian/refData/mouse_reference/gencode.vM37.chr_patch_hapl_scaff.annotation.bed
infer_experiment.py -r /data/shijian/refData/mouse_reference/gencode.vM37.chr_patch_hapl_scaff.annotation.bed -i /data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/5.star_bam/1/1Aligned.sortedByCoord.out.bam


#' @description  HTSEQ pipeline
#' @author  shijian
#' @param bamDir bam文件夹
#' @param pattern "sortedByCoord.out.bam$"
#' @param strand=c("yes","no","reverse")[1] 链特异性信息
#' @param gtf.path gtf注释文件路径
#' @param out.dir 输出结果路径
#' @param outfilepath 命令行文件路径
#' @usage htseq-count -f bam \ default: sam 指定bam格式，该参数的值可以是sam或bam。
#' -r name \  default: name 设置sam或bam文件的排序方式，该参数的值可以是name或pos。
#' -s yes \ default: yes 设置是否是链特异性测序。该参数的值可以是yes,no或reverse。
#' -a 10 \ default: 10 忽略比对质量低于此值的比对结果。在0.5.4版本以前该参数默认值是0。
#' -t exon \  default: exon 程序会对该指定的feature（gtf/gff文件第三列）进行表达量计算，而gtf/gff文件中其它的feature都会被忽略。
#' -i gene_id \ default: gene_id 设置feature ID是由gtf/gff文件第9列那个标签决定的；若gtf/gff文件多行具有相同的feature ID，则它们来自同一个feature，程序会计算这些features的表达量之和赋给相应的feature ID。
#' -m intersection-nonempty \ default: union 设置表达量计算模式。该参数的值可以有union, intersection-strict and intersection-nonempty。这三种模式的选择请见上面对这3种模式的示意图。从图中可知，对于原核生物，推荐使用intersection-strict模式；对于真核生物，推荐使用union模式。
#' yourfile_name.bam /Users/shijian/mydata/bulkRNA/0.reference/gencode.v38.chr_patch_hapl_scaff.annotation.gtf > counts.txt
runHTSEQ <- function(bamDir,
                     pattern = "sortedByCoord.out.bam$",
                     strand=c("yes","no","reverse")[1],
                     gtf.path,
                     out.dir,
                     outfilepath){
  library(stringr)
  bamfiles <- list.files(bamDir,full.names = T,recursive = T,pattern = pattern)
  #bamfiles <- list.files(bamDir,full.names = T,recursive = T,pattern = ".bam$")
  patients <- str_split(dirname(bamfiles),"/",simplify = T)[,dim(str_split(dirname(bamfiles),"/",simplify = T))[2]]
  if(!file.exists(out.dir))
    dir.create(out.dir)
  command <- paste0("htseq-count -f bam -r name -s ",strand," -m intersection-nonempty ")
  commands <- c()
  for(i in 1:length(bamfiles)){
    tempcommand <- paste0(command,bamfiles[i]," ",gtf.path," > ",file.path(out.dir,paste0(patients[i],".txt")))
    commands <- c(commands,tempcommand)
  }
  writeLines(commands,con = outfilepath)
  return(commands)
}
# htseq_res <- runHTSEQ(bamDir="/Users/shijian/mydata/bulkRNA/3.bam",
#                     #sampleInfo,
#                     strand=c("yes","no","reverse")[1],
#                     gtf.path="/Users/shijian/mydata/bulkRNA/0.reference/gencode.v38.chr_patch_hapl_scaff.annotation.gtf",
#                     out.dir="/Users/shijian/mydata/bulkRNA/4.counts",
#                     outfilepath="/Users/shijian/mydata/bulkRNA/htseq.sh")
# list.files("/Users/shijian/mydata/bulkRNA/3.bam",full.names = T,recursive = T,pattern = ".bam$")
#bash /Users/shijian/mydata/bulkRNA/htseq.sh > log.txt &

#' @description
#' htseq-count 聚合每个样本的表达成为最终的表达矩阵
#' @author shi jian
AggregateHTseq <- function(inputDir = "G:/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/6.counts1"){
  library(data.table)
  library(dplyr)
  filepaths <- list.files(path=inputDir,recursive = T,full.names = T) 
  sample_names <- stringr::str_split(basename(filepaths),pattern = "\\.",simplify = T)[,1]
  temp <- fread(filepaths[1])
  temp2 <- fread(filepaths[2])
  temp3 <- temp %>% bind_cols(.,temp2$V2)
  temp3 <- temp
  temp <- c()
  for(i in 1:length(filepaths)){
    file_i <- fread(filepaths[i])
    if( i == 1){
      temp <- c(temp,file_i)
    }else{
      temp <- c(temp,file_i[,2])
    }
  }
  temp <- as.data.frame(temp)
  temp %>% tibble::column_to_rownames("V1") %>% data.table::setnames(.,old = c("V2","V2.1","V2.2","V2.3","V2.4","V2.5"),
                                                                     new = sample_names) -> temp1
  temp2 <- head(temp1,-5)
  return( temp2 )
}
