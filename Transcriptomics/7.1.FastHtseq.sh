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
