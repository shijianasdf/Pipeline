### 运行featureCounts计数
#featureCounts -T 5 \
#              -p \
#              -t exon \ 
#              -g gene_id \
#              -s 0 \
#              -a /teach/database/gtf/gencode.v25.annotation.gtf.gz \
#              -o ~/6.featureCounts/all.id.txt \
#              *sort.bam

#-a：基因组注释文件（GTF/GFF格式），定义基因位置结构。
#-o：输出结果文件名。
#-t：特征类型（Feature type），默认为 exon（外显子）。
#-g：属性类型（Attribute type），默认为 gene_id，按基因级别汇总。
#-T：线程数，用于加速比对结果的统计。
#-p：针对双端测序（Paired-end）数据，必须添加此参数。
#-s：链特异性参数（0=未特异性, 1=正链, 2=负链）。
runFeatureCount <- function(bamDir,pattern = "sortedByCoord.out.bam$",process=10,paired=T,
                            strand=c(0,1,2)[1],type=c("exon","CDS")[1],
                            gtf.path,out.dir,outfilepath){
  library(stringr)
  bamfiles <- list.files(bamDir,full.names = T,recursive = T,pattern = pattern)
  patients <- str_split(dirname(bamfiles),"/",simplify = T)[,dim(str_split(dirname(bamfiles),"/",simplify = T))[2]]
  if(!file.exists(out.dir))
    dir.create(out.dir,recursive = T)
  if(paired == T){
    command <- paste0("featureCounts -T ",process," -p -t ",type, " -g gene_id -s ",strand," -a ",gtf.path)
  }else{
    command <- paste0("featureCounts -T ",process," -t ",type, " -g gene_id -s ",strand," -a ",gtf.path)
  }
  commands <- c()
  for(i in 1:length(bamfiles)){
    tempcommand <- paste0(command," -o ",file.path(out.dir,paste0(patients[i],".txt"))," ",bamfiles[i])
    commands <- c(commands,tempcommand)
  }
  writeLines(commands,con = outfilepath)
  return(commands)
}

runFeatureCount(bamDir="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0155_RNA_9hsa/5.star_bam",
                pattern = "sortedByCoord.out.bam$",process=10,paired=T,
                strand=c(0,1,2)[1],type=c("exon","CDS")[1],
                gtf.path="/data/shijian/refData/human_reference/gencode.v48.chr_patch_hapl_scaff.annotation.gtf",
                out.dir="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0155_RNA_9hsa/5.featurecounts",
                outfilepath="/data/shijian/project/NGSCommand/featurecounts.sh")
#nohup bash /data/shijian/project/NGSCommand/featurecounts.sh > featurecounts.log 2>&1 &

#' @description 
#' 对多个样本featurecount结果提取整理成表达矩阵
#' @param  inputDir featurecount定量的count文件
#' @param  pattern  选取文件名后缀
AggregateFeatureCount <- function(inputDir = "/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0155_RNA_9hsa/5.featurecounts",
                                  pattern=".txt$"){
  library(data.table)
  library(dplyr)
  library(magrittr)
  filepaths <- list.files(path=inputDir,recursive = T,full.names = T,pattern=pattern) 
  sample_names <- stringr::str_split(basename(filepaths),pattern = "\\.",simplify = T)[,1]
  #file_i <- fread(filepaths[1])
  temp <- c()
  for(i in 1:length(filepaths)){
    file_i <- fread(filepaths[i])
    file_i <- file_i[,c(1,7)]
    if( i == 1){
      temp <- c(temp,file_i)
    }else{
      temp <- c(temp,file_i[,2])
    }
  }
  temp <- as.data.frame(temp)
  temp %<>% tibble::column_to_rownames("Geneid") 
  colnames(temp) <- sample_names
  return( temp )
}

