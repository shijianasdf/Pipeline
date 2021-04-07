#' @description  MACS2生成的bdg文件转为bigwig文件，基于UCSC的bedGraphToBigWig
#' @author shi jian
#' @usage 
#' http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/
#' conda install ucsc-fetchchromSizes
#' conda install ucsc-bedGraphtobigwig
#' fetchChromSizes hg38 > /pub6/Temp/sj/GSE77737/hg38.chrom.sizes
#' bedGraphToBigWig /pub6/Temp/sj/GSE77737/DNase-seq/MACS2/SRS1282560/V576_Dnase_I_Hypersensitiviy_treat_pileup.bdg \
#'                 /pub6/Temp/sj/GSE77737/hg38.chrom.sizes \
#'                 /pub6/Temp/sj/GSE77737/DNase-seq/MACS2/SRS1282560/V576_Dnase_I_Hypersensitiviy_treat_pileup.bw	
#'  目前有bug，因为bowtie2比对参考基因组是NCBI的参考基因组，MACS2得到的bdg文件也是以NCBI为准，
#'  而bedGraphToBigWig基于UCSC参考基因组去转，所以出现染色体chrEBV找不到的情况，因为NCBI和UCSC参考基因组名字有不一样的情况
#'  所以，以后只使用UCSC的参考基因组就可以了
#' 
#' https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chromToUcsc 
#' wget -c --no-check-certificate https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz  -O /pub6/Temp/sj/GSE77737/GRCh38_latest_genomic.fna.gz              
#' samtools faidx /pub6/Temp/sj/GSE77737/GRCh38_latest_genomic.fa | cut -f 1,2 /pub6/Temp/sj/GSE77737/GRCh38_latest_genomic.fa.fai > /pub6/Temp/sj/GSE77737/GRCH38.chromsizes
#' conda install pyfaidx
#' faidx /pub6/Temp/sj/GSE77737/GRCh38_latest_genomic.fa -i chromsizes > /pub6/Temp/sj/GSE77737/sizes.genome
#' @param bdgDir bdg所在的文件目录               
#' @param fetchChromSizes hg38.chrom.sizes文件地址
#' @param outfilepath sh文件输出地址
#' @param outDir 生成的hg38.chrom.sizes输出目录
#' @param ref 参考基因组类型
#' @param pattern  文件匹配模式              
#' @export 
FastbedGraph2BigWig <- function(bdgDir,
                                fetchChromSizes = NULL,
                                outDir,
                                ref=c("hg19","hg38")[2],
                                outfilepath,
                                pattern=".bdg$"){
  bdgfiles <- list.files(bdgDir,pattern = pattern,full.names = T,recursive = T)
  filename <- basename(list.files(bdgDir,pattern = pattern,full.names = F,recursive = T))
  library(stringr)
  filename <- str_split(filename,"\\.",simplify=T)[,1]
  
  commands <- c()
  if(is.null(fetchChromSizes)){
    fetchCMD <- paste( "fetchChromSizes",ref, ">", file.path(outDir,paste(ref,"chrom.sizes",sep=".")) )
    fetchChromSizes <- file.path(outDir,paste(ref,"chrom.sizes",sep="."))
    commands <- c(commands,fetchCMD)
    #system(fetchCMD)
  }
  for(i in 1:length(bdgfiles)){
    runCMD <- paste("bedGraphToBigWig",bdgfiles[i],fetchChromSizes,file.path(dirname(bdgfiles[i]),paste(filename[i],"bw",sep=".")))
    commands <- c(commands,runCMD)
    #system(runCMD)
  }
  writeLines(commands,con = outfilepath)
  return(commands)
}
FastbedGraph2BigWig(bdgDir="/pub6/Temp/sj/GSE77737/DNase-seq/MACS2",
                    fetchChromSizes = NULL,
                    outDir="/pub6/Temp/sj/GSE77737",
                    ref=c("hg19","hg38")[2],
                    outfilepath="/pub6/Temp/sj/GSE77737/NGScommands/Dase_bdg2bigwig.sh",
                    pattern=".bdg$")
FastbedGraph2BigWig(bdgDir="/pub6/Temp/sj/GSE77737/DNase-seq/MACS2",
                    fetchChromSizes = "/pub6/Temp/sj/GSE77737/GRCH38.chromsizes",
                    outDir="/pub6/Temp/sj/GSE77737",
                    ref=c("hg19","hg38")[2],
                    outfilepath="/pub6/Temp/sj/GSE77737/NGScommands/Dase_bdg2bigwig.sh",
                    pattern=".bdg$")
./Dase_bdg2bigwig.sh 2>> Dase_bdg2bigwig.log
