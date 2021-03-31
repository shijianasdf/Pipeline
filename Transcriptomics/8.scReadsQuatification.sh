#' 单细胞read count量化基于salmon
#函数一：make.index
#' @@ 建立salmon要求的index文件
#' @param transcriptomeFile 参考转录组文件绝对地址，fq或fq.gz等gz压缩格式均可，暂不支持其他压缩形式。推荐使用genecode，打开 https://www.gencodegenes.org/human/ 网址，选择Fasta files--> Transcript sequences下载
#' salmon和kallisto采用的是参考转录组，不建议使用参考基因组。如"/pub6/temp/xlw/kallisto_linux-v0.44.0/gencode.v29.transcripts.fa.gz"。
#' @param outfilepath 输出.sh文件
#' @param outDir index文件输出目录，需要事先建好。
#' @param type 软件使用类型，当前支持"kallisto"和"salmon"，注意输入格式要保持完全正确。
#' @param k 有效匹配的最小可接受长度，当前默认为31。注意这个参数需要用户确认read的长度范围，当fastq文件中大量read长度小于31时，软件处理会出现问题。
#' kallisto软件默认为31，也是支持的最大值。
#' salmon官方文档说明，k设置为31对于read长度为75bp或者更长效果要好，如果用户需要处理更短的read，可以考虑设置一个更小的k值，但并没有给出具体设置方案。

#' @return index文件。kallisto index是一个idx文件结尾的文件，而salmon index是一个文件夹

#' @note 
#' 当前使用Salmon-v0.9.1，kallisto 0.44.0版本
#' 运行前需要确认这两个软件是否存在在给定的运行路径下:salmon.path = "/pub5/xiaoyun/BioSoftware/Salmon-v0.9.1/bin/salmon";kallisto.path = "/pub5/xiaoyun/BioSoftware/kallisto_linux-v0.44.0/kallisto"
#' 对于salmon，支持两种方式建立index（Quasi-index and FMD-index-based modes），当前使用更优的Quasi-index模式
#' salmon建立index，可能会出现大片黄色警告，一般情况都可以忽略：通常是表示有些转录本长度不满足31的要求；salmon还会删除定义为重复的转录本，记录在index文件夹下的duplicate_clusters.tsv

#' example
#' make.index(transcriptomeFile = "/pub6/temp/xlw/kallisto_linux-v0.44.0/gencode.v29.transcripts.fa.gz", 
#'            outfilepath = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/salmon.sh",
#'            outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test", 
#'            type = "salmon", k = 31)
#' make.index(transcriptomeFile = "/pub6/temp/xlw/kallisto_linux-v0.44.0/gencode.v29.transcripts.fa.gz", 
#'            outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test", 
#'            type = "kallisto", k = 31)

make.index <- function(transcriptomeFile,
                       outfilepath,
                       outDir, 
					   type = c("kallisto", "salmon"), 
					   k = 31){   
        #salmon.path <- "/pub5/xiaoyun/BioSoftware/Salmon-v0.9.1/bin/salmon"
        #kallisto.path <- "/pub5/xiaoyun/BioSoftware/kallisto_linux-v0.44.0/kallisto"
        salmon.path <- "salmon"
		kallisto.path <- "kallisto"
        if(!dir.exists(outDir)){
		    dir.create(outDir)
            #stop(outDir, "目录不存在，请确认...")
        }
        # 获取参考转录本名字用于命名index文件
        ##格式通常为fa/fa.gz/fasta/fasta.gz
        gz <- grep("gz$", transcriptomeFile)
        if(length(gz)>0){
                faName <- gsub("\\.f[^.]*\\.gz$", "", transcriptomeFile, ignore.case = T)
        }else{
                faName <- gsub("\\.f[^.]*$", "", transcriptomeFile, ignore.case = T)
        }
        faName <- gsub(".*/", "", faName)
        
        if(type == "kallisto"){
                outIndexFile <- paste0(faName, ".idx")
                outIndexFile <- file.path(outDir, outIndexFile)
                index.command <- paste(kallisto.path, "index", "-i", outIndexFile, "-k", k, transcriptomeFile, sep = " ")
        }
        if(type == "salmon"){
                outIndexFile <- faName
                outIndexFile <- file.path(outDir, outIndexFile)
                index.command <- paste(salmon.path, "index -t", transcriptomeFile, "-i", outIndexFile, "--type puff", "-k", k, sep = " ")
        }
		writeLines(index.command,con=outfilepath)
        return(index.command)  
}
make.index(transcriptomeFile="/pub6/temp/shijian/gencode.v37.transcripts.fa.gz",
           outfilepath="/pub6/temp/shijian/NGSCommand/salmonIndex.sh",
           outDir="/pub6/temp/shijian/salmon_result", 
		   type = c("kallisto", "salmon")[2], 
		   k = 31)
./salmonIndex.sh
#salmon index -t /pub6/temp/shijian/gencode.v37.transcripts.fa.gz -i /pub6/temp/shijian/salmon_result/gencode.v37.transcripts --type puff -k 31

#salmon alevin -i /pub6/temp/shijian/salmon_result/gencode.v37.transcripts \
#    -l ISR  \
#    -1 /pub6/temp/shijian/trimGalore_result/SRS2476507_1_val_1.fq.gz \
#    -2 /pub6/temp/shijian/trimGalore_result/SRS2476507_2_val_2.fq.gz \
#    --chromium  \
#    -p 10 \
#    -o /pub6/temp/shijian/salmon_result/alevin_result \
#    --tgMap /pub6/temp/shijian/txTogene.tsv \
#	–dumpMtx \
#    --dumpFeatures
	
#' 构建salmon单细胞量化所需的 --tgMap对应的转录本对基因tsv文件	
#'quant.id              gene.id
#'ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|processed_transcript|  ENSG00000223972.5                         
#'ENST00000450305.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000002844.2|DDX11L1-201|DDX11L1|632|transcribed_unprocessed_pseudogene|  ENSG00000223972.5
	
single.ks.quantitation <- function(fastqDir,
                            outfilepath, 
                            is.TrimGalore = T, 
							pattern = c("fastq.gz$","fq.gz$")[2], 
							type = c("kallisto", "salmon")[2],
							protocol=c("--dropseq","--chromiumV3","--chromium","--gemcode","--citeseq","--celseq","--celseq2","--quartzseq2")[3],
							index.path, 
							outDir, 
							tgMap,
							pairEND=T, 
							extraParameter = NULL, 
							threads = 20){
        # 软件运行路径
        salmon.path <- "salmon"
        
        # 检查文件是否存在
        if(!dir.exists(fastqDir)){
             stop(fastqDir, " not exist!please check...")
        }
        
        # 如果outDir不存在且不为空，则自动建立存储目录
        if(!dir.exists(outDir)){
                if(!is.null(outDir)){
                    dir.create(outDir)
                }
        }
        
        # 用于指示type参数输出是否有误
        if(!(type %in% c("kallisto", "salmon"))){
                stop("请确认输入的type参数是否正确!")
        }
		
        commands <- c()
        # 双末端
        if(pairEND){
                paired.fastqfile1 <- list.files(fastqDir, pattern=paste0("_1.", pattern), full.names=TRUE, recursive=F)
                paired.fastqfile2 <- list.files(fastqDir, pattern=paste0("_2.", pattern), full.names=TRUE, recursive=F)
                paired.sampleName <- list.files(fastqDir, pattern=paste0("_1.", pattern), full.names=F, recursive=F)  
				
                # 衔接trim galore处理结果，获得更准确的样本名
                if(is.TrimGalore){
                        paired.sampleName <- gsub(paste0("_1_val_1.", pattern), "", paired.sampleName)
                }else{
                        paired.sampleName <- gsub(paste0("_1.", pattern), "", paired.sampleName)
                }
                
                if(type == "salmon"){
                        # 检查index是否存在且格式是否正确
                        if(!file.exists(index.path) || (length(grep(".idx$", index.path)) > 0)){
                                stop(index.path, " 不存在或格式不符合!please check...")
                        }
                        
                        salmon.command <- paste(salmon.path, "alevin -i", index.path, "-p", threads, "-l ISR", "--tgMap",tgMap,protocol,sep = " ")
                        
                        # 加入扩展参数
                        if(!is.null(extraParameter)){
                                salmon.command <- paste(salmon.command, extraParameter, sep = " ")
                        }
                        
                        for(i in 1:length(paired.sampleName)){
                                # 构建输出目录
                                sampleDir <- file.path(outDir, paired.sampleName[i])
                                if(!dir.exists(sampleDir)){
                                        dir.create(sampleDir)                                
                                }
                                paired.command <- paste(salmon.command, "-1", paired.fastqfile1[i], "-2", paired.fastqfile2[i], "-o", sampleDir, sep = " ")
                                print(paired.command)
                                commands <- c(commands,paired.command)
                        }
                }
                
        }else{
        # 单末端
                single.fastqfile <- list.files(fastqDir, pattern=pattern, full.names=TRUE, recursive=F)
                single.sampleName <- list.files(fastqDir, pattern=pattern, full.names=F, recursive=F)
                
                # 衔接trim galore处理结果，获得更准确的样本名
                if(is.TrimGalore){
                        single.sampleName <- gsub(paste0("_trimmed.", pattern), "", single.sampleName) 
                }else{
                        single.sampleName <- gsub(paste0(".", pattern), "", single.sampleName) 
                }
                
                # salmon定量
                if(type == "salmon"){
                        # 检查index是否存在且格式是否正确
                        if(!file.exists(index.path) || (length(grep(".idx$", index.path)) > 0)){
                                stop(index.path, " 不存在或格式不符合!please check...")
                        }
                        # 平均片段长度
                        if(is.null(fldMean)){
                                salmon.command <- paste(salmon.path, "quant -i", index.path, "-p", threads, "-l A", sep = " ")
                        }else{
                                salmon.command <- paste(salmon.path, "quant -i", index.path, "-p", threads, "-l A --fldMean", fldMean, sep = " ")
                        }
                        # 片段长度方差
                        if(!is.null(fldSD)){
                                salmon.command <- paste(salmon.command, "--fldSD", fldSD, sep = " ")
                        }
                        # 是否进行偏性矫正
                        if(is.Bias){
                                salmon.command <- paste0(salmon.command, " --seqBias --gcBias") 
                        }
                        # 加入扩展参数
                        if(!is.null(extraParameter)){
                                salmon.command <- paste(salmon.command, extraParameter, sep = " ")
                        }
                        for(i in 1:length(single.sampleName)){
                                # 构建输出目录
                                sampleDir <- file.path(outDir, single.sampleName[i])
                                if(!dir.exists(sampleDir)){
                                        dir.create(sampleDir)                                
                                }
                                single.command <- paste(salmon.command, "-r", single.fastqfile[i], "-o", sampleDir, sep = " ")
                                
                                print(single.command)                        
                                commands <- c(commands,single.command)
                        }
                }
        }
		writeLines(commands,con=outfilepath)
		return(commands) 
}	
single.ks.quantitation(fastqDir="/pub6/temp/shijian/trimGalore_result",
                       outfilepath="/pub6/temp/shijian/NGSCommand/salmonAlevin.sh", 
                       is.TrimGalore = T, 
					   pattern = c("fastq.gz$","fq.gz$")[2], 
					   type = c("kallisto", "salmon")[2],
					   protocol=c("--dropseq","--chromiumV3","--chromium","--gemcode","--citeseq","--celseq","--celseq2","--quartzseq2")[3],
					   index.path="/pub6/temp/shijian/salmon_result/gencode.v37.transcripts", 
					   outDir="/pub6/temp/shijian/salmon_result/alevin_result", 
					   tgMap="/pub6/temp/shijian/txTogene.tsv",
					   pairEND=T, 
					   extraParameter = NULL, 
					   threads = 20)
./salmonAlevin.sh	
