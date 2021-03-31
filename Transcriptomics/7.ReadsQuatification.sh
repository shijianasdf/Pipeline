####################################################################################
# 整体功能简介：
# 使用kallisto/salmon进行bulk转录组定量表达
###################################################################################
# 重要使用函数简介
# 代码中共包括三个函数：
# 1）make.index，用于建立kallisto/salmon的index文件
# 2）ks.quantitation，实现kallisto/salmon的定量转录本表达
# 3）transcript2gene.quant，使用tximport包将kallisto/salmon定量的转录本表达转换为基因表达水平
###################################################################################


###################################################################################
# 创建作者: 龙志林
# 日期：2018/08/28
# 论坛网址：http://210.46.85.145/showtopic-2279.aspx
###################################################################################

# 使用前推荐先看一下论坛上salmon和kallisto的总结
# kallisto: http://210.46.85.145/showtopic-2243.aspx
# salmon: http://210.46.85.145/showtopic-2272.aspx
###################################################################################

###################################################################################
# 修改日志
###################################################################################
# 1) 日期：               修改者：
#    简要修改内容：
# 
# 

###################################################################################
#函数一：make.index
#' @@ 建立kallisto或salmon要求的index文件
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

###################################################################################
#函数二：ks.quantitation
#' @param fastqDir fastq文件目录，如/IData/CancerOMICS/SkinCancer/Riaz_Cell_2017/RNA_Seq/trim_galore
#' @param is.TrimGalore TRUE/FALSE，是否是经过trim galore处理过的fastq数据
#' @param pattern fastq文件的结尾模式，用于匹配fastqDir目录下的fastq文件，如"fastq/fq.gz"或"fastq/fq"，不能为".fastq/fq.gz"或".fastq/fq"
#' @param type 软件使用类型，当前支持kallisto和salmon。注意输入要保持正确
#' @param index.path 转录组index文件或目录，
#' salmon要求index目录，如/home/longzl2017/dataScience/PipeLine/RNA-seq/test/gencode_v28_transcripts_index
#' kallisto要求index文件，如/home/longzl2017/dataScience/PipeLine/RNA-seq/test/Homo_gencodeV29_transcripts.idx
#' @param outDir 输出目录，要求具有可写权限，如/IData/CancerOMICS/SkinCancer/Riaz_Cell_2017/RNA_Seq/quantitation_result。如果不存在，则会自动建立
#' @param pairEND: TRUE/FALSE, 是否是双末端测序, 这个一定要仔细查证之后再给定，因为kallisto或salmon处理单双末端参数不同
#' @param fldMean 估计的片段平均长度(mean fragment lenth of the sequencing library)，只适用单末端模式，且是kallisto必须指定的重要参数，会直接影响定量的准确性。
#' 当前采用kallisto官网的默认值100，实际使用时最好进一步查证，根据自己数据的实际情况进行修改设定。软件fastqc可以查看read长度分布情况
#' salmon官网上并没有强制要求单末端时必须指定该参数，所以type参数为salmon时，允许设置为NULL，但是明确说明了该参数很重要，会直接影响定量的准确性，使不使需要自己斟酌
#' @param fldSD 估计的片段长度的标准差(standard deviation of the fragment lenth distribution of the sequencing library)，只适用单末端模式，且是kallisto必须指定的重要参数，会直接影响定量的准确性。
#' salmon官网上并没有强制要求单末端时必须指定该参数，所以type参数为salmon时，允许设置为NULL，但是明确说明了该参数很重要，会直接影响定量的准确性，使不使需要自己斟酌
#' @param strand.specific 是否使用链特异模式运行kallisto，可选参数，默认为 "" ;支持"fr-stranded"或"rf-stranded"输入，"fr-stranded"指用链特异模式运行kallisto
#' "rf-stranded" 也指链特异，但first read 比对到转录本的反义链上。只支持kallisto
#' @param is.Bias 是否进行偏性矫正，默认为FALSE，不进行偏性矫正。kallisto支持基于sequence的bias矫正；salmon支持矫正序列特异的bias和fragment-level GC的bias。个人推荐使用，用户可以进一步查证自行选择。
#' @param is.BAMoutput 是否需要BAM格式的输出，默认FALSE，不输出BAM文件，只支持kallisto
#' @param extraParameter kallisto或salmon的自定义扩展参数字符串，可以加入不在函数设定参数内的运行参数，默认为NULL。输入需要完全按照参数组合方式排列，如"--numBootstraps 40"
#' @param threads 线程数目，默认为20 

#' @return 转录本定量结果文件。
#### kallisto会为每个样本3个文件： 详细说明见 --- http://210.46.85.145/showtopic-2243.aspx
#' abundances.tsv 是转录本表达丰度定量文件，不包含bootstrap estimates. 使用--plaintext 参数输出文本格式的abundance estimates。或者用h5dump命令把输出的HDF5文件转成文本格式。第一行是每列的列名，包含estimated counts, TPM, effective length.
#' run_info.json 是一个 json 文件，包含程序运行参数log的信息。

#### salmon会为每个样本产生一个quant.sf(定量结果)和很多其他记录文件：详细说明见 --- http://210.46.85.145/showtopic-2272.aspx
#' quant.sf 是转录本表达丰度定量文件，包含Name(转录本名字)、Length(length of the target transcript in nucleotides)、EffectiveLength(估计转录本有效长度)、TPM(估计的转录本TPM表达丰度)、NumReads(估计的比对到转录本上的read数目)
#' cmd_info.json 是程序运行参数记录文件
#' lib_format_counts.json记录了salmon估计的样本中所有read的文库协议分布情况
#' 其他文件可从 https://salmon.readthedocs.io/en/latest/file_formats.html#fileformats 进一步查阅

#' @note 
#' 内嵌的软件执行路径： 运行前需要确认一下
#' Salmon: /pub5/xiaoyun/BioSoftware/Salmon-v0.9.1/bin/salmon
#' kallisto: /pub5/xiaoyun/BioSoftware/kallisto_linux-v0.44.0/kallisto
#' Samtools: /pub5/xiaoyun/BioSoftware/Samtools-1.3/bin
#' 
#' 当前使用Salmon-v0.9.1，kallisto 0.44.0版本。salmon相对比kallisto运行速度更快
#' kallisto和salmon定量结果是转录本水平的，所以我们需要额外的一步转换(transcript2gene.quant.R)，将转录本表达转换为基因水平的表达
#' 对于多个run的样本，需要在定量前将同一个样本多个run合并
#' 
#### salmon
#' salmon定量默认使用参数：-l A 自动确定文库类型；
#' 使用salmon定量时输出控制台会显示出一段黄色输出：显示了有多少比例的read小于设定的K值和多少read比对上了，这段信息对于参数是否设置正确具有很好的指导意义，需要着重注意!!!
#### kallisto
#' kallisto定量默认使用参数：--plaintext 使用普通文本替代HDF5文件输出

#' example
#' single-end kallisto
#' ks.quantitation(fastqDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/fastq", is.TrimGalore = F, pattern = "fastq.gz", type = "kallisto", index.path = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/Homo_gencodeV29_transcripts.idx", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/fastq/kresult", pairEND = F, fldMean = 100, fldSD = 20)
#' ks.quantitation(fastqDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/trimGalore_result", is.TrimGalore = T, type = "kallisto", pattern = "fq.gz$", index.path = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/Homo_gencodeV29_transcripts.idx", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/trimGalore_result/kresult", pairEND = F, fldMean = 100, fldSD = 20)
#' 
#' single-end salmon
#' ks.quantitation(fastqDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/fastq", is.TrimGalore = F, pattern = "fastq.gz", type = "salmon", index.path = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/gencode_v29_transcripts_index", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/fastq/sresult", pairEND = F, fldMean = 100, fldSD = 20)
#' ks.quantitation(fastqDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/trimGalore_result", is.TrimGalore = T, type = "salmon", pattern = "fq.gz$", index.path = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/gencode_v29_transcripts_index", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/trimGalore_result/sresult", pairEND = F, fldMean = 100, fldSD = 20)
#' 
#' paired-end kallisto
#' ks.quantitation(fastqDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/paired/fastq", is.TrimGalore = F, type = "kallisto", pattern = "fastq.gz", index.path = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/Homo_gencodeV29_transcripts.idx", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/paired/fastq/kresult", pairEND = T)
#' ks.quantitation(fastqDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/paired/trimGalore_result", is.TrimGalore = T, type = "kallisto", pattern = "fq.gz$", index.path = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/Homo_gencodeV29_transcripts.idx", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/paired/trimGalore_result/kresult", pairEND = T)
#' 
#' paired-end salmon
#' ks.quantitation(fastqDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/paired/fastq", is.TrimGalore = F, type = "salmon", pattern = "fastq.gz", index.path = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/gencode_v29_transcripts_index", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/paired/fastq/sresult", pairEND = T)
#' ks.quantitation(fastqDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/paired/trimGalore_result", is.TrimGalore = T, type = "salmon", pattern = "fq.gz$", index.path = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/gencode_v29_transcripts_index", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/paired/trimGalore_result/sresult", pairEND = T)

ks.quantitation <- function(fastqDir,
                            outfilepath, 
                            is.TrimGalore = T, 
							pattern = c("fastq.gz$","fq.gz$")[2], 
							type = c("kallisto", "salmon")[2], 
							index.path, 
							outDir, 
							pairEND=T, 
							fldMean = 100, 
							fldSD = 20, 
							is.Bias = FALSE, 
							is.BAMoutput = FALSE, 
							strand.specific = "", 
							extraParameter = NULL, 
							threads = 20){
        
        # 软件运行路径
        salmon.path <- "salmon"
        kallisto.path <- "kallisto"
        
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
                
                if(type == "kallisto"){
                        # 检查index是否存在且格式是否正确
                        if(!file.exists(index.path) || (length(grep(".idx$", index.path)) == 0)){
                                stop(index.path, " 不存在或格式不符合!please check...")
                        }
                        
                        kallisto.command <- paste(kallisto.path, "quant -i", index.path, "-t", threads, "--plaintext", sep = " ")
                        
                        # 是否进行偏性矫正
                        if(is.Bias){
                                kallisto.command <- paste0(kallisto.command, " --bias") 
                        }
                        
                        # 是否是链特异
                        if(nchar(strand.specific) > 0){
                                kallisto.command <- paste0(kallisto.command, " --", strand.specific)
                        }
                        
                        # 加入扩展参数
                        if(!is.null(extraParameter)){
                                kallisto.command <- paste(kallisto.command, extraParameter, sep = " ") 
                        }
                        
                        for(i in 1:length(paired.sampleName)){
                                # 构建输出目录
                                sampleDir <- file.path(outDir, paired.sampleName[i])
                                if(!dir.exists(sampleDir)){
                                        dir.create(sampleDir)                                
                                }
                        
                                ##---------如果要输出BAM文件
                                if(is.BAMoutput){
                                        paired.command <- paste0(kallisto.command, " -o ", sampleDir, " --pseudobam ", paired.fastqfile1[i], " ", paired.fastqfile2[i], " | samtools view -Sb - > ", sampleDir, "/out.bam")
                                }else{
                                        paired.command <- paste(kallisto.command, "-o", sampleDir, paired.fastqfile1[i], paired.fastqfile2[i], sep = " ")
                                } 
                                
                                print(paired.command)
                                commands <- c(commands,paired.command)
                        }
                }
                
                if(type == "salmon"){
                        # 检查index是否存在且格式是否正确
                        if(!file.exists(index.path) || (length(grep(".idx$", index.path)) > 0)){
                                stop(index.path, " 不存在或格式不符合!please check...")
                        }
                        
                        salmon.command <- paste(salmon.path, "quant -i", index.path, "-p", threads, "-l A", sep = " ")
                        
                        # 是否进行偏性矫正
                        if(is.Bias){
                                salmon.command <- paste0(salmon.command, " --seqBias --gcBias") 
                        }
                        
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
                        
                # kallisto定量
                if(type == "kallisto"){
                        # 检查index是否存在且格式是否正确
                        if(!file.exists(index.path) || (length(grep(".idx$", index.path)) == 0)){
                                stop(index.path, " 不存在或格式不符合!please check...")
                        }
                        
                        kallisto.command <- paste(kallisto.path, "quant -i", index.path, "-t", threads, "--single -l", fldMean, "-s", fldSD, "--plaintext", sep = " ")
                        
                        # 是否进行偏性矫正
                        if(is.Bias){
                                kallisto.command <- paste0(kallisto.command, " --bias") 
                        }
                        
                        # 是否是链特异
                        if(nchar(strand.specific) > 0){
                                kallisto.command <- paste0(kallisto.command, " --", strand.specific)
                        }
                        
                        # 加入扩展参数
                        if(!is.null(extraParameter)){
                                kallisto.command <- paste(kallisto.command, extraParameter, sep = " ") 
                        }
                        
                        for(i in 1:length(single.sampleName)){
                                # 构建输出目录
                                sampleDir <- file.path(outDir, single.sampleName[i])
                                if(!dir.exists(sampleDir)){
                                        dir.create(sampleDir)                                
                                }
                                ##---------如果要输出BAM文件
                                if(is.BAMoutput){
                                        single.command <- paste0(kallisto.command, " -o ", sampleDir, " --pseudobam ", single.fastqfile[i], " | samtools view -Sb - > ", sampleDir, "/out.bam")
                                }else{
                                        single.command <- paste(kallisto.command, "-o", sampleDir, single.fastqfile[i], sep = " ")
                                }
                                
                                print(single.command)                        
                                commands <- c(commands,single.command)
                        }
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
ks.quantitation(fastqDir="/pub6/temp/shijian/trimGalore_result",
                outfilepath="/pub6/temp/shijian/NGSCommand/quantitation.sh", 
                is.TrimGalore = T, 
				pattern = "fq.gz$", 
				type = c("kallisto", "salmon")[2], 
				index.path="/pub6/temp/shijian/salmon_result/gencode.v37.transcripts", 
				outDir="/pub6/temp/shijian/salmon_result", 
				pairEND=T, 
				fldMean = 100, 
				fldSD = 20, 
				is.Bias = FALSE, 
				is.BAMoutput = FALSE, 
				strand.specific = "", 
				extraParameter = NULL, 
				threads = 20)
./quantitation.sh
salmon quant -i /pub6/temp/shijian/salmon_result/gencode.v37.transcripts -p 20 -l A -1 /pub6/temp/shijian/trimGalore_result/SRS2476507_1_val_1.fq.gz -2 /pub6/temp/shijian/t
rimGalore_result/SRS2476507_2_val_2.fq.gz -o /pub6/temp/shijian/salmon_result/SRS2476507
###################################################################################
#函数三：transcript2gene.quant
####原理
#' tximport package支持salmon和kallisto进行转录本和基因表达水平的转换
#' tximport只是导入kallisto/salmon定量的转录本TPM和count值，将同一个基因转录本异构体加和，计算出基因水平的TPM和count

####支持证据
#' salmon官网上有直接说明可以使用tximport包将转录本水平表达转换为基因水平表达
#' kallisto并没有直接的说明可以直接用tximport进行转录本和基因水平表达的转换，kallisto的开发人员推荐使用他们开发的专门适用kallisto下游分析sleuth软件进行处理
#' 但是确实有许多发表的文章使用了kallisto+tximport衔接运行

#' @param quantDir kallisto/salmon定量结果目录，如"/home/longzl2017/dataScience/PipeLine/RNA-seq/test/fastq/kresult"， "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/fastq/sresult"
#' @param outDir 结果输出目录，如"/pub6/temp/immunophenoscore/Data/28552987" 
#' @param type 使用软件类型。目前只支持salmon和kallisto 
#' @param is.RefGencode TRUE/FALSE，kallisto/salmon建立index时使用的参考转录组是否是gencode，推荐使用gencode转录组。如果为真，则tx2gene可以为NULL。但如果用户给定了tx2gene，即使是gencode注释类型，程序也会直接采用用户给定的tx2gene
#' @param tx2gene data.frame，转录本和基因的对应关系，要求第一列为转录本id，第二列为基因id，要求和参考转录组使用的版本一致，否则可能会出现转录本和基因的版本不一致情况，导致无法匹配。
#' 如果gencodeAnnotation为TRUE时，则可以为NULL
#' @param countsFromAbundance 设置从abundance估计count的方式。支持三种方式"no"、"scaledTPM"和"lengthScaledTPM"，通常情况下选择"no"。
#' tximport默认为no，即直接从定量软件中加和同一个基因转录本count值；
#' scaledTPM，基于library size和定量软件计算出TPM重新估计count；
#' lengthScaledTPM，基于library size、所有样本中平均转录本长度和定量软件计算出TPM重新估计count

#' @returnType list
#' @return gene.expression.matrix 包含3个matrix对象：TPM, counts, length。
#' TPM为基因表达TPM矩阵，counts为基因表达count矩阵，length为每个基因的平均转录本长度

#' @note 
#' 当前使用tximport package版本为1.6.0
#' scaledTPM、lengthScaledTPM两种方式，基于TPM重新估计每个基因的count，而不是直接加和软件估计的转录本count值，这种方式能够帮助矫正文库差异和转录本长度

#' example
# 非gencode注释 --- kallisto
#' load(file = "/pub6/temp/xlw/kallisto_linux-v0.44.0/Homo_sapiens.ENSG2ENST2.RData")
#' tx2gene <- data.frame(transcript_id = Homo_sapiens.ENSG2ENST$transcript_id, gene_id = Homo_sapiens.ENSG2ENST$Ensembl_id)，格式可以参考 https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
#' transcript2gene.quant(quantDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/fastq/kresult", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/fastq/kresult", type = "kallisto", tx2gene = tx2gene, countsFromAbundance = "no")

# gencode注释 --- kallisto
#' transcript2gene.quant(quantDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/fastq/kresult", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/fastq/kresult", type = "kallisto", is.RefGencode = T, tx2gene = NULL, countsFromAbundance = "no")

# gencode注释 --- salmon
#' transcript2gene.quant(quantDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/fastq/sresult", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/fastq/sresult", type = "salmon", is.RefGencode = T, tx2gene = NULL, countsFromAbundance = "no")

transcript2gene.quant <- function(quantDir, 
                                  outDir, 
								  type = c("salmon", "kallisto")[1], 
								  is.RefGencode = FALSE, 
								  tx2gene, 
								  countsFromAbundance = c("no", "scaledTPM", "lengthScaledTPM")){
        
        if(!dir.exists(quantDir)){
                stop(quantDir, "not exist!please check...")
        }
        
        # 判断type参数输入是否正确
        if(!(type %in% c("salmon", "kallisto"))){
                stop("type参数输入格式不符!")
        }
        
        if(type == "kallisto"){
                pattern <- "abundance.tsv"
        }
        if(type == "salmon"){
                pattern <- "quant.sf"
        }
        
        quantFiles <- list.files(quantDir, pattern=paste0("^", pattern, "$"), full.names=TRUE, recursive=T)
        if(length(quantFiles) == 0){
              stop(quantDir, "没有匹配相应类型的文件，请确认目录下是否有文件或输入的type参数和定量文件软件来源是否一致!")
        }
        
        sampleName <- list.files(quantDir, pattern=paste0("^", pattern, "$"), full.names=F, recursive=T)
        sampleName <- gsub(paste0("/", pattern), "", sampleName)
        names(quantFiles) <- sampleName
        
        #由于使用gencode中定量出quant.sf或abundance.tsv文件包含转录本和基因名字信息，且样本间是一致的，所以可以利用第一列信息来构建tx2gene
        if(is.RefGencode){
                # 如果用户给定了tx2gene，即使为gencode注释类型也采用用户给定的tx2gene
                if(is.null(tx2gene)){
                        geneInfo <- read.table(file = quantFiles[1], sep = "\t", stringsAsFactors = F, header = T)
                        geneInfo <- geneInfo[,1]
                        infoMatrix <- sapply(geneInfo, function(x){
                                attr <- unlist(strsplit(x, "\\|"))
                                return(attr)
                        })
                        infoMatrix <- t(infoMatrix)
                        tx2gene <- data.frame(quant.id = rownames(infoMatrix), gene.id = infoMatrix[, 2])  
                }
        }else{
                if(is.null(tx2gene)){
                        stop("如果不是gencode注释体系，需要给定tx2gene文件")
                }
        }
        library(tximport)
        # 判断countsFromAbundance参数输入是否正确
        if(!(countsFromAbundance %in% c("no", "scaledTPM", "lengthScaledTPM"))){
            stop("countsFromAbundance参数输入格式不符!")
        }
        # txOut = FALSE表示估计基因水平的表达
        gene.expression.matrix <- tximport(files = quantFiles, type = type, tx2gene = tx2gene, txOut = FALSE, countsFromAbundance = countsFromAbundance)
        gene.TPM.matrix <- gene.expression.matrix$abundance
        gene.count.matrix <- gene.expression.matrix$counts
        gene.length.matrix <- gene.expression.matrix$length
        
        gene.expression.matrix <- list(TPM = gene.TPM.matrix, count = gene.count.matrix, length = gene.length.matrix)
        save(gene.expression.matrix, file = file.path(outDir, "gene.expression.matrix.RData"))
}
transcript2gene.quant(quantDir = "/pub6/temp/shijian/salmon_result/Expressionquantification", 
                      outDir = "/pub6/temp/shijian/salmon_result/Expressionquantification", 
					  type = "salmon", is.RefGencode = T, tx2gene = NULL, countsFromAbundance = "no")

