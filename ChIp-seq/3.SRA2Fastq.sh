# 解压sra文件为fastq.gz文件
fastq-dump --gzip --split-files SRR6232298.sra -O /pub6/temp/shijian/SRP116382

#'@describeIn 拼接从SRA提取fastq的批量命令行
#'@param inputDir sra文件目录
#'@param outpath 命令行输出地址
#'@return 返回所有fastq-dump命令行
# fastq-dump --gzip --split-files SRR6232298.sra -O /pub6/temp/shijian/SRP116382
BatchSRA2fastq<- function(inputDir,
                          outDir,
						  outpath,
						  pattern=".sra$"){
  files <- list.files(inputDir,pattern = pattern, full.names = TRUE, recursive = TRUE)
  pre_command <- "fastq-dump --gzip --split-files"
  if(!file.exists(outDir)){
   dir.create(outDir)
  }
  commands <- c()
  for(i in files){
    temp <- paste(pre_command,i,outDir)
    commands <- c(commands,temp)
  }
  writeLines(commands,con=outpath)  
}
BatchSRA2fastq(inputDir="/home/shijian2015/ncbi/public/sra",outDir="/home/shijian2015/ncbi/public/fastq",outpath="/home/shijian2015/ncbi/public/sra/fastqdump.sh")
./fastqdump.sh

# 利用parallel-fastq-dump解压SRA
parallel-fastq-dump -s D:/Rsources/Project/StudySingleCell/data/SRA/SRR5988127.sra -t 10 --split-3 --gzip -O D:/Rsources/Project/StudySingleCell/data/fastq

#' @describeIn 拼接从SRA提取fastq的批量命令行（基于parallel-fastq-dump）
#' @param filePath: 待转换的.sra文件存放路径及文件名（仅用于每次处理单个SRA文件转换时）如：/IData/SingleCell/scRNAseq/Brain/SRP092584/SRR4919521/SRR4919521.sra
#' @param outfilepath: .sh文件存储地址
#' @param outDir: fastq.gz文件存储地址
# 在outDir为NULL情况下，对于批量sra文件和单个sra文件处理，会有不同的处理：
# 多sra文件处理，会自动在dirPath路径下建立名为fastq文件夹存储转换后的结果文件
# 单sra文件处理，会自动在sra文件所在目录下生成转换的fastq文件
#' @param dirPath: 待转换的.sra文件集的文件路径（仅用于SRA文件批量转换时）
#' @param pattern: 模式匹配sra文件，默认为".sra$"
#' @param pairEND: TURE/FALSE，TURE则表示为双末端测序，FALSE则表示为单末端测序（单个sra文件可能存储了pair-end的两个配对的read，因此，该参数可以从SRA文件中分别提取出两个_1与_2的FASTQ文件）
#' @param threads: 并行线程数目。一般默认为20，适合批量sra文件转换
#' @return 返回parallel-fastq-dump命令行
convertSRA2FASTQ <- function(filePath = NULL, 
                             outfilepath = NULL,
							 outDir = NULL, 
							 dirPath = NULL, 
							 pattern = ".sra$", 
							 threads = 20, 
							 pairEND = FALSE){
	# 运行命令
	SRA.dump.fastq <- "parallel-fastq-dump"
	
	if(!file.exists(outDir)){
     dir.create(outDir)
    }
	
	commands <- c()
	# 批量处理sra文件
	if(!is.null(dirPath)){
		# 匹配所有的sra文件
		t.files <- list.files(dirPath, pattern = pattern, full.names = TRUE, recursive = TRUE)
		# 确定样本名字，设定.sra前面的字段为样本名字 
		sampleNames <- list.files(dirPath, pattern = pattern, full.names = TRUE, recursive = TRUE)
		sampleNames <- gsub(".*/", "", sampleNames)
		sampleNames <- gsub(".sra$", "", sampleNames)
		for(i in 1:length(t.files)){
			if(!pairEND){
				t.command <- paste(SRA.dump.fastq, "-s", t.files[i], "-t", threads, "--gzip -O", outDir, sep = " ")
			    commands <- c(commands,t.command)
			}else{
				# split-3：针对pairend 测序结果，将原文件转换成3个fastq文件，文件名其中一个是*_1.fastq,一个是*_2.fastq,这两个文件里面的reads是成对的，而*.fastq里面的read是剩余不成对的
				t.command <- paste(SRA.dump.fastq, "-s", t.files[i], "-t", threads, "--split-3 --gzip -O", outDir, sep = " ")
			    commands <- c(commands,t.command)
			}
		}
	}else{
	#单样本处理
		if(pairEND){
			commands <- paste(SRA.dump.fastq, "-s", filePath, "-t", threads, "--split-files --gzip -O", outDir, sep = " ")
		}else{
			commands <- paste(SRA.dump.fastq, "-s", filePath, "-t", threads, "--gzip -O", outDir, sep = " ")
		}
	}
	writeLines(commands,con=outfilepath)
	return(commands)
}
convertSRA2FASTQ(outDir = "/pub6/temp/shijian/Fastq", 
                 outfilepath="/pub6/temp/shijian/NGSCommand/parallelfastqdump.sh",
				 dirPath = "/pub6/temp/shijian/SRP116382", 
				 pairEND = T, threads = 10)
./parrallelfastqdump.sh
