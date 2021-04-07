# FastQC检查fastq.gz文件质量，生成质量报告,每一个fastq.gz生成一个html文件和fastqc.zip文件
fastqc -o /pub6/temp/shijian/Fastq -t 10 /pub6/temp/shijian/Fastq/SRR5988127_1.fastq.gz /pub6/temp/shijian/Fastq/SRR5988127_2.fastq.gz

#' @description fastq质量评估基于fastQC
#' @param inputDir fastq.gz所在目录
#' @param outDir fastqc结果输出目录（包含每个fastq.gz对应的fastqc.zip）
#' @param outfilepath .sh文件
#' @param threads 线程数
#' @param pattern 文件名
qualityFASTQ <- function(inputDir,
						 outDir,
						 outfilepath,
						 threads,
						 pattern=".fastq.gz$"){
	files <- list.files(inputDir,pattern = pattern, full.names = TRUE, recursive = TRUE)
	if(!file.exists(outDir)){
		dir.create(outDir)
	}
	t.files <- paste(files,collapse=" ")
	command <- paste("fastqc -t",threads,"-o",outDir,t.files)
	writeLines(command,con=outfilepath)
	return(command)
}
qualityFASTQ(inputDir="/pub6/temp/shijian/Fastq",
			 outDir="/pub6/temp/shijian/Fastqc",
			 outfilepath="/pub6/temp/shijian/NGSCommand/Fastqc.sh",
		     threads=10,
             pattern=".fastq.gz$")
./Fastqc.sh

# multiQC整合fastQC生成的fastqc.zip文件 
multiqc test_7942raw_1_fastqc.zip test_7942raw_2_fastqc.zip -O /pub6/temp/shijian/Multiqc


#' @description multiqc质量评估基于fastQC结果
#' @param inputDir fastqc.zip所在目录
#' @param outDir multiqc结果输出目录
#' @param outfilepath 输出的.sh文件
#' @param pattern 文件名后缀
doMultiQC <- function(inputDir,
						 outDir,
						 outfilepath,
						 threads,
						 pattern=".fastqc.zip$"){
		files <- list.files(inputDir,pattern = pattern, full.names = TRUE, recursive = TRUE)
		if(!file.exists(outDir)){
			dir.create(outDir)
		}
		t.files <- paste(files,collapse=" ")
		command <- paste("multiqc",t.files,"-o",outDir)
		writeLines(command,con=outfilepath)
		return(command)				 						 
}
doMultiQC(inputDir="/pub6/temp/shijian/Fastqc",
			 outDir="/pub6/temp/shijian/Multiqc",
			 outfilepath="/pub6/temp/shijian/NGSCommand/multiqc.sh",
		     threads=10,
             pattern=".fastqc.zip$")
./multiqc.sh
