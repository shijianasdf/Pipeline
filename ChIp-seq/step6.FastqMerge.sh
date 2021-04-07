#' @description 合并一个样本多个run的fastq文件
#' @author dragon 2017.12.25
#' @param trim_galore_resultDir 字符串，数据分析结果存储路径
#' @param paired TRUE or FALSE，数据测序类型，单末端还是双末端
#' @return 在数据存储路径下的trim_galore产生合并多run的fq.gz文件

#' @note 要求basePath下存在SRPSampleInfo.RData文件
mergeFastq <- function(trim_galore_resultDir,
                       sampleInfoFile, 
					   paired = TRUE){
	if(!file.exists(sampleInfoFile)){
		return(cat("ERROR:", sampleInfoFile, "not exist!\n"))
	}
	load(file = sampleInfoFile)
	#检查fastq文件是否存在
	if(!dir.exists(trim_galore_resultDir)){
		return(cat(trim_galore_resultDir, "no exist!\n"))
	}
	if(paired){
		fastqFiles1 <- list.files(trim_galore_resultDir, pattern = "_1_val_1.fq.gz$", full.names = TRUE)
		fastqFiles2 <- list.files(trim_galore_resultDir, pattern = "_2_val_2.fq.gz$", full.names = TRUE)
		sampleNames <- list.files(trim_galore_resultDir, pattern = "_1_val_1.fq.gz$", full.names = F)
		sampleNames <- gsub("_1_val_1.fq.gz$", "", sampleNames)
	}else{
		fastqFiles <- list.files(trim_galore_resultDir, pattern = "_trimmed.fq.gz$", full.names = TRUE)
		sampleNames <- list.files(trim_galore_resultDir, pattern = "_trimmed.fq.gz$", full.names = F)
		sampleNames <- gsub("_trimmed.fq.gz$", "", sampleNames)
	}
	##基于sample_accession分辨具有多个run的样本
	sample_accession <- SRPSampleInfo$sample_accession
	run_table <- table(sample_accession)
	index <- which(run_table > 1 )
	if(length(index) == 0){
		return(cat("no sample of multi-run!\n"))
	}else{
		sample_uni_accession <- unique(sample_accession)
		a <- sapply(sample_uni_accession, function(x){
			index <- which(sample_accession == x)
			run_id <- SRPSampleInfo$run_accession[index]
			if(paired){
				index <- match(run_id, sampleNames)
				#由于sra转fastq时，是选择了部分样本，所以fastq文件可能会少于SRPSampleInfo中的样本数，所以会存在对应NA的情况
				#由于不可能出现多个run，部分run对应不同的属性，所以全为NA时不进行处理
				#存在多run样本，部分run由于转换错误而被删除，从而导致出现部分NA
				index1 <- which(is.na(index))
				if(length(index1) == 0 ){
					fastq1 <- paste(fastqFiles1[index], collapse = " ")
					fastq2 <- paste(fastqFiles2[index], collapse = " ")
					merge.command1 <- paste0("cat ", fastq1, " > ", trim_galore_resultDir, "/", x, "_1_val_1.fq.gz")
					system(merge.command1)
					merge.command2 <- paste0("cat ", fastq2, " > ", trim_galore_resultDir, "/", x, "_2_val_2.fq.gz")
					system(merge.command2)

					#必须删除单独的run的fastq文件，否则会影响后面结果
					delete.command1 <- paste("rm", fastq1, sep = " ")
					system(delete.command1)
					delete.command2 <- paste("rm", fastq2, sep = " ")
					system(delete.command2)
				}else{
					index2 <- na.omit(index)
					if(length(index2) > 0){
						#由于部分run缺失，导致样本无法进行run合并，删除样本对应剩余的run
						fastq1 <- paste(fastqFiles1[index2], collapse = " ")
						fastq2 <- paste(fastqFiles2[index2], collapse = " ")
						system(paste("rm", fastq1, sep = " "))
						system(paste("rm", fastq2, sep = " "))
					}
				}
			}else{
				index <- match(run_id, sampleNames)
				index1 <- which(is.na(index))
				if(length(index1) == 0 ){
					fastq <- paste(fastqFiles[index], collapse = " ")
					merge.command <- paste0("cat ", fastq, " > ", trim_galore_resultDir, "/", x, "_trimmed.fq.gz")
					system(merge.command)

					#必须删除单独的run的fastq文件
					delete.command <- paste("rm", fastq, sep = " ")
					system(delete.command)
				}else{
					index2 <- na.omit(index)
					if(length(index2) > 0){
						#由于部分run缺失，导致样本无法进行run合并，删除样本对应剩余的run
						fastq <- paste(fastqFiles[index2], collapse = " ")
						system(paste("rm", fastq, sep = " "))
					}
				}
			}
		})
	}
}
mergeFastq(trim_galore_resultDir="/pub6/temp/shijian/trimGalore_result", 
           sampleInfoFile="/pub6/temp/shijian/SRPSampleInfo.RData",
		   paired = TRUE)
