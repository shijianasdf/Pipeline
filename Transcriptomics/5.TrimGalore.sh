#函数名称：trimGalore
#' @param fastqDir: fastq文件所在目录，如/IData/CancerOMICS/SkinCancer/Hugo_Cell_2017/fastq
#' @param outfilepath: 输出命令行文件地址
#' @param outDir: 结果输出目录，要求具有可写权限，如/IData/CancerOMICS/SkinCancer/Hugo_Cell_2017。程序会自动在输出目录下建立trimGalore_result文件夹，如会自动在Hugo_Cell_2017文件夹下建立trimGalore_result文件夹用于存放结果，如果存在则会写入该目录。
#' @param pairEND: 逻辑值，指明是单末端还是双末端
#' @param is.phred33Encoding: 逻辑值，指示碱基编码格式是否为phred33，FALSE则为phred64。一般使用fastqc软件能够判断出来，sanger/illumina 1.9为phred33，illumina 1.3/1.5为phred64。 
#' @param pattern: fastq文件的末尾模式，通常为"fastq.gz"、"fastq"、"fa.gz"，不能带入点号，如".fastq.gz"。
#' @param q: read质量阈值，Trim Galore软件默认要求read最小质量得分为20
#' @param length: read长度阈值，Trim Galore软件默认要求read最小长度为20。对于paired-end数据，则要求两个配对的read长度都满足该要求才会保留。
#' @param threads: 并行线程数目。一般默认为20
#' @param extraParameter 自定义扩展参数字符串，可以加入不在函数设定参数内的运行参数，默认为NULL。输入需要完全按照Trim Galore参数组合方式排列，如"--suppress_warn"

#' @return trim galore质量处理后的文件。如果为双末端测序，则输出为1_val_1.fq.gz/2_val_2.fq.gz结尾文件；如果为单末端测序，则输出为trimmed.fq.gz结尾文件

#' @note 
#' Trim Galore当前使用版本为0.4.4，为2017年3月24发布
#' 函数内部内嵌了cutadapt/trim galore/parallel软件的路径，都在/pub5/xiaoyun/BioSoftware 目录下，运行时注意检查是否存在
#### 接头处理
#' Trim Galore能够处理3种接头类型：Illumina(AGATCGGAAGAGC)、Small RNA(TGGAATTCTCGG)、Nextera(CTGTCTCTTATA)，默认为illumina类型接头。
#' 如果没有指定接头序列，Trim Galore能够通过检测前1million read的前12~13bp来自动识别接头类型。当前采用的是自动检测接头模型。

#### 环境变量冲突
#' 由于是采用R来执行linux命令，所以会出现由于环境变量问题导致程序出现问题：
#' 如[1] "Trim Galore start..."
#  sh: cat: command not found
#  /usr/bin/env: perl: No such file or directory
#  sh: rm: command not found
#' 通常需要退出R，重新进入R即可

#' example
#' 单末端
#' trimGalore(fastqDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/fastq", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single", pairEND = F)
#' 双末端
#' trimGalore(fastqDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/paired/fastq", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/paired", pairEND = T)
trimGalore <- function(fastqDir,
						outfilepath,
						outDir, 
						pairEND, 
						is.phred33Encoding = TRUE, 
						pattern = "fastq.gz", 
						q = 20, 
						length = 20, 
						threads = 20, 
						extraParameter = NULL){
		path_to_cutadapt <- "cutadapt"				
        # 使用GNU Parallel工具来实现并行
        # -j 表示并行jobs的数量
        parallel_trimGalore <- paste("parallel", "-j", threads, "trim_galore", sep = " ")
        
        # 建立输出目录
        outDir <- file.path(outDir, "trimGalore_result")
        if(!dir.exists(outDir)){
            dir.create(outDir)
        }
		
        # 主要参数设置
        if(is.phred33Encoding){
                parameterSet <- paste(" --phred33", "-q", q, "--length", length, "--no_report_file --path_to_cutadapt", path_to_cutadapt, sep = " ")
        }else{
                parameterSet <- paste(" --phred64", "-q", q, "--length", length, "--no_report_file --path_to_cutadapt", path_to_cutadapt, sep = " ")
        }
        
        # 加入扩展参数
        if(!is.null(extraParameter)){
                parameterSet <- paste(parameterSet, extraParameter, sep = " ") 
        }
        
        fq.files <- list.files(fastqDir, pattern = pattern, full.names = TRUE, recursive = TRUE)
        fq.files <- gsub(paste0(".", pattern, "$"), "", fq.files)
        fq.files <- unique(gsub(paste0("_[12]", "$"), "", fq.files))
        # 输出文件便于并行程序调用
        tmpFile <- file.path(outDir, "sampleInfo.txt")
        write.table(fq.files, file = file.path(outDir, "sampleInfo.txt"), col.names = F, row.names = F, quote = F, sep = "\n")
        
        if(!file.exists(tmpFile)){
            stop("无法建立可使用的sampleInfo文件，请检查", outDir, "是否有写入权限或fastq文件是否是", pattern, "模式")
        }
        # 双末端
        if(pairEND){
                parameterSet <- paste(" --paired", parameterSet, sep = " ")
                tg.command <- paste("cat ", tmpFile, " | ", parallel_trimGalore, parameterSet, " -o ", outDir, " {}\\_1.", pattern, " {}\\_2.", pattern, sep = "")
        }else{
                tg.command <- paste("cat ", tmpFile, " | ", parallel_trimGalore, parameterSet, " -o ", outDir, " {}\\.", pattern, sep = "")
        }
        rm.command <- paste("rm", tmpFile, sep = " ")
		commands <- c(tg.command,rm.command)
        writeLines(commands,con=outfilepath)
		return(commands)
}
trimGalore(fastqDir = "/pub6/temp/shijian/Fastq", 
           outfilepath = "/pub6/temp/shijian/NGSCommand/trimFastq.sh",
           outDir = "/pub6/temp/shijian/", pairEND = T)
./trimFastq.sh		   
