# RNA-seq pipeline
# 1. fastqc 是一个样本一个样本的跑，每个样本的fastq文件在一个目录里，输出也是一个样本一个目录
qualityFASTQ(inputDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/Rawdata/1", #一个样本的fastq文件所在的文件夹
             outDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/1.Fastqc/1",  #一个样本的输出文件夹
             outfilepath="/data/shijian/project/NGSCommand/Fastqc.sh",
             threads=20,
             pattern=".fq.gz$")
qualityFASTQ(inputDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/Rawdata/2",
             outDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/1.Fastqc/2",
             outfilepath="/data/shijian/project/NGSCommand/Fastqc.sh",
             threads=20,
             pattern=".fq.gz$")
qualityFASTQ(inputDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/Rawdata/3",
             outDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/1.Fastqc/3",
             outfilepath="/data/shijian/project/NGSCommand/Fastqc.sh",
             threads=20,
             pattern=".fq.gz$")
qualityFASTQ(inputDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/Rawdata/CON1",
             outDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/1.Fastqc/CON1",
             outfilepath="/data/shijian/project/NGSCommand/Fastqc.sh",
             threads=20,
             pattern=".fq.gz$")
qualityFASTQ(inputDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/Rawdata/CON2",
             outDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/1.Fastqc/CON2",
             outfilepath="/data/shijian/project/NGSCommand/Fastqc.sh",
             threads=20,
             pattern=".fq.gz$")
qualityFASTQ(inputDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/Rawdata/CON3",
             outDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/1.Fastqc/CON3",
             outfilepath="/data/shijian/project/NGSCommand/Fastqc.sh",
             threads=20,
             pattern=".fq.gz$")
# 2.multiqc 是一个样本一个样本的跑，每个样本的fastqc文件在一个目录里，输出也是一个样本一个目录
doMultiQC(inputDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/1.Fastqc/1", #一个样本的fastqc文件所在的文件夹
          outDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/1.Multiqc/1",  #一个样本的输出文件夹
          outfilepath="/data/shijian/project/NGSCommand/multiqc.sh",
          threads=20,
          pattern=".fastqc.zip$")
doMultiQC(inputDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/1.Fastqc/2",
          outDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/1.Multiqc/2",
          outfilepath="/data/shijian/project/NGSCommand/multiqc.sh",
          threads=20,
          pattern=".fastqc.zip$")
doMultiQC(inputDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/1.Fastqc/3",
          outDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/1.Multiqc/3",
          outfilepath="/data/shijian/project/NGSCommand/multiqc.sh",
          threads=20,
          pattern=".fastqc.zip$")
doMultiQC(inputDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/1.Fastqc/CON1",
          outDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/1.Multiqc/CON1",
          outfilepath="/data/shijian/project/NGSCommand/multiqc.sh",
          threads=20,
          pattern=".fastqc.zip$")
doMultiQC(inputDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/1.Fastqc/CON2",
          outDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/1.Multiqc/CON2",
          outfilepath="/data/shijian/project/NGSCommand/multiqc.sh",
          threads=20,
          pattern=".fastqc.zip$")
doMultiQC(inputDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/1.Fastqc/CON3",
          outDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/1.Multiqc/CON3",
          outfilepath="/data/shijian/project/NGSCommand/multiqc.sh",
          threads=20,
          pattern=".fastqc.zip$")

# 3. trimGalore 需要考虑接头信息，管公司要接头序列
trimGalore(fastqDir = "/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/Rawdata",
           outfilepath = "/data/shijian/project/NGSCommand/trimGalore.sh",
           outDir = "/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/3.trimGalore",pairEND = T, 
           is.phred33Encoding = TRUE, pattern = "fq.gz", q = 20, 
           length = 20, threads = 20, extraParameter = NULL)

# 4. salmon  链特异性和非链特异性的是不一样的，需要管公司要建库方式
# Salmon 通过 -l 或 --libType 参数来识别您的测序数据是链特异性的哪种类型，从而正确地将reads定位到它们来自的转录本链上。
# 测序平台/建库套装说明书（比如 Illumina TruSeq Stranded → 1st-strand）
# RSeQC 工具检测
ks.quantitation(fastqDir = "/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/3.trimGalore/trimGalore_result",
                outfilepath = "/data/shijian/project/NGSCommand/salmon.sh", 
                is.TrimGalore = T, 
                pattern = c("fastq.gz$","fq.gz$")[2], 
                type = c("kallisto", "salmon")[2], 
                index.path = "/data/shijian/refData/salmon_result/gencode.vM37.transcripts", 
                outDir = "/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/4.salmon_result", 
                pairEND=T, 
                fldMean = 100, 
                fldSD = 20, 
                is.Bias = FALSE, 
                is.BAMoutput = FALSE, 
                strand.specific = "", 
                extraParameter = NULL, 
                threads = 20)
transcript2gene.quant(quantDir = "G:/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/4.salmon_result", 
                      outDir = "G:/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/4.salmon_result/Expressionquantification", 
                      type = "salmon", is.RefGencode = T, tx2gene = NULL, countsFromAbundance = "no")

# 4.1 star比对 ht-seq量化
runStar(fastqDir = "/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/3.trimGalore/trimGalore_result",
                   index="/data/shijian/refData/mouse_reference/mouse_index/star.index",
                   outfilepath="/data/shijian/project/NGSCommand/star.sh",
                   outDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/5.star_bam")
# nohup bash /data/shijian/project/NGSCommand/star.sh > star.log 2>&1 &
runSamtools("/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/5.star_bam",
            outfilepath="/data/shijian/project/NGSCommand/samtools.sh")
# nohup bash /data/shijian/project/NGSCommand/samtools.sh > samtools.log 2>&1 &
runHTSEQ(bamDir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/5.star_bam",
                      pattern = "sortedByCoord.out.bam$",
                      strand=c("yes","no","reverse")[1],
                      gtf.path="/data/shijian/refData/mouse_reference/gencode.vM37.chr_patch_hapl_scaff.annotation.gtf",
                      out.dir="/data/shijian/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/6.counts",
                      outfilepath="/data/shijian/project/NGSCommand/htseq.sh")
# nohup bash /data/shijian/project/NGSCommand/htseq.sh > htseq.log 2>&1 &

# 5.salmon差异基因识别
load("G:/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/4.salmon_result/Expressionquantification/gene.expression.matrix.RData")
gene.expression.matrix$count[1:4,1:4]
#差异分析
exp <- as.matrix(gene.expression.matrix$count)
gene_ensemble <- stringr::str_split(rownames(exp),"\\.",simplify = T)[,1]
rownames(exp) <- gene_ensemble
gene_annotation <- getAnnotation(gtf.path="G:/ANNO_XS01KF2023120019_PM-XS01KF2023120019-07_20250721/gencode.vM37.annotation.gtf.gz",select.col=NULL,
                                   Newcolnames=NULL,control.id=NULL,organism=c("human","mouse")[2])
gene_annotation <- getAnnotation2(gene_annotation,organism=c("human","mouse")[2],geo.type=F)
  
design <- data.frame(sample=c("1","2","3","CON1","CON2","CON3"),
                     condition=c("case","case","case","control","control","control"))
rownames(design) <- design$sample
setwd("G:/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/diffGenes")
dl <- DESeq1(counts=round(as.matrix(exp)),design=design,contrast.col= "condition",contrast.level =  c("control","case"),contrast.control = "control",
              count.filter=10,cutoff.lFC = 1,cutoff.padj = 0.05,save.file = F,db=gene_annotation,names = "wangyifan",report = T)
saveRDS(dl,file="G:/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/diffGenes/diff_genes.rds")
dl <- readRDS(file="G:/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/diffGenes/diff_genes.rds")
library(openxlsx)
write.xlsx(dl1$diffSig, "G:/ANNO_XS01KF2023120019_PM-XS01KF2023120019-12/diffGenes/diff_genes.xlsx", rowNames = TRUE)
dl1$diffSig[dl1$diffSig$lab != "NO",] %>% dim()

# 6.star+htseq表达分析






