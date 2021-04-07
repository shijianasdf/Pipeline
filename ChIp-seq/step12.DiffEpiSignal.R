#' @description  merge peak的差异表观信号分析 {\BioninforR<DESeq1>}
#' @author shi jian

# 导入样本注释信息和表观信号矩阵
load("D:/Rsources/Project/StudySingleCell/data/SampleInfo.rda")
load("D:/Rsources/Project/StudySingleCell/data/mergeBedReadCount.rda")
# 矩阵对齐
intersect.sample <- intersect( colnames(mergeBedReadCount) , SampleInfo$sample_accession)
SampleInfo <- SampleInfo[SampleInfo$sample_accession %in% intersect.sample,] 
rc.matrix <- assay(mergeBedReadCount)
rownames(rc.matrix) <- rowRanges(mergeBedReadCount)$X
rc.matrix <- rc.matrix[,SampleInfo$sample_accession]

# 观察三种表观信号矩阵
table(SampleInfo$chromState)
# H3K27ac ChIP Seq H3K27me3 ChIP Seq  H3K4me1 ChIP Seq        input DNA  
# 36                 2                35                30

# H3K27ac信号矩阵和样本矩阵
pos <- which(SampleInfo$chromState == "H3K27ac ChIP Seq")
rc.ac.matrix <- rc.matrix[,pos]
H3k27AC.SampleInfo <- SampleInfo[pos,]
H3k27AC.SampleInfo$diseaseStatus[1:4] <- "normal"
H3k27AC.SampleInfo$diseaseStatus[5:36] <- "tumor"
rownames(H3k27AC.SampleInfo) <- H3k27AC.SampleInfo$sample_accession
# 疾病相对于正常样本的差异分析
library(BioinforR)
de_re <- DESeq1(counts = rc.ac.matrix,
                design = H3k27AC.SampleInfo,
                contrast.col = "diseaseStatus",
                contrast.level = c("normal","tumor"),
                contrast.control = "normal",
                count.filter=10,
                cutoff.lFC = 1,
                cutoff.padj = 0.05,
                report = F,
                save.file = T,
                names = "love")
