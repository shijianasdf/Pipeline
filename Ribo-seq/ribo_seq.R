#' @description
#'ribo-seq测序质量评估 做之前问接头序列，链特异性测序
#’安装riboWaltz
#‘ribo-seq整合RNA-seq分析
#运行riboseq处理流程
ws=/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa
strand=1
nohup bash /data/shijian/project/NGSCommand/ribo_flow.sh $ws $strand > $ws/ribo.log 2>&1 &

#将STAR比对结果toTranscriptome.sorted.bam 复制到一个指定目录 (因为riboWaltz需要)
mkdir -p /data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/3.bams/Transcriptome_bams
find /data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/3.bams -name "*.toTranscriptome.sorted.bam" | while read it
do
  echo $it
  echo "cp ${it} /data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/3.bams/Transcriptome_bams"
  cp ${it} /data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/3.bams/Transcriptome_bams
done
# read质量评估, p-site offset, 3-nt 
{
  library(riboWaltz)
  # 获取基因组注释
  #setMethod("seqinfo", "GRanges", GenomeInfoDb::seqinfo)
  annotation_dt <- create_annotation(gtfpath = "/data/shijian/refData/human_reference/gencode.v48.chr_patch_hapl_scaff.annotation.gtf")
  bamfiles <- list.files("/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/3.bams/Transcriptome_bams",
                         recursive = T,pattern = "toTranscriptome.sorted.bam$",full.names = T)
  reads_list <- bamtolist(bamfolder = dirname(bamfiles)[1],
                          transcript_align = F,annotation = annotation_dt)
  reads_list1 <- bamtolist(bamfolder = dirname(bamfiles)[1],
                          transcript_align = T,annotation = annotation_dt)
  # 1. 过滤掉重复 read（可选择仅基于 5' 端或两端）
  # filtered_reads <- duplicates_filter(data = reads_list, extremity = "both")
  # 2. 筛选特定长度的 Reads（例如保留 25 到 35 nt 之间的读段）
  filtered_reads <- length_filter(data = reads_list1,  length_filter_mode = "custom", length_range = 25:35)
  # P-site offset (P-site偏移量)
  psite_offset <- psite(filtered_reads, flanking = 6, extremity = "auto")
  reads_psite_list <- psite_info(filtered_reads, psite_offset)
  codon_coverage <- codon_coverage(reads_psite_list, annotation_dt, psite = FALSE)
  cds_coverage <- cds_coverage(reads_psite_list, annotation_dt)
  saveRDS(reads_psite_list,file="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/results/reads_psite_list.rds")
  # Read length distribution
  input_samples <- stringr::str_split(names(reads_psite_list),pattern="_",simplify = T)[,1]
  names(reads_psite_list) <- input_samples
  length_dist <- rlength_distr(reads_psite_list,
                                       sample = input_samples,
                                       multisamples = "independent",
                                       plot_style = "facet",
                                       cl = 99, colour = paletteDiscrete(values=input_samples,set = "stallion2"))
  
  pdf("/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/results/read_length_distribution.pdf")
  print(length_dist[["plot"]])
  dev.off()
  #P-sites per region
  input_samples1 <- list("Sh-2" = input_samples[1:3],
                        "Sh-3" = input_samples[4:6],
                        "control" = input_samples[7:9])
  psite_per_region <- region_psite(reads_psite_list, annotation_dt,
                                           sample = input_samples1,
                                           multisamples = "average",
                                           plot_style = "dodge",
                                           cl = 85,
                                           colour = c("#333f50", "gray70", "#39827c"))
  pdf("/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/results/psite_per_region.pdf")
  psite_per_region[["plot"]]
  dev.off()
  
  # Trinucleotide periodicity
  frames <- frame_psite(reads_psite_list, annotation_dt,
                                sample = input_samples1,
                                multisamples = "average",
                                plot_style = "facet",
                                region = "cds",
                                colour = c("#333f50", "gray70", "#39827c"))
  pdf("/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/results/trinucleotide_periodicity.pdf")
  frames[["plot"]]
  dev.off()
  
  # metaplot
  metaprofile <- metaprofile_psite(reads_psite_list, annotation_dt,
                                           sample = input_samples1,
                                           multisamples = "average",
                                           utr5l = 20, cdsl = 40, utr3l = 20,
                                           colour = c("#333f50", "gray70", "#39827c"))
  
  pdf("/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/results/metplot.pdf",width=12)
  metaprofile[["plot_Sh-2"]]
  metaprofile[["plot_Sh-3"]]
  metaprofile[["plot_control"]]
  dev.off()
}
#得到count矩阵
ribo_counts <- AggregateFeatureCount_ribo(inputDir = "/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/4.counts",
                                          pattern="CDS_counts.txt$")
gene_ensemble <- stringr::str_split(rownames(ribo_counts),"\\.",simplify = T)[,1]
rownames(ribo_counts) <- gene_ensemble
#count to TPM
geneExonLength <- getExonLength_gtf(gtfpath="/data/shijian/refData/human_reference/gencode.v48.chr_patch_hapl_scaff.annotation.gtf")
names(geneExonLength) <- stringr::str_split(names(geneExonLength),"\\.",simplify = T)[,1]
ribo_TPM <- ConvertToTPM(ribo_counts, matrices = "count", GeneExonLen = geneExonLength, log = F, logbase = NULL)
save(ribo_counts,ribo_TPM,file="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/4.counts/ribo_expr_expression.rda")
load("/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/4.counts/ribo_expr_expression.rda")

design <- data.frame(sample=colnames(ribo_counts),
                     condition=c("case","case","case","case","case","case","control","control","control"),
                     condition1=c("Sh-2","Sh-2","Sh-2","Sh-3","Sh-3","Sh-3","control","control","control"))
rownames(design) <- design$sample
# Sh-2 Sh-3 vs control
dir.create("/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh2_3",recursive = T)
setwd("/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh2_3")
sh2_3_dl <- DESeq1(counts=as.matrix(ribo_counts),design=design,contrast.col= "condition",contrast.level =  c("control","case"),contrast.control = "control",
             count.filter=10,cutoff.lFC = log2(1.2),cutoff.padj = 0.05,save.file = F, db=common.annot,names = "SPC24_Sh2_3",report = F)
saveRDS(sh2_3_dl,file="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh2_3/diff_genes.rds")
sh2_3_dl <- readRDS(file="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh2_3/diff_genes.rds")
library(openxlsx)
write.xlsx(sh2_3_dl$diffSig, "/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh2_3/diff_genes.xlsx", rowNames = TRUE)
sh2_3_dl$diffSig[sh2_3_dl$diffSig$lab != "NO",] %>% dim()

#Sh2 vs control FC=1.5
dir.create("/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh2",recursive = T)
setwd("/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh2")
dl <- DESeq1(counts=as.matrix(ribo_counts)[,c(1:3,7:9)],design=design[c(1:3,7:9),],contrast.col= "condition1",contrast.level =  c("control","Sh-2"),contrast.control = "control",
             count.filter=10,cutoff.lFC = log2(1.5),cutoff.padj = 0.05,save.file = F,db=common.annot,names = "SPC24_Sh2",report = F)
saveRDS(dl,file="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh2/diff_genes.rds")
library(openxlsx)
write.xlsx(dl$diffSig, "/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh2/diff_genes.xlsx", rowNames = TRUE)
dl$diffSig[dl$diffSig$lab != "NO",] %>% dim()
#Sh3 vs control FC=1.5
dir.create("/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh3",recursive = T)
setwd("/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh3")
dl1 <- DESeq1(counts=as.matrix(ribo_counts)[,c(4:6,7:9)],design=design[c(4:6,7:9),],contrast.col= "condition1",contrast.level =  c("control","Sh-3"),contrast.control = "control",
             count.filter=10,cutoff.lFC = log2(1.5),cutoff.padj = 0.05,save.file = F,db=common.annot,names = "SPC24_Sh3",report = F)
saveRDS(dl1,file="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh3/diff_genes.rds")
library(openxlsx)
write.xlsx(dl1$diffSig, "/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh3/diff_genes.xlsx", rowNames = TRUE)
dl1$diffSig[dl1$diffSig$lab != "NO",] %>% dim()
#Sh2 vs control FC=1.2
dir.create("/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh2_FC1.2",recursive = T)
setwd("/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh2_FC1.2")
ribo_dl <- DESeq1(counts=as.matrix(ribo_counts)[,c(1:3,7:9)],design=design[c(1:3,7:9),],contrast.col= "condition1",contrast.level =  c("control","Sh-2"),contrast.control = "control",
             count.filter=10,cutoff.lFC = log2(1.2),cutoff.padj = 0.05,save.file = F,db=common.annot,names = "SPC24_Sh2",report = F)
saveRDS(ribo_dl,file="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh2_FC1.2/diff_genes.rds")
library(openxlsx)
write.xlsx(ribo_dl$diffSig, "/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh2_FC1.2/diff_genes.xlsx", rowNames = TRUE)
ribo_dl$diffSig[ribo_dl$diffSig$lab != "NO",] %>% dim()
#Sh3 vs control  FC=1.2
dir.create("/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh3_FC1.2",recursive = T)
setwd("/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh3_FC1.2")
ribo_dl1 <- DESeq1(counts=as.matrix(ribo_counts)[,c(4:6,7:9)],design=design[c(4:6,7:9),],contrast.col= "condition1",contrast.level =  c("control","Sh-3"),contrast.control = "control",
              count.filter=10,cutoff.lFC = log2(1.2),cutoff.padj = 0.05,save.file = F,db=common.annot,names = "SPC24_Sh3",report = F)
saveRDS(ribo_dl1,file="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh3_FC1.2/diff_genes.rds")
library(openxlsx)
write.xlsx(ribo_dl1$diffSig, "/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh3_FC1.2/diff_genes.xlsx", rowNames = TRUE)
ribo_dl1$diffSig[ribo_dl1$diffSig$lab != "NO",] %>% dim()
# RNA-seq和Ribo-seq整合
# 相关性计算
load("/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/4.counts/ribo_expr_expression.rda")
load("/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0155_RNA_9hsa/5.featurecounts_CDS/RNA_expr_expression.rda")
#过滤低表达基因
log_rna_tpm <- log2eset(rna_TPM)
log_ribo_tpm <- log2eset(ribo_TPM)
keep1 <- apply(log_ribo_tpm,1,var) > 0   
keep2 <- apply(log_rna_tpm,1,var) > 0 
intersect_genes <- intersect(rownames(log_rna_tpm[keep2,]),rownames(log_ribo_tpm[keep1,]))
cor_res <- MakeCorrelationMatrix(x=t(log_rna_tpm[intersect_genes,]), y= t(log_ribo_tpm[intersect_genes,]), size= 500, adjust= "BH",use=c("everything", "all.obs", "complete.obs", "na.or.complete","pairwise.complete.obs")[1],
                                 method = c("pearson", "kendall", "spearman")[1])
save(cor_res,file="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/cor_res.rda")
load(file="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/cor_res.rda")
cor_res$r[is.na(cor_res$r)] <- 0
heatmap_res <- FastComplexHeatmap(cor_res$r,df1=NULL,df1.col="gene",show_rowcol1=NULL,
                                   df2=NULL,df2.col="cell",show_rowcol2=NULL,scale=T,
                                   legend_name="correlation",cluster_columns=F,cluster_rows=F,
                                   show_column_names=F,show_row_names=F,
                                   row_split=NULL,column_split=NULL,
                                   palletes=c("horizon","blueYellow","horizonExtra","sambaNight","solarExtra","whitePurple","whiteBlue","whiteRed","comet","greenBlue","beach","coolwarm","fireworks","greyMagma","fireworks2","purpleOrange", #连续性颜色
                                              "stallion","stallion2","calm","kelly","bear","ironMan","circus","paired","grove","summerNight","zissou","darjeeling","rushmore","captain","yanru")[5], #离散型颜色
                                   color=NULL,color.type=c("continues","discrete")[1],
                                   path="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/results",
                                   picture.name="ribo_rna_cor",picture.type=c("png","pdf")[1])

#sh2 差异基因overlap 功能注释
rna_diff <- readRDS(file="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0155_RNA_9hsa/5.featureCount_Des/Sh2/diff_genes.rds")
ribo_diff <- readRDS(file="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh2_FC1.2/diff_genes.rds")
up_list <- list(rna_seq = rownames(rna_diff$diffSig)[rna_diff$diffSig$lab == "UP"],
                ribo_seq= rownames(ribo_diff$diffSig)[ribo_diff$diffSig$lab == "UP"])
up_p <- venn_plot(up_list, title = "upgenes venn diagram", colors = c("#4575b4", "#fee090", "#d73027"), # 默认渐变色 palettes(category = "continue", palette = "rdbu", alpha = 1, counts = 50, show_col = TRUE, show_message = FALSE)[c(5,6,8)]
               label_type = "both",      # "count", "percent", "both", "none"
               save_path = NULL,         # 文件路径，如 "venn.pdf" 或 "venn.png"
               show_edge = F,         # 是否显示边框
               edge_lty = "solid",       # 边框线条类型
               alpha = 0.2,              # 填充透明度
               label_size = 6,           # 内部数字大小
               set_size = 3              # 组名文字大小
          )
down_list <- list(rna_seq = rownames(rna_diff$diffSig)[rna_diff$diffSig$lab == "DOWN"],
                  ribo_seq= rownames(ribo_diff$diffSig)[ribo_diff$diffSig$lab == "DOWN"])
down_p <- venn_plot(down_list, title = "downgenes venn diagram", colors = c("#4575b4", "#fee090", "#d73027"), # 默认渐变色 palettes(category = "continue", palette = "rdbu", alpha = 1, counts = 50, show_col = TRUE, show_message = FALSE)[c(5,6,8)]
                  label_type = "both",      # "count", "percent", "both", "none"
                  save_path = NULL,         # 文件路径，如 "venn.pdf" 或 "venn.png"
                  show_edge = F,         # 是否显示边框
                  edge_lty = "solid",       # 边框线条类型
                  alpha = 0.2,              # 填充透明度
                  label_size = 6,           # 内部数字大小
                  set_size = 3              # 组名文字大小
)
ggsave(up_p,filename="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/results/upgenes_overlap.pdf",width=10)
ggsave(down_p,filename="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/results/downgenes_overlap.pdf",width=10)

#sh2差异基因，功能注释
setwd("/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa")
#up_genes <- intersect(up_list$rna_seq,up_list$ribo_seq)
up_genes <- up_list$ribo_seq
geneList <- log_ribo_tpm[,1]
up_kegg <- FastKEGG(genes=up_genes,geneList = geneList,default.universe = F,organism = c('hsa',"mmu")[1],db=common.annot,
                    pvalueCutoff = list(enrichKEGG = 0.05,enrichMKEGG = 0.05),qvalueCutoff = 0.05,
                    cnet.showCategory = 10,verbose = T,save.path = "up_KEGG",names = "6.EnrichmentRes")
up_go <- FastGO(genes=up_genes,geneList=geneList,organism = c("human","mouse")[1],default.universe = F,
                classlevel = 2:2,OrgDb  = NULL,  #c("org.Mm.eg.db","org.Hs.eg.db")
                keyType = NULL,pAdjustMethod = "BH",pvalueCutoff = 0.05,
                qvalueCutoff  = 0.05,cnet.showCategory = 5,verbose = TRUE,save.path = "up_GO",names = "6.EnrichmentRes")
up_hallmark_anno <- FastMS(genes = up_genes, geneList = geneList, default.universe = F,run_gsea = F,
                           pAdjustMethod = "BH", db = common.annot, organism = c("human","mouse")[1],
                           pvalueCutoff = list(enrichMS=0.05,
                                               gseMS = 0.05), qvalueCutoff  = 0.05,
                           cnet.showCategory = 5, verbose = TRUE, save.path = "up_HALLMARK",names = "6.EnrichmentRes")

#down_genes <- intersect(down_list$rna_seq,down_list$ribo_seq)
down_genes <- down_list$ribo_seq
down_kegg <- FastKEGG(genes=down_genes,geneList = geneList,default.universe = F,organism = c('hsa',"mmu")[1],db=common.annot,
                      pvalueCutoff = list(enrichKEGG = 0.05,enrichMKEGG = 0.05),qvalueCutoff = 0.05,
                      cnet.showCategory = 5,verbose = T,save.path = "down_KEGG",names = "6.EnrichmentRes")
down_go <- FastGO(genes=down_genes,geneList=geneList,organism = c("human","mouse")[1],default.universe = F,
                  classlevel = 2:2,OrgDb  = NULL,  #c("org.Mm.eg.db","org.Hs.eg.db")
                  keyType = NULL,pAdjustMethod = "BH",pvalueCutoff = 0.05,
                  qvalueCutoff  = 0.05,cnet.showCategory = 5,verbose = TRUE,save.path = "down_GO",names = "6.EnrichmentRes")
down_hallmark_anno <- FastMS(genes = down_genes, geneList = geneList, default.universe = F,run_gsea = F,
                             pAdjustMethod = "BH", db = common.annot, organism = c("human","mouse")[1],
                             pvalueCutoff = list(enrichMS=0.05,
                                                 gseMS = 0.05), qvalueCutoff  = 0.05,
                             cnet.showCategory = 5, verbose = TRUE, save.path = "down_HALLMARK",names = "6.EnrichmentRes")
save(up_kegg,up_go,up_hallmark_anno,down_go,down_hallmark_anno,file="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/6.EnrichmentRes/enrichRes1.rda")
load(file="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/6.EnrichmentRes/enrichRes1.rda")

#sh3 差异基因overlap 
sh3_rna_diff <- readRDS(file="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0155_RNA_9hsa/5.featureCount_Des/Sh3/diff_genes.rds")
sh3_ribo_diff <- readRDS(file="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh3_FC1.2/diff_genes.rds")
up_list <- list(rna_seq = rownames(sh3_rna_diff$diffSig)[sh3_rna_diff$diffSig$lab == "UP"],
                ribo_seq= rownames(sh3_ribo_diff$diffSig)[sh3_ribo_diff$diffSig$lab == "UP"])
up_p <- venn_plot(up_list, title = "sh3 upgenes venn diagram", colors = c("#4575b4", "#fee090", "#d73027"), # 默认渐变色 palettes(category = "continue", palette = "rdbu", alpha = 1, counts = 50, show_col = TRUE, show_message = FALSE)[c(5,6,8)]
                  label_type = "both",      # "count", "percent", "both", "none"
                  save_path = NULL,         # 文件路径，如 "venn.pdf" 或 "venn.png"
                  show_edge = F,         # 是否显示边框
                  edge_lty = "solid",       # 边框线条类型
                  alpha = 0.2,              # 填充透明度
                  label_size = 6,           # 内部数字大小
                  set_size = 3              # 组名文字大小
)
down_list <- list(rna_seq = rownames(sh3_rna_diff$diffSig)[sh3_rna_diff$diffSig$lab == "DOWN"],
                  ribo_seq= rownames(sh3_ribo_diff$diffSig)[sh3_ribo_diff$diffSig$lab == "DOWN"])
down_p <- venn_plot(down_list, title = "sh3 downgenes venn diagram", colors = c("#4575b4", "#fee090", "#d73027"), # 默认渐变色 palettes(category = "continue", palette = "rdbu", alpha = 1, counts = 50, show_col = TRUE, show_message = FALSE)[c(5,6,8)]
                    label_type = "both",      # "count", "percent", "both", "none"
                    save_path = NULL,         # 文件路径，如 "venn.pdf" 或 "venn.png"
                    show_edge = F,         # 是否显示边框
                    edge_lty = "solid",       # 边框线条类型
                    alpha = 0.2,              # 填充透明度
                    label_size = 6,           # 内部数字大小
                    set_size = 3              # 组名文字大小
)
ggsave(up_p,filename="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/results/sh3_upgenes_overlap.pdf",width=10)
ggsave(down_p,filename="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/results/sh3_downgenes_overlap.pdf",width=10)

# sh-2和sh-3敲低方式差异基因overlap
{
  sh3_ribo_diff <- readRDS(file="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh3_FC1.2/diff_genes.rds")
  sh2_ribo_diff <- readRDS(file="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/5.Des/Sh2_FC1.2/diff_genes.rds")
  ribo_intersect_up <- list(sh3_ribo=rownames(sh3_ribo_diff$diffSig)[sh3_ribo_diff$diffSig$lab == "UP"],
                            sh2_ribo=rownames(sh2_ribo_diff$diffSig)[sh2_ribo_diff$diffSig$lab == "UP"])
  up_p <- venn_plot(ribo_intersect_up, title = "sh3-sh2 ribo upgenes venn diagram", colors = c("#4575b4", "#fee090", "#d73027"), # 默认渐变色 palettes(category = "continue", palette = "rdbu", alpha = 1, counts = 50, show_col = TRUE, show_message = FALSE)[c(5,6,8)]
                    label_type = "both",      # "count", "percent", "both", "none"
                    save_path = NULL,         # 文件路径，如 "venn.pdf" 或 "venn.png"
                    show_edge = F,         # 是否显示边框
                    edge_lty = "solid",       # 边框线条类型
                    alpha = 0.2,              # 填充透明度
                    label_size = 6,           # 内部数字大小
                    set_size = 3              # 组名文字大小
  )
  ribo_intersect_down <- list(sh3_ribo=rownames(sh3_ribo_diff$diffSig)[sh3_ribo_diff$diffSig$lab == "DOWN"],
                            sh2_ribo=rownames(sh2_ribo_diff$diffSig)[sh2_ribo_diff$diffSig$lab == "DOWN"])
  down_p <- venn_plot(ribo_intersect_down, title = "sh3-sh2 ribo downgenes venn diagram", colors = c("#4575b4", "#fee090", "#d73027"), # 默认渐变色 palettes(category = "continue", palette = "rdbu", alpha = 1, counts = 50, show_col = TRUE, show_message = FALSE)[c(5,6,8)]
                      label_type = "both",      # "count", "percent", "both", "none"
                      save_path = NULL,         # 文件路径，如 "venn.pdf" 或 "venn.png"
                      show_edge = F,         # 是否显示边框
                      edge_lty = "solid",       # 边框线条类型
                      alpha = 0.2,              # 填充透明度
                      label_size = 6,           # 内部数字大小
                      set_size = 3              # 组名文字大小
  )
  ggsave(up_p,filename="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/results/sh3-sh2_ribo_upgenes_overlap.pdf",width=10)
  ggsave(down_p,filename="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/results/sh3-sh2_ribo_downgenes_overlap.pdf",width=10)
}

#sh3差异基因功能注释
setwd("/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa")
#up_genes <- intersect(up_list$rna_seq,up_list$ribo_seq)
up_genes <- up_list$ribo_seq
geneList <- log_ribo_tpm[,1]
up_kegg <- FastKEGG(genes=up_genes,geneList = geneList,default.universe = F,organism = c('hsa',"mmu")[1],db=common.annot,
                    pvalueCutoff = list(enrichKEGG = 0.05,enrichMKEGG = 0.05),qvalueCutoff = 0.05,
                    cnet.showCategory = 10,verbose = T,save.path = "sh3_up_KEGG",names = "6.EnrichmentRes")
up_go <- FastGO(genes=up_genes,geneList=geneList,organism = c("human","mouse")[1],default.universe = F,
                classlevel = 2:2,OrgDb  = NULL,  #c("org.Mm.eg.db","org.Hs.eg.db")
                keyType = NULL,pAdjustMethod = "BH",pvalueCutoff = 0.05,
                qvalueCutoff  = 0.05,cnet.showCategory = 5,verbose = TRUE,save.path = "sh3_up_GO",names = "6.EnrichmentRes")
up_hallmark_anno <- FastMS(genes = up_genes, geneList = geneList, default.universe = F,run_gsea = F,
                           pAdjustMethod = "BH", db = common.annot, organism = c("human","mouse")[1],
                           pvalueCutoff = list(enrichMS=0.05,
                                               gseMS = 0.05), qvalueCutoff  = 0.05,
                           cnet.showCategory = 5, verbose = TRUE, save.path = "sh3_up_HALLMARK",names = "6.EnrichmentRes")
#down_genes <- intersect(down_list$rna_seq,down_list$ribo_seq)
down_genes <- down_list$ribo_seq
down_kegg <- FastKEGG(genes=down_genes,geneList = geneList,default.universe = F,organism = c('hsa',"mmu")[1],db=common.annot,
                      pvalueCutoff = list(enrichKEGG = 0.05,enrichMKEGG = 0.05),qvalueCutoff = 0.05,
                      cnet.showCategory = 5,verbose = T,save.path = "sh3_down_KEGG",names = "6.EnrichmentRes")
down_go <- FastGO(genes=down_genes,geneList=geneList,organism = c("human","mouse")[1],default.universe = F,
                  classlevel = 2:2,OrgDb  = NULL,  #c("org.Mm.eg.db","org.Hs.eg.db")
                  keyType = NULL,pAdjustMethod = "BH",pvalueCutoff = 0.05,
                  qvalueCutoff  = 0.05,cnet.showCategory = 5,verbose = TRUE,save.path = "sh3_down_GO",names = "6.EnrichmentRes")
down_hallmark_anno <- FastMS(genes = down_genes, geneList = geneList, default.universe = F,run_gsea = F,
                             pAdjustMethod = "BH", db = common.annot, organism = c("human","mouse")[1],
                             pvalueCutoff = list(enrichMS=0.05,
                                                 gseMS = 0.05), qvalueCutoff  = 0.05,
                             cnet.showCategory = 5, verbose = TRUE, save.path = "sh3_down_HALLMARK",names = "6.EnrichmentRes")
save(up_kegg,up_go,up_hallmark_anno,down_kegg,down_go,down_hallmark_anno,file="/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/6.EnrichmentRes/sh3_enrichRes1.rda")


#差异翻译效率计算
{
  library(xtail)
  library(ggrepel)
  #导入ribo 和 RNA counts矩阵
  load("/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/4.counts/ribo_expr_expression.rda")
  load("/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0155_RNA_9hsa/5.featurecounts_CDS/RNA_expr_expression.rda")
  dim(ribo_counts)
  dim(rna_counts)
  intersect_genes <- intersect(rownames(ribo_counts),rownames(rna_counts))
  xtail_ribo_counts <- ribo_counts[intersect_genes,]
  xtail_rna_counts <- rna_counts[intersect_genes,]
  condition <- c("control","control","control","Sh_2","Sh_2","Sh_2")
  condition1 <- c("control","control","control","Sh_3","Sh_3","Sh_3")
  #
  sh2.results <- xtail(xtail_rna_counts[,c(7:9,1:3)],xtail_ribo_counts[,c(7:9,1:3)],condition,bins=1000)
  sh2.tab <- resultsTable(sh2.results,sort.by = "pvalue_final", log2FCs = TRUE,log2Rs = TRUE)
  sh3.results <- xtail(xtail_rna_counts[,c(7:9,4:6)],xtail_ribo_counts[,c(7:9,4:6)],condition1,bins=1000)
  sh3.tab <- resultsTable(sh3.results,sort.by = "pvalue_final", log2FCs = TRUE,log2Rs = TRUE)
  sh2.tab$SYMBOL <- convert(rownames(sh2.tab),fromtype="ENSEMBL",totype="SYMBOL")
  sh3.tab$SYMBOL <- convert(rownames(sh3.tab),fromtype="ENSEMBL",totype="SYMBOL")
  
  pdf("/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/results/sh2_xtail.pdf")
  sh2.tab <- sh2.tab[complete.cases(as.data.frame(sh2.tab)),]
  categories <- c("transcription only" = "#0D5FFF",
                  "translation only" = "#FF2E06",
                  "homodirectional" = "#75E805",
                  "opposite change" = "#FFDE13",
                  "stable" = "gray90")
  .add_category_label <- function(data, x_cat, y_cat, 
                                  cutoff, categories,add){
              data[,add] <- categories[5L]
              x_only <- abs(data[,x_cat]) >= cutoff &
                abs(data[,y_cat]) < cutoff
              data[,add][x_only] <- categories[1L]
              y_only <- abs(data[,y_cat]) >= cutoff &
                abs(data[,x_cat]) < cutoff
              data[,add][y_only] <- categories[2L]
              x_and_y <- abs(data[,x_cat]) >= cutoff &
                abs(data[,y_cat]) >= cutoff
              data[,add][x_and_y] <- categories[3L]
              x_opposite_y <- abs(data[,x_cat]) >= cutoff &
                abs(data[,y_cat]) >= cutoff &
                sign(data[,x_cat]) != sign(data[,y_cat])
              data[,add][x_opposite_y] <- categories[4L]
              data
  }
  sh2.tab <- .add_category_label(sh2.tab,x_cat="mRNA_log2FC",y_cat="RPF_log2FC",
                                   cutoff=log2(1.5),categories=names(categories), add="Category")
  cond <- paste0(conditions(sh2.results),"_log2TE")
  categories1 <- c("condition1" = "#0D5FFF",
                  "condition2" = "#FF2E06",
                  "homodirectional" = "#75E805",
                  "opposite change" = "#FFDE13",
                  "stable" = "gray90")
  names(categories1)[seq_len(2L)] <- paste0(conditions(sh2.results), " only")
  sh2.tab <- .add_category_label(data=sh2.tab,x_cat=cond[1L],y_cat=cond[2L],
                                 cutoff=log2(1.5),categories=names(categories1), add="Category1")
  sh2.tab <- as.data.frame(sh2.tab) %>% 
                     dplyr::mutate(lab=if_else(abs(log2FC_TE_final) > log2(1.2) & pvalue_final < 0.05,
                                               if_else(log2FC_TE_final > log2(1.2),"UP","DOWN"),"NO"))
  # 1. 筛选最显著的 10 个上调 + 10 个下调基因
  top_up <- sh2.tab %>%
    filter(pvalue_final < 0.05, log2FC_TE_final > log2(1.2)) %>%   # 显著上调
    arrange(pvalue_final) %>%                             # 按 padj 从小到大排序
    head(10)
  top_down <- sh2.tab %>%
    filter(pvalue_final < 0.05, log2FC_TE_final < log2(1/1.2)) %>%   # 显著下调
    arrange(pvalue_final) %>%
    head(10)
  sh2_top_genes <- bind_rows(top_up, top_down)
  
  p1 <- ggplot(data=as.data.frame(sh2.tab),mapping=aes(x = mRNA_log2FC, y = RPF_log2FC,color=Category)) +
    geom_hline(yintercept = log2(1.5), colour = "gray") +
    geom_hline(yintercept = -log2(1.5), colour = "gray") +
    geom_vline(xintercept = log2(1.5), colour = "gray") +
    geom_vline(xintercept = -log2(1.5), colour = "gray") +
    geom_point(shape = 19L, alpha = 0.65) +
    #geom_text_repel(data=sh2_top_genes,mapping=aes(x = mRNA_log2FC,y = RPF_log2FC,label=SYMBOL)) +
    scale_x_continuous(limits = NULL) +
    scale_y_continuous(limits = NULL) +
    scale_colour_manual(values = categories) +
    theme_bw()
  p2 <- ggplot(data=as.data.frame(sh2.tab),mapping=aes(x = Sh_2_log2TE, y = control_log2TE,color=Category1)) +
    geom_hline(yintercept = log2(1.5), colour = "gray") +
    geom_hline(yintercept = -log2(1.5), colour = "gray") +
    geom_vline(xintercept = log2(1.5), colour = "gray") +
    geom_vline(xintercept = -log2(1.5), colour = "gray") +
    geom_point(shape = 19L, alpha = 0.65) +
    #geom_text_repel(data=sh2_top_genes,aes(x = Sh.3_log2TE,y = control_log2TE,label=SYMBOL)) +
    scale_x_continuous(limits = NULL) +
    scale_y_continuous(limits = NULL) +
    scale_colour_manual(values = categories1) +
    theme_bw()
  p3 <- ggplot(data = as.data.frame(sh2.tab), mapping=aes(x = log2FC_TE_final,y = -log10(pvalue_final), color=lab)) +
             geom_point(alpha = 0.65, shape = 19L) +
             geom_text_repel(data=sh2_top_genes,mapping=aes(x = log2FC_TE_final,y = -log10(pvalue_final),label=SYMBOL)) +
             scale_colour_manual(values=paletteDiscrete(values=c("UP","NO","DOWN"),set = "ribowaltz")) +
             theme_bw()
  print(p1)
  print(p2)
  print(p3)
  dev.off()
  
  pdf("/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/results/sh3_xtail.pdf")
  sh3.tab <- sh3.tab[complete.cases(as.data.frame(sh3.tab)),]
  categories <- c("transcription only" = "#0D5FFF",
                  "translation only" = "#FF2E06",
                  "homodirectional" = "#75E805",
                  "opposite change" = "#FFDE13",
                  "stable" = "gray90")
  .add_category_label <- function(data, x_cat, y_cat, 
                                  cutoff, categories,add){
    data[,add] <- categories[5L]
    x_only <- abs(data[,x_cat]) >= cutoff &
      abs(data[,y_cat]) < cutoff
    data[,add][x_only] <- categories[1L]
    y_only <- abs(data[,y_cat]) >= cutoff &
      abs(data[,x_cat]) < cutoff
    data[,add][y_only] <- categories[2L]
    x_and_y <- abs(data[,x_cat]) >= cutoff &
      abs(data[,y_cat]) >= cutoff
    data[,add][x_and_y] <- categories[3L]
    x_opposite_y <- abs(data[,x_cat]) >= cutoff &
      abs(data[,y_cat]) >= cutoff &
      sign(data[,x_cat]) != sign(data[,y_cat])
    data[,add][x_opposite_y] <- categories[4L]
    data
  }
  sh3.tab <- .add_category_label(sh3.tab,x_cat="mRNA_log2FC",y_cat="RPF_log2FC",
                                 cutoff=log2(1.5),categories=names(categories), add="Category")
  cond <- paste0(conditions(sh3.results),"_log2TE")
  categories1 <- c("condition1" = "#0D5FFF",
                   "condition2" = "#FF2E06",
                   "homodirectional" = "#75E805",
                   "opposite change" = "#FFDE13",
                   "stable" = "gray90")
  names(categories1)[seq_len(2L)] <- paste0(conditions(sh3.results), " only")
  sh3.tab <- .add_category_label(data=sh3.tab,x_cat=cond[1L],y_cat=cond[2L],
                                 cutoff=log2(1.5),categories=names(categories1), add="Category1")
  sh3.tab <- as.data.frame(sh3.tab) %>% 
                    dplyr::mutate(lab=if_else(abs(log2FC_TE_final) > log2(1.2) & pvalue_final < 0.05,
                                              if_else(log2FC_TE_final > log2(1.2),"UP","DOWN"),"NO"))
  # 1. 筛选最显著的 10 个上调 + 10 个下调基因
  top_up <- sh3.tab %>%
    filter(pvalue_final < 0.05, log2FC_TE_final > log2(1.2)) %>%   # 显著上调
    arrange(pvalue_final) %>%                             # 按 padj 从小到大排序
    head(10)
  top_down <- sh3.tab %>%
    filter(pvalue_final < 0.05, log2FC_TE_final < log2(1/1.2)) %>%   # 显著下调
    arrange(pvalue_final) %>%
    head(10)
  sh3_top_genes <- bind_rows(top_up, top_down)
  
  p1 <- ggplot(data=as.data.frame(sh3.tab),mapping=aes(x = mRNA_log2FC, y = RPF_log2FC,color=Category)) +
          geom_hline(yintercept = log2(1.5), colour = "gray") +
          geom_hline(yintercept = -log2(1.5), colour = "gray") +
          geom_vline(xintercept = log2(1.5), colour = "gray") +
          geom_vline(xintercept = -log2(1.5), colour = "gray") +
          geom_point(shape = 19L, alpha = 0.65) +
          #geom_text_repel(sh3_top_genes,aes(x = mRNA_log2FC,y = RPF_log2FC,label=SYMBOL)) +
          scale_x_continuous(limits = NULL) +
          scale_y_continuous(limits = NULL) +
          scale_colour_manual(values = categories) +
          theme_bw()
  p2 <- ggplot(data=as.data.frame(sh3.tab),mapping=aes(x = Sh_3_log2TE, y = control_log2TE,color=Category1)) +
          geom_hline(yintercept = log2(1.5), colour = "gray") +
          geom_hline(yintercept = -log2(1.5), colour = "gray") +
          geom_vline(xintercept = log2(1.5), colour = "gray") +
          geom_vline(xintercept = -log2(1.5), colour = "gray") +
          geom_point(shape = 19L, alpha = 0.65) +
          #geom_text_repel(sh3_top_genes,aes(x = Sh.3_log2TE,y = control_log2TE,label=SYMBOL))+
          scale_x_continuous(limits = NULL) +
          scale_y_continuous(limits = NULL) +
          scale_colour_manual(values = categories1) +
          theme_bw()
  p3 <- ggplot(data = as.data.frame(sh3.tab), mapping=aes(x = log2FC_TE_final,y = -log10(pvalue_final), color=lab)) +
          geom_point(alpha = 0.65, shape = 19L) +
          geom_text_repel(data=sh3_top_genes,mapping=aes(x = log2FC_TE_final,y = -log10(pvalue_final),label=SYMBOL)) +
          scale_colour_manual(values=paletteDiscrete(values=c("UP","NO","DOWN"),set = "ribowaltz")) +
          theme_bw()
  print(p1)
  print(p2)
  print(p3)
  dev.off()
  
  library(openxlsx)
  #sh2.tab$SYMBOL <- convert(rownames(sh2.tab),fromtype="ENSEMBL",totype="SYMBOL")
  #sh3.tab$SYMBOL <- convert(rownames(sh3.tab),fromtype="ENSEMBL",totype="SYMBOL")
  write.xlsx(sh2.tab, "/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/results/sh2_diff_TE_genes.xlsx", rowNames = TRUE)
  write.xlsx(sh3.tab, "/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/results/sh3_diff_TE_genes.xlsx", rowNames = TRUE)
  
}





