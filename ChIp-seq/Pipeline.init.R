#' @description pipeline执行文件
#' @author shi jian

#DNase-seq
convertSRA2FASTQ(outDir = "/pub6/Temp/sj/GSE77737/DNase-seq/FASTQ", 
                 outfilepath="/pub6/Temp/sj/GSE77737/NGScommands/DNase_parallelfastqdump.sh",
                 dirPath = "/pub6/Temp/sj/GSE77737/DNase-seq/SRA", 
                 pairEND = F, threads = 10)
qualityFASTQ(inputDir="/pub6/Temp/sj/GSE77737/DNase-seq/FASTQ",
             outDir="/pub6/Temp/sj/GSE77737/DNase-seq/FastQC",
             outfilepath="/pub6/Temp/sj/GSE77737/NGScommands/Dnase_Fastqc.sh",
             threads=10,
             pattern=".fastq.gz$")
doMultiQC(inputDir="/pub6/Temp/sj/GSE77737/DNase-seq/FastQC",
          outDir="/pub6/Temp/sj/GSE77737/DNase-seq/Multiqc",
          outfilepath="/pub6/Temp/sj/GSE77737/NGScommands/Dnase_multiqc.sh",
          threads=10,
          pattern=".fastqc.zip$")
trimGalore(fastqDir = "/pub6/Temp/sj/GSE77737/DNase-seq/FASTQ", 
           outfilepath = "/pub6/Temp/sj/GSE77737/NGScommands/Dnase_trimFastq.sh",
           outDir = "/pub6/Temp/sj/GSE77737/DNase-seq/", pairEND = F)
mergeFastq(trim_galore_resultDir="/pub6/Temp/sj/GSE77737/DNase-seq/trimGalore_result", 
           sampleInfoFile="/pub6/Temp/sj/GSE77737/SampleInfo.rda",
           paired = F)
Fastbowtie2(FastqDir="/pub6/Temp/sj/GSE77737/DNase-seq/trimGalore_result",
            index="/pub6/Temp/sj/GSE77737/GRCh38_noalt_as/GRCh38_noalt_as",
            outDir="/pub6/Temp/sj/GSE77737/DNase-seq/Bowtie2",
            outfilepath = "/pub6/Temp/sj/GSE77737/NGScommands/bowtie2.sh",
            pairEND = F,
            is.trimGalore = T,
            isBAM=TRUE,
            threads=24)
FastMACS2(bamDir="/pub6/Temp/sj/GSE77737/DNase-seq/Bowtie2",
          SampleInfo="/pub6/Temp/sj/GSE77737/SampleInfo.rda",
          outfilepath="/pub6/Temp/sj/GSE77737/NGScommands/Dnase_MACS2.sh",
          outDir="/pub6/Temp/sj/GSE77737/DNase-seq/MACS2",
          fileType=c("AUTO","BAM","SAM","BED")[2],
          isBroad=FALSE, 
          broad.cutoff = 0.1,
          genome=c("hs","mm","ce","dm")[1])	
FastbedGraph2BigWig(bdgDir="/pub6/Temp/sj/GSE77737/DNase-seq/MACS2",
                    fetchChromSizes = NULL,
                    outDir="/pub6/Temp/sj/GSE77737",
                    ref=c("hg19","hg38")[2],
                    outfilepath="/pub6/Temp/sj/GSE77737/NGScommands/Dase_bdg2bigwig.sh",
                    pattern=".bdg$")
Fastbam2bw(bamDir="/pub6/Temp/sj/GSE77737/DNase-seq/Bowtie2",
                   SampleInfo="/pub6/Temp/sj/GSE77737/SampleInfo.rda",
                   outfilepath="/pub6/Temp/sj/GSE77737/NGScommands/DNasebam2bw.sh",
                   outDir="/pub6/Temp/sj/GSE77737/DNase-seq/GenomeBrowser",
                   bin=25,
                   type=c("RPKM","CPM","BPM","RPGC","None")[1],
                   threads="max",
				   extraParameters=NULL,
                   pattern=".sort.bam$")					
#Chip-seq
convertSRA2FASTQ(outDir = "/pub6/Temp/sj/GSE77737/Chip-seq/FASTQ", 
                 outfilepath="/pub6/Temp/sj/GSE77737/NGScommands/parallelfastqdump.sh",
                 dirPath = "/pub6/Temp/sj/GSE77737/Chip-seq/SRA", 
                 pairEND = F, threads = 10)
qualityFASTQ(inputDir="/pub6/Temp/sj/GSE77737/Chip-seq/FASTQ",
             outDir="/pub6/Temp/sj/GSE77737/Chip-seq/FastQC",
             outfilepath="/pub6/Temp/sj/GSE77737/NGScommands/Chip_Fastqc.sh",
             threads=10,
             pattern=".fastq.gz$")
doMultiQC(inputDir="/pub6/Temp/sj/GSE77737/Chip-seq/FastQC",
          outDir="/pub6/Temp/sj/GSE77737/Chip-seq/Multiqc",
          outfilepath="/pub6/Temp/sj/GSE77737/NGScommands/Chip_multiqc.sh",
          threads=10,
          pattern=".fastqc.zip$")
trimGalore(fastqDir = "/pub6/Temp/sj/GSE77737/Chip-seq/FASTQ", 
           outfilepath = "/pub6/Temp/sj/GSE77737/NGScommands/Chip_trimFastq.sh",
           outDir = "/pub6/Temp/sj/GSE77737/Chip-seq/", pairEND = F)
mergeFastq(trim_galore_resultDir="/pub6/Temp/sj/GSE77737/Chip-seq/trimGalore_result", 
           sampleInfoFile="/pub6/Temp/sj/GSE77737/SampleInfo.rda",
           paired = F)
Fastbowtie2(FastqDir="/pub6/Temp/sj/GSE77737/Chip-seq/trimGalore_result",
            index="/pub6/Temp/sj/GSE77737/GRCh38_noalt_as/GRCh38_noalt_as",
            outDir="/pub6/Temp/sj/GSE77737/Chip-seq/Bowtie2",
            outfilepath = "/pub6/Temp/sj/GSE77737/NGScommands/Chip_bowtie2.sh",
            pairEND = F,
            is.trimGalore = T,
            isBAM=T,
            threads=24)
FastMACS2(bamDir="/pub6/Temp/sj/GSE77737/Chip-seq/Bowtie2",
          SampleInfo="/pub6/Temp/sj/GSE77737/SampleInfo.rda",
          outfilepath="/pub6/Temp/sj/GSE77737/NGScommands/Dnase_MACS2.sh",
          outDir="/pub6/Temp/sj/GSE77737/Chip-seq/MACS2",
          fileType=c("AUTO","BAM","SAM","BED")[2],
          isBroad=FALSE, 
          broad.cutoff = 0.1,
          genome=c("hs","mm","ce","dm")[1])	
FastbedGraph2BigWig(bdgDir="/pub6/Temp/sj/GSE77737/Chip-seq/MACS2",
                    fetchChromSizes = "/pub6/Temp/sj/GSE77737/hg38.chrom.sizes",
                    outDir=NULL,
                    ref=c("hg19","hg38")[2],
                    outfilepath="/pub6/Temp/sj/GSE77737/NGScommands/Chip_bdg2bigwig.sh",
                    pattern=".bdg$")
Fastbam2bw(bamDir="/pub6/Temp/sj/GSE77737/Chip-seq/Bowtie2",
                   SampleInfo="/pub6/Temp/sj/GSE77737/SampleInfo.rda",
                   outfilepath="/pub6/Temp/sj/GSE77737/NGScommands/Chipbam2bw.sh",
                   outDir="/pub6/Temp/sj/GSE77737/Chip-seq/GenomeBrowser",
                   bin=25,
                   type=c("RPKM","CPM","BPM","RPGC","None")[1],
                   threads="max",
				   extraParameters=NULL,
                   pattern=".sort.bam$")					

#RNA-seq
convertSRA2FASTQ(outDir = "/pub6/Temp/sj/GSE77737/RNA-seq/FASTQ", 
                 outfilepath="/pub6/Temp/sj/GSE77737/NGScommands/RNA_parallelfastqdump.sh",
                 dirPath = "/pub6/Temp/sj/GSE77737/RNA-seq/SRA", 
                 pairEND = T, threads = 10)

qualityFASTQ(inputDir="/pub6/Temp/sj/GSE77737/RNA-seq/FASTQ",
             outDir="/pub6/Temp/sj/GSE77737/RNA-seq/FastQC",
             outfilepath="/pub6/Temp/sj/GSE77737/NGScommands/RNA_Fastqc.sh",
             threads=10,
             pattern=".fastq.gz$")

doMultiQC(inputDir="/pub6/Temp/sj/GSE77737/RNA-seq/FastQC",
          outDir="/pub6/Temp/sj/GSE77737/RNA-seq/Multiqc",
          outfilepath="/pub6/Temp/sj/GSE77737/NGScommands/RNA_multiqc.sh",
          threads=10,
          pattern=".fastqc.zip$")
trimGalore(fastqDir = "/pub6/Temp/sj/GSE77737/RNA-seq/FASTQ", 
           outfilepath = "/pub6/Temp/sj/GSE77737/NGScommands/RNA_trimFastq.sh",
           outDir = "/pub6/Temp/sj/GSE77737/RNA-seq/", pairEND = T)

mergeFastq(trim_galore_resultDir="/pub6/Temp/sj/GSE77737/RNA-seq/trimGalore_result", 
           sampleInfoFile="/pub6/Temp/sj/GSE77737/SampleInfo.rda",
           paired = TRUE)

make.index(transcriptomeFile="/pub6/Temp/sj/GSE77737/RNA-seq/gencode.v37.transcripts.fa.gz",
           outfilepath="/pub6/Temp/sj/GSE77737/NGScommands/salmonIndex.sh",
           outDir="/pub6/Temp/sj/GSE77737/RNA-seq/salmon_result", 
           type = c("kallisto", "salmon")[2], 
           k = 31)

ks.quantitation(fastqDir="/pub6/Temp/sj/GSE77737/RNA-seq/trimGalore_result",
                outfilepath="/pub6/Temp/sj/GSE77737/NGScommands/quantitation.sh", 
                is.TrimGalore = T, 
                pattern = "fq.gz$", 
                type = c("kallisto", "salmon")[2], 
                index.path="/pub6/Temp/sj/GSE77737/RNA-seq/salmon_result/gencode.v37.transcripts", 
                outDir="/pub6/Temp/sj/GSE77737/RNA-seq/salmon_result", 
                pairEND=T, 
                fldMean = 100, 
                fldSD = 20, 
                is.Bias = FALSE, 
                is.BAMoutput = FALSE, 
                strand.specific = "", 
                extraParameter = NULL, 
                threads = 20)

transcript2gene.quant(quantDir = "/pub6/Temp/sj/GSE77737/RNA-seq/salmon_result/ExpressionQuatification", 
                      outDir = "/pub6/Temp/sj/GSE77737/RNA-seq/salmon_result/ExpressionQuatification", 
                      type = "salmon", is.RefGencode = T, tx2gene = NULL, countsFromAbundance = "no")

