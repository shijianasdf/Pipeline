# cut&tag分析流程
# @author shi jian
# @date 2025.6.23
# 

#检查阿里云下载数据是不是完全，linux通过md5sum命令行验证
#cd /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25 
#md5sum -c md5.txt

##== linux 命令 ==##
#创建conda cut&tag环境
conda create -y --name epigenetics 
projPath="/data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25"
histName="12"
mkdir -p ${projPath}/fastqFileQC/${histName}

# 1.fastqc检测数据质量结果
fastqc -o ${projPath}/fastqFileQC/${histName} -f fastq ${projPath}/Rawdata/${histName}/${histName}_R1.fq.gz ${projPath}/Rawdata/${histName}/${histName}_R2.fq.gz

# 2.multiqc结果
mkdir -p ${projPath}/multiFileQC/${histName}
multiqc ${projPath}/fastqFileQC/${histName}/${histName}_R1_fastqc.zip ${projPath}/fastqFileQC/${histName}/${histName}_R2_fastqc.zip  -o ${projPath}/multiFileQC/${histName}

# 3.cutadapt去接头（Illuminal双端测序，P5 adapter P7 adaptor这两个不用去除，index也不用管，
#                主要是read1-SP和read2-SP需要去除，公司给的测序流程查到 read1-SP和read2-SP的序列，然后去除）
mkdir -p ${projPath}/cutadapt/${histName}
cutadapt \
  -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG \
  -A GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
  -o ${projPath}/cutadapt/${histName}/${histName}_trimmed_R1.fastq.gz \
  -p ${projPath}/cutadapt/${histName}/${histName}_trimmed_R2.fastq.gz \
  -j 10 \
  ${projPath}/Rawdata/${histName}/${histName}_R1.fq.gz ${projPath}/Rawdata/${histName}/${histName}_R2.fq.gz

# 4.bowtie2比对（no spike-in）
cores=12
ref="/home/shijian/refData/bowtie2Index/GRCm39.genome.mm"
mkdir -p ${projPath}/alignment/sam/bowtie2_summary
mkdir -p ${projPath}/alignment/bam
mkdir -p ${projPath}/alignment/bed
mkdir -p ${projPath}/alignment/bedgraph
mkdir -p ${projPath}/alignment/bigwig
## Build the bowtie2 reference genome index if needed:
mkdir -p ~/refData/bowtie2Index
bowtie2-build ~/refData/GRCm39.genome.fa.gz ~/refData/bowtie2Index/GRCm39.genome.mm
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${ref} -1 ${projPath}/cutadapt/${histName}/${histName}_trimmed_R1.fastq.gz -2 ${projPath}/cutadapt/${histName}/${histName}_trimmed_R2.fastq.gz -S ${projPath}/alignment/sam/${histName}_bowtie2.sam &> ${projPath}/alignment/sam/bowtie2_summary/${histName}_bowtie2.txt
# 5.对 spike-in 基因组进行校准，以进行 spike-in 校准（可选 / 推荐）
spikeInRef="/home/shijian/refData/bowtie2Index/Ecoli"
chromSize="/home/shijian/refData/mm39.chrom.sizes"
bowtie2-build ~/refData/GCF_000005845.2_ASM584v2_genomic.fa ~/refData/bowtie2Index/Ecoli
spikeInRef="/home/shijian/refData/spike_in_genome/spike_in.fa"
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --no-overlap --no-dovetail --phred33 -I 10 -X 700 -p ${cores} -x ${spikeInRef} -1 ${projPath}/cutadapt/${histName}/${histName}_trimmed_R1.fastq.gz -2 ${projPath}/cutadapt/${histName}/${histName}_trimmed_R2.fastq.gz -S $projPath/alignment/sam/${histName}_bowtie2_spikeIn.sam &> $projPath/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.txt
seqDepthDouble=$(samtools view -F 0x04 $projPath/alignment/sam/${histName}_bowtie2_spikeIn.sam | wc -l)
seqDepth=$((seqDepthDouble/2))
echo $seqDepth > $projPath/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.seqDepth
#spike_inrc=$(samtools view -c -F 4 -q 30 $projPath/alignment/sam/${histName}_bowtie2_spikeIn.sam)
#echo $spike_inrc >  $projPath/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.mappedread
# 6.去除重复 reads（可选 / 需要） 不建议去重
##== linux 命令 ==##
## depending on how you load picard and your server environment, the picardCMD can be different. Adjust accordingly.
#picardCMD="java -jar picard.jar"
#mkdir -p $projPath/alignment/rmDuplicate/picard_summary
## Sort by coordinate
## 按照坐标排序
#$picardCMD SortSam \
#I=$projPath/alignment/sam/${histName}_bowtie2.sam \
#O=$projPath/alignment/sam/${histName}_bowtie2.sorted.sam \
#SORT_ORDER=coordinate
## mark duplicates
## 标记重复
#$picardCMD MarkDuplicates \
#I=$projPath/alignment/sam/${histName}_bowtie2.sorted.sam \
#O=$projPath/alignment/removeDuplicate/${histName}_bowtie2.sorted.dupMarked.sam \
#METRICS_FILE=$projPath/alignment/removeDuplicate/picard_summary/${histName}_picard.dupMark.txt
## remove duplicates
## 去除重复 reads 
#picardCMD \
#MarkDuplicates I=$projPath/alignment/sam/${histName}_bowtie2.sorted.sam \
#O=$projPath/alignment/removeDuplicate/${histName}_bowtie2.sorted.rmDup.sam \
#REMOVE_DUPLICATES=true \
#METRICS_FILE=$projPath/alignment/removeDuplicate/picard_summary/${histName}_picard.rmDup.txt

#7. sam to bam  bam to bed and bigwig
## Filter and keep the mapped read pairs
## 筛选和保留比对上的双端 reads 
samtools view -bS -F 0x04 $projPath/alignment/sam/${histName}_bowtie2.sam -o $projPath/alignment/bam/${histName}_bowtie2.bam
samtools sort $projPath/alignment/bam/${histName}_bowtie2.bam -o $projPath/alignment/bam/${histName}_bowtie2.sort.bam
samtools index $projPath/alignment/bam/${histName}_bowtie2.sort.bam
## Convert into bed file format
## 将 BAM 文件转换为 bed 文件格式
bedtools bamtobed -i $projPath/alignment/bam/${histName}_bowtie2.bam -bedpe > $projPath/alignment/bed/${histName}_bowtie2.bed
## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
awk '$1==$4 && $6-$2 < 1000 {print $0}' $projPath/alignment/bed/${histName}_bowtie2.bed > $projPath/alignment/bed/${histName}_bowtie2.clean.bed
## Only extract the fragment related columns
cut -f 1,2,6 $projPath/alignment/bed/${histName}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  >$projPath/alignment/bed/${histName}_bowtie2.fragments.bed

## Spike-in calibration convert to bedgraph
if [[ "$seqDepth" -gt "1" ]]; then
    mkdir -p $projPath/alignment/bedgraph
    scale_factor=`echo "10000 / $seqDepth" | bc -l`
    echo "Scaling factor for $histName is: $scale_factor!"
    bedtools genomecov -bg -scale $scale_factor -i $projPath/alignment/bed/${histName}_bowtie2.fragments.bed -g $chromSize > $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph
fi

# convert to bigwig文件
bamCoverage -bs 25 -p max --ignoreDuplicates --normalizeUsing RPKM -b $projPath/alignment/bam/${histName}_bowtie2.sort.bam -o $projPath/alignment/bigwig/${histName}_bowtie2.bw

# 8.merge bigwigs and bedgraphs (wiggletools deeptools bedtools)
bigwigCompare -b1 /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bigwig/10_bowtie2.bw -b2 /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bigwig/11_bowtie2.bw \
  --operation mean \ #log2,ratio,subtract,add,mean,reciprocal_ratio,first,second
  -o /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bigwig/temp_merged_mean.bw
bigwigCompare -b1 /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bigwig/10_bowtie2.bw -b2 /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bigwig/11_bowtie2.bw \
  --operation mean \ #log2,ratio,subtract,add,mean,reciprocal_ratio,first,second
  -o /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bigwig/temp_merged_mean.bw
  
wiggletools write merged.bedGraph mean file1.bw file2.bw file3.bw  
bedGraphToBigWig merged.bedGraph chrom.sizes merged.bw

bedtools unionbedg -i /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bedgraph/10_bowtie2.fragments.normalized.bedgraph /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bedgraph/11_bowtie2.fragments.normalized.bedgraph /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bedgraph/12_bowtie2.fragments.normalized.bedgraph > /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bedgraph/merged.txt 
awk '{sum=0; for(i=4;i<=NF;i++) sum+=$i; $4=sum/(NF-3); print $1, $2, $3, $4}' /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bedgraph/merged.txt > /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bedgraph/H3K9me3_control_merged_mean.bedGraph                     

bedtools unionbedg -i /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bedgraph/13_bowtie2.fragments.normalized.bedgraph /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bedgraph/14_bowtie2.fragments.normalized.bedgraph /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bedgraph/15_bowtie2.fragments.normalized.bedgraph > /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bedgraph/merged.txt 
awk '{sum=0; for(i=4;i<=NF;i++) sum+=$i; $4=sum/(NF-3); print $1, $2, $3, $4}' /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bedgraph/merged.txt > /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bedgraph/H3K9me3_case_merged_mean.bedGraph                     
                      
bedtools unionbedg -i /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bedgraph/16_bowtie2.fragments.normalized.bedgraph /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bedgraph/17_bowtie2.fragments.normalized.bedgraph /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bedgraph/18_bowtie2.fragments.normalized.bedgraph > /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bedgraph/merged.txt 
awk '{sum=0; for(i=4;i<=NF;i++) sum+=$i; $4=sum/(NF-3); print $1, $2, $3, $4}' /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bedgraph/merged.txt > /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bedgraph/H3K27me2_control_merged_mean.bedGraph                     
  
bedtools unionbedg -i /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bedgraph/19_bowtie2.fragments.normalized.bedgraph /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bedgraph/20_bowtie2.fragments.normalized.bedgraph /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bedgraph/21_bowtie2.fragments.normalized.bedgraph > /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bedgraph/merged.txt 
awk '{sum=0; for(i=4;i<=NF;i++) sum+=$i; $4=sum/(NF-3); print $1, $2, $3, $4}' /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bedgraph/merged.txt > /data/shijian/ANNO_XS01KF2024060289_PM-XS01KF2024060289-25/alignment/bedgraph/H3K27me2_case_merged_mean.bedGraph                      



# 9. SEACR call peak 
##== linux command ==##
seacr="/home/shijian/software/SEACR/SEACR_1.3.sh"
#histControl=$2
mkdir -p $projPath/peakCalling/SEACR
bash $seacr $projPath/alignment/bedgraph/H3K9me3_case_merged_mean.bedGraph \
     $projPath/alignment/bedgraph/H3K9me3_control_merged_mean.bedGraph \
     norm stringent $projPath/peakCalling/SEACR/H3K9me3_seacr.peaks
bash $seacr $projPath/alignment/bedgraph/H3K27me2_case_merged_mean.bedGraph \
     $projPath/alignment/bedgraph/H3K27me2_control_merged_mean.bedGraph \
     norm stringent $projPath/peakCalling/SEACR/H3K27me2_seacr.peaks 

bash $seacr $projPath/alignment/bedgraph/H3K9me3_case_merged_mean.bedGraph \
     $projPath/alignment/bedgraph/H3K9me3_control_merged_mean.bedGraph \
     norm relaxed $projPath/peakCalling/SEACR/H3K9me3_seacr.peaks
bash $seacr $projPath/alignment/bedgraph/H3K27me2_case_merged_mean.bedGraph \
     $projPath/alignment/bedgraph/H3K27me2_control_merged_mean.bedGraph \
     norm relaxed $projPath/peakCalling/SEACR/H3K27me2_seacr.peaks



