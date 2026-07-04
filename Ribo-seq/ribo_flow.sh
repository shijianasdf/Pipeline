## function: ribo_flow

#### description 
# ribo-seq处理流程
# 1）trim_galore去除接头序列
# 2）bowtie2将单端序列比对到rrna和trna序列参考基因组，得到未比对到的read
# 3）star将未比对到的read比对到人类参考基因组
# 4）featurecount量化

#输入 ws 

#输出 1.trimGalore  
#     2.filterReads
#     3.bams
#     4.counts



#parameters
ws=$1   
adapter=AACTGTAGGCACCATCAAT
rrna_trna_index=/data/shijian/refData/human_reference/riboseq/rrna_trna_index/rrna_trna_index  #bowtie2构建的rrna和trna参考基因组序列索引
star_index=/data/shijian/refData/human_reference/human_index/star.index
gtf=/data/shijian/refData/human_reference/gencode.v48.chr_patch_hapl_scaff.annotation.gtf
path_trimGalore=${ws}/1.trimGalore
path_filtreads=${ws}/2.filterReads
path_bams=${ws}/3.bams
path_counts=${ws}/4.counts
config_file=${ws}/config.txt #配置文件，两列如下：Sh-2-2Ri /data/shijian/Sh-2-2Ri.fq.gz

## Usage
# ws=/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa
# find $ws/raw -name "*.fq.gz" | awk -F'/' '{ full_path=$0; sample_name=$NF; gsub(".fq.gz", "", sample_name); print sample_name, full_path }' > $ws/config.txt
# nohup bash ribo_flow.sh $ws > $ws/ribo.log 2>&1 &
# tail -f $ws/ribo.log

# dir
mkdir -p $path_trimGalore $path_filtreads $path_bams $path_counts

# cycle for ribo-seq
cat $config_file | while read it
do
  arr=(${it})
  case=${arr[0]}
  fq=${arr[1]}
  output_trimGalore=${path_trimGalore}/${case}
  if [ -f $fq ]; then
    mkdir -p ${output_trimGalore}
    echo `date` " trim_galore --cores 16 -a ${adapter} --length 20 --max_length 35 --fastqc -o ${output_trimGalore} ${fq}"
    trim_galore --cores 16 -a ${adapter} --length 20 --max_length 35 --fastqc -o ${output_trimGalore} ${fq}
  fi
  
  input_trimread=${output_trimGalore}/${case}_trimmed.fq.gz
  output_bowtie2=${path_filtreads}/${case}
  if [ -f $input_trimread ]; then
    mkdir -p ${output_bowtie2}
    echo `date` " bowtie2 -p 16 -x ${rrna_trna_index} -U ${input_trimread} --very-sensitive -L 18 --un-gz ${output_bowtie2}/${case}_clean.fq.gz -S ${output_bowtie2}/${case}.sam"
    bowtie2 -p 16 -x ${rrna_trna_index} -U ${input_trimread} --very-sensitive -L 18 --un-gz ${output_bowtie2}/${case}_clean.fq.gz -S ${output_bowtie2}/${case}.sam
  fi
  
  if [ -f ${output_bowtie2}/${case}_clean.fq.gz ]; then
    mkdir -p ${path_bams}/${case}
    echo `date` " STAR --outFilterType BySJout --runThreadN 16 --outFilterMismatchNmax 3 --genomeDir ${star_index} --readFilesIn ${output_bowtie2}/${case}_clean.fq.gz --outFileNamePrefix ${path_bams}/${case}/${case}_ --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --quantMode TranscriptomeSAM GeneCounts --outFilterMultimapNmax 20 --outFilterMatchNmin 16 --alignEndsType EndToEnd"
    STAR --outFilterType BySJout --runThreadN 16 --outFilterMismatchNmax 3 --genomeDir ${star_index} --readFilesIn ${output_bowtie2}/${case}_clean.fq.gz --outFileNamePrefix ${path_bams}/${case}/${case}_ --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --quantMode TranscriptomeSAM GeneCounts --outFilterMultimapNmax 20 --outFilterMatchNmin 16 --alignEndsType EndToEnd
  fi
  
  #feature counts量化
  if [ -f ${path_bams}/${case}/${case}_Aligned.sortedByCoord.out.bam ]; then
    echo `date` " featureCounts -T 8 -t CDS -g gene_id -a ${gtf} -o ${path_counts}/${case}/counts.txt --primary -s 1 ${path_bams}/${case}/${case}_Aligned.sortedByCoord.out.bam"
    mkdir -p ${path_counts}/${case}
    featureCounts -T 8 -t CDS -g gene_id -a ${gtf} -o ${path_counts}/${case}/counts.txt --primary -s 1 ${path_bams}/${case}/${case}_Aligned.sortedByCoord.out.bam
  fi
done
# End ribo_flow
