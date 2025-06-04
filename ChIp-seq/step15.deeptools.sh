#deeptools pipeline
## run compute matrix to collect the data needed for plotting
computeMatrix scale-regions \
                       -S /home/shijian/CXJ--Cut-Run/GSM5482194_PFC_WT_H3K4me3_Rep1.bw \
                          /home/shijian/CXJ--Cut-Run/GSM5482195_PFC_WT_H3K4me3_Rep2.bw \
                          /home/shijian/CXJ--Cut-Run/GSM5482196_PFC_Setd1aHet_H3K4me3_Rep1.bw \
                          /home/shijian/CXJ--Cut-Run/GSM5482197_PFC_Setd1aHet_H3K4me3_Rep2.bw \ #bigwig文件
                       -R /home/shijian/CXJ--Cut-Run/gencode.v48.annotation.gtf.gz \  #基因组区域注释文件
                       -b 3000 \ #upstream（区域前）保留长度（单位：bp）
                       --regionBodyLength 5000 \  #中间主体区域拉伸到统一长度（常设为 5k 或 10k）
                       -a 3000 \ #downstream（区域后）保留长度
                       --binSize  50 \  #每个 bin 的大小，默认为 10bp
                       --averageTypeBins mean \  # alternatives: median, min, max, sum
                       -p 25 \  # number of threads
                       --missingDataAsZero \  #缺失数据用 0 补
                       --skipZeros \  #排除全 0 区域，避免影响均值
                       -out /home/shijian/CXJ--Cut-Run/scaled_matrix.gz \  #必选 to be used with plotHeatmap and plotProfile
                       --outFileNameMatrix /home/shijian/CXJ--Cut-Run/scaled_matrix.tab \  #可选 scaled_matrix.gz的人类可读版
                       --outFileSortedRegions /home/shijian/CXJ--Cut-Run/scaled_sorted_regions.bed #可选 write out the values where each row corresponds to one region in the BED file
#运行必须放在一行
computeMatrix scale-regions -S /home/shijian/CXJ--Cut-Run/GSM5482194_PFC_WT_H3K4me3_Rep1.bw /home/shijian/CXJ--Cut-Run/GSM5482195_PFC_WT_H3K4me3_Rep2.bw /home/shijian/CXJ--Cut-Run/GSM5482196_PFC_Setd1aHet_H3K4me3_Rep1.bw /home/shijian/CXJ--Cut-Run/GSM5482197_PFC_Setd1aHet_H3K4me3_Rep2.bw -R /home/shijian/CXJ--Cut-Run/gencode.v48.annotation.gtf.gz -b 3000 --regionBodyLength 5000 -a 3000 --binSize  50 --averageTypeBins mean -p 25 --missingDataAsZero --skipZeros -out /home/shijian/CXJ--Cut-Run/scaled_matrix.gz --outFileNameMatrix /home/shijian/CXJ--Cut-Run/scaled_matrix.tab --outFileSortedRegions /home/shijian/CXJ--Cut-Run/scaled_sorted_regions.bed 

#运行不好使，必须是第一行的好使
computeMatrix  reference-point  \
               --referencePoint TSS  \  # alternatives: TES, center
               -S  /home/shijian/CXJ--Cut-Run/GSM5482194_PFC_WT_H3K4me3_Rep1.bw  /home/shijian/CXJ--Cut-Run/GSM5482195_PFC_WT_H3K4me3_Rep2.bw  /home/shijian/CXJ--Cut-Run/GSM5482196_PFC_Setd1aHet_H3K4me3_Rep1.bw   /home/shijian/CXJ--Cut-Run/GSM5482197_PFC_Setd1aHet_H3K4me3_Rep2.bw \
               -R  /home/shijian/CXJ--Cut-Run/Mus_musculus.GRCm38.102.gtf.gz  \
               --skipZeros  \  #排除全 0 区域，避免影响均值
               --binSize 10  \  # bin size in bp
               --averageTypeBins mean  \  # alternatives: median, min, max, sum
               -p 25  \  # number of threads
               -b 5000  -a 5000  \  #define the region you are interested in
               --missingDataAsZero  \  #缺失数据用 0 补
               -out /home/shijian/CXJ--Cut-Run/TSS_matrix.gz  \  #必选 to be used with plotHeatmap and plotProfile
               --outFileNameMatrix /home/shijian/CXJ--Cut-Run/TSS_matrix.tab  \  #可选 TSS_matrix.gz的人类可读版
               --outFileSortedRegions /home/shijian/CXJ--Cut-Run/sorted_regions.bed #可选 write out the values where each row corresponds to one region in the BED file
#运行必须放在一行
computeMatrix  reference-point  --referencePoint TSS  -S  /home/shijian/CXJ--Cut-Run/GSM5482194_PFC_WT_H3K4me3_Rep1.bw  /home/shijian/CXJ--Cut-Run/GSM5482195_PFC_WT_H3K4me3_Rep2.bw  /home/shijian/CXJ--Cut-Run/GSM5482196_PFC_Setd1aHet_H3K4me3_Rep1.bw   /home/shijian/CXJ--Cut-Run/GSM5482197_PFC_Setd1aHet_H3K4me3_Rep2.bw  -R  /home/shijian/CXJ--Cut-Run/Mus_musculus.GRCm38.102.gtf.gz  --skipZeros  --binSize 10  --averageTypeBins mean  -p 25  -b 5000  -a 5000  --missingDataAsZero  -out /home/shijian/CXJ--Cut-Run/TSS_matrix.gz  --outFileNameMatrix /home/shijian/CXJ--Cut-Run/TSS_matrix.tab  --outFileSortedRegions /home/shijian/CXJ--Cut-Run/sorted_regions.bed

#绘制信号热图
plotHeatmap -h #查看参数说明
plotHeatmap -m /home/shijian/CXJ--Cut-Run/TSS_matrix.gz \  #computeMatrix计算的TSS_matrix.gz
        -out /home/shijian/CXJ--Cut-Run/signalHeatmap.png \ #输出的热图文件
        --sortUsing  mean \ # alternatives: max, min, sum, median, mean, none，region_length()
        --sortRegions descend \ #descend or ascend or no or keep
        --colorMap RdBu viridis seismic YlGnBu \ #设置颜色
        --heatmapWidth 4 \ #热图宽度
        --heatmapHeight 28 \ #热图高度
        --zMin -3 --zMax 3 \ #使用zMin zMax限定颜色范围
        --kmeans 2 \ #kmeans聚类数目
        --plotTitle 'H3K4me3 enrichment at TSS' \ #热图标题
        --refPointLabel 'TSS' \ #参考点标签
        --regionsLabel 'cluster1 genes' 'cluster2 gene' \ #左侧区域标签
        --samplesLabel 'PFC_WT_H3K4me3_Rep1' 'PFC_WT_H3K4me3_Rep2' 'PFC_Setd1aHet_H3K4me3_Rep1' 'PFC_Setd1aHet_H3K4me3_Rep2' \ #列样本标签  
        --dpi 600 #输出图片分辨率（dpi）
#运行必须放在一行
plotHeatmap -m /home/shijian/CXJ--Cut-Run/TSS_matrix.gz -out /home/shijian/CXJ--Cut-Run/signalHeatmap.png --sortUsing  mean --sortRegions descend --colorMap RdBu viridis seismic YlGnBu --heatmapWidth 4 --heatmapHeight 28 --zMin -3 --zMax 3 --kmeans 2 --plotTitle 'H3K4me3 enrichment at TSS' --refPointLabel 'TSS' --regionsLabel 'cluster1 genes' 'cluster2 gene' --samplesLabel 'PFC_WT_H3K4me3_Rep1' 'PFC_WT_H3K4me3_Rep2' 'PFC_Setd1aHet_H3K4me3_Rep1' 'PFC_Setd1aHet_H3K4me3_Rep2' --dpi 600 


#绘制信号曲线图
plotProfile -h #查看参数说明
plotProfile -m /home/shijian/CXJ--Cut-Run/TSS_matrix.gz \
            -out /home/shijian/CXJ--Cut-Run/signalProfile.png \ #输出图片路径（支持 .png/.pdf/.svg/.eps）
            --outFileNameData /home/shijian/CXJ--Cut-Run/profile_values.txt \ #输出 profile 的数值（可用于统计分析）
            --perGroup \
            --samplesLabel 'PFC_WT_H3K4me3_Rep1' 'PFC_WT_H3K4me3_Rep2' 'PFC_Setd1aHet_H3K4me3_Rep1' 'PFC_Setd1aHet_H3K4me3_Rep2' \ #样本标签
            --colors red blue green black \ #颜色设置
            --legendLocation  best \ 图例位置，如 best, upper-right, none
            --plotTitle 'H3K4me3 enrichment at TSS' \
            --averageType mean \ # alternatives: median, min, max, std,sum
            --refPointLabel 'TSS' \ #中心参考点的标签，如 “TSS”、“TES” 等
            --plotHeight  7 \
            --plotWidth  11 \
            --dpi 600

plotProfile -m /home/shijian/CXJ--Cut-Run/TSS_matrix.gz -out /home/shijian/CXJ--Cut-Run/signalProfile.png --outFileNameData /home/shijian/CXJ--Cut-Run/profile_values.txt --perGroup --samplesLabel 'PFC_WT_H3K4me3_Rep1' 'PFC_WT_H3K4me3_Rep2' 'PFC_Setd1aHet_H3K4me3_Rep1' 'PFC_Setd1aHet_H3K4me3_Rep2' --colors red blue green black --legendLocation  best --plotTitle 'H3K4me3 enrichment at TSS' --averageType mean --refPointLabel 'TSS' --plotHeight  7 --plotWidth  11 --dpi 600
