#去SRA或者ENA下载整个项目的SRR注释大表，一般ENA字段最全
# prefetch下载 
prefetch --ascp-path "/pub6/temp/shijian/miniconda3/pkgs/aspera-cli-3.7.7-0/bin/ascp|/pub6/temp/shijian/miniconda3/pkgs/aspera-cli-3.7.7-0/etc/asperaweb_id_dsa.openssh" --option-file /pub6/temp/shijian/SRR_Acc_List1.txt -O /pub6/temp/shijian/SRP116382
prefetch --option-file /pub6/temp/shijian/SRR_Acc_List1.txt -O /pub6/temp/shijian/SRP116382
prefetch --option-file /pub6/temp/shijian/SRR_Acc_List.txt --max-size 100000000 --ascp-path "/home/shijian/miniconda3/bin/ascp|/home/shijian/miniconda3/etc/asperaweb_id_dsa.openssh" 

# aspera下载(目前不稳定)
/pub5/xiaoyun/Software/aspera/connect/bin/ascp -i /pub5/xiaoyun/Software/aspera/connect/etc/asperaweb_id_dsa.openssh -QT -l300m anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR713/SRR7138443/SRR7138443.sra /pub6/temp/shijian/SRP116382 
#' @param SRRTable 所有SRR号
#' @param outDir sra下载地址
#' @param outpath 命令行输出地址
#' @return 返回所有aspera下载命令行
#' @author shi jian 
#/pub5/xiaoyun/Software/aspera/connect/bin/ascp -i /pub5/xiaoyun/Software/aspera/connect/etc/asperaweb_id_dsa.openssh -QT -l300m anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR713/SRR7138443/SRR7138443.sra /pub6/temp/shijian/SRP116382 
BatchDownloadAspera<- function(SRRTable,outDir,outpath){
  srrtable <- read.table(SRRTable,sep="\t",header=F,fill=T,quote=NULL,stringsAsFactors=F)
  command_pre <- "/pub5/xiaoyun/Software/aspera/connect/bin/ascp -i /pub5/xiaoyun/Software/aspera/connect/etc/asperaweb_id_dsa.openssh -QT -l300m anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/"
  if(!file.exists(outDir)){
   dir.create(outDir)
  }
  aspera_commands <- c()
  for(i in 1:nrow(srrtable)){
    temp <- file.path(substring(srrtable[i,1],1,6),srrtable[i,1],paste0(srrtable[i,1],".sra"))
    tt <- paste0(command_pre, temp)
    tt <- paste(tt,outDir)
    aspera_commands[i] <- tt
  }
  writeLines(aspera_commands,con=outpath)  
}
BatchDownloadAspera("D:/Rsources/Project/StudySingleCell/data/SRR_Acc_List.txt","/pub6/temp/shijian/SRP116382","D:/Rsources/Project/StudySingleCell/data/aspera.sh")
./aspera.sh
