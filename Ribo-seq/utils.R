AggregateFeatureCount <- function(inputDir = "/data/shijian/project/Coorperation/HanJuanGong/data_2025_cold_client/2026-03/ID26-0154_RIBO_9hsa/4.counts",
                                  pattern="counts.txt$"){
  library(data.table)
  library(dplyr)
  library(magrittr)
  filepaths <- list.files(path=inputDir,recursive = T,full.names = T,pattern=pattern) 
  sample_names <- basename(dirname(filepaths))
  #file_i <- fread(filepaths[1])
  temp <- c()
  for(i in 1:length(filepaths)){
    file_i <- fread(filepaths[i])
    file_i <- file_i[,c(1,7)]
    if( i == 1){
      temp <- c(temp,file_i)
    }else{
      temp <- c(temp,file_i[,2])
    }
  }
  temp <- as.data.frame(temp)
  temp %<>% tibble::column_to_rownames("Geneid") 
  colnames(temp) <- sample_names
  return( temp )
}
