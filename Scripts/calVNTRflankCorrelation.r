library(tidyverse)
library(data.table)
calVNTRflankCorrelation<-function(data){
  ### data object is a data.frame and must contain five columns
  ### Coord, sampleID, VNTR, Flank3, Flank5
  ### Coord is VNTR cordinates in format chr:start-end
  ### VNTR, Flank3 and Flank5 are normalized read depth (from mosdepth and GATK) for VNTR
  ###    3 prime flank and 5 prime flank respectively
  
  data %>% group_by(Coord) %>%
    summarize(
      cor_3_flank = cor(VNTR, Flank3),
      cor_5_flank = cor(VNTR, Flank5)
    ) %>%
    ungroup %>% as.data.frame
  
}