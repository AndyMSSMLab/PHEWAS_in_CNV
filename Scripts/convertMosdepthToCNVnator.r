library(tidyverse)
library(data.table)

convertMosdepthToCNVnator <- function(mos, cnvnator){
  ## mos is dataframe containg three columns: Loc, sampleID and normRD (normalized read depth from mosdepth, GATK)
  ## cnvnator contains three columns Loc, sampleID, cnvnator (copy nuber estimates from CNVnator)
  ## Loc can be gene names for multicopy genes or VNTR coordinates
  ## mos will contain data for all topmed samples
  ## cnvnator will contain data for ~200 randomly selected samples
  
  data = mos %>$ inner_join(cnvnator, by = c("Loc", "sampleID"))
  
  ## Generate linear regression model for each loci
  model = lapply(unique(data$Loc), function(loc){
    d = data %>% filter(Loc ==loc)
    lm(cnvnator~normRD, d)
  })
  names(model) = unique(data$Loc)
  
  ## Convert normalized read depth of all samples into relative copy number
  out = lapply(names(model), function(loc){
    d = mos %>% filter(geneNames ==loc)
    predict = predict(model[[loc]],newdata = d %>% select(normRD))
    d %>% mutate(cnvnator = predict)
  })
  out = do.call(rbind, out) 
  out
  
}
