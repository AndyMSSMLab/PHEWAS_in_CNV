library(tidyverse)
library(data.table)
filterSamplesVNTR<-function(data){
  ### data object is a data.frame and must contain five columns
  ### Coord, sampleID, VNTR, Flank3, Flank5
  ### Coord is VNTR cordinates in format chr:start-end
  ### VNTR, Flank3 and Flank5 are normalized read depth (from mosdepth and GATK) for VNTR
  ###    3 prime flank and 5 prime flank respectively
  
  ## Calculate upper and lower found for 3 prime and 5 prime flanks to filter samples
  ## select samples that are between 30% and 70% quantile for each VNTR
  ## calculate mean and stdanrd deviation usint filtered samples
  ## calculate upper and lower bound by mean +/- 7*std dev
  stat_5 <- data %>% group_by(Coord) %>%
    filter(Flank5 > quantile(Flank5, 0.3) & Flank5 < quantile(Flank5, 0.7)) %>%
    summarize(mean = mean(Flank5), sd = sd(Flank5)) %>% as.data.frame %>%
    mutate(upper = mean + 7*sd, lower = mean - 7*sd) %>% 
    select(Coord, upper, lower) %>%
    rename(upper_5 = upper, lower_5 = lower)
  
  stat_3 <- data %>% group_by(Coord) %>%
    filter(Flank3 > quantile(Flank3, 0.3) & Flank3 < quantile(Flank3, 0.7)) %>%
    summarize(mean = mean(Flank3), sd = sd(Flank3)) %>% as.data.frame %>%
    mutate(upper = mean + 7*sd, lower = mean - 7*sd) %>% 
    select(Coord, upper, lower) %>%
    rename(upper_3 = upper, lower_3 = lower)
  
  ##Remove samples with normalized read depth higher than upper bound in both 3' and 5' flanks
  ## or lower than lower bound in both flanks
  out <- data %>% inner_join(stat_5, by = "Coord") %>%
    inner_join(stat_3, by = "Coord") %>%
    group_by(Coord) %>%
    filter(!(
      (Flank5 > upper_5 & Flank5 > upper_3) |
        (Flank3 < lower_5 & Flank3 < lower_3)
    )) %>%
    ungroup %>% select(Coord, sampleID, VNTR, Flank3, Flank5)
  
  out
}