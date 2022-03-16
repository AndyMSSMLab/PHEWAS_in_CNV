library(tidyverse)
library(data.table)
getMax <- function(motif){
  nuc = c("A","T","G","C")
  seq = expand.grid(nuc, nuc) %>% as.data.frame %>%
    mutate(Var1 = as.character(Var1), Var2= as.character(Var2)) %>%
    filter(Var1!=Var2) %>% rowwise %>% mutate(seq  = paste(sort(c(Var1,Var2)), collapse = "")) %>%
    select(seq) %>% distinct %>% pull
  
  out = sapply(seq, function(s){str_count(motif,s)})
  out = round(out/sum(out),3)
  out = out[which.max(out)]
  name = names(out)
  names(out) = NULL
  data.frame(MaxDinucPerc = out, MaxDinuc = name)
}

data = fread("VNTR_motif.txt") %>% select(VNTR, motif) %>% distinct %>% 
  separate_rows(motif, sep = ",") %>%
  as.data.frame

frac = cbind(data,do.call(rbind,lapply(data[,2], getMax))) %>% as.data.frame  %>%
  group_by(VNTR) %>%
  summarize(
    MaxDinucPerc = paste(MaxDinucPerc,collapse=",",sep=","),
    MaxDinuc =paste(MaxDinuc, collapse = ",")
  ) %>%
  ungroup %>% as.data.frame

out = frac %>% rowwise %>% 
  mutate(flag = !any(as.numeric(unlist(strsplit(MaxDinucPerc,split = ","))) >0.8)) %>% 
  as.data.frame

write.table(out, "VNTR.filtered.lowcomplexity.txt", row.names = F, quote = F, sep = "\t")
