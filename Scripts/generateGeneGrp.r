library(tidyverse)
library(data.table)
library(reshape2)

getGrpName<-function(words){
  ##Get largest common substring from gene names in each gene group
  ##This script is modification from the one available on stackoverflow.com
  ##https://stackoverflow.com/questions/26285010/r-find-largest-common-substring-starting-at-the-beginning
  words.split <- strsplit(words, '')
  words.split <- lapply(words.split, `length<-`, max(nchar(words)))
  words.mat <- do.call(rbind, words.split)
  common.substr.length <- which.max(apply(words.mat, 2, function(col) !length(unique(col)) == 1)) - 1
  substr(words[1], 1, common.substr.length)
}

getGeneGrp<-function(data){
  ### data object is a data.frame and must contain three columns
  ### geneNames, sampleID and normRD
  ### normRD is read depth from mosdepth after GATK normalization
  data = data %>% select(geneNames, sampleID, normRD) %>% 
    spread(geneNames, normRD) %>% column_to_rownames("sampleID")
  
  ##calculating pairwise correlation 
  cor = cor(data, method = "pearson")
  cor[lower.tri(cor,diag=T)] = NA
  cor = as.data.frame(cor) %>% rownames_to_column("geneNames.x") %>%
    melt(id.vars = "geneNames.x", variable.name = "geneNames.y", value.name = "cor") %>%
    filter(!is.na(cor)) %>% 
    arrange(desc(cor))
  cor = cor %>% filter(abs(cor) >= 0.9) ##selecting correlation >0.9

  ##Refseq gene annotation from UCSC genome browser
  ##columns chr, start, end, name2
  coord = fread("Refseq_genes_hg38.bed", header = F, sep = "\t") %>% as.data.frame %>%
    rename(chr = 1, start =2 , end = 3, geneNames = 4) %>% select(chr, geneNames)
  
  ## adding chr columns to the gene pairs
  cor = cor %>% inner_join(coord, by = c(geneNames.x = "geneNames")) %>% rename(chr.x = "chr") %>%
    inner_join(coord, by = c(geneNames.y = "geneNames")) %>% rename(chr.y = "chr") %>%
    select(geneNames.x, chr.x, geneNames.y, chr.y, cor)
  
  ## keeping only those pairs that are on same chromosome as well as 
  ## have same first 3 characters in geneNames
  cor = cor %>% group_by(geneNames.x, geneNames.y) %>%
    summarize(flag = any(chr.x == chr.y)) %>%
    ungroup %>% as.data.frame %>% distinct %>%
    mutate(flag = flag &  substr(geneNames.x,1,3) == substr(geneNames.y,1,3)) %>%
    filter(flag)
  
  ### Grouping genes
  out = list("AMY1A") ###initialize output
  for(i in 1:nrow(z)){
    genex = z[i,"geneNames.x"]
    geney = z[i,"geneNames.y"]
    
    ###check if any one of the gene in gene pairs have been assigned to previous gene group
    index = unlist(sapply(1:length(out), function(l){ if((genex %in% out[[l]])| (geney %in% out[[l]])){l} else {NULL} }))
    if(is.null(index)){
      ##if none of the gene pair has bee assigned, create new group
      out[[length(out)+1]] =  c(genex, geney)
    }else if(length(index) == 1){
      ##if only one of the gene in gene pair has been assigned to previous gene group,
      ## add other gene in gene pair to that group
      out[[index]] = unique(c(out[[index]], genex, geney))
    }else{
      ##if both genes in gene pairs have been found but formed two different gene group
      ## merge the two gene group
      out[[index[1]]] = unique(c(out[[index[1]]], out[[index[2]]], genex, geney))
      out[[index[2]]] = NULL
    }
  }
  
  ## assigning names to gene groups
  grp_name = lapply(out, getGrpName)
  
  ##check if two gene groups have same name
  grp_name = do.call(rbind, grp_name) %>% as.data.frame %>% rename(prefix = 1) %>%
    mutate(prefix = sub("-$","",prefix)) %>%
    group_by(prefix)%>%
    mutate(total = n(), n = 1:n())  %>%
    ungroup %>%
    mutate(name = ifelse(total == 1, paste0(prefix,"_grp"), paste0(prefix,"_grp", toupper(letters)[n]) )) %>%
    as.data.frame
  names(out) = grp_name$name
  
  ##the output file will be RDS object containing list of gene groups
  ##the name of each gene group is provided by names(out)
  write_rds(out, "Genes_grp.rds")
  
}