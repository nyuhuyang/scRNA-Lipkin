library(Seurat)
library(readxl)
library(dplyr)
library(tidyr)
library(magrittr)
library(kableExtra)
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

(load(file = "data/LynchSyndrome_6_20190802.Rda"))


(df <- table(object$singler1main,object$conditions) %>% as.data.frame())
colnames(df) = c("Cell.type", "conditions", "Freq")
(df %<>% spread(conditions, Freq))
rownames(df) = df$Cell.type
df = df[,-1]
p_value = c()
ColSum <- colSums(df)

for(i in 1:nrow(df)){
        conting <- rbind(df[i,],ColSum-df[i,])
        FISH <- fisher.test(conting,conf.int = T)
        p_value[i] = FISH$p.value
        #CHI = chisq.test(conting, correct = T)
        #chisq_p_value[i] = CHI$p.value             
}

df$p_value = p_value
df$p_val_adj = p.adjust(p = df$p_value, method = "bonferroni", 
                               n = nrow(df))
df %>% kable %>% kable_styling()
write.csv(df,paste0(path,"cell_numbers.csv"))
