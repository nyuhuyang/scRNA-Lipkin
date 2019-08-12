library(Seurat)
library(readxl)
library(dplyr)
library(tidyr)
library(magrittr)
library(kableExtra)
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

(load(file = "data/LynchSyndrome_6_20190802.Rda"))

# Fisher's Exact Test for values in Seurat
FisherTestSeurat <- function(object,var = c("singler1main","orig.ident"),percentage = FALSE,
                             file.name = "cell_numbers",...){
        (df <- table(object@meta.data[,var[1]],object@meta.data[,var[2]]) %>% as.data.frame %>%
                 spread(Var2, Freq))
        colnames(df)[1] = "Cell.type"
        rownames(df) = df$Cell.type
        df = df[,-1]
        
        if(percentage) {
                (df1 <- table(object$singler1main,object$conditions) %>% prop.table(2) %>%
                         as.data.frame %>%
                         spread(Var2, Freq))
                colnames(df1)[1] = "Cell.type"
                rownames(df1) = df1$Cell.type
                df1 = df1[,-1]
        }        
        p_value = c()
        ColSum <- colSums(df)
        
        for(i in 1:nrow(df)){
                (conting <- rbind(df[i,],ColSum-df[i,]))
                FISH <- fisher.test(conting,workspace = 2000000,...)
                (p_value[i] = FISH$p.value)
                #CHI = chisq.test(conting, correct = T)
                #chisq_p_value[i] = CHI$p.value             
        }
        
        df$p_value = p_value
        df$p_val_adj = p.adjust(p = df$p_value, method = "bonferroni", 
                                n = nrow(df))
        if(percentage) df = cbind(df1,df)
        write.csv(df,paste0(path,file.name,".csv"))
        return(df)
}
FisherTestSeurat(object,var = c("singler1main","conditions"),
                 percentage= F, file.name = "cell_numbers") %>% kable %>% kable_styling
FisherTestSeurat(object,var = c("singler1main","orig.ident"),
                 percentage= F, hybrid=TRUE,file.name = "cell_numbers_6") %>% kable %>% kable_styling
#===
(df <- table(object$singler1main,object$conditions) %>% as.data.frame %>%
                spread(Var2, Freq))
colnames(df)[1] = "Cell.type"
rownames(df) = df$Cell.type
df = df[,-1]

#===


#===
(df2 <- table(object$singler1main,object$orig.ident) %>% as.data.frame %>%
                spread(Var2, Freq))
colnames(df2)[1] = "Cell.type"
rownames(df2) = df2$Cell.type
df2 = df2[,-1]

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
df = cbind(df,df1)
df %>% kable %>% kable_styling()
write.csv(df,paste0(path,"cell_numbers.csv"))

fisher.test(df2, hybrid=TRUE,simulate.p.value=TRUE)
