########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(gplots)
library(MAST)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 
# 3.1.1 load data
# Rename ident
(load(file = "data/LynchSyndrome_6_20190802.Rda"))

object@meta.data$celltype.conditions = paste0(object@meta.data$singler1main,"_",
                                                object@meta.data$conditions)
(df_table <- table(object$celltype.conditions) %>% as.data.frame() %>% 
                .[order(.$Freq, decreasing = T)[1:12],])

(celltypes <- df_table$Var1 %>% gsub("_.*","",.) %>% unique)
ident.1 <- paste0(celltypes,"_Contorl")
ident.2 <- paste0(celltypes,"_Naproxen")

Idents(object) = "celltype.conditions"
table(Idents(object))
subfolder <- paste0(path,"DEG/")
DefaultAssay(object) = "SCT"
gde.pair <- FindPairMarkers(object, ident.1 = ident.1, ident.2 = ident.2,
                            logfc.threshold = 0.01, min.cells.group =3,assay.type = "SCT",
                            return.thresh = 1, only.pos = FALSE, save.path = subfolder)
#gde.pair = gde.pair[gde.pair$p_val_adj< 0.25,]
write.csv(gde.pair, paste0(subfolder,"pairwise_comparision_SCT.csv"))
gde.pair = read.csv("output/20190803/DEG/pairwise_comparision.csv",row.names = 1)
head(gde.pair,10) %>% kable %>% kable_styling

(titles <- paste(ident.1, "vs.", ident.2))
# Volcano plot=========
(clusters <- unique(gde.pair$cluster1.vs.cluster2))
for(i in 1:length(clusters)){
        df <- gde.pair[gde.pair$cluster1.vs.cluster2 %in% clusters[i],]
        df$log10_p_val_adj = -log10(df$p_val_adj)
        df$log10_p_val_adj[df$log10_p_val_adj == "Inf"] = 400
        left = df[df$avg_logFC < -0.1,]
        right = df[df$avg_logFC > 0.1,]
        left = rownames(left)[left$log10_p_val_adj >= head(sort(left$log10_p_val_adj,decreasing = T),15) %>% tail(1)]
        right = rownames(right)[right$log10_p_val_adj >= head(sort(right$log10_p_val_adj,decreasing = T),15) %>% tail(1)]
        g <- ggplot(df,aes(avg_logFC,log10_p_val_adj)) + 
                geom_point() + 
                ggtitle(titles[i]) + 
                ylab("-log10(p_value_adj)")+
                theme_minimal()+
                theme(plot.title = element_text(size=20,hjust = 0.5))+
                ggrepel::geom_text_repel(aes(label = gene), 
                                         data=df[c(left,right),]) +
                geom_point(color = ifelse((df$avg_logFC > 0.1  & df$p_val_adj < 0.05) , "red",
                                          ifelse((df$avg_logFC < -0.1 & df$p_val_adj < 0.05), "blue","gray")))
        jpeg(paste0(path,"Volcano_plot",clusters[i],".jpeg"), units="in", width=10, height=7,res=600)
        print(g)
        dev.off()
}

subfolder <- paste0(path,"DEG/")
DefaultAssay(object) = "SCT"
object@meta.data$integrated_snn_res.0.6.conditions = paste0(object@meta.data$integrated_snn_res.0.6,"_",
                                              object@meta.data$conditions)

Idents(object) = "integrated_snn_res.0.6.conditions"
Foxp3<- AverageExpression(object,features = c("Foxp3","Lef1"), assays = "SCT")
Foxp3$SCT %>% kable %>% kable_styling()
