library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
library(dplyr)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file = "data/LynchSyndrome_6_20190802.Rda"))
(load(file = "output/singler_T_LynchSyndrome_6_20190805.Rda"))
# if singler didn't find all cell labels`
length(singler$singler[[1]]$SingleR.single$labels) == ncol(object@assays$RNA@data)
if(length(singler$singler[[1]]$SingleR.single$labels) < ncol(object@assays$RNA@data)){
        all.cell = colnames(object);length(all.cell)
        know.cell = rownames(singler$singler[[1]]$SingleR.single$labels);length(know.cell)
        object = object[,know.cell]
}

table(rownames(singler$singler[[1]]$SingleR.single$labels) %in% colnames(object))
singler$meta.data$orig.ident = object@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = object@reductions$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = Idents(object) # the Seurat clusters (if 'clusters' not provided)
save(singler,file="output/singler_T_LynchSyndrome_6_20190805.Rda")

##############################
# add singleR label to Seurat
###############################
singlerDF = data.frame("singler1sub" = singler$singler[[1]]$SingleR.single$labels,
                       #"singler1main" = singler$singler[[1]]$SingleR.single.main$labels,
                       "orig.ident" = object@meta.data$orig.ident,
                       row.names = rownames(object@meta.data))

table(rownames(singlerDF) %in% colnames(object))
head(singlerDF)
apply(singlerDF,2,function(x) length(unique(x)))

##############################
# check the spearman correlation
###############################
#Or by all cell types (showing the top 50 cell types):
jpeg(paste0(path,"DrawHeatmap_sub1.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single, top.n = 50,normalize = F))
dev.off()
jpeg(paste0(path,"DrawHeatmap_sub1_N.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single,top.n = 50,normalize = T))
dev.off()

#Finally, we can also view the labeling as a table compared to the original identities:

kable(table(singlerDF$singler1sub, singlerDF$orig.ident)) %>%
        kable_styling()
singlerDF$orig.ident %>% table() %>% kable() %>% kable_styling()
singlerDF$singler1sub %>% table() %>% kable() %>% kable_styling()

singlerDF$singler1main = gsub("Hepatocytes","Epithelial cells",singlerDF$singler1main)
singlerDF$singler1main = gsub("Microglia","Macrophages",singlerDF$singler1main)
object$singler1main = gsub("Microglia","Macrophages",object$singler1main)

singlerDF$singler1sub = gsub("Tcm|Tem","T-cells",singlerDF$singler1sub)
singlerDF$singler1sub = gsub("CLP|CMP|GMP|MEP|MPP","Progenitor cells",singlerDF$singler1sub)
singlerDF$singler1sub = gsub("mv Endothelial cells","Endothelial cells",singlerDF$singler1sub)
singlerDF$singler1sub = gsub("Macrophages M1","Macrophages",singlerDF$singler1sub)

Idents(object) = "integrated_snn_res.0.6"
##############################
# process color scheme
##############################
object <- AddMetaData(object = object,metadata = singlerDF)
object <- AddMetaColor(object = object, label= "singler1sub", colors = singler.colors)

table(object$singler1sub,object$orig.ident) %>% prop.table(margin = 2) %>%  kable() %>% kable_styling()
#
Idents(object) <- "singler1sub"
object %<>% sortIdent()
TSNEPlot.1(object, group.by="singler1sub",cols = ExtractMetaColor(object),
           label = T,pt.size = 1,no.legend = T,label.repel = T,
         label.size = 4, repel = T,do.return= T,do.print = T,alpha = 0.9,
         title = "All cell types identified by Human RNA-seq reference database")
##############################
# Adjust cell type manually
##############################
cluster_2_17 <- rownames(object@meta.data)[object$integrated_snn_res.0.6 %in% c(2,17)]
object@meta.data[cluster_2_17,"singler1main"] = "Epithelial cells"
object$singler1main = gsub("Microglia","Macrophages",object$singler1main)

singler_colors <- readxl::read_excel("doc/singler.colors.xlsx")
singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
singler_colors1[duplicated(singler_colors1)]
length(singler_colors1)
apply(object@meta.data[,c("singler1sub","singler1main")],2,function(x) length(unique(x)))
object <- AddMetaColor(object = object, label= "singler1main", colors = singler_colors1)
Idents(object) <- "singler1main"
object %<>% sortIdent()
TSNEPlot.1(object, group.by="singler1main",cols = ExtractMetaColor(object),
           label = T,pt.size = 1,no.legend = T,label.repel = T,
           label.size = 4, repel = T,do.return= T,do.print = T,alpha = 0.9,
           title = "All cell types in tSNE plot")
object$singler1sub %>% table() %>% kable() %>% kable_styling()

save(object,file="data/LynchSyndrome_6_20190802.Rda")
##############################
# split tsne plot
##############################
TSNEPlot.1(object, cols = ExtractMetaColor(object),label = T,pt.size = 1,label.repel = T,
         split.by = "conditions", group.by = "singler1main",label.size = 4, repel = T, 
         no.legend = T, do.print = T,border = T,
         ncol=3,title = "Compare cell types in both conditions")

TSNEPlot.1(object, cols = ExtractMetaColor(object),label = F,pt.size = 1,
           split.by = "orig.ident", group.by = "singler1main",label.size = 4, repel = T, 
           no.legend = T, do.print = T,border = T,
           ncol=3,title = "Compare cell types in all samples")

