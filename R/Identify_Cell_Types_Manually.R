library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
marker_path <- paste0(path,"markers/")
if(!dir.exists(marker_path))dir.create(marker_path, recursive = T)

# ======== 2.1 =========== test with known markers==================
(load(file="data/LynchSyndrome_6_20190802.Rda"))
DefaultAssay(object) <- "integrated"
#df_markers <- readxl::read_excel("doc/Lynch mouse model scRNAseq genes of interest 073119.xlsx")

df_markers <- readxl::read_excel("../seurat_resources/bio-rad-markers.xlsx")
colnames(df_markers) = gsub(" ","_",colnames(df_markers))
colnames(df_markers) = gsub(":|\\/","_",colnames(df_markers))
colnames(df_markers) = gsub("\\+","",colnames(df_markers))
markers = df_markers[,grep("Alias",colnames(df_markers),invert = T)]
marker.list <- df2list(markers)

marker.list %<>% lapply(function(x) x[1:18]) %>% 
     lapply(function(x) FilterGenes(object,x)) %>% 
     lapply(function(x) x[!is.na(x)]) %>% 
    lapply(function(x) x[1:min(length(x),12)])
marker.list <- marker.list[!is.na(marker.list)]
marker.list %>% list2df %>% t %>% kable() %>% kable_styling()

object@meta.data = object@meta.data[,grep("ident",colnames(object@meta.data),invert = T)]
Idents(object) = "integrated_snn_res.0.6"
for(i in 1:length(marker.list)){
    if(length(marker.list[[i]]) == 0) next
    p <- lapply(marker.list[[i]], function(marker) {
        FeaturePlot(object = object, feature = marker,pt.size = 0.5,
                    reduction="tsne", label = T)+
            NoLegend()+
            ggtitle(paste0(marker,Alias(df = df_markers,gene = marker)))+
            theme(plot.title = element_text(hjust = 0.5,size = 15,face = "plain"))
    })
    jpeg(paste0(path,"markers/",names(marker.list)[i],".jpeg"),units="in", width=10, height=7,res=600)
    print(do.call(cowplot::plot_grid, p)+ ggtitle(paste(names(marker.list)[i],"markers"))+
              theme(plot.title = element_text(hjust = 0.5,size = 20)))
    dev.off()
    print(paste0(i,":",length(marker.list)))
}

object_data <- object@assays$RNA@data
save(object_data, file = "data/LynchSyndrome_6_20190802.Rda")
