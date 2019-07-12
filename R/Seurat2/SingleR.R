library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
library(dplyr)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
lname1 = load(file = "./data/EC_20181022.Rda");lname1
singler = CreateSinglerObject(EC@data, annot = NULL, project.name=EC@project.name,
                              min.genes = 500,technology = "10X", species = "Mouse", citation = "",
                              normalize.gene.length = F, variable.genes = "de",
                              fine.tune = F, do.signatures = F, clusters = NULL)
# singler didn't find all cell labels
length(singler$singler[[1]]$SingleR.single$labels) == ncol(EC@data)
singler$meta.data$orig.ident = EC@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = EC@dr$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = EC@ident # the Seurat clusters (if 'clusters' not provided)
save(singler,file="./output/singler_EC_20181022.RData")
#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
lnames = load(file = "./output/singler_EC_20181022.RData")
lnames = load(file = "../SingleR/data/immgen.rda")
lnames
singler$seurat = EC
SingleR.DrawScatter(sc_data = singler$seurat@data,cell_id = 10, 
                    ref = immgen, sample_id = 232)

# Step 2: Multiple correlation coefficients per cell types are aggregated 
# to provide a single value per cell type per single-cell. 
# In the examples below we use the 80% percentile of correlation values.
# for visualization purposes we only present a subset of cell types (defined in labels.use)
out = SingleR.DrawBoxPlot(sc_data = singler$seurat@data,cell_id = 10, 
                          ref = immgen,main_types = T,
                          labels.use=c('B cells','T cells','DC','Macrophages','Monocytes','NK cells',
                                       'Mast cells','Neutrophils','Fibroblasts','Endothelial cells'))
print(out$plot)

##############################
# Human Primary Cell Atlas (HPCA)
###############################
SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single.main, top.n = Inf)
#Or by all cell types (showing the top 50 cell types):
jpeg(paste0(path,"DrawHeatmap.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single, top.n = 50))
dev.off()
jpeg(paste0(path,"DrawHeatmap_normF.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single,top.n = 50,normalize = F))
dev.off()
#Next, we can use the fine-tuned labels to color the t-SNE plot:

out1 = SingleR.PlotTsne.1(singler$singler[[1]]$SingleR.single,
                          singler$meta.data$xy,do.label=F,
                          do.letters = F,labels = singler$singler[[1]]$SingleR.single$labels,
                          label.size = 5, dot.size = 2,do.legend = T,alpha = 1,
                          label.repel = T,force=2)
jpeg(paste0(path,"/PlotTsne_sub.jpeg"), units="in", width=10, height=7,
     res=600)
out1+  ggtitle("Supervised sub cell type labeling by immgen")+
        theme(text = element_text(size=20),
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) #title in middle
dev.off()
# main types-------
out2 = SingleR.PlotTsne.1(singler$singler[[1]]$SingleR.single.main,
                          singler$meta.data$xy,do.label=T,
                          do.letters = F,labels = singler$singler[[1]]$SingleR.single.main$labels,
                          label.size = 5, dot.size = 2,do.legend = F,alpha = 1,
                          label.repel = T,force=2)
jpeg(paste0(path,"/PlotTsne_main.jpeg"), units="in", width=10, height=7,
     res=600)
out2+  ggtitle("Supervised main cell type labeling by immgen")+
        theme(text = element_text(size=20),     
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) #title in middle
dev.off()

#Finally, we can also view the labeling as a table compared to the original identities:

kable(table(singler$singler[[2]]$SingleR.single$labels,
            singler$meta.data$orig.ident)) %>%
        kable_styling()
kable(table(singler$meta.data$orig.ident,
            singler$singler[[1]]$SingleR.single.main$labels)) %>%
        kable_styling()

##############################
# add singleR label to Seurat
###############################
singlerDF = data.frame("singler1sub"=singler$singler[[1]]$SingleR.single$labels,
                       "singler1main"=singler$singler[[1]]$SingleR.single.main$labels,
                       "singler2sub"=singler$singler[[2]]$SingleR.single$labels,
                       "singler2main"=singler$singler[[2]]$SingleR.single.main$labels,
                       "cell.names" = rownames(singler$singler[[1]]$SingleR.single$labels))
rownames(singlerDF) = singlerDF[,"cell.names"]
singlerDF = singlerDF[,-which(colnames(singlerDF) == "cell.names")]
EC <- AddMetaData(object = EC,
                  metadata = as.data.frame(singlerDF))
EC <- SetAllIdent(object = EC, id = "singler1main")
##############################
# process singler.color
##############################
singler_colors <- readxl::read_excel("../SingleR/singler.colors.xlsx")
singler_colors = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
length(singler_colors)
length(unique(singler$singler[[1]]$SingleR.single.main$labels))

#' AddMetaColor: prepare meta.data data frame to store color code
#' @mat factor at column 1, rownames is cell.names
#' @colors vector of hex colors
#' @example 
# singlerDF = factor2color(mat = singler$singler[[2]]$SingleR.single$labels,
#                       colors = singler_colors[3:26], col.name = "singler2sub")
AddMetaColor <- function(mat, colors){
        if(class(mat) != "data.frame") mat = as.data.frame(mat)
        mat$cell.names = rownames(mat)
        mat$index <- as.numeric(as.factor(mat[,1]))
        df_colors = data.frame("colors" = colors[1:length(unique(mat$index))],
                               "index" = 1:length(unique(mat$index)))
        mat_colors <- full_join(mat, df_colors, by = "index")
        # remove NA color if there is any
        mat_colors <- mat_colors[(mat_colors$cell.names %in% rownames(mat)),]
        rownames(mat_colors) = mat_colors$cell.names
        df_colors <- data.frame("colors" = mat_colors$colors,
                                row.names = mat_colors$cell.names)
        
        return(df_colors)
}
singlerDF = AddMetaColor(mat = singler$singler[[1]]$SingleR.single.main$labels,
                         colors = singler_colors[4:14])
EC <- AddMetaData(object = EC, metadata = singlerDF)

##############################
# draw tsne plot
##############################

#' ExtractMetaColor: extract color code from meta.data
#' @object seurat object, meta.data slot must have "color"
#' @label corresponding label
#' @example 
# singlerDF = factor2color(mat = singler$singler[[2]]$SingleR.single$labels,
#                       colors = singler_colors[3:26], col.name = "singler2sub")
ExtractMetaColor <- function(object, label = "singler1main"){
        meta.data =object@meta.data
        meta.data = meta.data[,c(label,"colors")]
        meta.data$index <- as.numeric(as.factor(meta.data[,1]))
        df_colors = meta.data[!duplicated(meta.data$index),]
        df_colors = df_colors[order(df_colors$index),]
        
        return(as.character(df_colors$colors))
}

df_colors <- ExtractMetaColor(object = EC, label = "singler1main")
p3 <- DimPlot.1(object = EC, reduction.use = "tsne", dim.1 = 1, dim.2 = 2, 
                group.by = "ident", do.return = TRUE,
                cols.use = df_colors, no.legend = F,
                do.label =F,label.size=4, label.repel = T,force=1)+
        ggtitle("Supervised sub cell type labeling by immune gene")+
        theme(text = element_text(size=15),
              plot.title = element_text(hjust = 0.5))

jpeg(paste0(path,"PlotTsne_main.jpeg"), units="in", width=10, height=7,
     res=600)
print(p3)
dev.off()

save(EC,file="./output/EC_20181022.RData")
##############################
# subset Seurat
###############################
table(EC@meta.data$orig.ident)
table(EC@ident)

p4 <- SplitTSNEPlot(EC,group.by = "ident", split.by = "orig.ident",
                    select.plots = c(2,1),
                    no.legend = F,do.label =T,label.size=3,
                    return.plots =T, label.repel = T,force=2)
jpeg(paste0(path,"/splitetsne.jpeg"), units="in", width=10, height=7,
     res=600)
do.call(plot_grid, p4)
dev.off()