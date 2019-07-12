########################################################################
#
#  0 setup environment, install libraries if nLynchessary, load libraries
# 
# ######################################################################

library(Seurat)
library(magrittr)
library(harmony)
library(dplyr)
library(kableExtra)
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("./data/")) dir.create("data")
########################################################################
#
#  1 Data preprocessing
# 
# ######################################################################
#======1.1 Load the data files and Set up Seurat objLyncht =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/181128_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c("test",paste0("test",1:6)))
df_samples[sample_n,] %>% kable() %>% kable_styling()
table(df_samples$tests);nrow(df_samples)
samples <- df_samples$sample[sample_n]
sample.id <- df_samples$sample.id[sample_n]
projects <- df_samples$project[sample_n]
tests <- df_samples$tests[sample_n]
conditions <- df_samples$conditions[sample_n]

# check missing data
current <- list.files("data")[!grepl(".Rda|RData",list.files("data"))]
(new_samples <- sample.id[!(samples %in% current)])

# Move files from Download to ./data and rename them
current.folder <- "/Users/yah2014/Downloads"
species = "mm10"
for(new_sample in new_samples){
        old.pth  <- paste(current.folder, new_sample,"outs",
                          "filtered_gene_bc_matrices",species,
                          sep = "/")
        list.of.files <- list.files(old.pth)
        new.sample <- gsub("_","-",new_sample)
        new.folder <- paste("./data", new.sample,"outs",
                            "filtered_gene_bc_matrices",
                            species,sep = "/")
        if(!dir.exists(new.folder)) dir.create(new.folder, recursive = T)
        # copy the files to the new folder
        file.copy(paste(old.pth, list.of.files, sep = "/"), new.folder)
}

## Load the dataset
Lynch_raw <- list()
Lynch_Seurat <- list()
for(i in 1:length(samples)){
        Lynch_raw[[i]] <- Read10X(data.dir = paste0("./data/",
                                                   samples[i],"/outs/filtered_gene_bc_matrices/",species,"/"))
        colnames(Lynch_raw[[i]]) <- paste0(samples[i],"_",colnames(Lynch_raw[[i]]))
        Lynch_Seurat[[i]] <- CreateSeuratObject(Lynch_raw[[i]],
                                               min.cells = 3,
                                               min.genes = 100,
                                               names.delim = "_",
                                               projects = projects[i],)
        Lynch_Seurat[[i]]@meta.data$conditions <- conditions[i]
}

#======1.1.2 QC before merge =========================
cell.number <- sapply(Lynch_Seurat, function(x) length(x@cell.names))
QC_list <- lapply(Lynch_Seurat, function(x) as.matrix(x = x@raw.data))
median.nUMI <- sapply(QC_list, function(x) median(colSums(x)))
median.nGene <- sapply(QC_list, function(x) median(apply(x,2,function(y) sum(length(y[y>0])))))

min.nUMI <- sapply(QC_list, function(x) min(colSums(x)))
min.nGene <- sapply(QC_list, function(x) min(apply(x,2,function(y) sum(length(y[y>0])))))

QC.list <- cbind(df_samples,cell.number, median.nUMI,median.nGene,min.nUMI,min.nGene,
                 row.names = samples)
write.csv(QC.list,paste0(path,"QC_list.csv"))
QC.list %>% kable() %>% kable_styling()
remove(QC_list,median.nUMI,median.nGene,min.nUMI,min.nGene,QC.list);GC()

#========1.1.3 merge ===================================
Lynch <- Reduce(function(x, y) MergeSeurat(x, y, do.normalize = F), Lynch_Seurat)
remove(Lynch_raw,Lynch_Seurat);GC()

mito.genes <- grep(pattern = "^mt-", x = rownames(x = Lynch@data), value = TRUE)
percent.mito <- Matrix::colSums(Lynch@raw.data[mito.genes, ])/Matrix::colSums(Lynch@raw.data)
Lynch <- AddMetaData(object = Lynch, metadata = percent.mito, col.name = "percent.mito")

Lynch@ident = factor(Lynch@ident,levels = samples)

g1 <- VlnPlot(object = Lynch, features.plot = c("nGene", "nUMI", "percent.mito"),
              nCol = 1,point.size.use = 0.2,,size.x.use = 10, group.by = "ident",
              x.lab.rot = T, do.return = T,return.plotlist =T)

save(g1,file="./data/g1_14_20181207.Rda")

#======1.2 load  SingleCellExperiment =========================
(load(file = "./data/sce_14_20181207.Rda"))
names(sce_list)
Lynch_Seurat <- lapply(sce_list, as.seurat) %>%
        lapply(NormalizeData) %>%
        #lapply(ScaleData) %>%
        lapply(FindVariableGenes, do.plot = FALSE)

for(i in 1:length(samples)){
        Lynch_Seurat[[i]]@meta.data$tests <- tests[i]
        Lynch_Seurat[[i]]@meta.data$conditions <- conditions[i]
}
# we will take the union of the top 1k variable genes in each dataset for alignment
genes.use <- Lynch_Seurat %>% 
        lapply(function(object) head(rownames(object@hvg.info), 800)) %>%
        unlist %>% unique
length(genes.use)

#========1.3 merge ===================================
Lynch <- Reduce(function(x, y) MergeSeurat(x, y, do.normalize = F), Lynch_Seurat)
Lynch@var.genes = genes.use
remove(sce_list,Lynch_Seurat);GC()

Lynch = SetAllIdent(Lynch, id = "orig.ident")
#======1.4 mito, QC, filteration =========================
mito.genes <- grep(pattern = "^mt-", x = rownames(x = Lynch@data), value = TRUE)
percent.mito <- Matrix::colSums(Lynch@raw.data[mito.genes, ])/Matrix::colSums(Lynch@raw.data)
Lynch <- AddMetaData(object = Lynch, metadata = percent.mito, col.name = "percent.mito")

(load(file = "./data/g1_14_20181207.Rda"))

Lynch <- FilterCells(object = Lynch, subset.names = c("nGene","nUMI","percent.mito"),
                   low.thresholds = c(50,100, -Inf), 
                   high.thresholds = c(Inf,Inf, 0.5))

Lynch@ident = factor(Lynch@ident,levels = samples)
g2 <- VlnPlot(object = Lynch, features.plot = c("nGene", "nUMI", "percent.mito"), 
              nCol = 1,point.size.use = 0.2,size.x.use = 10, group.by = "ident",
              x.lab.rot = T, do.return = T,return.plotlist =T)
save(g2,file = "./data/g2_14_20181207.Rda")
jpeg(paste0(path,"/S1_nGene.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[1]]+ggtitle("nGene in raw data")+ 
                        scale_y_log10(limits = c(100,10000)),
                g2[[1]]+ggtitle("nGene after filteration")+ 
                        scale_y_log10(limits = c(100,10000))))
dev.off()
jpeg(paste0(path,"/S1_nUMI.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[2]]+ggtitle("nUMI in raw data")+ 
                        scale_y_log10(limits = c(500,100000)),
                g2[[2]]+ggtitle("nUMI after filteration")+ 
                        scale_y_log10(limits = c(500,100000))))
dev.off()
jpeg(paste0(path,"/S1_mito.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[3]]+ggtitle("mito % in raw data")+ 
                        ylim(c(0,0.5)),
                g2[[3]]+ggtitle("mito % after filteration")+ 
                        ylim(c(0,0.5))))
dev.off()

#======1.5 Add Cell-cycle score =========================
# Read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "./data/seurat_resources/regev_lab_cell_cycle_genes.txt")
s.genes <- HumanGenes(Lynch,cc.genes[1:43])
g2m.genes <- HumanGenes(Lynch,cc.genes[44:97])
Lynch <- CellCycleScoring(object = Lynch, s.genes = s.genes, g2m.genes = g2m.genes, 
                        set.ident = TRUE)
RidgePlot(object = Lynch, features.plot = HumanGenes(Lynch,c("CCND1","CDK4","CCND2","CDK6","CCND3","RB1")), 
          nCol = 2)
Lynch@meta.data$CC.Difference <- Lynch@meta.data$S.Score - Lynch@meta.data$G2M.Score
Lynch@meta.data$S.Score = Lynch@meta.data$S.Score - min(Lynch@meta.data$S.Score)
Lynch@meta.data$G2M.Score = Lynch@meta.data$G2M.Score - min(Lynch@meta.data$G2M.Score)
head(x = Lynch@meta.data)

#======1.5 loom pca=======================
Lynch <- NormalizeData(object = Lynch)
jpeg(paste0(path,"/S1_dispersion.jpeg"), units="in", width=10, height=7,res=600)
Lynch <- FindVariableGenes(object = Lynch, mean.function = ExpMean, 
                         dispersion.function = LogVMR, do.plot = T, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 1)
dev.off()
length(Lynch@var.genes)

# Convert from Seurat to loom Convert takes and object in 'from', a name of
# a class in 'to', and, for conversions to loom, a filename

Lynch %<>% ScaleData %>%
        RunPCA(pc.genes = Lynch@var.genes, pcs.compute = 50, do.print = F)

jpeg(paste0(path,"/S1_DimElbowPlot_pca.jpeg"), units="in", width=10, height=7,res=600)
DimElbowPlot(Lynch, reduction.type = "pca", dims.plot = 50)
dev.off()

jpeg(paste0(path,"/S1_PCHeatmap.jpeg"), units="in", width=10, height=7,res=600)
PCHeatmap(Lynch, pc.use = c(1:3, 45:50), cells.use = 500, do.balanced = TRUE)
dev.off()

saveRDS(Lynch@scale.data, file = "./data/Lynch.scale.data_23_20181205.Rda")

#Lynch@scale.data = readRDS("./data/Lynch.scale.data_23_20181205.Rda")
#Lynch <- RunPCA(Lynch, pc.genes = Lynch@var.genes, pcs.compute = 50, do.print = F)
GC()

#======1.6 RunHarmony=======================
jpeg(paste0(path,"/S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(Lynch %<>% RunHarmony("orig.ident", dims.use = 1:50,
                                theta = 2, plot_convergence = TRUE,
                                nclust = 50, max.iter.cluster = 100))
dev.off()

Lynch@ident %<>% factor(levels = samples)
p1 <- DimPlot(object = Lynch, reduction.use = "harmony", pt.size = .1, group.by = "orig.ident", do.return = T)
p2 <- VlnPlot(object = Lynch, features.plot = "Harmony1", group.by = "orig.ident", do.return = TRUE,
              x.lab.rot = T)
jpeg(paste0(path,"/S1_Harmony_vplot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1,p2)
dev.off()

jpeg(paste0(path,"/S1_Harmony_DimHeatmap.jpeg"), units="in", width=10, height=7,res=600)
DimHeatmap(object = Lynch, reduction.type = "harmony", cells.use = 500, 
           dim.use = c(1:3,48:50), do.balanced = TRUE)
dev.off()

#========1.6 Seurat tSNE Functions for Integrated Analysis Using Harmony Results=======
system.time({
        Lynch %<>% RunTSNE(reduction.use = "harmony", dims.use = 1:50, do.fast = TRUE)
        Lynch %<>% FindClusters(reduction.type = "harmony", resolution = 0.6, dims.use = 1:50,
                              save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                              force.recalc = TRUE, print.output = FALSE)
})

p3 <- TSNEPlot(Lynch, do.return = T, pt.size = 0.5, group.by = "orig.ident")
np4 <- TSNEPlot(Lynch, do.label = T, do.return = T, pt.size = 0.5)
jpeg(paste0(path,"/S1_Harmony_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p3, p4)
dev.off()

g_Harmony <- TSNEPlot.1(object = Lynch, do.label = F, group.by = "ident", 
                        do.return = TRUE, no.legend = T, 
                        colors.use = ExtractMetaColor(Lynch),
                        pt.size = 1,label.size = 6 )+
        ggtitle("Tsne plot of all cell types")+
        theme(text = element_text(size=15),							
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"/TSNEplot-Harmony.jpeg"), units="in", width=10, height=7,res=600)
print(g_Harmony)
dev.off()

save(Lynch, file = "./data/Lynch_Harmony_12_20181121.Rda")
