########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
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
#======1.1 Load the data files and Set up Seurat object =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/181128_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c("test",paste0("test",1:6)))
df_samples[sample_n,] %>% kable() %>% kable_styling()
table(df_samples$tests);nrow(df_samples)
samples <- df_samples$sample[sample_n]
sample.id <- df_samples$sample[sample_n]
projects <- df_samples$project[sample_n]
tests <- df_samples$tests[sample_n]
conditions <- df_samples$conditions[sample_n]

# check missing data
current <- list.files("data")[!grepl(".Rda|RData",list.files("data"))]
(new_samples <- sample.id[!(samples %in% current)])

# Move files from Download to ./data and rename them
current.folder <- "/Users/yah2014/Downloads"
new.folders <- "./data"
for(new_sample in new_samples){
    list.of.files <- list.files(paste(current.folder, new_sample,"outs","filtered_gene_bc_matrices","mm10",sep = "/")
    new.sample <- sub("_","-",new_sample)
    new.folder <- list.files(paste("./data", new.sample,"outs","filtered_gene_bc_matrices","mm10",sep = "/")
    if(!dir.exists(new.folder)) dir.create(new.folder, recursive = T)
    # copy the files to the new folder
    file.copy(list.of.files, new.folder)
}

## Load the dataset
EC_raw <- list()
EC_Seurat <- list()
for(i in 1:length(samples)){
        EC_raw[[i]] <- Read10X(data.dir = paste0("./data/",
                                                   samples[i],"/outs/filtered_gene_bc_matrices/mm10/"))
        colnames(EC_raw[[i]]) <- paste0(samples[i],"_",colnames(EC_raw[[i]]))
        EC_Seurat[[i]] <- CreateSeuratObject(EC_raw[[i]],
                                               min.cells = 3,
                                               min.genes = 200,
                                               names.delim = "_",
                                               project = projects[i],)
        EC_Seurat[[i]]@meta.data$conditions <- conditions[i]
}
EC <- Reduce(function(x, y) MergeSeurat(x, y, do.normalize = F), EC_Seurat)
remove(EC_raw,EC_Seurat);GC()
EC <- FilterCells(EC, subset.names = "nGene",
                    low.thresholds = 200,
                    high.thresholds = Inf) %>%
        NormalizeData() %>%
        ScaleData(display.progress = FALSE) %>%
        FindVariableGenes(do.plot = FALSE, display.progress = FALSE)
save(EC, file = "./data/EC_20181022.Rda")

#======1.2 QC, pre-processing and normalizing the data=========================
# 1.2.1 Calculate median UMI per cell
EC_raw_data <- as.matrix(x = EC@raw.data)
mean(colSums(EC_raw_data))
median(colSums(EC_raw_data))
min(colSums(EC_raw_data))
remove(EC_raw_data);GC()

# 1.2.3 calculate mitochondria percentage
mito.genes <- grep(pattern = "^mt-", x = rownames(x = EC@data), value = TRUE)
percent.mito <- Matrix::colSums(EC@raw.data[mito.genes, ])/Matrix::colSums(EC@raw.data)
EC <- AddMetaData(object = EC, metadata = percent.mito, col.name = "percent.mito")

g1 <- VlnPlot(object = EC, features.plot = c("nGene", "nUMI", "percent.mito"), 
              nCol = 1,point.size.use = 0.2,
              x.lab.rot = T, do.return = T,return.plotlist =T)

EC <- FilterCells(object = EC, subset.names = c("nGene", "nUMI", "percent.mito"),
                    low.thresholds = c(500,3000,0), 
                    high.thresholds = c(Inf,Inf,0.1))

g2 <- VlnPlot(object = EC, features.plot = c("nGene", "nUMI", "percent.mito"), 
              nCol = 1,point.size.use = 0.2,
              x.lab.rot = T, do.return = T,return.plotlist =T)

jpeg(paste0(path,"/S1_mito.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[3]]+ggtitle("mito % in raw data")+ 
                        ylim(c(0,0.2)),
                g2[[3]]+ggtitle("mito % after filteration")+ 
                        ylim(c(0,0.2))))
dev.off()

jpeg(paste0(path,"/S1_GenePlot.jpeg"), units="in", width=10, height=7,res=600)
par(mfrow = c(2, 1))
GenePlot(object = EC, gene1 = "nUMI", gene2 = "percent.mito",use.raw = T)
GenePlot(object = EC, gene1 = "nUMI", gene2 = "nGene",use.raw = T)
dev.off()


# After removing unwanted cells from the dataset, the next step is to normalize the data.
EC <- NormalizeData(object = EC, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
EC <- FindVariableGenes(object = EC, mean.function = ExpMean, 
                          dispersion.function = LogVMR, do.plot = FALSE, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(EC@var.genes)
#======1.3 1st run of pca-tsne  =========================
EC <- ScaleData(object = EC) %>%
        RunPCA() %>%
        FindClusters(dims.use = 1:20, force.recalc = T, print.output = FALSE) %>%
        RunTSNE()

p1 <- TSNEPlot(object = EC, do.label = F, group.by = "orig.ident", 
               do.return = TRUE, no.legend = F, #colors.use = singler.colors,
               pt.size = 1,label.size = 8 )+
        ggtitle("Original")+
        theme(text = element_text(size=15),							
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

save(EC, file = "./data/EC_20180825.Rda")
Iname = load("./data/EC_20180825.Rda")

#======1.4 Add Cell-cycle score =========================
# Read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "../R/seurat_resources/regev_lab_cell_cycle_genes.txt")
# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- MouseGenes(EC,cc.genes[1:43])
g2m.genes <- MouseGenes(EC,cc.genes[44:97])
# Assign Cell-Cycle Scores
EC <- CellCycleScoring(object = EC, s.genes = s.genes, g2m.genes = g2m.genes, 
                         set.ident = TRUE)
# Visualize the distribution of cell cycle markers across
RidgePlot(object = EC, features.plot = MouseGenes(EC,c("PCNA", "TOP2A", "MCM6", "MKI67")), 
          nCol = 2)
# regressing out the difference between the G2M and S phase scores
EC@meta.data$CC.Difference <- EC@meta.data$S.Score - EC@meta.data$G2M.Score
# view cell cycle scores and phase assignments
head(x = EC@meta.data)

#======1.5 Add batch id =========================
batchname = EC@meta.data$orig.ident
batch.effect = as.numeric(factor(batchname,levels = samples))
names(batch.effect) = rownames(EC@meta.data)
EC <- AddMetaData(object = EC, metadata = batch.effect, col.name = "batch.effect")
table(EC@meta.data$batch.effect)
head(x = EC@meta.data)

#======1.6 batch-correct using ComBat =========================
jpeg(paste0(path,"S1_nUMI~.jpeg"), units="in", width=10, height=7,res=600)

SingleFeaturePlot.1(EC,"nUMI",threshold=15000)
SingleFeaturePlot.1(EC,"batch.effect",threshold=2.0)
SingleFeaturePlot.1(EC,"percent.mito",threshold=0.05)
SingleFeaturePlot.1(EC,"CC.Difference",threshold=0.05)

EC <- ScaleData(object = EC, 
                model.use = "linear", do.par=T, do.center = T, do.scale = T,
                vars.to.regress = c("nUMI","percent.mito","batch.effect"),
                display.progress = T)
#======1.7 Performing MNN-based correction =========================
#https://bioconductor.org/packages/3.8/workflows/vignettes/simpleSingleCell/inst/doc/work-5-mnn.html#4_performing_mnn-based_correction
set.seed(100)
original <- lapply(samples, function(x) EC@scale.data[EC@var.genes,
                                                      (EC@meta.data$orig.ident %in% x)])
mnn.out <- do.call(fastMNN, c(original, list(k=20, d=50, auto.order=T,
                                             approximate=TRUE)))
dim(mnn.out$corrected)
rownames(mnn.out$corrected) = EC@cell.names
colnames(mnn.out$corrected) = paste0("MNN_",1:ncol(mnn.out$corrected))
#Storing a new MNN
EC <- SetDimReduction(object = EC, reduction.type = "MNN", slot = "cell.embeddings",
                              new.data = mnn.out$corrected)
EC <- SetDimReduction(object = EC, reduction.type = "MNN", slot = "key", 
                              new.data = "MNN_")
remove(original);GC()
EC <- SetAllIdent(EC,id = "orig.ident")
DimPlot(object = EC, reduction.use = "MNN", pt.size = 0.5)

#======1.7 unsupervised clustering based on MNN =========================
EC <- RunPCA(object = EC, pc.genes = EC@var.genes, pcs.compute = 100, 
                     do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = EC)
PCElbowPlot(object = EC, num.pc = 100)
PCHeatmap(EC, pc.use = c(1:3, 25:30), cells.use = 500, do.balanced = TRUE)

DimElbowPlot.1(object = EC, reduction.type = "MNN", 
               dims.plot = 50,slot = "cell.embeddings")

EC <- RunTSNE(object = EC, reduction.use = "MNN", dims.use = 1:50, 
                      do.fast = TRUE, perplexity= 30)

EC <- FindClusters(object = EC, reduction.type = "MNN", 
                           dims.use = 1:50, resolution = 0.6, 
                           k.param = 30,force.recalc = T,
                           save.SNN = TRUE, n.start = 100, nn.eps = 0, print.output = FALSE)

p2 <- TSNEPlot.1(object = EC, do.label = F, group.by = "orig.ident", 
                 do.return = TRUE, no.legend = T, 
                 pt.size = 1,label.size = 4 )+
        ggtitle("Corrected")+
        theme(text = element_text(size=15),							
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"remove_batch.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1 + theme(text = element_text(size=15),
                legend.position="none",
                plot.title = element_text(hjust = 0.5,size = 18, face = "bold")),
          p2)
dev.off()

jpeg(paste0(path,"tsneplot.jpeg"), units="in", width=10, height=7,res=600)
p2
dev.off()

save(EC, file = "./data/EC_20181022.Rda")

EC.subset <- SplitSeurat(EC, split.by = "orig.ident")
gg1 <- SingleFeaturePlot.1(object = EC.subset[[1]], "nUMI",threshold=15000)
gg2 <- SingleFeaturePlot.1(object = EC.subset[[2]], "nUMI",threshold=15000)

jpeg(paste0(path,"S1_nUMI~.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(gg1,gg2)
dev.off()
