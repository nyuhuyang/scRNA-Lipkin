library(Seurat)
library(SingleR)
library(dplyr)
library(tidyr)
library(kableExtra)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")

#====== 2.1 identify phenotype for each cluster  ==========================================
lnames = load(file="./output/EC_20181022.RData")
lnames

test.markers <- unique(c(Myeloid,Lymphoid))
test.markers <- MouseGenes(EC,"PTPRC")
for(i in 1:length(test.markers)) {
        jpeg(paste0(path,"/markers/",test.markers[i],".jpeg"), units="in", width=10, height=7,
             res=600)
        p1 <- SingleFeaturePlot.1(object = EC, feature = test.markers[i])
        print(p1)
        print(paste0(i,":",length(test.markers)))
        dev.off()
}


#------Myeloid----
megakaryocytes <-  MouseGenes(EC,c("PPBP","GNG11"))
erythrocyte <-  MouseGenes(EC,c("HBA2","HBB"))
MastCells <- MouseGenes(EC,c("Cma1","Mcpt4","Tpsb2","Cpa3"))
Monocytes <-  MouseGenes(EC,c("LYZ","S100A9","CD14","CCL2","FCGR3A","MS4A7","VMO1"))
Macrophages <- MouseGenes(EC,c("LYZ","CD68","MARCO","Emr1"))
DendriticCells <- MouseGenes(EC,c("Itgax","GPR183","CST3","HLA-DQA1","FCER1A","TSPAN13",
                                     "IL3RA","IGJ"))
Myeloid <-  MouseGenes(EC,c(megakaryocytes,erythrocyte,MastCells,
                               Monocytes,Macrophages,DendriticCells))
#------Lymphoid----

# T cell
T_Cell <- MouseGenes(EC,c("CD3G","CD3D","CD2","CD8A","IL2RA",
                             "FOXP3","NCAM1","FCGR3A"))
CD4_Naive_T <- MouseGenes(EC,c("CD4","IL7R","GIMAP5","SELL","IL2RG"))
NK <- MouseGenes(EC,c("NKG7","CCL5"))
# B cell
Immature_B <- MouseGenes(EC,c("MME","CD19","MS4A1","CD34","CR2","CD40","CD24","IL4R"))


Lymphoid <- MouseGenes(EC,c(T_Cell,CD4_Naive_T,NK,Immature_B))
