########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(gplots)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)


# 3.1.1 load data
# Rename ident
(load(file = "data/LynchSyndrome_6_20190802.Rda"))
object@meta.data$orig.ident %<>% gsub("Aspirin","Naproxen",.)
object@meta.data$celltype.conditions %<>% gsub("Aspirin","Naproxen",.)
object@meta.data$conditions %<>% gsub("Aspirin","Naproxen",.)

NewNames = gsub("Aspirin","Naproxen",colnames(object))
object %<>% RenameCells(new.names = NewNames)
rownames(object@reductions$tsne@cell.embeddings) = colnames(object)
gsub("_.*","",rownames(object@reductions$tsne@cell.embeddings)) %>% table

Idents(object) = "integrated_snn_res.0.6"
table(Idents(object))
object %<>% sortIdent(numeric = T)

TSNEPlot.1(object)
save(object, file = "data/LynchSyndrome_6_20190802.Rda")
