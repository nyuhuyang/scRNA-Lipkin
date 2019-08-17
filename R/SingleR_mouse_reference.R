############################################
# combine mouse.rnaseq and 
############################################
library(SingleR)
library(genefilter)
library(dplyr)
library(magrittr)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
####functions===========

attach(mouse.rnaseq)
table(mouse.rnaseq$types)
dim(mouse.rnaseq$data)
head(mouse.rnaseq$data[,1:5])
anyNA(mouse.rnaseq$data)
testMMM(mouse.rnaseq$data)
boxplot(mouse.rnaseq$data, main="mouse.rnaseq")#slow!

# remove low quanlity mouse.rnaseq data
par(mfrow=c(2,1))
hist(colMeans(mouse.rnaseq$data),breaks=ncol(mouse.rnaseq$data))
quantile_75 <- apply(mouse.rnaseq$data,2,function(x) quantile(x,probs =0.75))
hist(quantile_75, breaks=ncol(mouse.rnaseq$data))
rm_samples <- names(quantile_75)[quantile_75<5]
(rm_index <- which(colnames(mouse.rnaseq$data) %in% rm_samples))
mouse.rnaseq_rm <- mouse.rnaseq$data[,-rm_index]
par(mfrow=c(1,1))
boxplot(mouse.rnaseq_rm,main="mouse.rnaseq_rm")#slow!
testMMM(mouse.rnaseq_rm)

remove(mouse.rnaseq_rm);GC()

sort(unique(mouse.rnaseq$main_types))
rm_T = grep("T cells", mouse.rnaseq$main_types)
(rm_index <- unique(c(rm_index,rm_T)))

mouse.rnaseq_rm = mouse.rnaseq$data[,-rm_index]
types = FineTune(as.character(mouse.rnaseq$types[-rm_index]),
                       main.type = FALSE)
main_types = FineTune(as.character(mouse.rnaseq$main_types[-rm_index]),
                            main.type = TRUE)
mouse.rnaseq = CreateSinglerReference(name = 'mouse.rnaseq',
                                          expr = mouse.rnaseq_rm,
                                          types = types, 
                                          main_types = main_types)
save(mouse.rnaseq,file='../SingleR/data/mouse.rnaseq.RData')

#===========================================
# check GSE109125 data
#===========================================

GSE109125 <- read.csv("data/GSE109125_Normalized_Gene_count_table.csv")
head(GSE109125[,1:3])
# remove NA rows
GSE109125 <- GSE109125[!apply(GSE109125,1, function(x) all(is.na(x))),]

#remove duplicate rownames with lower rowsumns
#' @param mat input as data.frame with gene name
#' @export mat matrix with gene as rownames, no duplicated genes
RemoveDup <- function(mat){
        gene_id <- as.matrix(mat[,1])
        mat <- mat[,-1]
        if(!is.matrix(mat)) mat <- sapply(mat,as.numeric)
        rownames(mat) <- 1:nrow(mat)
        mat[is.na(mat)] = 0
        mat <- cbind(mat, "rowSums" = rowSums(mat))
        mat <- mat[order(mat[,"rowSums"],decreasing = T),]
        gene_id <- gene_id[as.numeric(rownames(mat))]
        remove_index <- duplicated(gene_id)
        mat <- mat[!remove_index,]
        rownames(mat) <- gene_id[!remove_index]
        return(mat[,-ncol(mat)])
}
GSE109125 <- RemoveDup(GSE109125)
dim(GSE109125)
head(GSE109125[,1:3])
GSE109125 = log2(GSE109125+1)
#testMMM(GSE109125)

# GSE109125 cell types======

GSE109125_label <- read.delim2("data/GSE109125_series_matrix.txt",stringsAsFactors = F)
rownames(GSE109125_label) = GSE109125_label[,1]
GSE109125_label = GSE109125_label[,-1]
head(GSE109125_label[,1:3])
GSE109125_label = as.data.frame(t(GSE109125_label))
GSE109125_label[,"main_type"] = GSE109125_label[,"Sample_characteristics_cell type"]
GSE109125_label[,"main_type"] = gsub("gd T cells","T_cells:gd", GSE109125_label[,"main_type"])
GSE109125_label[,"type"] = GSE109125_label[,"main_type"]
GSE109125_label[,"main_type"] = gsub(".*T cells","T cells", GSE109125_label[,"main_type"])
GSE109125_label[,"rownames"] = rownames(GSE109125_label)
GSE109125_label = apply(GSE109125_label,2,as.character)
rownames(GSE109125_label) = GSE109125_label[,"rownames"]
GSE109125_label[,"rownames"] = sub('.[1-5]$', '', rownames(GSE109125_label))

GSE109125_label[grepl("NKT",GSE109125_label[,"rownames"]),"type"] = "T_cells:NKT"
GSE109125_label[grepl("T\\.4",GSE109125_label[,"rownames"]),"type"] = "T_cells:CD4+"
GSE109125_label[grepl("T\\.8|T8\\.",GSE109125_label[,"rownames"]),"type"] = "T_cells:CD8+"
GSE109125_label[grepl("Treg\\.4",GSE109125_label[,"rownames"]),"type"] = "T_cells:Treg"
GSE109125_label[,"type"] = gsub("ab T cells","T_cells:ab",GSE109125_label[,"type"])

GSE109125_label[,"main_type"] = gsub("T_cells\\:gd","T_cells",GSE109125_label[,"main_type"])

Th <- grepl("Th",GSE109125_label[,"rownames"]) & 
        GSE109125_label[,"type"]  %in% c("T_cells:ab")
GSE109125_label[Th,"type"] = "T_cells:Th"
GSE109125_label[grepl("preT",GSE109125_label[,"rownames"]),"type"] = "T_cells:preT_Th"

GSE109125_label[,"type"] = FineTune(GSE109125_label[,"type"], main.type = FALSE)
GSE109125_label[,"main_type"] = FineTune(GSE109125_label[,"main_type"],main.type = TRUE)

#---- test ------
(colnames <- sub('.[1-5]$', '', colnames(GSE109125)))
df_table <- table(rownames(GSE109125_label), GSE109125_label[,"type"]) %>%
        as.data.frame() %>% spread(Var2, Freq)
rownames(df_table) = df_table$Var1
df_table = df_table[,-1]
df_table[df_table$`ab T cells` >0,]

#================================
# merge MCL and mouse.rnaseq
#================================

dim(mouse.rnaseq$data)
GSE109125_mouse.rnaseq <- merge(GSE109125[,rownames(GSE109125_label)],mouse.rnaseq$data,
                         by="row.names",all=FALSE)
rownames(GSE109125_mouse.rnaseq) = GSE109125_mouse.rnaseq$Row.names
GSE109125_mouse.rnaseq <- GSE109125_mouse.rnaseq[-which(colnames(GSE109125_mouse.rnaseq)=="Row.names")]
testMMM(GSE109125_mouse.rnaseq)

colsum <- colSums(GSE109125_mouse.rnaseq)
scale_factor = median(colsum)
GSE109125_mouse.rnaseq = GSE109125_mouse.rnaseq/colsum * scale_factor
testMMM(GSE109125_mouse.rnaseq)


jpeg(paste0(path,"boxplot_GSE109125_mouse.rnaseq.jpeg"), units="in", width=10, height=7,res=600)
boxplot(GSE109125_mouse.rnaseq) #slow
title(main = "boxplot for GSE109125 + mouse.rnaseq")
dev.off()

# Create Singler Reference=============
ref = CreateSinglerReference(name = 'GSE109125_mouse.rnaseq',
                             expr = as.matrix(GSE109125_mouse.rnaseq), # the expression matrix
                             types = c(GSE109125_label[,"type"],mouse.rnaseq$types), 
                             main_types = c(GSE109125_label[,"main_type"], mouse.rnaseq$main_types))

save(ref,file='data/ref_GSE109125_mouse.rnaseq.RData') # it is best to name the object and the file with the same name.
