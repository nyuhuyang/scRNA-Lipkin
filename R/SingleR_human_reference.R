library(SingleR)
library(genefilter)
library(dplyr)
source("../R/Seurat3_functions.R")
source("../R/SingleR_functions.R")
################################################
# filter blueprint_encode data
################################################
attach(blueprint_encode)
unique(blueprint_encode$types)
dim(blueprint_encode$data)
head(blueprint_encode$data[,1:5])
table(is.na(blueprint_encode$data))
blueprint_encode$data[is.na(blueprint_encode$data)] = 0
head(colSums(blueprint_encode$data))
testMMM(blueprint_encode$data)
# remove low quanlity blueprint_encode data
par(mfrow=c(2,1))
hist(colMeans(blueprint_encode$data),breaks=ncol(blueprint_encode$data))
quantile_75 <- apply(blueprint_encode$data,2,function(x) quantile(x,probs =0.75))
hist(quantile_75, breaks=ncol(blueprint_encode$data))
rm_samples <- names(quantile_75)[quantile_75<1]
rm_index <- which(colnames(blueprint_encode$data) %in% rm_samples)
blueprint_encode_rm <- blueprint_encode$data[,-rm_index]
quantile_75_new <- apply(blueprint_encode_rm,2,function(x) quantile(x,probs =0.75))
hist(quantile_75, breaks=ncol(blueprint_encode$data))
hist(quantile_75_new, breaks=ncol(blueprint_encode_rm),xlim = c(0,4.1))

par(mfrow=c(1,1))
boxplot(blueprint_encode_rm)#slow!
title(main="blueprint_encode_rm")

types = FineTune(as.character(blueprint_encode$types[-rm_index]),
                       main.type = FALSE)
main_types = FineTune(as.character(blueprint_encode$main_types[-rm_index]),
                            main.type = TRUE)

Blueprint_encode = CreateSinglerReference(name = 'Blueprint_encode',
                                             expr = blueprint_encode_rm,
                                             types = types, 
                                             main_types = main_types)
remove(blueprint_encode);GC()
################################################
# check GSE118974 data
################################################
GSE118974 <- read.csv("data/GSE118974_Th0_Th17_RNASeq_72h_Processed_Data.csv")
head(GSE118974[,1:3])
# remove NA rows
GSE118974 <- GSE118974[!apply(GSE118974,1, function(x) all(is.na(x))),]

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
GSE118974 <- RemoveDup(GSE118974)
dim(GSE118974)
head(GSE118974[,1:3])
testMMM(GSE118974)

#================================
# merge MCL and blueprint_encode
#================================
GSE118974_types <- paste0("T_cells:",gsub("\\_.*","",colnames(GSE118974)))
GSE118974_blueprint_encode <- merge(GSE118974,blueprint_encode_rm,
                                by="row.names",all=FALSE)
rownames(GSE118974_blueprint_encode) = GSE118974_blueprint_encode$Row.names
GSE118974_blueprint_encode <- GSE118974_blueprint_encode[-which(colnames(GSE118974_blueprint_encode)=="Row.names")]
testMMM(GSE118974_blueprint_encode)

colsum <- colSums(GSE118974_blueprint_encode)
scale_factor = median(colsum)
GSE118974_blueprint_encode = GSE118974_blueprint_encode/colsum * scale_factor
testMMM(GSE118974_blueprint_encode)


jpeg(paste0(path,"boxplot_GSE118974_blueprint_encode.jpeg"), units="in", width=10, height=7,res=600)
boxplot(GSE118974_blueprint_encode) #slow
title(main = "boxplot for GSE118974 + blueprint_encode")
dev.off()

# Create Singler Reference=============
ref = CreateSinglerReference(name = 'GSE118974_blueprint_encode',
                             expr = as.matrix(GSE118974_blueprint_encode), # the expression matrix
                             types = c(GSE118974_types,Blueprint_encode$types), 
                             main_types = c(gsub("\\:.*","",GSE118974_types), Blueprint_encode$main_types))

save(ref,file='data/ref_GSE118974_blueprint_encode.RData') # it is best to name the object and the file with the same name.
