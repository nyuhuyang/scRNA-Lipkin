library(SingleR)
source("../R/SingleR_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file = "data/object_data_mm10_6_20190802.Rda"))
(load(file='data/ref_GSE118974_blueprint_encode.RData'))

singler <- CreateBigSingleRObject.1(object_data, annot = NULL, project.name="EC-OG-5903",
                                    N = 5000, min.genes = 500, technology = "10X",
                                    species = "Human", citation = "", 
                                    ref.list = list(ref),
                                    normalize.gene.length = F, variable.genes = "de", 
                                    fine.tune = F,
                                    reduce.file.size = T, do.signatures = F, do.main.types = T,
                                    temp.dir = getwd(), numCores = SingleR.numCores)
save(singler,file="output/singler_F_LynchSyndrome_6_20190816.Rda")
  
