library(SingleR)
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file = "data/object_data_mm10_6_20190802.Rda"))
attach(mouse.rnaseq)
attach(blueprint_encode)

singler = CreateSinglerObject(object_data, annot = NULL, project.name="EC-OG-5903",
                              min.genes = 500,technology = "10X", species = "Human", citation = "",
                              ref.list = list(blueprint_encode),
                              normalize.gene.length = F, variable.genes = "de",
                              fine.tune = T, do.signatures = F, clusters = NULL)
save(singler,file="output/singler_T_LynchSyndrome_6_20190805.Rda")
  
