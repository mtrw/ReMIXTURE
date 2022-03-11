setwd("\\\\filer-5.ipk-gatersleben.de\\agruppen\\GGR\\wallace\\workspace\\eggplant_remixture")

#devtools::install_github("mtrw/ReMIXTURE")
#library(ReMIXTURE)
#library(data.table)
#?ReMIXTURE



#Set up DM ############################################################
dm <- readRDS("data/ibs_TARGET_SNPs_melongena.rds")$ibs
dm <- 1-dm
ind_regions <- fread("./data/Melongena_origin_updated.txt",col.names=c("sample_id","region"))
tmp <- data.table(sample_id=readRDS("data/ibs_TARGET_SNPs_melongena.rds")$sample.id)
regions <- ind_regions[tmp,on="sample_id"]
#regions[,.N,by=.(region)]
colnames(dm) <- rownames(dm) <- regions$region
selectinds <- !is.na(rownames(dm)) & !rownames(dm) %in% c("Unknown","E_S_America","W_S_America","Meso_America-Caribbean","Oceania")
dm <- dm[selectinds,selectinds]
dm <- dm[order(rownames(dm)),order(rownames(dm))]
dm[1:5,1:5]
#######################################################################

#Set up region table  #################################################
rt <- fread("data/region_info.csv")
rt
############################################################

rx <- ReMIXTURE$new(distance_matrix = dm,region_table = rt)

rx$run(iterations = 1000,subsample_proportion = 0.8,h_cutoff = 0.13)

rx$plot_heatmap()


rx$plot_maps()
