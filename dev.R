
devtools::load_all()
devtools::document()
devtools::build_vignettes()
vignette("ReMIXTURE")

setwd("/data/gpfs/projects/punim1869/users/mrabanuswall/workspace/vitis_RM")
source("scripts/0_setup.R")
library(ReMIXTURE)


dm <- readRDS(dmFname)
regTab <- readRDS(rtFname)

rm <- ReMIXTURE$new(
  distance_matrix = dm,
  region_table = regTab
)


mds <- rm$plot_MDS(doPlot = TRUE)

# Check goot range of cutoffs to try
rm$plot_distance_densities()

rm$run(
  iterations = 1,
  h_cutoffs = seq( .05, 0.7, l=14),
  subsample_proportions = 0.8#,
  #diagnosticPlotMDSclusters = TRUE
)

rm$run_results

rm$plot_h_optimisation()

rm$region_table

rm$plot_maps(run=4,focalRegion = "Greece",range_lon = c(-10,60),range_lat = c(20,60),width_lims = c(.5,3),alpha_lims = c(.1,3))
rm$plot_maps(run=4,range_lon = c(-10,60),range_lat = c(20,60),width_lims = c(.5,3),alpha_lims = c(.1,3))


#










setwd("/data/gpfs/projects/punim1869/users/mrabanuswall/my_packages/ReMIXTURE/")
devtools::document()
#devtools::install_github("https://github.com/mtrw/ReMIXTURE")
library(ReMIXTURE)

?ReMIXTURE

rm$penis()

# Params
RmHumansRunFname <- "/data/gpfs/projects/punim1869/users/amadhusudans/Trial_genomes/filt_Trial_output/barley_remixture/rmbarley_Obj.rds"
tmp <- readRDS(RmHumansRunFname)
dm <- tmp$distance_matrix
diag(dm) <- 1.0
dm[1:5,1:5]

my_analysis <- ReMIXTURE$new(
  distance_matrix = 1-dm,
  region_table = tmp$region_table
)

my_analysis$run( iterations=5 , subsample_proportions=c(0.95) , h_cutoffs=seq(0,.14,l=30) )

my_analysis$plot_distance_densities(samePlot = TRUE)
par(mfrow=c(1,1))

my_analysis$run_results

my_analysis$plot_h_optimisation()

my_analysis$plot_heatmaps()

my_analysis$plot_clustercounts()

my_analysis$region_table

#debugonce(my_analysis$plot_maps)
my_analysis$region_table
my_analysis$plot_maps(
  focalRegion = "Ethiopia",
  run = 12,
  #range_lon = c(-23.0,60.0),
  #range_lat = c(20,60),
  width_lims = c(5,20),
  alpha_lims = c(.05,1),
  projection = winkelIII
  #curvature_matrix = cm
)

#



# DEV:
# auto make curve matrix
# vignette
# fix map edges in plot with 180 at edges










































setwd("/data/gpfs/projects/punim1869/users/mrabanuswall/workspace/vitis_RM")
source("scripts/0_setup.R")
library(ReMIXTURE)

dm <- readRDS(dmFname)
regTab <- readRDS(rtFname)

rm <- ReMIXTURE$new(
  distance_matrix = dm,
  region_table = regTab
)

debugonce(rm$plot_MDS)
mds <- rm$plot_MDS()

# Check goot range of cutoffs to try
rm$plot_distance_densities()
### RUN ####
debugonce(rm$run)
rm$run(
  iterations = 1,
  h_cutoffs = seq(.1,.7,l=30),
  subsample_proportions = 0.8,
  diagnosticPlotMDSclusters = TRUE
)
rm$plot_h_optimisation()


# #### Run Diagnostics ###
rm$plot_clustercounts()


rm$plot_heatmaps()
#### Curviness ####
{
  set.seed(2)
  cm <- matrix(abs(rnorm(nrow(regTab)**2))/2,nrow=nrow(regTab),dimnames = list(regTab$region,regTab$region))
}

# Plot maps
debugonce(rm$plot_maps)
rm$plot_maps(
  #focalRegion = "France",
  run = 7,
  range_lon = c(-23.0,60.0),
  range_lat = c(20,60),
  width_lims = c(.4,3),
  alpha_lims = c(0,1),
  curvature_matrix = cm,
  projection = eckertIV
)






#### Plot ####
#debugonce(rm$plot_maps)
dev.off()
par(mfrow=c(4,5),mar=c(0,0,1,0))
rm$plot_maps(
  #focalRegion = "France",
  run = 7,
  range_lon = c(-23.0,60.0),
  range_lat = c(20,60),
  width_lims = c(.4,3),
  alpha_lims = c(0,1),
  curvature_matrix = cm
)


setwd("/data/gpfs/projects/punim1869/shared_projects/ReMIXTURE")
devtools::document()
devtools::load_all()

