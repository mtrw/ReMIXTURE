setwd("/data/gpfs/projects/punim1869/users/mrabanuswall/my_packages/ReMIXTURE")
# usethis::use_package("data.table")
# usethis::use_package("magrittr")
# usethis::use_package("ggplot2")
# usethis::use_package("truncnorm")

unloadNamespace("ReMIXTURE")

devtools::document()
devtools::install()

?ReMIXTURE

require(ReMIXTURE)

rm <- ReMIXTURE$new(
  distance_matrix = ReMIXTURE_example_distance_matrix,
  region_table = ReMIXTURE_example_region_table
)

debugonce(rm$plot_MDS)
rm$plot_MDS()

debugonce(rm$run)
rm$run()
rm$run(diagnosticPlotMDSclusters=T)

# rm$plot_distance_densities(samePlot=T)
# par(mfrow=c(1,1)) # reset graphics layout
rm$plot_MDS(xlim=c(-0.01,0.01),ylim=c(-0.01,0.01),doPlot = F)
rm$plot_MDS(xlim=c(-0.01,0.01),ylim=c(-0.01,0.01))


rm$plot_h_optimisation()
rm$plot_heatmaps()
rm$plot_maps(run = 4,focalRegion = "Africa",width_max = 10)
rm$plot_heatmaps()

# Default run
# A run involves taking a subsample of all the samples in the distance matrix, such that every region has the same number of samples included. This subsample is heirarchically clustered (imagine creatting a dentrogram and cutting it at a certain height). Clusters are then counted. The number of clusters in which a region appears is a proxy for that region's total genetic diversity. The number of clusters in which *only* one region appears is a proxy for the diversity unique to that region. The number of clusters in which members of a pair of clusters appears is a proxy for the diversity that is overlapping between those two regions. the concept is similar to the way a Venn diagram might work. This is done `iterations` times, and all these various counts are averaged at the end.

# How is it decided where the tree is cut? This value 'H' can be set several ways. By default FINISH ME!!! (and eventually turn me into a vignette)
rm$run(diagnosticPlotMDSclusters = T)
rm$plot_maps()









private <- list()
private$results$runs <- rm$run_results[[1]]


