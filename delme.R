





#devtools::install_github("https://github.com/mtrw/ReMIXTURE")
library(ReMIXTURE)

#new analysis using example data provided
my_analysis <- ReMIXTURE$new(
  distance_matrix = ReMIXTURE::ReMIXTURE_example_distance_matrix,
  region_table = ReMIXTURE::ReMIXTURE_example_region_positions
)

#A good selection of h-cutoffs should span the lower end of all distance peaks

my_analysis$run(iterations = 80,subsample_proportions = c(0.8),h_cutoffs=seq(.008,.04,l=8))
#my_analysis$run(iterations = 250,subsample_proportions = c(0.8),h_cutoffs=0.013)

my_analysis$run_results$runs[[1]]$correlations %>% pheatmap(cluster_rows = F,cluster_cols = F)



#my_analysis$plot_distance_densities(set_xlims = c(0,0.2))

#Raw results
#my_analysis$run_results

#Tools to assess h-cutoff values
   #Good values are often found between where the median clustercounts are equal, and where the inter-region cluster counts maximise
my_analysis$plot_h_optimisation()
my_analysis$plot_h_optimisation(plot_entropy = T)
#A good value has a "balanced" heatmap, the main features in which are stable at nearby values
my_analysis$plot_heatmaps()
   #A good value does not incur many or any per-region cluster counts near one. A few may be inevitable in some datasets though.
my_analysis$plot_clustercounts()

#Plot the regions on the globe, and add raw plotting data to my_analysis$run_results
my_analysis$plot_maps(run = 1,alpha_norm_per_region = T)
my_analysis$plot_maps(run = 1,alpha_correlation = T ,alpha_norm_per_region = F)

