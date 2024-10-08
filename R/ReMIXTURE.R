#' ReMIXTURE: Regionwise similarity analysis using a resampling/heirarchical clustering based method.
#'
#'
#' @section Warning:
#' This is under development, in early stages.
#' This version is significantly updated and differs non-trivially from the version presented in the paper Tripodi & Rabanus-Wallace et. al. 2021 about peppers. If you wish to use this method, please contact me to talk about how to cite it. A publication is coming soon.
#'
#' @examples
#'
#' #installation
#' devtools::install_github("https://github.com/mtrw/ReMIXTURE")
#'
#' #new analysis using example data provided
#' my_analysis <- ReMIXTURE$new(
#'   distance_matrix = ReMIXTURE::ReMIXTURE_example_distance_matrix,
#'   region_table = ReMIXTURE::ReMIXTURE_example_region_positions
#' )
#'
#' #A good selection of h-cutoffs should span the lower end of all distance peaks
#' my_analysis$plot_distance_densities(set_xlims = c(0,0.2))
#'
#' #Run the analysis
#' my_analysis$run(iterations = 500,subsample_proportions = c(0.7,0.8,0.9),h_cutoffs=seq(.008,.02,l=10))
#'
#' #Raw results
#' my_analysis$run_results
#'
#' #Tools to assess h-cutoff values
#'    #Good values are often found between where the median clustercounts are equal, and where the inter-region cluster counts maximise
#' my_analysis$plot_clustercount_diag_nondiag_means()
#'    #A good value has a "balanced" heatmap, the main features in which are stable at nearby values
#' my_analysis$plot_heatmaps()
#'    #A good value does not incur many or any per-region cluster counts near one. A few may be inevitable in some datasets though.
#' my_analysis$plot_clustercounts()
#'
#' #Plot the regions on the globe, and add raw plotting data to my_analysis$run_results
#' my_analysis$plot_maps(run = 9,alpha_norm_per_region = T)
#'
#' @docType _PACKAGE
#' @name ReMIXTURE
NULL
