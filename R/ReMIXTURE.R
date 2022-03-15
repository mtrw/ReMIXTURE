#' ReMIXTURE: Regionwise similarity analysis using a resampling/heirarchical clustering based method.
#'
#'
#' @section Warning:
#' This is under development, in early stages.
#' This version is significantly updated and differs non-trivially from the version presented in the paper Tripodi & Rabanus-Wallace et. al. 2021 about peppers. If you wish to use this method, please contact me to talk about how to cite it. A publication is coming soon.
#'
#' @examples
#'
#' devtools::install_github("https://github.com/mtrw/ReMIXTURE")
#'
#' my_analysis <- ReMIXTURE$new(
#'   distance_matrix = ReMIXTURE::ReMIXTURE_example_distance_matrix,
#'   region_table = ReMIXTURE::ReMIXTURE_example_region_positions
#' )
#'
#' my_analysis$run(iterations = 300,subsample_proportions = c(0.7,0.8,0.9),h_cutoffs=seq(.008,.02,l=6))
#' my_analysis$run_results
#'
#' my_analysis$plot_clustercount_diag_nondiag_means()
#' my_analysis$plot_heatmaps()
#' my_analysis$plot_clustercounts()
#'
#' my_analysis$plot_maps(run = 9,alpha_norm_per_region = T)
#'
#' @docType package
#' @name ReMIXTURE
NULL
