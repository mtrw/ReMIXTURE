#' ReMIXTURE: Regionwise similarity analysis using a resampling/heirarchical clustering based method.
#'
#'
#' @section Warning:
#' This is under development, in early stages.
#' This version is significantly updated and differs non-trivially from the version presented in the paper Tripodi & Rabanus-Wallace et. al. 2021 about peppers. If you wish to use this method, please contact me to talk about how to cite it. A publication is coming soon.
#'
#' @examples
#'
#' my_analysis <- ReMIXTURE$new(
#' distance_matrix = ReMIXTURE::ReMIXTURE_example_distance_matrix,
#' region_table = ReMIXTURE::ReMIXTURE_example_region_positions
#' )
#'
#' my_analysis$run(iterations = 500,subsample_proportion = 0.8,h_cutoff = 0.013)
#'
#' my_analysis$plot_heatmap()
#'
#' my_analysis$plot_maps(alpha_norm_per_region = T)
#'
#' @docType package
#' @name ReMIXTURE
NULL
