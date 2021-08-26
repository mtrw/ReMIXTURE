#' Metadata for map plotting
#'
#' A dataset contains metadata, used for map ploting.
#'
#' @format A data table object with 11 rows and 4 variables:
#' \describe{
#'   \item{region}{A region, as stated in colnames or rownames in a distance matrix}
#'   \item{y}{y coordinates for a region center, i.e., standard latitude values in decimal format as you could read off Google maps}
#'   \item{x}{x coordinates for a region center, i.e., standard longitude values in decimal format as you could read off Google maps}
#'   \item{col}{Color to use in hexadecimal format}
#' }
#' @source {Rabanus-Wallace & Tripodi et al. (2021) }
"geocenters"
