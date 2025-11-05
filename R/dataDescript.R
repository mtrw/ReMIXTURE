#' Example region table
#'
#' The region table contains metadata, used for map plotting.
#'
#' @format A data table object with 11 rows and 3 variables:
#' \describe{
#'   \item{region}{A region, as stated in colnames or rownames in a distance matrix}
#'   \item{y}{y coordinates for a region center, i.e., standard latitude values in decimal format as you could read off Google maps}
#'   \item{x}{x coordinates for a region center, i.e., standard longitude values in decimal format as you could read off Google maps}
#' }
#' @source { Rabanus-Wallace & Tripodi et al. (2021) }
"ReMIXTURE_example_region_table"

#' Matrix of distances between different samples
#'
#' Dataset contains 972  chili hot pepper samples from each region
#' from Rabanus-Wallace & Tripodi et al. (2021) data.
#' This is a subset of the original data.
#'
#' @format A matrix with 972 rows and 972 columns. Matrix is diagonal.
#' Upper triangle and lower triangle are identical.
#' It is possible to populate one of the triangles with zeroes
#'
#' @source { Rabanus-Wallace & Tripodi et al. (2021) }
"ReMIXTURE_example_distance_matrix"

#' Map data, 10 metre resolution
#'
#' This is used to draw the maps, and can be passed to `$run()` via the `mapData` argument. In principle, a custom data.table with the same columns could be substituted to get fully customisable maps. The data come from NaturalEarth.
#'
#' @format A data.table giving the polygon outlines of map objects (countries, typically).
#' The columns are as follows:
#
#' lon (numeric) The longitude of a point
#' lat (numeric) The latitude of a point
#' sovereignty / name (character) Identifying information for each polygon from NaturalEarth. Not used by ReMIXTURE.
#' plotGrp (integer) Used to uniquely identify each polygon
#' ptIdx (integer) The order in which the points should be drawn/joined to draw the polygon correctly.
#'
#' @source { https://www.naturalearthdata.com/ }
"mapData10"

#' Map data, 110 metre resolution
#'
#' This is used to draw the maps, and is the default map data passed to `$run()` via the `mapData` argument. In principle, a custom data.table with the same columns could be substituted to get fully customisable maps. The data come from NaturalEarth.
#'
#' @format A data.table giving the polygon outlines of map objects (countries, typically).
#' The columns are as follows:
#' lon (numeric) The longitude of a point
#' lat (numeric) The latitude of a point
#' sovereignty / name (character) Identifying information for each polygon from NaturalEarth. Not used by ReMIXTURE.
#' plotGrp (integer) Used to uniquely identify each polygon
#' ptIdx (integer) The order in which the points should be drawn/joined to draw the polygon correctly.
#'
#' @source { https://www.naturalearthdata.com/ }
"mapData110"


