# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' sayname
#'
#' Says the name of the intended future of this package
#'
#' @section Warning:
#' Does nothing interesting.
#'
#' @param x Anything printable.
#' @return Nothing. Just prints nothing interesting.
#' @examples
#' sayname("Humboldt")
#' @export
sayname <- function(x="nothing") {
  print(paste0("This is ReMIXTURE: ",x))
  print(head(melt(iris)))
}

# melt <- function(x) {
#   cat("abc",x)
# } #as a distraction to check namespacing, will cause sayname to fail


#rmdata <- 1:10
#usethis::use_data(rmdata,overwrite=TRUE)
#.rmdata2 <- 10:1
#usethis::use_data(.rmdata2, internal = TRUE,overwrite=TRUE)
