.onLoad <- function(libname, pkgname) {
  # Your code here
  message("Thank you for using ", pkgname,"!")
  message("-------------------------------------------------")
  message("For help on all functions use ?ReMIXTURE.")
  message("For a quick tutorial/workflow to get started, follow the examples at the bottom.")
  message("-------------------------------------------------")
  message("Always remember:")
  message("\t1) Separate ReMIXTURE analyses are not comparable unless the same run settings are used for all datasets, and the distance matrix construction method was also identical.")
  message("\t2) The algorithm has changed a lot recently as the package is prepared for official release--cite the github along with the commit hash.")
  message("\t3) Please include your run settings in figure captions, to help build this as a convention, e.g. \"ReMIXTURE map plots, constructed using truncated random normal H (mean = ..., sd = ..., range=<data range>)\".")
  message("\t4) ReMIXTURE is all about overlapping diversity. If your regions do not have overlapping diversity, ReMIXTURE may not be the tool you need.")
  message("\t5) ReMIXTURE was designed to describe genetic data grouped by region. But the data needn't be genetic and the grouping needn't be geographic. Can you think of other good use cases for ReMIXTURE? Please email Tim!")
  message("-------------------------------------------------")
  message("Find any bugs, or have some suggestions? Contact Tim at tim.rabanuswallace@unimelb.edu.au.")
  message("-------------------------------------------------")
}
