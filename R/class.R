
# The fundamental class. Here, just data members.

#' @export
ReMIXTURE <- R6::R6Class(


  private = list(


    ################################ INITIAL INPUTS ##################################
    m = matrix(), # a distance matrix with rownames and colnames giving regions
    rt = data.table(), # region table, an info table with columns "region", "lat" , "long" , and optionally "colour"






    ################################ PARAMETERS ######################################

    #### ALGORITHM #####
    nits = NULL, # a record of the number of iterations used for the analysis
    hcut = NULL, # hcut used for the clustering step
    pr_samp = NULL, # proportion of dataset subsampled each iteration

    #### PLOTTING ######




    ################################ RESULTS #########################################
    diversity = NULL, # average number of clusters per region
    var_diversity = NULL, # variance of number of clusters per region
    overlap = NULL, # average number of region-overlapping / region-unique (diagonal) clusters per region
    var_overlap = NULL, # variance of number of region-overlapping / region-unique (diagonal) clusters per region

    plot_circles_diversity = NULL,
    plot_circles_uniqueness = NULL,
    plot_lines_overlap = NULL,


    runflag = FALSE, #flags a run has been done, results computed and saved
    plotflag = FALSE #flags a plot has been done, results computed and saved
  ),


  public = list(

    #' @description
    #' Create a new ReMIXTURE object.
    #' @param distance_matrix An all-vs-all, full numeric distance matrix, with rownames and
    #'      colnames giving the region of origin of the corresponding individual.
    #' @param region_table A data.table describing the lat(y)/long(x)s of each region (numeric or integer values), with columns named "region", "x" and "y".
    #' @return a new `ReMIXTURE` object.
    initialize = function(distance_matrix,region_table){

      ################################################################################################
      ce("------------------------------------------------")
      ce("Initialising ReMixture object ...")
      ce("------------------------------------------------\n")


      ################################################################################################
      ce("\tValidating input distance matrix ...")
      private$validate_m(distance_matrix)

      ################################################################################################
      ce("\tValidating input region table ...")
      private$validate_rt(region_table)
      ################################################################################################
      ce("\tChecking distance matrix and region table compatibility ...")
      private$validate_dm_rt(distance_matrix,region_table)

      ################################################################################################
      ce("\tModifying distance matrix if necessary ... ")

      if (all(is.na(diag(distance_matrix)))){
        ce("Diagonals of distance matrix all `NA` --- these will be replaced by zeroes.")
        diag(distance_matrix) <- 0
      } else if (all(is.infinite(diag(distance_matrix)))){
        ce("Diagonals of distance matrix are all +/- Inf --- these will be replaced by zeroes. Be sure to assure the matrix has larger values for more distant pairs.")
        diag(distance_matrix) <- 0
      }
      if(any(distance_matrix>1)){
        ce("Distance matrix has values > 1, and will now have all entries linearly scaled to fit the range [0,1]. If this is an issue, please provide a pre-scaled distance matrix.")
        distance_matrix <- distance_matrix %>% scale_between(0,1)
        stopifnot(all(distance_matrix[diag(distance_matrix)]==0))
      }

      ################################################################################################
      ce("\tSaving input to object ...")
      private$m <- distance_matrix
      diag(private$m) <- Inf
      private$rt <- region_table


      ################################################################################################
      ce("\n------------------------------------------------")
      ce("Initialisation complete. Have a very. safe. day.")
      ce("------------------------------------------------")
    }

  )
)

