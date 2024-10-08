
ReMIXTURE$set( "active" , "distance_matrix" ,
  function(){
    out_dm <- copy(private$m)
    diag(out_dm) <- 0.0
    return(out_dm)
  }
)



#' @field region_table Get the region table from a ReMIXTURE object
#' @name ReMIXTURE
#' @rdname ReMIXTURE
#' @usage ReMIXTURE_object$region_table
#' @format An active binding in an R6 class.
ReMIXTURE$set( "active" , "region_table" ,
  function(){
    rt <- copy(private$rt)
    return(rt)
  }
)







ReMIXTURE$set( "active" , "run_results" ,
  function(){
    if(private$runflag==FALSE){
      stop("Analysis has not been run. Perform using `$run()`")
    }
    return(copy(private$results))
  }
)
