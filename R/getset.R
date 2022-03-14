
ReMIXTURE$set( "active" , "distance_matrix" ,
  function(){
    out_dm <- copy(private$m)
    diag(out_dm) <- 0.0
    return(out_dm)
  }
)




#' @description
#' getset region table
#' @return table
ReMIXTURE$set( "active" , "region_table" ,
  function(){
    rt <- copy(private$rt)
    return(out_dm)
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
