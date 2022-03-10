
ReMIXTURE$set( "active" , "distance_matrix" ,
  function(){
    out_dm <- copy(private$m)
    diag(out_dm) <- 0.0
    return(out_dm)
  }
)





ReMIXTURE$set( "active" , "region_table" ,
  function(){
    rt <- copy(private$rt)
    return(out_dm)
  }
)







ReMIXTURE$set( "active" , "results" ,
  function(){
    if(private$runflag==FALSE){
      stop("Analysis has not been run. Perform using `$run()`")
    }
    return(
      list(
        parameters = list (
          n_iterations = private$nits,
          subsample_proportion = private$pr_samp,
          h_clustering_cutoff = private$hcut
        ),
        diversity_nclusts = private$diversity,
        var_diversity_nclusts = private$var_diversity,
        overlapping_nclusts = private$overlap,
        var_overlapping_nclusts = private$var_overlap
      )
    )
  }
)
