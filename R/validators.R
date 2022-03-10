
ReMIXTURE$set( "private" , "validate_rt" ,
   function(in_rt){
     if( !data.table::is.data.table(in_rt) ){
       stop("Region position table must be a data.table")
     }
     if( any(!c("region","x","y") %in% colnames(in_rt) ) ){
       stop("Region position table must include columns named \"x\", \"y\", and \"region\".")
     }
     #check types of cols
     if(!(class(in_rt$x) == "numeric" | class(in_rt$x) == "integer") ){
       stop("Column `x` must contain numeric or integer values.")
     }
     if(!(class(in_rt$y) == "numeric" | class(in_rt$y) == "integer") ){
       stop("Column `y` must contain numeric or integer values.")
     }
     if(!(class(in_rt$region) == "character") ){
       stop("Column `region` must be a character vector.")
     }
     if( !all(unique(colnames(private$dm)) %in% in_rt$region) ){
       stop("All regions present in distance matrix must have entries in the region position table.")
     }
     if( !all(in_rt$x %between% c(-180,180)) ){
       stop("All x (longitude) values must fall between +/- 180.")
     }
     if( !all(in_rt$y %between% c(-85,85)) ){
       stop("All y (latitude) values must fall between +/- 85.")
     }
   }
)






ReMIXTURE$set( "private" , "validate_m" ,
  function(in_dm){
    if( !is.matrix(in_dm) ){
      stop( paste0("Argument to distance_matrix must be a matrix.") )
    }
    if( ncol(in_dm) != nrow(in_dm) ){
      stop( paste0("Argument to distance_matrix must be a square matrix.") )
    }
    if( !all(in_dm[diag(in_dm)]==0) ){
      stop("Self-distance (i.e. distance matrix diagonals) should all be zero")
    }
    if ( !all(in_dm[upper.tri(in_dm)]==t(in_dm)[upper.tri(in_dm)]) ){
      stop("Distance matrix must be diagonal.")
    }
    if(!all(in_dm>=0)){
      stop("All distance matrix entries must be positive")
    }
    if (is.null(colnames(in_dm)) | is.null(rownames(in_dm))){
      stop( "Column and row names of input matrix must provide region information." )
    }
    if( !all(colnames(in_dm) == rownames(in_dm)) ) {
      stop( "Column and row names of input matrix must be the same." )
    }
    if( min(table(colnames(in_dm))) < 20 ) {
      warning( paste0("Some regions contain dangerously low numbers of samples (the smallest group has ",min(table(colnames(in_dm)))," members), which could lead to dodgy results. Careful!") )
    }
  }
)







ReMIXTURE$set( "private" , "validate_dm_rt" ,
               function(in_dm,in_rt){
                 if(! all(colnames(in_dm) %in% in_rt$region)){
                   stop("All regions described in the column / row names of the input distance matrix must correspond to entries in the region information table.")
                 }
               }
)


