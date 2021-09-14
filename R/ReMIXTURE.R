#' ReMIXTURE
#'
#' Regionwise similarity analysis using a resampled nearest-neighbour method.
#'
#' @section Warning:
#' Under development.
#'
#' @return A ReMIXTURE class object.
#' @examples
#' my_analysis <- ReMIXTURE$new(
#' distance_matrix = ReMIXTURE::distance_matrix,
#'    info_table = ReMIXTURE::geocenters
#' )
#'
#'#Takes about 5--10 minutes without increasing the parallelise argument (avilable on linux-based systems).
#' my_analysis$run(iterations = 500,parallelize = 1)
#'
#' my_analysis$plot_heatmap()
#'
#' m <- my_analysis$plot_maps()
#' dev.off()
#' pdf("MY_PLOTS.pdf",width=5,height=50)
#' print(m)
#' dev.off()
#' @export
ReMIXTURE <- R6::R6Class(
  ################# Public ################
  public = list(
    #' @description
    #' Create a new ReMIXTURE object.
    #' @param distance_matrix An all-vs-all, full numeric distance matrix, with rownames and
    #'      colnames giving the region of origin of the corresponding individual.
    #' @param info_table A data.table rescribing the lat(y)/long(x)s of each region, with columns named "region", "x", "y", and optionally "col", to give a HEX colour to each region.
    #' @return a new `ReMIXTURE`` object.
    initialize = function(distance_matrix,info_table=NULL){ #constructor, overrides self$new
      #browser()

      #call validators for dm and it if they exist
      #browser()
      if( ncol(distance_matrix) != nrow(distance_matrix)){
        stop("Distance matrix must be square")
      }
      if( #lower triangular dm --- fill
        all(distance_matrix[lower.tri(distance_matrix,diag=F)]==0) & !all(distance_matrix[upper.tri(distance_matrix,diag=F)]==0)
      ){
        warning("Detected a probable triangular distance matrix as input. Zero entries in lower triangle will be filled based on the upper triangle")
        distance_matrix <- fill_lower_from_upper(distance_matrix)
      }

      if( #upper triangular dm --- fill
        !all(distance_matrix[lower.tri(distance_matrix,diag=F)]==0) & all(distance_matrix[upper.tri(distance_matrix,diag=F)]==0)
      ){
        warning("Detected a probable triangular distance matrix as input. Zero entries in upper triangle will be filled based on the lower triangle")
        distance_matrix <- fill_upper_from_lower(distance_matrix)
      }




      #call validators for dm and it if they exist
      private$validate_dm(distance_matrix)
      if( !is.null(info_table) ){
        private$validate_it(info_table)
        private$it <- info_table
        #if colour not present, auto-fill
        if( is.null(info_table$col) ){ # No colours provided --- assign!
          warning("No colour column in info_table provided. Colour will be manually added.")
          info_table[ , col := replace_levels_with_colours(region) ]
        }
      } else {
        warning("No info table provided. Must be inputted manually with $info_table() before $run() can be called.")
      }

      private$dm <- distance_matrix

    },


    #' @description
    #' Run the ReMIXTURE analysis. Requires the information table to have been provided upon initialisation or later with $info_table().
    #' @param iterations The number of samplings requested.
    #' @param resample A length two numeric vector in the format c(`iterations`,`proportion`). If `iterations`>100 and 0.1<=`proportion`<=0.9, then $run_resampling() will be automatically run with `iterations` and `proportions` as arguments.
    #' @param parallelize If TRUE, then all available cores will be used for an analysis
    #' @return Nothing
    run = function(iterations=1000, resample=c(0,0), parallelize = 1){
      #run the method to fill private$counts (define this somewhere else for clarity and call it here)
      # if resample==T, then run the resampling stuff too
      dm <- private$dm
      gpcol <- colnames(dm)
      gplist <- data.table::data.table(region=colnames(dm))[,.N,by=.(region)]

      #index the positions of each region group
      gplist$offset <- c(0,rle(gpcol)$lengths) %>% `[`(-length(.)) %>% cumsum %>% `+`(1)
      gplist[,idx:=1:.N]

      sampsize <- (min(table(gpcol)) * (2/3)) %>% round #SET: how many samples per iteration (from each region)

      #set up some vectors to store info later
      blocksize <- sampsize * nrow(gplist)
      outsize <- iterations * blocksize
      raw_out <- data.table::data.table( #to store raw output each iteration
        p1 = character(length=outsize),
        p2 = character(length=outsize),
        dist = numeric(length=outsize)
      )

      #run the iterations
      l_ply( seq( 1 , outsize-blocksize+1 , by=blocksize ), function(top_row){ #one mclapply call here. should give the top row of raw_out to write
        #initialise and fill the `select` vector
        select <- vector(mode="integer",length=blocksize) #to store a list of the randomly selected samples each iteration
          #dev iteration = 1
        gplist[,{
          select[(sampsize*(idx-1)+1):((sampsize*(idx-1))+sampsize)] <<- sample(N,sampsize)-1+offset
        },by="idx"] %>% invisible

        #Find closest neighbours for the selected sample, store results in output table
        mclapply( mc.cores=parallelize , 1:blocksize,function(i){# dm[select,select],1,function(r){ #mclapply call here, must have prior knowledge in each round of the rownums of raw_out to which it will write (which will be the iterated-overt thing, and based on the top_row given in the outer loop)
          raw_out[top_row+i-1]$p1 <<- colnames(dm)[select[i]]
          raw_out$p2[top_row+i-1] <<- colnames(dm)[select][which(dm[select[i],select]==min(dm[select[i],select]))[1]]
          raw_out$dist[top_row+i-1] <<- min(dm[select[i],select])
        }) %>% invisible
      }) %>% invisible



      private$raw_out <- raw_out

      #summarise the output
      private$counts <- private$raw_out[ , .(count=.N) , by=.(p1,p2) ][ is.na(count) , count:=0 ]

      data.table::setorder(private$raw_out,p1,p2,-dist)
      private$raw_out[,idx:=1:.N,by=.(p1)]

      if (resample[1]>100 & 0.1<=resample[2] & resample[2]<=0.9){
        self$run_resampling(iterations=resample,qualifier=resample)
      }

      invisible(self)
    },
    #' @description
    #' Run the resampling of ReMIXTURE analysis results. Requires the information table to have been provided upon initialisation or later with $info_table().
    #' @param iterations The number of samplings requested.
    #' @param qualifier Float number between 0 and 1, which states how many samples to retrieve from every iteration.
    #' @return A sense of profound satisfaction.
    run_resampling = function(iterations=1000, qualifier = 0.9){

      samplesize <- nu(private$raw_out$iteration)*qualifier #SET: How many items to sample each time
      nrowsit <- (nu(private$raw_out$p1)**2)
      nrowsout <- nrowsit*iterations
      #to store output
      itcount <- data.table::data.table(
        p1=character(length=nrowsout),
        p2=character(length=nrowsout),
        count=numeric(length=nrowsout),
        resamp=numeric(length=nrowsout)
      )

      #perform resampling
      for(it in 1:iterations){
        #it <- 1
        ce("It: ",it)
        selectit <- sample(unique(private$raw_out$iteration),samplesize)

        fill <- data.table::setDT(expand.grid(p1=unique(private$raw_out$p1),p2=unique(private$raw_out$p2)))
        insert <- private$raw_out[ iteration %in% selectit , .(count=.N,resamp=it) , by=.(p1,p2) ]
        insert <- insert[fill,on=.(p1,p2)]
        insert[is.na(count),count:=0]
        insert[is.na(resamp),resamp:=it]

        itcount[(nrowsit*(it-1)+1):((nrowsit*(it-1))+nrow(insert))] <- insert
      }

      #summarise output
      itcount[, pct:=(count/sum(count))*100 , by=.(resamp,p1) ]
      itcount <- itcount[, .(sd_pct=sd(pct),mean_pct=mean(pct)) , by=.(p1,p2) ]
      itcount[, description:=paste0( round(mean_pct-(2*sd_pct),digits=2)," (",round(mean_pct,digits=2),") ",round(mean_pct+(2*sd_pct),digits=2)  )  ]
      private$resample <- itcount
    },

    #' @description
    #' Plot the heatmap of RGOs. Requires that $run() has been called.
    #' @return A ggplot2 plot object.
    plot_heatmap = function(){
      #produce plots
      cnormed <- data.table::copy(private$counts)[,prop:=count/sum(count),by=.(p1)]
      cnormed[p1!=p2][order(prop)]
      cm <- as.matrix(data.table::dcast(cnormed,formula=p1~p2,value.var="prop")[,-"p1"])
      rownames(cm) <- colnames(cm)
      hmplot <- pheatmap::pheatmap(cm,cluster_rows = F,cluster_cols = F)
      return(hmplot)
    },



    #' @description
    #' Plot all inter-region RGOs on a map. Requires that $run() has been called.
    #' @param lat_angle A value indicating how the globe should be rotated along the latitude axis, in standard latitude values in decimal format as you could read off Google maps
    #' @param lon_angle A value indicating how the globe should be rotated along the longitude axis, in standard longitude values in decimal format as you could read off Google maps
    plot_maps = function(lat_angle=30,lon_angle=40){
      #dev lat_angle=30;lon_angle=40
      #check
      #produce plots
      cnormed <- data.table::copy(private$counts)[,prop:=count/sum(count),by=.(p1)]
      cnormed<-cnormed[p1!=p2][order(prop)]

      coords<-data.table::copy(private$it)
      coords[,size:=private$counts[p1==region & p2==region]$count,by=region]
      coords[,size:=size %>% scale_between(2,7)]
      coords <- coords[!is.na(size)]

      cnormed[,id:=1:.N]
      cnormed <- coords[,.(p1=region,x1=x,y1=y)][cnormed,on="p1"]
      cnormed <- coords[,.(p2=region,x2=x,y2=y)][cnormed,on="p2"]

      cdat <- plyr::ldply(1:nrow(coords),function(i){
        c <- circle(coords$x[i],coords$y[i],coords$size[i])
        c[,region:=coords$region[i]]
      })
      data.table::setDT(cdat)
      cdat <- coords[,.(region,col)][cdat,on="region"]

      ldat <- cnormed[p1 != p2][, data.table(
        region = rep(p1,2),
        x = c(x1,x2),
        y = c(y1,y2),
        count = as.numeric(rep(prop,2))
      ) ,by="id"]
      ldat <- coords[,.(region,col)][ldat,on="region"]


      require(ggplot2)
      require(ggspatial)
      world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
      p <- ggplot(data = world) +
        geom_sf(lwd=0.05) + theme_bw() + theme(panel.border = element_blank()) +
        coord_sf(crs = paste0("+proj=laea +lat_0=",lat_angle," +lon_0=",lon_angle," +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "))  + xlab(NULL) + ylab(NULL)

      #p + geom_spatial_polygon(data=cdat,aes(x=x,y=y,fill=region,group=region),alpha=0.8)

      P <- p +
        geom_spatial_path(data=ldat,aes(x=y,y=x,alpha=count,size=count,colour=col),lineend="round") +
        geom_spatial_polygon(data=cdat,aes(x=x,y=y,group=region,fill=col),colour="black") +
        theme(legend.position = "none") +
        facet_grid(region~.)
      return(P)
    }
  ),




  ################# Private ################
  private = list(
    dm = matrix(), # a distance matrix with rownames and colnames giving regions
    it = data.table::data.table(), # an info table with columns "region", "lat" , "long" , and optionally "colour"
    iterations = NA_integer_, # a record of the number of iterations used for the
    validate_dm = function(in_dm){
      #check matrix is a legit distance matrix
      if( !is.matrix(in_dm) ){
        stop( paste0("Argument to distance_matrix must be a matrix (class(in_dm)==",class((in_dm)),")") )
      }
      if( ncol(in_dm) != nrow(in_dm) ){
        stop( paste0("Argument to distance_matrix must be a square matrix") )
      }

      #check if there is NAs or Inf on a diag. Convert to zeroes
      if (all(is.na(diag(in_dm)))){
        diag(in_dm) <- 0
      } else if (all(is.infinite(diag(in_dm)))){
        diag(in_dm) <- 0
      }
      #check zeroes on diagonal
      if( !all(in_dm[diag(in_dm)]==0) ){
        stop("Self-distance (i.e. distance matrix diagonals) should always be zero")
      }

      #check rows and columns are the same in_dm <- N
      if ( !sapply(1:nrow(in_dm),function(r) { all(in_dm[r,]==in_dm[,r]) }) %>% all ){
        stop("Distance matrix is not diagonal")
      }

      #check groups have decent numbers
      #check rowsnames/colnames exist and rownames==colnames
      if (is.null(colnames(in_dm)) | is.null(rownames(in_dm))){
        stop( "Column and row names of input matrix must provide region information" )
      }
      if( !all(colnames(in_dm) == colnames(in_dm)) ) {
        stop( "Column and row names of input matrix must be the same" )
      }

    },
    validate_it = function(in_it){
      #check all columns "region", "x"(longitude) , "y"(latitude) present and character/numeric/numeric
      if( !data.table::is.data.table(in_it) ){
        stop("Info table must be a data.table")
      }
      if( any(!c("region","x","y") %in% colnames(in_it) ) ){
        stop("Info table must have all( c(\"regions\",\"x\",\"y\") %in% colnames(.) )")
      }
      if( !all(unique(colnames(private$dm)) %in% in_it$region) ){
        stop("All regions present in distance matrix must have entries in the info table.")
      }
    },
    raw_out = data.table(), #raw output from sampling
    counts = data.table(), #(normalised, prefereably) count data from sampling
    resample = data.table() #(normalised, prefereably) count data from sampling
  ),

  ################# Active ################
  active = list( #functions that look like vars. mostly for getters and setters of privates, since they can perform checks
    #' @field distance_matrix Return the distance matrix as a data.table
    distance_matrix = function(in_dm){
      if(missing(in_dm)){
        private$dm
      } else {
        warning("Distance matrix cannot be set after initialisation")
      }
    },
    #' @field results Return the numerical results of $run() as a data.table
    results = function(res){
      if (missing(res)){
        private$counts
      } else {
        warning("Results matrix cannot be modified manually")
      }
    },
    #' @field resample_matrix Return the results of resampling conducted by $run_resampling() as a data.table
    resample_matrix = function(resample){
      if (missing(resample)){
        private$resample
      } else {
        warning("Resample matrix cannot be modified manually")
      }
    },
    #' @field info_table Return the information table that assigns latitudes, longitudes, and colours to the regions as a data.table
    info_table = function(in_it){
      if(missing(in_it)){
        return(private$it)
      } else { #validate and replace private$it
        private$validate_it(in_it)
        #if colour not present, auto-fill
        if( is.null(in_it$col) ){ # No colours provided --- assign!
          warning("No colour column in info_table provided. Colour will be manually added.")
          in_it[ , col := replace_levels_with_colours(region) ]
        }
        private$it <- in_it
      }
    }
  )
)














































































#non-parallel $run

# run = function(iterations=1000, resample=c(0,0), parallelize = F){
#   #run the method to fill private$counts (define this somewhere else for clarity and call it here)
#   # if resample==T, then run the resampling stuff too
#   dm <- private$dm
#   gpcol <- colnames(dm)
#   gplist <- data.table::data.table(region=colnames(dm))[,.N,by=.(region)]
#
#   #index the positions of each region group
#   gplist$offset <- c(0,rle(gpcol)$lengths) %>% `[`(-length(.)) %>% cumsum %>% `+`(1)
#   gplist[,idx:=1:.N]
#
#   sampsize <- (min(table(gpcol)) * (2/3)) %>% round #SET: how many samples per iteration (from each region)
#
#   #set up some vectors to store info later
#   outsize <- iterations * sampsize * nrow(gplist)
#   select <- vector(mode="integer",length=sampsize*nrow(gplist)) #to store a list of the randomly selected samples each iteration
#   raw_out <- data.table::data.table( #to store raw output each iteration
#     p1 = character(length=outsize),
#     p2 = character(length=outsize),
#     dist = numeric(length=outsize),
#     iteration = integer(length=outsize)
#   )
#   insert <- 1 #a flag
#
#   #run the iterations
#   if (parallelize){
#     ce("Parallelisation is not yet implemented. Sorry.")
#     # options(parallelly.makeNodePSOCK.setup_strategy = "sequential")
#     # cat("Setting up parallel architecture \n")
#     # future::plan(future::multisession)
#   } else {
#     # future::plan(future::sequential)
#   }
#   for(iteration in 1:iterations){
#     #fill the `select` vector
#     #dev iteration = 1
#     gplist[,{
#       select[(sampsize*(idx-1)+1):((sampsize*(idx-1))+sampsize)] <<- sample(N,sampsize)-1+offset
#     },by="idx"] %>% invisible
#
#     #Find closest neighbours for the selected sample, store results in output table
#     rnum <- 1
#     #r = dm[select,select][1,]
#     # future.apply::future_apply(dm[select,select],1,function(r){
#     apply(dm[select,select],1,function(r){
#       raw_out$p1[insert] <<- colnames(dm)[select][rnum]
#       raw_out$p2[insert] <<- colnames(dm)[select][which(r==min(r))[1]]
#       raw_out$dist[insert] <<- min(r)[1]
#       raw_out$iteration[insert] <<- iteration
#       rnum <<- rnum+1
#       insert <<- insert+1
#     }) %>% invisible
#
#     ce("% complete: ",round((iteration/iterations)*100, 4))
#   }
#   private$raw_out <- raw_out
#   #summarise the output
#   private$counts <- private$raw_out[ , .(count=.N) , by=.(p1,p2) ][ is.na(count) , count:=0 ]
#
#   data.table::setorder(private$raw_out,p1,p2,-dist)
#   private$raw_out[,idx:=1:.N,by=.(p1)]
#
#   if (resample[1]>100 & 0.1<=resample[2] & resample[2]<=0.9){
#     self$run_resampling(iterations=resample,qualifier=resample)
#   }
#
#   invisible(self)
# }
