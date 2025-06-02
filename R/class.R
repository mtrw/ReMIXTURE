#' ReMIXTURE
#'
#' @description
#' ReMIXTURE is a tool for visualising genetic (or any) diversity across space.
#'
#' It produces intuitive (and pretty!) plots that show the diversity in different geographic regions, and how it is shared between regions.
#' These plots are based on ReMIXTURE's diversity metric, which is designed to be maximally intuitive. Basically the plots are designed work like Venn diagrams: Diversity in a group is partitioned into diversity that is unique to that group, and that which overlaps other groups--and this is quantified for every focal--target group pair.
#'
#' ReMIXTURE is designed for data where many samples exist from many (or ... at least a few) locations. It requires just a distance matrix describing the difference between every sample pair, and a table describing where on Earth each sample is located (latitude and longitude).
#'
#' This is the first release version of the package, and development is ongoing.
#'
#' For a short lesson on how to run ReMIXTURE, follow the example below.
#'
#' @export
#' @examples
#' # Create ReMIXTURE analysis object from an inbuilt dataset (Tripodi & Rabanus-Wallace et al. 2021. PNAS).
#' rm <- ReMIXTURE$new(
#'   distance_matrix = ReMIXTURE_example_distance_matrix,
#'   region_table = ReMIXTURE_example_region_positions
#' )
#'
#' rm$plot_distance_densities(samePlot=T)
#' par(mfrow=c(1,1)) # reset graphics layout
#' rm$plot_MDS()
#'
#' # Default run
#' # A run involves taking a subsample of all the samples in the distance matrix, such that every region has the same number of samples included. This subsample is heirarchically clustered (imagine creatting a dentrogram and cutting it at a certain height). Clusters are then counted. The number of clusters in which a region appears is a proxy for that region's total genetic diversity. The number of clusters in which *only* one region appears is a proxy for the diversity unique to that region. The number of clusters in which members of a pair of clusters appears is a proxy for the diversity that is overlapping between those two regions. the concept is similar to the way a Venn diagram might work. This is done `iterations` times, and all these various counts are averaged at the end.
#'
#' # How is it decided where the tree is cut? This value 'H' can be set several ways. By default FINISH ME!!! (and eventually turn me into a vignette)
#' rm$run()
#' rm$plot_maps(focalRegion = "Central America")
#' rm$plot_maps(focalRegion = "Central America",curvature_matrix = "random") # Add curved joining lines (randomly generated)
#'
#'
#' # Try different fixed cutoffs (a trial run with low iterations)
#' rm$run(
#'   iterations=100,
#'   h_cutoff = seq(from=0.001,to=0.03,length.out=15)
#' )
#' rm$plot_heatmaps()
#' rm$plot_h_optimisation()
#'
#' rm$plot_maps(
#'   run=6, #Smaller H emphasises uniqueness
#'   width_max = 20.0,
#'   alpha_max = 1.0, # Adjust alpha scale. Since this is a low-overlap dataset, it can make it easier to see small differences in the mounts of overlap
#'   focalRegion = "Central America"
#' )
#' rm$plot_maps(
#'   run=8, #Larger H emphasises overlap
#'   width_max = 20.0,
#'   alpha_max = 1.0,
#'   focalRegion = "Central America"
#' )
#'
#' ##### Try with a selection of user-given random cutoffs
#' rm$run(
#'   iterations=500, # For random cutoffs we want a fair number of iterations, get a cup of tea
#'   h_cutoff = seq(from=0.001,to=0.03,length.out=15),
#'   h_cutoff_normal_sigma = rep(0.01,15) # Set the standard deviation for each cutoff tried
#' )
#' rm$plot_heatmaps()
#' rm$plot_h_optimisation()
#'
#' rm$plot_maps(
#'   run=12, # Strikes a good balance I would say
#'   width_max = 20.0,
#'   alpha_max = 1.0,
#'   focalRegion = "Central America"
#' )
#'
#'#### Run using random-empirical cutoffs (not recommended except for very very weirdly structured datasets maybe, and maybe even then ReMIXTURE is not the tool for you)
#' rm$run(
#'   iterations=500, # For random cutoffs we want a fair number of iterations, get a cup of tea
#'   h_cutoff="random-empirical"
#' )
#'
#' rm$plot_maps(
#'   width_max = 20.0,
#'   alpha_max = 1.0,
#'   focalRegion = "Central America"
#' )
ReMIXTURE <- R6::R6Class("ReMIXTURE",

  #### PRIVATE ####
  private = list(

    m = matrix(), # a distance matrix with rownames and colnames giving regions
    rt = data.table(), # region table, an info table with columns "region", "lat" , "long" , and optionally "colour"

    results = NULL,

    runflag = FALSE, #flags a run has been done, results computed and saved
    plotflag = FALSE, #flags a plot has been done, results computed and saved

    #### VALIDATORS ####
    validate_rt = function(in_rt){
      if( !data.table::is.data.table(in_rt) ){
        stop("Region position table must be a data.table")
      }
      if( any(!c("region","lon","lat") %in% colnames(in_rt) ) ){
        stop("Region position table must include columns named \"lon\", \"lat\", and \"region\".")
      }
      #check types of cols
      if(!(class(in_rt$lon) == "numeric" | class(in_rt$lon) == "integer") ){
        stop("Column `lon` must contain numeric or integer values.")
      }
      if(!(class(in_rt$lat) == "numeric" | class(in_rt$lat) == "integer") ){
        stop("Column `lat` must contain numeric or integer values.")
      }
      if(!(class(in_rt$region) == "character") ){
        stop("Column `region` must be a character vector.")
      }
      if( !all(unique(colnames(private$dm)) %in% in_rt$region) ){
        stop("All regions present in distance matrix must have entries in the region position table.")
      }
      if( !all(in_rt$lon %between% c(-180,180)) ){
        stop("All lon (longitude) values must fall between +/- 180.")
      }
      if( !all(in_rt$lat %between% c(-85,85)) ){
        stop("All lat (latitude) values must fall between +/- 85.")
      }
    },

    validate_m = function(in_dm){
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
    },

    validate_dm_rt = function(in_dm,in_rt){
      if(! all(colnames(in_dm) %in% in_rt$region)){
        stop("All regions described in the column / row names of the input distance matrix must correspond to entries in the region information table.")
      }
    }
  ),

  #### ACTIVE BINDINGS ####
  active = list(

    #' @field distance_matrix
    #' Return the distance matrix.
    #'
    #' @examples
    #' rm <- ReMIXTURE$new(
    #'   distance_matrix = ReMIXTURE_example_distance_matrix,
    #'   region_table = ReMIXTURE_example_region_positions
    #' )
    #' rm$distance_matrix[1:10,1:10]
    #' rm$region_table
    distance_matrix=function(){
      out_dm <- copy(private$m)
      diag(out_dm) <- 0.0
      return(out_dm)
    },


    #' @field region_table
    #' Return the region table with lat/lon information.
    #'
    #' @examples
    #' rm <- ReMIXTURE$new(
    #'   distance_matrix = ReMIXTURE_example_distance_matrix,
    #'   region_table = ReMIXTURE_example_region_positions
    #' )
    #' rm$region_table
    #' rm$distance_matrix[1:10,1:10]
    region_table=function(){
      rt <- copy(private$rt)
      return(rt)
    },

    #' @field run_results
    #' Return the results of a run (performed with `<ReMIXTURE_object>$run()`)
    #'
    #' @examples
    #' rm <- ReMIXTURE$new(
    #'   distance_matrix = ReMIXTURE_example_distance_matrix,
    #'   region_table = ReMIXTURE_example_region_positions
    #' )
    #' rm$run() # Default run settings
    #' rm$run_results
    run_results=function(){
      if(private$runflag==FALSE){
        stop("Analysis has not been run. Perform using `$run()`")
      }
      return(copy(private$results))
    }
  ),

  #### PUBLIC ####
  public = list(


    #### INTIALISER ####

    #' @description
    #' Create a new ReMIXTURE object.
    #' @param distance_matrix An all-vs-all, full numeric distance matrix, with rownames and colnames giving the region of origin of the corresponding individual.
    #' @param region_table A data.table describing the longitudes/latitudes of each region, with columns named "region" (character), and "lon" and "lat" (numeric or integer). The "region" column must have names corresponding to all the row/column names of the distance matrix.
    #' @return A new ReMIXTURE object.
    #' @examples
    #' # Using build-in ReMIXTURE test datasets
    #' rm <- ReMIXTURE$new(
    #'   distance_matrix = ReMIXTURE_example_distance_matrix,
    #'   region_table = ReMIXTURE_example_region_positions
    #' )
    #' rm$distance_matrix[1:10,1:10]
    #' rm$region_table
    initialize = function(distance_matrix,region_table){


      ce("------------------------------------------------")
      ce("Initialising ReMixture object ...")
      ce("------------------------------------------------\n")


      ce("\tValidating input distance matrix ...")
      private$validate_m(distance_matrix)

      ce("\tValidating input region table ...")
      private$validate_rt(region_table)

      ce("\tChecking distance matrix and region table compatibility ...")
      private$validate_dm_rt(distance_matrix,region_table)

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

      ce("\tSaving input to object ...")
      private$m <- distance_matrix
      diag(private$m) <- Inf
      private$rt <- region_table

      ce("\n------------------------------------------------")
      ce("Initialisation complete.")
      ce("------------------------------------------------")
    },

    #### RUN ####
    #' @description
    #' Run the ReMIXTURE algorithm and save the results in the object.
    #'
    #' @param iterations [1000] The number of times subsamples are drawn, clustered, and counted.
    #' @param subsample_proportions [0.8] size of subsample to draw each iteration. An equal number of samples will be selected from each region, the number being
    #'   \eqn{\text{subsample\_proportions} \times \text{\# of samples in region with the fewest samples}}
    #' (rounded to the nearest integer).
    #' @param h_cutoff ["auto"] Controls the clustering cutoff \eqn{H}. There are three modes possible:
    #' - (Default) Random \eqn{H}, truncated normal distribution. Each iteration, \eqn{H} will be drawn from a normal distribution, truncated to the values given by `h_cutoff_normal_truncRange`. To use this option, you must provide the means of the distributions you want to try in `h_cutoff`, and the standard deviations (one for each mean) in `h_cutoff_normal_sigma`.
    #'   - If `h_cutoff="auto"` (the default), then a single normal distribution will be used, with the same mean and sd as the (non-self i.e. off-diagonal) distances in the distance matrix.
    #' - Fixed \eqn{H}. To use this option, provide a list of values to try as \eqn{H}. A run will be conducted with each, using the same value at each iteration. If you use this mode, you should clearly indicate it in any figure or publication, e.g. "ReMIXTURE plot using truncated random normal \eqn{H} (mean = ..., sd = ..., range=[...,...])".
    #' - Random \eqn{H}, empirical distribution. This samples values of \eqn{H} from the distances in the distance matrix. This should probably never be used, because while it does help produce meaningful plots in cases of extremely structured populations, the results are just not very intuitive. If you use this mode, you should clearly indicate it in any figure or publication, e.g. "ReMIXTURE plot using empirically distributed \eqn{H}".
    #' @param h_cutoff_normal_sigma [NULL] See description for `h_cutoff`.
    #' @param h_cutoff_normal_truncRange [range(self$distance_matrix)] See description for `h_cutoff`.
    #' @param diagnosticPlotMDSclusters [FALSE] A diagnostic tool. If true, at each iteration, will draw the subsampled samples on an MDS plot, and draw hulls around the clusters.
    #'
    #' @return NULL
    #'
    #' @examples
    #' # See ?ReMIXTURE
    run = function(
      iterations=1000,
      subsample_proportions=c(0.8),
      h_cutoff="auto",
      h_cutoff_normal_sigma=NULL,
      h_cutoff_normal_truncRange=range(self$distance_matrix),
      diagnosticPlotMDSclusters=FALSE
    ){
      ce("------------------------------------------------")
      ce("Running ReMixture analysis ...")
      ce("------------------------------------------------\n")

      ce("\tResetting results fields and flags")
      private$results <- list()
      private$runflag <- FALSE
      private$plotflag <- FALSE


      ce("\tSetting up local variables and containers ...")
      #local params
      nits <- iterations
      ind_info <- data.table(
        gp=colnames(private$m),
        idx=1:ncol(private$m)
      )
      gp_vec <- colnames(private$m)
      gp_list <- sort(ind_info[,unique(gp)])
      gp_info <- data.table(
        gp = sort(gp_list),
        gp_idx = 1:length(gp_list)
      )
      ind_info <- gp_info[ind_info,on="gp"]
      nind <- ncol(private$m)
      ngp <- nrow(gp_info)
      param_test_insert <- 1:ngp
      results_insert <- 1

      h_cutoff_type <- if(h_cutoff[1]=="auto"){
        h_cutoff <- mean(private$m[upper.tri(private$m)])
        h_cutoff_normal_sigma <- sd(private$m[upper.tri(private$m)])
        ce("\tRandom H-cutoffs following an automatically chosen (truncated) normal distribution are requested.")
        ce("\t`h_cutoff`=",h_cutoff)
        ce("\t`h_cutoff_normal_sigma`=",h_cutoff_normal_sigma)
        "random-truncNormal"
      } else if(h_cutoff[1]=="random-empirical"){
        "random-empirical"
      } else if (!is.null(h_cutoff_normal_sigma)){
        ce("\tRandom H-cutoffs following user-supplied (truncated) normal distributions are requested.")
        if(
          !(is.numeric(h_cutoff) | is.integer(h_cutoff)) |
          !(is.numeric(h_cutoff_normal_sigma) | is.integer(h_cutoff_normal_sigma)) |
          !(is.numeric(h_cutoff_normal_truncRange) | is.integer(h_cutoff_normal_truncRange)) |
          (length(h_cutoff)!=length(h_cutoff_normal_sigma)) |
          (length(h_cutoff_normal_truncRange)!=2) |
          !(h_cutoff_normal_truncRange[1]<h_cutoff_normal_truncRange[2]) |
          any(h_cutoff_normal_sigma<=0.0)
        ){
          stop("For random (truncated normal) cutoffs, `h_cutoff` and `h_cutoff_normal_sigma` must be numerical or integers vectors of equal length, and `h_cutoff_normal_truncRange` must be a length two vector specifying a range of numbers > 0.")
        }
        "random-truncNormal"
      } else {
        ce("\tConstant h-cutoffs requested.")
        if( !(is.numeric(h_cutoff) | is.integer(h_cutoff)) ){
          stop("h-cutoffs must be numeric or integers.")
        }
        "constant"
      }


      # In case we are doing diagnostic MDS plots, get the coords
      if( diagnosticPlotMDSclusters==TRUE ){
        mds <- self$plot_MDS(doPlot = FALSE)
      }


      for( pr_samp in subsample_proportions ){
        #dev pr_samp <- 0.8
        nSelectPerGrp <- round(ind_info[,.N,by=.(gp)][,min(N)] * pr_samp)
        if(nSelectPerGrp<2){
          private$results <- list()
          private$runflag <- FALSE
          private$plotflag <- FALSE
          stop(paste0("Some regions do not have enough samples for this value of `subsample_proportions` (",subsample_proportions,"), which results in a subsampling number of ",nSelectPerGrp," samples. Increasing `subsample_proportions` could help but most likely you should think about omitting low-sample regions or combining them into other regions."))
        }
        ce("Each iteration, ",nSelectPerGrp," samples from each region will be used.\n")
        for(i_hcut in seq_along(h_cutoff)){
          #dev hcut="auto"; pr_samp=0.8
          sayHc <- if(h_cutoff_type=="random-empirical"){
            h_cutoff_type
          } else if (h_cutoff_type=="random-truncNormal"){
            paste0("truncNorm( mean=",h_cutoff[i_hcut],", sd=",h_cutoff_normal_sigma[i_hcut]," ,range=[",h_cutoff_normal_truncRange[1],",",h_cutoff_normal_truncRange[2],"] )")
          } else if (h_cutoff_type=="constant"){
            round(h_cutoff[i_hcut],digits = 4)
          } else {
            stop("Something is fundamentally wrong with the universe, email mtrw85@gmail.com and alert Tim.")
          }
          ce( "Begin analysis for h_cutoff==" , sayHc , " and subsample_proportions==" , pr_samp , " ..." )

          #container for parameter testing output
          param_test_out <- expand.grid( #innermost loops first
            gp_idx=gp_info$gp, #
            it=1:nits,
            hcut=h_cutoff,
            pr_samp=subsample_proportions,
            nclust=integer(1),
            run=integer(1)
          ) %>% setDT()

          #local result containers
          nclust_counts <- rep(0L,length(gp_list))
          nclust_counts_2 <- rep(0L,length(gp_list))
          counts_mat <- matrix(0L,nrow=ngp,ncol=ngp)
          colnames(counts_mat) <- rownames(counts_mat) <- gp_list
          counts_mat_accumulator_empty <- counts_mat
          counts_2_mat <- matrix(0.0,nrow=ngp,ncol=ngp)
          colnames(counts_2_mat) <- rownames(counts_2_mat) <- gp_list
          ind_info[,clust:=NA_integer_]

          ce("\tIterating ...")
          for(it in 1:nits){ #Begin iteration loop
            if(it %% 100 == 0) {
              ce("\t\tBegin iteration: ",it)
            }
            #Subsample
            ss_selector <- ind_info[,.(s=sample(idx,nSelectPerGrp)),by=.(gp_idx)]$s

            # set hc
            hc <- if(h_cutoff_type=="random-empirical"){
              sample(private$m[ss_selector,ss_selector],1)
            } else if (h_cutoff_type=="random-truncNormal"){
              rTruncNorm(1,h_cutoff[i_hcut],h_cutoff_normal_sigma[i_hcut],h_cutoff_normal_truncRange)
            } else if (h_cutoff_type=="constant"){
              h_cutoff[i_hcut]
            } else {
              stop("Something is fundamentally wrong with the universe, email mtrw85@gmail.com and alert Tim.")
            }



            #Cluster
            ind_info[ss_selector, clust:=cutree(hclust(as.dist(private$m[ss_selector,ss_selector])),h=hc)]
            # Diagnostic plots if requested
            if(diagnosticPlotMDSclusters==TRUE){
              plot(
                mds$mds[ss_selector,]$axisA,
                mds$mds[ss_selector,]$axisB,
                col=mds$legend[data.table(region=mds$mds[ss_selector,]$region),on=.(region)]$col,
                xlim=range(mds$mds$axisA),
                ylim=range(mds$mds$axisB),
                xlab="Axis 1",
                ylab="Axis 2",
                pch=20,
                cex=0.4,
                main=paste0( "Iteration " , it , "; Subsampling " , pr_samp*100,"%; H-cutoff (",h_cutoff_type,"): ", round(hc,2) )
              )
              l_ply(unique(ind_info[ss_selector,]$clust),function(cl){
                #browser()
                cl_selector <- which(ind_info$clust==cl)
                if(length(cl_selector)==1){
                  points(
                    mds$mds[cl_selector,]$axisA,
                    mds$mds[cl_selector,]$axisB,
                    col="#00000077",
                    pch=20,
                    cex=2
                  )
                } else if(length(cl_selector)==2){
                  lines(
                    mds$mds[cl_selector,]$axisA,
                    mds$mds[cl_selector,]$axisB,
                    col="#00000077",
                    lwd=11
                  )
                } else {
                  hull <- chull(mds$mds[cl_selector,]$axisA, mds$mds[cl_selector,]$axisB)
                  polygon(
                    mds$mds[cl_selector,][hull,]$axisA,
                    mds$mds[cl_selector,][hull,]$axisB,
                    col="#00000033"
                  )
                }
              })
              wait("Press any key to go plot hulls for the next iteration ...")
            }

            counts_mat_accumulator <- counts_mat_accumulator_empty
            nclust_counts <- nclust_counts + (t<-ind_info[ss_selector,.(add=nu(clust)),by=.(gp_idx)][gp_info,on=.(gp_idx)][is.na(add),add:=0L][]$add) # in how many unique clusters does each region occur
            nclust_counts_2 <- nclust_counts_2 + t**2

            param_test_out[param_test_insert]$nclust <- t # keep a record of each run
            param_test_out[param_test_insert]$run <- results_insert
            param_test_insert <- param_test_insert + ngp

            ind_info[ss_selector,{ #over clusters
              ugidx <- unique(gp_idx)
              if(length(ugidx) == 1){
                counts_mat_accumulator[gp_idx,gp_idx] <<- counts_mat_accumulator[gp_idx,gp_idx] + 1 #solo cluster--contributes to genetic 'uniqueness'
              } else {
                apply(combn(ugidx,2),2,function(c) { counts_mat_accumulator[c[1],c[2]] <<- counts_mat_accumulator[c[1],c[2]] + 1 } ) #over permutations of members of multigroup cluster
              }
            },by=.(clust)] %>% invisible

            counts_mat_accumulator <- fold_matrix(counts_mat_accumulator)
            counts_mat <- counts_mat + counts_mat_accumulator
            counts_2_mat <- counts_2_mat + counts_mat_accumulator**2

          } #end iteration loop
          #browser()
          ce("\tSummarising and saving results ...")
          private$results$runs[[results_insert]] <- list(
            subsample_proportion=pr_samp,
            h_cutoff_type=h_cutoff_type,
            h_cutoff=if(h_cutoff_type %in% c("constant","random-truncNormal") ){h_cutoff[i_hcut]}else{NA_real_},
            h_cutoff_normal_sigma=if(h_cutoff_type %in% c("random-truncNormal") ){h_cutoff_normal_sigma[i_hcut]}else{NA_real_},
            iterations=nits,
            overlap = counts_mat/nits,
            var_overlap = (counts_2_mat/nits) - (counts_mat/nits)**2,
            diversity = nclust_counts/nits,
            var_diversity = (nclust_counts_2/nits) - (nclust_counts/nits)**2
          )

          results_insert <- results_insert + 1
        } # end hcut loop
      } # end pr_samp loop

      #Save param test output
      private$results$parameter_selection_clustercounts <- param_test_out

      private$runflag <- TRUE
      ce("\n------------------------------------------------")
      ce("ReMixture analysis complete ...")
      ce("------------------------------------------------")
    },


    #### PLOTTING ####

    #' @description
    #' Plot heatmaps representing the unique and overlapped diversity recorded in a ReMIXTURE run. If multiple runs were done, then the user must press [ENTER] to get the heatmaps for each run.
    #'
    #' Heatmaps are extremely useful for assessing how different values of \eqn{H} affect the output.
    #'
    #' @param colPalette [c("#f2f5ff","#000570")] The colour palette for the heatmap. A vector of two or more colours (Hex-RGB or names are fine), which will be interpolated to make the map.
    #' @param ... Additional arguments passed to `stats::heatmap()`
    #'
    #' @return Nothing
    #'
    #' @examples
    #' # See ?ReMIXTURE
    plot_heatmaps = function(
      colPalette=c("#f2f5ff","#000570"),...
    ){
      if(private$runflag==FALSE){
        stop("Analysis has not been run. Perform using `$run()`")
      }
      for(i in 1:(length(private$results$runs))){
        ssp <- private$results$runs[[i]]$subsample_proportion
        hc <- private$results$runs[[i]]$h_cutoff
        its <- private$results$runs[[i]]$iterations
        if(length(private$results$runs)>1){
          wait( paste0("Please press [ENTER] to produce the next heatmap (Run: ",i,"; Subsample proportion: ",ssp,"; H cutoff: ",hc,"; Iterations: ",its,")") )
        }

        heatmap(
          private$results$runs[[i]]$overlap,
          Rowv=NA,
          Colv=NA,
          col=colPalette,
          cexRow=.8,cexCol=.8,
          main="",
          ...
        )
        title(main=paste0("Run ",i,": Subsample proportion: ",ssp,"; H cutoff: ",round(hc,digits = 2),"; Iterations: ",its),side=2,line=0,adj=0)
      }
    },

    #' @description
    #' Plot the cluster counts (ReMIXTURE's diversity metric) for each region in each run. This is useful for assessing how different values of \eqn{H} affect the output, in particular checking there aren't a lot of regions with very low cluster counts, which could be unreliable.
    #'
    #' @return Nothing
    #'
    #' @examples
    #' # See ?ReMIXTURE
    plot_clustercounts = function(){
      if(private$runflag==FALSE){
        stop("Analysis has not been run. Perform using `$run()`")
      }
      require(ggplot2)
      labeldat <- private$results$parameter_selection_clustercounts[,.(nclust=mean(nclust)),by=.(run,pr_samp,hcut)]
      ggplot(private$results$parameter_selection_clustercounts,aes( x=as.factor(hcut) , y=nclust , fill=as.factor(gp_idx) )) +
        geom_violin(scale = "width") +
        facet_grid(as.factor(pr_samp)~.) +
        geom_hline(aes(yintercept=0)) +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          labs(fill="Region / group" , x="h-cutoff" , y="Number of clusters") +
        ggtitle(paste0("Cluster counts by region (",private$results$runs[[1]]$iterations," iterations run)")) +
        geom_text(data=labeldat,aes(label=paste0("Run: ",run),fill=NULL))
    },

    #' @description
    #' These counts are excellent for determining a good value for \eqn{H}, in cases where the default choice is not desirable. A good value typically gives a nice balance of multi- and single-region clusters, typically between where the two values cross over and (as occurs in most datasets) a point where the number of multi-region clusters hits some maximum.
    #'
    #' @return Nothing
    #'
    #' @examples
    #' # See ?ReMIXTURE
    plot_h_optimisation = function(){
      if(private$runflag==FALSE){
        stop("Analysis has not been run. Perform using `$run()`")
      }
      #if(private$)
      d <- setDT(ldply(private$results$runs,function(r){
        t <- c(get_upper_tri(r$overlap),diag(r$overlap))
        tt <- t[t>0]
        p <- tt/sum(tt)
        data.table(
          subsample_proportion = r$subsample_proportion,
          h_cutoff = r$h_cutoff,
          iterations = r$iterations,
          median_clustercount_diag = median(diag(r$overlap)),
          median_clustercount_nondiag = median(get_upper_tri(r$overlap)), #or lower tri
          total_clustercount_diag = sum(diag(r$overlap)),
          total_clustercount_nondiag = sum(get_upper_tri(r$overlap)), #or lower tri
          entropy = -sum( p * log2(p) )
        )
      }))[,run:=paste0("Run: ",1:.N)][]
      pdat_s <- melt(d,id.vars=c("subsample_proportion","h_cutoff","iterations","run"),measure.vars=c("total_clustercount_diag","total_clustercount_nondiag"))
      pdat_e <- melt(d,id.vars=c("subsample_proportion","h_cutoff","iterations","run"),measure.vars=c("entropy"))

      print(
        ggplot(pdat_s,aes(x=h_cutoff,y=value,colour=variable)) +
          geom_line() +
          geom_point() +
          facet_grid(subsample_proportion~.) +
          theme_classic() +
          geom_hline(aes(yintercept=0)) +
          labs( colour = "Cluster counts:" ) +
          scale_color_manual(labels = c("Single-region clusters", "Multi-region clusters"), values = c("#11888a", "#c94d4d")) +
          ylab("Count") +
          xlab("h-cutoff") +
          geom_text(aes(label=run))
      )
    },

    #' @description
    #' Produces the canonical ReMIXTURE plots.
    #'
    #' @param run [NULL] If multiple runs were done, choose which run to use. This is best assessed using `<ReMIXTURE Object>$plot_heatmaps()` and `<ReMIXTURE Object>$plot_h_optimisation()`.
    #' @param focalRegion [NULL] A character string naming one of the regions as it occurs in the region table, chosen as the focal region (from which lines showing inter-region overlap will eminate). If left NULL then a map will be plotted with each region as the focus in turn--this is best used after setting up a multi-panel plot with some variant of `par(mfrow=c(<number of rows>,<number of columns>))`.
    #' @param range_lon [c(-179.0,179.0)] Limit the map to some range of longitudes. Must be a vector of 2 numbers. For awkward reasons, using 180 or -180 can cause graphical bugs.
    #' @param range_lat [c(-85.0,85.0)] Limit the map to some range of latitudes. Must be a vector of 2 numbers.
    #' @param width_max [10.0] The maximum width of the circles/lines, in units of lat/lons.
    #' @param alpha_max [] Normally the alpha of the connecting lines is scaled so
    #' @param projection []
    #' @param curvature_matrix []
    #'
    #' @return Nothing
    #'
    #' @examples
    #' # See ?ReMIXTURE
    plot_maps = function(
      run=NULL,
      focalRegion=NULL,
      range_lon=c(-179.0,179.0),
      range_lat=c(-85.0,85.0),
      width_max=10.0,
      alpha_max=NULL,
      projection=eckertIV,
      curvature_matrix=NULL
    ){
      if(private$runflag==FALSE){
        stop("Analysis has not been run. Perform using `$run()`")
      }
      if(is.null(run) & length(private$results$runs)>1){
        stop("Please provide a run number to plot from (consider using plot_heatmaps() and plot_clustercounts() to assess which run parameters are appropriate).")
      }
      if(is.null(run) & length(private$results$runs)==1){
        run <- 1
      }

      st <- private$results$runs[[run]]$overlap
      rt <- data.table(
        region=colnames(st),
        totDiv=private$results$runs[[run]]$diversity
      )

      rt[,uniqueDiv:=private$results$runs[[run]]$overlap[r,r],by=.(r=region)] # could just pull out diagonal but this seems safer
      rt <- private$rt[rt,on=.(region)]
      maxDiv <- max(rt$totDiv)
      rt[,wTotDiv   := totDiv/maxDiv    * width_max]
      rt[,wUniqueDiv:= uniqueDiv/maxDiv * width_max]
      wst <- st/maxDiv*width_max


      ct <- if(is.null(curvature_matrix)){
        matrix(0.0,ncol=nrow(rt),nrow=nrow(rt),dimnames=list(rt$region,rt$region))
      } else if (curvature_matrix[1]=="random") {
        cm <- matrix(rnorm(nrow(private$rt)**2,0,0.3),nrow=nrow(private$rt))
        rownames(cm) <- colnames(cm) <- private$rt$region
        cm
      } else {
        curvature_matrix
      }

      at <- copy(st)
      for(j in 1:ncol(at)){
        at[,j] <- st[,j]/rt[j]$totDiv
      }
      plotMiddle <- findCentreLL(range_lon,range_lat)
      trt <- copy(rt) %>% rotateLatLonDtLL(-plotMiddle[1],-plotMiddle[2],splitPlotGrps=F)
      rIdxList <- if(!is.null(focalRegion)){
        which(rt$region %in% focalRegion)
      } else {
        1:nrow(trt)
      }

      for(i in rIdxList){
        #i <- 1
        pe <- plotEmptyMap( range_lon, range_lat, projFun=projection )
        #plotMapItem( makeBorder() , range_lon , range_lat, projFun=projection , plotFun=lines,   col="#00000022" , lwd=.4 )
        plotMapItem( makeMapDataLatLonLines() , range_lon , range_lat, projFun=projection , plotFun=lines,   col="#00000022" , lwd=.4 )
        plotMapItem( mapData110               , range_lon , range_lat, projFun=projection , plotFun=polygon, col="#f7bf25" , lwd=0.2 )
        plotMapBorder(                          range_lon , range_lat, projFun=projection , plotEdges=pe ,lwd=4 )

        for(j in 1:nrow(trt)){
          #j <- 3
          if(i==j){next}
          ldt <- curved_rounded_line(
            x1=trt$lon[i],y1 = trt$lat[i],
            x2=trt$lon[j],y2 = trt$lat[j],
            width =     wst[trt$region[i],trt$region[j]],
            curvature = ct[trt$region[i],trt$region[j]]
          ) %>% mat2dtLL()
          plotMapItem(ldt,projFun=projection,plotFun=polygon,col=alpha("black",at[trt$region[i],trt$region[j]]),border="#000000",lwd=0.15)
          cdt <- circle_seg(trt[j]$lon,trt[j]$lat,radius=trt$wTotDiv[j]/2   ) %>% mat2dtLL()
          udt <- circle_seg(trt[j]$lon,trt[j]$lat,radius=trt$wUniqueDiv[j]/2) %>% mat2dtLL()
          plotMapItem(cdt,projFun=projection,plotFun=polygon,col="#00000055") # 'shadow' effect
          plotMapItem(cdt,projFun=projection,plotFun=polygon,col="#000000FF")
          plotMapItem(udt,projFun=projection,plotFun=polygon,col="#FFFFFF")
        }
        cdt <- circle_seg(trt[i]$lon,trt[i]$lat,radius=trt$wTotDiv[i]/2   ) %>% mat2dtLL()
        udt <- circle_seg(trt[i]$lon,trt[i]$lat,radius=trt$wUniqueDiv[i]/2) %>% mat2dtLL()
        plotMapItem(cdt,projFun=projection,plotFun=polygon,col="#00000055") # 'shadow' effect
        plotMapItem(cdt,projFun=projection,plotFun=polygon,col="#000000FF")
        plotMapItem(udt,projFun=projection,plotFun=polygon,col="#FFFFFF")
        title(main=trt$region[i])
      }
    },
    #' @description
    #' Plot inter-sample distances as a density plot over the distance matrix--basically, a convenient way to judge how far apart samples tend to be, and thus how the clustering might behave when various cutoff values (\eqn{H}) are used.
    #'
    #' Each plot tells the story from the perspective of a focal region. Each region is assigned a colour that remains constant across plots. A region's samples' distances to each other are marked with a thick line allowing you to see what colour is for what region.
    #'
    #' @param set_bw [NULL] The bandwidth parameter 'bw' passed to `density()`. By default it will be set automatically.
    #' @param set_xlims [NULL] Limits on the x-axis. Default is set by `density()`.
    #' @param samePlot [FALSE] If TRUE, it will put all the plots on one page.
    #' @param H [NULL] Visualise your values of \eqn{H} on the plot. If only this argument is used, a single line will indicate the cutoff.
    #' @param truncNorm_sd [NULL] Visualise your values of \eqn{H} on the plot. If this argument is used, a truncated normal distribution will be shown.
    #' @param truncNorm_lims [The range of the ReMIXTURE object's distance matrix] Visualise your values of \eqn{H} on the plot. If this argument is used, it will set the truncation limits of the normal distribution.
    #'
    #' @return Nothing
    #'
    #' @examples
    #' # See ?ReMIXTURE
    plot_distance_densities = function(
      set_bw=0.001,
      set_xlims=range(private$m[upper.tri(private$m)]),
      samePlot=FALSE,
      H="auto",
      truncNorm_sd=NULL,
      truncNorm_lims=range(private$m[upper.tri(private$m)])
    ){
      if(H=="auto"){
        H <- mean(private$m[upper.tri(private$m)])
        truncNorm_sd <- sd(private$m[upper.tri(private$m)])
        truncNorm_lims <- range(private$m[upper.tri(private$m)])
      }

      dt1 <- ldply(unique(colnames(private$m)),function(r){ #dev r = "Africa"
        selr <- rownames(private$m)==r
        ldply(unique(colnames(private$m)),function(c){ #dev r = "Africa"
          selc <- colnames(private$m)==c
          data.table(
            x=(0:5000)/5000,
            y=density(private$m[selr,selc],bw=set_bw,from=0,to=1,n=5001)$y,
            region1=c,
            region2=r
          )
        })
      }) %>% setDT

      colTable <- data.table(
        region2 = unique(colnames(private$m)),
        col     = rgb( t(col2rgb(hsv(seq(0, 1, length.out = nu(colnames(private$m))), 1, 1)) / 255) )
      )
      dt1 <- colTable[dt1,on=.(region2)]
      dt1[,col:=rgb(t(col2rgb(col))/255,alpha=fifelse(region1==region2,1,0.4))]
      dt1[,lwd:=fifelse(region1==region2,3,1)]

      if(samePlot==TRUE){
        warning("Altering mfrow parameters for multi-plot graphics. It will be reset to c(1,1).")
        par( mfrow=c(ceiling(sqrt(nrow(private$rt))),ceiling(sqrt(nrow(private$rt)))) )
      }

      rangeY <- range(dt1$y)
      for(r1 in unique(colnames(private$m))){
        sel1 <- rownames(private$m)==r1
        null_plot(set_xlims,rangeY,main=paste0(r1))
        for(r2 in unique(colnames(private$m))){
          sel2 <- rownames(private$m)==r2
          dt1[ region1==r1 & region2==r2 , lines(x,y,col=col,lwd=lwd) ]
        }
        if(!is.null(H)){
          if(!is.null(truncNorm_sd)){
            polygon(
              x=c(
                seq(set_xlims[1],set_xlims[2],l=1000L),
                set_xlims[2],
                set_xlims[1]
              ),
              y=c(
                ((t<-truncnorm::dtruncnorm(seq(set_xlims[1],set_xlims[2],l=1000L),mean=H,sd=truncNorm_sd,a=truncNorm_lims[1],b=truncNorm_lims[2]))/max(t)) * rangeY[2],
                0,
                0
              ),
              border="#00000000",
              col="#00000022"
            )
          } else {
          }
        }
      }
      if(samePlot==TRUE){
        par( mfrow=c(1,1) )
      }
    },
    #' @description
    #' An MDS plot showing relationships among samples
    #'
    #' An MDS ordination plot (in the same flavour as a PCA--in fact in some circumstances they will be the same depending on what your distances represent), useful as a first-pass assessment of the relationships among your samples--what overlaps what, what is over/underrepresented in the sample set, etc.
    #'
    #' These can take a very long time if the number of samples is high. In these cases, I suggest taking a smaller sample of the whole dataset.
    #'
    #' @param axes [c(1,2)] A vector of the MDS axes that should be outputted.
    #' @param colPalette [c("#DD000088","#DDDD0088","#00DD0088","#0000DD88","#DD00DD88")] A vector collection of colours that will be interpolated and applied to the different
    #' @param showLegend [TRUE] Add a legend?
    #' @param doPlot [TRUE] Produce a plot? If not, will still return the data for the MDS.
    #' @param pch [20] The shape of the points. See `?points`.
    #' @param cex [0.8] The size of the points. See `?points`.
    #' @param ... Other arguments passed to `plot()`.
    #'
    #' @return Silently returns a 2-list containing the (1) results of the MDS in a data.table and (2) information about the legend (another data.table).
    #'
    #' @examples
    #' # See ?ReMixture
    plot_MDS = function(
      axes=c(1L,2L),
      colPalette=c("#DD000088","#DDDD0088","#00DD0088","#0000DD88","#DD00DD88"),
      showLegend=TRUE,
      doPlot=TRUE,
      pch=20,
      cex=0.8,
      ...
    ){
      if( length(axes)!=2L | !(is.numeric(axes) | is.integer(axes)) ){ stop("`axes` must be a numeric or integer vector of length 2") }
      if(!is.null(distanceMatrix)){m<-distanceMatrix}
      colTable <- data.table(
        region = unique(colnames(private$m))
      )[,col:=colorRampPalette(colPalette)(.N)]
      dim <- max(axes)
      mds <- cmdscale( as.dist(private$m) , k=dim )  # Also has function of ignoring the Inf diagonals (done for reasons that were presumably good at the time)

      if(doPlot==TRUE){
        plot(
          mds[,axes[1]],
          mds[,axes[2]],
          pch=pch,
          col=colTable[data.table(region=rownames(mds)),on=.(region)]$col,
          cex=cex,
          xlab=paste0("Axis ",axes[1]),
          ylab=paste0("Axis ",axes[2]),
          ...
        )
        if(showLegend==TRUE){
          legend(min(mds[,axes[1]]),max(mds[,axes[2]]),colTable$region,colTable$col,bg="#FFFFFFAA")
        }
      }

      invisible(return(list(axes=axes,mds=data.table(region=rownames(mds),axisA=mds[,1],axisB=mds[,2]),legend=colTable)))
    }



  )
)

