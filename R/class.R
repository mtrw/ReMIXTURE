#' ReMIXTURE
#'
#' @description
#' The ReMIXTURE package
#'
#' \figure{rm_logo.png}{options: width="50\%" alt="ReMIXTURE logo"}
#'
#' ## *GETTING STARTED?*
#' Your first stop is the vignette, a quick and easy tutorial that will teach you what you need to know about how ReMIXTURE works and how to run it properly.
#'
#' `vignettes("ReMIXTURE")`
#'
#' See Examples for a quick reminder of how a normal workflow might look
#'
#' @examples
#'
#'
#' rm <- ReMIXTURE$new(
#'   distance_matrix = ReMIXTURE_example_distance_matrix,
#'   region_table = ReMIXTURE_example_region_table
#' )
#'
#' rm$plot_MDS()
#'
#' # A collection of 16 test runs with automatically-selected parameters
#' rm$run(iterations = 5)
#'
#' rm$plot_h_optimisation()
#'
#' rm$plot_distance_densities(HdistFromRun=5,samePlot = T) # Overlay H distribution as used in run 5
#'
#' rm$plot_results_grid(runs = 5)
#'
#' # Run 5 seems good ... do a proper run with those parameters
#' (params <- rm$get_run_parameters(run = 5))
#' rm$run(
#'   iterations=500,
#'   H=params$H,
#'   H_truncNorm_sd=params$H_truncNorm_sd
#' )
#'
#' #undebug(rm$plot_maps)
#' rm$plot_maps(
#'   # implicit: run=1 (since only 1 run is available)
#'   focalRegion = "Africa"
#' )
#'
#'  Explore some plot options
#' rm$plot_maps(
#'   focalRegion = "Africa",
#'   width_max = 30, # Adjust point/line widths
#'   curvature = 0.5 # Add curvature for lines
#' )
#'
#' par(mfrow=c(2,2),mar=c(1,1,4,1))
#' rm$plot_maps(
#'   focalRegion = rm$region_table$region[1:4],
#'   width_max = 30,
#'   curvature = {set.seed(42.42);"random"},   # Random curvatures
#'   range_lon = c(-120,140),                  # Crop longitude
#'   range_lat = c(-30,60),                    # Crop latitude
#'   projection = equirectangular,             # Set map projection. Also try `winkelIII` or `eckertIV` (default)
#'   diversityCirclesFocalOnly = TRUE          # Omit diversity circles for all regions except the focal region of the subplot
#' )
#' dev.off() # to reset plot parameters!
#'
#' rm$plot_maps(
#'   allOnOneMap = TRUE, # All overlap lines on one map ... sometimes great, mostly too busy.
#'   width_max = 20,
#'   curvature = c(0.1,0.75),
#'   range_lat = c(-20,90),
#'   projection = winkelIII
#' )
#'
#' # Extract info from the object
#' rm$get_run_parameters() # With multiple runs, a list will be returned
#' rm$distance_matrix[1:5,1:5] %>% round(digits=3) # the distance matrix
#' rm$region_table # the region table
#' rm$mds # the mds data (only if $plot_MDS() has been run)
#'
#' @export
ReMIXTURE <- R6::R6Class("ReMIXTURE",

  #### PRIVATE ####
  private = list(

    m = matrix(), # a distance matrix with rownames and colnames giving regions
    rt = data.table(), # region table, an info table with columns "region", "lat" , "long" , and optionally "colour"

    results = NULL,
    mdsPlot = NULL,

    runFlag = FALSE, #flags a run has been done, results computed and saved

    #### VALIDATORS ####
    validate_rt = function(in_rt){
      if( !(data.table::is.data.table(in_rt) || is.data.frame(in_rt)) ){
        stop("Region position table must be a data.frame or data.table")
      }
      in_rt <- data.table::as.data.table(in_rt)
      if( any(!c("region","lon","lat") %in% colnames(in_rt) ) ){
        stop("Region position table must include columns named \"lon\", \"lat\", and \"region\".")
      }
      #check types of cols
      if(!(is.numeric(in_rt$lon) | is.integer(in_rt$lon))){
        stop("Column `lon` must contain numeric or integer values.")
      }
      if(!(is.numeric(in_rt$lat) | is.integer(in_rt$lat))){
        stop("Column `lat` must contain numeric or integer values.")
      }
      if(!(is.character(in_rt$region))){
        stop("Column `region` must be a character vector.")
      }
      if( !all(in_rt$lon %between% c(-180,180)) ){
        stop("All lon (longitude) values must fall between +/- 180.")
      }
      if( !all(in_rt$lat %between% c(-85,85)) ){
        stop("All lat (latitude) values must fall between +/- 85.")
      }
      if( nrow(in_rt)<2 ){
        stop("At least two regions are required ...")
      }
      in_rt
    },

    validate_m = function(in_dm){
      if( !is.matrix(in_dm) ){
        stop( paste0("Argument to distance_matrix must be a matrix.") )
      }
      if( ncol(in_dm) != nrow(in_dm) ){
        stop( paste0("Argument to distance_matrix must be a square matrix.") )
      }
      if(any(!is.finite(in_dm))){
        stop("All distance matrix entries must be finite")
      }
      if( !all(diag(in_dm)==0) ){
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
    },

    validate_variant_sample_table = function(in_sample_table){
      if(is.null(in_sample_table)){
        stop("`vcf_sample_table` must be provided when using `vcf_file`.")
      }
      if(!data.table::is.data.table(in_sample_table) && !is.data.frame(in_sample_table)){
        stop("`vcf_sample_table` must be a data.frame or data.table.")
      }

      in_sample_table <- data.table::as.data.table(in_sample_table)

      if(!all(c("sample", "region") %in% colnames(in_sample_table))){
        stop("`vcf_sample_table` must contain columns named `sample` and `region`.")
      }
      if(!is.character(in_sample_table$sample)){
        stop("Column `sample` in `vcf_sample_table` must be character.")
      }
      if(!is.character(in_sample_table$region)){
        stop("Column `region` in `vcf_sample_table` must be character.")
      }
      if(any(is.na(in_sample_table$sample)) || any(in_sample_table$sample == "")){
        stop("Column `sample` in `vcf_sample_table` must not contain missing or empty values.")
      }
      if(any(is.na(in_sample_table$region)) || any(in_sample_table$region == "")){
        stop("Column `region` in `vcf_sample_table` must not contain missing or empty values.")
      }
      if(anyDuplicated(in_sample_table$sample)){
        stop("`vcf_sample_table$sample` must not contain duplicates.")
      }

      in_sample_table
    },

    filter_invalid_distance_samples = function(dm){
      invalid_message <- paste0(
        "Fewer than 2 samples remain after filtering invalid distances. ",
        "This often means the VCF has too much missing data. ",
        "Consider building the distance matrix yourself, performing missing-data handling or imputation first, or using data with fewer missing values."
      )

      if(nrow(dm) < 2){
        stop(invalid_message)
      }

      #browser()
      invalid <- !isBehaved(dm)
      drop_all_invalid <- apply(invalid, 1, all)
      if(any(drop_all_invalid)){
        dm <- dm[!drop_all_invalid, !drop_all_invalid, drop = FALSE]
      }
      if(nrow(dm) < 2){
        stop(invalid_message)
      }
      invalid <- !isBehaved(dm)
      #diag(invalid) <- FALSE
      drop_any_invalid <- apply(invalid, 1, any)
      drop_all_invalid[1:10]
      if(any(drop_any_invalid)){
        dm <- dm[!drop_any_invalid, !drop_any_invalid, drop = FALSE]
      }

      if(nrow(dm) < 2){
        stop(invalid_message)
      }

      dm
    },

    distance_matrix_from_vcf_inputs = function(vcf_file, vcf_sample_table){
      elapsed_string <- function(start_time){
        sprintf("%.2f sec", proc.time()[["elapsed"]] - start_time)
      }

      total_start <- proc.time()[["elapsed"]]
      snprelate_verbose <- interactive()

      if(!requireNamespace("SNPRelate", quietly = TRUE)){
        stop("Reading VCF input requires the 'SNPRelate' package. Install it from Bioconductor and try again.")
      }
      if(!is.character(vcf_file) || length(vcf_file) != 1L || is.na(vcf_file) || vcf_file == ""){
        stop("`vcf_file` must be a length-1 character path.")
      }
      if(!file.exists(vcf_file)){
        stop("`vcf_file` does not exist: ", vcf_file)
      }
      if(grepl("\\.bcf$", vcf_file, ignore.case = TRUE)){
        stop("Only .vcf and .vcf.gz input files are supported. Please convert the BCF file to VCF or VCF.GZ before using ReMIXTURE.")
      }
      if(!grepl("\\.vcf(\\.gz)?$", vcf_file, ignore.case = TRUE)){
        stop("Only .vcf and .vcf.gz input files are supported.")
      }

      stage_start <- proc.time()[["elapsed"]]
      ce("\tValidating sample table and ...")
      st <- private$validate_variant_sample_table(vcf_sample_table)
      ce("\tValidating sample table ... done (", elapsed_string(stage_start), ")")

      gds_tmp_fname <- tempfile(fileext = ".gds")
      gds_handle <- NULL

      on.exit({
        if(!is.null(gds_handle)){
          tryCatch(
            SNPRelate::snpgdsClose(gds_handle),
            error = function(e) NULL
          )
        }
        if(file.exists(gds_tmp_fname)){
          unlink(gds_tmp_fname)
        }
      }, add = TRUE)

      stage_start <- proc.time()[["elapsed"]]
      ce("\tConverting VCF to temporary GDS ...")
      tryCatch(
        SNPRelate::snpgdsVCF2GDS(
          vcf_file,
          out.fn = gds_tmp_fname,
          method = "biallelic.only",
          verbose = snprelate_verbose
        ),
        error = function(e){
          stop("VCF could not be converted to GDS: ", conditionMessage(e))
        }
      )
      ce("\tConverting VCF to temporary GDS ... done (", elapsed_string(stage_start), ")")

      gds_handle <- tryCatch(
        SNPRelate::snpgdsOpen(gds_tmp_fname),
        error = function(e){
          stop("VCF was converted to GDS, but the temporary GDS file could not be opened: ", conditionMessage(e))
        }
      )

      stage_start <- proc.time()[["elapsed"]]
      ce("\tComputing IBS matrix ...")
      ibs <- tryCatch(
        SNPRelate::snpgdsIBS(
          gds_handle,
          sample.id = st$sample,
          missing.rate = NaN,
          remove.monosnp = FALSE,
          num.thread = NA_real_,
          verbose = snprelate_verbose
        ),
        error = function(e){
          stop("IBS calculation failed for the provided VCF: ", conditionMessage(e))
        }
      )
      ce("\tComputing IBS matrix ... done (", elapsed_string(stage_start), ")")

      stage_start <- proc.time()[["elapsed"]]
      ce("\tConverting IBS to distance matrix ...")
      sample_list <- ibs$sample.id
      missing_samples <- setdiff(st$sample, sample_list)
      if(length(missing_samples) > 0){
        stop(
          "`vcf_sample_table` contains sample IDs that could not be mapped to VCF samples. Examples: ",
          paste(utils::head(missing_samples, 5), collapse = ", ")
        )
      }

      region_labels <- st$region[match(sample_list, st$sample)]

      if(anyNA(region_labels)){
        stop(
          "Internal mapping error: VCF samples were returned by SNPRelate, but some could not be mapped back to `vcf_sample_table$region`."
        )
      }

      dm <- 1 - ibs$ibs
      colnames(dm) <- rownames(dm) <- region_labels
      ce("\tConverting IBS to distance matrix ... done (", elapsed_string(stage_start), ")")

      stage_start <- proc.time()[["elapsed"]]
      ce("\tFiltering invalid distances ...")
      dm <- private$filter_invalid_distance_samples(dm)
      ce("\tFiltering invalid distances ... done (", elapsed_string(stage_start), ")")

      ce("\tDone. Total elapsed time: ", elapsed_string(total_start))

      dm
    }
  ),

  #### ACTIVE BINDINGS ####
  active = list(

    #' @field distance_matrix
    #' Return a copy of the distance matrix.
    distance_matrix=function(){
      out_dm <- copy(private$m)
      diag(out_dm) <- 0.0
      return(out_dm)
    },


    #' @field region_table
    #' Return a copy of the region table with lat/lon information.
    region_table=function(){
      rt <- copy(private$rt)
      return(rt)
    },

    #' @field run_results
    #' Return the raw results of ReMIXTURE runs (runs must be first performed with `<ReMIXTURE_object>$run()`)
    run_results=function(){
      if(private$runFlag==FALSE){
        stop("Analysis has not been run. Perform using `$run()`")
      }
      return(copy(private$results))
    },

    #' @field mds
    #' Return the mds data (created when `<ReMIXTURE_object>$plot_MDS()` is run)
    mds=function(){
      if(is.null(private$mdsPlot)){
        stop("MDS has not been produced. Make it using e.g. `<ReMIXTURE_object>$plot_MDS()`")
      }
      return(copy(private$mdsPlot))
    }
  ),

  #### PUBLIC ####
  public = list(

#' @description
#' Print information about a ReMIXTURE object. To probe the object further use `str(<ReMIXTURE object>)`
#'
#' @returns NULL
#' @usage
#' print(<ReMIXTURE object>)
#' # or just interactively run
#' <ReMIXTURE object>
#' @export
    print = function(){
      #browser()
      cat("---------------- A ReMIXTURE object ----------------\n\n")
      cat(dim(private$m)[1],"samples in",nrow(private$rt),"regions:\n")
      l_ply(1:nrow(private$rt),function(i){cat("\t",private$rt$region[i]," (n=",private$rt$N[i],")\n",sep="")})
      cat("\n")
      cat(if(private$runFlag==FALSE){
        "No runs performed\n"
      } else {
          cat(length(private$results),"run(s) performed with the following parameters:\n")
          l_ply(1:length(private$results),function(i){
            #i<-1
            cat("\tRun:",i)
            cat("\n\t\t")
            cat("Iterations:",private$results[[i]]$iterations)
            cat("\n\t\t")
            cat("Subsample proportion:",private$results[[i]]$subsample_proportion)
            cat("\n\t\t")
            cat("Cutoff (H) method:",private$results[[i]]$H_type)
            cat("\n\t\t")
            cat("H=",private$results[[i]]$H)
            cat("\n\t\t")
            if(is.null(private$results[[i]]$H_truncNorm_sd)){
              cat("Hsd=<NULL>")
            } else {
              cat("Hsd=",private$results[[i]]$H_truncNorm_sd)
            }
            cat("\n\t\t")
            if(is.null(private$results[[i]]$H_truncNorm_sd)){
              cat("Hrange=<NULL>")
            } else {
              cat("Hrange=","[",private$results[[i]]$H_truncNorm_range[1],"--",private$results[[i]]$H_truncNorm_range[2],"]",sep="")
            }
            cat("\n")
        })
        cat("\n")
        if(!is.null(private$mdsPlot)){
          cat("MDS coordinates calculated\n")
        } else {
          cat("MDS coordinates not calculated\n")
        }
        cat("\n")
        #browser()
      })
      cat("----------------------------------------------------")
    },
    #### INTIALISER ####

    #' @description
    #' Create a new ReMIXTURE object.
    #' @param distance_matrix \[NULL\] An all-vs-all, full numeric distance matrix, with rownames and colnames giving the region of origin of the corresponding individual.
    #' @param region_table \[no default\] A data.table or data.frame describing the longitudes/latitudes of each region, with columns named "region" (character), and "lon" and "lat" (numeric or integer). The "region" column must have names corresponding to all the row/column names of the distance matrix.
    #' @param vcf_file \[NULL\] Path to a `.vcf` or `.vcf.gz` file. Provide this instead of `distance_matrix` to build the matrix via `SNPRelate`.
    #' @param vcf_sample_table \[NULL\] A data.frame or data.table with character columns named `sample` and `region`, mapping VCF sample IDs to ReMIXTURE regions.
    #' @return A new ReMIXTURE object.
    initialize = function(
      distance_matrix = NULL,
      region_table,
      vcf_file = NULL,
      vcf_sample_table = NULL
    ){


      ce("------------------------------------------------")
      ce("Initialising ReMixture object ...")
      ce("------------------------------------------------\n")

      if(!is.null(distance_matrix) && !is.null(vcf_file)){
        stop("Provide either `distance_matrix` or `vcf_file`, not both.")
      }
      if(is.null(distance_matrix) && is.null(vcf_file)){
        stop("Provide either `distance_matrix` or `vcf_file`.")
      }
      if(is.null(distance_matrix)){
        ce("\tConstructing distance matrix from VCF input ...")
        distance_matrix <- private$distance_matrix_from_vcf_inputs(
          vcf_file = vcf_file,
          vcf_sample_table = vcf_sample_table
        )
      }

      ce("\tModifying distance matrix if necessary ... ")
      bad_diagonal <- is.na(diag(distance_matrix)) | is.infinite(diag(distance_matrix))
      if(any(bad_diagonal)){
        ce("Diagonals of distance matrix contain `NA` or +/-Inf --- these will be replaced by zeroes.")
        diag_values <- diag(distance_matrix)
        diag_values[bad_diagonal] <- 0
        diag(distance_matrix) <- diag_values
      }

      ce("\tValidating input distance matrix ...")
      private$validate_m(distance_matrix)

      ce("\tValidating input region table ...")
      region_table <- private$validate_rt(region_table)

      ce("\tChecking distance matrix and region table compatibility ...")
      private$validate_dm_rt(distance_matrix,region_table)

      if(any(distance_matrix>1)){
        ce("Distance matrix has values > 1, and will now have all entries linearly scaled to fit the range [0,1]. If this is an issue, please provide a pre-scaled distance matrix.")
        distance_matrix <- distance_matrix %>% scale_between(0,1)
        stopifnot(all(diag(distance_matrix)==0))
      }

      ce("\tTrimming region table if necessary ... ")
      region_table <- region_table[region %in% colnames(distance_matrix)]

      ce("\tSaving input to object ...")
      private$m <- distance_matrix
      diag(private$m) <- Inf
      private$rt <- region_table

      ce("\tAdding sample counts to internal region table ...")
      if(!is.null(private$rt$N)){ warning("Column named 'N' discovered in region table. This will be overwritten. To preserve it, please rename it and initialise again.") }
      tmp <- as.data.table(table(colnames(private$m))) %>% setnames(c("V1"),c("region"))
      private$rt <- tmp[private$rt,on=.(region)]

      ce("\n------------------------------------------------")
      ce("Initialisation complete.")
      ce("------------------------------------------------")
    },

    #' @description
    #' Return the parameters of a ReMIXTURE run.
    #' @param run \[NULL\] If multiple runs were done, choose which run to use. This is best assessed using `<ReMIXTURE Object>$plot_h_optimisation()`, `<ReMIXTURE Object>$plot_results_grid()`, and `<ReMIXTURE Object>$plot_clustercounts()`.
    #' @return An appropriately-labelled list
    get_run_parameters = function(run=NULL){
      if(is.null(run) & length(private$results)!=1){
        stop("Please provide a valid run number to plot from (consider using `<ReMIXTURE Object>$plot_h_optimisation()`, `<ReMIXTURE Object>$plot_results_grid()`, and `<ReMIXTURE Object>$plot_clustercounts()` to assess which run parameters are appropriate).")
      }
      if(is.null(run) & length(private$results)==1){
        run <- 1
      }
      if(private$runFlag==FALSE){
        stop("Analysis has not been run. Perform using `$run()`")
      }
      return(copy(private$results[[run]][c("subsample_proportion","H_type","H","H_truncNorm_sd","H_truncNorm_range","iterations")]))
    },


    #' @description
    #' Build a pairwise distance matrix directly from a `.vcf` or `.vcf.gz` file using `SNPRelate`.
    #' Direct `.bcf` input is not supported by this helper.
    #'
    #' The VCF is converted to a temporary GDS file with `snpgdsVCF2GDS()`, IBS is computed with `snpgdsIBS()`, and the distance matrix is returned as `1 - ibs`.
    #'
    #' @param vcf_file \[no default\] Path to a `.vcf` or `.vcf.gz` file.
    #' @param vcf_sample_table \[no default\] A data.frame or data.table with character columns named `sample` and `region`.
    #' @return A full square numeric distance matrix with region labels in the row and column names.
    distance_matrix_from_vcf = function(vcf_file, vcf_sample_table){
      private$distance_matrix_from_vcf_inputs(
        vcf_file = vcf_file,
        vcf_sample_table = vcf_sample_table
      )
    },

    #### RUN ####
    #' @description
    #' Run the ReMIXTURE algorithm and save the results in the object. Multiple runs can be requested to aid parameter selection (see description for `H`), which is strongly recommended.
    #'
    #' @param iterations \[1000\] The number of times subsamples are drawn, clustered, and counted.
    #' @param subsample_proportions \[0.8\] Size of subsample to draw each iteration. An equal number of samples will be selected from each region, the number being
    #'   \eqn{\text{subsample\_proportions} \times \text{\# of samples in region with the fewest samples}}
    #' (rounded to the nearest integer).
    #' @param H \["auto"\] Controls the clustering cut heights \eqn{H} used in each run. There are three modes possible:
    #' - (Default) Random \eqn{H}, truncated normal distributions. Each iteration, \eqn{H} will be drawn from a truncated normal distribution, truncated to the values given by `H_truncNorm_range`. To use this option, you must provide the means of the distributions you want to try in `H`, and the standard deviations (one corresponding to each mean) in `H_truncNorm_sd`.
    #'   - If `H="auto"` (the default), then a selection of 16 normal distributions will be chosen based on a simple heuristic that typically captures some very good values. Briefly, the values of `H` fully span the interguartile ranges of inter-sample distances of each region, `H_truncNorm_sd` is half the standard deviation of inter-sample distances in the region with the smallest such value (i.e. it is given the same value in all 16 distributions), and `H_truncNorm_range` is just the range of all inter-sample distances in the whole dataset.
    #' - Fixed \eqn{H}. To use this option, provide a list of values to try as \eqn{H}. A run will be conducted with each, using the same value at each iteration.
    #' - Random \eqn{H}, empirical distribution. This samples values of \eqn{H} from the distances in the distance matrix. This should probably never be used, because while it does help produce meaningful plots in cases of extremely structured populations, the results are just not very intuitive. If you use this mode, you should clearly indicate it in any figure or publication, e.g. "ReMIXTURE plot using empirically distributed \eqn{H}". But just don't. This option will be removed.
    #' @param H_truncNorm_sd \[NULL\] See description for `H`.
    #' @param H_truncNorm_range \[upper_tri_ply(private$m,range)\] See description for `H`. By default it is set to the range of values in the distance matrix.
    #' @param diagnosticPlotMDSclusters \[FALSE\] A very useful diagnostic tool, especially for weirdly-structured datasets. If true, at each iteration, will draw the subsampled samples on an MDS plot, and draw hulls around the clusters.
    #'
    #' @return Nothing
    run = function(
      iterations=1000,
      subsample_proportions=c(0.8),
      H="auto",
      H_truncNorm_sd=NULL,
      H_truncNorm_range=upper_tri_ply(private$m,range),
      diagnosticPlotMDSclusters=FALSE
    ){
      ce("------------------------------------------------")
      ce("Running ReMixture analysis ...")
      ce("------------------------------------------------\n")

      if( diagnosticPlotMDSclusters==TRUE & interactive()==FALSE ){
        ce("Diagnostic MDS clustering plots are only available in interactive sessions--disabling this feature.")
        diagnosticPlotMDSclusters <- FALSE
      }

      ce("\tResetting results fields and flags")
      private$results <- list()
      private$runFlag <- FALSE

      ce("\tSetting up local variables and containers ...")
      #local params
      nits <- iterations
      ind_info <- data.table(
        gp=colnames(private$m),
        idx=1:ncol(private$m)
      )
      gp_vec <- colnames(private$m)
      gp_list <- sort(ind_info[,unique(gp)]) # check this lines up
      gp_info <- data.table(
        gp = sort(gp_list),
        gp_idx = 1:length(gp_list)
      )
      ind_info <- gp_info[ind_info,on="gp"]
      nind <- ncol(private$m)
      ngp <- nrow(gp_info)
      # param_test_insert <- 1:ngp
      results_insert <- 1

      H_type <- if(H[1]=="auto"){
        ce("\t\tCalculating some (hopefully) sensible parameters.")
        tmp <- apply(as.matrix(gp_list),1,function(gp){
          # gp<-"Africa"
          selRegionIdx <- which(colnames(private$m)==gp)
          utIdx <- upper.tri(private$m[selRegionIdx,selRegionIdx])
          matrix(
            c(
              quantile(private$m[selRegionIdx,selRegionIdx][utIdx],0.25),
              quantile(private$m[selRegionIdx,selRegionIdx][utIdx],0.75),
              # median(private$m[selRegionIdx,selRegionIdx][utIdx]),
              sd(private$m[selRegionIdx,selRegionIdx][utIdx])/2
            ),
            ncol=1
          )
        })

        h_min <- tmp[1,which(tmp[1,]==min(tmp[1,]))][1]
        h_max <- tmp[2,which(tmp[2,]==max(tmp[2,]))][1]
        H <- seq(h_min,h_max,length.out=16)
        H_truncNorm_sd <- min(tmp[3,which(tmp[3,]==min(tmp[3,]))][1])/2 %>% rep(16)

        ce("\tRandom H values following a selection of automatically chosen (truncated) normal distributions are requested.")
        ce("\t`H`=c(",paste0(signif(H,4),collapse=","),")")
        ce("\t`H_truncNorm_sd`=c(",paste0(signif(H_truncNorm_sd,4),collapse=","),")")
        "random-truncNormal"
      } else if(H[1]=="random-empirical"){
        "random-empirical"
      } else if (!is.null(H_truncNorm_sd)){
        ce("\tRandom H values following user-supplied (truncated) normal distributions are requested.")
        if(
          !(is.numeric(H) | is.integer(H)) |
          !(is.numeric(H_truncNorm_sd) | is.integer(H_truncNorm_sd)) |
          !(is.numeric(H_truncNorm_range) | is.integer(H_truncNorm_range)) |
          (length(H)!=length(H_truncNorm_sd)) |
          (length(H_truncNorm_range)!=2) |
          !(H_truncNorm_range[1]<H_truncNorm_range[2]) |
          any(H_truncNorm_sd<=0.0)
        ){
          stop("For random (truncated normal) H values, `H` and `H_truncNorm_sd` must be numerical or integers vectors of equal length, and `H_truncNorm_range` must be a length two vector specifying a range of numbers > 0.")
        }
        "random-truncNormal"
      } else {
        ce("\tConstant H requested.")
        if( !(is.numeric(H) | is.integer(H)) ){
          stop("H values must be numeric or integers.")
        }
        "constant"
      }


      # In case we are doing diagnostic MDS plots, get the coords
      if( diagnosticPlotMDSclusters==TRUE ){
        makePlot <- if(!is.null(private$mdsPlot)){
          if (all(private$mdsPlot$axes==c(1,2))){ FALSE } else { TRUE }
        } else { TRUE }

        if(makePlot==TRUE){
          colTable <- data.table(
            region = unique(colnames(private$m))
          )[,col:=colorRampPalette(c("#DD000088","#DDDD0088","#00DD0088","#0000DD88","#DD00DD88"))(.N)]
          mds <- cmdscale( as.dist(private$m) , k=2 )
          private$mdsPlot <- list(axes=c(1,2),mds=data.table(region=rownames(mds),axisA=mds[,1],axisB=mds[,2]),legend=colTable)
        }
      }

      for( pr_samp in subsample_proportions ){
        #dev pr_samp <- 0.8
        nSelectPerGrp <- round(ind_info[,.N,by=.(gp)][,min(N)] * pr_samp)
        if(nSelectPerGrp<2){
          private$results <- list()
          private$runFlag <- FALSE
          stop(paste0("Some regions do not have enough samples for this value of `subsample_proportions` (",subsample_proportions,"), which results in a subsampling number of ",nSelectPerGrp," samples. Increasing `subsample_proportions` could help but most likely you should think about omitting low-sample regions or combining them into other regions."))
        }
        ce("Each iteration, ",nSelectPerGrp," samples from each region will be used.\n")
        for(i_hcut in seq_along(H)){
          sayHc <- if(H_type=="random-empirical"){
            H_type
          } else if (H_type=="random-truncNormal"){
            paste0("truncNorm( mean=",H[i_hcut],", sd=",H_truncNorm_sd[i_hcut]," ,range=[",H_truncNorm_range[1],",",H_truncNorm_range[2],"] )")
          } else if (H_type=="constant"){
            round(H[i_hcut],digits = 4)
          } else {
            stop("Something is fundamentally wrong with the universe, email mtrw85@gmail.com and alert Tim.")
          }
          ce( "Begin analysis for H==" , sayHc , " and subsample_proportions==" , pr_samp , " ..." )

          #container for parameter testing output
          # param_test_out <- expand.grid( #innermost loops first
          #   gp_idx=gp_info$gp, #
          #   it=1:nits,
          #   hcut=H,
          #   pr_samp=subsample_proportions,
          #   nclust=integer(1),
          #   run=integer(1)
          # ) %>% setDT()

          #local result containers
          nclust_counts <- rep(0L,length(gp_list))
          nclust_counts_2 <- rep(0L,length(gp_list))
          counts_mat <- matrix(0L,nrow=ngp,ncol=ngp)
          colnames(counts_mat) <- rownames(counts_mat) <- gp_list
          counts_mat_accumulator_empty <- counts_mat
          total_population_diversity <- 0.0
          counts_2_mat <- matrix(0.0,nrow=ngp,ncol=ngp)
          colnames(counts_2_mat) <- rownames(counts_2_mat) <- gp_list
          ind_info[,clust:=NA_integer_]

          skipPlotUntil <- 1L
          ce("\tIterating ...")
          for(it in 1:nits){ #Begin iteration loop
            if(it %% 100 == 0) {
              ce("\t\tBegin iteration: ",it)
            }
            #Subsample
            ss_selector <- ind_info[,.(s=sample(idx,nSelectPerGrp)),by=.(gp_idx)]$s

            # set hc
            hc <- if(H_type=="random-empirical"){
              sample(private$m[ss_selector,ss_selector],1)
            } else if (H_type=="random-truncNormal"){
              truncnorm::rtruncnorm(1,H_truncNorm_range[1],H_truncNorm_range[2],H[i_hcut],H_truncNorm_sd[i_hcut])
            } else if (H_type=="constant"){
              H[i_hcut]
            } else {
              stop("Something is fundamentally wrong with the universe, email mtrw85@gmail.com and alert Tim.")
            }

            #Cluster
              #reset
            ind_info[,clust:=NULL]
            ind_info[ss_selector, clust:=cutree(hclust(as.dist(private$m[ss_selector,ss_selector])),h=hc)]


            # Diagnostic plots if requested
            if( diagnosticPlotMDSclusters==TRUE & skipPlotUntil==it)
            {
              plot(
                private$mdsPlot$mds[ss_selector,]$axisA,
                private$mdsPlot$mds[ss_selector,]$axisB,
                col=private$mdsPlot$legend[data.table(region=private$mdsPlot$mds[ss_selector,]$region),on=.(region)]$col,
                xlim=range(private$mdsPlot$mds$axisA),
                ylim=range(private$mdsPlot$mds$axisB),
                xlab="Axis 1",
                ylab="Axis 2",
                pch=20,
                cex=0.4,
                main=paste0( "Iteration " , it , "; Subsampling " , pr_samp*100,"%; H (",H_type,"): ", round(hc,2) )
              )

              for(cl in unique(ind_info[ss_selector,]$clust)){
                cl_selector <- which(ind_info$clust==cl)
                if(length(cl_selector)==1){
                  points(
                    private$mdsPlot$mds[cl_selector,]$axisA,
                    private$mdsPlot$mds[cl_selector,]$axisB,
                    col="#00000077",
                    pch=20,
                    cex=2
                  )
                } else if(length(cl_selector)==2){
                  lines(
                    private$mdsPlot$mds[cl_selector,]$axisA,
                    private$mdsPlot$mds[cl_selector,]$axisB,
                    col="#00000077",
                    lwd=11
                  )
                } else {
                  hull <- chull(private$mdsPlot$mds[cl_selector,]$axisA, private$mdsPlot$mds[cl_selector,]$axisB)
                  polygon(
                    private$mdsPlot$mds[cl_selector,][hull,]$axisA,
                    private$mdsPlot$mds[cl_selector,][hull,]$axisB,
                    col="#00000033"
                  )
                }
              }
              skipPlotUntil <- skipPlotUntil+1L
              ans <- ask(paste0("[Round ",it,"] Press <return> for the next round, enter an integer N to skip to the Nth round, or 'f' to stop showing clusters and finish all rounds."),YN=FALSE)
              if(!is.na(as.integer(ans)%>%suppressWarnings())){ skipPlotUntil<-(as.integer(ans)%>%suppressWarnings()); ce("You have requested another plot at round ",skipPlotUntil,"; You are currently at round ",it,".")}
              if(ans=="f"){ skipPlotUntil <- -1L; ce("No more plots will be drawn.") }
            }

            counts_mat_accumulator <- counts_mat_accumulator_empty
            nclust_counts <- nclust_counts + (t<-ind_info[ss_selector,.(add=nu(clust)),by=.(gp_idx)][gp_info,on=.(gp_idx)][is.na(add),add:=0L][]$add) # in how many unique clusters does each region occur
            nclust_counts_2 <- nclust_counts_2 + t**2

            # param_test_out[param_test_insert]$run <- results_insert
            # param_test_insert <- param_test_insert + ngp

            ind_info[ss_selector,{ #over clusters
              ugidx <- unique(gp_idx)
              if(length(ugidx) == 1){
                counts_mat_accumulator[gp_idx,gp_idx] <<- counts_mat_accumulator[gp_idx,gp_idx] + 1 #solo cluster--contributes to genetic 'uniqueness'
              } else {
                apply(combn(ugidx,2),2,function(c) { counts_mat_accumulator[c[1],c[2]] <<- counts_mat_accumulator[c[1],c[2]] + 1 } ) #over permutations of members of multigroup cluster
              }
            },by=.(clust)] %>% invisible

            total_population_diversity <- total_population_diversity + ind_info[ss_selector,nu(clust)]

            counts_mat_accumulator <- fold_matrix(counts_mat_accumulator)
            counts_mat <- counts_mat + counts_mat_accumulator
            counts_2_mat <- counts_2_mat + counts_mat_accumulator**2

          } #end iteration loop
          #browser()
          ce("\tSummarising and saving results ...")
          private$results[[results_insert]] <- list(
            subsample_proportion=pr_samp,
            H_type=H_type,
            H=if(H_type %in% c("constant","random-truncNormal") ){H[i_hcut]}else{NULL},
            H_truncNorm_sd=if(H_type %in% c("random-truncNormal") ){H_truncNorm_sd[i_hcut]}else{NULL},
            H_truncNorm_range=if(H_type %in% c("random-truncNormal") ){H_truncNorm_range}else{NULL},
            iterations=nits,
            total_population_diversity = total_population_diversity/nits,
            overlap = counts_mat/nits,
            var_overlap = (counts_2_mat/nits) - (counts_mat/nits)**2,
            diversity = nclust_counts/nits,
            var_diversity = (nclust_counts_2/nits) - (nclust_counts/nits)**2
          )

          results_insert <- results_insert + 1
        } # end hcut loop
      } # end pr_samp loop

      #browser()
      private$runFlag <- TRUE
      ce("\n------------------------------------------------")
      ce("ReMixture analysis complete ...")
      ce("------------------------------------------------")
    },


    #### PLOTTING ####

    #' @description
    #' Plot the unique and overlapped diversity recorded in a ReMIXTURE run, in something like a heatmap format but using the circle-area conventions as per `plot_maps()`.
    #'
    #' @param runs \[NULL\] For which runs would you like a plot? An integer vector. By default, all runs.
    #' @param ... Additional arguments ultimately passed to `base::plot()`
    #'
    #' @return Nothing
    plot_results_grid = function( runs=NULL ,...){
      if(private$runFlag==FALSE){
        stop("Analysis has not been run. Perform using `$run()`")
      }
      if(is.null(runs)){
        ce("Plotting a results grids for all runs.")
        runs <- 1:length(private$results)
      }
      for(i in runs){
        ssp <- private$results[[i]]$subsample_proportion
        hc <- private$results[[i]]$H
        hcsd <- private$results[[i]]$H_truncNorm_sd
        its <- private$results[[i]]$iterations

        ce(paste0("Plotting run ",i,":\n\tSubsample proportion: ",ssp,"\n\tH mean: ",round(hc,digits = 2),"\n\tH sd:",hcsd,"\n\tIterations:",its,"\n------------------------------\n"))

        td <- private$results[[i]]$diversity
        ol_ud <- private$results[[i]]$overlap
        maxDiv <- max(c(td,ol_ud))
        tdScale <- td/maxDiv
        ol_udScale <- ol_ud/maxDiv
        n <- m <- nrow(ol_ud)


        null_plot(1:(n+1),1:(m+1),xaxt="n",yaxt="n",main=paste0("Results grid for run ",i,"\nTotal/unique diversity (diagonal) &\noverlapped diversity (off-diagonal)"))#,...)
        axis(2,at=(1:n)+0.5,labels=rownames(ol_ud),las=2,cex.axis=0.6)
        axis(1,at=(1:n)+0.5,labels=colnames(ol_ud),las=3,cex.axis=0.6)
        mtext(paste0("(Maximum circle size=",maxDiv," clusters)"), side = 4,cex=.6)
        for(i in 1:n){
          for(j in 1:m){
            drawUnitSquareTopRight(i,j,col="#00000000",lwd=0.2)
            if(i==j){
              drawUnitCircleTopRight(i,j,scale=tdScale[i]     ,col="#000000")
              drawUnitCircleTopRight(i,j,scale=ol_udScale[i,j],col="#FFFFFF")
            } else {
              drawUnitCircleTopRight(i,j,scale=ol_udScale[i,j],col="#000000")
            }
          }
        }
      }
    },


    #' @description
    #' Plot the cluster counts for each region in each run. This is useful for assessing how different values of \eqn{H} affect the output, in particular checking there aren't regions with very low cluster counts, or that they are not "maxxed out", either of which will probably give dodgy results unreliable.
    #' @param plotOverlaps \[FALSE\] Also draw (faint) lines to record the cluster counts contributing to inter-region overlaps? Typically hard to see but useful if you fear some significant major observation may not be robust to parameter settings ...
    #' @return Nothing
    plot_clustercounts = function(
      plotOverlaps=FALSE
    ){
      #browser()
      if(private$runFlag==FALSE){
        stop("Analysis has not been run. Perform using `$run()`")
      }

      pd <- setDT(ldply(1:length(private$results),function(i){
        data.table(
          run=i,
          region=private$results[[i]]$overlap %>% rownames,
          cluster_count=private$results[[i]]$diversity
        )
      }))

      setkey(pd,run,region)
      pd[,xIdx:=1:.N,by=run]
      pd[,colRun:=rgb(run/length(private$results),0.0,0.0)]

      #i=1
      null_plot(pd$xIdx,pd$cluster_count,xaxt="n",ylab="Cluster count (run #)")
      axis(1,1:nu(pd$region),pd[run==1,region],las=2,cex.axis=0.6)
      #browser()
      for(i in unique(pd$run)){
        if(plotOverlaps==TRUE){
          for(r2 in pd[run==i]$region){
            lines(x=pd[run==i]$xIdx,y=private$results[[i]]$overlap[r2,],col=pd[run==i]$colRun,lwd=.4)
          }
        }
        lines(x=pd[run==i]$xIdx,y=pd[run==i]$cluster_count,col=pd[run==i]$colRun,type="b",pch=20)
      }
      pd[region==region[1],text(xIdx,cluster_count,paste0("(",run,")"),cex=0.5,pos=2)]

      invisible(NULL)
    },

    #' @description
    #' These counts are excellent for determining a good value for \eqn{H}, in cases where the default choice is not desirable. A good value typically gives a nice balance of multi- and single-region clusters, typically between where the two values cross over and (as occurs in most datasets) a point where the number of multi-region clusters hits some maximum.
    #'
    #' @return Nothing
    plot_h_optimisation = function(){
      if(private$runFlag==FALSE){
        stop("Analysis has not been run. Perform using `$run()`")
      }
      d <- setDT(ldply(private$results,function(r){
        data.table(
          subsample_proportion = r$subsample_proportion,
          H = r$H,
          H_type = r$H_type,
          iterations = r$iterations,
          clustercount_uniq = weighted.mean(diag(r$overlap),w=1/private$rt$N),
          clustercount_shared = weighted.mean(r$diversity-diag(r$overlap),w=1/private$rt$N) #or lower tri
        )
      }))[,run:=1:.N][]

      d[,{
        null_plot(run,c(clustercount_uniq,clustercount_shared),xaxt="n",ylab="Aggregated cluster count",xlab=paste0("Run\n(H) [mean; type=\"",H_type[1],"\"]"),)
        abline(v=run,lty=2,col="#00000044")
        axis(1,run,paste0(run,"\n(",signif(H,3),")."),padj = -0.2,cex.axis=0.6)
        lines(run,clustercount_uniq,type="b",col="#880000",pch=20)
        lines(run,clustercount_shared,type="b",col="#ffcb00",pch=20)
        legend(1,max(c(clustercount_uniq,clustercount_shared)),c("Single-region clusters","Multi-region clusters"),fill=c("#880000","#ffcb00"))
      }]
      invisible(NULL)
    },

    #' @description
    #' Produces the canonical ReMIXTURE plots.
    #'
    #' @param run \[NULL\] If multiple runs were done, choose which run to use. This is best assessed using `<ReMIXTURE Object>$plot_h_optimisation()`, `<ReMIXTURE Object>$plot_results_grid()`, and `<ReMIXTURE Object>$plot_clustercounts()`.
    #' @param focalRegion \[NULL\] A character vector naming regions to plot as the focal region (from which lines showing inter-region overlap will emanate). If left NULL then a map will be plotted with each region as the focus in turn--this is best used after setting up a multi-panel plot with some variant of `par(mfrow=c(<number of rows>,<number of columns>))`. By default, all regions will be plotted sequentially.
    #' @param range_lon \[c(-179.0,179.0)\] Limit the map to some range of longitudes. Must be a vector of 2 numbers. For awkward reasons, using 180 or -180 can cause graphical bugs and I am deeply sorry.
    #' @param range_lat \[c(-85.0,85.0)\] Limit the map to some range of latitudes. Must be a vector of 2 numbers.
    #' @param width_max \[15.0\] The maximum width of the circles/lines, in units of lat/lons. This width will correspond to the highest cluster count and the others will be scaled such that zero clusters <=> zero width.
    #' @param alpha_max \[1.0\] As per width_max but controlling the alpha of connecting lines. Setting to NULL will disable alpha and make all lines solid. Can be set above 1.0, with weird results--probably don't do this.
    #' @param diversityCirclesFocalOnly \[FALSE\] Plot the circle representing a region's total/unique/shared diversity at only the focal region. Otherwise, all plots will show the diversities at all regions. This is good when you have multiple plots on the page, see description for `focalRegion`.
    #' @param allOnOneMap \[FALSE\] If TRUE, plot all regional diversity circles and overlap lines on the same map. This is usually far too busy to be very good, but can be good for datasets with not too many regions.
    #' @param projection \[EckertIV\] A function the performs projection of the map coordinates. Included in the package are `eckertIV`, `winkelIII`, and `equirectangular`. A custom function can be given--see details.
    #' @param curvature \[NULL\] Describes how the lines emanating from the focal region should bend to aid visualisation.
    #' One of the following options:
    #' * To automatically deflect all the angles a certain number of radians clockwise from the simple straight line path (that is, from a straight line on the equirectangular projection), provide a single integer/numeric value. This creates nice 'pinwheel' type plots.
    #' * If the above still has pesky overlaps, add another bit of deflection proportional to the (equirectangular euclidean) geographic distance between regions, provide a length two numeric/integer vector. The first element provides the baseline deflection, the second a factor by which the distance should be multiplied to provide the subsequent deflection, where the longest distance in the dataset is given distance 1 and the shortest is given distance 0. Try, e.g. `c(1,-1)`.
    #' * If these fail to please (or if you are just chaotic-aligned), the argument "random" will auto-generate a matrix filled with entries ~ normal(0,0.3), which usually works well.
    #' * If you really want full control, provide a numeric matrix where the (i,j)th entry gives the deflection of the line beginning at region i and ending at region j (indexed as per the order in the region table).
    #' @param mapData \[mapData110\] The polygon data from which the map is made. For higher resolution, use mapData10. See e.g. ?mapData110.
    #'
    #' @details
    #' A custom projection function should follow these rules: its first argument, named 'dtLL' should take a data.table with columns 'lon' and 'lat'. The second argument, named 'projColNames', should take a length-2 character vector. The output should be a data.table with columns containing the transformed longitude and latitude values. These columns should be named as per the entries in the input 'projColNames'.
    #'
    #'
    #' @return Nothing
    plot_maps = function(
      run=NULL,
      focalRegion=NULL,
      range_lon=c(-179.0,179.0),
      range_lat=c(-85.0,85.0),
      width_max=10.0,
      alpha_max=1.0,
      diversityCirclesFocalOnly=FALSE,
      allOnOneMap = FALSE,
      projection=eckertIV,
      curvature=NULL,
      mapData=mapData110
    ){
      if(private$runFlag==FALSE){
        stop("Analysis has not been run. Perform using `$run()`")
      }
      if(is.null(run) & length(private$results)!=1){
        stop("Please provide a valid run number to plot from (consider using `<ReMIXTURE Object>$plot_h_optimisation()`, `<ReMIXTURE Object>$plot_results_grid()`, and `<ReMIXTURE Object>$plot_clustercounts()` to assess which run parameters are appropriate).")
      }
      if(is.null(run) & length(private$results)==1){
        run <- 1
      }
      if(!is.null(focalRegion) && any(!focalRegion %in% private$rt$region)){
        stop("Entries given to `focalRegion` must be in the region table")
      }

      st <- private$results[[run]]$overlap
      rt <- data.table(
        region=colnames(st),
        totDiv=private$results[[run]]$diversity
      )

      rt[,uniqueDiv:=private$results[[run]]$overlap[r,r],by=.(r=region)] # could just pull out diagonal but this seems safer
      rt <- private$rt[rt,on=.(region)]
      width_scaling <- compute_plot_width_scaling(
        totDiv = rt$totDiv,
        uniqueDiv = rt$uniqueDiv,
        overlap = st,
        width_max = width_max
      )
      rt[,wTotDiv := width_scaling$wTotDiv]
      rt[,wUniqueDiv := width_scaling$wUniqueDiv]
      wst <- width_scaling$wst

      ct <- if(is.null(curvature)){
        matrix(0.0,ncol=nrow(private$rt),nrow=nrow(private$rt),dimnames=list(private$rt$region,private$rt$region))
      } else if (curvature[1]=="random"){
        matrix(rnorm(nrow(private$rt)**2,0,0.3),nrow=nrow(private$rt),dimnames=list(private$rt$region,private$rt$region))
      } else if (length(curvature)==1 && (is.numeric(curvature) | is.integer(curvature))){
        matrix(curvature,ncol=nrow(rt),nrow=nrow(rt),dimnames=list(private$rt$region,private$rt$region))
      } else if (length(curvature)==2 && (is.numeric(curvature) | is.integer(curvature))){
        geoDistMat <- outer(1:nrow(rt),1:nrow(rt),function(i,j){
          euc_dist(rt[i,]$lon,rt[i,]$lat,rt[j,]$lon,rt[j,]$lat)
        })
        geoDistMat[,] <- geoDistMat %>% scale_between(0,1)
        matrix(curvature[1],ncol=nrow(private$rt),nrow=nrow(private$rt),dimnames=list(private$rt$region,private$rt$region)) + curvature[2]*geoDistMat
      } else {
        curvature
      }

      at <- copy(st)
      diag(at) <- 0.0
      if(is.null(alpha_max)){
        at[,] <- 1.0
      } else {
        maxAlphaSource <- suppressWarnings(max(at, na.rm = TRUE))
        if(!is.finite(maxAlphaSource) || maxAlphaSource <= 0){
          at[,] <- 0.0
        } else {
          at[,] <- at/maxAlphaSource*alpha_max
        }
      }

      plotMiddle <- findCentreLL(range_lon,range_lat)
      trt <- copy(rt) %>% rotateLatLonDtLL(-plotMiddle[1],-plotMiddle[2],splitPlotGrps=F)

      rIdxList <- if(is.null(focalRegion)){
        1:nrow(trt)
      } else {
        which(rt$region %in% focalRegion)
      }

      if(allOnOneMap!=TRUE){
        ### SEPARATE MAPS PER FOCAL REGIONS CASE ###

        for(i in rIdxList){ # Over focal regions

          pe <- plotEmptyMap( range_lon, range_lat, projFun=projection )
          plotMapItem( makeMapDataLatLonLines() , range_lon , range_lat, projFun=projection , plotFun=lines,   col="#00000022" , lwd=.4 )
          plotMapItem( mapData                  , range_lon , range_lat, projFun=projection , plotFun=polygon, col="#f7bf25" , lwd=0.2 )
          plotMapBorder(                          range_lon , range_lat, projFun=projection , plotEdges=pe ,lwd=4 )

          for(j in 1:nrow(trt)){ # Over target regions
            if(i==j){next}
            ldt <- curved_rounded_line(
              x1=trt$lon[i],y1 = trt$lat[i],
              x2=trt$lon[j],y2 = trt$lat[j],
              width =     wst[trt$region[i],trt$region[j]],
              curvature = ct[trt$region[i],trt$region[j]]
            ) %>% mat2dtLL()
            plotMapItem(ldt,projFun=projection,plotFun=polygon,col=alpha("black",at[trt$region[i],trt$region[j]]),border="#000000",lwd=0.15)
            if(diversityCirclesFocalOnly==FALSE){
              cdt <- circle_seg(trt[j]$lon,trt[j]$lat,radius=trt$wTotDiv[j]/2   ) %>% mat2dtLL()
              udt <- circle_seg(trt[j]$lon,trt[j]$lat,radius=trt$wUniqueDiv[j]/2) %>% mat2dtLL()
              plotMapItem(cdt,projFun=projection,plotFun=polygon,col="#00000055") # 'shadow' effect
              plotMapItem(cdt,projFun=projection,plotFun=polygon,col="#000000FF")
              plotMapItem(udt,projFun=projection,plotFun=polygon,col="#FFFFFF")
            }
          }
          cdt <- circle_seg(trt[i]$lon,trt[i]$lat,radius=trt$wTotDiv[i]/2   ) %>% mat2dtLL()
          udt <- circle_seg(trt[i]$lon,trt[i]$lat,radius=trt$wUniqueDiv[i]/2) %>% mat2dtLL()
          plotMapItem(cdt,projFun=projection,plotFun=polygon,col="#000000FF")
          plotMapItem(udt,projFun=projection,plotFun=polygon,col="#FFFFFF")
          title(main=trt$region[i],line=0)
        }
      } else {
        ### ALL OVERLAPS ON ONE MAP CASE ###

        pe <- plotEmptyMap( range_lon, range_lat, projFun=projection )
        plotMapItem( makeMapDataLatLonLines() , range_lon , range_lat, projFun=projection , plotFun=lines,   col="#00000022" , lwd=.4 )
        plotMapItem( mapData                  , range_lon , range_lat, projFun=projection , plotFun=polygon, col="#f7bf25" , lwd=0.2 )
        plotMapBorder(                          range_lon , range_lat, projFun=projection , plotEdges=pe ,lwd=4 )

        # Draw connecting lines
        for(i in 1:(nrow(trt)-1)){ # From region
          for(j in (i+1):nrow(trt)){ # To region
            # dev i<-2;j<-7;
            ldt <- curved_rounded_line(
              x1=trt$lon[i],y1 = trt$lat[i],
              x2=trt$lon[j],y2 = trt$lat[j],
              width =     wst[trt$region[i],trt$region[j]],
              curvature = ct[trt$region[i],trt$region[j]]
            ) %>% mat2dtLL()
            plotMapItem(ldt,projFun=projection,plotFun=polygon,col=alpha("black",at[trt$region[i],trt$region[j]]),border="#000000",lwd=0.15)
          }
        }
        # Draw circles
        for(i in 1:nrow(trt)){
          cdt <- circle_seg(trt[i]$lon,trt[i]$lat,radius=trt$wTotDiv[i]/2   ) %>% mat2dtLL()
          udt <- circle_seg(trt[i]$lon,trt[i]$lat,radius=trt$wUniqueDiv[i]/2) %>% mat2dtLL()
          plotMapItem(cdt,projFun=projection,plotFun=polygon,col="#00000055") # 'shadow' effect
          plotMapItem(cdt,projFun=projection,plotFun=polygon,col="#000000FF")
          plotMapItem(udt,projFun=projection,plotFun=polygon,col="#FFFFFF")
        }
      }
    },
    #' @description
    #' Plot inter-sample distances as a density plot over the distance matrix--basically, a convenient way to judge how far apart samples tend to be, and thus how the clustering might behave when various H values (\eqn{H}) are used.
    #'
    #' Each plot tells the story from the perspective of a focal region. Each region is assigned a colour that remains constant across plots. A region's samples' distances to each other are marked with a thick line allowing you to see what colour is for what region.
    #'
    #' @param set_bw \[NULL\] The bandwidth parameter 'bw' passed to `density()`. By default it will be set automatically.
    #' @param set_xlims \[NULL\] Limits on the x-axis. Default is set by `density()`.
    #' @param samePlot \[FALSE\] If TRUE, it will put all the plots on one page.
    #' @param H \[NULL\] Visualise chosen values of (or distribution over) \eqn{H} on the plot. If only this argument is used, a single line will indicate that H value.
    #' @param H_truncNorm_sd \[NULL\] Visualise a given distribution of \eqn{H} on the plot. If this argument is used (in conjunction with `H`), a truncated normal distribution will be shown.
    #' @param H_truncNorm_lims \[upper_tri_ply(private$m,range)\] Visualise your values of \eqn{H} on the plot. If this argument is used, it will set the truncation limits of the normal distribution. By default it will be the range of inter-sample distances in the distance matrix.
    #' @param HdistFromRun \[NULL\] Show the H distribution for a particular run. Obviously, this requires that a run has been done.
    #' @param plotLegend \[FALSE\] If true, the last plot produced will be the colour legend.
    #' @param whichRegions \[private$rt$region\] Plot only chosen regions. A character vector giving the names.
    #'
    #' @return Invisibly returns a data.table giving the legend colours.
    plot_distance_densities = function(
      set_bw=0.001,
      set_xlims=upper_tri_ply(private$m,range),
      samePlot=FALSE,
      H=NULL,
      H_truncNorm_sd=NULL,
      H_truncNorm_lims=NULL,
      HdistFromRun=NULL,
      plotLegend=FALSE,
      whichRegions=private$rt$region
    ){

      if(!is.null(HdistFromRun)){
        if(private$runFlag==FALSE){
          stop("Analysis has not been run. Perform using `$run()`")
        }
        if(HdistFromRun==TRUE){ HdistFromRun <- 1L }
        ce("Plotting H distribution from run ",HdistFromRun)
        H <- private$results[[HdistFromRun]]$H
        H_truncNorm_sd <- private$results[[HdistFromRun]]$H_truncNorm_sd
        H_truncNorm_range <- private$results[[HdistFromRun]]$H_truncNorm_range
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
        col     = rgb( t(col2rgb(hsv(seq(0, 0.8, length.out = nu(colnames(private$m))), 1, 1)) / 255) )
      )
      dt1 <- colTable[dt1,on=.(region2)]
      dt1[,col:=rgb(t(col2rgb(col))/255,alpha=fifelse(region1==region2,1,0.4))]
      dt1[,lwd:=fifelse(region1==region2,3,1)]

      if(samePlot==TRUE){
        warning("Altering plot parameters for multi-plot graphics. They will be reset to the defaults, i.e. `mfrow=c(1,1)` and `mar=c(5.1,4.1,4.1,2.1)`.")
        par( mfrow=c(ceiling(sqrt(nrow(private$rt))),ceiling(sqrt(nrow(private$rt)))), mar=c(2,2,4,1) )
      }

      rangeY <- range(dt1$y)
      for(r1 in whichRegions){
        sel1 <- rownames(private$m)==r1
        null_plot(set_xlims,rangeY,main=paste0(r1))
        for(r2 in unique(colnames(private$m))){
          sel2 <- rownames(private$m)==r2
          dt1[ region1==r1 & region2==r2 , lines(x,y,col=col,lwd=lwd) ]
        }
        if(!is.null(H)){
          if(!is.null(H_truncNorm_sd)){
            if(is.null(H_truncNorm_lims)){
              ce("Argument `H_truncNorm_lims` not provided. Will be set to +/- Inf.")
              H_truncNorm_lims <- c(-Inf,Inf)
            }
            #require(truncnorm)
            polygon(
              x=c(
                seq(set_xlims[1],set_xlims[2],l=1000L),
                set_xlims[2],
                set_xlims[1]
              ),
              y=c(
                ((t<-truncnorm::dtruncnorm(seq(set_xlims[1],set_xlims[2],l=1000L),mean=H,sd=H_truncNorm_sd,a=H_truncNorm_lims[1],b=H_truncNorm_lims[2]))/max(t)) * rangeY[2],
                0,
                0
              ),
              border="#00000000",
              col="#00000022"
            )
          } else {
            abline(v=H,lty=2,col="#00000033")
          }
        }
      }
      if(samePlot==TRUE){
        par( mfrow=c(1,1), mar=c(5.1,4.1,4.1,2.1) )
      }

      if(plotLegend==TRUE){
        null_plot(0:2,0:nrow(colTable),xaxt="n",yaxt="n")
        for(i in nrow(colTable):1){
          rect(0,i-1,1,i,col=colTable[i]$col)
          text(x=1.0,y=i-0.5,label=colTable[i]$region2,pos=4,cex=0.8)
        }
      }
      invisible(return(colTable))
    },
    #' @description
    #' An MDS plot showing relationships among samples
    #'
    #' An MDS ordination plot (in the same flavour as a PCA--in fact in some circumstances they will be the same depending on what your distances represent), useful as a first-pass assessment of the relationships among your samples--what overlaps what, what is over/underrepresented in the sample set, etc.
    #'
    #' These can take a very long time if the number of samples is high. In these cases, I suggest taking a smaller sample of the whole dataset.
    #'
    #' @param axes \[c(1,2)\] A vector of the MDS axes that should be outputted.
    #' @param colPalette \[c("#DD000088","#DDDD0088","#00DD0088","#0000DD88","#DD00DD88")\] A vector collection of colours that will be interpolated and applied to the different
    #' @param showLegend \[TRUE\] Add a legend?
    #' @param doPlot \[TRUE\] Produce a plot? If not, will still return the data for the MDS.
    #' @param pch \[20\] The shape of the points. See `?points`.
    #' @param cex \[0.8\] The size of the points. See `?points`.
    #' @param ... Other arguments passed to `plot()`.
    #'
    #' @return Invisibly returns a 2-list containing the (1) results of the MDS in a data.table and (2) information about the legend (another data.table).
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
      colTable <- data.table(
        region = unique(colnames(private$m))
      )[,col:=colorRampPalette(colPalette)(.N)]
      dim <- max(axes)

      makePlot <- if(!is.null(private$mdsPlot)){
        if (all(private$mdsPlot$axes==axes)){ FALSE } else { TRUE }
      } else { TRUE }

      if(makePlot==TRUE){
        mds <- cmdscale( as.dist(private$m) , k=dim )  # Also has function of ignoring the Inf diagonals (done for reasons to make the main algorithm work)
        private$mdsPlot <- list(axes=axes,mds=data.table(region=rownames(mds),axisA=mds[,axes[1]],axisB=mds[,axes[2]]),legend=colTable)
      }

      if(doPlot==TRUE){
        plot(
          private$mdsPlot$mds$axisA,
          private$mdsPlot$mds$axisB,
          pch=pch,
          col=colTable[data.table(region=private$mdsPlot$mds$region),on=.(region)]$col,
          cex=cex,
          xlab=paste0("Axis ",axes[1]),
          ylab=paste0("Axis ",axes[2]),
          ...
        )
        if(showLegend==TRUE){
          legend(min(private$mdsPlot$mds$axisA),max(private$mdsPlot$mds$axisB),colTable$region,colTable$col,bg="#FFFFFFAA")
        }
      }

      invisible(return(private$mdsPlot))
    }



  )
)

