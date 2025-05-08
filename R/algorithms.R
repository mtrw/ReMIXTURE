
ReMIXTURE$set( "public" , "run" ,
#' run
#' @name run
#' @description
#' Run the ReMIXTURE algorithm and save the results in the object.
#'
#' @param iterations \[1000\] The number of times subsamples are drawn, clustered, and counted.
#' @param subsample_proportions \[0.8\] size of subsample to draw each iteration. An equal number of samples will be selected from each region, the number being
#'   \eqn{\text{subsample\_proportions} \times \text{\# of samples in region with the fewest samples}}
#' (rounded to the nearest integer).
#' @param h_cutoff \["auto"\] Controls the clustering cutoff \eqn{H}. There are three modes possible:
#' - (Default) Random \eqn{H}, truncated normal distribution. Each iteration, \eqn{H} will be drawn from a normal distribution, truncated to the values given by `h_cutoff_normal_truncRange`. To use this option, you must provide the means of the distributions you want to try in `h_cutoff`, and the standard deviations (one for each mean) in `h_cutoff_normal_sigma`.
#'   - If `h_cutoff="auto"` (the default), then a single normal distribution will be used, with the same mean and sd as the distances in the distance matrix, equivalent to `remixture_obj$run(h_cutoff=mean(remixture_obj$distance_matrix),h_cutoff_normal_sigma=sd(remixture_obj$distance_matrix))`
#' - Fixed \eqn{H}. To use this option, provide a list of values to try as \eqn{H}. A run will be conducted with each, using the same value at each iteration. If you use this mode, you should clearly indicate it in any figure or publication, e.g. "ReMIXTURE plot using truncated random normal \eqn{H} (mean = ..., sd = ..., range=\[...,...\])".
#' - Random \eqn{H}, empirical distribution. This samples values of \eqn{H} from the distances in the distance matrix. This should probably never be used, because while it does help produce meaningful plots in cases of extremely structured populations, the results are just not very intuitive. If you use this mode, you should clearly indicate it in any figure or publication, e.g. "ReMIXTURE plot using empirically distributed \eqn{H}".
#' @param h_cutoff_normal_sigma \[NULL\] See description for `h_cutoff`.
#' @param h_cutoff_normal_truncRange \[range(self$distance_matrix)\] See description for `h_cutoff`.
#' @param diagnosticPlotMDSclusters \[FALSE\] A diagnostic tool. If true, at each iteration, will draw the subsampled samples on an MDS plot, and draw hulls around the clusters.
#'
#' @return NULL
#'
#' @examples
#'(1+1)==2
               function(
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
                 if(length(h_cutoff)==1){
                   warning("Only one h_cutoff value being tried. It is highly recommended to experiment with this parameter and select a value that gives meaningful results, i.e., where the cluster size captures a good balance of region-unique and region-overlapping clusters. Functions to help with this include plot_clustercount_diag_nondiag_means(), plot_heatmaps(), plot_clustercounts()")
                 }
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

                 # If random h-cutoffs are being used, alter the name for the output HHHHHHERE
                 h_cutoff_type <- if(h_cutoff[1]=="auto"){
                   h_cutoff <- mean(self$distance_matrix)
                   h_cutoff_normal_sigma <- sd(self$distance_matrix)
                   ce("Random H-cutoffs following an automatically chosen (truncated) normal distribution are requested.\n
                      `h_cutoff`=",h_cutoff,"; `h_cutoff_normal_sigma`=",h_cutoff_normal_sigma)
                   "random-truncNormal"
                 } else if(h_cutoff[1]=="random-empirical"){
                   "random-empirical"
                 } else if (!is.null(h_cutoff_normal_sigma)){
                   ce("Random H-cutoffs following user-supplied (truncated) normal distributions are requested.")
                   if(
                      !(is.numeric(h_cutoff) | is.integer(h_cutoff)) |
                      !(is.numeric(h_cutoff_normal_sigma) | is.integer(h_cutoff_normal_sigma)) |
                      !(is.numeric(h_cutoff_normal_truncRange) | is.integer(h_cutoff_normal_truncRange)) |
                      (length(h_cutoff)!=length(h_cutoff_normal_sigma) | length(h_cutoff)!=1) |
                      (length(h_cutoff_normal_truncRange)!=2) |
                      !(h_cutoff_normal_truncRange[1]<h_cutoff_normal_truncRange[2]) |
                      any(h_cutoff_normal_sigma<=0.0)
                     ){
                     stop("For random (truncated normal) cutoffs, `h_cutoff` and `h_cutoff_normal_sigma` must be numerical or integers vectors of equal length, and `h_cutoff_normal_truncRange` must be a length two vector specifying a range of numbers > 0.")
                   }
                   "random-truncNormal"
                 } else {
                   ce("Preset constant h-cutoffs requested.")
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
                   nSelectPerGrp <- round(ind_info[,.N,by=.(gp)][,min(N)] * pr_samp)
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

                       #browser("HELLO! 1")
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
               }
)
