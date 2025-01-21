


ReMIXTURE$set( "public" , "run" ,
  function(
    iterations=2500,
    subsample_proportions=c(1.0),
    h_cutoffs,
    diagnosticPlotMDSclusters=FALSE,
    ...
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
    if(length(h_cutoffs)==1){
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

    # In case we are doing diagnostic MDS plots, get the coords
    if( diagnosticPlotMDSclusters==TRUE ){
      mds <- self$plot_MDS(doPlot = FALSE)
    }

    #container for parameter testing output
    param_test_out <- expand.grid( #innermost loops first
      gp_idx=gp_info$gp, #
      it=1:nits,
      hcut=h_cutoffs,
      pr_samp=subsample_proportions,
      nclust=integer(1),
      run=integer(1)
    ) %>% setDT()


    param_test_insert <- 1:ngp
    results_insert <- 1
    for( pr_samp in subsample_proportions ){
      nSelectPerGrp <- round(ind_info[,.N,by=.(gp)][,min(N)] * pr_samp)
      for(hcut in h_cutoffs){
        #dev hcut=.013; pr_samp=0.8
        ce( "Begin analysis for h_cutoff==" , round(hcut,digits=4) , " and subsample_proportions==" , pr_samp , " ..." )

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
          #Cluster

          ind_info[ss_selector, clust:=cutree(hclust(as.dist(private$m[ss_selector,ss_selector])),h=hcut)]
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
              main=paste0( "Iteration " , it , "; Subsampling " , pr_samp*100,"%; H-cutoff ", round(hcut,2) )
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
          nclust_counts <- nclust_counts + (t<-ind_info[ss_selector,.(add=nu(clust)),by=.(gp_idx)][order(gp_idx),]$add) # in how many unique clusters does each region occur
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

          # gp_presabs <- ind_info[ss_selector,.(gp_idx=1:ngp,presence=1:ngp %in% unique(gp_idx)),by=.(clust)]
          # setkey(gp_presabs,clust,gp_idx)
          rowSums(counts_mat_accumulator) - nclust_counts


        } #end iteration loop

        ce("\tSummarising and saving results ...")

        private$results$runs[[results_insert]] <- list(
          subsample_proportion=pr_samp,
          h_cutoff=hcut,
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

