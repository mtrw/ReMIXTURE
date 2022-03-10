






ReMIXTURE$set( "public" , "run" ,
  function(iterations=2500,subsample_proportion=0.8,h_cutoff=private$hcut){
    ce("------------------------------------------------")
    ce("Running ReMixture analysis ...")
    ce("------------------------------------------------\n")

    ce("\tSetting up local variables and containers ...")
    #local params
    nits <- iterations
    pr_samp <- subsample_proportion
    hcut <- h_cutoff

    if(is.null(hcut)){
      stop("h_cutoff not provided or set up in advance. Please provide one, preferebly based on the insights of `$test_h_cutoffs()`")
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

    #local result containers
    nclust_counts <- rep(0L,length(gp_list))
    nclust_counts_2 <- rep(0L,length(gp_list))
    counts_mat <- matrix(0L,nrow=ngp,ncol=ngp)
    colnames(counts_mat) <- rownames(counts_mat) <- gp_list
    counts_mat_accumulator_empty <- counts_mat
    counts_2_mat <- matrix(0.0,nrow=ngp,ncol=ngp)
    colnames(counts_2_mat) <- rownames(counts_2_mat) <- gp_list

    ce("\tIterating ...")
    for(it in 1:nits){ #dev it = 1
      if(it %% 100 == 0) {
        ce("\t\tBegin iteration: ",it)
      }
      ss_selector <- ind_info[sample(idx,round(pr_samp*.N))]$idx
      ind_info[ss_selector, clust:=cutree(hclust(as.dist(private$m[ss_selector,ss_selector]),),h=hcut)]
      counts_mat_accumulator <- counts_mat_accumulator_empty
      nclust_counts <- nclust_counts + (t<-ind_info[ss_selector,.(add=nu(clust)),by=.(gp_idx)][order(gp_idx),]$add)
      nclust_counts_2 <- nclust_counts_2 + t**2
      rm(t)

      ind_info[ss_selector,{ #over clusters
        ugidx <- unique(gp_idx)
        if(length(ugidx) == 1){
          counts_mat_accumulator[gp_idx,gp_idx] <<- counts_mat_accumulator[gp_idx,gp_idx] + 1 #solo cluster
        } else {

          apply(combn(ugidx,2),2,function(c) { counts_mat_accumulator[c[1],c[2]] <<- counts_mat_accumulator[c[1],c[2]] + 1 } ) #over permutations of members of multigroup cluster
        }
      },by=.(clust)] %>% invisible

      counts_mat_accumulator <- fold_matrix(counts_mat_accumulator)
      counts_mat <- counts_mat + counts_mat_accumulator
      counts_2_mat <- counts_2_mat + counts_mat_accumulator**2
    }
    ce("Summarising and saving results ...")

    private$overlap <- counts_mat / nits
    private$var_overlap <- (counts_2_mat/nits) - (counts_mat/nits)**2
    private$diversity <- nclust_counts/nits
    private$var_diversity <- (nclust_counts_2/nits) - (nclust_counts/nits)**2
    private$runflag <- TRUE

    ce("Updating parameter records ...")
    private$nits <- nits
    private$hcut <- hcut
    private$pr_samp <- pr_samp

    ce("\n------------------------------------------------")
    ce("ReMixture analysis complete ...")
    ce("------------------------------------------------")
  }
)



























# ReMIXTURE$set( "public" , "select_h_cutoff" ,
#                function(){
#                  out_dm <- copy(private$m)
#                  diag(out_dm) <- 0.0
#                  return(out_dm)
#                }
# )
