


ReMIXTURE$set( "public" , "plot_heatmaps" ,
  function(colpalette=colorRampPalette(c("#f2f5ff","#214feb","#001261"))(100),...){
    if(private$runflag==FALSE){
      stop("Analysis has not been run. Perform using `$run()`")
    }
    for(i in 1:(length(private$results$runs))){
      ssp <- private$results$runs[[i]]$subsample_proportion
      hc <- private$results$runs[[i]]$h_cutoff
      its <- private$results$runs[[i]]$iterations
      wait( paste0("Please press [ENTER] to produce the next heatmap (Run: ",i,"; Subsample proportion: ",ssp,"; H cutoff: ",hc,"; Iterations: ",its,")") )
      pheatmap::pheatmap(private$results$runs[[i]]$overlap,cluster_rows=F,cluster_cols=F,color=colpalette,main=paste0("Run: ",i,"\nSubsample proportion: ",ssp,"; H cutoff: ",hc,"; Iterations: ",its),...)
    }
  }
)








ReMIXTURE$set( "public" , "plot_clustercounts" ,
  function(colpalette=colorRampPalette(c("#f2f5ff","#214feb","#001261"))(100),...){
    if(private$runflag==FALSE){
      stop("Analysis has not been run. Perform using `$run()`")
    }
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
  }
)













ReMIXTURE$set( "public" , "plot_h_optimisation" ,
  function(colpalette=colorRampPalette(c("#f2f5ff","#214feb","#001261"))(100),plot_entropy=FALSE,...){
    if(private$runflag==FALSE){
      stop("Analysis has not been run. Perform using `$run()`")
    }
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

    if(!plot_entropy){
      print(ggplot(pdat_s,aes(x=h_cutoff,y=value,colour=variable)) +
        geom_line() +
        geom_point() +
        facet_grid(subsample_proportion~.) +
        theme_classic() +
        geom_hline(aes(yintercept=0)) +
        labs( colour = "Cluster counts:" ) +
        scale_color_manual(labels = c("Single-region clusters", "Multi-region clusters"), values = c("#11888a", "#c94d4d")) +
        ylab("Count") +
        xlab("h-cutoff") +
        geom_text(aes(label=run)))
    } else {
      print(ggplot(pdat_e,aes(x=h_cutoff,y=value)) +
        geom_line() +
        geom_point() +
        facet_grid(subsample_proportion~.) +
        theme_classic() +
        ylab("Overlap cluster count matrix entropy") +
        xlab("h-cutoff") +
        geom_text(aes(label=run)))
    }
  }
)









ReMIXTURE$set( "public" , "plot_maps" ,
  function(run=NULL,focalRegion=NULL,range_lon=c(-179,179),range_lat=c(-85,85),width_lims=c(5,13),alpha_lims=c(0.05,0.99),projection=eckertIV,curvature_matrix=NULL){
    if(private$runflag==FALSE){
      stop("Analysis has not been run. Perform using `$run()`")
    }
    if(is.null(run) & length(private$results$runs)>1){
      stop("Please provide a run number to plot from (consider using plot_heatmaps() and plot_clustercounts() to assess which run parameters are appropriate).")
    }
    if(is.null(run) & length(private$results$runs)==1){
     warning("Only one run was performed, and will be plotted. If you haven't already, please be sure to try multiple runs with a good range of parameter values, and assess their appropriateness with plot_heatmaps() and plot_clustercounts()")
     run <- 1
    }

    st <- private$results$runs[[run]]$overlap
    rt <- data.table(
      region=colnames(st),
      totDiv=private$results$runs[[run]]$diversity
    )
    rt[,uniqueDiv:=private$results$runs[[run]]$overlap[r,r],by=.(r=region)]
    rt <- private$rt[rt,on=.(region)]
    divRange <- range(rt$totDiv,rt$uniqueDiv,st)

    rt[,wTotDiv:=(c(totDiv,divRange) %>% scale_between(width_lims[1],width_lims[2]))[1:.N]]
    rt[,wUniqueDiv:=(c(uniqueDiv,divRange) %>% scale_between(width_lims[1],width_lims[2]))[1:.N]]
    wst <- st
    wst[,] <- c(st,range(st)) %>% scale_between(width_lims[1],width_lims[2]) %>% head(-2)

    ct <- if(is.null(curvature_matrix)){
      matrix(0.0,ncol=nrow(rt),nrow=nrow(rt),dimnames=list(rt$region,rt$region))
    } else {
      curvature_matrix
    }

    at <- st
    diag(at) <- mean(at) # just any value not at the extreme, so it doesn't affect the scaling of the shared div values
    at[,] <- at %>% scale_between(alpha_lims[1],alpha_lims[2])

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
        plotMapItem(ldt,projFun=projection,plotFun=polygon,col=alpha("black",at[trt$region[i],trt$region[j]]),border="black",lwd=1)

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
    }


  }
)






































ReMIXTURE$set( "public" , "plot_distance_densities" ,
  function(set_bw=0.001,set_xlims=c(0,1)){
    dt1 <- ldply(unique(colnames(private$m)),function(r){ #dev r = "Africa"
      selr <- rownames(private$m)==r
      ldply(unique(colnames(private$m)),function(c){ #dev r = "Africa"
        selc <- colnames(private$m)==c
        data.table(
          x=(0:1000)/1000,
          y=density(private$m[selr,selc],bw=set_bw,from=0,to=1,n=1001)$y,
          region1=c,
          region2=r
        )
      })
    }) %>% setDT

    wait("Press [ENTER] for next plot ...")
    print(ggplot(dt1,aes(x=x,y=y)) +
      geom_line() +
      xlim(set_xlims) +
      theme_classic() +
      facet_grid(region1~region2) +
      labs(x=NULL,y=NULL))



    dt2 <- ldply(unique(colnames(private$m)),function(r){
      selr <- rownames(private$m)==r
      data.table(
        x=(0:1000)/1000,
        y=density(private$m[selr,],bw=set_bw,from=0,to=1,n=1001)$y,
        region=r
      )
    }) %>% setDT
    wait("Press [ENTER] for next plot ...")
    print(ggplot(dt2,aes(x=x,y=y,colour=region)) +
      geom_line() +
      xlim(set_xlims) +
      theme_classic() +
      labs(colour="Region" , x=NULL, y=NULL))
  }
)

