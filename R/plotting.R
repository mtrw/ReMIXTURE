


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
function(run=NULL,alpha_norm_per_region=NULL,alpha_correlation=F,lon_angle=0,lat_angle=0,width_lims=c(5,35),alpha_lims=c(0.05,0.99)){
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
  if(alpha_correlation==TRUE & is.null(alpha_norm_per_region)){
    alpha_norm_per_region=FALSE
  }
  if(alpha_correlation==TRUE & alpha_norm_per_region==TRUE){
    stop("alpha_correlation and alpha_norm_per_region cannot both be TRUE")
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

  regionlist <- private$rt

  ol_dt <- as.data.table(private$results$runs[[run]]$overlap)[,region := gp_list] %>% melt(id.vars="region") %>% setnames(c("region","variable","value"),c("p1","p2","count"))
  ol_dt <- regionlist[,.(p1=region,x1=x,y1=y)][ol_dt,on=.(p1)]
  ol_dt <- regionlist[,.(p2=region,x2=x,y2=y)][ol_dt,on=.(p2)]
  ol_dt <- ol_dt[p1!=p2,]
  ol_dt[,idx:=1:.N]
  ol_dt[,count_norm := count/sum(count),by=.(p1)]

  cor_dt <- as.data.table(private$results$runs[[run]]$correlation)[,region := gp_list] %>% melt(id.vars="region") %>% setnames(c("region","variable","value"),c("p1","p2","correlation"))
  ol_dt <- cor_dt[ol_dt,on=.(p1,p2)]


  coords <- regionlist[,.(region,x,y)][gp_info[,.(region=gp)],on=.(region)]
  coords[,diversity:=private$results$runs[[run]]$diversity]
  coords[,self:=diag(private$results$runs[[run]]$overlap)]

  scale_between_f(c(private$results$runs[[run]]$overlap,private$results$runs[[run]]$diversity),width_lims[1],width_lims[2]) -> scaler_width
  scale_between_f(ol_dt[p1!=p2]$count,alpha_lims[1],alpha_lims[2]) -> scaler_alpha
  scale_between_f(ol_dt[p1!=p2]$count_norm,alpha_lims[1],alpha_lims[2]) -> scaler_alpha_norm

  cdat_div <- plyr::ldply(1:ngp, function(i) { #dev i=2
    c <- circle_seg_dt(coords$x[i], coords$y[i], coords[i,scaler_width(diversity)/2])
    c[, `:=`(region, coords$region[i])]
    c
  })
  data.table::setDT(cdat_div)
  cdat_self <- plyr::ldply(1:ngp, function(i) {
    c <- circle_seg_dt(coords$x[i], coords$y[i], coords[i,(scaler_width(diversity)/2)*(self/diversity)])
    c[, `:=`(region, coords$region[i])]
    c
  })
  data.table::setDT(cdat_self)

  cdat_div <- coords[, .(region)][cdat_div, on = "region"]
  cdat_self <- coords[, .(region)][cdat_self, on = "region"]


  #avert your eyes
  ldat_in  <- ol_dt[p1 != p2][, data.table(
    region = rep(p1, 2),
    region_tgt = rep(p2, 2),
    x = c(x1, x2),
    y = c(y1, y2),
    count = as.numeric(rep(count, 2)),
    count_norm = as.numeric(rep(count_norm, 2)),
    correlation=rep(correlation,2)
  ), by = "idx"]
  ldat <- ldply(unique(ldat_in$idx),function(i){ # dev i = 1
    ldat_in[ idx==i , {
      x1_ <<- x[1]
      x2_ <<- x[2]
      y1_ <<- y[1]
      y2_ <<- y[2]
      w_  <<- scaler_width(count[1])
      c_  <<- if(alpha_norm_per_region){alpha("#000000",scaler_alpha_norm(count_norm[1]))} else if(alpha_correlation){alpha("#000000",scaler_alpha_norm(correlation[1]))} else {alpha("#000000",scaler_alpha(count[1]))}
      a_  <<- scaler_alpha(count[1])
      r_  <<- region[1]
    }] %>% invisible
    rounded_line(x1_,y1_,x2_,y2_,w_)[,col:=c_][,region:=r_][,idx:=i][,alpha:=a_]#add gp later
  }) %>% setDT
  #\avert


  private$results$plot_data$lines <- ldat
  private$results$plot_data$circles_regiondiversity <- cdat_div
  private$results$plot_data$circles_regionunique <- cdat_self

  ce("Plotting. Raw plot data is now accessible via `$run_results`")


  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  p <- ggplot(data = world) + geom_sf(lwd = 0.05) + theme_bw() +
    theme(panel.border = element_blank()) +
    coord_sf(crs = paste0("+proj=laea +lat_0=", lat_angle, " +lon_0=", lon_angle, " +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs ")) +
    xlab(NULL) + ylab(NULL)
  pp <- p + geom_spatial_polygon(data = ldat, aes(x = x, y = y , fill=col , group=region)) +
    scale_fill_identity() +
    facet_wrap("region",ncol=4)
  ppp <- pp +
    geom_spatial_polygon(data = cdat_div, aes(x = x, y = y, group = region, fill = "#000000")) +
    geom_spatial_polygon(data = cdat_self, aes(x = x, y = y, group = region) , fill = "#FFFFFF") +
    theme(legend.position = "none")#
  ppp

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

