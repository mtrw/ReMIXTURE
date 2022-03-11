


ReMIXTURE$set( "public" , "plot_heatmap" ,
  function(colpalette=colorRampPalette(c("#f2f5ff","#214feb","#001261"))(100),...){
    if(private$runflag==FALSE){
      stop("Analysis has not been run. Perform using `$run()`")
    }
    pheatmap::pheatmap(private$overlap,cluster_rows = F,cluster_cols = F,color = colpalette,...)
  }
)


















ReMIXTURE$set( "public" , "plot_maps" ,
function(alpha_norm_per_region,lon_angle=0,lat_angle=0,width_lims=c(5,35),alpha_lims=c(0.05,0.99)){
  if(private$runflag==FALSE){
    stop("Analysis has not been run. Perform using `$run()`")
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



  ol_dt <- as.data.table(private$overlap)[,region := gp_list] %>% melt(id.vars="region") %>% setnames(c("region","variable","value"),c("p1","p2","count"))
  ol_dt <- regionlist[,.(p1=region,x1=x,y1=y)][ol_dt,on=.(p1)]
  ol_dt <- regionlist[,.(p2=region,x2=x,y2=y)][ol_dt,on=.(p2)]
  ol_dt <- ol_dt[p1!=p2,]
  ol_dt[,idx:=1:.N]


  ol_dt[,count_norm := count/sum(count),by=.(p1)]

  coords <- regionlist[,.(region,x,y)][gp_info[,.(region=gp)],on=.(region)]
  coords[,diversity:=private$diversity]
  coords[,self:=diag(private$overlap)]

  scale_between_f(c(private$overlap,private$diversity),width_lims[1],width_lims[2]) -> scaler_width
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

  ldat_in  <- ol_dt[p1 != p2][, data.table(region = rep(p1, 2),region_tgt = rep(p2, 2), x = c(x1, x2), y = c(y1, y2), count = as.numeric(rep(count, 2)), count_norm = as.numeric(rep(count_norm, 2))), by = "idx"]
  #NORMAL COLS
  #ldat_in <- coords[, .(region)][ldat_in, on = "region"] #kill

  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  p <- ggplot(data = world) + geom_sf(lwd = 0.05) + theme_bw() +
    theme(panel.border = element_blank()) +
    coord_sf(crs = paste0("+proj=laea +lat_0=", lat_angle, " +lon_0=", lon_angle, " +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs ")) +
    xlab(NULL) + ylab(NULL)

  ldat <- ldply(unique(ldat_in$idx),function(i){ # dev i = 1
    ldat_in[ idx==i , {
      x1_ <<- x[1]
      x2_ <<- x[2]
      y1_ <<- y[1]
      y2_ <<- y[2]
      w_  <<- scaler_width(count[1])
      c_  <<- alpha("#000000",scaler_alpha_norm(count_norm[1]))
      #c_  <<- alpha("#000000",scaler_alpha(count[1])) # OR un-normalised
      a_  <<- scaler_alpha(count[1])
      r_  <<- region[1]
    }] %>% invisible
    rounded_line(x1_,y1_,x2_,y2_,w_)[,col:=c_][,region:=r_][,idx:=i][,alpha:=a_]#add gp later
  }) %>% setDT

  ce("Plotting. Raw plot data is now accessible via `$results`")

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



#if(alpha_norm_per_region==TRUE)
