
ce <- function(...){   cat(paste0(...,"\n"), sep='', file=stderr()) %>% eval(envir = globalenv() ) %>% invisible() }

argGiven <- function(x){
  !is.null(x)
}
argNotGiven <- function(x){
  is.null(x)
}

clamp <- function(x,lower=min(x,na.rm=TRUE),upper=max(x,na.rm=TRUE)){
  x[x<lower] <- lower
  x[x>upper] <- upper
  x
}

diffLonRight <- function(range_lon){
  #ll=c(170,-10)
  lon1 <- wrap2lon(range_lon[1])
  lon2 <- wrap2lon(range_lon[2])
  if(length(lon1)==1){lon1 <- rep(lon1,length(lon2))}
  if(length(lon2)==1){lon2 <- rep(lon2,length(lon1))}
  sapply(1:length(lon1),function(i){
    #i=1
    d <- diff(rotateLon(c(lon1[i],lon2[i]),-lon1[i]))
    if(d<0){360+d}else{d}
  })
}

diffLat <- function(range_lat){
  range_lat[2]-range_lat[1]
}

clampNearestLon <- function(lon,range_lon,sideFlag=FALSE){
  f <- abs(diffLonRight(c(lon,range_lon[1])))<abs(diffLonRight(c(lon,range_lon[2])))
  c <- fifelse(f,range_lon[1],range_lon[2])
  if(sideFlag==TRUE){
    return(list(lon=c,sideFlag=f))
  } else {
    return(c)
  }
}

mat2dtLL <- function(m){ data.table(lon=m[1,],lat=m[2,]) }
mat2dtXY <- function(m){ data.table(x=m[1,],y=m[2,]) }
# matrix(1:4,byr=TRUE,nc=2) %T>% print %>% mat2dtXY

#dtXY2mat <- function(dt){ ... do we need it?  }

null_plot <- function(x,y,xlab=NA,ylab=NA,revx=F,revy=F,...){
  xl<-range(x,na.rm=TRUE)
  yl<-range(y,na.rm=TRUE)
  if(revx==TRUE){ xl <- rev(xl) }
  if(revy==TRUE){ yl <- rev(yl) }
  plot(NULL,xlim=xl,ylim=yl,xlab=xlab,ylab=ylab,...)
}

scale_between <- function(x,lower,upper){
  if(all(x==mean(x,na.rm=TRUE))) return(rep(mean(c(lower,upper),na.rm=TRUE),length(x)))
  ( x - min(x,na.rm=TRUE) ) / (max(x,na.rm=TRUE)-min(x,na.rm=TRUE)) * (upper-lower) + lower
}
# 1:10 %>% scale_between(-100,1)

wrapValsRange <- function(x,lower,upper){ # CANNOT be used for both lat and long, long always goes upwards as we circle, lat bounces between up wards and downwards
  #browser()
  keepSame <- (x==lower) | (x==upper)
  rangeLength <- upper-lower
  shift <- -lower
  w <- ((x+shift) %% rangeLength) - shift
  w[keepSame] <- x[keepSame]
  w
}
# plot(1:1000,1:1000 %>% wrapValsRange(-20,50))


findCentreLL <- function(range_lon,range_lat){
  clo <- if(diff(range_lon)<0){ (mean(range_lon)+180) %>% wrap2lon } else { mean(range_lon) }
  cla <- mean(range_lat)
  c(clo,cla)
}
# findCentreLL(c(20,-10),c(-90,90))

fold_matrix <- function(m){ #add lower and upper so the result is symmetrical
  m[lower.tri(m)] <- m[lower.tri(m)] + t(m)[lower.tri(m)] #lower = lower + upper
  m[upper.tri(m)] <- t(m)[upper.tri(m)] #upper = lower
  m
}
# matrix(1:4,byr=TRUE,nc=2) %T>% print %>% fold_matrix()


rev_mat_cols <- function(m){
  apply(m,1,rev) %>% t
}
# matrix(1:6,byr=TRUE,nc=3) %T>% print %>% rev_mat_j()

rev_mat_rows <- function(m){
  apply(m,2,rev)
}
# matrix(1:6,byr=TRUE,nc=3) %T>% print %>% rev_mat_i()

round2dp <- function(x){
  round(x,digits = 2)
}
# round2dp(pi)


# From BioDT
alpha <- function(colChain,setAlpha=1L){
  parseColChain(colChain,setAlpha) %>%
    t %>%
    apply( 1 , function(r) {grDevices::rgb(r[1],r[2],r[3],r[4])} )
}

# From BioDT
parseColChain <- function(colChain,setAlpha=NULL){
  # Parse and wrangle chain to rgb
  chain <- grDevices::col2rgb(colChain,a=T)/255

  # Set alpha as requested
  if(argGiven(setAlpha)){
    if(!all(setAlpha %between% 0:1)){ stop("'setAlpha' values should be scaled between 0 and 1") }
    chain[4,] <- if( length(setAlpha)==1 | length(setAlpha)==length(colChain) ){
      setAlpha
    } else {
      stats::approx((1:length(setAlpha)) %>% scale_between(1,length(colChain)),setAlpha,1:length(colChain))$y
    }
  }
  return(chain)
}


euc_dist<-function(x1,y1,x2,y2){
  sqrt( ((x1-x2)**2) + ((y1-y2)**2) )
}
# euc_dist(0,0,1,1)

circle_seg <- function(
    x=0,
    y=0,
    radius=1,
    start_radians=pi,
    end_radians=pi,
    n=200,
    ends=TRUE
){
  n <- n+1
  start_radians <- start_radians %% (2*pi)
  end_radians <- end_radians %% (2*pi)
  if(start_radians!=end_radians){ #not a full circle
    if(start_radians > end_radians){
      #n <- round( n*((2*pi)+end_radians-start_radians)/(2*pi))
      s <- seq(start_radians,(2*pi)+end_radians,length.out=n) %% (2*pi)
    } else {
      #n <- n * ((end_radians-start_radians)/(2*pi))
      s <- seq(start_radians,end_radians,length.out=n) %% (2*pi)
    }
  } else { #full circle
    s <- seq(start_radians,start_radians+(2*pi),length.out=n) %% (2*pi)
  }
  m <- matrix(
    c(x+radius*sin(s),y+radius*cos(s)),
    ncol=length(s),
    byrow=T
  )
  if(ends==T){
    m
  } else {
    m[,2:(ncol(m)-1)] # This is good
  }
}
# circle_seg(10,-5,4,-pi/2,pi) %>% mat2dtXY() %>% `[`(,plot(x,y))
# circle_seg(10,-5,3.9,-pi/2,pi,ends=F) %>% mat2dtXY() %>% `[`(,points(x,y,col="red"))


translate <- function(m,by_x=0,by_y=0){
  T_ <- matrix(
    c(
      1,0,by_x,
      0,1,by_y,
      0,0,1
    ),nrow=3,byrow=T
  )
  T_ %*% rbind(m,1) %>% `[`(1:2,)
}
# circle_seg(n=5) %>% translate(by_x = -2, by_y= 10) %>% mat2dtXY() %>% `[`(,plot(x,y))

rotate <- function(m,theta=pi){
  R_ <- matrix(
    c(
      cos(-theta) , -sin(-theta),
      sin(-theta)          , cos(-theta)
    ),nrow=2,byrow=T
  )
  R_ %*% m
}
# circle_seg(n=5) %>% rotate(pi/32) %>% mat2dtXY() %>% `[`(,plot(x,y))
# circle_seg(n=5) %>% rotate(-pi/16) %>% mat2dtXY() %>% `[`(,points(x,y,col="red"))

scale <- function(m,by_x,by_y){
  S_ <- matrix(
    c(
      by_x , 0,
      0    , by_y
    ),nrow=2,byrow=T
  )
  S_ %*% m
}
# circle_seg(n=50) %>% mat2dtXY() %>% `[`(,plot(x,y))
# circle_seg(n=50) %>% scale(by_x=.9,by_y=.2) %>% mat2dtXY() %>% `[`(,points(x,y,col="red"))

reflect <- function(m,about_x=T,about_y=F){
  S_ <- matrix(
    c(
      1 - 2*about_x , 0,
      0    , 1 - 2*about_y
    ),nrow=2,byrow=T
  )
  S_ %*% m
}
# circle_seg(.03,.03,n=50) %>% mat2dtXY() %>% `[`(,plot(x,y))
# circle_seg(.03,.03,n=50) %>% reflect(T,T) %>% mat2dtXY() %>% `[`(,points(x,y,col="red"))

angle <- function(x1,y1,x2,y2){
  x <- x2 - x1
  y <- y2 - y1
  l <- euc_dist(x1,y1,x2,y2)
  ifelse(x<0,2*pi - acos(y/l),acos(y/l))
}
# angle(0,0,1,0)/pi
# angle(0,0,-1,0)/pi
# angle(0,0,1,0)/pi
# angle(1,0,0,0)/pi

filled_line <- function(x1,y1,x2,y2,n=200,ends=T){
  m <- matrix(
    c(seq(x1,x2,length.out=n),
      seq(y1,y2,length.out=n)),
    nrow=2,byrow=T
  )
  if(ends==T){
    m
  } else {
    m[,2:(ncol(m)-1)] # This is good
  }
}
# filled_line(0,0,1,1,10) %>% mat2dtXY() %>% `[`(,plot(x,y))
# filled_line(0,0,1,.9,10,ends = F) %>% mat2dtXY() %>% `[`(,points(x,y,col="red"))

rounded_line <- function(x1,y1,x2,y2,width,n=200){
  theta <- angle(x1,y1,x2,y2)
  l <- euc_dist(x1,y1,x2,y2)
  cbind(
    circle_seg(0,0,radius=width/2,start_radians=pi*(1/2),end_radians=pi*(3/2),n=n),
    filled_line(-width/2,0,-width/2,l,n=n,ends=F),
    circle_seg(0,l,radius=width/2,start_radians=pi*(3/2),end_radians=pi/2,n=n),
    filled_line(width/2,l,width/2,0,n=n,ends=F)
  ) %>% rotate(theta) %>% translate(x1,y1)
}
# rounded_line(0,0,1,1,2) %>% mat2dtXY() %>% `[`(,plot(x,y))
# rounded_line(1,-1,0,0,1) %>% mat2dtXY() %>% `[`(,points(x,y,col="red"))


curved_rounded_line <- function(x1,y1,x2,y2,width=1,curvature=0,n=200){
  # x1<- -1
  # y1<- 2
  # x2<- 10
  # y2<- -14
  # curvature <- -pi*.6
  # n <- 30
  # width <- 2
  if(curvature==0){ return(rounded_line(x1,y1,x2,y2,width,n)) }
  if(! (curvature>-pi & curvature<pi) ){ stop("Curvatures must be in [-pi,+pi]") }
  reflectX <- sign(curvature)==1
  curvature <- abs(curvature)
  l <- euc_dist(x1,y1,x2,y2)
  x_c <- (l/2)/tan(curvature) # x coord of the circle used to generate the curved 'long' lines
  h_c <- (l/2)/sin(pi-curvature) # hypotenuse -- dist of circle centre from origin
  cbind(
    circle_seg(0,0,width/2,(pi*(1/2))-curvature,(pi*(3/2))-curvature), #el
    circle_seg(x_c,l/2,h_c+width/2,(pi*(3/2))-curvature,(pi*(3/2))+curvature), #lu
    circle_seg(0,l,width/2,(pi*(3/2))+curvature,(pi*(1/2))+curvature), #er
    circle_seg(x_c,l/2,h_c-width/2,(pi*(3/2))-curvature,(pi*(3/2))+curvature) %>% rev_mat_cols() #ll
  ) %>% reflect(about_x=reflectX,about_y=FALSE) %>% rotate(angle(x1,y1,x2,y2)) %>% translate(by_x=x1,by_y=y1)
}
# null_plot(-50:50,-50:50)
# data.table(i=seq(-pi+.0001,pi-.0001,l=200))[,{
#   curved_rounded_line(-1,2,10,-14,abs(i)/4,i) %>% mat2dtXY() %>% `[`(,polygon(x,y,col = "#99000066"))
# },by=.I]

wrap2lon <- function(x){
  wrapValsRange(x,-180,180)
}

wrap2lat <- function(x){
  isTop <- (x==90)
  X <- x + 90
  cycle <- (floor(X/180) %% 2)
  progress <- X %% 180
  w <- (((-1)**cycle)*progress + cycle*180) - 90
  w[isTop] <- 90
  w
}

ll2rad <- function(x){
  x*(360/(2*pi))
}

deg2rad <- function(x){
  x*((2*pi)/360)
}

rad2deg <- function(x){
  x*(360/(2*pi))
}

nu <-function(x){
  unique(x) %>% length
}

pd <- function(x,add=F,...){
  if(!add){
    x %>% stats::density(na.rm=TRUE,...) %>% plot(main=NA)
  } else {
    x %>% stats::density(na.rm=TRUE,...) %>% graphics::lines(main=NA)
  }
}

get_upper_tri <- function(x,...){
  x[upper.tri(x,...)]
}

#' Apply a function FUN to a vector made from the entries in the upper triangle of the matrix of m.
#'
#' @param m \[no default\] A matrix.
#' @param FUN \[no default\] The function to apply.
#' @param ... Extra arguments passed to FUN( ... )
#'
#' @return The output as described
#'
#' @export
upper_tri_ply <- function(m,FUN,...){
  FUN(m[upper.tri(m,...)],...)
}

ask <- function(q,YN = TRUE){
  if(YN){
    qn <-  paste0(q," [Y/N]: ")
  } else {
    qn <- paste0(q,": ")
  }
  ce(qn)
  r <- readline()
  if(YN){
    if(r=="Y"){
      return(TRUE)
    } else if (r=="N") {
      return(FALSE)
    } else {
      stop("Answer `Y` or `N` required.")
    }
  } else {
    return(r)
  }
}



rotateLon <- function(lon,deg,flipFlag=FALSE){
  isTop <- (lon+deg)==180 # This accounts for the case that -180 is distinct from 180 for plotting purposes. wrapValsRange will send 180 -> -180 since the wrapping is done to range [min,max)
  l <- wrapValsRange(lon+deg,-180,180)
  l[isTop] <- (lon+deg)[isTop]
  if(flipFlag==TRUE){
    return(
      list(
        wrapVals=l,
        flipped=! ( (lon+deg) %between% c(-180,180) )
      )
    )
  } else {
    return(l)
  }
}

rotateLat <- function(lat,deg){
  l <- lat+deg
  if(any(!l %between% c(-90,90))){
    #warning("After rotation, latitudes exceeded [-90,90]. These have been clamped to the range.")
    l <- clamp(l,-90,90)
  }
  l
}

rotateLatLonDtLL <- function(dtLL,lon_deg,lat_deg,splitPlotGrps=TRUE,rotatedColnames=c("lon","lat"),flipColname="lon_flipped",flipFlags=FALSE){
  # replace flip flags with a full polygon/line splitter that detects any straddling line and splits it into plotgroups, adding new points for each that sit on the line
  # lon_deg <- 140
  # lat_deg <- 0
  # flipColname="lon_flipped"
  # rotatedColnames=c("lon","lat")
  dt <- copy(dtLL)

  rLo <- rotateLon(dt$lon,lon_deg,flipFlag = TRUE)
  rLa <- rotateLat(dt$lat,lat_deg)

  dt[,c(rotatedColnames) := .( rLo$wrapVals , rLa )]

  if(splitPlotGrps==TRUE){
    if(is.null(dt$plotGrp)){ dt$plotGrp <- 1L }
    dt[,c(flipColname) := .( rLo$flipped )]
    dt[, nextLon :=shift(lon,-1), by=.(plotGrp) ]
    dt[, nextLat :=shift(lat,-1), by=.(plotGrp) ]
    dt[,crossLon:=lon_flipped!=shift(lon_flipped,-1),by=.(plotGrp)]
    dt[is.na(crossLon),crossLon:=FALSE]
    dt <- dt[,{
      if(crossLon==FALSE){
        .SD
      } else {
        replacer <- .SD[rep(1,3)]
        clampLon1 <- sign(lon)*180
        clampLon2 <- -clampLon1
        intLat <- stats::approx(rotateLon(c(lon,nextLon),-lon),c(lat,nextLat),clampLon1-lon)$y
        replacer[2:3,]$lon <- c(clampLon1,clampLon2)
        replacer[2:3,]$lat <- intLat
        replacer[3,]$lon_flipped <- replacer[3,!lon_flipped]
        replacer
      }
    },by=.I]

    dt[,plotGrp:=.GRP,by=.(plotGrp,lon_flipped)]
    dt[,ptIdx:=1:.N,by=.(plotGrp)]
    dt[,I:=NULL][,nextLon:=NULL][,nextLat:=NULL][,crossLon:=NULL][,c(flipColname):=NULL][]


  }

  if (flipFlags==TRUE) {
    dt[,c(flipColname) := .( rLo$flipped )]
  }

  dt[]
}

fillPathDtLL <- function(dtLL,maxGap=0.5){
  dt <- copy(dtLL)
  dt[,nextLat:=shift(lat,-1,type="cyclic"),by=.(plotGrp)]
  dt[,nextLon:=shift(lon,-1,type="cyclic"),by=.(plotGrp)]
  dt[,toNextLat:=abs(nextLat-lat),by=.(plotGrp)]
  dt[,toNextLon:=abs(nextLon-lon),by=.(plotGrp)]
  dt[,fillNext:=(toNextLat>maxGap) | (toNextLon>maxGap)]
  dt[is.na(fillNext),fillNext:=FALSE]
  dt <- dt[,{
    if(fillNext==TRUE){
      nFill <- floor(max(toNextLat,toNextLon)/maxGap)
      replacer <- .SD[rep(1,nFill+1),]
      replacer$lon <- seq(lon,nextLon,length.out=nFill+1)
      replacer$lat <- seq(lat,nextLat,length.out=nFill+1)
      replacer
    } else {
      .SD
    }
  },by=.I][,I:=NULL][,nextLat:=NULL][,nextLon:=NULL][,toNextLat:=NULL][,toNextLon:=NULL][,fillNext:=NULL][]
  dt
}

fillPathDtXY <- function(dtXY,maxGap=0.5){
  dt <- copy(dtXY)
  dt[,nextY:=shift(Y,-1,type="cyclic"),by=.(plotGrp)]
  dt[,nextX:=shift(X,-1,type="cyclic"),by=.(plotGrp)]
  dt[,toNextY:=abs(nextY-Y),by=.(plotGrp)]
  dt[,toNextX:=abs(nextX-X),by=.(plotGrp)]
  dt[,fillNext:=(toNextY>maxGap) | (toNextX>maxGap)]
  dt[is.na(fillNext),fillNext:=FALSE]
  dt <- dt[,{
    if(fillNext==TRUE){
      nFill <- floor(max(toNextY,toNextX)/maxGap)
      replacer <- .SD[rep(1,nFill+1),]
      replacer$X <- seq(X,nextX,length.out=nFill+1)
      replacer$Y <- seq(Y,nextY,length.out=nFill+1)
      replacer
    } else {
      .SD
    }
  },by=.I][,I:=NULL][,nextY:=NULL][,nextX:=NULL][,toNextY:=NULL][,toNextX:=NULL][,fillNext:=NULL][]
  dt
}

sinc_un <- function(x){ fifelse(x==0,1,sin(x)/x) }


makeBorder <- function(range_lon,range_lat,n=500,clippingMask=FALSE){

  llo <- diffLonRight(range_lon)
  lla <- diffLat(range_lat)
  left   <- -llo/2
  right  <- llo/2
  bottom <- -lla/2
  top    <- lla/2

  if( clippingMask==FALSE ){

    out <- data.table(
      lon = c( seq(left,by=llo/(n-1),length.out=n)   , rep(right,n)                           , seq(right,by=-llo/(n-1),length.out=n) , rep(left,n) ) ,
      lat = c( rep(bottom,n)                         , seq(bottom,by=lla/(n-1),length.out=n)  , rep(top,n)                            , seq(top,by=-lla/(n-1),length.out=n) )
    )[,plotGrp:=1L][]

  } else {

    out <- rbind(
      data.table( #bottom from right
        lon=c( seq(left,right,length.out=n) , 180 , -180 , left ),
        lat=c( rep(bottom,n) , -90 , -90 , bottom ),
        plotGrp=1,
        TB=c(rep(NA_character_,n),"B","B",NA_character_),
        LR=c(rep(NA_character_,n),"R","L",NA_character_)
      ),
      data.table( # top from right
        lon=c( seq(left,right,length.out=n) , 180 , -180 , left ),
        lat=c( rep(top,n) , 90 , 90 , top ),
        plotGrp=2,
        TB=c(rep(NA_character_,n),"T","T",NA_character_),
        LR=c(rep(NA_character_,n),"R","L",NA_character_)
      ),
      data.table( # left from bottom
        lon=c( rep(left,n) , -180 , -180 , left ),
        lat=c( seq(bottom,top,length.out=n) , 90 , -90 , bottom ),
        plotGrp=3,
        TB=c(rep(NA_character_,n),"T","B",NA_character_),
        LR=c(rep(NA_character_,n),"L","L",NA_character_)
      ),
      data.table( #right from bottom
        lon=c( rep(right,n) , 180 , 180 , right ),
        lat=c( seq(bottom,top,length.out=n) , 90 , -90 , bottom ),
        plotGrp=4,
        TB=c(rep(NA_character_,n),"T","B",NA_character_),
        LR=c(rep(NA_character_,n),"R","R",NA_character_)
      )
    )
  }
  out[,ptIdx:=1:.N,by=.(plotGrp)]
  out[]
}

makeMapDataLatLonLines <- function(range_lon=c(-180,180),range_lat=c(-90,90),n=100){
  plotGrp <- 0L
  dt <- c(
    lapply(seq(-90,90,by=10),function(la){
      plotGrp <<- plotGrp+1
      data.table(
        lon=seq(-180,180,l=n),
        lat=la,
        plotGrp=plotGrp
      )
    }),
    lapply(seq(-180,180,by=10),function(lo){
      plotGrp <<- plotGrp+1
      data.table(
        lon=lo,
        lat=seq(-90,90,l=round(n/2)),
        plotGrp=plotGrp
      )
    })
  ) %>% rbindlist()
  dt[(lon %betweenLon% range_lon) & (lat %betweenLat% range_lat)][]
}

#' Winkel III map projection
#'
#' @description
#' One of ReMIXTURE's three inbuilt map projection functions, the Winkel Tripel aka Winkel III. See (https://en.wikipedia.org/wiki/Winkel_tripel_projection)
#'
#' @param dtLL A data.table with at minimum two numerical columns, one named 'lat' and the other 'lon', containing latitude and longitude points (in degrees) to be transformed.
#' @param projColNames A length-2 character vector. In the returned data.table, what should the transformed columns be named?
#'
#' @export
winkelIII <- function(dtLL,projColNames=c("x_W3","y_W3")){
  dt <- copy(dtLL)
  radLon <- dt$lon %>% deg2rad()
  radLat <- dt$lat %>% deg2rad()
  phi <- acos(2/pi)
  alpha <- acos(cos(radLat) * cos(radLon/2))
  lon_W3 <- ((1/2) * ((radLon*cos(phi)) + ((2*cos(radLat)*sin((radLon/2)))/(sinc_un(alpha))))) %>% rad2deg()
  lat_W3 <- ((1/2) * ( radLat + (sin(radLat)/sinc_un(alpha)) )) %>% rad2deg()
  dt[,  c(projColNames):=.( lon_W3 , lat_W3 )  ][]
}

#' Eckert IV Map Projection
#'
#' @description
#' One of ReMIXTURE's three inbuilt map projection functions, the Eckert IV. See (https://en.wikipedia.org/wiki/Eckert_IV_projection)
#'
#' @param dtLL A data.table with at minimum two numerical columns, one named 'lat' and the other 'lon', containing latitude and longitude points (in degrees) to be transformed.
#' @param precision The precision with which the projected values are calculated (the calculation is numerical using Newton's Method).
#' @param projColNames A length-2 character vector. In the returned data.table, what should the transformed columns be named?
#'
#' @export
eckertIV <- function(dtLL,precision=0.001,projColNames=c("x_EIV","y_EIV")){
  dt <- copy(dtLL)
  radLon <- dt$lon %>% deg2rad()
  radLat <- dt$lat %>% deg2rad()
  th  <- radLat/2
  dTh <- rep(1.0+precision,length(radLat))
  while(any(dTh>precision)){
    dTh <- -( th + sin(th)*cos(th) + 2*sin(th) - (2+(pi/2))*sin(radLat) ) / ( 2*cos(th)*(1+cos(th)) )
    th <- th + dTh
  }
  lon_W3 <- (2 / ((pi*(4+pi))**(1/2)) ) * radLon * (1+cos(th))
  lat_W3 <- (2 * ((pi*(4+pi))**(1/2)) ) * sin(th)
  dt[,  c(projColNames):=.( lon_W3 , lat_W3 )  ][]
}

#' Equirectangular Map Projection
#'
#' @description
#' One of ReMIXTURE's three inbuilt map projection functions, the simple equirectangular projection. See (https://en.wikipedia.org/wiki/Equirectangular_projection)
#'
#' @param dtLL A data.table with at minimum two numerical columns, one named 'lat' and the other 'lon', containing latitude and longitude points (in degrees) to be transformed.
#' @param projColNames A length-2 character vector. In the returned data.table, what should the transformed columns be named?
#'
#' @export
equirectangular <- function(dtLL,projColNames=c("x_ER","y_ER")){
  copy(dtLL)[,c(projColNames):=.(lon,lat)][]
}


emptyPlot <- function(rangeX,rangeY,...){
  plot(
    x=NULL,
    y=NULL,
    xlab=NA,
    ylab=NA,
    xlim=rangeX,
    ylim=rangeY,
    axes=0,
    ...
  )
}

plotEmptyMap <- function( range_lon=c(-180,180), range_lat=c(-90,90), projFun=equirectangular, ... ){
  p <- makeBorder(range_lon,range_lat) %>%
    projFun(projColNames=c("x","y"))
  emptyPlot(range(p$x),range(p$y),...)
  return(data.table(x=range(p$x),y=range(p$y)) %>% invisible)
}



`%betweenLon%` <- function(lon,range_lon){
  if((range_lon[2]-range_lon[1])<0){
    (lon>=range_lon[1]) | (lon<=range_lon[2])
  } else {
    lon %between% range_lon
  }
}

`%betweenLat%` <- function(lat,range_lat){
  lat %between% range_lat
}

truncateLatLonDtLL <- function(dtLL,range_lon,range_lat,updatePlotgrp=TRUE,plotGrpClamp){
  dtLL <- copy(dtLL)
  # dtLL <- copy(mapData110)
  # range_lon = c(50,-50)
  # range_lat = c(-90,-75)

  dtLL[ ,inBox:=lon %betweenLon% range_lon & lat %between% range_lat ]

  if(plotGrpClamp==TRUE){
    dtLL[,clampGrps:=nu(inBox)==2,by=.(plotGrp)]
    dtLL[clampGrps==TRUE & !lon %betweenLon% range_lon,lon:=clampNearestLon(lon,range_lon)]
    dtLL[clampGrps==TRUE & !lat %betweenLat% range_lat,lat:=clampNearestLat(lat,range_lat)]
    dtLL[inBox==FALSE,outRunIdx:=1:.N,by=.(rleid(inBox))]
    dtLL[,inPlot:=clampGrps==TRUE | inBox==TRUE, by=.()]
    dtLL[,clampGrps:=NULL]
  }

  if(updatePlotgrp==TRUE){
    if(is.null(dtLL$plotGrp)){ dtLL$plotGrp <- NA }
    dtLL[,plotGrp:=rleid(plotGrp,inBox)]
    dtLL[inBox==TRUE,plotGrp:=rleid(plotGrp)]
  }

  dtLL[inBox==TRUE,][,inBox:=NULL][]
}


plotMapItem <- function(dtLL,range_lon=NULL,range_lat=NULL,projFun=equirectangular,plotFun=graphics::lines,splitPlotGrps=TRUE,...){
  centre_lonLat <- if(is.null(range_lon) & is.null(range_lat)){
    c(0.0,0.0)
  } else {
    findCentreLL(range_lon,range_lat)
  }
  dt <- copy(dtLL) %>% rotateLatLonDtLL( -centre_lonLat[1], -centre_lonLat[2] , splitPlotGrps ) %>% projFun(projColNames=c("x","y"))

  for( pg in unique(dt$plotGrp) ){
    plotFun( dt[plotGrp==pg]$x, dt[plotGrp==pg]$y, ... )
    #polygon( dt[plotGrp==pg]$x, dt[plotGrp==pg]$y, col="green")
  }
}

plotMapBorder <- function(range_lon,range_lat,projFun=equirectangular,type=c("both","border","mask"),plotEdges=NULL,extraEdge=40L,maskCol="#FFFFFF",...){
  dt1 <- dt2 <- NULL
  if(type[1]=="mask" | type[1]=="both"){
    dt2 <- makeBorder(range_lon,range_lat,clippingMask = T) %>% projFun(projColNames=c("x","y"))
    if(!is.null(plotEdges)){
      dt2[LR=="L",   x:=plotEdges$x[1]-extraEdge ]
      dt2[LR=="R",  x:=plotEdges$x[2]+extraEdge ]
      dt2[TB=="B", y:=plotEdges$y[1]-extraEdge ]
      dt2[TB=="T",    y:=plotEdges$y[2]+extraEdge ]
    }
    dt2 %<>% fillPathDtLL
    for( pg in unique(dt2$plotGrp) ){
      graphics::polygon( dt2[plotGrp==pg]$x, dt2[plotGrp==pg]$y, col=maskCol , border = NA )
    }
  }

  if(type[1]=="border" | type[1]=="both"){
    dt1 <- makeBorder(range_lon,range_lat) %>% fillPathDtLL %>% projFun(projColNames=c("x","y"))
    dt1[,ptIdx:=1:.N,by=.(plotGrp)]
    for( pg in unique(dt1$plotGrp) ){
      graphics::lines( dt1[plotGrp==pg]$x, dt1[plotGrp==pg]$y, ... )
    }
  }

  list(dt1,dt2) %>% invisible
}

drawUnitSquareTopRight <- function(x,y,border="#000000",col="#AA0000",xAdj=0.0,yAdj=0.0,scale=1.0,...){
  graphics::rect(x+xAdj,y+yAdj,x+xAdj+scale,y+yAdj+scale,border=border,col=col,...)
}

drawUnitCircleTopRight <- function(x,y,border="#000000",col="#AA0000",xAdj=0.0,yAdj=0.0,scale=1.0,n = 180,...){
  c <- circle_seg(x+0.5+xAdj,y+0.5+yAdj,radius=scale/2,start_radians = 0.0,end_radians = 0.0)
  graphics::polygon(c[1,],c[2,],border=border,col=col,...)
}
