ce <- function(...){   cat(paste0(...,"\n"), sep='', file=stderr()) %>% eval(envir = globalenv() ) %>% invisible() }

wait <- function(message="Press [enter] to continue"){
  invisible(readline(prompt=message))
}

null_plot <- function(x,y,xlab=NA,ylab=NA,revx=F,revy=F,...){
  xl<-range(x,na.rm=T)
  yl<-range(y,na.rm=T)
  if(revx==T){ xl <- rev(xl) }
  if(revy==T){ yl <- rev(yl) }
  plot(NULL,xlim=xl,ylim=yl,xlab=xlab,ylab=ylab,...)
}

scale_between <- function(x,lower,upper){
  if(all(x==mean(x,na.rm=T))) return(rep(mean(c(lower,upper),na.rm=T),length(x)))
  ( x - min(x,na.rm=T) ) / (max(x,na.rm=T)-min(x,na.rm=T)) * (upper-lower) + lower
}


fold_matrix <- function(m){ #add lower and upper so the result is symmetrical
  m[lower.tri(m)] <- m[lower.tri(m)] + t(m)[lower.tri(m)] #lower = lower + upper
  m[upper.tri(m)] <- t(m)[upper.tri(m)] #upper = lower
  m
}

round2dp <- function(x){
  round(x,digits = 2)
}


euc_dist<-function(x1,x2,y1,y2){
  sqrt( ((x1-x2)**2) + ((y1-y2)**2) )
}


circle_seg <- function(
  x=0,
  y=0,
  radius=1,
  start_radians=pi,
  end_radians=pi,
  n=200
){
  start_radians <- start_radians %% (2*pi)
  end_radians <- end_radians %% (2*pi)
  if(start_radians!=end_radians){ #not a full circle
    if(start_radians > end_radians){
      n <- round( n*((2*pi)+end_radians-start_radians)/(2*pi))
      s <- seq(start_radians,(2*pi)+end_radians,length.out=n) %% (2*pi)
    } else {
      n <- n * ((end_radians-start_radians)/(2*pi))
      s <- seq(start_radians,end_radians,length.out=n) %% (2*pi)
    }
  } else { #full circle
    n <- round(n)+1
    s <- seq(start_radians,start_radians+(2*pi),length.out=n) %% (2*pi)
  }
  matrix(
    c(x+radius*sin(s),y+radius*cos(s)),
    ncol=length(s),
    byrow=T
  )
}

circle_seg_dt <- function(x=0,y=0,radius=1,start_radians=pi,end_radians=pi,n=200){
  d <- circle_seg(x=x,y=y,radius=radius,start_radians=start_radians,end_radians=end_radians,n=n)
  data.table(
    x=d[1,],
    y=d[2,]
  )
}

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

rotate <- function(m,theta=pi){
  R_ <- matrix(
    c(
      cos(-theta) , -sin(-theta),
      sin(-theta)          , cos(-theta)
    ),nrow=2,byrow=T
  )
  R_ %*% m
}

scale <- function(m,by_x,by_y){
  S_ <- matrix(
    c(
      by_x , 0,
      0    , by_y
    ),nrow=2,byrow=T
  )
  S_ %*% m
}

reflect <- function(m,about_x=T,about_y=F){
  S_ <- matrix(
    c(
      1 - 2*about_x , 0,
      0    , 1 - 2*about_y
    ),nrow=2,byrow=T
  )
  S_ %*% m
}

angle <- function(x1,y1,x2,y2){
  x <- x2 - x1
  y <- y2 - y1
  l <- euc_dist(x1,x2,y1,y2)
  ifelse(x<=0,2*pi - acos(y/l),acos(y/l))
}

filled_line <- function(x1,y1,x2,y2,n=200,ends=T){
  m <- matrix(
    c(seq(x1,x2,length.out=n),
      seq(y1,y2,length.out=n)),
    nrow=2,byrow=T
  )
  if(ends==T){
    m
  } else {
    m[,2:(ncol(m)-1)]
  }
}

rounded_line <- function(x1,y1,x2,y2,width,return_mat=F,n=200){
  theta <- angle(x1,y1,x2,y2)
  end1 <- circle_seg(0,0,radius=width/2,start_radians=pi*(1/2),end_radians=pi*(3/2),n=n)
  end2 <- end1 %>% rotate(pi) %>% translate(by_y = euc_dist(x1,x2,y1,y2))
  mid1 <- filled_line(last(t(end1))[,1],last(t(end1))[,2],first(t(end2))[,1] , first(t(end2))[,2], n=n, ends = F )
  mid2 <- filled_line(last(t(end2))[,1],last(t(end2))[,2],first(t(end1))[,1] , first(t(end1))[,2], n=n, ends = F )
  #ce("RLine, theta is: ",theta %>% round2dp)
  shape <- cbind(
    end1,
    mid1,
    end2,
    mid2
  ) %>%
    rotate(theta) %>%
    translate(x1,y1)

  if(return_mat==TRUE){
    shape
  } else {
    data.table(
      x=shape[1,],
      y=shape[2,]
    )
  }

}

scale_between_f <- function(x,lower,upper){
  function(y){
    if(all(y==mean(x,na.rm=T))) return(rep(mean(c(lower,upper),na.rm=T),length(y)))
    ( y - min(x,na.rm=T) ) / (max(x,na.rm=T)-min(x,na.rm=T)) * (upper-lower) + lower
  }
}

deg2rad <- function(x){
  x*(360/(2*pi))
}

rad2deg <- function(x){
  (x/360)*(2*pi)
}

nu <-function(x){
  unique(x) %>% length
}

pd <- function(x,add=F,...){
  if(!add){
    x %>% density(na.rm=TRUE,...) %>% plot(main=NA)
  } else {
    x %>% density(na.rm=TRUE,...) %>% lines(main=NA)
  }
}

get_upper_tri <- function(x,...){
  x[upper.tri(x,...)]
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



