#' @importFrom data.table data.table
NULL
#' @importFrom magrittr `%>%`
NULL
#' @import R6
NULL

ce <- function(...){   cat(paste0(...,"\n"), sep='', file=stderr()) %>% eval(envir = globalenv() ) %>% invisible() }
nu <-function(x){
  unique(x) %>% length
}
scale_between <- function(x,lower,upper){
  if(all(x==mean(x,na.rm=T))) return(rep(mean(c(lower,upper),na.rm=T),length(x)))
  ( x - min(x,na.rm=T) ) / (max(x,na.rm=T)-min(x,na.rm=T)) * (upper-lower) + lower
}

replace_levels_with_colours <- function(x,palette="Berlin",alpha=1,fun="diverge_hcl",plot=FALSE,newplot=TRUE){
  #require(colorspace)
  n <- nu(x[!is.na(x)])
  cols <- match.fun(fun)(n,palette = palette,alpha = alpha)
  colvec <- swap( x , unique(x[!is.na(x)]) , cols , na.replacement = NA )
  if(plot==FALSE) {
    return(colvec)
  } else {
    # null_plot(y=1:length(cols),x=rep(1,length(cols)),xaxt="n",yaxt="n")
    # text(y=1:length(cols),x=rep(1,length(cols)),labels=unique(x),col=cols)
    if(newplot) {null_plot(x=0,y=0,xaxt="n",yaxt="n",bty="n")}
    legend(x="topleft",legend=unique(x[!is.na(x)]),fill=cols,text.col=cols)
  }
}
swap <- function(vec,matches,names,na.replacement=NA){
  orig_vec <- vec
  #if(sum(! matches %in% names ) > 0 ) { stop("Couldn't find all matches in names") }
  if(length(matches) != length(names)) { stop("Lengths of `matches` and `names` vectors don't match, you old bison!") }
  if(is.factor(vec)) { levels(vec) <- c(levels(vec),names,na.replacement) }
  vec[is.na(orig_vec)] <- na.replacement
  plyr::l_ply( 1:length(matches) , function(n){
    vec[orig_vec==matches[n]] <<- names[n]
  })
  vec
}
null_plot <- function(x,y,xlab=NA,ylab=NA,...){
  plot(NULL,xlim=range(x,na.rm=T),ylim=range(y,na.rm=T),xlab=xlab,ylab=ylab,...)
}

fill_upper_from_lower <- function(M){
  M[upper.tri(M)] <- t(M)[upper.tri(M)]
  M
}

fill_lower_from_upper <- function(M){
  M[lower.tri(M)] <- t(M)[lower.tri(M)]
  M
}

circle <- function(x=0,y=0,rad=1,n_pts=200){
  theta <- seq(from=0,to=((2*pi)-(2*pi)/n_pts),length.out=n_pts)
  data.table(
    x = y+sin(theta)*rad,
    y = x+cos(theta)*rad
  )
}
