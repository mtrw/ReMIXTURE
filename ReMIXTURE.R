

setwd("") #set as appropriate
source("https://raw.githubusercontent.com/mtrw/tim_r_functions/master/tim_functions.R") #will also load some required packages and convenience functions.
library(ggspatial)
library(rnaturalearth)
library(pheatmap)


dm <- readRDS("MATRIX.Rds") #a matrix of similarities (e.g., IBS scores, where higher=more similar) between individuals, where the row and column names give the regions from which the sample comes. Naturally, the individuals should be in the same order along rows and columns. Also, importantly, the regions MUST be grouped together, e.g. Regions(rows) = A,A,A,C,C,D,D,D,D,D,B,B and NOT A,B,D,A,C,C,A, ...
dm[lower.tri(dm)] <- dm[upper.tri(dm)] #assure it's symmetrical

#index the positions of each region group
gplist$offset <- c(0,rle(gpcol)$lengths) %>% `[`(-length(.)) %>% cumsum %>% `+`(1)
gplist[,idx:=1:.N]

nits <- 1000 #SET: how many iterations?
sampsize <- (min(table(gpcol)) * (2/3)) %>% round #SET: how many samples per iteration (from each region)

#set up some vectors to store info later
outsize <- nits * sampsize * nrow(gplist)
select <- vector(mode="integer",length=sampsize*nrow(gplist)) #to store a list of the randomly selected samples each iteration
rawoutput <- data.table( #to store raw output each iteration
  p1 = character(length=outsize),
  p2 = character(length=outsize),
  dist = numeric(length=outsize),
  iteration = integer(length=outsize)
  
)
insert <- 1 #a flag

#run the iterations
for(iteration in 1:nits){
  #fill the `select` vector
  gplist[,{
    select[(sampsize*(idx-1)+1):((sampsize*(idx-1))+sampsize)] <<- sample(N,sampsize)-1+offset
  },by="idx"] %>% invisible
  
  #Find closest neighbours for the selected sample, store results in output table
  rnum <- 1
  #r = dm[select,select][1,]
  apply(dm[select,select],1,function(r){
    rawoutput$p1[insert] <<- colnames(dm)[select][rnum]
    rawoutput$p2[insert] <<- colnames(dm)[select][which(r==min(r))[1]]
    rawoutput$dist[insert] <<- min(r)[1]
    rawoutput$iteration[insert] <<- iteration
    rnum <<- rnum+1
    insert <<- insert+1
  }) %>% invisible
  
  ce("% complete: ",(insert/outsize)*100)
}

#summarise the output
counts <- rawoutput[ , .(count=.N) , by=.(p1,p2) ]

#look at the raw counts ("summary") matrix to check it all seems to have gone ok
dcast(counts,formula=p1~p2,value.var="count")

setorder(rawoutput,p1,p2,-dist)
rawoutput[,idx:=1:.N,by=.(p1)]

#saveRDS(counts,"COUNTS.Rds") #it's a good idea to save the output, these runs take a while
#saveRDS(rawoutput,"RAW_OUTPUT.Rds") #it's a good idea to save the output, these runs take a while


#heatmap
cnormed <- copy(counts)[,prop:=count/sum(count),by=.(p1)]
cnormed[p1!=p2][order(prop)]
cm <- as.matrix(dcast(cnormed,formula=p1~p2,value.var="prop")[,-"p1"])
dim(cm)
rownames(cm) <- colnames(cm)
hmplot <- pheatmap(cm,cluster_rows = F,cluster_cols = F)
hmplot

#Run me to save it
# dev.off()
# pdf("ReMIXTURE_heatmap.pdf",height=7,width=7,onefile = TRUE)
# print(hmplot)
# dev.off()

#significance testing by resampling from the raw output
nits <- 1000 #SET: How many samples
samplesize <- nu(rawoutput$iteration)*0.1 #SET: How many items to sample each time
nrowsit <- (nu(rawoutput$p1)**2)
nrowsout <- nrowsit*nits
#to store output
itcount <- data.table(
  p1=character(length=nrowsout),
  p2=character(length=nrowsout),
  count=numeric(length=nrowsout),
  resamp=numeric(length=nrowsout)
)

#perform resampling
for(it in 1:nits){
  #it <- 1
  ce("It: ",it)
  selectit <- sample(unique(rawoutput$iteration),samplesize)
  
  fill <- setDT(expand.grid(p1=unique(rawoutput$p1),p2=unique(rawoutput$p2)))
  insert <- rawoutput[ iteration %in% selectit , .(count=.N,resamp=it) , by=.(p1,p2) ]
  insert <- insert[fill,on=.(p1,p2)]
  insert[is.na(count),count:=0]
  insert[is.na(resamp),resamp:=it]
  
  itcount[(nrowsit*(it-1)+1):((nrowsit*(it-1))+nrow(insert))] <- insert
}

#summarise output
itcount[, pct:=(count/sum(count))*100 , by=.(resamp,p1) ]
itcount <- itcount[, .(sd_pct=sd(pct),mean_pct=mean(pct)) , by=.(p1,p2) ]
itcount[, description:=paste0( round(mean_pct-(2*sd_pct),digits=2)," (",round(mean_pct,digits=2),") ",round(mean_pct+(2*sd_pct),digits=2)  )  ]
itcount
#To save the result ...
#write.csv(itcount,"MIX_data_stdev.csv")




#Geospatial plots

#Table of lat("y")/long("x") positions for each region, and their chosen colour for the plot (recommend hex format [https://www.google.com/search?q=color+picker]). Format (example):

centres <- fread("REGIONS.csv",col.names=c("region","x","y","col"))
world <- ne_countries(scale = "medium", returnclass = "sf")
circle <- function(x=0,y=0,rad=1,n_pts=200){
  theta <- seq(from=0,to=((2*pi)-(2*pi)/n_pts),length.out=n_pts)
  data.table(
    x = y+sin(theta)*rad,
    y = x+cos(theta)*rad
  )
}

#Plot the view of the globe. Use me to optimise the angle your globe will take (see geom_sf() help for details)
p <- ggplot(data = world) +
  geom_sf(lwd=0.05) +
  coord_sf(crs = "+proj=laea +lat_0=30 +lon_0=0 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs ")
p

#create circle data
centres[,size:=counts[p1==region & p2==region]$count,by=region]
centres[,size:=size %>% scale_between(2,7)]
centres <- centres[!is.na(size)]

cdat <- ldply(1:nrow(centres),function(i){
  c <- circle(centres$x[i],centres$y[i],centres$size[i])
  c[,region:=centres$region[i]]
}) %>% setDT

#Check they plot ok. Use to adjust colours and sizes
p + geom_spatial_polygon(data=cdat,aes(x=x,y=y,fill=region,group=region),alpha=0.8) #May generate a harmless warning for leaving a stat_spatial_identity() argument default

#Normalise the nearest-neighbour counts (not strictly necessary but, for sanity purposes)
cnormed[,id:=1:.N]
cnormed <- centres[,.(p1=region,x1=x,y1=y)][cnormed,on="p1"]
cnormed <- centres[,.(p2=region,x2=x,y2=y)][cnormed,on="p2"]

#compile data for lines joining regions
ldat <- cnormed[p1 != p2][, data.table(
  region = rep(p1,2),
  x = c(x1,x2),
  y = c(y1,y2),
  count = as.numeric(rep(prop,2))
) ,by="id"]

#Do the plots
for(i in unique(ldat$region)){
  P <- p +
    geom_spatial_path(data=ldat[region==i],aes(x=y,y=x,alpha=count,size=count),colour=centres[region==i]$col,lineend="round") +
    geom_spatial_polygon(data=cdat[region==i],aes(x=x,y=y,group=region),fill=centres[region==i]$col) +
    theme(legend.position = "none")
  print(P)
}

#... or in case you'd like to save them
# dev.off()
# pdf("ReMIXTURE_geospatial.pdf",height=5,width=5,onefile = TRUE)
# for(i in unique(ldat$region)){
#   P <- p +
#     geom_spatial_path(data=ldat[region==i],aes(x=y,y=x,alpha=count,size=count),colour=centres[region==i]$col,lineend="round") +
#     geom_spatial_polygon(data=cdat[region==i],aes(x=x,y=y,group=region),fill=centres[region==i]$col) +
#     theme(legend.position = "none")
#   print(P)
# }
# dev.off()


#Plot superimposition of top two (or more) RGOs for each region
topn <- 2 #SET: How many top RGOs?
ldat_top <- ldat[,.SD[order(-count)][1:(2*topn)],by="region"]

P <- p
for(i in unique(ldat_top$region)){
  P <- P +
    geom_spatial_path(data=ldat_top[region==i],aes(x=y,y=x,alpha=count,size=count),lineend="round") +
    #geom_spatial_polygon(data=cdat[region==i],aes(x=x,y=y,group=region),fill="") +
    theme(legend.position = "none")
}
print(P)
#To plot :
# dev.off()
# pdf("ReMIXTURE_geospatial_jux.pdf",height=5,width=5,onefile = TRUE)
# print(P)
# dev.off()


# Thank you for running ReMIXTURE. If you are interested in expanding, improving, re-writing, or properly producing a package from this early implementation, please email Tim at mtrw85@gmail.com.