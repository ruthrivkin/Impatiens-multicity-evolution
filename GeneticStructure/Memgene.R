setwd("/Users/ruthrivkin/Dropbox/Grad School 2015-2020/Impatiens capensis/Data files/9CitiesSeq/Population Genetics/NewAssembly/")


#impute missing genotypes
dist_populations<-read.table("ld.city.mdist",header=F)

#memgene analysis of gen dist date

memgenedist <- mgQuick(dist_populations, coord, longlat = TRUE, truncation = NULL,
        transformation = NULL, forwardPerm = 100, forwardAlpha = 0.05,
        finalPerm = NULL, doPlot = NULL, verbose = TRUE)

## Find the proportional variation explained by each MEMGENE variable
MEMGENEProp <- memgenedist$sdev/sum(memgenedist$sdev)
format(signif(MEMGENEProp, 3)[1:3], scientific=FALSE)

pdf("memGENE_axis1.pdf", width = 6.78, height = 5.3)
plot(coord, type="n",  xlab = "Longitude (°W)", ylab = "Latitude (°N)", axes=T)
mgMap(coord, memgenedist$memgene[, 1], add.plot=TRUE,
      legend=F)
dev.off()

pdf("memGENE_axis2.pdf", width = 6.78, height = 5.3)
plot(coord, type="n",  xlab = "Longitude (°W)", ylab = "Latitude (°N)", axes=T)
mgMap(coord, memgenedist$memgene[, 2], add.plot=TRUE,
      legend=F)
dev.off()

pdf("memGENE_axis3.pdf", width = 6.78, height = 5.3)
plot(coord, type="n",  xlab = "Longitude (°W)", ylab = "Latitude (°N)", axes=T)
mgMap(coord, memgenedist$memgene[, 3], add.plot=TRUE,
      legend=F)
dev.off()

memgenedist$RsqAdj # =  0.0287916
#Note that this low value should be interpreted not as an inadequacy of the regression to explain variation, but rather that there is only a small proportion of all genetic variation that can be attributed to spatial patterns

## Find the proportional variation explained by each MEMGENE variable
memgeneProp <- memgenedist$sdev/sum(memgenedist$sdev)

## Neatly print proportions for the first three MEMGENE variables
format(signif(memgeneProp, 3)[1:3], scientific=FALSE)

#MEMGENE1 MEMGENE2 MEMGENE3 
#""0.381"  "0.228"  "0.215" 


