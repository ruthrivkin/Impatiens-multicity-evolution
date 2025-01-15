setwd("/Users/ruthrivkin/Dropbox/Grad School 2015-2020/Impatiens capensis/Data files/9CitiesSeq/Population Genetics/NewAssembly/")

ng <- theme(aspect.ratio=0.7,panel.background = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border=element_rect(fill = NA, color = "black"),
            axis.line = element_line(linewidth=1), 
            axis.ticks=element_line(color="black"),
            axis.text=element_text(color="black",size=15, margin = 0.5), 
            axis.title=element_text(color="black",size=1), 
            axis.title.y=element_text(vjust=1,face="bold",size=17),
            axis.title.x=element_text(vjust=0.1,face="bold",size=17),
            axis.text.x=element_text(size=12),
            axis.text.y=element_text(size=12),
            # legend.position = "top", legend.direction="horizontal", 
            legend.text=element_text(size=12), legend.key = element_rect(fill = "white"), 
            legend.title = element_text(size=12),legend.key.size = unit(1, "cm"))


#impute missing genotypes
dist_populations<-read.table("PCA/ld.city.mdist",header=F)
sample.info <- read.csv("../../Datasheets/24.10.21_MasterData.csv", header = T)

coord <- dplyr::select(sample.info, "Longitude", "Latitude")

# Extract environmental vars
env <- dplyr::select(sample.info,"City", "Habitat", "City.Area.km2", "PercentImpervious", "Pop.change.rate", "NDVI", "Temp", "Prec")


#memgene analysis of gen dist date
library(memgene)
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


#Plot memgene variables through space

# Set map boundary (xmin, xmax, ymin, ymax)
library(sf)
library(raster)
boundary = extent(-80.7,-79, 43, 44.1) 
boundary

# Get map outlines from rworldmap package
library(rworldmap)
library(rworldxtra)
library(maps)
library(ggsn)
library(maptools)
library(grid)
library(miscTools)
library(stringr)
library(ggpubr)
library(ggmap)

map.outline = getMap(resolution = "high")

# Crop to boundary and convert to dataframe
map.outline = raster::crop(map.outline, y = boundary) %>% fortify()

bbox <- c(left = -80.7, bottom = 43, right = 79, top = 44.1)
# Plot basemap
data.all <- cbind(env, coord,memgenedist$memgene[, 1], memgenedist$memgene[, 2], memgenedist$memgene[, 3])
head(data.all)


(memgene1 = data.all %>%
    mutate(value = ifelse(memgenedist$memgene[, 1] < 0, "negative", "positive")) %>%
    ggplot()+
    ggmap(get_stamenmap(bbox, maptype = "stamen_terrain_background", zoom = 2)) +
    geom_point(aes(x=Longitude, y=Latitude, size = abs(memgenedist$memgene[, 1]), shape = Habitat, col=value)) +
    scale_size("MEMGENE-1") +
    coord_quickmap(expand=F)+
    ggsn::north(map.outline, symbol = 10, scale = 0.05, location = "topleft")+
    xlab("Longitude (°W)")+
    ylab("Latitude (°N)") + 
    ggsn::scalebar(data = map.outline, dist = 10, dist_unit = "km", height = 0.02,
                   transform = TRUE, model = "WGS84", 
                   location = "bottomleft", anchor = c(x = -80.65, y = 43.01),
                   st.bottom = FALSE, st.size = 3.5, st.dist = 0.025)
)
ggsave("Figures/memgene1.pdf", width = 6.78, height = 5.3, dpi = 300)


(memgene2 = data.all %>%
    mutate(value = ifelse(memgenedist$memgene[, 2] < 0, "negative", "positive")) %>%
    ggplot()+
    geom_polygon(data=map.outline, aes(x=long, y=lat, group=group), fill="white",
                 colour="black")+
    geom_point(aes(x=Longitude, y=Latitude, size = abs(memgenedist$memgene[, 2]), shape = Habitat, col=value)) +
    scale_size("MEMGENE-2") +
    coord_quickmap(expand=F)+
    ggsn::north(map.outline, symbol = 10, scale = 0.05, location = "topleft")+
    xlab("Longitude (°W)")+
    ylab("Latitude (°N)")+ 
    ggsn::scalebar(data = map.outline, dist = 10, dist_unit = "km", height = 0.02,
                   transform = TRUE, model = "WGS84", 
                   location = "bottomleft", anchor = c(x = -80.65, y = 43.01),
                   st.bottom = FALSE, st.size = 3.5, st.dist = 0.025)
)
ggsave("Figures/memgene2.pdf", width = 6.78, height = 5.3, dpi = 300)

(memgene3 = data.all %>%
    mutate(value = ifelse(memgenedist$memgene[, 3] < 0, "negative", "positive")) %>%
    ggplot()+
    geom_polygon(data=map.outline, aes(x=long, y=lat, group=group), fill="white",
                 colour="black")+
    geom_point(aes(x=Longitude, y=Latitude, size = abs(memgenedist$memgene[, 3]), shape = Habitat, col=value)) +
    scale_size("MEMGENE-3") +
    coord_quickmap(expand=F)+
    ggsn::north(map.outline, symbol = 10, scale = 0.05, location = "topleft")+
    xlab("Longitude (°W)")+
    ylab("Latitude (°N)")+ 
    ggsn::scalebar(data = map.outline, dist = 10, dist_unit = "km", height = 0.02,
                   transform = TRUE, model = "WGS84", 
                   location = "bottomleft", anchor = c(x = -80.65, y = 43.01),
                   st.bottom = FALSE, st.size = 3.5, st.dist = 0.025)
)
ggsave("Figures/memgene3.pdf", width = 6.78, height = 5.3, dpi = 300)

