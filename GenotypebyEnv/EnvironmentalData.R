
setwd("/Users/ruthrivkin/Dropbox/Grad School 2015-2020/Impatiens capensis/Data files/9CitiesSeq/EnvironmentData/")

#Start by extracting urban data -- impervious surface and built area
# loads the packages used in this guide

library(raster)
library(ggplot2)
library(dplyr)
library(tidyverse)

#load layer

#SEDAC: Brown de Colstoun, E. C., C. Huang, P. Wang, J. C. Tilton, B. Tan, J. Phillips, S. Niemczura, P.-Y. Ling, and R. E. Wolfe. 2017. Global Man-made Impervious Surface (GMIS) Dataset From Landsat. Palisades, New York: NASA Socioeconomic Data and Applications Center (SEDAC). https://doi.org/10.7927/H4P55KKF. Accessed DAY MONTH YEAR.
#https://sedac.ciesin.columbia.edu/data/set/ulandsat-gmis-v1
#Global Man-made Impervious Surface (GMIS) Dataset V1 from 2010
impsurf.rast <- raster("Rasters/CAN_gmis_impervious_surface_percentage_geographic_30m/CAN_gmis_impervious_surface_percentage_geographic_30m.tif")
names(impsurf.rast) <- "PercentImpervious"

#Get worldclim data for temp and preci
wclim <- getData("worldclim",var="bio", res=0.5, lon=-80, lat=44)

wclim <- wclim[[c(1,12)]]
names(wclim) <- c("Temp","Prec")

#Get modis NDVI data for May-Sept 2019 (sampling period)
ndvi.may <- raster("Rasters/MODIS_NDVI/MOD13A3.061__1_km_monthly_NDVI_doy2019121_aid0001.tif")
ndvi.jun <- raster("Rasters/MODIS_NDVI/MOD13A3.061__1_km_monthly_NDVI_doy2019152_aid0001.tif")
ndvi.jul <- raster("Rasters/MODIS_NDVI/MOD13A3.061__1_km_monthly_NDVI_doy2019182_aid0001.tif")
ndvi.aug <- raster("Rasters/MODIS_NDVI/MOD13A3.061__1_km_monthly_NDVI_doy2019213_aid0001.tif")
ndvi.sept <- raster("Rasters/MODIS_NDVI/MOD13A3.061__1_km_monthly_NDVI_doy2019244_aid0001.tif")

ndvi.monthly <- stack(ndvi.may,ndvi.jun,ndvi.jul,ndvi.aug,ndvi.sept)
plot(ndvi.monthly)
ndvi <- mean(ndvi.monthly) #average for annual growing season value
ndvi <- ndvi*0.0001 #multiply by conversion factor

ndvi.latlon <- projectRaster(ndvi, crs='+proj=longlat +datum=WGS84 +no_def')

names(ndvi.latlon) <- c("NDVI")
plot(ndvi.latlon)

#Stack layers


#Change extent to match study area
extent <- c(-81,-79, 43, 44.5) 

impsurf.rast1 <- crop(impsurf.rast, extent)
ndvi1 <- crop(ndvi.latlon, extent)
wclim1 <- crop(wclim, extent)

#Change values of impsurf raster to set 250 and 255 to 0

impsurf.rast1[impsurf.rast1 == 255] <- 0 #Change 255 values to 0
impsurf.rast1[impsurf.rast1 == 200] <- 0 #Change 200 values to 0

plot(impsurf.rast1)
plot(ndvi1)
plot(wclim1)

#Merge layers
library(terra)

ndvi2<- projectRaster(ndvi1, impsurf.rast1) #project the same resolution
wclim2 <-projectRaster(wclim1, impsurf.rast1)

#convert farenheight to celcies
temp <- ((wclim2$Temp)/(10)) #average for annual growing season value, need to fix in all.raster because incorrect
prec <- wclim2$Prec
names(ndvi2) <- c("NDVI")

plot(temp)

complete.raster <- stack(impsurf.rast1,ndvi2,temp, prec)
nlayers(complete.raster)
plot(complete.raster)

writeRaster(complete.raster, "Rasters/isa.ndvi.wclim.tif", overwrite=TRUE)


#4. Extract data from raster for current conditions
# Load site locations (needs coordinates for each sampling/target site)
sites <- read.csv("../Datasheets/24.02.26_SeedSamples.csv")


# set x and y for lat and long
x <- sites$Longitude
y <- sites$Latitude

# set up coordinates as a tibble
pts <- cbind(x, y) %>% as_tibble()
pts <- SpatialPoints(pts, proj4string = CRS("+proj=longlat +datum=WGS84"))
pts


#Check points CRS matches raster CRS

projection(pts) == projection(complete.raster)

#Create a tibble or data.frame to store data for each point.
urban.data = tibble(ID = 1:nrow(pts@coords),
                    Lon = pts$x,
                    Lat = pts$y
)
urban.data



#Create raster list with  buffer around each sample
store_data = list()
for (i in 1:nlayers(complete.raster)){
  store_data[[i]] = raster::extract(complete.raster[[i]], pts, buffer=20, fun = mean)
}

names(store_data) = names(complete.raster)
all.data = bind_cols(urban.data, as_tibble(store_data))
all.data

# Check each column for NA values
na.check = map_int(all.data, ~sum(is.na(.)))
summary(na.check > 0) #FALSE so should be good?

# Remove NA records
all.data.nona = all.data %>% drop_na

# Update column values to be biologically meaningful
#ImpSurface: 0-100 percent; non-HBASE 200 = 0 impervious surface (aka no urban built area)
str(all.data.nona)
(all.data.final <- all.data.nona %>%
    rename(Longitude = Lon, Latitude = Lat) %>%
    mutate(PercentImpervious = as.numeric(PercentImpervious) / 100)
)

str(all.data.final)

#Export data to a csv file and individual id coordinates
head(sites)
sites <- rename(sites, "ID" = "X")

urbandata <- merge(sites, all.data.final, by = "ID")
head(urbandata)
str(urbandata)
#write to csv
write_csv(urbandata, "24.03.25_FinalData.csv") 



#Environmental differences between pops
library(terra)
library(maps)
library(RStoolbox)

sites <- read.csv("../Datasheets/24.03.25_FinalData.csv")

rast <- rast("Rasters/isa.ndvi.wclim.tif")
names(rast) <- c("ISA", "NDVI", "Temp", "Prec")
plot(rast)

#PCA
pcamap<-rasterPCA(rast,spca=TRUE)

knitr::kable(round(pcamap$model$loadings[,1:3],3))
#|     | Comp.1| Comp.2| Comp.3|
#  |:----|------:|------:|------:|
#  |ISA  |  0.352|  0.933|  0.052|
#  |NDVI | -0.509|  0.136|  0.840|
#  |Temp |  0.558| -0.274|  0.270|
#  |Prec | -0.553|  0.191| -0.467|

summary(pcamap$model) # Prop var: PC1- 0.5680118 PC2-0.2058721 PC3-0.1385219 PC-40.08759421

c10 <- c(
  "dodgerblue2", "maroon", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6")# lt purple

pdf("../Figures/PCA1.Map.Environment.pdf")
plot(pcamap$map,1)
points(sites$Longitude, sites$Latitude, col = c10[factor(sites$City)])
legend("topleft",
       legend = levels(factor(sites$City)),
       pch = 19,
       col = c10)
dev.off()

pdf("../Figures/PCA2.Map.Environment.pdf")
plot(pcamap$map,2)
points(sites$Longitude, sites$Latitude, col = c10[factor(sites$City)])
legend("topleft",
       legend = levels(factor(sites$City)),
       pch = 19,
       col = c10)
dev.off()

#Look at differences between sites

cor.test(log(sites$City.Area.km2), sites$PercentImpervious) #-0.2407659 
cor.test(log(sites$City.Area.km2), sites$NDVI) #0.1917156
cor.test(log(sites$City.Area.km2), sites$Temp) #0.3702436
cor.test(log(sites$City.Area.km2), sites$Prec) #-0.5374937 




library(car)
(summary <- sites %>%
    group_by(PopName, City, Habitat, Latitude, Longitude,
             City.Census.2021,City.Area.km2,Pop.change.rate, City.Pca, 
             PercentImpervious, BuiltArea, LandUseType, Prec, Temp, NDVI,
             Pi, FIS) %>%
    reframe(MeanHo = mean(Gl.Ho))
)

Size <- lm(log(City.Area.km2) ~ PercentImpervious + NDVI + Temp + Prec, 
           data = summary)
summary(Size)
Anova(Size, type = 3)

ISA <- lm(PercentImpervious ~ log(City.Area.km2)*Habitat + Pop.change.rate, 
          data = summary)
summary(ISA)
Anova(ISA, type = 3)

temp <- lm(Temp ~ log(City.Area.km2)*Habitat + Pop.change.rate, 
           data = summary)
summary(temp)
Anova(temp, type = 3)

Prec <- lm(Prec ~ log(City.Area.km2)*Habitat + Pop.change.rate, 
           data = summary)
summary(Prec)
Anova(Prec, type = 3)

plot(ISA)

##Add Land use type from MODIS (Friedl, M., Sulla-Menashe, D. (2022). MODIS/Terra+Aqua Land Cover Type Yearly L3 Global 500m SIN Grid V061. NASA EOSDIS Land Processes Distributed Active Archive Center. Accessed 2023-08-22 from https://doi.org/10.5067/MODIS/MCD12Q1.061. Accessed August 22, 2023)

lc.raster <- raster("Rasters/MODIS/MCD12Q1.061_LC_Type1_doy2010001_aid0001.tif")
names(lc.raster) <- "LandUseType"

#Check points CRS matches raster CRS

projection(pts) == projection(lc.raster)

#Create raster list with 50 km buffer around each sample
store_data = list()
for (i in 1:nlayers(lc.raster)){
  store_data[[i]] = raster::extract(lc.raster[[i]], pts)
}

names(store_data) = names(lc.raster)
lc.data = bind_cols(urban.data, as_tibble(store_data))
lc.data

# Check each column for NA values
na.check = map_int(lc.data, ~sum(is.na(.)))
summary(na.check > 0) #FALSE so should be good?

# Remove NA records
lc.data.nona = lc.data %>% drop_na

# Update column values to be biologically meaningful from user guide
str(lc.data.nona)
unique(lc.data.nona$LandUseType)

(lc.data.final <- lc.data.nona %>%
    rename(Longitude = Lon, Latitude = Lat) %>%
    mutate(LandUseType = str_replace(LandUseType, "12", "Cropland"),
           LandUseType = str_replace(LandUseType, "13", "Urban"),
           LandUseType = str_replace(LandUseType, "14", "CropVegMosaic"),
           LandUseType = str_replace(LandUseType, "9", "Savanna"),
           LandUseType = str_replace(LandUseType, "11", "Wetland"),
           LandUseType = str_replace(LandUseType, "5", "MixedForest"),
           LandUseType = str_replace(LandUseType, "10", "Grassland"))
)

str(lc.data.final)
unique(lc.data.final$LandUseType)

#Export data to a csv file and individual id coordinates
alldata <- read.csv("23.08.22_OffspringDataEnv.csv")
lcdata <- merge(alldata, lc.data.final)
head(lcdata)
str(lcdata)
#write to csv
write_csv(lcdata, "23.08.22_OffspringDataEnv.csv")




##Add Land use type from MODIS (Friedl, M., Sulla-Menashe, D. (2022). MODIS/Terra+Aqua Land Cover Type Yearly L3 Global 500m SIN Grid V061. NASA EOSDIS Land Processes Distributed Active Archive Center. Accessed 2023-08-22 from https://doi.org/10.5067/MODIS/MCD12Q1.061. Accessed August 22, 2023)

lc.raster <- raster("Rasters/MODIS/MCD12Q1.061_LC_Type1_doy2010001_aid0001.tif")
names(lc.raster) <- "LandUseType"

#Check points CRS matches raster CRS

projection(pts) == projection(lc.raster)

#Create raster list with 50 km buffer around each sample
store_data = list()
for (i in 1:nlayers(lc.raster)){
  store_data[[i]] = raster::extract(lc.raster[[i]], pts)
}

names(store_data) = names(lc.raster)
lc.data = bind_cols(urban.data, as_tibble(store_data))
lc.data

# Check each column for NA values
na.check = map_int(lc.data, ~sum(is.na(.)))
summary(na.check > 0) #FALSE so should be good?

# Remove NA records
lc.data.nona = lc.data %>% drop_na

# Update column values to be biologically meaningful from user guide
str(lc.data.nona)
unique(lc.data.nona$LandUseType)

(lc.data.final <- lc.data.nona %>%
    rename(Longitude = Lon, Latitude = Lat) %>%
    mutate(LandUseType = str_replace(LandUseType, "12", "Cropland"),
           LandUseType = str_replace(LandUseType, "13", "Urban"),
           LandUseType = str_replace(LandUseType, "14", "CropVegMosaic"),
           LandUseType = str_replace(LandUseType, "9", "Savanna"),
           LandUseType = str_replace(LandUseType, "11", "Wetland"),
           LandUseType = str_replace(LandUseType, "5", "MixedForest"),
           LandUseType = str_replace(LandUseType, "10", "Grassland"))
)

str(lc.data.final)
unique(lc.data.final$LandUseType)

#Export data to a csv file and individual id coordinates
alldata <- read.csv("23.08.22_OffspringDataEnv.csv")
lcdata <- merge(alldata, lc.data.final)
head(lcdata)
str(lcdata)
#write to csv
write_csv(lcdata, "23.08.22_OffspringDataEnv.csv")



#Get future data from wclim
library(geodata)
extent <- c(-81,-79, 43, 44.5) 

clim <- geodata::worldclim_global(var = 'bio', res=2.5, lon=-80, lat=44 , download = F, 
                                path = "/Volumes/OneTouch/Grad School 2015-2021/9CitiesSeq/EnvironmentData/geodata")
crop.clim <- crop(clim,extent) #crop rasters
writeRaster(crop.clim, "/Volumes/OneTouch/Grad School 2015-2021/9CitiesSeq/EnvironmentData/Rasters/wclim.present.tif", overwrite=TRUE)

clim_fut245 <- geodata::cmip6_world(model='ACCESS-ESM1-5', ssp='245', time='2061-2080', 
                                   var='bioc', download=F,
                                   res=2.5, lon=-80, lat=44,
                                   path = "/Volumes/OneTouch/Grad School 2015-2021/9CitiesSeq/EnvironmentData/geodata")
crop.245 <- crop(clim_fut245,extent) #crop rasters
writeRaster(crop.245, "/Volumes/OneTouch/Grad School 2015-2021/9CitiesSeq/EnvironmentData/Rasters/clim.future.245.tif", overwrite=TRUE)

clim_fut370<- geodata::cmip6_world(model='ACCESS-ESM1-5', ssp='370', time='2061-2080', 
                                    var='bioc', download=F,
                                    res=2.5, lon=-80, lat=44,
                                    path = "/Volumes/OneTouch/Grad School 2015-2021/9CitiesSeq/EnvironmentData/geodata")
crop.370 <- crop(clim_fut370,extent) #crop rasters
writeRaster(crop.370, "/Volumes/OneTouch/Grad School 2015-2021/9CitiesSeq/EnvironmentData/Rasters/clim.future.370.tif", overwrite=TRUE)

clim_fut585<- geodata::cmip6_world(model='ACCESS-ESM1-5', ssp='585', time='2061-2080', 
                                    var='bioc', download=F,
                                    res=2.5, lon=-80, lat=44,
                                    path = "/Volumes/OneTouch/Grad School 2015-2021/9CitiesSeq/EnvironmentData/geodata")
crop.585 <- crop(clim_fut585,extent) #crop rasters
writeRaster(crop.585, "/Volumes/OneTouch/Grad School 2015-2021/9CitiesSeq/EnvironmentData/Rasters/clim.future.585.tif", overwrite=TRUE)



#Extract data
#4. Extract data from raster for current conditions
# Load site locations (needs coordinates for each sampling/target site)
sites <- read.csv("../Datasheets/24.05.10_FinalData.csv")


# set x and y for lat and long
x <- sites$Longitude
y <- sites$Latitude

# set up coordinates as a tibble
pts <- cbind(x, y) %>% as_tibble()

#Create raster list with  buffer around each sample 
data <- terra::extract(clim, pts)
data1 <- cbind(pts, data)

head(data1)
names(data1)[1] = "Lon"
names(data1)[2] = "Lat"
names(data1)[4:22] = c(paste0("bio", 1:19))

# Check each column for NA values
na.check = map_int(data1, ~sum(is.na(.)))
summary(na.check > 0) #FALSE so should be good?

# Remove NA records
all.data.nona = data1 %>% drop_na
head(all.data.nona)
#Export data to a csv file and individual id coordinates
head(sites)
sites <- dplyr::select(sites, c("ID", "SampleName", "Order", "PopName"))

climdata <- merge(sites, all.data.nona, by = "ID")
head(climdata)
str(climdata)

#write to csv
write_csv(climdata, "24.06.10_Climate_Present.csv") 

#Repeat for future data
#Create raster list with  buffer around each sample 
fut245 <- terra::extract(clim_fut245, pts)
fut245 <- cbind(pts, fut245)

head(fut245)
names(fut245)[1] = "Lon"
names(fut245)[2] = "Lat"
names(fut245)[4:22] = c(paste0("bio", 1:19))

# Check each column for NA values
na.check = map_int(fut245, ~sum(is.na(.)))
summary(na.check > 0) #FALSE so should be good?

#Export data to a csv file and individual id coordinates
head(sites)
sites <- dplyr::select(sites, c("ID", "SampleName", "Order", "PopName"))

fut245.all <- merge(sites, fut245, by = "ID")
head(fut245.all)
str(fut245.all)

#write to csv
write_csv(fut245.all, "24.06.10_Climate_Fut245.csv") 
