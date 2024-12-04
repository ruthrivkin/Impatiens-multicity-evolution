setwd("/Users/ruthrivkin/Dropbox/Grad School 2015-2020/Impatiens capensis/Data files/9CitiesSeq/Population Genetics/NewAssembly/")

library(ggplot2)

ng <- theme(aspect.ratio=0.7,panel.background = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border=element_rect(fill = NA, color = "NA"),
            axis.line = element_line(linewidth=0.5), 
            axis.ticks=element_line(color="black"),
            axis.text=element_text(color="black",size=10, margin = 0.5), 
            axis.title=element_text(color="black",size=1), 
            axis.title.y=element_text(vjust=1,size=12),
            axis.title.x=element_text(vjust=0.1,size=12),
            axis.text.x=element_text(size=10),
            axis.text.y=element_text(size=10))


#snmf
library(LEA)
library(tidyverse)


# Create geno file
system("~/plink/plink --file variants_hardfilter.dp10.geno.maf.ld_city --allow-extra-chr --recode vcf --out snmf/variants_hardfilter.dp10.geno.maf.ld_city")

#Reorder samples so cities are in order
system("~/plink/plink --file variants_hardfilter.dp10.geno.maf.ld_city --indiv-sort f snmf/Reordersamples.txt --allow-extra-chr -make-bed --out snmf/variants_hardfilter.dp10.geno.maf.ld_cityreorder")

system("~/plink/plink --bfile snmf/variants_hardfilter.dp10.geno.maf.ld_cityreorder --allow-extra-chr --recode vcf --out snmf/variants_hardfilter.dp10.geno.maf.ld_cityreoder")



vcf2lfmm("snmf/variants_hardfilter.dp10.geno.maf.ld_cityreoder.vcf", "snmf/ld.cityreorder")

#snmf
#Load geno
pop.info <- read.csv("../../Datasheets/24.10.21_Master.Data.csv", header = T)
str(pop.info)
pop.info <- dplyr::arrange(pop.info, City)
#pop.info <- subset(pop.info, City != "Toronto"), no effect of removing toronto

genoin <- "snmf/variants_hardfilter.dp10.geno.maf.ld.hwe_cityreoder.geno"


##identify best clusters
library(reshape2)
library(dplyr)

snmf.reorder = snmf(genoin, K = 1:15, ploidy = 2, entropy = T,
                alpha = 100, project = "new", repetitions = 10, seed = 42)
summary(snmf.reorder)

snmf = load.snmfProject("snmf/variants_hardfilter.dp10.geno.maf.ld.hwe_city.snmfProject")

K <- summary(snmf.reorder)$crossEntropy %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("temp") %>%
  mutate(K = as.numeric(str_extract(temp, "(\\d)+"))) %>%
  dplyr::select(-temp)

# choose omptimal (K = 6)
ggplot(K, aes(x = K, y = mean)) +
  geom_line(color = "black", size = 0.25 ) +
  geom_segment(aes(x = K, y = min, xend = K, yend = max)) +
  geom_point(shape = 21, size = 4, color = "black", fill = "blue") +
  scale_x_continuous(breaks = seq(0, 25, by = 1)) +
  labs(x = "Number of ancestral populations", y = "Cross-entropy criterion") +
  ng
ggsave("../../Manuscript/Figures/Cross_Entropy.pdf", width = 6.78, height = 5.3, dpi = 300)

ce = cross.entropy(snmf.reorder, K = 6)
ce
lowest.ce = which.min(ce)
lowest.ce



#plot ancestry matrix
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
   "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown", "black")
  

pdf("snmf/Ind Admixture K=6.pdf", width = 6.78, height = 5.3)
barchart(snmf.reorder, K = 6, run = lowest.ce, sort.by.Q = FALSE,
         border = NA, space = 0,
         col = c25,
         xlab = "Individuals",
         ylab = "Ancestry proportions K = 6",
         main = "Ancestry matrix") -> bp

bp$city <- pop.info$City
axis(1, at = 1:length(bp$city), 
     labels = bp$city, las = 3, 
     cex.axis = .4)

dev.off()



qmatrix = as.data.frame(Q(snmf.reorder, K = 6, run = lowest.ce))
head(qmatrix)

# Label column names of qmatrix
ncol(qmatrix)
cluster_names = c()
for (i in 1:ncol(qmatrix)){
  cluster_names[i] = paste("Cluster", i)
}
cluster_names
colnames(qmatrix) = cluster_names
head(qmatrix)

# Add individual IDs
qmatrix$Ind = pop.info$SampleName

# Add site IDs
qmatrix$Site = pop.info$City
head(qmatrix)

# Convert dataframe to long format

qlong = reshape2::melt(qmatrix, id.vars=c("Ind","Site"))
head(qlong)

# Change order of sites by using the factor function
#site.order = c("BB","DS", "FB","GB","KB","LS","MC","NB","NW","SB","SH", "VM","WH")
qlong$Site= as.factor(qlong$Site)
levels(qlong$Site)

# Define colour palette
#cols = cols(length(unique(qlong$variable)))

admix.bar = ggplot(data=qlong, aes(x=Ind, y=value, fill=variable))+
  geom_bar(stat = "identity")+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = c25)+
  ylab("Admixture proportion")+
  # xlab("Individual")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(colour="black", size=12),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
admix.bar

# Plot admixture barplot 
admix.bar = ggplot(data=qlong, aes(x=Ind, y=value, fill=variable))+
  geom_bar(stat = "identity")+
  scale_y_continuous(expand = c(0,0))+
  facet_wrap(~Site, scales = "free", ncol = 2)+
  scale_fill_manual(values = c25)+
  ylab("Admixture proportion")+
  # xlab("Individual")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(colour="black", size=12),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
admix.bar

ggsave("City_admixture_barplot_offspring.pdf", width = 6.78, height = 5.3, dpi=300)


# ----------------- #
#
# Prepare pie charts
#
# ----------------- #

# Calculate mean admixture proportions for each site
head(qmatrix)
clusters = grep("Cluster", names(qmatrix)) # indexes of cluster columns
avg_admix = aggregate(qmatrix[, clusters], list(qmatrix$Site), mean)

# Order alphabetically by site
avg_admix = avg_admix[order(as.character(avg_admix$Group.1)), ]
avg_admix

# Convert dataframe from wide to long format
avg_admix = reshape2::melt(avg_admix, id.vars = "Group.1")
head(avg_admix)

# Define a function to plot pie charts using ggplot for each site
pie_charts = function(admix_df, site, cols){
  # admix_df = dataframe in long format of admixture proportions per site 
  # site = string 
  # cols = vector of colours of length(clusters)
  ggplot(data = subset(admix_df, Group.1 == site),
         aes(x = "", y = value, fill = variable))+
    geom_bar(width = 1, stat = "identity", colour = "black", show.legend = FALSE)+
    coord_polar(theta = "y")+
    scale_fill_manual(values = c25)+
    theme_void()
}

# Test function on one site
pie_charts(avg_admix, site = "Georgetown", cols = c25)

# Apply function to all sites using for loop
site.ordered = sort(unique(qlong$Site))

pies = list()
for (i in site.ordered){
  pies[[i]] = pie_charts(admix_df = avg_admix, site = i, cols = c25) 
}


# ----------------- #
#
# Prepare basemap
#
# ----------------- #

# Estimate city coordinates
city.coords <- pop.info  %>% 
  group_by(City)  %>% 
  dplyr::summarize(N = n(),
                   Mean.Lon = mean(Longitude),
                   Mean.Lat = mean(Latitude)
  )
city.coords


# Order alphabetically by site]
coords = city.coords[order(city.coords$City), ] 
coords

# Check order matches coords order
as.character(avg_admix$Group.1) == as.character(coords$City)

# Set map boundary (xmin, xmax, ymin, ymax)
library(sf)
library(raster)
boundary = extent(-81,-79, 43, 44) 
boundary

# Get map outlines from rworldmap package
library(rworldmap)
library(rworldxtra)
library(rgeos)
library(maps)
library(ggsn)
library(maptools)
library(grid)
library(miscTools)
library(stringr)
library(ggpubr)


map.outline = getMap(resolution = "high")

# Crop to boundary and convert to dataframe
map.outline = raster::crop(map.outline, y = boundary) %>% fortify()


# Plot basemap
basemap = ggplot()+
  geom_polygon(data=map.outline, aes(x=long, y=lat, group=group), fill="white",
               colour="black", size=0.5)+
  coord_quickmap(expand=F)+
  ggsn::north(map.outline, symbol = 10, scale = 0.06, location = "topleft")+
  ggsn::scalebar(data = map.outline, dist = 10, dist_unit = "km", height = 0.01,
                 transform = TRUE, model = "WGS84", 
                 location = "bottomleft", anchor = c(x = -79.3, y = 43.01),
                 st.bottom = FALSE, st.size = 3, st.dist = 0.015)+
  xlab("Longitude")+
  ylab("Latitude")+
  ng
basemap

# Add pie charts to basemap
#
# ----------------- #

# Extract coordinates for each site
coord.list = list()

for (i in site.ordered){
  coord.list[[i]] = c(subset(coords, City == i)$Mean.Lon, subset(coords, City == i)$Mean.Lat)
}
coord.list

# Define pie chart sizes
radius = 0.1

# Convert ggplot pie charts to annotation_custom layers
pies.ac = list()
for (i in 1:length(site.ordered)){
  pies.ac[[i]] = annotation_custom(grob = ggplotGrob(pies[[i]]),
                                   xmin = coord.list[[i]][[1]] - radius,
                                   xmax = coord.list[[i]][[1]] + radius,
                                   ymin = coord.list[[i]][[2]] - radius,
                                   ymax = coord.list[[i]][[2]] + radius)
}


pie.map = basemap + pies.ac
pie.map
ggsave("snmf/Admixture_pie_charts_map_city.pdf", width = 6.78, height = 5.3, dpi = 300)


#Assign Individual cluster based on 75% ancestry proportion
names(qmatrix) <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5", "Cluster6", "Site", "Ind")
cluster.assign <- qmatrix %>% 
  mutate(
    Cluster = case_when(
      Cluster1 > .75 ~ "Cluster 1",
      Cluster2 > .75 ~ "Cluster 2",
      Cluster3 > .75 ~ "Cluster 3",
      Cluster4 > .75 ~ "Cluster 4",
      Cluster5 > .75 ~ "Cluster 5",
      Cluster6 > .75 ~ "Cluster 6",
          .default = "Admixed"
    )
  )

write.csv(cluster.assign, "snmf/24.10.08_snmfcluster.csv", row.names = F)

snmf.cluster <- read.csv("snmf/24.10.08_snmfcluster.csv")

#append to master data
pop.info$Cluster <- snmf.cluster$Cluster
data.all$Cluster <- pop.info$Cluster

