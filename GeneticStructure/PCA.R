
setwd("/Users/ruthrivkin/Dropbox/Grad School 2015-2020/Impatiens capensis/Data files/9CitiesSeq/Population Genetics/Genetics Files/")

library(ggplot2)
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

ng <- theme(aspect.ratio=0.7,panel.background = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border=element_rect(fill = NA, color = "NA"),
            axis.line = element_line(linewidth = 0.5), 
            axis.ticks=element_line(color="black"),
            axis.text=element_text(color="black",size=10, margin = 0.5), 
            axis.title=element_text(color="black",size=1), 
            axis.title.y=element_text(vjust=1,size=12),
            axis.title.x=element_text(vjust=0.1,size=12),
            axis.text.x=element_text(size=10),
            axis.text.y=element_text(size=10))

#Generate pca and mdist files for pca
system("~/plink/plink --file variants_hardfilter.dp10.geno.maf.ld_city --allow-extra-chr --distance-matrix --pca --out PCA/ld.city")

## Load data
dist_populations<-read.table("PCA/ld.city.mdist",header=F)
### Extract breed names
fam <- data.frame(City=read.table("PCA/ld.city.mdist.id")[,1])
### Extract individual names 
famInd <- data.frame(SampleName=read.table("PCA/ld.city.mdist.id")[,2])

## Perform PCA using the cmdscale function 
# Time intensive step - takes a few minutes with the 4.5K animals
mds_populations <- cmdscale(dist_populations,eig=T,3)

## Extract the eigen vectors
eigenvec_populations <- cbind(fam,famInd,mds_populations$points)

## Proportion of variation captured by each eigen vector
eigen_percent <- round(((mds_populations$eig)/sum(mds_populations$eig))*100,2)


# Visualize PCA in tidyverse
#load sample info and extract relevant pop data
sample.info <- read.csv("../../Datasheets/24.03.25_FinalData.csv", header = T)
pops <- dplyr::select(sample.info, "SampleName", "Habitat")

eigen_pop.city <- merge(pops, eigenvec_populations)


# PCA plot
(pca <- ggplot(data = eigen_pop.city) +
  geom_point(mapping = aes(x = `1`, y = `2`,color = City, shape = Habitat), show.legend = TRUE , size = 2.5) + 
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(x = paste0("Principal component 1 (",eigen_percent[1]," %)"),
       y = paste0("Principal component 2 (",eigen_percent[2]," %)")) + 
  scale_color_manual(values = c25) +
  ng
)

ggsave("PCA_plink_mdist.pdf", width = 6.78, height = 5.3, dpi = 300)



#PLink PCA
# read in data
library(readr)
pca <- read_table("PCA/ld.city.eigenvec", col_names = FALSE)
eigenval <- scan("PCA/ld.city.eigenval")


# sort out the pca data

# set names
names(pca)[1] <- "City"
names(pca)[2] <- "SampleName"
names(pca)[3:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-2))

pca_pop.city <- merge(pops, pca)

(pcaplot <- ggplot(data = pca_pop.city) +
    geom_point(mapping = aes(x = PC1, y = PC2,color = City, shape = Habitat), show.legend = TRUE, size = 2.5 ) + 
    geom_hline(yintercept = 0, linetype="dotted") + 
    geom_vline(xintercept = 0, linetype="dotted") +
    labs(x = paste0("Principal component 1 (",round(eigenval[1],2),"%)"),
         y = paste0("Principal component 2 (",round(eigenval[2],2),"%)")) + 
    scale_color_manual(values = c25) +
    ng
)
ggsave("PCA/PCA_plink_pca.pdf", width = 6.78, height = 5.3, dpi = 300)


