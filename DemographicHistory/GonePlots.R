setwd("~/Dropbox/Grad School 2015-2020/Impatiens capensis/Data files/9CitiesSeq/Population Genetics/NewAssembly/")
library(dartRverse)

#Load gl file that has been filtered for HWE and LD with city level info
gl <- gl.load("Fst/Cityld.hw.gl")
gl <- gl.compliance.check(gl)
gl@pop

ne <- dartR.popgen::gl.LDNe(gl, outfile="cityLD.txt", 
              outpath=getwd(),
              neest.path ="~/Downloads/NeEstimator/", 
              critical=c(0,0.05), singleton.rm=TRUE, mating='random',
              Waples.correction = "nChromosomes", Waples.correction.value = 20,
              naive = TRUE)


#Compare with strataG
options(repos = c(
  zkamvar = 'https://zkamvar.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'))

install.packages('strataG')

library(vcfR)
library(adegenet)
library(strataG)


# convert genpop
gl <- gl.load("Fst/Popld.hw.gl")

genind <- gl2gi(gl)

# convert genind to gtypes
snps_gtypes <- genind2gtypes(genind)
#class(snps_gtypes)

# Estimating Ne using ldNe (from https://github.com/jdalapicolla/Ne_StrataG.R/blob/master/Ne_Estimation.R)
Ne <- ldNe(snps_gtypes, 
           maf.threshold=0, 
           by.strata=TRUE, 
           ci=0.95, 
           drop.missing=TRUE,
           num.cores=1)
Ne


#Investigate demographic history
#Plot GONE results
library(tidyverse)

#Map gone ne results
# Use the 100 iterations to make a 95% confidence interval (this is based off of Kardos et al. 2023's script: https://github.com/martykardos/KillerWhaleInbreeding/blob/main/FigureCode/rCode_Fig_ED_1.R)

library(scales)
library(matrixStats)

# BB first (high arctic)
# load all the iteration files and put it in a matrix
files <- paste("/Volumes/OneTouch/Grad School 2015-2021/9CitiesSeq/NewAssembly/GONE/TEMPORARY_FILES_All/outfileLD_TEMP/outfileLD_",1:100,"_GONE_Nebest",sep="")
NeMat <- NULL
for(i in 1:100){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}

# Cet CI for the recent 200 generations
NeCI <- matrix(NA,nrow=200,ncol=2)
for(i in 1:200){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.05,0.95))
}

# Set up data to get ready for ggplot
NeCI_dat <- as.data.frame(NeCI)
NeCI_dat$Generation <- 1:nrow(NeCI_dat)

# Median
Ne_med <- as.data.frame(rowMedians(NeMat[1:200,]))
colnames(Ne_med) <- "median"
Ne_med$Generation <- 1:nrow(Ne_med)


# Define generation time b/c we want to plot in years
gen=1 

(all <- ggplot()+
    # High Arctic
    geom_ribbon(data=NeCI_dat, aes(x=Generation*gen, ymin=V1, ymax=V2), fill="#92C5DE", alpha=0.5)+
    geom_line(data=Ne_med, aes(x=(Generation*gen),y=median), color="#4D9BD3", lwd=1.5)+
    # Cutting off at 150 generations
    xlim(0,200*gen)+
#    ylim(0,15000)+
    theme_bw()+
    xlab("Years ago")+
    ylab(expression(paste("Effective population size (",italic(""*N*"")[e],")",sep="")))+
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
    ng
)

ggsave("Figures/Gone_NE_all_withCI.pdf", width = 6.78, height = 5.3, dpi = 300)



# Load data for individual cities
BF <- read.delim("GONE/Output_Ne_BF", header=T, skip = 1)
head(BF)
BT <- read.delim("GONE/Output_Ne_BT", header=T, skip = 1)
head(BT)
CM <- read.delim("GONE/Output_Ne_CM", header=T, skip = 1)
head(CM)
GT <- read.delim("GONE/Output_Ne_GT", header=T, skip = 1)
head(GT)
GU <- read.delim("GONE/Output_Ne_GU", header=T, skip = 1)
head(GU)
HM <- read.delim("GONE/Output_Ne_HM", header=T, skip = 1)
head(HM)
KT <- read.delim("GONE/Output_Ne_KT", header=T, skip = 1)
head(KT)
MT <- read.delim("GONE/Output_Ne_MT", header=T, skip = 1)
head(MT)
OJ <- read.delim("GONE/Output_Ne_OJ", header=T, skip = 1)
head(OJ)
TO <- read.delim("GONE/Output_Ne_TO", header=T, skip = 1)
head(TO)
c10 <- c(
  "dodgerblue2", "maroon", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6")# lt purple
ng <- theme(aspect.ratio=0.7,panel.background = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_line(linewidth=0.5), 
            axis.ticks=element_line(color="black"),
            axis.text=element_text(color="black",size=12, margin = 0.5), 
            axis.title=element_text(color="black",size=1), 
            axis.title.y=element_text(vjust=1,size=14),
            axis.title.x=element_text(vjust=0.1,size=14),
            axis.text.x=element_text(size=12),
            axis.text.y=element_text(size=12),
            # legend.position = "top", legend.direction="horizontal", 
            legend.text=element_text(size=12), legend.key = element_rect(fill = "white"), 
            legend.title = element_text(size=12),legend.key.size = unit(1, "cm"))



# Quick plot of the means (cut off at 200 generation time)
ggplot()+
  geom_line(data=BF, aes(x=(Generation),y=Geometric_mean), color="#2166AC", lwd=1)+
  geom_line(data=BT, aes(x=(Generation),y=Geometric_mean), color="maroon", lwd=1)+
  geom_line(data=CM, aes(x=(Generation),y=Geometric_mean), color="green4", lwd=1)+
  geom_line(data=GT, aes(x=(Generation),y=Geometric_mean), color="#6A3D9A", lwd=1)+
  geom_line(data=GU, aes(x=(Generation),y=Geometric_mean), color="#FF7F00", lwd=1)+
  geom_line(data=HM, aes(x=(Generation),y=Geometric_mean), color="gold1", lwd=1)+
  geom_line(data=KT, aes(x=(Generation),y=Geometric_mean), color="skyblue2", lwd=1)+
  geom_line(data=MT, aes(x=(Generation),y=Geometric_mean), color="#FB9A99", lwd=1)+
  geom_line(data=OJ, aes(x=(Generation),y=Geometric_mean), color="palegreen2", lwd=1)+
  geom_line(data=TO, aes(x=(Generation),y=Geometric_mean), color="#CAB2D6", lwd=1)+
  ng +
  labs(y = "Effective Population Size", x = "Generation") +
  xlim(0,400)
ggsave("Figures/GONE_Ne.pdf")

