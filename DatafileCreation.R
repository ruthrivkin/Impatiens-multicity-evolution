library(dplyr)
library(tidyverse)
library(ggplot2)
library(RColorBrewer) # needed for some extra colours in one of the graphs
library(ggthemes)
library(ggridges)
library(readxl)

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

setwd("/Users/ruthrivkin/Dropbox/Grad School 2015-2020/Impatiens capensis/Data files/9CitiesSeq/Datasheets/")

#Create new dataset for individuals that passed qc
#Full dataset (458 samples)
data <- read.csv("RefGenome_SampleSummary.csv")
str(data)

data.filtered <- data %>%
  select(SampleName, PopName, Subpop, Population, Plant, Family, Type, ReadsMapped, MappingRate)
str(data.filtered)


#Summarize pops/Cities to see if data is still evenly distributed
(summary.pop <- data %>% 
    group_by(PopName) %>% 
    reframe(N = n())
)

#Change name of toronto pops to match other samples
data$PopName[data$PopName == "TO_KSR"] <- "TO_R1"
data$PopName[data$PopName == "TO_WLN"] <- "TO_R2"
data$PopName[data$PopName == "TO_ORT"] <- "TO_R3"
data$PopName[data$PopName == "TO_MUD"] <- "TO_U1"
data$PopName[data$PopName == "TO_YCP"] <- "TO_U2"
data$PopName[data$PopName == "TO_TCS"] <- "TO_U3"

#Pull population information (environment, city size, etc)
env <- read.csv("23.08.22_OffspringDataEnv.csv")
(env.pop <- env %>% 
    group_by(Latitude, Longitude, PopName, City, City.Census.2021, Pop.change.rate, 
             City.Area.km2, Habitat, PercentImpervious, BuiltArea, LandUseType) %>% 
    reframe(N = n())
)
str(env.pop)

#Create new data set that has all population level data for all samples
samples.env <- inner_join(data, env.pop, by = "PopName")
str(samples.env)

write.csv(samples.env, "23.12.12_RefGenome_SampleSummary.csv")



#Create Seed Dataset

full <- read.csv("24.01.05_RefGenome_SampleSummary.csv")
full.filtered <- full %>%
  dplyr::select(c(1,5:22))
str(full.filtered)


seed.sample <- read_excel("SampleSummary.xlsx")
seed.filtered <- seed.sample %>%
  filter(SeedSampleRun == "Yes")  %>%
  dplyr::select(SampleName)

str(seed.filtered)

seedsamples <- inner_join(full.filtered, seed.filtered, by = "SampleName")
str(seedsamples)



#Create heterozygosity files and population pi
het <- read.delim("../Population Genetics/Genetics Files/ld.hwe.city.het.txt", na.strings = " ")
het <- het %>% rename(SampleName= IID)
head(het)


het.formated <- het %>%
  mutate( Ho = 1-(O.HOM./N.NM.),
          He = 1-(E.HOM./N.NM.)) %>%
  dplyr::select(SampleName, F, Ho, He)

#Take absolute value of F
het.formated$F <- abs(het.formated$F)   
head(het.formated)

seedsamples.het <- inner_join(seedsamples, het.formated, by = "SampleName")
str(seedsamples.het)

write.csv(seedsamples.het, "24.02.26_SeedSamples.csv")

#Create pi strata
strata <- seedsamples.het %>%
  select(PopName, SampleName) 
colnames(strata) <- c("POP_ID", "INDIVIDUALS")
write_delim(strata, "../Population Genetics/Genetics Files/SeedSamples.Strata.tsv")




#Now add pi
pi <- read.delim("../Population Genetics/Genetics Files/pi_20240226@0656/pi.populations.tsv", na.strings = "")
pi <- pi[-c(54),] #remove overall pi
colnames(pi) <- c("PopName", "Pi")

#now join with main dataframe
popgen <- inner_join(seedsamples.het, pi)
head(popgen)
dim(popgen)

write.csv(popgen, "24.02.26_RefGenome_SampleSummary.csv")


#Write new tsv strata file for radiator
samples.env <- read.csv("23.12.12_RefGenome_SampleSummary.csv")

strata <- samples.env %>%
  select(PopName, SampleName, City, Habitat) 
colnames(strata) <- c("POP_ID", "INDIVIDUALS", "STRATA", "STRATA")

write_delim(strata, "../Population Genetics/strata.all.tsv")
#Had to update pop names because they were slightly different in vcf import in radiator. No idea why


#Add pixy pi

pi <- read.delim("../Population Genetics/Genetics Files/pi_20240509@1135/all.sites.pi.populations.tsv", na.strings = "")

data <- read.csv("24.03.25_FinalData.csv") #Complete dataset (I hope)

head(data)


data.pi <- inner_join(data, pi)
data.pi

write_csv(data.pi, "24.05.10_FinalData.csv")
