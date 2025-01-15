setwd("/Users/ruthrivkin/Dropbox/Grad School 2015-2020/Impatiens capensis/Data files/9CitiesSeq/Population Genetics/NewAssembly/")

#Memgene

#load sample info and extract relevant pop data
sample.info <- read.csv("../../Datasheets/24.10.21_MasterData.csv", header = T)
pop.info <- sample.info[order(sample.info$ID),]


#import files - offspring analysis
library("adegenet")
library("dartR")

infile.city<- "variants_hardfilter.dp10.geno.maf.ld.hwe_city.raw"
myData.city <- read.PLINK(infile.city)
myData.city@pop
gl.save(myData.city, "Fst/Cityld.hw.gl")

city.gl <- gl.load("Fst/Cityld.hw.gl")


#Calculate basic stats
library(StAMPP)

#load data
city.gl@pop

fst <- stamppFst(city.gl, nboots = 1000, percent = 95, nclusters = 1)
fst.matrix <- as.data.frame(fst$Fsts)
p.matrix <- as.data.frame(fst$Pvalues)

write.csv(fst.matrix, "fst/city.fst.matrix.csv")
write.csv(p.matrix, "fst/city.pfst.matrix.csv")

fst.boots <- as.data.frame(fst$Bootstraps) # p=0 for all


#Repeat for population differentiation

infile.pop<- "variants_hardfilter.dp10.geno.maf.ld.hwe_pop.raw"
myData.pop <- read.PLINK(infile.pop)
myData.pop@pop
gl.save(myData.pop, "Fst/Popld.hw.gl")

pop.gl <- gl.load("Fst/Popld.hw.gl")


#load data
pop.gl@pop

pop.fst <- stamppFst(pop.gl, nboots = 100, percent = 95, nclusters = 1)
pop.fst.matrix <- as.data.frame(pop.fst$Fsts)
pop.p.matrix <- as.data.frame(pop.fst$Pvalues)

write.csv(pop.fst.matrix, "fst/pop.fst.matrix.csv")
write.csv(pop.p.matrix, "fst/pop.pfst.matrix.csv")

pop.fst.boots <- as.data.frame(pop.fst$Bootstraps) # p=0 for all

