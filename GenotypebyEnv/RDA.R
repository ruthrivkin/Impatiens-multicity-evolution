library(dplyr)
library(raster)
library(vegan)

setwd("/Users/ruthrivkin/Dropbox/Grad School 2015-2020/Impatiens capensis/Data files/9CitiesSeq/Population Genetics/NewAssembly/")

#Load genetic data
library(LEA)

#rerun snmf for ld snps only

vcf2lfmm("snmf/variants_hardfilter.dp10.geno.maf.ld_cityreoder.vcf", "snmf/ld.cityreorder")
genoin <- "snmf/variants_hardfilter.dp10.geno.maf.ld_cityreoder.geno"

snmf.reorder = snmf(genoin, K = 1:15, ploidy = 2, entropy = T,
                    alpha = 100, project = "new", repetitions = 10, seed = 42)
summary(snmf.reorder)
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


project=load.snmfProject("snmf/variants_hardfilter.dp10.geno.maf.ld_cityreoder.snmfProject")
ce = cross.entropy(project, K = 6)
lowest.ce = which.min(ce)

#Impute missing data
impute(project, "snmf/variants_hardfilter.dp10.geno.maf.ld_cityreoder.lfmm", method = 'mode', K = 6, run = lowest.ce)

lfmm_data <- read.lfmm("snmf/variants_hardfilter.dp10.geno.maf.ld_cityreoder.lfmm_imputed.lfmm")

#> Loading required namespace: adegenet
sample.info <- read.csv("../../Datasheets/24.10.21_MasterData.csv", header = T)

coord <- dplyr::select(sample.info, "Longitude", "Latitude")

# Extract environmental vars
env <- dplyr::select(sample.info, "City.Area.km2", "PercentImpervious", "Pop.change.rate", "NDVI", "Temp", "Prec")



#Control for pop structure
#RDA
pca <- read.table("PCA/ld.city.eigenvec")
names(pca)[1] <- "City"
names(pca)[2] <- "SampleName"
names(pca)[3:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-2))

env$PC1 <- pca$PC1
env$PC2 <- pca$PC2

#Check correlation between PC1 and environmental predictors
round(cor(env, env$PC1), 2)
round(cor(env, env$PC2), 2)

# correlation tests for whole dataset
library(Hmisc)
res <- rcorr(as.matrix(env)) # rcorr() accepts matrices only

# display p-values (rounded to 3 decimals)
round(res$P, 3)
#Conclusion most are sign correlated

#Run rda with PC1 alone (only correlations less than 0.5)
env <- subset(env, select = -PC2)
rda <- rda(lfmm_data ~ ., data=env, coords = coord, scale=T)

rda
RsquareAdj(rda) # 0.0336776
summary(rda)$concont
vif.cca(rda) #all below 10

plot(rda, scaling=3)


#make some pretty plots
sample.info$City <- as.factor(sample.info$City)
city <- sample.info$City
c10 <- c(
  "dodgerblue2", "maroon", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6")# lt purple

pdf("Figures/RDA-Axis12-popstructure.pdf")
plot(rda, type="n", scaling=3)
points(rda, display="species", pch=20, cex=0.7, col="lightgrey", scaling=3)           # the SNPs
points(rda, display="sites", pch=21, cex=1.3, col="grey", scaling=3, bg=c10[city]) # the wolves
text(rda, scaling=3, display="bp", col="blue", cex=1)                           # the predictors
legend("bottomleft", legend=levels(city), bty="n", col="gray32", pch=21, cex=1, pt.bg=c10)
dev.off()

pdf("Figures/RDA-Axis13.pdf")
plot(rda, type="n", scaling=3, choices=c(1,3))
points(rda, display="species", pch=20, cex=0.7, col="lightgrey", scaling=3)           # the SNPs
points(rda, display="sites", pch=21, cex=1.3, col="grey", scaling=3, bg=c10[city]) # the wolves
text(rda, scaling=3, display="bp", col="blue", cex=1)                           # the predictors
legend("bottomright", legend=levels(city), bty="n", col="gray32", pch=21, cex=1, pt.bg=c10)
dev.off()

pdf("Figures/RDA-Axis23.pdf")
plot(rda, type="n", scaling=3, choices=c(2,3))
points(rda, display="species", pch=20, cex=0.7, col="lightgrey", scaling=3)           # the SNPs
points(rda, display="sites", pch=21, cex=1.3, col="grey", scaling=3, bg=c10[city]) # the wolves
text(rda, scaling=3, display="bp", col="blue", cex=1)                           # the predictors
legend("bottomright", legend=levels(city), bty="n", col="gray32", pch=21, cex=1, pt.bg=c10)
dev.off()

#Model significance
signif.full <- anova.cca(rda, parallel=getOption("mc.cores"), permutations = how(nperm=10000)) # default is permutation=999
signif.full #p < 0.001


#Check sign axes
signif.axis <- anova.cca(rda, by="axis")
signif.axis #RDA 1,2,3

#FInd outlier snps

load.rda <- scores(rda, choices=c(1:3), display="species")  # Species scores for the first three constrained axes

#RDA axes have pretty normal distn
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA2")

#Find outliers with 3 standard deviation cutoff (two-tailed p-value = 0.012 (FDR = 1%)
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}


cand1 <- outliers(load.rda[,1],3) # 110
cand2 <- outliers(load.rda[,2],3) # 39
cand3 <- outliers(load.rda[,3],3) # 26

ncand <- length(cand1) + length(cand2) + length(cand3)
ncand #175 outliers total

#Make dataframe of outlier snps
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

#Add in the correlations of each candidate SNP with the  environmental predictors:

foo <- matrix(nrow=(ncand), ncol=7)  # 6 columns for 6 predictors
colnames(foo) <- c("City.Area.km2", "PercentImpervious", "Pop.change.rate", "NDVI", "Temp", "Prec", "PC1")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- lfmm_data[,nam]
  foo[i,] <- apply(env,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

#Investigate candidates
#Check dups
length(cand$snp[duplicated(cand$snp)]) #4

foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2]) #none on RDA1
table(foo[foo[,1]==2,2]) #0 on RDA2
table(foo[foo[,1]==3,2]) #4 on RDA2

#remove
cand <- cand[!duplicated(cand$snp),] 
length(cand$snp) #171

#Count which snps correlated most strongly with which predictors
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,11] <- names(which.max(abs(bar[4:10]))) # gives the variable
  cand[i,12] <- max(abs(bar[4:10]))              # gives the correlation
}

colnames(cand)[11] <- "predictor"
colnames(cand)[12] <- "correlation"

table(cand$predictor) 
#Most snps most strongly correlated with City area and temp

#City.Area.km2          NDVI           PC1          Prec          Temp 
#50                      3           112             2             4        

#Plot outliers

sel <- cand$snp
env <- cand$predictor
env[env=="City.Area.km2"] <- '#ffff33'
env[env=="Pop.change.rate"] <- '#6a3d9a'
env[env=="PercentImpervious"] <- 'orange'
env[env=="NDVI"] <- '#33a02c'
env[env=="Temp"] <- '#e31a1c'
env[env=="Prec"] <- '#1f78b4'
env[env=="PC1"] <- 'pink'

# color by predictor:
col.pred <- rownames(rda$CCA$v) # pull the SNP names

for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("V",col.pred)] <- '#f1eef6' # non-candidate SNPs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#ffff33','#6a3d9a','orange','#33a02c','#e31a1c','#1f78b4', "pink")


# axes 1 & 2
pdf("Figures/RDA-outlier-Axis12-popstructure.pdf")
plot(rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(rda, display="species", pch=21, cex=1, col="grey", bg=col.pred, scaling=3)
points(rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(rda, scaling=3, display="bp", col="black", cex=1)
legend("bottomright", legend=c("City.Area.km2", "PercentImpervious", "Pop.change.rate", "NDVI", "Temp", "Prec", "PC1"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
dev.off()

#axes 1 and 3
pdf("Figures/RDA-outlier-Axis13-popstructure.pdf")
plot(rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(1,3))
points(rda, display="species", pch=21, cex=1, col="grey", bg=col.pred, scaling=3, choices=c(1,3))
points(rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
text(rda, scaling=3, display="bp", col="black", cex=1, choices=c(1,3))
legend("bottomright", legend=c("City.Area.km2", "PercentImpervious", "Pop.change.rate", "NDVI", "Temp", "Prec", "PC1"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
dev.off()

#axes 2 and 3
pdf("Figures/RDA-outlier-Axis23-popstructure.pdf")
plot(rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(2,3))
points(rda, display="species", pch=21, cex=1, col="grey", bg=col.pred, scaling=3, choices=c(2,3))
points(rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(2,3))
text(rda, scaling=3, display="bp", col="black", cex=1, choices=c(2,3))
legend("bottomright", legend=c("City.Area.km2", "PercentImpervious", "Pop.change.rate", "NDVI", "Temp", "Prec", "PC1"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
dev.off()


#Which environments are most strongly cor
intersetcor(rda)[,1:3]

#PCAAdapt
library(pcadapt)

lfmm_data.pca <- "snmf/variants_hardfilter.dp10.geno.maf.ld_cityreoder.lfmm"
pca.ad <- read.pcadapt(lfmm_data.pca, type = "lfmm")

#Find pop structure first, find K
findk <- pcadapt(input = pca.ad, K = 20) #k should be large to start
plot(findk, option = "screeplot") #looks like k=3 is the best fit

#Plot scores
plot(findk, option = "scores", pop = pop.info$City) #looks like Toronto split, consistent with plink pca. Going to stick with K = 2

plot(findk, option = "scores", i = 3, j = 4, pop = sample.info$City) #just a blur, no pop structure here

#Detect outliers
outliers <- pcadapt(pca.ad, K = 3)
summary(outliers)

#definitely outliers
plot(outliers , option = "manhattan")
plot(outliers, option = "qqplot")
hist(outliers$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

#Plot test stat D
plot(outliers, option = "stat.distribution")


#Choose cutoff point
library(qvalue)
#FDR 10%
qval <- qvalue(outliers$pvalues)$qvalues
alpha <- 0.01
outliers.qv <- which(qval < alpha)
length(outliers.qv) #126


#Benjamini-Hochberg Procedure 
padj <- p.adjust(outliers$pvalues, method="BH")
alpha <- 0.01
outliers.BH <- which(padj < alpha)
length(outliers.BH) #126, the same snps are id'd


#snmf outlier detection
project=load.snmfProject("snmf/variants_hardfilter.dp10.geno.maf.ld_cityreoder.snmfProject")

p <- snmf.pvalues(project, entropy=TRUE, ploidy=2, K=6)

padj.snmf <- p.adjust(p$pvalues, method="BH")
alpha <- 0.01
snmf.BH <- which(padj.snmf < alpha)
length(snmf.BH) #459

#create dataframe with contig information
# load map file
snp_info <- read.table("snmf/variants_hardfilter.dp10.geno.maf.ld_cityreorder.bim")

# pull out only CHROM and POS
snp_info1 <- snp_info %>%
  rename(Contig = V1, 
        ID = V2,
        Poscm = V3,
        Pos = V4,
        All1 = V5,
        All2 = V6) %>%
  mutate(Order = c(1:3965))

snp_info1 

#Merge with candidate snps

#snmf
snmf.BH <- as.data.frame(snmf.BH)
names(snmf.BH) <- "Order"

snmf.outliers <- left_join(snmf.BH, snp_info1)

#pcadapt
outliers.BH <- as.data.frame(outliers.BH)
names(outliers.BH) <- "Order"

pcadapt.outliers <- left_join(outliers.BH, snp_info1)

#RDA
library(tidyverse)
cand.merge <- dplyr::select(cand,c("snp", "axis", "predictor"))
cand.merge$snp<-gsub("V","", as.character(cand.merge$snp))
colnames(cand.merge)[1] <- "Order"

cand.merge$Order<- as.integer(cand.merge$Order)

rda.outliers <- left_join(cand.merge, snp_info1)

overlap.rda.snmf <- merge(rda.outliers, snmf.outliers)
length(overlap.rda.snmf$Order) #130

overlap.pca.snmf <- merge(pcadapt.outliers, snmf.outliers)
length(overlap.pca.snmf$Order) #117

overlap.pca.rda<- merge(pcadapt.outliers, rda.outliers)
length(overlap.pca.rda$Order) #13


overlap.all <- merge(overlap.pca.snmf, rda.outliers) #0
length(overlap.all$Order) #13

#Get list of outlier loci regions for blast
write_csv(overlap.all, "RDA/all_outliers.txt")

#Create bed for searching genome (formate chr start end)
overlap.rda.snmf <- read.csv("RDA/RDA_snmf_outliers.txt")

bed <- dplyr::select(overlap.all, c("Contig", "Pos"))

#Search withing 100 bp of each location
bed100bp <- bed %>%
  mutate(Start = Pos - 100,
         End = Pos + 100) %>%
  dplyr::select(Contig, Start, End)
names(bed100bp) <- NULL


bed1kb <- bed %>%
  mutate(Start = Pos - 1000,
         End = Pos + 1000) %>%
  dplyr::select(Contig, Start, End)
names(bed1kb) <- NULL

bed10kb <- bed %>%
  mutate(Start = Pos - 10000,
         End = Pos + 10000) %>%
  dplyr::select(Contig, Start, End)
names(bed10kb) <- NULL

readr::write_delim(bed100bp, "RDA/outlierlocs100bp.txt", delim = " ", col_names = F)
readr::write_delim(bed1kb, "RDA/outlierlocs1kb.txt", delim = " ", col_names = F)
readr::write_delim(bed10kb, "RDA/outlierlocs10kb.txt", delim = " ", , col_names = F)





