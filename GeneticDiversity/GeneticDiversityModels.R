library(ggplot2)
library(dplyr)
library(tidyverse)
library(RColorBrewer) # needed for some extra colours in one of the graphs
library(ggthemes)
library(ggridges)
library(car)

c25 <- c(
  "dodgerblue2", "maroon", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "#E31A1C",  "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown", "black")

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

setwd("/Users/ruthrivkin/Dropbox/Grad School 2015-2020/Impatiens capensis/Data files/9CitiesSeq/Population Genetics/NewAssembly/")
#options(contrasts=c("contr.sum", "contr.poly"))

#Estimate Heterozygosity and Nucleotide diversity

library(adegenet)
library(dartR)

#Calculate Ho and pop Fis using in genlight
pop.gl <- gl.load("Fst/Popld.hw.gl")

#Calculate genetic diversity at pop level
gl.compliance.check(pop.gl)
pop.gl.nona <- gl.impute(pop.gl, method="neighbour")
gl.compliance.check(pop.gl.nona)

#Ho
pop.het <- gl.report.heterozygosity(pop.gl, method = "ind")
colnames(pop.het)[1:2] =c("SampleName", "ind.Ho")
pop.het <- dplyr::select(pop.het, "SampleName", "ind.Ho")

#Fis
pop.f <-  gl.report.heterozygosity(pop.gl.nona, method = "pop")
colnames(pop.f)[1] ="PopName"
pop.f <- dplyr::select(pop.f, "PopName", "FIS", "Ho")

#Private alleles
pa <- gl.report.pa(pop.gl, method = "one2rest")
colnames(pa)[3] = "PopName"
colnames(pa)[10] = "PrivateAlleles"

pa <- dplyr::select(pa, "PopName", "PrivateAlleles", "AFD")


#Append to complete dataset
data <- read.csv("../../Datasheets/24.10.21_Master.Data.csv") #Complete dataset (I hope)

data.all <- merge(data, pop.het)
data.all <- left_join(data.all, pop.f)
data.all <- left_join(data, pa)
str(data.all)





#Update sample sheet
write.csv(data.all, "../../Datasheets/24.10.21_MasterData.csv")  



#Run popgen analyses
data <- read.csv("../../Datasheets/24.10.21_MasterData.csv") #Complete dataset (I hope)
data$LandUseType <- as.factor(data$LandUseType)
data$Habitat <- as.factor(data$Habitat)
data$City <- as.factor(data$City)
data$Type <- as.factor(data$Type)
data$BuiltArea <- as.factor(data$BuiltArea)
data$Cluster <- as.factor(data$Cluster)



#Plot
hist(data$Ho)
hist(data$ind.Ho)
hist(data$FIS)
hist(data$Pi)
hist(data$AFD)
hist(data$PrivateAlleles)

cor(data$Temp, data$PercentImpervious)
cor.test(data$Ho, data$PrivateAlleles)


#plot data
c10 <- c(
  "dodgerblue2", "maroon", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6")# lt purple





#Look at population level data
(summary <- data %>%
    group_by(PopName, City, Habitat, Latitude, Longitude,
             City.Census.2021,City.Area.km2,Pop.change.rate, 
             PercentImpervious, BuiltArea, LandUseType, Prec, Temp, NDVI,
             Ho, FIS, Pi, PrivateAlleles) %>%
    reframe(count = n())
)

cor.test(summary$City.Census.2021, summary$City.Area.km2) #0.9963754 
cor.test(summary$PercentImpervious, summary$NDVI) #-0.6980012 
cor.test(summary$Prec, summary$Temp) #-0.6000745 
cor.test(summary$PercentImpervious, summary$City.Area.km2) #-0.2033429 
cor.test(summary$Pi, summary$PrivateAlleles)



#Get summaries
summary(summary$Pi)
summary(summary$Ho)
summary(summary$PrivateAlleles)


#Models
Pi <- lm(Pi ~ PercentImpervious + NDVI + Temp + Prec + log(City.Area.km2) + Pop.change.rate,
         data = summary)
summary(Pi)
Anova(Pi, type = 3)

Ho <- lm(Ho ~ PercentImpervious + NDVI + Temp + Prec +log(City.Area.km2) + Pop.change.rate, 
         data = summary)
summary(Ho)
Anova(Ho, type = 3)

coef(Ho)["log(City.Area.km2)"]/100 #1.462766e-05

Pa <- lm(PrivateAlleles ~ PercentImpervious + NDVI + Temp + Prec + log(City.Area.km2) + Pop.change.rate,
          data = summary)
summary(Pa)
Anova(Pa, type = 3)


cor.test(summary$Pi, summary$Ho)
#Check residuals in space

res <- Pi$residuals
cor.test(res, summary$Latitude) #ns cor
cor.test(res, summary$Longitude) #ns cor
plot(summary$Latitude, res) #no trend


res <- Ho$residuals
cor.test(res, summary$Latitude) #ns cor
cor.test(res, summary$Longitude) #ns cor
plot(summary$Latitude, res) #no trend


res <- Pa$residuals
cor.test(res, summary$Latitude) #ns cor
cor.test(res, summary$Longitude) #ns cor
plot(summary$Latitude, res) #no trend

#Check autocorrelation
library(lmtest)
dwtest(Ho) #ns
dwtest(Pi) #ns
dwtest(Pi) #ns


#All looks okay
vif(Pi)
vif(Ho)
vif(Pa)

plot(Pi)
plot(Ho)
plot(Pa)

#Plots 
(ggplot(data  = summary,
        aes(x = Prec,
            y = FIS ,
            col = Habitat))+ #to add the colours for different classes
    geom_point(size     = 1.2,
               alpha    = .8,
               position = "jitter")+ #to add some random noise for plotting purposes
    geom_smooth(method = lm,
                se     = TRUE, 
                col    = "black",
                size   = .5, 
                alpha  = .8)+ # to add regression line
    theme(legend.position = "right") + 
    labs(y = "FIS", x = "Precipitation (mm)") +
    scale_color_manual(values = c25) +
    ng
)

ggsave("Fis_by_Prec.pdf", width = 6.78, height = 5.3, dpi = 300)

(ggplot(data  = summary,
        aes(x = log(City.Area.km2),
            y =Pi_Pixy,
            col = City)) + #to add the colours for different classes
    geom_point(size     = 1.2,
               alpha    = .8,
               position = "jitter")+ #to add some random noise for plotting purposes
    geom_smooth(method = lm,
                se     = TRUE, 
                col    = "black",
                size   = .5, 
                alpha  = .8)+ # to add regression line
    theme(legend.position = "right") + 
    labs(y = "Pi", x = "log(City.Area.km2)") +
    scale_color_manual(values = c10) +
    ng
)
ggsave("../../Figures/Pi_pixy.cityarea.pdf", width = 6.78, height = 5.3, dpi = 300)




#Old analyses
#Create a pca variable summarizing city size
library(ggfortify)

city.size<-data[,c(13:15)]  #this is just the city size variables
city.pca<-prcomp(city.size,scale=TRUE) # result is quite intuitive. nice separation on PC1 and PC2 
summary(city.pca)
#                         PC1    PC2     PC3
# Proportion of Variance 0.7195 0.2797 0.00074
autoplot(city.pca, data = data, colour = "City", loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3)

#extract variables to make plots of the eigenvectors and scores
loadings<-city.pca$rotation # "rotation" is what R calls the PCA loadings
rownames(loadings)<-colnames(city.size)
scores<-city.pca$x # "x" is what R calls the species scores

#Append city pc1 score to dataset
data$City.Pca <- (scores[,1])
#Save as variable




#H 
library(lme4)
library(car)

ind.Ho <- lmer(ind.Ho~ PercentImpervious + NDVI + Temp + Prec + LandUseType + log(City.Area.km2) + Pop.change.rate + Cluster+ (1|PopName), 
               data = data)
summary(ind.Ho)
Anova(ind.Ho)

(ggplot(data  = data,
        aes(x = PercentImpervious,
            y = ind.Ho ,
            col = City))+ #to add the colours for different classes
    geom_point(size     = 1.2,
               alpha    = .8,
               position = "jitter")+ #to add some random noise for plotting purposes
    geom_smooth(method = lm,
                se     = TRUE, 
                col    = "black",
                size   = .5, 
                alpha  = .8)+ # to add regression line
    theme(legend.position = "right") + 
    labs(y = "Individual Ho", x = "Percent Impervious Surface") +
    scale_color_manual(values = c10) +
    ng
)
ggsave("Figures/Ho_by_PercImp.pdf", width = 6.78, height = 5.3, dpi = 300)

(ggplot(data  = data,
        aes(x = log(City.Area.km2),
            y = ind.Ho ,
            col = City))+ #to add the colours for different classes
    geom_point(size     = 1.2,
               alpha    = .8,
               position = "jitter")+ #to add some random noise for plotting purposes
    geom_smooth(method = lm,
                se     = TRUE, 
                col    = "black",
                size   = .5, 
                alpha  = .8)+ # to add regression line
    theme(legend.position = "right") + 
    labs(y = "Individual Ho", x = "log(City Size)") +
    scale_color_manual(values = c10) +
    ng
)
ggsave("Figures/Ho_by_PercImp.pdf", width = 6.78, height = 5.3, dpi = 300)

(ggplot(data  = data,
        aes(x = Prec,
            y = ind.Ho ,
            col = City))+ #to add the colours for different classes
    geom_point(size     = 1.2,
               alpha    = .8,
               position = "jitter")+ #to add some random noise for plotting purposes
    geom_smooth(method = lm,
                se     = TRUE, 
                col    = "black",
                size   = .5, 
                alpha  = .8)+ # to add regression line
    theme(legend.position = "right") + 
    labs(y = "Individual Ho", x = "Precipitation (mm)") +
    scale_color_manual(values = c10) +
    ng
)
ggsave("Figures/Ho_by_Precimp.pdf", width = 6.78, height = 5.3, dpi = 300)

(Ho.plot <- ggplot(data, aes(x=City, y=ind.Ho)) + 
    geom_boxplot(outlier.size = 1, colour = c10) + 
    geom_dotplot(binaxis='y', stackdir ='center', dotsize=0.01) +
    geom_jitter(shape=1, position=position_jitter(0.2)) + 
    xlab("City") +
    ylab("Ho")+
    ng
)

ggsave("Figures/Ind_Ho_by_City.pdf", width = 6.78, height = 5.3, dpi = 300)

(Ho.plot <- ggplot(data, aes(x=City, y=Ho)) + 
    geom_boxplot(outlier.size = 1, colour = c10) + 
    geom_dotplot(binaxis='y', stackdir ='center', dotsize=0.01) +
    geom_jitter(shape=1, position=position_jitter(0.2)) + 
    xlab("City") +
    ylab("Ho")+
    ng
)
ggsave("Pop_Ho_by_City.pdf", width = 6.78, height = 5.3, dpi = 300)

(Fis.plot <- ggplot(data, aes(x=City, y=FIS) + 
                      geom_boxplot(outlier.size = 1, colour = c10)) + 
    geom_dotplot(binaxis='y', stackdir ='center', dotsize=0.01) +
    geom_jitter(shape=1, position=position_jitter(0.2)) + 
    xlab("City") +
    ylab("FIS")+
    ng
)

ggsave("Figures/Fis_by_City.pdf", width = 6.78, height = 5.3, dpi = 300)

#BRMS analysis
#Heterozygosity
library(brms)
library(report)

#Imp surface and Habitat and city size
priorsm1 <- c(
  prior(normal(0, 1), class = b, coef = PercentImpervious),
  prior(normal(0, 1), class = b, coef = NDVI),
  prior(normal(0, 1), class = b, coef = Prec),
  prior(normal(0, 1), class = b, coef = Temp),
  prior(normal(0, 1), class = b, coef = City.Area.km2),      
  prior(normal(0, 1), class = b, coef = Pop.change.rate)
)

model1 <- brm(Pi ~  PercentImpervious + NDVI + Temp + Prec + City.Area.km2 + Pop.change.rate, 
              data = summary, 
              prior = priorsm1,
              sample_prior = TRUE,
              warmup = 3000, iter = 10000,
              cores = 4, chains = 4) #to run the model
summary(model1)
report(model1)


model2 <- brm(FIS ~  PercentImpervious*Habitat + City.Pca,  
              data = summary, 
              prior = priorsm1,
              sample_prior = TRUE,
              warmup = 3000, iter = 10000, 
              cores = 4, chains = 4) #to run the model

summary(model2)
report(model2)

model3 <- brm(Gl.Ho ~  PercentImpervious + Habitat + City + (1|PopName),  
              data = data, 
              prior = priorsm1,
              sample_prior = TRUE,
              warmup = 3000, iter = 10000, max_treedepth = 10,
              cores = 4, chains = 4) #to run the model
summary(model3)
report(model3)

#Check model output

pairs(model1)

#relationship of the predictors with the response
plot(conditional_effects(model1, effects = "PercentImpervious"))
plot(conditional_effects(model1, effects = "Habitat"))
plot(conditional_effects(model1, effects = "City.Pca"))

bayesplot::mcmc_intervals(as.array(model1), pars=c("b_PercentImpervious",
                                                   "b_HabitatUrban",
                                                   "b_City.Pca",
                                                   "b_PercentImpervious:HabitatUrban")
)



#Extract posterior values
post_model1 <-  as_draws_array(model1)
post_sum_model1 <- posterior::summarize_draws(post_model1)

#Hypothesis testing
h <- c(" PercentImpervious > 0", "LandUseTypeGrassland < 0", "LandUseTypeUrban < 0")
(hyp1 <- hypothesis(model1, h))
plot(hyp1, ignore_prior = TRUE)


model1 <- brm(prob.obs.het ~ Outcrossing  + PercentImpervious + LandUseType + (1|City:Pop),  
              data = full, 
              prior = priorsm1,
              sample_prior = TRUE,
              warmup = 1000, iter = 3000, 
              cores = 3, chains = 3) #to run the model

summary(model1)
prior_summary(model1)

#Check model output
pairs(model1)

#relationship of the predictors with the response
plot(conditional_effects(model1, effects = "PercentImpervious"))
plot(conditional_effects(model1, effects = "Outcrossing"))

bayesplot::mcmc_intervals(as.array(model1), pars=c("b_Outcrossing",
                                                   "b_PercentImpervious",
                                                   "b_LandUseTypeGrassland",
                                                   "b_LandUseTypeCropVegMosaic", 
                                                   "b_LandUseTypeMixedForest",
                                                   "b_LandUseTypeSavanna",
                                                   "b_LandUseTypeUrban",
                                                   "b_LandUseTypeWetland")
)

bayesplot::mcmc_areas(as.array(model1), pars=c("b_Outcrossing", "b_PercentImpervious"), prob=.95)


#Extract posterior values
post_model1 <-  as_draws_array(model1)
post_sum_model1 <- posterior::summarize_draws(post_model1)

#Hypothesis testing
h <- c("Outcrossing > 0", " PercentImpervious > 0", "LandUseTypeGrassland < 0", "LandUseTypeUrban < 0")
(hyp1 <- hypothesis(model1, h))
plot(hyp1, ignore_prior = TRUE)


#Get allelic richness
library(hierfstat)

