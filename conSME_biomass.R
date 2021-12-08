################################################################################
##  conSME_biomass.R: Analysis of aboveground biomass responses in conSME experiment.
##
##  Author: Kimberly Komatsu
##  Date created: December 8, 2021
################################################################################

library(PerformanceAnalytics)
library(nlme)
library(emmeans)
library(tidyverse)

setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\konza projects\\conSME\\data\\biomass')


#set options
options(contrasts=c('contr.sum','contr.poly'))

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

#homemade functions
#barGraphStats(data=, variable="", byFactorNames=c(""))

barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  


##### data #####
trt <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\konza projects\\conSME\\data\\conSME_treatments.csv')

biomass2019 <- read.csv('conSME_biomass_2019.csv')%>%
  mutate(pdead=0)
biomass2020 <- read.csv('conSME_biomass_2020.csv')%>%
  mutate(drop=ifelse(block=='I'&plot==6&strip==1, 1, 0))%>% #drop missing sample
  filter(drop!=1)%>%
  select(-drop)

biomass <- rbind(biomass2019, biomass2020)%>%
  rename(project_name=experiment)%>%
  left_join(trt)%>%
  filter(forb!='< 1g', forb!='<1 g', forb!='<1g', woody!='<1g')%>%
  mutate(forb=as.numeric(forb), woody=as.numeric(woody))%>%
  mutate_all(~replace(., is.na(.), 0))%>%
  mutate(total=gram+forb+woody)

biomassMean <- biomass%>%
  group_by(watershed, block, plot, year, bison, small_mammal, invertebrates, trt)%>%
  summarise(gram=mean(gram)*10, forb=mean(forb)*10, woody=mean(woody)*10, pdead=mean(pdead)*10, total=mean(total)*10)%>%
  ungroup()
  

##### checking for outliers #####
dataVis <- biomass%>%
  select(gram, forb, woody, pdead, total) #make visualization dataframe
chart.Correlation(dataVis, histogram=T, pch=19)

chart.Correlation(subset(biomass, trt=='BSI')%>%
                    select(gram, forb, woody, pdead, total), histogram=T, pch=19)
chart.Correlation(subset(biomass, trt=='BSX')%>%
                    select(gram, forb, woody, pdead, total), histogram=T, pch=19)
chart.Correlation(subset(biomass, trt=='XSI')%>%
                    select(gram, forb, woody, pdead, total), histogram=T, pch=19)
chart.Correlation(subset(biomass, trt=='XSX')%>%
                    select(gram, forb, woody, pdead, total), histogram=T, pch=19)
chart.Correlation(subset(biomass, trt=='XXI')%>%
                    select(gram, forb, woody, pdead, total), histogram=T, pch=19)
chart.Correlation(subset(biomass, trt=='XXX')%>%
                    select(gram, forb, woody, pdead, total), histogram=T, pch=19)



##### ANOVAs #####
#total biomass
summary(totBiomassModel <- lme(total~watershed*trt,
                               data=biomassMean,
                               random=~1|block/plot,
                               correlation=corCompSymm(form=~year|block/plot), 
                               control=lmeControl(returnObject=T)))
anova.lme(totBiomassModel, type='sequential') 
emmeans(totBiomassModel, pairwise~watershed*trt, adjust="tukey")

#N1A total biomass
summary(totBiomassN1AModel <- lme(total~trt,
                               data=subset(biomassMean, watershed=='N1A'),
                               random=~1|block/plot,
                               correlation=corCompSymm(form=~year|block/plot), 
                               control=lmeControl(returnObject=T)))
anova.lme(totBiomassN1AModel, type='sequential') 
emmeans(totBiomassN1AModel, pairwise~trt, adjust="tukey")

#N4B total biomass
summary(totBiomassN4BModel <- lme(total~trt,
                                  data=subset(biomassMean, watershed=='N4B'),
                                  random=~1|block/plot,
                                  correlation=corCompSymm(form=~year|block/plot), 
                                  control=lmeControl(returnObject=T)))
anova.lme(totBiomassN4BModel, type='sequential') 
emmeans(totBiomassN4BModel, pairwise~trt, adjust="tukey")

#figure - total biomass
ggplot(data=barGraphStats(data=biomassMean, variable="total", byFactorNames=c("trt", "watershed")), aes(x=trt, y=mean)) +
  geom_bar(position=position_dodge(0.1), size=5, stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1)) +
  facet_wrap(~watershed)
#export at 500x600

#figure - graminoid biomass
ggplot(data=barGraphStats(data=biomassMean, variable="gram", byFactorNames=c("trt", "watershed")), aes(x=trt, y=mean)) +
  geom_bar(position=position_dodge(0.1), size=5, stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Graminoid Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1)) +
  facet_wrap(~watershed)
#export at 500x600

#figure - forb biomass
ggplot(data=barGraphStats(data=biomassMean, variable="forb", byFactorNames=c("trt", "watershed")), aes(x=trt, y=mean)) +
  geom_bar(position=position_dodge(0.1), size=5, stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Forb Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1)) +
  facet_wrap(~watershed)
#export at 500x600