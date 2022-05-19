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

# setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\konza projects\\conSME\\data\\biomass') #laptop
# setwd('C:\\Users\\komatsuk\\Dropbox (Smithsonian)\\konza projects\\conSME\\data\\biomass') #desktop


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
trt <- read.csv('C:\\Users\\komatsuk\\Dropbox (Smithsonian)\\konza projects\\conSME\\data\\conSME_treatments.csv')

biomass2019 <- read.csv('conSME_biomass_2019.csv')%>%
  mutate(pdead=0)
biomass2020 <- read.csv('conSME_biomass_2020.csv')%>%
  mutate(drop=ifelse(block=='I'&plot==6&strip==1, 1, 0))%>% #drop missing sample
  filter(drop!=1)%>%
  select(-drop)
biomass2021 <- read.csv('conSME_biomass_2021.csv')%>%
  mutate(experiment='conSME')%>%
  rename(gram=grass)%>%
  select(-date)%>%
  filter(notes!='no biomass in any bag') #filter out missing sample

biomass <- rbind(biomass2019, biomass2020, biomass2021)%>%
  rename(project_name=experiment)%>%
  left_join(trt)%>%
  select(-notes)%>%
  mutate_all(~replace(., is.na(.), 0))%>% #don't be scared by the errors, they just don't add NAs to the factor columns (watershed, trt, etc)
  mutate(total=gram+forb+woody)

biomassMean <- biomass%>%
  group_by(watershed, block, plot, year, bison, small_mammal, invertebrates, trt)%>%
  summarise(gram=mean(gram)*10, forb=mean(forb)*10, woody=mean(woody)*10, pdead=mean(pdead)*10, total=mean(total)*10)%>%
  ungroup()

#subsetting out the first year of trts, which is different in patterns from all subsequent years
biomassLater <- biomassMean%>%
  filter(year!=2019)
  

##### checking for outliers #####
#outliers have been confirmed as true values for 2019-2021
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


##### ANOVA - total biomass #####
summary(totBiomassModel <- lme(total~watershed*trt*year,
                               data=biomassLater,
                               random=~1|block/trt,
                               correlation=corCompSymm(form=~year|block/trt), 
                               control=lmeControl(returnObject=T)))
anova.lme(totBiomassModel, type='sequential') 
emmeans(totBiomassModel, pairwise~year*trt, adjust="tukey")

#N1A total biomass
summary(totBiomassN1AModel <- lme(total~trt*year,
                               data=subset(biomassLater, watershed=='N1A'),
                               random=~1|block/trt,
                               correlation=corCompSymm(form=~year|block/trt), 
                               control=lmeControl(returnObject=T)))
anova.lme(totBiomassN1AModel, type='sequential') 
emmeans(totBiomassN1AModel, pairwise~trt, adjust="tukey")

#N4B total biomass
summary(totBiomassN4BModel <- lme(total~trt*year,
                                  data=subset(biomassLater, watershed=='N4B'),
                                  random=~1|block/trt,
                                  correlation=corCompSymm(form=~year|block/trt), 
                                  control=lmeControl(returnObject=T)))
anova.lme(totBiomassN4BModel, type='sequential') 
emmeans(totBiomassN4BModel, pairwise~trt, adjust="tukey")

#figure - total biomass by year
ggplot(data=barGraphStats(data=biomassLater, variable="total", byFactorNames=c("trt", "year")), aes(x=trt, y=mean)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black', fill='white') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
  scale_x_discrete(limits=c('BSI', 'BSX', 'XSI', 'XXI', 'XSX', 'XXX')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  coord_cartesian(ylim=c(0,750)) +
  facet_grid(cols=vars(year))
#export at 1400x600


##### ANOVA - graminoid biomass #####
summary(gramBiomassModel <- lme(gram~watershed*trt*year,
                               data=biomassLater,
                               random=~1|block/trt,
                               correlation=corCompSymm(form=~year|block/trt), 
                               control=lmeControl(returnObject=T)))
anova.lme(gramBiomassModel, type='sequential') 
emmeans(gramBiomassModel, pairwise~watershed*trt*year, adjust="tukey")

#N1A graminoid biomass
summary(gramBiomassN1AModel <- lme(gram~trt*year,
                                  data=subset(biomassLater, watershed=='N1A'),
                                  random=~1|block/trt,
                                  correlation=corCompSymm(form=~year|block/trt), 
                                  control=lmeControl(returnObject=T)))
anova.lme(gramBiomassN1AModel, type='sequential') 
emmeans(gramBiomassN1AModel, pairwise~trt*year, adjust="tukey")

#N4B graminoid biomass
summary(gramBiomassN4BModel <- lme(gram~trt*year,
                                  data=subset(biomassLater, watershed=='N4B'),
                                  random=~1|block/trt,
                                  correlation=corCompSymm(form=~year|block/trt), 
                                  control=lmeControl(returnObject=T)))
anova.lme(gramBiomassN4BModel, type='sequential') 
emmeans(gramBiomassN4BModel, pairwise~trt*year, adjust="tukey")

#figure - graminoid biomass
ggplot(data=barGraphStats(data=biomassLater, variable="gram", byFactorNames=c("trt", "watershed", "year")), aes(x=trt, y=mean)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black', fill='white') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Graminoid Biomass (g m'^'-2',')'))) +
  scale_x_discrete(limits=c('BSI', 'BSX', 'XSI', 'XXI', 'XSX', 'XXX')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  coord_cartesian(ylim=c(0,700)) +
  facet_grid(rows=vars(watershed), cols=vars(year))
#export at 1400x1200



##### ANOVA - forb biomass #####
summary(forbBiomassModel <- lme(forb~watershed*trt*year,
                                data=biomassLater,
                                random=~1|block/trt,
                                correlation=corCompSymm(form=~year|block/trt), 
                                control=lmeControl(returnObject=T)))
anova.lme(forbBiomassModel, type='sequential') 
emmeans(forbBiomassModel, pairwise~trt, adjust="tukey")

#N1A forb biomass
summary(forbBiomassN1AModel <- lme(forb~trt*year,
                                   data=subset(biomassLater, watershed=='N1A'),
                                   random=~1|block/trt,
                                   correlation=corCompSymm(form=~year|block/trt), 
                                   control=lmeControl(returnObject=T)))
anova.lme(forbBiomassN1AModel, type='sequential') 
emmeans(forbBiomassN1AModel, pairwise~trt, adjust="tukey")

#N4B forb biomass
summary(forbBiomassN4BModel <- lme(forb~trt*year,
                                   data=subset(biomassLater, watershed=='N4B'),
                                   random=~1|block/trt,
                                   correlation=corCompSymm(form=~year|block/trt), 
                                   control=lmeControl(returnObject=T)))
anova.lme(forbBiomassN4BModel, type='sequential') 
emmeans(forbBiomassN4BModel, pairwise~trt, adjust="tukey")

#figure - forb biomass across years and watersheds
ggplot(data=barGraphStats(data=biomassLater, variable="forb", byFactorNames=c("trt")), aes(x=trt, y=mean)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black', fill='white') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Forb Biomass (g m'^'-2',')'))) +
  scale_x_discrete(limits=c('BSI', 'BSX', 'XSI', 'XXI', 'XSX', 'XXX')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  annotate('text', x=1, y=170, label='a', size=10) +
  annotate('text', x=2, y=166, label='a', size=10) +
  annotate('text', x=3, y=70, label='b', size=10) +
  annotate('text', x=4, y=95, label='b', size=10) +
  annotate('text', x=5, y=130, label='ab', size=10) +
  annotate('text', x=6, y=150, label='ab', size=10)
#export at 500x600


##### ANOVA - woody biomass #####
summary(woodyBiomassModel <- lme(woody~trt*year*watershed,
                                 data=biomassLater,
                                 random=~1|block/trt,
                                 correlation=corCompSymm(form=~year|block/trt), 
                                 control=lmeControl(returnObject=T)))
anova.lme(woodyBiomassModel, type='sequential') 
emmeans(woodyBiomassModel, pairwise~trt*year, adjust="tukey")


##### ANOVA - pdead biomass #####
#N4B pdead biomass
summary(pdeadBiomassN4BModel <- lme(pdead~trt*year,
                                  data=subset(biomassLater, watershed=='N4B'),
                                  random=~1|block/trt,
                                  correlation=corCompSymm(form=~year|block/trt), 
                                  control=lmeControl(returnObject=T)))
anova.lme(pdeadBiomassN4BModel, type='sequential') 
emmeans(pdeadBiomassN4BModel, pairwise~trt*year, adjust="tukey")

#figure - pdead biomass all years
ggplot(data=barGraphStats(data=subset(biomassLater, watershed=='N4B'), variable="pdead", byFactorNames=c("trt", "year", "watershed")), aes(x=trt, y=mean)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Litter Biomass (g m'^'-2',')'))) +
  scale_x_discrete(limits=c('BSI', 'BSX', 'XSI', 'XXI', 'XSX', 'XXX')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  facet_grid(cols=vars(year), rows=vars(watershed))
#export at 1400x600


# ##### ANCOVA - total biomass coupled with graminoid biomass? #####
# summary(biomassRegressions <- lme(total~trt*gram,
#                                   data=biomassLater,
#                                   random=~1|block/trt))
# anova.lme(biomassRegressions)
# 
# #figure - regression
# ggplot(data=biomassLater, aes(x=gram, y=total, color=trt)) +
#   geom_point(size=5) +
#   geom_smooth(method='lm', se=F)
# 
#   ylab(expression(paste('Litter Biomass (g m'^'-2',')'))) +
#   scale_x_discrete(limits=c('BSI', 'BSX', 'XSI', 'XXI', 'XSX', 'XXX')) +
#   theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
#   facet_grid(cols=vars(year), rows=vars(watershed))
# #export at 1400x600

