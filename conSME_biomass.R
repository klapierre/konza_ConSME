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

setwd('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\konza projects\\conSME\\data') #desktop


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
trt <- read.csv('conSME_treatments.csv')

biomass2019 <- read.csv('biomass\\conSME_biomass_2019.csv') %>%
  mutate(pdead=0)
biomass2020 <- read.csv('biomass\\conSME_biomass_2020.csv') %>%
  mutate(drop=ifelse(block=='I'&plot==6&strip==1, 1, 0)) %>% #drop missing sample
  filter(drop!=1) %>%
  select(-drop)
biomass2021 <- read.csv('biomass\\conSME_biomass_2021.csv') %>%
  mutate(experiment='conSME') %>%
  rename(gram=grass) %>%
  select(-date) %>%
  filter(notes!='no biomass in any bag') #filter out missing sample
biomass2022 <- read.csv('biomass\\conSME_biomass_2022.csv') %>%
  mutate(experiment='conSME') %>%
  rename(gram=grass) %>%
  filter(!is.na(gram)) #filter out missing samples


biomass <- rbind(biomass2019, biomass2020, biomass2021, biomass2022) %>%
  rename(project_name=experiment) %>%
  left_join(trt) %>%
  select(-notes) %>%
  mutate_all(~replace(., is.na(.), 0)) %>% #don't be scared by the errors, they just don't add NAs to the factor columns (watershed, trt, etc)
  mutate(total=gram+forb+woody) %>%
  ungroup() %>% mutate(ws_label=ifelse(watershed=='N1A', 'Annual', '4 Year'))

biomassMean <- biomass %>%
  group_by(watershed, ws_label, block, plot, year, bison, small_mammal, invertebrates, trt) %>%
  summarise(gram=mean(gram)*10, forb=mean(forb)*10, woody=mean(woody)*10, pdead=mean(pdead)*10, total=mean(total)*10)

#subsetting out the first year of trts, which is different in patterns from all subsequent years
biomassLater <- biomassMean %>%
  filter(year!=2019)
  

# ##### checking for outliers #####
# #outliers have been confirmed as true values for 2019-2021
# dataVis <- biomass%>%
#   select(gram, forb, woody, pdead, total) #make visualization dataframe
# chart.Correlation(dataVis, histogram=T, pch=19)
# 
# chart.Correlation(subset(biomass, trt=='BSI') %>%
#                     select(gram, forb, woody, pdead, total), histogram=T, pch=19)
# chart.Correlation(subset(biomass, trt=='BSX') %>%
#                     select(gram, forb, woody, pdead, total), histogram=T, pch=19)
# chart.Correlation(subset(biomass, trt=='XSI') %>%
#                     select(gram, forb, woody, pdead, total), histogram=T, pch=19)
# chart.Correlation(subset(biomass, trt=='XSX') %>%
#                     select(gram, forb, woody, pdead, total), histogram=T, pch=19)
# chart.Correlation(subset(biomass, trt=='XXI') %>%
#                     select(gram, forb, woody, pdead, total), histogram=T, pch=19)
# chart.Correlation(subset(biomass, trt=='XXX') %>%
#                     select(gram, forb, woody, pdead, total), histogram=T, pch=19)


##### ANOVA - total biomass #####

summary(totBiomassModel <- lme(total~watershed*year*invertebrates*bison + watershed*year*invertebrates*small_mammal,
                               data=biomass,
                               random=~1|block/trt,
                               correlation=corCompSymm(form=~year|block/trt), 
                               control=lmeControl(returnObject=T)))
anova.lme(totBiomassModel, type='sequential') 
emmeans(totBiomassModel, pairwise~invertebrates*bison, adjust="tukey") #ws*year*small_mammal and ws*year*bison effects, invert*bison interaction

#figure - total biomass bison*ws*year
ggplot(data=barGraphStats(data=biomass, variable="total", byFactorNames=c("bison", "ws_label", "year")), aes(x=bison, y=mean)) +
  geom_bar(position=position_dodge(), size=2, stat="identity", color="black", fill="white") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  facet_grid(cols=vars(year), rows=vars(ws_label))
#export at 1400x1200

#figure - total biomass small_mammal*ws*year
ggplot(data=barGraphStats(data=biomass, variable="total", byFactorNames=c("small_mammal", "ws_label", "year")), aes(x=small_mammal, y=mean)) +
  geom_bar(position=position_dodge(), size=2, stat="identity", color="black", fill="white") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  coord_cartesian(ylim=c(0,90)) +
  facet_grid(cols=vars(year), rows=vars(ws_label))
#export at 1400x1200

#bison*invertebrate interaction
temp <- biomass%>%
  group_by(bison, invertebrates) %>%
  summarise(mean=mean(total), sd=sd(total), N=length(total)) %>%
  ungroup() %>%
  mutate(trt=ifelse(bison=='B'&invertebrates=='I', 'BI', ifelse(bison=='B'&invertebrates=='X', 'B', ifelse(bison=='X'&invertebrates=='I', 'I', 'X')))) %>%
  mutate(se=sd/sqrt(N))

ggplot(data=temp, aes(x=trt, y=mean)) +
  geom_bar(stat='identity', color='black', fill='white', size=2) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  scale_x_discrete(limits=c('BI', 'B', 'I', 'X')) +
  xlab('') + ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  coord_cartesian(ylim=c(0,60))


# #figure - total biomass by year
# ggplot(data=barGraphStats(data=biomass, variable="total", byFactorNames=c("trt", "year")), aes(x=trt, y=mean)) +
#   geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black', fill='white') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
#   ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
#   scale_x_discrete(limits=c('BSI', 'BSX', 'XSI', 'XXI', 'XSX', 'XXX')) +
#   theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
#   # coord_cartesian(ylim=c(0,750)) +
#   facet_grid(cols=vars(year))
# #export at 1400x600


##### total biomass - response ratios #####
###bison effect
bisonBiomassResponse <- biomass%>%
  filter(trt %in% c('XSI', 'BSI', 'XSX', 'BSX')) %>%
  mutate(comparison=ifelse(trt %in% c('XSI', 'BSI'), 'with inverts', 'without inverts'))

bisonBiomassResponse2 <- bisonBiomassResponse%>%
  filter(watershed=='N1A' & year>2019) %>%
  group_by(bison) %>%
  summarise(mean=mean(total)) %>%
  ungroup() 
bisonBiomassResponse3 <- bisonBiomassResponse%>%
  filter(watershed=='N4B') %>%
  group_by(bison) %>%
  summarise(mean=mean(total)) %>%
  ungroup() 
#bison consumption: 32% of biomass in annual and 48% of biomass in 4yr

# #with inverts
# summary(bisonBiomassModel <- lme(total~watershed*trt*year,
#                                 data=subset(bisonBiomassResponse, comparison='with inverts'),
#                                 random=~1|block/trt,
#                                 correlation=corCompSymm(form=~year|block/trt), 
#                                 control=lmeControl(returnObject=T)))
# anova.lme(bisonBiomassModel, type='sequential') 
# emmeans(bisonBiomassModel, pairwise~watershed*trt*year, adjust="tukey") #watershed*trt*year effect
# 
# ggplot(data=barGraphStats(data=subset(bisonBiomassResponse, comparison='with inverts'), variable="total", byFactorNames=c("bison", "watershed", "year")), aes(x=bison, y=mean)) +
#   geom_bar(position=position_dodge(), size=2, stat="identity", color="black", fill="white") +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
#   ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
#   theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
#   coord_cartesian(ylim=c(0,90)) +
#   facet_grid(cols=vars(year), rows=vars(watershed))
# #export at 1400x1200
# 
# #without inverts
# summary(bisonBiomassModel <- lme(total~watershed*trt*year,
#                                 data=subset(bisonBiomassResponse, comparison='withput inverts'),
#                                 random=~1|block/trt,
#                                 correlation=corCompSymm(form=~year|block/trt), 
#                                 control=lmeControl(returnObject=T)))
# anova.lme(bisonBiomassModel, type='sequential') 
# emmeans(bisonBiomassModel, pairwise~watershed*trt*year, adjust="tukey") #watershed*trt*year effect
# 
# ggplot(data=barGraphStats(data=subset(bisonBiomassResponse, comparison='without inverts'), variable="total", byFactorNames=c("bison", "watershed", "year")), aes(x=bison, y=mean)) +
#   geom_bar(position=position_dodge(), size=2, stat="identity", color="black", fill="white") +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
#   ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
#   theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
#   coord_cartesian(ylim=c(0,90)) +
#   facet_grid(cols=vars(year), rows=vars(watershed))
# #export at 1400x1200

#interaction
summary(bisonBiomassModel <- lme(total~watershed*bison*invertebrates*year,
                                data=bisonBiomassResponse,
                                random=~1|block/trt,
                                correlation=corCompSymm(form=~year|block/trt), 
                                control=lmeControl(returnObject=T)))
anova.lme(bisonBiomassModel, type='sequential') 
emmeans(bisonBiomassModel, pairwise~bison, adjust="tukey") #no interaction - can plot just bison effect (above)

ggplot(data=barGraphStats(data=bisonBiomassResponse, variable="total", byFactorNames=c("bison", "ws_label", "year")), aes(x=bison, y=mean)) +
  geom_bar(position=position_dodge(), size=2, stat="identity", color="black", fill="white") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  coord_cartesian(ylim=c(0,80)) +
  facet_grid(cols=vars(year), rows=vars(ws_label))
#export at 1400x1200



###small mammal effect
smallMammalBiomassResponse <- biomass%>%
  filter(trt %in% c('XSI', 'XXI', 'XXX', 'XSX')) %>%
  mutate(comparison=ifelse(trt %in% c('XSI', 'XXI'), 'with inverts', 'without inverts'))

smallMammalBiomassResponse2 <- smallMammalBiomassResponse%>%
  filter(watershed=='N1A' & year<2021) %>%
  group_by(small_mammal) %>%
  summarise(mean=mean(total)) %>%
  ungroup() 
smallMammalBiomassResponse3 <- smallMammalBiomassResponse%>%
  filter(watershed=='N1A' & year==2021) %>%
  group_by(small_mammal) %>%
  summarise(mean=mean(total)) %>%
  ungroup() 
smallMammalBiomassResponse4 <- smallMammalBiomassResponse%>%
  filter(watershed=='N4B') %>%
  group_by(small_mammal) %>%
  summarise(mean=mean(total)) %>%
  ungroup() 
#small mammal consumption: 17% of biomass in annual and no effect in 4yr

#interaction
summary(smallMammalBiomassModel <- lme(total~watershed*small_mammal*invertebrates*year,
                                data=smallMammalBiomassResponse,
                                random=~1|block/trt,
                                correlation=corCompSymm(form=~year|block/trt), 
                                control=lmeControl(returnObject=T)))
anova.lme(smallMammalBiomassModel, type='sequential') 
emmeans(smallMammalBiomassModel, pairwise~watershed*small_mammal*year, adjust="tukey") #no interaction - can plot just small mammal effect

ggplot(data=barGraphStats(data=smallMammalBiomassResponse, variable="total", byFactorNames=c("small_mammal", "watershed", "year")), aes(x=small_mammal, y=mean)) +
  geom_bar(position=position_dodge(), size=2, stat="identity", color="black", fill="white") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  coord_cartesian(ylim=c(0,90)) +
  facet_grid(cols=vars(year), rows=vars(watershed))
#export at 1400x1200




###invertebrate effect
invertebrateBiomassResponse <- biomass%>%
  filter(trt %in% c('XSX', 'XSI', 'BSX', 'BSI')) %>%
  mutate(comparison=ifelse(trt %in% c('BSX', 'BSI'), 'with bison', 'without bison'))

invertebrateBiomassResponse2 <- invertebrateBiomassResponse%>%
  filter(comparison=='without bison') %>%
  group_by(invertebrates) %>%
  summarise(mean=mean(total)) %>%
  ungroup() 
#invert consumption: -4% of biomass in absence of bison
invertebrateBiomassResponse3 <- invertebrateBiomassResponse%>%
  filter(comparison=='with bison') %>%
  group_by(invertebrates) %>%
  summarise(mean=mean(total)) %>%
  ungroup() 
#invert consumption: 9% of biomass in presence of bison

ggplot(data=barGraphStats(data=invertebrateBiomassResponse, variable="total", byFactorNames=c("invertebrates")), aes(x=invertebrates, y=mean, fill=invertebrates)) +
  geom_bar(position=position_dodge(), size=2, stat="identity", color="black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30))

ggplot(data=barGraphStats(data=invertebrateBiomassResponse, variable="gram", byFactorNames=c("invertebrates", "watershed")), aes(x=invertebrates, y=mean, fill=invertebrates)) +
  geom_bar(position=position_dodge(), size=2, stat="identity", color="black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  (facet_wrap(~watershed))

ggplot(data=barGraphStats(data=invertebrateBiomassResponse, variable="forb", byFactorNames=c("invertebrates", "watershed")), aes(x=invertebrates, y=mean, fill=invertebrates)) +
  geom_bar(position=position_dodge(), size=2, stat="identity", color="black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  facet_wrap(~watershed)


#interaction
summary(invertebrateBiomassModel <- lme(total~watershed*bison*invertebrates*year,
                                       data=invertebrateBiomassResponse,
                                       random=~1|block/trt,
                                       correlation=corCompSymm(form=~year|block/trt), 
                                       control=lmeControl(returnObject=T)))
anova.lme(invertebrateBiomassModel, type='sequential') 
emmeans(invertebrateBiomassModel, pairwise~watershed*invertebrates*year, adjust="tukey") #bison interaction

ggplot(data=barGraphStats(data=invertebrateBiomassResponse, variable="total", byFactorNames=c("invertebrates", "bison", "watershed", "year")), aes(x=bison, y=mean, fill=invertebrates)) +
  geom_bar(position=position_dodge(), size=2, stat="identity", color="black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  coord_cartesian(ylim=c(0,90)) +
  facet_grid(cols=vars(year), rows=vars(watershed))
#export at 1400x1200

#with bison
summary(invertebrateBiomassModel <- lme(total~watershed*invertebrates*year,
                                        data=subset(invertebrateBiomassResponse, comparison=='with bison'),
                                        random=~1|block/trt,
                                        correlation=corCompSymm(form=~year|block/trt), 
                                        control=lmeControl(returnObject=T)))
anova.lme(invertebrateBiomassModel, type='sequential') 
emmeans(invertebrateBiomassModel, pairwise~invertebrates, adjust="tukey") #no interaction - can plot just invertebrate effect

ggplot(data=barGraphStats(data=subset(invertebrateBiomassResponse, comparison=='with bison'), variable="total", byFactorNames=c("invertebrates")), aes(x=invertebrates, y=mean)) +
  geom_bar(position=position_dodge(), size=2, stat="identity", color="black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30))
#export at 1400x1200

#without bison
summary(invertebrateBiomassModel <- lme(total~watershed*invertebrates*year,
                                        data=subset(invertebrateBiomassResponse, comparison=='without bison'),
                                        random=~1|block/trt,
                                        correlation=corCompSymm(form=~year|block/trt), 
                                        control=lmeControl(returnObject=T)))
anova.lme(invertebrateBiomassModel, type='sequential') 
emmeans(invertebrateBiomassModel, pairwise~invertebrates, adjust="tukey") #watershed*year*invert effect

ggplot(data=barGraphStats(data=subset(invertebrateBiomassResponse, comparison=='without bison'), variable="total", byFactorNames=c("invertebrates", "watershed", "year")), aes(x=invertebrates, y=mean)) +
  geom_bar(position=position_dodge(), size=2, stat="identity", color="black", fill='white') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  coord_cartesian(ylim=c(0,90)) +
  facet_grid(cols=vars(year), rows=vars(watershed))
#export at 1400x1200




##### ANOVA - graminoid biomass #####
summary(gramBiomassModel <- lme(gram~watershed*year*invertebrates*bison + watershed*year*invertebrates*small_mammal,
                               data=biomass,
                               random=~1|block/trt,
                               correlation=corCompSymm(form=~year|block/trt), 
                               control=lmeControl(returnObject=T)))
anova.lme(gramBiomassModel, type='sequential') 
emmeans(gramBiomassModel, pairwise~watershed*bison*year, adjust="tukey")

#figure - graminoid biomass by bison
ggplot(data=barGraphStats(data=subset(biomass, trt %in% c('XSI', 'BSI')), variable="gram", byFactorNames=c("trt", "watershed", "year")), aes(x=trt, y=mean)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black', fill='white') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Graminoid Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  facet_grid(rows=vars(watershed), cols=vars(year))
#export at 1400x1200

#figure - graminoid biomass by bison
ggplot(data=barGraphStats(data=subset(biomass, trt %in% c('XSI', 'BSI')), variable="gram", byFactorNames=c("bison")), aes(x=bison, y=mean)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black', fill='white') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Graminoid Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  coord_cartesian(ylim=c(0,45))
#export at 1400x1200

#figure - graminoid biomass by small mammal
ggplot(data=barGraphStats(data=subset(biomass, trt %in% c('XXI', 'XSI')), variable="gram", byFactorNames=c("trt", "watershed", "year")), aes(x=trt, y=mean)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black', fill='white') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Graminoid Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  # coord_cartesian(ylim=c(0,700)) +
  facet_grid(rows=vars(watershed), cols=vars(year))
#export at 1400x1200

#figure - graminoid biomass by small mammal
ggplot(data=barGraphStats(data=subset(biomass, trt %in% c('XXI', 'XSI')), variable="gram", byFactorNames=c("trt")), aes(x=trt, y=mean)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black', fill='white') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Graminoid Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  coord_cartesian(ylim=c(0,45))
#export at 1400x1200


##### ANOVA - forb biomass #####
summary(forbBiomassModel <- lme(forb~watershed*year*invertebrates*bison + watershed*year*invertebrates*small_mammal,
                                data=biomass,
                                random=~1|block/trt,
                                correlation=corCompSymm(form=~year|block/trt), 
                                control=lmeControl(returnObject=T)))
anova.lme(forbBiomassModel, type='sequential') 
emmeans(forbBiomassModel, pairwise~trt, adjust="tukey")

invertebrateBiomassResponse <- biomass%>%
  mutate(comparison=ifelse(trt %in% c('BSX', 'BSI'), 'with bison', 'without bison'))

invertebrateBiomassResponse2 <- invertebrateBiomassResponse%>%
  filter(comparison=='without bison') %>%
  group_by(invertebrates) %>%
  summarise(mean=mean(forb)) %>%
  ungroup() 

#without bison
summary(forbBiomassModel <- lme(forb~watershed*year*invertebrates*small_mammal,
                                data=subset(biomass, !(trt %in% c('BSX', 'BSI'))),
                                random=~1|block/trt,
                                correlation=corCompSymm(form=~year|block/trt), 
                                control=lmeControl(returnObject=T)))
anova.lme(forbBiomassModel, type='sequential') 
emmeans(forbBiomassModel, pairwise~watershed*year*invertebrates*small_mammal, adjust="tukey")

#figure - forb biomass with invert and small mammal
ggplot(data=barGraphStats(data=subset(biomass, !(trt %in% c('BSX', 'BSI'))), variable="forb", byFactorNames=c("trt")), aes(x=trt, y=mean)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black', fill='white') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Forb Biomass (g m'^'-2',')'))) +
  scale_x_discrete(limits=c('XSI', 'XXI', 'XSX', 'XXX'), labels=c('SI', 'I', 'S', 'X')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  coord_cartesian(ylim=c(0,22))
#export at 700x600

#figure - forb biomass with invert without bison
ggplot(data=barGraphStats(data=subset(biomass, !(trt %in% c('BSX', 'BSI'))), variable="forb", byFactorNames=c("invertebrates")), aes(x=invertebrates, y=mean)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black', fill='white') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Forb Biomass (g m'^'-2',')'))) +
  # scale_x_discrete(limits=c('BSI', 'BSX', 'XSI', 'XXI', 'XSX', 'XXX')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  coord_cartesian(ylim=c(0,18))
#export at 500x700

#figure - forb biomass with small mammal without bison
ggplot(data=barGraphStats(data=subset(biomass, !(trt %in% c('BSX', 'BSI'))), variable="forb", byFactorNames=c("small_mammal")), aes(x=small_mammal, y=mean)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black', fill='white') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Forb Biomass (g m'^'-2',')'))) +
  # scale_x_discrete(limits=c('BSI', 'BSX', 'XSI', 'XXI', 'XSX', 'XXX')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  coord_cartesian(ylim=c(0,18))
#export at 500x700


##### ANOVA - woody biomass #####
summary(woodyBiomassModel <- lme(woody~watershed*year*invertebrates*bison + watershed*year*invertebrates*small_mammal,
                                 data=biomassLater,
                                 random=~1|block/trt,
                                 correlation=corCompSymm(form=~year|block/trt), 
                                 control=lmeControl(returnObject=T)))
anova.lme(woodyBiomassModel, type='sequential') 
emmeans(woodyBiomassModel, pairwise~trt*year, adjust="tukey")


##### ANOVA - pdead biomass #####
#N4B pdead biomass
summary(pdeadBiomassN4BModel <- lme(pdead~year*invertebrates*bison + year*invertebrates*small_mammal,
                                  data=subset(biomass, watershed=='N4B' & year>2019), #burned in 2019
                                  random=~1|block/trt,
                                  correlation=corCompSymm(form=~year|block/trt), 
                                  control=lmeControl(returnObject=T)))
anova.lme(pdeadBiomassN4BModel, type='sequential') 
emmeans(pdeadBiomassN4BModel, pairwise~trt*year, adjust="tukey")

#figure - pdead biomass all years
ggplot(data=barGraphStats(data=subset(biomass, watershed=='N4B' & year>2019), variable="pdead", byFactorNames=c("trt", "year")), aes(x=trt, y=mean)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Litter Biomass (g m'^'-2',')'))) +
  # scale_x_discrete(limits=c('XSI', 'XXI', 'XSX', 'XXX')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  facet_grid(cols=vars(year))
#export at 1400x600

#figure - pdead biomass small mammal
ggplot(data=barGraphStats(data=subset(biomass, watershed=='N4B' & year>2019 & !(trt %in% c('BSI', 'BSX'))), variable="pdead", byFactorNames=c("trt", "year")), aes(x=trt, y=mean)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Litter Biomass (g m'^'-2',')'))) +
  # scale_x_discrete(limits=c('XSI', 'XXI', 'XSX', 'XXX')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  facet_grid(cols=vars(year))
#export at 1400x600


# ##### ANCOVA - total biomass coupled with graminoid biomass? #####
# summary(biomassRegressions <- lme(total~trt*gram,
#                                   data=biomass,
#                                   random=~1|block/trt))
# anova.lme(biomassRegressions)
# 
# #figure - regression
# ggplot(data=biomass, aes(x=gram, y=total, color=trt)) +
#   geom_point(size=5) +
#   geom_smooth(method='lm', se=F) +
#   ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
#   theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
#   facet_grid(cols=vars(year), rows=vars(watershed))
# #export at 1400x600
# 
# #with bison only
# summary(biomassRegressions <- lme(total~trt*gram,
#                                   data=subset(biomass, trt %in% c('BSI', 'BSX')),
#                                   random=~1|block/trt))
# anova.lme(biomassRegressions)
# 
# #figure - regression
# ggplot(data=subset(biomass), aes(x=gram, y=total, color=invertebrates)) +
#   geom_point(size=5) +
#   geom_smooth(method='lm', se=F) +
#   ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
#   theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
#   facet_grid(cols=vars(year), rows=vars(watershed))
# #export at 1400x600
