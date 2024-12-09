################################################################################
##  conSME_biomass.R: Analysis of aboveground biomass responses in conSME experiment.
##
##  Author: Kimberly Komatsu
##  Date created: December 8, 2021
################################################################################

library(PerformanceAnalytics)
library(nlme)
library(emmeans)
library(gridExtra)
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
biomass2022 <- read.csv('biomass\\conSME_biomass_2022_corrected.csv') %>%
  mutate(experiment='conSME') %>%
  rename(gram=grass) %>%
  filter(!is.na(strip)) %>% 
  mutate(year=2022) %>% 
  select(year, experiment, watershed, block, plot, strip, gram, forb, woody, pdead, notes)
biomass2023 <- read.csv('biomass\\conSME_biomass_2023.csv') %>% 
  mutate(experiment='conSME') %>% 
  select(year, experiment, watershed, block, plot, strip, gram, forb, woody, pdead, notes)
  
biomass <- rbind(biomass2019, biomass2020, biomass2021, biomass2022, biomass2023) %>%
  rename(project_name=experiment) %>%
  left_join(trt) %>%
  select(-notes) %>%
  mutate_all(~replace(., is.na(.), 0)) %>% #don't be scared by the errors, they just don't add NAs to the factor columns (watershed, trt, etc)
  filter(year>0) %>% 
  mutate(gram=10*gram, forb=10*forb, woody=10*woody) %>% 
  mutate(total=gram+forb+woody) %>%
  ungroup() %>% 
  mutate(ws_label=ifelse(watershed=='N1A', 'Annual', '4 Year'),
         experiment_year=year-2018) %>% 
  filter(total>0) #removes two missing plots of data from 2022 :(

biomassMean <- biomass %>%
  group_by(watershed, ws_label, block, plot, year, experiment_year, bison, small_mammal, invertebrates, trt) %>%
  summarise(gram=mean(gram), forb=mean(forb), woody=mean(woody), 
            pdead=mean(pdead), total=mean(total)) %>% 
  ungroup()

# #subsetting out the first year of trts, which is different in patterns from all subsequent years
# biomassLater <- biomassMean %>%
#   filter(year!=2019)
  

# ##### checking for outliers #####
# #outliers have been confirmed as true values for 2019-2023
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

summary(totBiomassModel <- lme(total ~ watershed*experiment_year*invertebrates*bison + 
                                       watershed*experiment_year*invertebrates*small_mammal,
                               data=subset(biomassMean, year<2023),
                               random=~1|block/trt,
                               correlation=corCompSymm(form=~year|block/trt), 
                               control=lmeControl(returnObject=T)))
anova.lme(totBiomassModel, type='sequential') 
contrast(emmeans(totBiomassModel, pairwise~as.factor(experiment_year)*bison*invertebrates), "consec", simple = "each", combine = TRUE, adjust = "mvt") #ws*year*bison, year*small mammal, year*invert*bison and ws*year*bison effects, marginal ws*year*invert

#figure - total biomass ws*bison*year
ggplot(data=barGraphStats(data=biomassMean, variable="total", byFactorNames=c("bison", "ws_label", "year")), aes(x=bison, y=mean, fill=bison)) +
  geom_bar(position=position_dodge(), size=2, stat="identity", color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
  scale_fill_manual(values=c('lightgrey','limegreen')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), 
        axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), 
        axis.text.y=element_text(size=26), 
        legend.position='none', strip.text=element_text(size=30)) +
  facet_grid(cols=vars(year), rows=vars(ws_label))
# ggsave(file='C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\conSME\\figures\\2025\\Fig A_biomass_wsYearBison.png', width=12, height=10, units='in', dpi=300, bg='white')


#figure - total biomass ws*small_mammal*year
ggplot(data=barGraphStats(data=biomassMean, variable="total", byFactorNames=c("small_mammal", "ws_label", "year")), aes(x=small_mammal, y=mean, fill=small_mammal)) +
  geom_bar(position=position_dodge(), size=2, stat="identity", color="black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
  scale_fill_manual(values=c('lightgrey','orange')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), 
        axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), 
        axis.text.y=element_text(size=26), legend.position='none', 
        strip.text=element_text(size=30)) +
  # coord_cartesian(ylim=c(0,90)) +
  facet_grid(cols=vars(year), rows=vars(ws_label))
# ggsave(file='C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\conSME\\figures\\2025\\Fig B_biomass_wsYearSmallMammal.png', width=12, height=10, units='in', dpi=300, bg='white')

#figure - total biomass year*bison*invert
ggplot(data=barGraphStats(data=biomassMean, variable="total", 
                          byFactorNames=c("invertebrates", "bison", "year")), 
       aes(x=interaction(bison,invertebrates), y=mean,  
           fill=interaction(bison, invertebrates))) +
  geom_bar(position=position_dodge(), size=2, stat="identity", color="black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
  scale_fill_manual(values=c('lightgrey','limegreen','azure4','darkgreen')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), 
        axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), 
        axis.text.y=element_text(size=26), 
        legend.position='none', strip.text=element_text(size=30)) +
  # coord_cartesian(ylim=c(0,90)) +
  facet_grid(cols=vars(year)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# ggsave(file='C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\conSME\\figures\\2025\\Fig Ca_biomass_yearBisonInvert.png', width=12, height=6, units='in', dpi=300, bg='white')

#figure - total biomass year*watershed*bison*invert
ggplot(data=barGraphStats(data=biomassMean, variable="total", 
                          byFactorNames=c("invertebrates", "bison", "year", 'ws_label')), 
       aes(x=interaction(bison,invertebrates), y=mean,  
           fill=interaction(bison, invertebrates))) +
  geom_bar(position=position_dodge(), size=2, stat="identity", color="black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
  scale_fill_manual(values=c('lightgrey','limegreen','azure4','darkgreen')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), 
        axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), 
        axis.text.y=element_text(size=26), 
        legend.position='none', strip.text=element_text(size=30)) +
  # coord_cartesian(ylim=c(0,90)) +
  facet_grid(cols=vars(year), rows=vars(ws_label)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# ggsave(file='C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\conSME\\figures\\2025\\Fig Cb_biomass_yearWatershedBisonInvert.png', width=12, height=10, units='in', dpi=300, bg='white')



##### total biomass - response ratios #####
###bison effect
bisonBiomassResponse <- biomassMean %>%
  filter(trt %in% c('XSI', 'BSI', 'XSX', 'BSX')) %>%
  mutate(comparison=ifelse(trt %in% c('XSI', 'BSI'), 'with inverts', 'without inverts'))

bisonBiomassResponse2 <- bisonBiomassResponse %>%
  filter(watershed=='N1A') %>%
  group_by(bison) %>%
  summarise(mean=mean(total)) %>%
  ungroup() 
bisonBiomassResponse3 <- bisonBiomassResponse %>%
  filter(watershed=='N4B') %>%
  group_by(bison) %>%
  summarise(mean=mean(total)) %>%
  ungroup() 
#bison consumption: 26% of biomass in annual and 44% of biomass in 4yr

###small mammal effect
smallMammalBiomassResponse <- biomassMean %>%
  filter(trt %in% c('XSI', 'XXI', 'XXX', 'XSX')) %>%
  mutate(comparison=ifelse(trt %in% c('XSI', 'XXI'), 'with inverts', 'without inverts'))

smallMammalBiomassResponse2 <- smallMammalBiomassResponse %>%
  filter(watershed=='N1A' & year>2018) %>%
  group_by(small_mammal) %>%
  summarise(mean=mean(total)) %>%
  ungroup() 
#8% diff in biomass
smallMammalBiomassResponse3 <- smallMammalBiomassResponse %>%
  filter(watershed=='N1A' & year %in% c(2019, 2020, 2022)) %>%
  group_by(small_mammal) %>%
  summarise(mean=mean(total)) %>%
  ungroup() 
smallMammalBiomassResponse4 <- smallMammalBiomassResponse%>%
  filter(watershed=='N4B' & year %in% c(2019, 2020, 2023)) %>%
  group_by(small_mammal) %>%
  summarise(mean=mean(total)) %>%
  ungroup() 
#small mammal consumption: 14% of biomass in annual and 2% in 4yr IN YEARS WHEN THEY HAVE AN EFFECT



##### invertebrate effect --  only look at plots with small mammals present! (drops plots where small mammals have been excluded) #####
invertebrateBiomassResponse <- biomassMean %>%
  filter(trt %in% c('XSX', 'XSI', 'BSX', 'BSI')) %>%
  mutate(comparison=ifelse(trt %in% c('BSX', 'BSI'), 'with bison', 'without bison'))

invertebrateBiomassResponse2 <- invertebrateBiomassResponse %>%
  filter(comparison=='without bison' & year>2018) %>%
  group_by(invertebrates, watershed) %>%
  summarise(mean=mean(total)) %>%
  ungroup() 
#invert consumption: 3% of biomass in absence of bison in annual burn, 1% in 4 yr
invertebrateBiomassResponse3 <- invertebrateBiomassResponse %>%
  filter(comparison=='with bison' & year>2018) %>%
  group_by(invertebrates, watershed) %>%
  summarise(mean=mean(total)) %>%
  ungroup() 
#invert consumption: -8% of biomass in presence of bison in annual burn, -1% in 4 yr

# ggplot(data=barGraphStats(data=invertebrateBiomassResponse, variable="total", byFactorNames=c("invertebrates", "watershed", "year")), aes(x=invertebrates, y=mean, fill=invertebrates)) +
#   geom_bar(position=position_dodge(), size=2, stat="identity", color="black") +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
#   ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
#   scale_fill_manual(values=c('lightgrey', 'lightgoldenrod')) +
#   theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position='none', legend.justification=c(0,1), strip.text=element_text(size=30)) +
#   facet_grid(cols=vars(year), rows=vars(watershed))
# 
# ggplot(data=barGraphStats(data=invertebrateBiomassResponse, variable="gram", byFactorNames=c("invertebrates", "watershed")), aes(x=invertebrates, y=mean, fill=invertebrates)) +
#   geom_bar(position=position_dodge(), size=2, stat="identity", color="black") +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
#   ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
#   theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
#   (facet_wrap(~watershed))
# 
# ggplot(data=barGraphStats(data=invertebrateBiomassResponse, variable="forb", byFactorNames=c("invertebrates", "watershed")), aes(x=invertebrates, y=mean, fill=invertebrates)) +
#   geom_bar(position=position_dodge(), size=2, stat="identity", color="black") +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
#   ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
#   theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
#   facet_wrap(~watershed)
#  
# #interaction
# summary(invertebrateBiomassModel <- lme(total~watershed*as.factor(year)*bison*invertebrates,
#                                        data=invertebrateBiomassResponse,
#                                        random=~1|block/trt,
#                                        correlation=corCompSymm(form=~year|block/trt), 
#                                        control=lmeControl(returnObject=T)))
# anova.lme(invertebrateBiomassModel, type='sequential') 
# contrast(emmeans(invertebrateBiomassModel, pairwise~watershed*as.factor(year)*bison*invertebrates), "consec", simple = "each", combine = TRUE, adjust = "mvt") 
# 
# ggplot(data=barGraphStats(data=invertebrateBiomassResponse, variable="total", byFactorNames=c("invertebrates", "watershed", "year")), aes(x=invertebrates, y=mean, fill=invertebrates)) +
#   geom_bar(position=position_dodge(), size=2, stat="identity", color="black") +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
#   ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
#   scale_fill_manual(values=c('lightgrey','skyblue')) +
#   theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), 
#         axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), 
#         axis.text.y=element_text(size=26), legend.position='none', 
#         strip.text=element_text(size=30)) +
#   # coord_cartesian(ylim=c(0,90)) +
#   facet_grid(cols=vars(year), rows=vars(watershed))
# # ggsave(file='C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\conSME\\figures\\2025\\Fig D_biomass_yearWatershedInvert.png', width=12, height=10, units='in', dpi=300, bg='white')

# #with bison
# summary(invertebrateBiomassModel <- lme(total~watershed*invertebrates*year,
#                                         data=subset(invertebrateBiomassResponse, comparison=='with bison'),
#                                         random=~1|block/trt,
#                                         correlation=corCompSymm(form=~year|block/trt), 
#                                         control=lmeControl(returnObject=T)))
# anova.lme(invertebrateBiomassModel, type='sequential') 
# emmeans(invertebrateBiomassModel, pairwise~invertebrates, adjust="tukey") #no interaction - can plot just invertebrate effect
# 
# ggplot(data=barGraphStats(data=subset(invertebrateBiomassResponse, comparison=='with bison'), variable="total", byFactorNames=c("invertebrates")), aes(x=invertebrates, y=mean)) +
#   geom_bar(position=position_dodge(), size=2, stat="identity", color="black") +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
#   ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
#   theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30))
# #export at 1400x1200
# 
# #without bison
# summary(invertebrateBiomassModel <- lme(total~watershed*invertebrates*year,
#                                         data=subset(invertebrateBiomassResponse, 
#                                                     comparison=='without bison'),
#                                         random=~1|block/trt,
#                                         correlation=corCompSymm(form=~year|block/trt), 
#                                         control=lmeControl(returnObject=T)))
# anova.lme(invertebrateBiomassModel, type='sequential') 
# emmeans(invertebrateBiomassModel, pairwise~invertebrates, adjust="tukey") #watershed*year*invert effect
# 
# ggplot(data=barGraphStats(data=subset(invertebrateBiomassResponse, comparison=='without bison'), variable="total", byFactorNames=c("invertebrates")), aes(x=invertebrates, y=mean, fill=invertebrates)) +
#   geom_bar(position=position_dodge(), size=2, stat="identity", color="black") +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
#   ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
#   scale_fill_manual(values=c('lightgrey', 'lightgoldenrod')) +
#   theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position='none', legend.justification=c(0,1), strip.text=element_text(size=30)) +
#   annotate("text", x=1, y=550, label='a', size=9) +
#   annotate("text", x=2, y=590, label='b', size=9)
# # ggsave('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\konza projects\\conSME\\figures\\2023\\invert_total.png', width=5, height=5, units='in', dpi=600, bg='white')

# ggplot(data=barGraphStats(data=subset(invertebrateBiomassResponse), variable="total", byFactorNames=c("invertebrates")), aes(x=invertebrates, y=mean)) +
#   geom_bar(position=position_dodge(), size=2, stat="identity", color="black", fill='white') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
#   ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
#   theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30))
# #export at 1400x1200
# 
# ggplot(data=barGraphStats(data=subset(invertebrateBiomassResponse, comparison=='with bison'), variable="total", byFactorNames=c("invertebrates")), aes(x=invertebrates, y=mean)) +
#   geom_bar(position=position_dodge(), size=2, stat="identity", color="black", fill='white') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
#   ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
#   theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30))


##### ANOVA - graminoid biomass #####
summary(gramBiomassModel <- lme(gram~watershed*year*invertebrates*bison + 
                                watershed*year*invertebrates*small_mammal,
                               data=biomassMean,
                               random=~1|block/trt,
                               correlation=corCompSymm(form=~year|block/trt), 
                               control=lmeControl(returnObject=T)))
anova.lme(gramBiomassModel, type='sequential') 
contrast(emmeans(gramBiomassModel, pairwise~watershed*as.factor(year)*bison), "consec", simple = "each", combine = TRUE, adjust = "mvt")

#figure - total biomass ws*bison*year
ggplot(data=barGraphStats(data=biomassMean, variable="gram", byFactorNames=c("bison", "ws_label", "year")), aes(x=bison, y=mean, fill=bison)) +
  geom_bar(position=position_dodge(), size=2, stat="identity", color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Graminoid Biomass (g m'^'-2',')'))) +
  scale_fill_manual(values=c('lightgrey','limegreen')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), 
        axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), 
        axis.text.y=element_text(size=26), 
        legend.position='none', strip.text=element_text(size=30)) +
  facet_grid(cols=vars(year), rows=vars(ws_label))
# ggsave(file='C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\conSME\\figures\\2025\\Fig D_graminoid_wsYearBison.png', width=12, height=10, units='in', dpi=300, bg='white')

#figure - graminoid biomass by small mammal
ggplot(data=barGraphStats(data=biomassMean, variable="gram", byFactorNames=c("small_mammal", "year")), aes(x=small_mammal, y=mean, fill=small_mammal)) +
  geom_bar(position=position_dodge(), size=2, stat="identity", color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Graminoid Biomass (g m'^'-2',')'))) +
  scale_fill_manual(values=c('lightgrey','orange')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), 
        axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), 
        axis.text.y=element_text(size=26), 
        legend.position='none', strip.text=element_text(size=30)) +
  facet_grid(cols=vars(year))
# ggsave(file='C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\conSME\\figures\\2025\\Fig Ea_graminoid_yearSmallMammal.png', width=12, height=6, units='in', dpi=300, bg='white')

# ggplot(data=barGraphStats(data=biomassMean, variable="gram", byFactorNames=c("small_mammal", "ws_label", "year")), aes(x=small_mammal, y=mean, fill=small_mammal)) +
#   geom_bar(position=position_dodge(), size=2, stat="identity", color='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
#   ylab(expression(paste('Graminoid Biomass (g m'^'-2',')'))) +
#   scale_fill_manual(values=c('lightgrey','orange')) +
#   theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), 
#         axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), 
#         axis.text.y=element_text(size=26), 
#         legend.position='none', strip.text=element_text(size=30)) +
#   facet_grid(cols=vars(year), rows=vars(ws_label))
# # ggsave(file='C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\conSME\\figures\\2025\\Fig Eb_graminoid_watershedYearSmallMammal.png', width=12, height=10, units='in', dpi=300, bg='white')


##### ANOVA - forb biomass #####
summary(forbBiomassModel <- lme(forb~watershed*year*invertebrates*bison + 
                                watershed*year*invertebrates*small_mammal,
                                data=biomassMean,
                                random=~1|block/trt,
                                correlation=corCompSymm(form=~year|block/trt), 
                                control=lmeControl(returnObject=T)))
anova.lme(forbBiomassModel, type='sequential') 
emmeans(forbBiomassModel, pairwise~trt, adjust="tukey")

bisonBiomassResponse2 <- biomassMean %>%
  group_by(bison, watershed) %>%
  summarise(mean=mean(forb)) %>%
  ungroup() 
#bison reduce forb biomass by 21% overall, 25% in annual and 9% in 4-year

invertebrateBiomassResponse <- biomassMean %>%
  mutate(comparison=ifelse(trt %in% c('BSX', 'BSI'), 'with bison', 'without bison'))

invertebrateBiomassResponse2 <- invertebrateBiomassResponse %>%
  # filter(comparison=='without bison') %>%
  group_by(invertebrates, comparison) %>%
  summarise(mean=mean(forb)) %>%
  ungroup() 
#invertebrates remove 29% of forb biomass in the absence of bison (10% higher with invert removal in presence of bison)

smallMammalBiomassResponse2 <- invertebrateBiomassResponse %>%
  # filter(comparison=='without bison') %>%
  group_by(small_mammal, comparison) %>%
  summarise(mean=mean(forb)) %>%
  ungroup() 
#invertebrates remove 22% of forb biomass in the absence of bison (can't say about with bison present, due to fencing nestedness)

#without bison
summary(forbBiomassModel <- lme(forb~watershed*year*invertebrates*small_mammal,
                                data=subset(biomassMean, !(trt %in% c('BSX', 'BSI'))),
                                random=~1|block/trt,
                                correlation=corCompSymm(form=~year|block/trt), 
                                control=lmeControl(returnObject=T)))
anova.lme(forbBiomassModel, type='sequential') 
emmeans(forbBiomassModel, pairwise~invertebrates*small_mammal, adjust="tukey")

#figure - forb biomass with invert and small mammal
ggplot(data=barGraphStats(data=subset(biomassMean, !(trt %in% c('BSX', 'BSI'))), variable="forb", byFactorNames=c("trt")), aes(x=trt, y=mean)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black', fill='white') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Forb Biomass (g m'^'-2',')'))) +
  scale_x_discrete(limits=c('XSI', 'XXI', 'XSX', 'XXX'), labels=c('SI', 'I', 'S', 'X')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  annotate("text", x=1, y=100, label='a', size=9) +
  annotate("text", x=2, y=130, label='ab', size=9) +
  annotate("text", x=3, y=140, label='ab', size=9) +
  annotate("text", x=4, y=170, label='b', size=9) 
  # coord_cartesian(ylim=c(0,22))
# ggsave('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\konza projects\\conSME\\figures\\2023\\forb_invertSmallMammal.png', width=5, height=5, units='in', dpi=600, bg='white')


#forb biomass -- bison effect
forbBiomassFigA <- ggplot(data=barGraphStats(data=subset(biomassMean), variable="forb", byFactorNames=c("bison")), aes(x=bison, y=mean, fill=bison)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Forb Biomass (g m'^'-2',')'))) +
  scale_fill_manual(values=c('lightgrey','limegreen')) +
  coord_cartesian(ylim=c(0,150)) +
  # scale_x_discrete(limits=c('BSI', 'BSX', 'XSI', 'XXI', 'XSX', 'XXX')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), 
        axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), 
        axis.text.y=element_text(size=26), legend.position='none', 
        strip.text=element_text(size=30)) 

forbBiomassFigAa <- ggplot(data=barGraphStats(data=subset(biomassMean), variable="forb", byFactorNames=c("bison", "ws_label", "year")), aes(x=bison, y=mean, fill=bison)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Forb Biomass (g m'^'-2',')'))) +
  scale_fill_manual(values=c('lightgrey','limegreen')) +
  # coord_cartesian(ylim=c(0,150)) +
  # scale_x_discrete(limits=c('BSI', 'BSX', 'XSI', 'XXI', 'XSX', 'XXX')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), 
        axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), 
        axis.text.y=element_text(size=26), legend.position='none', 
        strip.text=element_text(size=30)) +
  facet_grid(cols=vars(year),rows=vars(ws_label))

#figure - forb biomass with invert without bison
forbBiomassFigB <- ggplot(data=barGraphStats(data=subset(biomassMean, !(trt %in% c('BSX', 'BSI'))), variable="forb", byFactorNames=c("invertebrates")), aes(x=invertebrates, y=mean, fill=invertebrates)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab('') +
  scale_fill_manual(values=c('limegreen','skyblue')) +
  coord_cartesian(ylim=c(0,150)) +
  # scale_x_discrete(limits=c('BSI', 'BSX', 'XSI', 'XXI', 'XSX', 'XXX')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), 
        axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), 
        axis.text.y=element_text(size=26), legend.position='none', 
        strip.text=element_text(size=30)) 

forbBiomassFigBb <- ggplot(data=barGraphStats(data=subset(biomassMean, !(trt %in% c('BSX', 'BSI'))), variable="forb", byFactorNames=c("invertebrates", "year", "ws_label")), aes(x=invertebrates, y=mean, fill=invertebrates)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  # ylab(expression(paste('Forb Biomass (g m'^'-2',')'))) +
  ylab('') +
  scale_fill_manual(values=c('limegreen','skyblue')) +
  # coord_cartesian(ylim=c(0,150)) +
  # scale_x_discrete(limits=c('BSI', 'BSX', 'XSI', 'XXI', 'XSX', 'XXX')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), 
        axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), 
        axis.text.y=element_text(size=26), legend.position='none', 
        strip.text=element_text(size=30)) +
  facet_grid(cols=vars(year),rows=vars(ws_label))

#figure - forb biomass with small mammal without bison
forbBiomassFigC <- ggplot(data=barGraphStats(data=subset(biomassMean, !(trt %in% c('BSX', 'BSI'))), variable="forb", byFactorNames=c("small_mammal")), aes(x=small_mammal, y=mean, fill=small_mammal)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab('') +
  scale_fill_manual(values=c('limegreen','orange')) +
  coord_cartesian(ylim=c(0,150)) +
  # scale_x_discrete(limits=c('BSI', 'BSX', 'XSI', 'XXI', 'XSX', 'XXX')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), 
        axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), 
        axis.text.y=element_text(size=26), legend.position='none', 
        strip.text=element_text(size=30)) 

forbBiomassFigCc <- ggplot(data=barGraphStats(data=subset(biomassMean, !(trt %in% c('BSX', 'BSI'))), variable="forb", byFactorNames=c("small_mammal", "ws_label", "year")), aes(x=small_mammal, y=mean, fill=small_mammal)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Forb Biomass (g m'^'-2',')'))) +
  scale_fill_manual(values=c('limegreen','orange')) +
  # coord_cartesian(ylim=c(0,150)) +
  # scale_x_discrete(limits=c('BSI', 'BSX', 'XSI', 'XXI', 'XSX', 'XXX')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), 
        axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), 
        axis.text.y=element_text(size=26), legend.position='none', 
        strip.text=element_text(size=30))  +
  facet_grid(cols=vars(year),rows=vars(ws_label))

grid.arrange(forbBiomassFigA, forbBiomassFigC, forbBiomassFigB, ncol=3)
# export to file Fig F_forb_mainEffects.png at 1600x800

grid.arrange(forbBiomassFigAa, forbBiomassFigCc, forbBiomassFigB, nrow=3)
# export to file Fig F_forb_interactiveEffects.png at 1600x1600


##### ANOVA - woody biomass #####
summary(woodyBiomassModel <- lme(woody~watershed*year*invertebrates*bison + 
                                 watershed*year*invertebrates*small_mammal,
                                 data=biomassMean,
                                 random=~1|block/trt,
                                 correlation=corCompSymm(form=~year|block/trt), 
                                 control=lmeControl(returnObject=T)))
anova.lme(woodyBiomassModel, type='sequential') 


##### ANOVA - pdead biomass #####
#N4B pdead biomass
summary(pdeadBiomassN4BModel <- lme(pdead~year*invertebrates*bison + 
                                    year*invertebrates*small_mammal,
                                  data=subset(biomassMean, watershed=='N4B' & 
                                              year>2019 & year<2023), #burned in 2019 and 2023
                                  random=~1|block/trt,
                                  correlation=corCompSymm(form=~year|block/trt), 
                                  control=lmeControl(returnObject=T)))
anova.lme(pdeadBiomassN4BModel, type='sequential') 