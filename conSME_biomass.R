################################################################################
##  conSME_biomass.R: Analysis of aboveground biomass responses in conSME experiment.
##
##  Author: Kimberly Komatsu
##  Date created: December 8, 2021
##  Modified by: Allison Louthan January 22, 2025
################################################################################

library(PerformanceAnalytics)
library(nlme)
library(emmeans)
library(gridExtra)
library(cowplot)
library(tidyverse)

user <- "AL" # change based on your initials to deal with directory issues


if (user== "KK"){setwd('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\konza projects\\conSME\\data') }#desktop


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
trt <- read.csv('data/conSME_treatments.csv')
if (user== "KK") {
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
}

if (user== "AL") {
  biomass2019 <- read.csv('data/biomass/conSME_biomass_2019.csv') %>%
    mutate(pdead=0)
  biomass2020 <- read.csv('data/biomass/conSME_biomass_2020.csv') %>%
    mutate(drop=ifelse(block=='I'&plot==6&strip==1, 1, 0)) %>% #drop missing sample
    filter(drop!=1) %>%
    select(-drop)
  biomass2021 <- read.csv('data/biomass/conSME_biomass_2021.csv') %>%
    mutate(experiment='conSME') %>%
    rename(gram=grass) %>%
    select(-date) %>%
    filter(notes!='no biomass in any bag') #filter out missing sample
  biomass2022 <- read.csv('data/biomass/conSME_biomass_2022_corrected.csv') %>%
    mutate(experiment='conSME') %>%
    rename(gram=grass) %>%
    filter(!is.na(strip)) %>% 
    mutate(year=2022) %>% 
    select(year, experiment, watershed, block, plot, strip, gram, forb, woody, pdead, notes)
  biomass2023 <- read.csv('data/biomass/conSME_biomass_2023.csv') %>% 
    mutate(experiment='conSME') %>% 
    select(year, experiment, watershed, block, plot, strip, gram, forb, woody, pdead, notes)
}
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

if (user== "AL") {write.csv(biomassMean, file='derived_data/01biomassMean.csv')}
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
summary(totBiomassModel <- lme(total ~ watershed*year*invertebrates*bison + 
                                       watershed*year*invertebrates*small_mammal,
                               data=subset(biomassMean, year<2023),
                               random=~1|block/trt,
                               correlation=corCompSymm(form=~year|block/trt), 
                               control=lmeControl(returnObject=T)))
anova.lme(totBiomassModel, type='sequential') 
emmeans(totBiomassModel, ~bison*invertebrates*year*watershed)
contrast(emmeans(totBiomassModel, pairwise~year*watershed*bison*invertebrates), "consec", simple = "each", combine = TRUE, adjust = "mvt") #ws*year*bison, year*small mammal, year*invert*bison and ws*year*bison effects, marginal ws*year*invert

#figure - total biomass ws*bison*year
biomassFig1a <- ggplot(data=barGraphStats(data=biomassMean, variable="total", byFactorNames=c("bison", "ws_label", "year")), aes(x=bison, y=mean, fill=bison)) +
  geom_bar(position=position_dodge(), size=2, stat="identity", color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
  scale_fill_manual(values=c('lightgrey','limegreen')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), 
        axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), 
        axis.text.y=element_text(size=26), 
        legend.position='none', strip.text=element_text(size=30)) +
  coord_cartesian(ylim=c(0,800)) +
  facet_grid(cols=vars(year), rows=vars(ws_label))
# if (user= "KK") {ggsave(file='C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\conSME\\figures\\2025\\Fig A_biomass_wsYearBison.png', width=12, height=10, units='in', dpi=300, bg='white')}

#figure - total biomass small_mammal*year
biomassSmallMammalTemp <- biomassMean %>% mutate(ws2='ws')
biomassFig1b <- ggplot(data=barGraphStats(data=biomassSmallMammalTemp, variable="total", byFactorNames=c("small_mammal", "year", "ws2")), aes(x=small_mammal, y=mean, fill=small_mammal)) +
  geom_bar(position=position_dodge(), size=2, stat="identity", color="black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
  scale_fill_manual(values=c('lightgrey','orange')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), 
        axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), 
        axis.text.y=element_text(size=26), legend.position='none', 
        strip.text=element_text(size=30)) +
  coord_cartesian(ylim=c(0,800)) +
  facet_grid(cols=vars(year), rows=vars(ws2))
# if (user= "KK") {ggsave(file='C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\conSME\\figures\\2025\\Fig B_biomass_wsYearSmallMammal.png', width=12, height=10, units='in', dpi=300, bg='white')}

#figure - total biomass year*bison*invert
biomassInvertTemp <- biomassMean %>% mutate(ws2='ws')
biomassFig1c <- ggplot(data=barGraphStats(data=biomassInvertTemp, variable="total", 
                          byFactorNames=c("invertebrates", "bison", "year", 'ws2')), 
       aes(x=interaction(bison,invertebrates), y=mean,  
           fill=interaction(bison, invertebrates))) +
  geom_bar(position=position_dodge(), size=2, stat="identity", color="black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
  scale_fill_manual(values=c('lightgrey','limegreen','skyblue','turquoise4')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), 
        axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), 
        axis.text.y=element_text(size=26), 
        legend.position='none', strip.text=element_text(size=30)) +
  coord_cartesian(ylim=c(0,800)) +
  facet_grid(cols=vars(year), rows=vars(ws2)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# if (user= "KK") {ggsave(file='C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\conSME\\figures\\2025\\Fig Ca_biomass_yearBisonInvert.png', width=12, height=6, units='in', dpi=300, bg='white')}

plot_grid(biomassFig1a, biomassFig1b, biomassFig1c, 
          align = "v", nrow = 3, rel_heights = c(5, 2.8, 3))
# save as Fig 1_biomass_all.png at 1600 x 1600 


#figure - total biomass year*watershed*bison*invert
ggplot(data=barGraphStats(data=biomassMean, variable="total", 
                          byFactorNames=c("invertebrates", "bison", "year", 'ws_label')), 
       aes(x=interaction(bison,invertebrates), y=mean,  
           fill=interaction(bison, invertebrates))) +
  geom_bar(position=position_dodge(), size=2, stat="identity", color="black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
  scale_fill_manual(values=c('lightgrey','limegreen','skyblue','darkgreen')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), 
        axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), 
        axis.text.y=element_text(size=26), 
        legend.position='none', strip.text=element_text(size=30)) +
  # coord_cartesian(ylim=c(0,90)) +
  facet_grid(cols=vars(year), rows=vars(ws_label)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# if (user= "KK") {ggsave(file='C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\conSME\\figures\\2025\\Fig Cb_biomass_yearWatershedBisonInvert.png', width=12, height=10, units='in', dpi=300, bg='white')}



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
  filter(watershed=='N1A' & year %in% c(2019, 2020, 2022)) %>%
  group_by(small_mammal) %>%
  summarise(mean=mean(total)) %>%
  ungroup() 
#5% diff in biomass over all years, 14% of biomass in years with significant response only
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

### invertebrate effect
invertebrateBiomassResponse <- biomassMean %>%
  # filter(trt %in% c('XSX', 'XSI', 'BSX', 'BSI')) %>%
  mutate(comparison=ifelse(trt %in% c('BSX', 'BSI'), 'with bison', 'without bison'))

invertebrateBiomassResponse2 <- invertebrateBiomassResponse %>%
  filter(comparison=='without bison' & year>2018) %>%
  group_by(invertebrates, watershed) %>%
  summarise(mean=mean(total)) %>%
  ungroup() 
#invert consumption: 7% of biomass in absence of bison across both watersheds; 5% in annual burn, 9% in 4 yr
invertebrateBiomassResponse3 <- invertebrateBiomassResponse %>%
  filter(comparison=='with bison' & year>2018) %>%
  group_by(invertebrates, watershed) %>%
  summarise(mean=mean(total)) %>%
  ungroup() 
#invert consumption: -8% of biomass in presence of bison in annual burn, -1% in 4 yr

#with bison
summary(invertebrateBiomassModel <- lme(total~watershed*year*invertebrates,
                                        data=subset(invertebrateBiomassResponse, comparison=='with bison'),
                                        random=~1|block/trt,
                                        correlation=corCompSymm(form=~year|block/trt),
                                        control=lmeControl(returnObject=T)))
anova.lme(invertebrateBiomassModel, type='sequential')
emmeans(invertebrateBiomassModel, pairwise~invertebrates, adjust="tukey") #no effect

#without bison
summary(invertebrateBiomassModel <- lme(total~watershed*small_mammal*invertebrates*year,
                                        data=subset(invertebrateBiomassResponse,
                                                    comparison=='without bison'),
                                        random=~1|block/trt,
                                        correlation=corCompSymm(form=~year|block/trt),
                                        control=lmeControl(returnObject=T)))
anova.lme(invertebrateBiomassModel, type='sequential')
emmeans(invertebrateBiomassModel, pairwise~invertebrates, adjust="tukey") #marignally significant invert effect

ggplot(data=barGraphStats(data=subset(invertebrateBiomassResponse, comparison=='without bison'), variable="total", byFactorNames=c("invertebrates")), aes(x=invertebrates, y=mean, fill=invertebrates)) +
  geom_bar(position=position_dodge(), size=2, stat="identity", color="black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
  scale_fill_manual(values=c('limegreen', 'skyblue')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position='none', legend.justification=c(0,1), strip.text=element_text(size=30)) +
  annotate("text", x=1, y=550, label='a', size=9) +
  annotate("text", x=2, y=590, label='b', size=9)
# ggsave('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\konza projects\\conSME\\figures\\2025\\Fig 2_biomass_invert.png', width=5, height=5, units='in', dpi=600, bg='white')


##### ANOVA - graminoid biomass #####
summary(gramBiomassModel <- lme(gram~watershed*year*invertebrates*bison + 
                                watershed*year*invertebrates*small_mammal,
                               data=biomassMean,
                               random=~1|block/trt,
                               correlation=corCompSymm(form=~year|block/trt), 
                               control=lmeControl(returnObject=T)))
anova.lme(gramBiomassModel, type='sequential') 
contrast(emmeans(gramBiomassModel, pairwise~watershed*year*bison), "consec", simple = "each", combine = TRUE, adjust = "mvt")

#figure - total biomass ws*bison*year
gramFig3a <- ggplot(data=barGraphStats(data=biomassMean, variable="gram", byFactorNames=c("bison", "ws_label", "year")), aes(x=bison, y=mean, fill=bison)) +
  geom_bar(position=position_dodge(), size=2, stat="identity", color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Graminoid\nBiomass (g m'^'-2',')'))) +
  scale_fill_manual(values=c('lightgrey','limegreen')) +
  coord_cartesian(ylim=c(0,500)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), 
        axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), 
        axis.text.y=element_text(size=26), 
        legend.position='none', strip.text=element_text(size=30)) +
  facet_grid(cols=vars(year), rows=vars(ws_label))
# ggsave(file='C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\conSME\\figures\\2025\\Fig D_graminoid_wsYearBison.png', width=12, height=10, units='in', dpi=300, bg='white')

#figure - graminoid biomass by small mammal
smallMammalTemp <- biomassMean %>% mutate(ws2='ws')
gramFig3b <- ggplot(data=barGraphStats(data=smallMammalTemp, variable="gram", byFactorNames=c("small_mammal", "year", "ws2")), aes(x=small_mammal, y=mean, fill=small_mammal)) +
  geom_bar(position=position_dodge(), size=2, stat="identity", color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Graminoid\nBiomass (g m'^'-2',')'))) +
  scale_fill_manual(values=c('lightgrey','orange')) +
  coord_cartesian(ylim=c(0,500)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), 
        axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), 
        axis.text.y=element_text(size=26), 
        legend.position='none', strip.text=element_text(size=30)) +
  facet_grid(cols=vars(year), rows=vars(ws2))
# ggsave(file='C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\conSME\\figures\\2025\\Fig Ea_graminoid_yearSmallMammal.png', width=12, height=6, units='in', dpi=300, bg='white')

plot_grid(gramFig3a, gramFig3b,  nrow = 2, rel_heights = c(5, 2.8))
# save as Fig S1_graminoid_all.png at 1600 x 1067 



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
  # ylab(expression(paste('Forb Biomass (g m'^'-2',')'))) +
  ylab('') +
  scale_fill_manual(values=c('limegreen','skyblue')) +
  coord_cartesian(ylim=c(0,150)) +
  # scale_x_discrete(limits=c('BSI', 'BSX', 'XSI', 'XXI', 'XSX', 'XXX')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), 
        axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), 
        axis.text.y=element_text(size=26), legend.position='none', 
        strip.text=element_text(size=30))

invertTemp <- biomassMean %>% mutate(year2=2020, ws2='ws')
forbBiomassFigBb <- ggplot(data=barGraphStats(data=subset(invertTemp, !(trt %in% c('BSX', 'BSI'))), variable="forb", byFactorNames=c("invertebrates", 'year2', 'ws2')), aes(x=invertebrates, y=mean, fill=invertebrates)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Forb Biomass (g m'^'-2',')'))) +
  # ylab('') +
  scale_fill_manual(values=c('limegreen','skyblue')) +
  coord_cartesian(ylim=c(0,150)) +
  # scale_x_discrete(limits=c('BSI', 'BSX', 'XSI', 'XXI', 'XSX', 'XXX')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), 
        axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), 
        axis.text.y=element_text(size=26), legend.position='none', 
        strip.text=element_text(size=30))  +
  facet_grid(cols=vars(year2),rows=vars(ws2))

forbBiomassFigBb2 <- ggplot(data=barGraphStats(data=subset(biomassMean, !(trt %in% c('BSX', 'BSI'))), variable="forb", byFactorNames=c("invertebrates", "year", "ws_label")), aes(x=invertebrates, y=mean, fill=invertebrates)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Forb Biomass (g m'^'-2',')'))) +
  # ylab('') +
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
# export to file Fig 3_forb_mainEffects.png at 1600x800

plot_grid(forbBiomassFigAa, forbBiomassFigCc, arrangeGrob(forbBiomassFigBb, nullGrob(), ncol=2), 
          nrow=3, rel_widths=c(3,3,1), rel_heights=c(2,2,1))
# export to file Fig 3_forb_interactiveEffects.png at 1600x1600


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

if (user== "AL") {save(totBiomassModel, gramBiomassModel, forbBiomassModel, woodyBiomassModel, 
                       file='derived_data/01biomass_models.csv')}
