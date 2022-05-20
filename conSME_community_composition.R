################################################################################
##  conSME_community_composition.R: Analysis of plant community composition responses in conSME experiment.
##
##  Author: Kimberly Komatsu
##  Date created: December 8, 2021
################################################################################

library(codyn)
library(nlme)
library(emmeans)
library(tidyverse)

setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\konza projects\\conSME\\data\\species composition') #laptop
setwd('C:\\Users\\komatsuk\\Dropbox (Smithsonian)\\konza projects\\conSME\\data\\species composition') #desktop


##### functions #####
`%notin%` <- Negate(`%in%`)

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

sp2018 <- read.csv('ConSME_species composition_2018.csv')
sp2019 <- read.csv('ConSME_species composition_2019.csv')
sp2020 <- read.csv('ConSME_species composition_2020.csv')%>%
  select(-X, -taxa)
sp2021 <- read.csv('ConSME_species composition_2021.csv')%>%
  select(-taxa, -flw_cover, -flw_number)

spAll <- rbind(sp2018,sp2019,sp2020,sp2021)%>%
  group_by(year, watershed, block, plot, sppnum)%>%
  summarise(max_cover=max(cover))%>%
  ungroup()%>%
  left_join(read.csv('PPS011_new KNZ spp list.csv'))%>%
  filter(gen %notin% c('litter', 'rock', 'dung', 'bare_ground', 'bison_trail'))


##### relative cover #####
totCover <- spAll%>%
  group_by(year, watershed, block, plot)%>%
  summarise(total_cover=sum(max_cover))%>% #calculate total cover
  ungroup()

relCover <- spAll%>%
  left_join(totCover)%>%
  mutate(rel_cover=100*(max_cover/total_cover))%>% #calculate relative cover
  select(-total_cover)%>%
  mutate(replicate=paste(watershed, block, plot, sep='::'))


##### community metrics #####
commMetrics <- community_structure(relCover, time.var='year', abundance.var='rel_cover', replicate.var='replicate')%>%
  separate(replicate, into=c('watershed', 'block', 'plot'), sep='::')%>%
  mutate(plot=as.integer(plot))%>%
  left_join(trt)

#######START HERE
summary(richModel <- lme(richness~watershed*trt*year,
                               data=biomassLater,
                               random=~1|block/trt,
                               correlation=corCompSymm(form=~year|block/trt), 
                               control=lmeControl(returnObject=T)))
anova.lme(totBiomassModel, type='sequential') 
emmeans(totBiomassModel, pairwise~year*trt, adjust="tukey")




##### community difference #####





#to do: richness, evenness, comp diff, spp turnover?, trends for dominant spp, simper for which forbs come to dominate (same in bison vs invert areas?)

























# #finding the most frequent species across all plots and years for updating the datasheets
# freq <- spAll%>%
#   mutate(genus_species=paste(genus,species, sep='_'))%>%
#   filter(max_cover>0, genus_species!='NA_NA')%>%
#   group_by(watershed, block, plot, sppnum, genus_species, lifeform)%>%
#   summarize(cover=mean(max_cover))%>%
#   ungroup()%>%
#   group_by(sppnum, genus_species, lifeform, watershed)%>%
#   summarise(freq=length(genus_species))%>%
#   ungroup()%>%
#   spread(key=watershed, value=freq)
# 
# # write.csv(freq, 'species_frequency_2018-2019.csv', row.names=F)
# 
# #finding the most frequent species across all plots and years for updating the datasheets
# abund <- spAll%>%
#   mutate(genus_species=paste(genus,species, sep='_'))%>%
#   group_by(watershed, block, plot, sppnum, genus_species, lifeform)%>%
#   summarize(cover=mean(max_cover))%>%
#   ungroup()%>%
#   filter(cover>0, genus_species!='NA_NA')%>%
#   group_by(sppnum, genus_species, lifeform, watershed)%>%
#   summarise(freq=length(genus_species), avg_cover=mean(cover))%>%
#   ungroup()
# 
# # write.csv(abund, 'conSME_species_dominance_2018-2019.csv', row.names=F)