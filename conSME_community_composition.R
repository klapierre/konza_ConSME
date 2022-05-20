################################################################################
##  conSME_community_composition.R: Analysis of plant community composition responses in conSME experiment.
##
##  Author: Kimberly Komatsu
##  Date created: December 8, 2021
################################################################################

library(codyn)
library(nlme)
library(emmeans)
library(vegan)
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
  filter(gen %notin% c('litter', 'rock', 'dung', 'bare_ground', 'bison_trail'))%>%
  mutate(genus_species=paste(genus, species, sep='_'))


##### relative cover #####
totCover <- spAll%>%
  group_by(year, watershed, block, plot)%>%
  summarise(total_cover=sum(max_cover))%>% #calculate total cover
  ungroup()

relCover <- spAll%>%
  left_join(totCover)%>%
  mutate(rel_cover=100*(max_cover/total_cover))%>% #calculate relative cover
  select(-total_cover)%>%
  mutate(replicate=paste(watershed, block, plot, sep='::'))%>%
  left_join(trt)%>%
  filter(!is.na(sppnum)) #remove 5 entries that were unknowns


##### community metrics #####
commMetrics <- community_structure(relCover, time.var='year', abundance.var='rel_cover', replicate.var='replicate')%>%
  separate(replicate, into=c('watershed', 'block', 'plot'), sep='::')%>%
  mutate(plot=as.integer(plot))%>%
  left_join(trt)

##### richness response #####
summary(richModel <- lme(richness~watershed*trt*year,
                               data=commMetrics,
                               random=~1|block/trt,
                               correlation=corCompSymm(form=~year|block/trt), 
                               control=lmeControl(returnObject=T)))
anova.lme(richModel, type='sequential') 
emmeans(richModel, pairwise~year*trt, adjust="tukey")

#figure - richness by year
ggplot(data=barGraphStats(data=commMetrics, variable="richness", byFactorNames=c("trt", "year")), aes(x=trt, y=mean)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black', fill='white') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Total Biomass (g m'^'-2',')'))) +
  scale_x_discrete(limits=c('BSI', 'BSX', 'XSI', 'XXI', 'XSX', 'XXX')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  # coord_cartesian(ylim=c(0,750)) +
  facet_grid(cols=vars(year))
#export at 1400x600



##### evenness response #####
summary(evarModel <- lme(Evar~watershed*trt*year,
                         data=subset(commMetrics, year>2018),
                         random=~1|block/trt,
                         correlation=corCompSymm(form=~year|block/trt), 
                         control=lmeControl(returnObject=T)))
anova.lme(evarModel, type='sequential') 
emmeans(evarModel, pairwise~year*trt*watershed, adjust="tukey")

anova.lme(evarModel <- lme(Evar~trt,
                           data=subset(commMetrics, year==2021 & watershed=='N4B'),
                           random=~1|block/trt))
emmeans(evarModel, pairwise~trt, adjust="tukey")

#figure - evenness by year
ggplot(data=barGraphStats(data=subset(commMetrics, year>2018), variable="Evar", byFactorNames=c("trt", "year", "watershed")), aes(x=trt, y=mean)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black', fill='white') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Evenness'))) +
  scale_x_discrete(limits=c('BSI', 'BSX', 'XSI', 'XXI', 'XSX', 'XXX')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  coord_cartesian(ylim=c(0,0.55)) +
  facet_grid(cols=vars(year), rows=vars(watershed))
#export at 1400x1200



##### community difference #####
commDiff <- multivariate_difference(df=relCover, time.var='year', species.var='genus_species', abundance.var='rel_cover', replicate.var='replicate', treatment.var='trt', reference.treatment='BSI')

ggplot(data=subset(commDiff, year>2018), aes(x=year, y=composition_diff, color=trt2)) +
  geom_point() +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F)



##### PERMANOVA #####
relCover2021 <- relCover%>%
  select(year, replicate, trt, genus_species, rel_cover)%>%
  pivot_wider(names_from='genus_species', values_from='rel_cover', values_fill=list(rel_cover=0))%>%
  filter(year==2021)

print(permanova <- adonis2(formula = relCover2021[,4:183]~trt, data=relCover2021, permutations=999, method="bray"))
#F=3.487, df=5,104, p=0.001

#betadisper
veg <- vegdist(relCover2021[,4:183], method = "bray")
dispersion <- betadisper(veg, relCover2021$trt)
permutest(dispersion, pairwise=TRUE, permutations=999) 
#F=3.982, df=5,99, p=0.005

sppBC <- metaMDS(relCover2021[,4:183])

plotData <- relCover2021[,2:3]

#Use the vegan ellipse function to make ellipses
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100)
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

BC_NMDS = data.frame(MDS1 = sppBC$points[,1], MDS2 = sppBC$points[,2],group=relCover2021$trt)
BC_NMDS_Graph <- cbind(plotData,BC_NMDS)
BC_Ord_Ellipses<-ordiellipse(sppBC, plotData$trt, display = "sites",
                             kind = "se", conf = 0.95, label = T)               

ord3 <- data.frame(plotData,scores(sppBC,display="sites"))%>%
  group_by(trt)

BC_Ord_Ellipses<-ordiellipse(sppBC, plotData$trt, display = "sites",
                             kind = "se", conf = 0.95, label = T)
BC_Ellipses <- data.frame() #Make a new empty data frame called BC_Ellipses  
for(g in unique(BC_NMDS$group)){
  BC_Ellipses <- rbind(BC_Ellipses, cbind(as.data.frame(with(BC_NMDS[BC_NMDS$group==g,], 
                                                             veganCovEllipse(BC_Ord_Ellipses[[g]]$cov,BC_Ord_Ellipses[[g]]$center,BC_Ord_Ellipses[[g]]$scale)))
                                          ,group=g))
} #Generate ellipses points

ggplot(BC_NMDS_Graph, aes(x=MDS1, y=MDS2, color=group,linetype = group, shape = group)) +
  geom_point(size=6)+ 
  geom_path(data = BC_Ellipses, aes(x = NMDS1, y = NMDS2), size = 3) +
  labs(color="", linetype = "", shape = "") +
  scale_colour_manual(values=c("brown", "brown", "dark green", "dark green", "dark green", "dark green"), name = "") +
  scale_linetype_manual(values = c("twodash", "solid", "twodash", "solid", "twodash", "solid"), name = "") +
  xlab("NMDS1")+ 
  ylab("NMDS2")+ 
  theme(axis.text.x=element_text(size=24, color = "black"), axis.text.y = element_text(size = 24, color = "black"), legend.text = element_text(size = 24))



##### simper #####
print(sim <- with(relCover2021, simper(relCover2021[,4:183], trt)))


##### trends for dominant species #####
ggplot(barGraphStats(data=subset(relCover, genus_species=='andropogon_gerardii' & year>2018), variable="max_cover", byFactorNames=c("year", "trt", "watershed")), aes(x=year, y=mean, color=trt)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  ylab('Andropogon gerardii cover') +
  facet_wrap(~watershed)

ggplot(barGraphStats(data=subset(relCover, genus_species=='bouteloua_curtipendula' & year>2018), variable="max_cover", byFactorNames=c("year", "trt", "watershed")), aes(x=year, y=mean, color=trt)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  ylab('Bouteloua curtipendula cover') +
  facet_wrap(~watershed)

ggplot(barGraphStats(data=subset(relCover, genus_species=='lespedeza_violacea' & year>2018), variable="max_cover", byFactorNames=c("year", "trt", "watershed")), aes(x=year, y=mean, color=trt)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  ylab('Lespedeza violacea cover') +
  facet_wrap(~watershed)

ggplot(barGraphStats(data=subset(relCover, genus_species=='ambrosia_psilostachya' & year>2018), variable="max_cover", byFactorNames=c("year", "trt", "watershed")), aes(x=year, y=mean, color=trt)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  ylab('Ambrosia psilostachya cover') +
  facet_wrap(~watershed)

ggplot(barGraphStats(data=subset(relCover, genus_species=='oxalis_violacea' & year>2018), variable="max_cover", byFactorNames=c("year", "trt", "watershed")), aes(x=year, y=mean, color=trt)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  ylab('Oxalis violacea cover') +
  facet_wrap(~watershed)



#spp turnover?

























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