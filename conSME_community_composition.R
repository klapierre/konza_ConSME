################################################################################
##  conSME_community_composition.R: Analysis of plant community composition responses in conSME experiment.
##
##  Author: Kimberly Komatsu
##  Date created: December 8, 2021
##  Modified by: Allison Louthan Feb 26, 2025
################################################################################

library(codyn)
library(nlme)
library(emmeans)
library(vegan)
library(tidyverse)

user <- "AL" # change based on your initials to deal with directory issues. options= "KK", "AL


if (user== "KK"){setwd('C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\conSME\\data')}


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


if (user== "KK") {
  trt <- read.csv('conSME_treatments.csv')
  sp2018 <- read.csv('species composition\\ConSME_species composition_2018.csv')
  sp2019 <- read.csv('species composition\\ConSME_species composition_2019.csv')
  sp2020 <- read.csv('species composition\\ConSME_species composition_2020.csv') %>%
            select(-X, -taxa)
  sp2021 <- read.csv('species composition\\ConSME_species composition_2021.csv') %>%
            select(-taxa, -flw_cover, -flw_number)
  sp2022 <- read.csv('species composition\\ConSME_species composition_2022.csv') %>% 
            select(-taxa) %>% 
            filter(!(is.na(cover)))
  sp2023 <- read.csv('species composition\\ConSME_species composition_2023.csv') %>% 
            select(-taxa, -flowernum) %>% 
            filter(!(is.na(cover)), cover>0)

spAll <- rbind(sp2018,sp2019,sp2020,sp2021,sp2022,sp2023) %>%
  group_by(year, watershed, block, plot, sppnum) %>%
  summarise(max_cover=max(cover)) %>%
  ungroup() %>%
  left_join(read.csv('species composition\\PPS011_new KNZ spp list.csv')) %>%
  filter(gen %notin% c('litter', 'rock', 'dung', 'bare_ground', 'bison_trail')) %>%
  mutate(genus_species=paste(genus, species, sep='_'))
}
if (user == "AL"){
  trt <- read.csv('data/conSME_treatments.csv')
  sp2018 <- read.csv('data/species composition/ConSME_species composition_2018.csv')
  sp2019 <- read.csv('data/species composition/ConSME_species composition_2019.csv')
  sp2020 <- read.csv('data/species composition/ConSME_species composition_2020.csv') %>%
    select(-X, -taxa)
  sp2021 <- read.csv('data/species composition/ConSME_species composition_2021.csv') %>%
    select(-taxa, -flw_cover, -flw_number)
  sp2022 <- read.csv('data/species composition/ConSME_species composition_2022.csv') %>% 
    select(-taxa) %>% 
    filter(!(is.na(cover)))
  sp2023 <- read.csv('data/species composition/ConSME_species composition_2023.csv') %>% 
    select(-taxa, -flowernum) %>% 
    filter(!(is.na(cover)), cover>0)
  
  spAll <- rbind(sp2018,sp2019,sp2020,sp2021,sp2022,sp2023) %>%
    group_by(year, watershed, block, plot, sppnum) %>%
    summarise(max_cover=max(cover)) %>%
    ungroup() %>%
    left_join(read.csv('data/species composition/PPS011_new KNZ spp list.csv')) %>%
    filter(gen %notin% c('litter', 'rock', 'dung', 'bare_ground', 'bison_trail')) %>%
    mutate(genus_species=paste(genus, species, sep='_'))
}

##### relative cover #####
totCover <- spAll %>%
  group_by(year, watershed, block, plot) %>%
  summarise(total_cover=sum(max_cover)) %>% #calculate total cover
  ungroup()

relCover <- spAll %>%
  left_join(totCover) %>%
  mutate(rel_cover=100*(max_cover/total_cover)) %>% #calculate relative cover
  select(-total_cover) %>%
  mutate(replicate=paste(watershed, block, plot, sep='::')) %>%
  left_join(trt)


##### community metrics #####
commMetrics <- community_structure(relCover, time.var='year', abundance.var='rel_cover', replicate.var='replicate') %>%
  separate(replicate, into=c('watershed', 'block', 'plot'), sep='::') %>%
  mutate(plot=as.integer(plot)) %>%
  left_join(trt) %>% 
  mutate(ws_label=ifelse(watershed=='N1A', 'Annual', '4 Year'),
         experiment_year=year-2018)

hist(log(commMetrics$richness)) #looks as good as we can get!
shapiro.test(log(commMetrics$richness))
# W = 0.98961, p-value = 0.0001609

hist(log(commMetrics$Evar))
shapiro.test(log(commMetrics$Evar))
# W = 0.99854, p-value = 0.881


##### richness response #####
summary(richModel <- lme(log(richness)~watershed*as.factor(year)*invertebrates*bison + watershed*as.factor(year)*invertebrates*small_mammal,
                               data=subset(commMetrics, year>2018),
                               random=~1|block/trt,
                               correlation=corAR1(form=~year|block/trt), 
                               control=lmeControl(returnObject=T)))
anova.lme(richModel, type='sequential') 
emmeans(richModel, pairwise~year*watershed*bison, adjust="tukey")
emmeans(richModel, pairwise~invertebrates, adjust="tukey")

#figure - richness bison by year*ws, either put in supplement or don't include at all
richnessBisonSupplementalFig <- ggplot(data=barGraphStats(data=subset(commMetrics, year>2018), variable="richness", byFactorNames=c("bison", "year", "ws_label")), aes(x=bison, y=mean, fill=bison)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  # coord_cartesian(ylim=c(0,23)) +
  ylab(expression(paste('Plant Species\nRichness'))) +
  scale_fill_manual(values=c('lightgrey', 'limegreen')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=35), axis.title.y=element_text(size=35, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=35), legend.position='none', legend.justification=c(0,1), strip.text=element_text(size=35)) +
  facet_grid(cols=vars(year), rows=vars(ws_label))

# ggsave('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\konza projects\\conSME\\figures\\2023\\invert_richness.png', width=4, height=7, units='in', dpi=600, bg='white')

#figure - richness by bison and ws (include in main text and show year interaction in a supplement?)
richnessBisonFig <- ggplot(data=barGraphStats(data=subset(commMetrics, year>2018), variable="richness", byFactorNames=c("bison", 'watershed')), aes(x=watershed, y=mean, fill=bison)) +
  geom_bar(position=position_dodge(0.9), size=2, stat="identity", color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Plant Species\nRichness'))) +
  # ylab('') +
  coord_cartesian(ylim=c(0,30)) +
  scale_fill_manual(values=c('lightgrey', 'limegreen')) +
  scale_x_discrete(labels=c('Annual', '4 Year')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=35), axis.title.y=element_text(size=35, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=35), legend.position=c(0.98, 0.98), legend.justification=c(1,1), strip.text=element_text(size=35), legend.text=element_text(size=35))
# ggsave('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\konza projects\\conSME\\figures\\2023\\bison_richness.png', width=7, height=7, units='in', dpi=600, bg='white')

# #figure - richness by year, bison - don't include figure (it's obv from stats)
# ggplot(data=barGraphStats(data=subset(commMetrics, year>2018), variable="richness", byFactorNames=c("bison", 'year')), aes(x=as.factor(year), y=mean, fill=bison)) +
#   geom_bar(position=position_dodge(0.9), size=2, stat="identity", color='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
#   ylab(expression(paste('Plant Species Richness'))) +
#   coord_cartesian(ylim=c(0,30)) +
#   scale_fill_manual(values=c('lightgrey', 'limegreen')) +
#   # scale_x_discrete(labels=c('Annual', '4 Year')) +
#   theme(axis.title.x=element_blank(), axis.text.x=element_text(size=35), axis.title.y=element_text(size=35, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=35), legend.position=c(0.98, 0.98), legend.justification=c(1,1), strip.text=element_text(size=35), legend.text=element_text(size=35)) 
# # ggsave('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\konza projects\\conSME\\figures\\2023\\bison_year_richness.png', width=7, height=7, units='in', dpi=600, bg='white')

#figure - richness by invertebrates
richnessInvertFig <- ggplot(data=barGraphStats(data=subset(commMetrics, year>2018), variable="richness", byFactorNames=c("invertebrates")), aes(x=invertebrates, y=mean, fill=invertebrates)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  coord_cartesian(ylim=c(0,23)) +
  # ylab(expression(paste('Plant Species\nRichness'))) +
  ylab('') +
  scale_fill_manual(values=c('lightgrey', 'skyblue')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=35), axis.title.y=element_text(size=35, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=35), legend.position='none', legend.justification=c(0,1), strip.text=element_text(size=35))
# ggsave('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\konza projects\\conSME\\figures\\2023\\invert_richness.png', width=4, height=7, units='in', dpi=600, bg='white')

#figure - richness by small mammals
richnessSmallMammalFig <- ggplot(data=barGraphStats(data=subset(commMetrics, year>2018), variable="richness", byFactorNames=c("small_mammal", "year", "ws_label")), aes(x=small_mammal, y=mean, fill=small_mammal)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  # coord_cartesian(ylim=c(0,23)) +
  ylab(expression(paste('Plant Species\nRichness'))) +
  scale_fill_manual(values=c('lightgrey', 'orange')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=35), axis.title.y=element_text(size=35, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=35), legend.position='none', legend.justification=c(0,1), strip.text=element_text(size=35)) +
  facet_grid(cols=vars(year), rows=vars(ws_label))
# ggsave('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\konza projects\\conSME\\figures\\2023\\invert_richness.png', width=4, height=7, units='in', dpi=600, bg='white')

plot_grid(arrangeGrob(richnessBisonFig, richnessInvertFig, ncol=2), 
          richnessSmallMammalFig,
          nrow=2, ncol=1, rel_heights = c(0.5, 1))
# export to file Fig 4_richness_interactiveEffects.png at 1600x1600

##### response ratios #####
###bison
commMetricsMeans <- commMetrics %>%
  group_by(watershed, bison) %>%
  summarise(richness_mean=mean(richness), sd=sd(richness), N=length(richness)) %>%
  ungroup() %>%
  mutate(se=sd/sqrt(N))

###small mammal
commMetricsMeans <- commMetrics %>%
  filter(year %in% c(2020, 2023)) %>% 
  group_by(watershed, small_mammal) %>%
  summarise(richness_mean=mean(richness), sd=sd(richness), N=length(richness)) %>%
  ungroup() %>%
  mutate(se=sd/sqrt(N))

###invert
commMetricsMeans <- commMetrics %>%
  group_by(invertebrates) %>%
  summarise(richness_mean=mean(richness), sd=sd(richness), N=length(richness)) %>%
  ungroup() %>%
  mutate(se=sd/sqrt(N))

  
##### evenness response #####
summary(evarModel <- lme(log(Evar)~watershed*as.factor(year)*invertebrates*bison + watershed*as.factor(year)*invertebrates*small_mammal,
                         data=subset(commMetrics, year>2018),
                         random=~1|block/trt,
                         correlation=corCompSymm(form=~year|block/trt), 
                         control=lmeControl(returnObject=T)))
anova.lme(evarModel, type='sequential') 
emmeans(evarModel, pairwise~year*trt*watershed, adjust="tukey")
 if (user== "AL"){save(richModel,evarModel, file='derived_data/02comp_models.RData')}
#figure - evenness by watershed, year, bison
evennessBisonFig <- ggplot(data=barGraphStats(data=subset(commMetrics, year>2018), variable="Evar", byFactorNames=c("bison", 'ws_label', "year")), aes(x=bison, y=mean, fill=bison)) +
  geom_bar(position=position_dodge(0.9), size=2, stat="identity", color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Plant Species Evenness'))) +
  # ylab('') +
  coord_cartesian(ylim=c(0,0.6)) +
  scale_fill_manual(values=c('lightgrey', 'limegreen')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=35), 
        axis.title.y=element_text(size=35, angle=90, vjust=1, margin=margin(r=15)), 
        axis.text.y=element_text(size=35), 
        legend.position='none', strip.text=element_text(size=35), 
        legend.text=element_text(size=35)) +
  facet_grid(cols=vars(year), rows=vars(ws_label))
# ggsave('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\konza projects\\conSME\\figures\\2023\\bison_evenness.png', width=7, height=7, units='in', dpi=600, bg='white')

#figure - evenness by watershed, year, small mammals
evennessSmallMammalFig <- ggplot(data=barGraphStats(data=subset(commMetrics, year>2018), variable="Evar", byFactorNames=c("small_mammal", "year", "ws_label")), aes(x=small_mammal, y=mean, fill=small_mammal)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  coord_cartesian(ylim=c(0,0.6)) +
  ylab(expression(paste('Plant Species Evenness'))) +
  scale_fill_manual(values=c('lightgrey', 'orange')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=35), axis.title.y=element_text(size=35, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=35), legend.position='none', legend.justification=c(0,1), strip.text=element_text(size=35)) +
  facet_grid(cols=vars(year), rows=vars(ws_label))
# ggsave('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\konza projects\\conSME\\figures\\2023\\invert_richness.png', width=4, height=7, units='in', dpi=600, bg='white')

plot_grid(evennessBisonFig, evennessSmallMammalFig, nrow=2, ncol=1)
# export to file Fig 5_evenness_interactiveEffects.png at 1600x1400


##### community difference #####
commDiff <- multivariate_difference(df=relCover, time.var='year', species.var='genus_species', abundance.var='rel_cover', replicate.var='replicate', treatment.var='trt', reference.treatment='BSI')

ggplot(data=subset(commDiff, year>2018), aes(x=year, y=composition_diff, color=trt2)) +
  geom_point() +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F)
# bison removal overwhelms all community change


#looking at just insect effects with bison
commDiffBison <- multivariate_difference(df=subset(relCover, trt %in% c('BSX', 'BSI')), time.var='year', species.var='genus_species', abundance.var='rel_cover', replicate.var='replicate', treatment.var='trt', reference.treatment='BSI')

ggplot(data=subset(commDiffBison, year>2018), aes(x=year, y=composition_diff, color=trt2)) +
  geom_point() +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F)
# increasing difference through time (except 2022, where it goes way down)


#looking at effects without bison
commDiffNoBison <- multivariate_difference(df=subset(relCover, !(trt %in% c('BSX', 'BSI'))), time.var='year', species.var='genus_species', abundance.var='rel_cover', replicate.var='replicate', treatment.var='trt', reference.treatment='XSI')

ggplot(data=subset(commDiffNoBison, year>2018), aes(x=year, y=composition_diff, color=trt2)) +
  geom_point() +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F)
# insects create community change


##### PERMANOVA #####
relCover2023 <- relCover %>%
  # filter(bison=='X') %>% 
  mutate(bison_ws=paste(bison, watershed, sep='::')) %>% 
  select(year, watershed, block, replicate, trt, bison_ws, bison, small_mammal, invertebrates, genus_species, rel_cover) %>%
  pivot_wider(names_from='genus_species', values_from='rel_cover', values_fill=list(rel_cover=0)) %>%
  filter(year==2023)

print(permanova <- adonis(formula = relCover2023[,10:204]~
                            year*watershed*bison*invertebrates + 
                            year*watershed*small_mammal*invertebrates, 
                          data=relCover2023, permutations=999, method="bray"))
result_table <- as.data.frame(permanova$aov.tab)
#watershed*bison F=1.88, df=1,104, p=0.067; bison F=14.47, df=1,104, p=0.001


#betadisper
veg <- vegdist(relCover2023[,10:204], method = "bray")
dispersion <- betadisper(veg, relCover2023$bison)
permutest(dispersion, pairwise=TRUE, permutations=999) 
#F=1.4561, p=0.226

sppBC <- metaMDS(relCover2023[,10:204])

plotData <- relCover2023[,1:9]

#bison and watershed NMDS
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100)
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

BC_NMDS = data.frame(MDS1 = sppBC$points[,1], MDS2 = sppBC$points[,2],group=relCover2023$bison_ws)
BC_NMDS_Graph <- cbind(plotData,BC_NMDS)
BC_Ord_Ellipses<-ordiellipse(sppBC, plotData$bison_ws, display = "sites",
                             kind = "se", conf = 0.95, label = T)               

ord3 <- data.frame(plotData,scores(sppBC,display="sites")) %>%
  group_by(bison_ws)

BC_Ord_Ellipses<-ordiellipse(sppBC, plotData$bison_ws, display = "sites",
                             kind = "se", conf = 0.95, label = T)
BC_Ellipses <- data.frame() #Make a new empty data frame called BC_Ellipses  
for(g in unique(BC_NMDS$group)){
  BC_Ellipses <- rbind(BC_Ellipses, cbind(as.data.frame(with(BC_NMDS[BC_NMDS$group==g,], 
                                                             veganCovEllipse(BC_Ord_Ellipses[[g]]$cov,BC_Ord_Ellipses[[g]]$center,BC_Ord_Ellipses[[g]]$scale)))
                                          ,group=g))
} #Generate ellipses points

ggplot(BC_NMDS_Graph, aes(x=MDS1, y=MDS2, color=group,linetype = group, shape = group)) +
  geom_point(size=3)+ 
  geom_path(data = BC_Ellipses, aes(x = NMDS1, y = NMDS2), size = 1) +
  labs(color="", linetype = "", shape = "") +
  scale_colour_manual(values=c("brown", "brown", "dark green", "dark green", "dark green", "dark green"), name = "") +
  scale_linetype_manual(values = c("twodash", "solid", "twodash", "solid", "twodash", "solid"), name = "") +
  xlab("NMDS1")+ 
  ylab("NMDS2")+ 
  theme(axis.text.x=element_text(size=24, color = "black"), axis.text.y = element_text(size = 24, color = "black"), legend.text = element_text(size = 24))
# ggsave(file='C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\conSME\\figures\\2025\\Fig S2_NMDSbison2023.png', width=10, height=8, units='in', dpi=300, bg='white')


#invertebrate NMDS
BC_NMDS = data.frame(MDS1 = sppBC$points[,1], MDS2 = sppBC$points[,2],group=relCover2023$invertebrates)
BC_NMDS_Graph <- cbind(plotData,BC_NMDS)
BC_Ord_Ellipses<-ordiellipse(sppBC, plotData$invertebrates, display = "sites",
                             kind = "se", conf = 0.95, label = T)               

ord3 <- data.frame(plotData,scores(sppBC,display="sites")) %>%
  group_by(invertebrates)

BC_Ord_Ellipses<-ordiellipse(sppBC, plotData$invertebrates, display = "sites",
                             kind = "se", conf = 0.95, label = T)
BC_Ellipses <- data.frame() #Make a new empty data frame called BC_Ellipses  
for(g in unique(BC_NMDS$group)){
  BC_Ellipses <- rbind(BC_Ellipses, cbind(as.data.frame(with(BC_NMDS[BC_NMDS$group==g,], 
                                                             veganCovEllipse(BC_Ord_Ellipses[[g]]$cov,BC_Ord_Ellipses[[g]]$center,BC_Ord_Ellipses[[g]]$scale)))
                                          ,group=g))
} #Generate ellipses points

ggplot(BC_NMDS_Graph, aes(x=MDS1, y=MDS2, color=group,linetype = group, shape = group)) +
  geom_point(size=3)+ 
  geom_path(data = BC_Ellipses, aes(x = NMDS1, y = NMDS2), size = 1) +
  labs(color="", linetype = "", shape = "") +
  scale_colour_manual(values=c("lightgrey", "skyblue"), name = "") +
  scale_linetype_manual(values = c("twodash", "solid", "twodash", "solid", "twodash", "solid"), name = "") +
  xlab("NMDS1")+ 
  ylab("NMDS2")+ 
  theme(axis.text.x=element_text(size=24, color = "black"), axis.text.y = element_text(size = 24, color = "black"), legend.text = element_text(size = 24))
# if (user= "KK"){ggsave(file='C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\conSME\\figures\\2025\\Fig S2_NMDSbison2023.png', width=10, height=8, units='in', dpi=300, bg='white')}


##### simper #####
#bison and watershed
summary(sim <- with(relCover2023, simper(relCover2023[,10:204], bison_ws)))

#invertebrates
summary(sim <- with(relCover2023, simper(relCover2023[,10:204], invertebrates)))


##### trends for dominant species #####

###bison effect
#andropogon gerardii
ggplot(barGraphStats(data=subset(relCover, genus_species=='andropogon_gerardii' & year>2018), variable="max_cover", byFactorNames=c("year", "bison", "watershed")), aes(x=year, y=mean, color=bison)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  facet_wrap(~watershed)

#lespedeza violacea
ggplot(barGraphStats(data=subset(relCover, genus_species=='lespedeza_violacea' & year>2018), variable="max_cover", byFactorNames=c("year", "bison", "watershed")), aes(x=year, y=mean, color=bison)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  facet_wrap(~watershed)

#bouteloua dactyloides
ggplot(barGraphStats(data=subset(relCover, genus_species=='bouteloua_dactyloides' & year>2018), variable="max_cover", byFactorNames=c("year", "bison", "watershed")), aes(x=year, y=mean, color=bison)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  facet_wrap(~watershed)

#vernonia baldwinii
ggplot(barGraphStats(data=subset(relCover, genus_species=='vernonia_baldwinii' & year>2018), variable="max_cover", byFactorNames=c("year", "bison", "watershed")), aes(x=year, y=mean, color=bison)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  facet_wrap(~watershed)

#dalea multiflora
ggplot(barGraphStats(data=subset(relCover, genus_species=='dalea_multiflora' & year>2018), variable="max_cover", byFactorNames=c("year", "bison", "watershed")), aes(x=year, y=mean, color=bison)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  facet_wrap(~watershed)

#schizachyrium scoparium
ggplot(barGraphStats(data=subset(relCover, genus_species=='schizachyrium_scoparium' & year>2018), variable="max_cover", byFactorNames=c("year", "bison", "watershed")), aes(x=year, y=mean, color=bison)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  facet_wrap(~watershed)

#bouteloua curtipendula
ggplot(barGraphStats(data=subset(relCover, genus_species=='bouteloua_curtipendula' & year>2018), variable="max_cover", byFactorNames=c("year", "bison", "watershed")), aes(x=year, y=mean, color=bison)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  facet_wrap(~watershed)

#sporobolus compositus
ggplot(barGraphStats(data=subset(relCover, genus_species=='sporobolus_compositus' & year>2018), variable="max_cover", byFactorNames=c("year", "bison", "watershed")), aes(x=year, y=mean, color=bison)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  facet_wrap(~watershed)

#bouteloua hirsuta
ggplot(barGraphStats(data=subset(relCover, genus_species=='bouteloua_hirsuta' & year>2018), variable="max_cover", byFactorNames=c("year", "bison", "watershed")), aes(x=year, y=mean, color=bison)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  facet_wrap(~watershed)

###invert effect
#symphotrichum ericoides
ggplot(barGraphStats(data=subset(relCover, genus_species=='symphyotrichum_ericoides' & year>2018), variable="max_cover", byFactorNames=c("year", "invertebrates", "watershed")), aes(x=year, y=mean, color=invertebrates)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  facet_wrap(~watershed)

#dichanthelium oligosanthes
ggplot(barGraphStats(data=subset(relCover, genus_species=='dichanthelium_oligosanthes' & year>2018), variable="max_cover", byFactorNames=c("year", "invertebrates", "watershed")), aes(x=year, y=mean, color=invertebrates)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  facet_wrap(~watershed)

#echinacia angustifolia
ggplot(barGraphStats(data=subset(relCover, genus_species=='echinacea_angustifolia' & year>2018), variable="max_cover", byFactorNames=c("year", "invertebrates", "watershed")), aes(x=year, y=mean, color=invertebrates)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  facet_wrap(~watershed)

#liatris punctata
ggplot(barGraphStats(data=subset(relCover, genus_species=='liatris_punctata' & year>2018), variable="max_cover", byFactorNames=c("year", "invertebrates", "watershed")), aes(x=year, y=mean, color=invertebrates)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  facet_wrap(~watershed)


# ##### RACs #####
# 
# #bison by watershed
# rankAbundance <- relCover %>%
#   filter(year==2023) %>%
#   mutate(spp_name=str_to_sentence(paste(genus, species, sep=' '))) %>%
#   group_by(watershed, bison, spp_name, growthform, lifeform) %>%
#   summarize(avg_cover=mean(rel_cover)) %>%
#   ungroup() %>%
#   arrange(bison, watershed, -avg_cover) %>%
#   mutate(bison_ws=paste(bison, watershed, sep='::')) %>%
#   group_by(bison_ws) %>%
#   mutate(rank=seq_along(bison_ws)) %>%
#   ungroup() %>%
#   mutate(lifeform2=ifelse(spp_name=='Sisyrinchium campestre', 'f', ifelse(lifeform=='o', 'f', ifelse(lifeform=='s', 'g', as.character(lifeform))))) %>%
#   filter(spp_name!='NA NA')
# 
# ggplot(data=subset(rankAbundance, bison_ws=='B::N1A', avg_cover>0), aes(x=rank, y=avg_cover)) +
#   geom_line() +
#   geom_point(aes(colour=lifeform2, shape=growthform), size=3) +
#   scale_color_manual(labels=c("Forb", "Graminoid", "Woody"), values=c('#DE8C00', '#009CC0', '#83431E')) +
#   scale_shape_discrete(labels=c("Annual", "Biennial", "Perennial"))+
#   xlab('') +
#   ylab('N1A Bison\nRelative Percent Cover\n') +
#   # scale_x_continuous(expand=c(0,0), limits=c(0.5,17), breaks=seq(0,17,5)) +
#   # scale_y_continuous(expand=c(0,0), limits=c(0,60), breaks=seq(0,60,10)) +
#   # geom_text(aes(y=avg_cover+1.2, x=rank+0.1, label=spp_name), hjust='left', vjust='center', angle=90, size=4) +
#   expand_limits(y=30, x=90)
# #export at 1400x400
# 
# ggplot(data=subset(rankAbundance, bison_ws=='X::N1A', avg_cover>0), aes(x=rank, y=avg_cover)) +
#   geom_line() +
#   geom_point(aes(colour=lifeform2, shape=growthform), size=3) +
#   scale_color_manual(labels=c("Forb", "Graminoid", "Woody"), values=c('#DE8C00', '#009CC0', '#83431E')) +
#   scale_shape_discrete(labels=c("Annual", "Biennial", "Perennial"))+
#   xlab('') +
#   ylab('N1A Bison Removed\nRelative Percent Cover\n') +
#   # scale_x_continuous(expand=c(0,0), limits=c(0.5,17), breaks=seq(0,17,5)) +
#   # scale_y_continuous(expand=c(0,0), limits=c(0,60), breaks=seq(0,60,10)) +
#   # geom_text(aes(y=avg_cover+1.2, x=rank+0.1, label=spp_name), hjust='left', vjust='center', angle=90, size=4) +
#   expand_limits(y=30, x=90)
# #export at 1400x400
# 
# ggplot(data=subset(rankAbundance, bison_ws=='B::N4B', avg_cover>0), aes(x=rank, y=avg_cover)) +
#   geom_line() +
#   geom_point(aes(colour=lifeform2, shape=growthform), size=3) +
#   scale_color_manual(labels=c("Forb", "Graminoid", "Woody"), values=c('#DE8C00', '#009CC0', '#83431E')) +
#   scale_shape_discrete(labels=c("Annual", "Biennial", "Perennial"))+
#   xlab('') +
#   ylab('N4B Bison\nRelative Percent Cover\n') +
#   # scale_x_continuous(expand=c(0,0), limits=c(0.5,17), breaks=seq(0,17,5)) +
#   # scale_y_continuous(expand=c(0,0), limits=c(0,60), breaks=seq(0,60,10)) +
#   # geom_text(aes(y=avg_cover+1.2, x=rank+0.1, label=spp_name), hjust='left', vjust='center', angle=90, size=4) +
#   expand_limits(y=40, x=90)
# #export at 1400x400
# 
# ggplot(data=subset(rankAbundance, bison_ws=='X::N4B', avg_cover>0), aes(x=rank, y=avg_cover)) +
#   geom_line() +
#   geom_point(aes(colour=lifeform2, shape=growthform), size=3) +
#   scale_color_manual(labels=c("Forb", "Graminoid", "Woody"), values=c('#DE8C00', '#009CC0', '#83431E')) +
#   scale_shape_discrete(labels=c("Annual", "Biennial", "Perennial"))+
#   xlab('') +
#   ylab('N4B Bison Removed\nRelative Percent Cover\n') +
#   # scale_x_continuous(expand=c(0,0), limits=c(0.5,17), breaks=seq(0,17,5)) +
#   # scale_y_continuous(expand=c(0,0), limits=c(0,60), breaks=seq(0,60,10)) +
#   # geom_text(aes(y=avg_cover+1.2, x=rank+0.1, label=spp_name), hjust='left', vjust='center', angle=90, size=4) +
#   expand_limits(y=40, x=90)
# #export at 1400x400
# 
# 
# #invertebrate
# rankAbundance <- relCover %>%
#   filter(year==2023) %>%
#   mutate(spp_name=str_to_sentence(paste(genus, species, sep=' '))) %>%
#   group_by(invertebrates, spp_name, growthform, lifeform) %>%
#   summarize(avg_cover=mean(rel_cover)) %>%
#   ungroup() %>%
#   arrange(invertebrates, -avg_cover) %>%
#   group_by(invertebrates) %>%
#   mutate(rank=seq_along(invertebrates)) %>%
#   ungroup() %>%
#   mutate(lifeform2=ifelse(spp_name=='Sisyrinchium campestre', 'f', ifelse(lifeform=='o', 'f', ifelse(lifeform=='s', 'g', as.character(lifeform))))) %>%
#   filter(spp_name!='NA NA')
# 
# ggplot(data=subset(rankAbundance, invertebrates=='I', avg_cover>0), aes(x=rank, y=avg_cover)) +
#   geom_line() +
#   geom_point(aes(colour=lifeform2, shape=growthform), size=3) +
#   scale_color_manual(labels=c("Forb", "Graminoid", "Woody"), values=c('#DE8C00', '#009CC0', '#83431E')) +
#   scale_shape_discrete(labels=c("Annual", "Biennial", "Perennial"))+
#   xlab('') +
#   ylab('Invertebrates\nRelative Percent Cover\n') +
#   # scale_x_continuous(expand=c(0,0), limits=c(0.5,17), breaks=seq(0,17,5)) +
#   # scale_y_continuous(expand=c(0,0), limits=c(0,60), breaks=seq(0,60,10)) +
#   geom_text(aes(y=avg_cover+1.2, x=rank+0.1, label=spp_name), hjust='left', vjust='center', angle=90, size=4) +
#   expand_limits(y=30, x=90)
# #export at 1400x400
# 
# ggplot(data=subset(rankAbundance, invertebrates=='X', avg_cover>0), aes(x=rank, y=avg_cover)) +
#   geom_line() +
#   geom_point(aes(colour=lifeform2, shape=growthform), size=3) +
#   scale_color_manual(labels=c("Forb", "Graminoid", "Woody"), values=c('#DE8C00', '#009CC0', '#83431E')) +
#   scale_shape_discrete(labels=c("Annual", "Biennial", "Perennial"))+
#   xlab('') +
#   ylab('Invertebrates Removed\nRelative Percent Cover\n') +
#   # scale_x_continuous(expand=c(0,0), limits=c(0.5,17), breaks=seq(0,17,5)) +
#   # scale_y_continuous(expand=c(0,0), limits=c(0,60), breaks=seq(0,60,10)) +
#   geom_text(aes(y=avg_cover+1.2, x=rank+0.1, label=spp_name), hjust='left', vjust='center', angle=90, size=4) +
#   expand_limits(y=30, x=90)
# #export at 1400x400



# ##### finding the most frequent species across all plots and years for updating the datasheets #####
# freq <- spAll%>%
#   mutate(genus_species=paste(genus,species, sep='_')) %>%
#   filter(max_cover>0, genus_species!='NA_NA') %>%
#   group_by(watershed, block, plot, sppnum, genus_species, lifeform) %>%
#   summarize(cover=mean(max_cover)) %>%
#   ungroup() %>%
#   group_by(sppnum, genus_species, lifeform, watershed) %>%
#   summarise(freq=length(genus_species)) %>%
#   ungroup() %>%
#   spread(key=watershed, value=freq)
# 
# # write.csv(freq, 'species_frequency_2018-2019.csv', row.names=F)
# 
# #finding the most frequent species across all plots and years for updating the datasheets
# abund <- spAll%>%
#   mutate(genus_species=paste(genus,species, sep='_')) %>%
#   group_by(watershed, block, plot, sppnum, genus_species, lifeform) %>%
#   summarize(cover=mean(max_cover)) %>%
#   ungroup() %>%
#   filter(cover>0, genus_species!='NA_NA') %>%
#   group_by(sppnum, genus_species, lifeform, watershed) %>%
#   summarise(freq=length(genus_species), avg_cover=mean(cover)) %>%
#   ungroup()
# 
# # write.csv(abund, 'conSME_species_dominance_2018-2019.csv', row.names=F)