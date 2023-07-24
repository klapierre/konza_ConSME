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

setwd('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\konza projects\\conSME\\data')


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

spAll <- rbind(sp2018,sp2019,sp2020,sp2021,sp2022) %>%
  group_by(year, watershed, block, plot, sppnum) %>%
  summarise(max_cover=max(cover)) %>%
  ungroup() %>%
  left_join(read.csv('species composition\\PPS011_new KNZ spp list.csv')) %>%
  filter(gen %notin% c('litter', 'rock', 'dung', 'bare_ground', 'bison_trail')) %>%
  mutate(genus_species=paste(genus, species, sep='_'))


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
  left_join(trt) %>%
  filter(!is.na(sppnum)) #remove 5 entries that were unknowns


##### community metrics #####
commMetrics <- community_structure(relCover, time.var='year', abundance.var='rel_cover', replicate.var='replicate') %>%
  separate(replicate, into=c('watershed', 'block', 'plot'), sep='::') %>%
  mutate(plot=as.integer(plot)) %>%
  left_join(trt)

hist(log(commMetrics$richness))
shapiro.test(log(commMetrics$richness))
# W = 0.99191, p-value = 0.00514

hist(log(commMetrics$Evar))
shapiro.test(log(commMetrics$Evar))
# W = 0.99836, p-value = 0.9


##### richness response #####
summary(richModel <- lme(log(richness)~watershed*year*invertebrates*bison + watershed*year*invertebrates*small_mammal,
                               data=subset(commMetrics, year>2018),
                               random=~1|block/trt,
                               correlation=corCompSymm(form=~year|block/trt), 
                               control=lmeControl(returnObject=T)))
anova.lme(richModel, type='sequential') 
emmeans(richModel, pairwise~year*watershed*bison, adjust="tukey")
emmeans(richModel, pairwise~invertebrates, adjust="tukey")

#figure - richness by watershed, bison
ggplot(data=barGraphStats(data=subset(commMetrics, year>2018), variable="richness", byFactorNames=c("bison", 'watershed')), aes(x=watershed, y=mean, fill=bison)) +
  geom_bar(position=position_dodge(0.9), size=2, stat="identity", color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9), size=2) +
  ylab(expression(paste('Plant Species Richness'))) +
  coord_cartesian(ylim=c(0,30)) +
  scale_fill_manual(values=c('lightgrey', 'limegreen')) +
  scale_x_discrete(labels=c('Annual', '4 Year')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=35), axis.title.y=element_text(size=35, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=35), legend.position=c(0.98, 0.98), legend.justification=c(1,1), strip.text=element_text(size=35), legend.text=element_text(size=35)) 
# ggsave('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\konza projects\\conSME\\figures\\2023\\bison_richness.png', width=7, height=7, units='in', dpi=600, bg='white')


#figure - richness by invertebrates
ggplot(data=barGraphStats(data=subset(commMetrics, year>2018), variable="richness", byFactorNames=c("invertebrates")), aes(x=invertebrates, y=mean, fill=invertebrates)) +
  geom_bar(position=position_dodge(0.1), size=2, stat="identity", color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  coord_cartesian(ylim=c(0,23)) +
  ylab(expression(paste('Plant Species Richness'))) +
  scale_fill_manual(values=c('lightgrey', 'lightgoldenrod')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=35), axis.title.y=element_text(size=35, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=35), legend.position='none', legend.justification=c(0,1), strip.text=element_text(size=35))
# ggsave('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\konza projects\\conSME\\figures\\2023\\invert_richness.png', width=4, height=7, units='in', dpi=600, bg='white')


#response ratios 
commMetricsMeans <- commMetrics %>%
  group_by(watershed, year, bison) %>%
  summarise(richness_mean=mean(richness), sd=sd(richness), N=length(richness)) %>%
  ungroup() %>%
  mutate(se=sd/sqrt(N))

#figure - richness by watershed, year, bison
ggplot(data=subset(commMetricsMeans, year>2018), aes(x=as.factor(year), y=richness_mean, color=bison)) +
  geom_point(size=2) +
  geom_smooth(method='lm') +
  geom_errorbar(aes(ymin=richness_mean-se, ymax=richness_mean+se), width=.1) +
  ylab(expression(paste('Plant Species Richness'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  # coord_cartesian(ylim=c(0,750)) +
  facet_grid(rows=vars(watershed))
#export at 1400x600


commMetricsRR <- commMetrics %>%
  group_by(watershed, block, year, bison) %>%
  summarise(richness_mean=mean(richness)) %>%
  ungroup() %>%
  pivot_wider(names_from='bison', values_from='richness_mean') %>%
  mutate(richness_percent_loss=100*(B-X)/B)

temp <- commMetricsRR %>%
  group_by(watershed, year) %>%
  summarise(richness_mean_bison=mean(B), richness_mean_X=mean(X)) %>%
  ungroup()

ggplot(data=barGraphStats(data=subset(commMetricsRR, year>2018), variable="richness_percent_loss", byFactorNames=c("year", 'watershed')), aes(x=as.factor(year), y=-(mean))) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=(-mean-se), ymax=(-mean+se)), width=.1, size=2) +
  ylab('Plant Species Richness Response (%)') +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
  # coord_cartesian(ylim=c(0,50)) +
  geom_hline(yintercept=0) +
  facet_grid(rows=vars(watershed))
#export at 1400x600
  
  
  
##### evenness response #####
summary(evarModel <- lme(log(Evar)~watershed*year*invertebrates*bison + watershed*year*invertebrates*small_mammal,
                         data=subset(commMetrics, year>2018),
                         random=~1|block/trt,
                         correlation=corCompSymm(form=~year|block/trt), 
                         control=lmeControl(returnObject=T)))
anova.lme(evarModel, type='sequential') 
emmeans(evarModel, pairwise~year*trt*watershed, adjust="tukey")

#figure - evenness by small mammal and inverts
ggplot(data=barGraphStats(data=subset(commMetrics, year>2018 & trt %in% c('XSI', 'XSX', 'XXI', 'XXX')), variable="Evar", byFactorNames=c("small_mammal")), aes(x=small_mammal, y=mean)) +
  geom_point(position=position_dodge(0.1), size=5, stat="identity", color='black', fill='white') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Evenness'))) +
  # scale_x_discrete(limits=c('XSI', 'XXI', 'XSX', 'XXX')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) #+
  # coord_cartesian(ylim=c(0,0.55)) +
  # facet_grid(cols=vars(year), rows=vars(watershed))
#export at 1400x1200

# #figure - evenness by bison
# ggplot(data=barGraphStats(data=subset(commMetrics, year>2018), variable="Evar", byFactorNames=c("bison", "watershed", "year")), aes(x=year, y=mean, color=bison)) +
#   geom_point(position=position_dodge(0.1), size=5, stat="identity", color='black', fill='white') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
#   ylab(expression(paste('Evenness'))) +
#   # scale_x_discrete(limits=c('XSI', 'XXI', 'XSX', 'XXX')) +
#   theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1), strip.text=element_text(size=30)) +
# # coord_cartesian(ylim=c(0,0.55)) +
# facet_grid(rows=vars(watershed))
# #export at 1400x1200


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
# effect in 2021, but goes away in 2022


#looking at effects without bison
commDiffNoBison <- multivariate_difference(df=subset(relCover, !(trt %in% c('BSX', 'BSI'))), time.var='year', species.var='genus_species', abundance.var='rel_cover', replicate.var='replicate', treatment.var='trt', reference.treatment='XSI')

ggplot(data=subset(commDiffNoBison, year>2018), aes(x=year, y=composition_diff, color=trt2)) +
  geom_point() +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F)
# insects create community change


##### PERMANOVA #####
relCover2022 <- relCover %>%
  # filter(bison=='X') %>% 
  mutate(bison_ws=paste(bison, watershed, sep='::')) %>% 
  select(year, watershed, replicate, trt, bison_ws, bison, small_mammal, invertebrates, genus_species, rel_cover) %>%
  pivot_wider(names_from='genus_species', values_from='rel_cover', values_fill=list(rel_cover=0)) %>%
  filter(year==2022)

print(permanova <- adonis2(formula = relCover2022[,9:194]~watershed*bison*small_mammal*invertebrates, data=relCover2022, permutations=999, method="bray"))
#watershed*bison F=1.88, df=1,104, p=0.067; bison F=14.47, df=1,104, p=0.001

#betadisper
veg <- vegdist(relCover2022[,9:194], method = "bray")
dispersion <- betadisper(veg, relCover2022$invertebrates)
permutest(dispersion, pairwise=TRUE, permutations=999) 
#F=1.4561, p=0.226

sppBC <- metaMDS(relCover2022[,9:194])

plotData <- relCover2022[,1:8]

#Use the vegan ellipse function to make ellipses
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100)
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

BC_NMDS = data.frame(MDS1 = sppBC$points[,1], MDS2 = sppBC$points[,2],group=relCover2022$bison_ws)
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
  geom_point(size=6)+ 
  geom_path(data = BC_Ellipses, aes(x = NMDS1, y = NMDS2), size = 3) +
  labs(color="", linetype = "", shape = "") +
  scale_colour_manual(values=c("brown", "brown", "dark green", "dark green", "dark green", "dark green"), name = "") +
  scale_linetype_manual(values = c("twodash", "solid", "twodash", "solid", "twodash", "solid"), name = "") +
  xlab("NMDS1")+ 
  ylab("NMDS2")+ 
  theme(axis.text.x=element_text(size=24, color = "black"), axis.text.y = element_text(size = 24, color = "black"), legend.text = element_text(size = 24))



##### simper #####
print(sim <- with(relCover2022, simper(relCover2022[,7:172], trt)))


##### trends for dominant species #####
#bison effect
ggplot(barGraphStats(data=subset(relCover, genus_species=='andropogon_gerardii' & year>2018), variable="max_cover", byFactorNames=c("year", "bison", "watershed")), aes(x=year, y=mean, color=bison)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  ylab('Andropogon gerardii cover') +
  facet_wrap(~watershed)

#bison effect in 4 yr
ggplot(barGraphStats(data=subset(relCover, genus_species=='bouteloua_curtipendula' & year>2018), variable="max_cover", byFactorNames=c("year", "bison", "watershed")), aes(x=year, y=mean, color=bison)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  ylab('Bouteloua curtipendula cover') +
  facet_wrap(~watershed)

#invert effect
ggplot(barGraphStats(data=subset(relCover, genus_species=='lespedeza_violacea' & year>2018), variable="max_cover", byFactorNames=c("year", "invertebrates", "watershed")), aes(x=year, y=mean, color=invertebrates)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  ylab('Lespedeza violacea cover') +
  facet_wrap(~watershed)

#bison effect
ggplot(barGraphStats(data=subset(relCover, genus_species=='ambrosia_psilostachya' & year>2018), variable="max_cover", byFactorNames=c("year", "bison", "watershed")), aes(x=year, y=mean, color=bison)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  ylab('Ambrosia psilostachya cover') +
  facet_wrap(~watershed)

# ??? who is driving this sp?
ggplot(barGraphStats(data=subset(relCover, genus_species=='oxalis_violacea' & year>2018), variable="max_cover", byFactorNames=c("year", "trt", "watershed")), aes(x=year, y=mean, color=trt)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula=y~poly(x,2), se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  ylab('Oxalis violacea cover') +
  facet_wrap(~watershed)



#RACs
rankAbundance <- relCover %>%
  filter(year==2021) %>%
  mutate(spp_name=str_to_sentence(paste(genus, species, sep=' '))) %>%
  group_by(watershed, bison, spp_name, growthform, lifeform) %>%
  summarize(avg_cover=mean(rel_cover)) %>%
  ungroup() %>%
  arrange(bison, watershed, -avg_cover) %>%
  mutate(bison_ws=paste(bison, watershed, sep='::')) %>%
  group_by(bison_ws) %>%
  mutate(rank=seq_along(bison_ws)) %>%
  ungroup() %>%
  mutate(lifeform2=ifelse(spp_name=='Sisyrinchium campestre', 'f', ifelse(lifeform=='o', 'f', ifelse(lifeform=='s', 'g', as.character(lifeform))))) %>%
  filter(spp_name!='NA NA')

ggplot(data=subset(rankAbundance, bison_ws=='B::N1A', avg_cover>0), aes(x=rank, y=avg_cover)) +
  geom_line() +
  geom_point(aes(colour=lifeform2, shape=growthform), size=3) +
  scale_color_manual(labels=c("Forb", "Graminoid", "Woody"), values=c('#DE8C00', '#009CC0', '#83431E')) +
  scale_shape_discrete(labels=c("Annual", "Biennial", "Perennial"))+
  xlab('') +
  ylab('N1A Bison\nRelative Percent Cover\n') +
  # scale_x_continuous(expand=c(0,0), limits=c(0.5,17), breaks=seq(0,17,5)) +
  # scale_y_continuous(expand=c(0,0), limits=c(0,60), breaks=seq(0,60,10)) +
  # geom_text(aes(y=avg_cover+1.2, x=rank+0.1, label=spp_name), hjust='left', vjust='center', angle=90, size=4) +
  expand_limits(y=30, x=90)
#export at 1400x400

ggplot(data=subset(rankAbundance, bison_ws=='X::N1A', avg_cover>0), aes(x=rank, y=avg_cover)) +
  geom_line() +
  geom_point(aes(colour=lifeform2, shape=growthform), size=3) +
  scale_color_manual(labels=c("Forb", "Graminoid", "Woody"), values=c('#DE8C00', '#009CC0', '#83431E')) +
  scale_shape_discrete(labels=c("Annual", "Biennial", "Perennial"))+
  xlab('') +
  ylab('N1A Bison Removed\nRelative Percent Cover\n') +
  # scale_x_continuous(expand=c(0,0), limits=c(0.5,17), breaks=seq(0,17,5)) +
  # scale_y_continuous(expand=c(0,0), limits=c(0,60), breaks=seq(0,60,10)) +
  # geom_text(aes(y=avg_cover+1.2, x=rank+0.1, label=spp_name), hjust='left', vjust='center', angle=90, size=4) +
  expand_limits(y=30, x=90)
#export at 1400x400

ggplot(data=subset(rankAbundance, bison_ws=='B::N4B', avg_cover>0), aes(x=rank, y=avg_cover)) +
  geom_line() +
  geom_point(aes(colour=lifeform2, shape=growthform), size=3) +
  scale_color_manual(labels=c("Forb", "Graminoid", "Woody"), values=c('#DE8C00', '#009CC0', '#83431E')) +
  scale_shape_discrete(labels=c("Annual", "Biennial", "Perennial"))+
  xlab('') +
  ylab('N4B Bison\nRelative Percent Cover\n') +
  # scale_x_continuous(expand=c(0,0), limits=c(0.5,17), breaks=seq(0,17,5)) +
  # scale_y_continuous(expand=c(0,0), limits=c(0,60), breaks=seq(0,60,10)) +
  # geom_text(aes(y=avg_cover+1.2, x=rank+0.1, label=spp_name), hjust='left', vjust='center', angle=90, size=4) +
  expand_limits(y=40, x=90)
#export at 1400x400

ggplot(data=subset(rankAbundance, bison_ws=='X::N4B', avg_cover>0), aes(x=rank, y=avg_cover)) +
  geom_line() +
  geom_point(aes(colour=lifeform2, shape=growthform), size=3) +
  scale_color_manual(labels=c("Forb", "Graminoid", "Woody"), values=c('#DE8C00', '#009CC0', '#83431E')) +
  scale_shape_discrete(labels=c("Annual", "Biennial", "Perennial"))+
  xlab('') +
  ylab('N4B Bison Removed\nRelative Percent Cover\n') +
  # scale_x_continuous(expand=c(0,0), limits=c(0.5,17), breaks=seq(0,17,5)) +
  # scale_y_continuous(expand=c(0,0), limits=c(0,60), breaks=seq(0,60,10)) +
  # geom_text(aes(y=avg_cover+1.2, x=rank+0.1, label=spp_name), hjust='left', vjust='center', angle=90, size=4) +
  expand_limits(y=40, x=90)
#export at 1400x400






















# #finding the most frequent species across all plots and years for updating the datasheets
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