library(tidyverse)

setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\konza projects\\conSME\\data\\species composition') #laptop
setwd('C:\\Users\\komatsuk\\Dropbox (Smithsonian)\\konza projects\\conSME\\data\\species composition') #desktop


###functions
`%notin%` <- Negate(`%in%`)

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