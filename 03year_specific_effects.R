################################################################################
##  03year_specific_effects.R: tests whether annual variation in herbivore effect size differs with rainfall
##
##  Author: Allison Louthan
##  Date created: January 16, 2025
##  Modified : Feb 26, 2025 
################################################################################

# ok, this whole code is meaningless because the year effect should not be a factor. 
# load in data----
rm(list = ls())
load("derived_data/01biomass_models.RData") # biomass models from Kim's conSME_biomass code
load('derived_data/02comp_models.RData') # composition models from Kim's community metrics code

# temp data-- daily
clim_dat_raw <- read_csv("./data/AWE012.csv", col_types = cols(.default = "?",  TMIN = "n",TAVE = "n",DHUMID = "n",
                                                               DSRAD = "n",DPPT = "n",SMAX = "n",SMIN = "n",S_AVE = "n",
                                                               WAVE = "n"), na = c(".", "NA")) %>%  rename_with(tolower) # standardize col name formats

# Konza HQ precipitation. More accurate than the AWE01 tipping bucket method (per the AWE01 metadata)
precip_dat_raw <- read_csv("./data/APT011.csv", col_types = c(.default = "?",
                                                              ppt = "n"), na = c(".", "NA")) # These should be replaced as NA when converting to numeric


# calc annual weather----

# Sarah's code to make SPEI: 

# get drought index
# average climate data by month
clim_by_month <- clim_dat_raw %>% 
  group_by(recyear, recmonth) %>% 
  summarise(#pr_mean = mean(dppt, na.rm = T),
    #tmax_mean = mean(tmax, na.rm = T),
    #tmin_mean = mean(tmin, na.rm = T),
    tmean_mean = mean(tave, na.rm = T)) %>% # temp means
filter(recyear %in% 2015:2023)  # trim data
clim_by_month <- 
  precip_dat_raw %>%  # precip means
  mutate(RecDate = lubridate::ymd(RecDate), # convert to date format
         recyear = year(RecDate),
         recmonth = month(RecDate)) %>% 
  group_by(recyear, recmonth) %>% 
  summarise(pr_mean = mean(ppt, na.rm = T)) %>% 
  ungroup() %>% 
  add_row(recyear = c(1986,1989, 2009, 2016), recmonth = c(1, 11, 1, 1), pr_mean = 0) %>%   # adding some dummy data for months when precip was too low to register
  filter(recyear %in% 2015:2023) %>%  # trim data
  full_join(clim_by_month) %>% # add precip and temp dat together
  
  arrange(recyear) #cosmetic

# calculate potential evapotranspiration (PET)
clim_by_month$pet <- SPEI::thornthwaite(Tave = clim_by_month$tmean_mean, lat = 39.09306) # konza latitude
# calculate climatic water balance (BAL)
clim_by_month$bal <- clim_by_month$pr_mean - clim_by_month$pet
# convert to time series (ts) for convenience
clim_by_month <- ts(clim_by_month[, -c(1, 2)], end = c(2023, 12), frequency = 12) #end = c(2020, 12)
plot(clim_by_month)
# three and twelve-months SPEI
spei3 <- SPEI::spei(clim_by_month[, "bal"], 3, na.rm = T) 
spei12 <- SPEI::spei(clim_by_month[, "bal"], 12, na.rm = T)

# join
spei_res <- full_join(data.frame(spei3.z=as.matrix(spei3$fitted), date=zoo::as.Date(time(spei3$fitted))),
                      data.frame(spei12.z=as.matrix(spei12$fitted), date=zoo::as.Date(time(spei12$fitted))))
spei_res <- spei_res[which(lubridate::month(spei_res$date)==8), ] # here, the spei values correspond to the spei over the n months prior to Aug of th eindicated year
spei_res$year <- lubridate::year(spei_res$date)
spei_res$yearprior <- spei_res$year-1

# # SPEI values between -1 and 1 are considered near normal for a given area, whereas values below -1 signify drought and values above 1 signify unusually moist conditions. Because drought conditions fluctuate naturally, it is helpful to look at average conditions over several years to explore how drought is connected to long-term climate change.5


# biomass----
# reornanzing the coefficients form totBiomassModel into year-specific main effect and interaction terms
intrxt_terms1 <- stack(totBiomassModel$coefficients$fixed[grepl("as.factor(year)1" , names(totBiomassModel$coefficients$fixed), fixed=TRUE)] )
intrxt_terms1$ind <- gsub("\\s*\\([^\\)]+\\)", "", intrxt_terms1$ind)
intrxt_terms1$ind <- gsub("as.factor1", "year", intrxt_terms1$ind)
intrxt_terms1$expt_year <- 2019 
  intrxt_terms2 <- stack(totBiomassModel$coefficients$fixed[grepl("as.factor(year)2" , names(totBiomassModel$coefficients$fixed), fixed=TRUE)] )
intrxt_terms2$ind <- gsub("\\s*\\([^\\)]+\\)", "", intrxt_terms2$ind)
intrxt_terms2$ind <- gsub("as.factor2", "year", intrxt_terms2$ind)
intrxt_terms2$expt_year <- 2020 
intrxt_terms3 <- stack(totBiomassModel$coefficients$fixed[grepl("as.factor(year)3" , names(totBiomassModel$coefficients$fixed), fixed=TRUE)]) # extract the terms that include a year effect
intrxt_terms3$ind <- gsub("\\s*\\([^\\)]+\\)", "", intrxt_terms3$ind)
intrxt_terms3$ind <- gsub("as.factor3", "year", intrxt_terms3$ind)
intrxt_terms3$expt_year <- 2021 
intrxt_terms4 <- as.data.frame(rep(0, 12))
names(intrxt_terms4) <- "values"
intrxt_terms4$ind <- intrxt_terms3$ind
intrxt_terms4$expt_year <- 2022

intrxt_terms <- rbind(
  cbind(pivot_wider(intrxt_terms1, values_from = values, names_from = c(ind)), year = 1), 
cbind(pivot_wider(intrxt_terms2, values_from = values, names_from = c(ind)), year = 2), 
cbind(pivot_wider(intrxt_terms3, values_from = values, names_from = c(ind)), year = 3), 
cbind(pivot_wider(intrxt_terms4, values_from = values, names_from = c(ind)), year = 4)) 

# combine climate data with interaction term data---
# first metric of annual weather: SPEI over the entire year prior
intrxt_termsSPEI <-intrxt_terms %>% select(-year) %>%
left_join(spei_res, by= join_by(expt_year== year),) %>%
  rename(last3SPEI = spei3.z, last12SPEI = spei12.z) %>%
  left_join(spei_res[,c("yearprior","spei3.z" )], by= join_by(expt_year== yearprior)) %>%
  rename(priorgrowingseasonSPEI= spei3.z)
# (NB this last12SPEI and last3SPEI can include weather after plant data were collected; e.g., biomass data were collected in June but 
# this metric includes weather from JUly)

# what I would do is say: bison x year and small mammal x year and insect x year effects are strongest (i.e.,more negative) in low SPEI years, 
# so compare across the differnet herbivore guilds I guess
# but what about the interaction term between herbivores... bot sure what to do there


covariance_intrxt_terms_SPEI <- as.data.frame(cov(select(intrxt_termsSPEI, -c(date, expt_year, yearprior))))
covariance_intrxt_terms_SPEI <- select(covariance_intrxt_terms_SPEI, c(last3SPEI, last12SPEI, priorgrowingseasonSPEI))
covariance_intrxt_terms_SPEI <- covariance_intrxt_terms_SPEI[-which(rownames(covariance_intrxt_terms_SPEI) %in% c("watershed1:year", "last3SPEI"   ,   "priorgrowingseasonSPEI")), ]  
# what to add to paper:
# we watned to understand why the effect of bison, small mammals differed across years.
# here is what these numbers in covariance_intrxt_terms_SPEI mean: priorgrowingseasonSPEI is positively correlated with most year-specific interaction terms that involve herbivores
# in both watersheds!
# by contrast, last12SPEI is sometimes positively, sometimes negatively correlated with year-specific interaction terms
# we can't do formal statistical tests on these results becuase each correlation is only represnted by 4 data points (and because they are non-sig-- see below)
# thus, when SPEI is high, the negative effect of herbivores on biomass is weaker-- particualrly for small mammals AND for bison
# and, most other interaction terms between two herbivore guilds and year are also positive, 
# , such that presence of one herbivore tends to weaken the effect of other herbivores' presence in high SPEI years
# importantly, there is a non-sig correlation between year of experiment and any SPEI metric, indicating we can 
# disentangle climate year effects (though note that there are only 4 data points here) (p val for prior growing season SPEI is 0.7681, for last 12 SPEI is 0.1351, for last 3 SPEI is 0.3622)
# NB that only some of the year x herbivore treatmetn effects were significatn according to an ANOVA 

year_clim_corrs <- intrxt_termsSPEI %>% select(c(expt_year, last3SPEI, last12SPEI,priorgrowingseasonSPEI ))
cor.test(year_clim_corrs$expt_year, year_clim_corrs$last3SPEI)
cor.test(year_clim_corrs$expt_year, year_clim_corrs$last12SPEI)
cor.test(year_clim_corrs$expt_year, year_clim_corrs$priorgrowingseasonSPEI)

warning("repeat this analysis for grass and forb biomass because grass biomass has sig year x herb effects but forb and woody biomass do not")

# richness----
# reornanzing the coefficients form totBiomassModel into year-specific main effect and interaction terms
intrxt_terms1 <- stack(richModel$coefficients$fixed[grepl("as.factor(year)1" , names(richModel$coefficients$fixed), fixed=TRUE)] )
intrxt_terms1$ind <- gsub("\\s*\\([^\\)]+\\)", "", intrxt_terms1$ind)
intrxt_terms1$ind <- gsub("as.factor1", "year", intrxt_terms1$ind)
intrxt_terms1$expt_year <- 2019 
intrxt_terms2 <- stack(richModel$coefficients$fixed[grepl("as.factor(year)2" , names(richModel$coefficients$fixed), fixed=TRUE)] )
intrxt_terms2$ind <- gsub("\\s*\\([^\\)]+\\)", "", intrxt_terms2$ind)
intrxt_terms2$ind <- gsub("as.factor2", "year", intrxt_terms2$ind)
intrxt_terms2$expt_year <- 2020 
intrxt_terms3 <- stack(richModel$coefficients$fixed[grepl("as.factor(year)3" , names(richModel$coefficients$fixed), fixed=TRUE)]) # extract the terms that include a year effect
intrxt_terms3$ind <- gsub("\\s*\\([^\\)]+\\)", "", intrxt_terms3$ind)
intrxt_terms3$ind <- gsub("as.factor3", "year", intrxt_terms3$ind)
intrxt_terms3$expt_year <- 2021 
intrxt_terms4 <- as.data.frame(rep(0, 12))
names(intrxt_terms4) <- "values"
intrxt_terms4$ind <- intrxt_terms3$ind
intrxt_terms4$expt_year <- 2022

intrxt_terms <- rbind(
  cbind(pivot_wider(intrxt_terms1, values_from = values, names_from = c(ind)), year = 1), 
  cbind(pivot_wider(intrxt_terms2, values_from = values, names_from = c(ind)), year = 2), 
  cbind(pivot_wider(intrxt_terms3, values_from = values, names_from = c(ind)), year = 3), 
  cbind(pivot_wider(intrxt_terms4, values_from = values, names_from = c(ind)), year = 4)) 

# combine climate data with interaction term data
# first metric of annual weather: SPEI over the entire year prior
intrxt_termsSPEI <-intrxt_terms %>% select(-year) %>%
  left_join(spei_res, by= join_by(expt_year== year),) %>%
  rename(last3SPEI = spei3.z, last12SPEI = spei12.z) %>%
  left_join(spei_res[,c("yearprior","spei3.z" )], by= join_by(expt_year== yearprior)) %>%
  rename(priorgrowingseasonSPEI= spei3.z)
# (NB this last12SPEI and last3SPEI can include weather after plant data were collected; e.g., biomass data were collected in June but 
# this metric includes weather from JUly)


covariance_intrxt_terms_SPEI <- as.data.frame(cov(select(intrxt_termsSPEI, -c(date, expt_year, yearprior))))
covariance_intrxt_terms_SPEI <- select(covariance_intrxt_terms_SPEI, c(last3SPEI, last12SPEI, priorgrowingseasonSPEI))
covariance_intrxt_terms_SPEI <- covariance_intrxt_terms_SPEI[-which(rownames(covariance_intrxt_terms_SPEI) %in% c("watershed1:year", "last3SPEI"   ,   "priorgrowingseasonSPEI")), ]  

# here is what these numbers in covariance_intrxt_terms_SPEI mean: because the sign of the correlation between annual climate and year x herbivore interaction terms is inconsistent
# we watned to understand why the effect of bison, small mammals differed across years. (i.e., why there was a significatn year x bison and year x small mammal terms)
# it does not seem like the (sig?) year interaction term on richness is driven by annual variation in climate effects 

# composition/ evar----
# reornanzing the coefficients form totBiomassModel into year-specific main effect and interaction terms
intrxt_terms1 <- stack(evarModel$coefficients$fixed[grepl("as.factor(year)1" , names(evarModel$coefficients$fixed), fixed=TRUE)] )
intrxt_terms1$ind <- gsub("\\s*\\([^\\)]+\\)", "", intrxt_terms1$ind)
intrxt_terms1$ind <- gsub("as.factor1", "year", intrxt_terms1$ind)
intrxt_terms1$expt_year <- 2019 
intrxt_terms2 <- stack(evarModel$coefficients$fixed[grepl("as.factor(year)2" , names(evarModel$coefficients$fixed), fixed=TRUE)] )
intrxt_terms2$ind <- gsub("\\s*\\([^\\)]+\\)", "", intrxt_terms2$ind)
intrxt_terms2$ind <- gsub("as.factor2", "year", intrxt_terms2$ind)
intrxt_terms2$expt_year <- 2020 
intrxt_terms3 <- stack(evarModel$coefficients$fixed[grepl("as.factor(year)3" , names(evarModel$coefficients$fixed), fixed=TRUE)]) # extract the terms that include a year effect
intrxt_terms3$ind <- gsub("\\s*\\([^\\)]+\\)", "", intrxt_terms3$ind)
intrxt_terms3$ind <- gsub("as.factor3", "year", intrxt_terms3$ind)
intrxt_terms3$expt_year <- 2021 
intrxt_terms4 <- as.data.frame(rep(0, 12))
names(intrxt_terms4) <- "values"
intrxt_terms4$ind <- intrxt_terms3$ind
intrxt_terms4$expt_year <- 2022

intrxt_terms <- rbind(
  cbind(pivot_wider(intrxt_terms1, values_from = values, names_from = c(ind)), year = 1), 
  cbind(pivot_wider(intrxt_terms2, values_from = values, names_from = c(ind)), year = 2), 
  cbind(pivot_wider(intrxt_terms3, values_from = values, names_from = c(ind)), year = 3), 
  cbind(pivot_wider(intrxt_terms4, values_from = values, names_from = c(ind)), year = 4)) 

# combine climate data with interaction term data
# first metric of annual weather: SPEI over the entire year prior
intrxt_termsSPEI <-intrxt_terms %>% select(-year) %>%
  left_join(spei_res, by= join_by(expt_year== year),) %>%
  rename(last3SPEI = spei3.z, last12SPEI = spei12.z) %>%
  left_join(spei_res[,c("yearprior","spei3.z" )], by= join_by(expt_year== yearprior)) %>%
  rename(priorgrowingseasonSPEI= spei3.z)
# (NB this last12SPEI and last3SPEI can include weather after plant data were collected; e.g., biomass data were collected in June but 
# this metric includes weather from JUly)


covariance_intrxt_terms_SPEI <- as.data.frame(cov(select(intrxt_termsSPEI, -c(date, expt_year, yearprior))))
covariance_intrxt_terms_SPEI <- select(covariance_intrxt_terms_SPEI, c(last3SPEI, last12SPEI, priorgrowingseasonSPEI))
covariance_intrxt_terms_SPEI <- covariance_intrxt_terms_SPEI[-which(rownames(covariance_intrxt_terms_SPEI) %in% c("watershed1:year", "last3SPEI"   ,   "priorgrowingseasonSPEI")), ]  
# we watned to understand why the effect of bison differed across years. (i.e., why there was a significatn year x bison term)

# here is what these numbers in covariance_intrxt_terms_SPEI mean: because the sign of the correlation between annual climate and year x herbivore interaction terms is inconsistent
# it does not seem like the (sig?) year interaction term on compoisition is driven by annual variation in climate effects 


