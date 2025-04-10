################################################################################
##  04 specieslvl_effects.R: how individual species responded to these drivers
##
##  Author: Allison Louthan
##  Date created: January 16, 2025
##  Last modified : March 11, 2025
################################################################################

# libraries etc--- 
rm(list = ls())
library("data.table")
library("dplyr")
library("tidyr")
user <- "AL"

# NB that the performance::check_singularity command doesn't appear to be working. next steps: chagne the glmer_converged function to both glmer converged and singular, then replace all uses of the performance library
glmer_converged <- function(model_function){
  !any(grepl("failed to converge", model_function@optinfo$conv$lme4$messages))
}


# data preparation----
if (user == "AL"){
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
    filter(!(gen %in% c('litter', 'rock', 'dung', 'bare_ground', 'bison_trail'))) %>%
    mutate(genus_species=paste(genus, species, sep='_'))
  trt <- read.csv('data/conSME_treatments.csv')
}


# add zeros if that species was not included in that year x watershed x block x plot combination:
species_list <- unique(spAll$genus_species)
plotyear_list <- as.vector(t(unite(as.data.frame(unique(cbind(spAll$watershed, spAll$block, spAll$plot, spAll$year ))), plotyear_list)))

for (i in 1:length(species_list)){
  spAll_i <- spAll[which(spAll$genus_species== species_list[i]), ]
  
  missing_plot_years <- plotyear_list[
    which(!(plotyear_list %in% as.vector(t(unite(as.data.frame(cbind(spAll_i$watershed, spAll_i$block, spAll_i$plot, spAll_i$year )), plotyear_i )))))]
  zeros_data <- data.frame(matrix(data=NA, nrow=length(missing_plot_years), ncol= ncol(spAll)))
  names(zeros_data) <- names(spAll)
  zeros_data[,c("sppnum", "gen", "spec", "genus", "species", "family", "growthform","lifeform",  "origin", "photo", "genus_species" )] <- 
    spAll_i[1,c("sppnum", "gen", "spec", "genus", "species", "family", "growthform","lifeform",  "origin", "photo", "genus_species" )]
  zeros_data[,c("watershed", "block", "plot", "year")] <- stringr::str_split_fixed(missing_plot_years, "_", n= 4)
  zeros_data$max_cover <- 0
  spAll <- rbind(spAll, zeros_data)
  
}
# add on previous years' data
spAll <- spAll %>% mutate(prioryear = as.numeric(year)-1, numeric_year = as.numeric(year))
all_transitions <- spAll  %>% 
  left_join(spAll, by= join_by(prioryear== numeric_year, 
                               watershed, block, plot,sppnum,gen   , 
                               spec  ,genus     ,    species ,family ,growthform ,lifeform, origin ,photo ,genus_species), 
            suffix= c("", "_prioryear")) %>% select(-prioryear_prioryear, -year_prioryear)
all_transitions$year <- as.numeric(all_transitions$year)
all_transitions$plot <- as.numeric(all_transitions$plot)

all_transitions <-all_transitions %>% 
  left_join(trt) %>% 
mutate(ws_label=ifelse(watershed=='N1A', 'Annual', '4 Year'),
       experiment_year=year-2018)


sp_list <- unique(all_transitions$genus_species)
glm.nb_problems <- NULL # these are the species who threw problems on the glm.nb command in the loop; you will have to run these models seperately
sp_list_no_glm.nb_problems <- sp_list[!sp_list %in% glm.nb_problems]

coefficients_allsp <- array(NA, dim=c(length(sp_list), 2,12 ), dimnames= list(sp_list, c("coefficient", "P-val"),
   c("(Intercept)", "watershedN4B"   ,     "trtBSI"      ,        "trtBSX"    ,          "trtXSI"   ,           "trtXSX"    ,         
 "trtXXI"        ,      "watershedN4B:trtBSI", "watershedN4B:trtBSX", "watershedN4B:trtXSI", "watershedN4B:trtXSX" ,"watershedN4B:trtXXI"
)))

# start here at i= 90, there is a weird wraning there. 
for (i in 1:length(sp_list_no_glm.nb_problems)){
  options (warn=0)
data_i <- all_transitions[which(all_transitions$genus_species== sp_list_no_glm.nb_problems[i]), ] 
if (dim(data_i)[1]<260){                             
  coefficients_allsp[i, 1, ] <- NA
  coefficients_allsp[i, 2, ] <- NA} else {
# NB for consistency with Kims analysis I am regressing cover (not change in cover!) against predictors
# and using her same error structure
    data_i$trt <- as.factor(data_i$trt)
  data_i$trt <- relevel(data_i$trt, ref= "XXX")
  data_i$ year <- as.factor(data_i$year)
  data_i$ block <- as.factor(data_i$block)
my_mod <- try(lme4::glmer.nb(max_cover~watershed*trt + (1|block) + (1|year), 
                      #  (1|block/trt), # ideally, you would have this RE structure, plus year nested in there somehow, but the model does not converge even for Andropogon ger
                                data=data_i ), silent=TRUE)
if (class(my_mod)[1]== "try-error" || performance::check_singularity(my_mod) || !performance::check_convergence(my_mod)) { # if the most-complex model is a try-error, is singular or DN converge
  my_mod <- try(lme4::glmer.nb(max_cover~watershed*trt + year + (1|block) , data=data_i ), silent= TRUE) # then fit a simpler model
if (class(my_mod)[1]== "try-error" || performance::check_singularity(my_mod) || !performance::check_convergence(my_mod)) {  # if that model is singular or DN converge
  options(warn= 2) # make sure warnings are converted to errors, and then
  my_mod <- try(MASS::glm.nb(max_cover~watershed*trt + year + block , data=data_i ), silent= TRUE) # use an even-simpler model-- take out RE
  if (class(my_mod)[1] == "try-error"){# if that model DN fit 
    my_mod <- try(MASS::glm.nb(max_cover~watershed*trt + block , data=data_i) , silent= TRUE)# use an even-simpler model that removes year effect
   if (class(my_mod)[1] == "try-error"){# if that model DN fit 
    try(my_mod <- MASS::glm.nb(max_cover~watershed*trt  , data=data_i), silent= TRUE)
    # use an even-simpler model that removes year & block effects. code will add i to glm.nb problems vector if this model failts to fit
  }}
  options (warn=0) # make sure warnings stay as warnings
  }
  if (class(my_mod)[1] == "try-error") {glm.nb_problems <- c(glm.nb_problems, i)} else {
    
if(class(my_mod)[1] == "negbin") {my_coefs <- coef(my_mod)} else  {my_coefs <-  lme4::fixef(my_mod)} 
  
  coef_present_wNA <- match(names(my_coefs), dimnames(coefficients_allsp)[[3]] )# which coefficients are actually present, put them in the correct order-- some might be missing
  coef_present <- na.omit(match(names(my_coefs), dimnames(coefficients_allsp)[[3]] ))# which coefficients are actually present, put them in the correct order-- some might be missing
  
      coefficients_allsp[i, "coefficient",coef_present ] <-  my_coefs[which(!is.na(coef_present_wNA))]
  coefficients_allsp[i, "P-val",coef_present ] <- summary(my_mod)$coefficients[which(!is.na(coef_present_wNA)),"Pr(>|z|)"]}
}
                                                          



print(i)}
} # can safely ignore non-convergence or singularity warnings; they are handled within the loop; the loop will break if any singularity or convergene issues are problematic


# do the glm.nb problems individuallY: https://stackoverflow.com/questions/67360883/how-to-fix-fitted-probabilities-numerically-0-or-1-occurred-warning-in-r

glm.nb_problems
