################################################################################
##  04 specieslvl_effects.R: how individual species responded to these drivers
##
##  Author: Allison Louthan
##  Date created: January 16, 2025
##  Last modified : April 22, 2025
################################################################################

# libraries etc--- 
rm(list = ls())
library("data.table")
library("dplyr")
library("tidyr")
user <- "AL"


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

# fitting the models for each species----
sp_list <- unique(all_transitions$genus_species)
complete_sep_problems <- # these are the species who had all zeros in one or more treatments. you'll have to usepenalized regression 
  glm_problems <- NULL # these are the species who threw problems on the glm.nb command in the loop; you will have to run these models seperately

coefficients_allsp <- array(NA, dim=c(length(sp_list), 2,12 ), dimnames= list(sp_list, c("coefficient", "P-val"),
   c("(Intercept)", "watershedN4B"   ,     "trtBSI"      ,        "trtBSX"    ,          "trtXSI"   ,           "trtXSX"    ,         
 "trtXXI"        ,      "watershedN4B:trtBSI", "watershedN4B:trtBSX", "watershedN4B:trtXSI", "watershedN4B:trtXSX" ,"watershedN4B:trtXXI"
)))


for (i in 1:length(sp_list)){
  options (warn=0)
data_i <- all_transitions[which(all_transitions$genus_species== sp_list[i]), ] 
if (dim(data_i[which(data_i$max_cover >0), ])[1]<10){                             
  coefficients_allsp[i, 1, ] <- NA
  coefficients_allsp[i, 2, ] <- NA
  } else {
# NB for consistency with Kims analysis I am regressing cover (not change in cover!) against predictors
# and using her same error structure
    
    data_i$trt<- factor(data_i$trt, levels= c( "XXX", "XXI", "XSX", "BSX", "XSI", "BSI")) %>% droplevels()
  data_i$ year <- as.factor(data_i$year)
  data_i$ block <- as.factor(data_i$block)
  options (warn=2) # converts warnings to errors in code below

  if (var(data_i$max_cover)/ mean(data_i$max_cover) >1){
my_mod <- try( # try command allows us to try to fit the model, but allow for code to continue if it fails
  lme4::glmer.nb(max_cover~watershed*trt + (1|block) + (1|year), 
                      #  (1|block/trt), # ideally, you would have this RE structure, plus year nested in there somehow, but the model does not converge even for Andropogon ger
                                data=data_i ), 
  silent=TRUE) # suppress printing of errors (and warnings converted to errors)
options (warn=2)
if (class(my_mod)[1]== "try-error" || lme4::isSingular(my_mod)) { # if the most-complex model is a try-error, is singular or DN converge
   my_mod <- try(lme4::glmer.nb(max_cover~watershed*trt + year + (1|block) , data=data_i ), silent= TRUE) # then fit a simpler model
if (class(my_mod)[1]== "try-error" || lme4::isSingular(my_mod)) {  # if that model is singular or DN converge
  my_mod <- try(
    pkgcond::suppress_warnings( # if the model throws a warning of "step size truncated", we don't want that warning converted to an error, b/c its not fatal & the model is still valid. only applies to GLM's
    MASS::glm.nb(max_cover~watershed*trt + year + block , data=data_i , control=glm.control(maxit=2500)),  pattern = "step size truncated due to divergence", class = "warning") ,
    silent= TRUE) # use an even-simpler model-- take out RE. NB that for all glm models, complete separation sometimes occurs. if it does, we take out the predictor variable that's causing complete separation (one of watershed, block, or year-- never treatment-- b/c complete separation reflects that a species is not present in a block, watershed, or year, for example)
  if (class(my_mod)[1] == "try-error"){# if that model DN fit 
    my_mod <- try(
      pkgcond::suppress_warnings(
        MASS::glm.nb(max_cover~watershed*trt + block , data=data_i, control=glm.control(maxit=2500)),  pattern = "step size truncated due to divergence", class = "warning") ,
      silent= TRUE)# use an even-simpler model that removes year effect
   if (class(my_mod)[1] == "try-error"){# if that model DN fit 
    try(
      pkgcond::suppress_warnings(
        my_mod <- MASS::glm.nb(max_cover~watershed*trt  , data=data_i, control=glm.control(maxit=2500)),  pattern = "step size truncated due to divergence", class = "warning") , 
      silent= TRUE)
    # use an even-simpler model that removes year & block effects. code will add i to glm.nb problems vector if this model failts to fit
     if (class(my_mod)[1] == "try-error"){# if that model DN fit 
       my_mod <- try(
         pkgcond::suppress_warnings(
           MASS::glm.nb(max_cover~watershed+trt  , data=data_i, control=glm.control(maxit=2500)),  
           pattern = "step size truncated due to divergence", class = "warning")  ,
         silent= TRUE)}# use an even-simpler model that removes interaction term
     if (class(my_mod)[1] == "try-error"){# if that model DN fit 
       print(i) ### take out
       my_mod <- try(
         pkgcond::suppress_warnings(
           MASS::glm.nb(max_cover~trt, data=data_i, control=glm.control(maxit=2500)) , 
           pattern = "step size truncated due to divergence", class = "warning") , 
         silent= FALSE)}# use an even-simpler model in case the species is only present in one watershed
  }} }} # ends the outer try-error loop
options (warn=0) # make sure warnings stay as warnings 
  } else { # this else is var/ mean <1
    
    my_mod <- try(lme4::glmer(max_cover~watershed*trt + (1|block) + (1|year), 
                                 #  (1|block/trt), # ideally, you would have this RE structure, plus year nested in there somehow, but the model does not converge even for Andropogon ger
                                 data=data_i, family= "poisson" ), silent=TRUE)
    options (warn=2)
    if (class(my_mod)[1]== "try-error" || lme4::isSingular(my_mod)) { # if the most-complex model is a try-error, is singular or DN converge
      my_mod <- try(lme4::glmer(max_cover~watershed*trt + year + (1|block) , data=data_i, family= "poisson" ), silent= TRUE) # then fit a simpler model
      if (class(my_mod)[1]== "try-error" || lme4::isSingular(my_mod)) {  # if that model is singular or DN converge
        my_mod <- try(
          pkgcond::suppress_warnings(
            glm(max_cover~watershed*trt + year + block , data=data_i, family= "poisson" ),  pattern = "step size truncated due to divergence", class = "warning") , silent= TRUE) # use an even-simpler model-- take out RE
        if (class(my_mod)[1] == "try-error"){# if that model DN fit 
          my_mod <- try(
            pkgcond::suppress_warnings(
              glm(max_cover~watershed*trt + block , data=data_i, family= "poisson") ,  pattern = "step size truncated due to divergence", class = "warning") , silent= TRUE)# use an even-simpler model that removes year effect
          if (class(my_mod)[1] == "try-error"){# if that model DN fit 
            try(
              pkgcond::suppress_warnings(
                my_mod <- glm(max_cover~watershed*trt  , data=data_i, family= "poisson"),  pattern = "step size truncated due to divergence", class = "warning") , silent= TRUE)
            # use an even-simpler model that removes year & block effects. code will add i to glm.nb problems vector if this model failts to fit
            if (class(my_mod)[1] == "try-error"){# if that model DN fit 
              my_mod <- try(
                pkgcond::suppress_warnings(
                  glm(max_cover~watershed+trt, data=data_i, family= "poisson"),  pattern = "step size truncated due to divergence", class = "warning")  , silent= TRUE)}# use an even-simpler model that removes interaction term
            if (class(my_mod)[1] == "try-error"){# if that model DN fit 
              print(i) ### take out
              my_mod <- try(
                pkgcond::suppress_warnings(
                  glm(max_cover~trt, data=data_i, family= "poisson"),  pattern = "step size truncated due to divergence", class = "warning") , 
                silent= FALSE)}# use an even-simpler model in case the species is only present in one watershed
          }} }} # ends the outer try-error loop
    options (warn=0) # make sure warnings stay as warnings 
  } 
  
if (class(my_mod)[1] == "try-error") {
 # if the error message contains 0/1 separation term, then what has happened is there has been complete separation along the #trt condition-- you cannot quantify trt effects using a glm approach so you'll need to use a penalized regression approach
  if (grepl( "fitted rates numerically 0 or 1", attr(my_mod, "condition") , fixed= TRUE)) { # if the error message in the try-error object indicates complete separation
    complete_sep_problems <- c(complete_sep_problems, i)} else {
  glm_problems <- c(glm_problems, i)}
  } else {
 
my_coefs <- summary(my_mod)$coefficient[, "Estimate"]
names(my_coefs)[which(names(my_coefs) %in% c("trtBSI:watershedN4B", "trtBSX:watershedN4B", "trtXSI:watershedN4B", "trtXSX:watershedN4B", "trtXXI:watershedN4B"  ))] <- c("watershedN4B:trtBSI","watershedN4B:trtBSX", "watershedN4B:trtXSI", "watershedN4B:trtXSX", "watershedN4B:trtXXI"  ) # if the interaction term was written in the opposite order in the coefficient list, replace with standardized order
  coef_present_wNA <- match(names(my_coefs), dimnames(coefficients_allsp)[[3]] )# which coefficients are actually present, put them in the correct order-- some might be missing
  coef_present <- na.omit(match(names(my_coefs), dimnames(coefficients_allsp)[[3]] ))# which coefficients are actually present, put them in the correct order-- some might be missing
  
  coefficients_allsp[i, "coefficient",coef_present ] <-  my_coefs[which(!is.na(coef_present_wNA))]
  coefficients_allsp[i, "P-val",coef_present ] <- summary(my_mod)$coefficients[which(!is.na(coef_present_wNA)),"Pr(>|z|)"]} # ends the loop that says: if nothing works to fit this model, then add to problems list
} # ends if dim<260 loop
print(i)
} # can safely ignore non-convergence or singularity warnings; they are handled within the loop; the loop will break if any singularity or convergene issues are problematic

# figure----- 

#75/195 species did not have enough data to fit the model. 

# for the remaining 120 species, which cofficients are significant and what is their sign?

coef_names <- dimnames(coefficients_allsp)[[3]][-1]
coef_summary <- as.data.frame(coef_names)
coef_summary$no.sig <- NA
coef_summary$no.sig.pos <- NA
which_sig <- as.data.frame(matrix(NA, nrow= 100, ncol= length(coef_names)))
names(which_sig) <- coef_names
which_sig_pos <- as.data.frame(matrix(NA, nrow= 100, ncol= length(coef_names)))
names(which_sig_pos) <- coef_names
  
for (i in 1:length(coef_names)){
  
  coef_summary$no.sig[i] <- length(which(coefficients_allsp[,"P-val", coef_names[i]] <= 0.05)) 
  names_sig_i <- names(which(coefficients_allsp[,"P-val", coef_names[i]] <= 0.05))
  names_sig_neg_i <- names(which(coefficients_allsp[,"P-val", coef_names[i]] <= 0.05 & coefficients_allsp[,"coefficient", coef_names[i]]<0))
  nos_sig_neg_i <- which(names_sig_i %in% names_sig_neg_i)
  names_sig_i[nos_sig_neg_i] <- paste(names_sig_i[nos_sig_neg_i] , "*", sep= "") # adds an asterisk to those species whose coefficients are negative 
  names_sig_i <- sub("_", " ",names_sig_i )
  which_sig[1:length(names_sig_i),i] <- sort(stringr::str_to_sentence(names_sig_i))
  coef_summary$no.sig.pos[i] <- length(which(coefficients_allsp[,"P-val", coef_names[i]] <= 0.05 & coefficients_allsp[,"coefficient", coef_names[i]]>0))

  
  }

species_names <- dimnames(coefficients_allsp)[[1]]
species_summary <- as.data.frame((species_names))
species_summary$no.sig <- NA

for (i in 1:length(species_list)){
  if (all(is.na(coefficients_allsp[i,"P-val",]))) { # if all the p-vals are NA, that means that species model wasn't fit
    species_summary$no.sig[i] <- NA} else {
      species_summary$no.sig[i] <-  sum(coefficients_allsp[i,"P-val", 3:12]<=0.05, na.rm=TRUE) # how many sig trmt effects did that species have?
    } }
species_summary$any.sig <- NA
species_summary$any.sig[which(species_summary$no.sig>0)] <- 1
species_summary$any.sig[which(species_summary$no.sig==0)] <- 0


write.csv(coef_summary,file= "derived_data/04_coef_sig_summary.csv" )
write.csv(which_sig, file= "derived_data/04_coef_sig_species_names.csv", na="")
write.csv(species_summary,file= "derived_data/04_species_sig_summary.csv" )

# summary statistics for paper---- 
sum(species_summary$any.sig, na.rm=TRUE) # how many species had any significant treatment effects
sig.species <- species_summary[which(species_summary$any.sig== 1),]
sp_unique <- unique(spAll[, c("genus_species","family"   ,     "growthform"  ,  "lifeform"   ,   "origin"  )])
sig.species <- left_join(sig.species, sp_unique, by= c("X.species_names."= "genus_species"))
sig.species.summary <- sig.species %>% count(lifeform) # how many species had sig effects of each growth form
nsig.species <- species_summary[which(species_summary$any.sig== 0),] # which species did NOT have sig effects-- given that we could fit a treatment effect
nsig.species <- left_join(nsig.species, sp_unique, by= c("X.species_names."= "genus_species"))
nsig.species.summary <- nsig.species %>% count(lifeform) # how many species had sig effects of each growth form
cont.table <- left_join(nsig.species.summary, sig.species.summary, by= "lifeform")
colnames(cont.table) <- c("lifeform", "no.nsig", "no.sig")
rownames(cont.table) <- cont.table$lifeform
cont.table$no.sig[which(is.na(cont.table$no.sig))] <- 0

chisq.test(cont.table[,-1], simulate.p.value = TRUE)

coef_summary
