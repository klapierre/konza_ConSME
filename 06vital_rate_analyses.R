##  06vital_rate_analyses.R: tests for herbivore effects on vital rates of 3 forbs
##
##  Author: Allison Louthan
##  Date created: April 14, 2025
##  Modified : 
################################################################################

start here, with this list of problems:
    1. even though I wrote the catch, some of the models are still try-error; WTF
    2. once you figure that out, fix this code, the 04 code, and email Johannes
    3. combine pf and cf into a neg bin response to help with convergence
    
rm(list = ls())

# functions and libraries----
library("lme4") #NB this is DIFFERENT than what Kim used 

# load in data----

amca_data <- read.csv("derived_data/05_amca_clean_demo_data.RData")
ecan_data <- read.csv("derived_data/05_ecan_clean_demo_data.RData")
kueu_data <- read.csv("derived_data/05_kueu_clean_demo_data.RData")

trt <- read.csv('data/conSME_treatments.csv')

# fitting vital rates----

sur <- gr <- pf <- cf <- list(NA, NA, NA)

for (i in 1:3){
  if (i == 1) data_i <- amca_data
  if (i == 2) data_i <- ecan_data
  if (i == 3) data_i <- kueu_data
  
  data_i <- left_join(data_i, trt)
  data_i$watershed <- factor(data_i$watershed) %>% droplevels()
  data_i$bison <- factor(data_i$bison) %>% droplevels()
  data_i$small_mammal <- factor(data_i$small_mammal) %>% droplevels()
  data_i$invertebrates <- factor( data_i$invertebrates) %>% droplevels()
  data_i$trt<- factor(data_i$trt) %>% droplevels()
  data_i$block <- factor(data_i$block) %>% droplevels()
  
  data_i$log_biom_1 <- log(data_i$biom_1)
  data_i$log_biom_2 <- log(data_i$biom_2)
  data_i$log_biom_1 <- log(data_i$biom_1)
  data_i$log_biom_2 <- log(data_i$biom_2)
  
  data_i$prf_2 <- NA
  data_i$prf_2[which(data_i$f_2 >0)] <- 1
  data_i$prf_2[which(data_i$f_2 ==0)] <- 0
  data_i$ log_cf_2 <- NA
  data_i$ log_cf_2[which(data_i$f_2 >0)] <- log(data_i$f_2[which(data_i$f_2 >0)])
  #   watershed*trt + (1|block) + (1|year) # for the counts of species
  #   watershed*as.factor(year)*invertebrates*bison + watershed*as.factor(year)*invertebrates*small_mammal # Kim's analysis 
  contrasts_present <- data_i %>% 
    group_by(watershed) %>%
    summarise(n=sum(!is.na(sur_1_2)))
  options(warn= 2)
  
  if (length(which(contrasts_present$n> 0))>1){ # if there are data in both watersheds, keep the watershed effect
# survival
    sur[[i]] <- try(
      lme4::glmer(sur_1_2 ~ log_biom_1 + watershed*trt  + (1|block)+ (1|startyear), data= data_i, family= "binomial")
                    , silent=TRUE)
    if (class(sur[[i]])[1]== "try-error" || lme4::isSingular(sur[[i]])) { # if the most-complex model is a try-error, is singular or DN converge
      sur[[i]] <- try(lme4::glmer(sur_1_2 ~ log_biom_1 + watershed*trt  + block + (1|startyear), data= data_i, family= "binomial"), silent= TRUE) # then fit a simpler model
      if (class(sur[[i]])[1] == "try-error" || lme4::isSingular(sur[[i]])) {  # if that model is singular or DN converge
        sur[[i]] <- try(glm(sur_1_2 ~ log_biom_1 + watershed*trt  + block + startyear, data= data_i, family= "binomial"), silent=TRUE)
        if (class(sur[[i]])[1] == "try-error") {print(paste("species", i, "sur DN converge"))}}}
    
# growth (mean only, no variatnce)
    gr[[i]] <-try(
      glmer(log_biom_2 ~ log_biom_1  + watershed*trt + (1|block)+ (1|startyear), data= data_i), 
      silent=TRUE)
    if (class(gr[[i]])[1]== "try-error" || lme4::isSingular(gr[[i]])) { # if the most-complex model is a try-error, is singular or DN converge
      gr[[i]] <- try(lme4::glmer(log_biom_2 ~ log_biom_1 + watershed*trt  + block + (1|startyear), data= data_i), silent= TRUE) # then fit a simpler model
      if (class(gr[[i]])[1] == "try-error" || lme4::isSingular(gr[[i]])) {  # if that model is singular or DN converge
        gr[[i]] <- try(glm(log_biom_1 ~ log_biom_2 + watershed*trt  + block + startyear, data= data_i), silent=TRUE)
        if (class(gr[[i]])[1] == "try-error") {print(paste("species", i, "gr DN converge"))}}}

  # prob of fruiting    
    pf[[i]] <-try(glmer(prf_2 ~ log_biom_1  + watershed*trt + (1|block)+ (1|startyear), data= data_i, family= "binomial"), silent=TRUE)
    if (class(pf[[i]])[1]== "try-error" || lme4::isSingular(pf[[i]])) { # if the most-complex model is a try-error, is singular or DN converge
      pf[[i]]  <- try(lme4::glmer(prf_2 ~ log_biom_1 + watershed*trt  + block + (1|startyear), data= data_i, family= "binomial"), silent= TRUE) # then fit a simpler model
      if (class(pf[[i]])[1] == "try-error" || lme4::isSingular(pf[[i]])) {  # if that model is singular or DN converge
        pf[[i]]  <- try(glm(prf_2 ~ log_biom_1 + watershed*trt  + block + startyear, data= data_i, family= "binomial"), silent=TRUE)
      
        if (class(pf[[i]])[1] == "try-error") {print(paste("species", i, "pf DN converge"))}}}
    
# number of fruits
    cf[[i]]<-try(glmer(log_cf_2 ~ log_biom_1 + watershed*trt + (1|block)+ (1|startyear), data= data_i)  , silent=TRUE) 
    if (class(cf[[i]])[1]== "try-error" || lme4::isSingular(cf[[i]])) { # if the most-complex model is a try-error, is singular or DN converge
      cf[[i]] <- try(lme4::glmer(log_cf_2 ~ log_biom_1 + watershed*trt  + block + (1|startyear), data= data_i), silent= TRUE) # then fit a simpler model
      if (class(cf[[i]])[1] == "try-error" || lme4::isSingular(cf[[i]])) {  # if that model is singular or DN converge
        cf[[i]] <- try(glm(log_cf_2 ~ log_biom_1 + watershed*trt  + block + startyear, data= data_i), silent=TRUE)
        if (class(cf[[i]])[1] == "try-error") {print(paste("species", i, "cf DN converge"))}}}
    
    
    } else { # if there is data in only one watershed, remove the watershed effect
      sur[[i]] <- try(
        lme4::glmer(sur_1_2 ~ log_biom_1 + trt  + (1|block)+ (1|startyear), data= data_i, family= "binomial")
        , silent=TRUE)
      if (class(sur[[i]])[1]== "try-error" || lme4::isSingular(sur[[i]])) { # if the most-complex model is a try-error, is singular or DN converge
        sur[[i]]  <- try(lme4::glmer(sur_1_2 ~ log_biom_1 + trt  + block + (1|startyear), data= data_i, family= "binomial"), silent= TRUE) # then fit a simpler model
        if (class(sur[[i]])[1] == "try-error" || lme4::isSingular(sur[[i]])) {  # if that model is singular or DN converge
          sur[[i]]  <- try(glm(sur_1_2 ~ log_biom_1 + trt  + block + startyear, data= data_i, family= "binomial"), silent=TRUE)
          if (class(sur[[i]])[1] == "try-error") {print(paste("species", i, "sur DN converge- trt only"))}}}
      
      # growth (mean only, no variatnce)
      gr[[i]] <-try(
        glmer(log_biom_2 ~ log_biom_1  + trt + (1|block)+ (1|startyear), data= data_i), 
        silent=TRUE)
      if (class(gr[[i]])[1]== "try-error" || lme4::isSingular(gr[[i]])) { # if the most-complex model is a try-error, is singular or DN converge
        gr[[i]] <- try(lme4::glmer(log_biom_2 ~ log_biom_1 + trt  + block + (1|startyear), data= data_i), silent= TRUE) # then fit a simpler model
        if (class(gr[[i]])[1] == "try-error" || lme4::isSingular(gr[[i]])) {  # if that model is singular or DN converge
          gr[[i]] <- try(glm(log_biom_1 ~ log_biom_2 + trt  + block + startyear, data= data_i), silent=TRUE)
          if (class(gr[[i]])[1] == "try-error") {print(paste("species", i, "gr DN converge- trt only"))}}}
      
      # prob of fruiting    
      pf[[i]] <-try(glmer(prf_2 ~ log_biom_1  + trt + (1|block)+ (1|startyear), data= data_i, family= "binomial"), silent=TRUE)
      if (class(pf[[i]])[1]== "try-error" || lme4::isSingular(pf[[i]])) { # if the most-complex model is a try-error, is singular or DN converge
        pf[[i]] <- try(lme4::glmer(prf_2 ~ log_biom_1 + trt  + block + (1|startyear), data= data_i, family= "binomial"), silent= TRUE) # then fit a simpler model
        if (class(pf[[i]])[1] == "try-error" || lme4::isSingular(pf[[i]])) {  # if that model is singular or DN converge
          pf[[i]] <- try(glm(prf_2 ~ log_biom_1 + trt  + block + startyear, data= data_i, family= "binomial"), silent=TRUE)
          if (class(pf[[i]])[1] == "try-error") {print(paste("species", i, "pf DN converge- trt only"))}}}
      
      # number of fruits
      cf[[i]] <- try(glmer(log_cf_2 ~ log_biom_1 + trt + (1|block)+ (1|startyear), data= data_i)  , silent=TRUE) 
      if (class(cf[[i]])[1]== "try-error" || lme4::isSingular(cf[[i]])) { # if the most-complex model is a try-error, is singular or DN converge
        cf[[i]] <- try(lme4::glmer(log_cf_2 ~ log_biom_1 + trt  + block + (1|startyear), data= data_i), silent= TRUE) # then fit a simpler model
        if (class(cf[[i]])[1] == "try-error" || lme4::isSingular(cf[[i]])) {  # if that model is singular or DN converge
          cf[[i]] <- try(glm(log_cf_2 ~ log_biom_1 + trt  + block + startyear, data= data_i), silent=TRUE)
          if (class(cf[[i]])[1] == "try-error") {print(paste("species", i, "cf DN converge- trt only"))}}}    }
print(i)}


# species 2 pf & cf DN converge; fix by hand-----
{data_i <- ecan_data
data_i <- left_join(data_i, trt)
data_i$watershed <- factor(data_i$watershed) %>% droplevels()
data_i$bison <- factor(data_i$bison) %>% droplevels()
data_i$small_mammal <- factor(data_i$small_mammal) %>% droplevels()
data_i$invertebrates <- factor( data_i$invertebrates) %>% droplevels()
data_i$trt<- factor(data_i$trt) %>% droplevels()
data_i$block <- factor(data_i$block) %>% droplevels()

data_i$log_biom_1 <- log(data_i$biom_1)
data_i$log_biom_2 <- log(data_i$biom_2)
data_i$log_biom_1 <- log(data_i$biom_1)
data_i$log_biom_2 <- log(data_i$biom_2)

data_i$prf_2 <- NA
data_i$prf_2[which(data_i$f_2 >0)] <- 1
data_i$prf_2[which(data_i$f_2 ==0)] <- 0
data_i$ log_cf_2 <- NA
data_i$ log_cf_2[which(data_i$f_2 >0)] <- log(data_i$f_2[which(data_i$f_2 >0)])}
contrasts_present <- data_i %>% 
  group_by(watershed) %>%
  summarise(n=sum(!is.na(sur_1_2)))


pf[[i]] <- glm(prf_2 ~ log_biom_1 +  trt  + watershed + startyear, data= data_i, family= "binomial") # block will not converge due to ecan's patchy distribution & limited occurrence in N1A
cf[[i]] <- glm(log_cf_2 ~ log_biom_1 + trt  + block + startyear, data= data_i) # ditto here-- watershed will not converge


mods <- list("A. canescens survival"= sur[[1]] , 
             "A. canescens growth"= gr[[1]] , 
             "A. canescens prob. rep."= pf[[1]] , 
             "A. canescens amt. rep."= cf[[1]] , 
             "E. angustifolia survival"= sur[[2]] , 
             "E. angustifolia growth"= gr[[2]] , 
             "E. angustifolia prob. rep."= pf[[2]] , 
             "E. angustifolia amt. rep."= cf[[2]] , 
             "B. eupatorioides survival"= sur[[2]] , 
             "B. eupatorioides growth"= gr[[2]] , 
             "B. eupatorioides prob. rep."= pf[[2]] , 
             "B. eupatorioides amt. rep."= cf[[2]] )

texreg::wordreg(mods, digits=3, stars=numeric(0) ,custom.coef.map = list(
  "(Intercept)" = "Intercept", 
  "GrowthFormW" = "Woody growth form"),
  custom.gof.rows = list("No. obs."= lapply(mods, nobs)
  ), # add row for AICc weight & no obs,
  file = "data/06_all_demo_rates.docx")
