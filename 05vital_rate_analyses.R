##  06vital_rate_analyses.R: tests for herbivore effects on vital rates of 3 forbs
##
##  Author: Allison Louthan
##  Date created: April 14, 2025
##  Modified : 
################################################################################


rm(list = ls())

# functions and libraries----
library("lme4") #NB this is DIFFERENT than what Kim used 
library("dplyr")
library("ggplot2")
library("lmerTest")
# load in data----

amca_data <- read.csv("derived_data/05_amca_clean_demo_data.RData")
ecan_data <- read.csv("derived_data/05_ecan_clean_demo_data.RData")
kueu_data <- read.csv("derived_data/05_kueu_clean_demo_data.RData")


amca_rec <- read.csv("derived_data/05_amca_clean_rec_data.RData")
ecan_rec <- read.csv("derived_data/05_ecan_clean_rec_data.RData")
kueu_rec <- read.csv("derived_data/05_kueu_clean_rec_data.RData")


trt <- read.csv('data/conSME_treatments.csv')

# fitting vital rates----
#warning("there were ZERO new recruits for any of the species in any of the yeras")
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
  data_i$trt<- factor(data_i$trt, levels= c("BSI","BSX", "XSI",  "XSX", "XXI", "XXX")) %>% droplevels()
  data_i$block <- factor(data_i$block) %>% droplevels()
  data_i$ startyear <- as.factor(data_i$startyear)
  if (length(which( data_i$trt == "BSI"))==0){stop("BSI reference treatment not present!")}
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
      glmer(sur_1_2 ~ log_biom_1 + watershed*trt  + (1|block)+ (1|startyear), data= data_i, family= "binomial")
                    , silent=TRUE)
    if (class(sur[[i]])[1]== "try-error" || lme4::isSingular(sur[[i]])) { # if the most-complex model is a try-error, is singular or DN converge
      sur[[i]] <- try(glmer(sur_1_2 ~ log_biom_1 + watershed*trt  + block + (1|startyear), data= data_i, family= "binomial"), silent= TRUE) # then fit a simpler model
      if (class(sur[[i]])[1] == "try-error" || lme4::isSingular(sur[[i]])) {  # if that model is singular or DN converge
        sur[[i]] <- try(glm(sur_1_2 ~ log_biom_1 + watershed*trt  + block + startyear, data= data_i, family= "binomial"), silent=TRUE)
        if (class(sur[[i]])[1] == "try-error") {print(paste("species", i, "sur DN converge"))}}}
    
# growth (mean only, no variatnce)
    gr[[i]] <-try(
      lmer(log_biom_2 ~ log_biom_1  + watershed*trt + (1|block)+ (1|startyear), data= data_i), 
      silent=TRUE)
    if (class(gr[[i]])[1]== "try-error" || lme4::isSingular(gr[[i]])) { # if the most-complex model is a try-error, is singular or DN converge
      gr[[i]] <- try(lmer(log_biom_2 ~ log_biom_1 + watershed*trt  + block + (1|startyear), data= data_i), silent= TRUE) # then fit a simpler model
      if (class(gr[[i]])[1] == "try-error" || lme4::isSingular(gr[[i]])) {  # if that model is singular or DN converge
        gr[[i]] <- try(glm(log_biom_2 ~ log_biom_1 + watershed*trt  + block + startyear, data= data_i), silent=TRUE)
        if (class(gr[[i]])[1] == "try-error") {print(paste("species", i, "gr DN converge"))}}}

  # prob of fruiting    
    pf[[i]] <-try(glmer(prf_2 ~ log_biom_1  + watershed*trt + (1|block)+ (1|startyear), data= data_i, family= "binomial"), silent=TRUE)
    if (class(pf[[i]])[1]== "try-error" || lme4::isSingular(pf[[i]])) { # if the most-complex model is a try-error, is singular or DN converge
      pf[[i]]  <- try(glmer(prf_2 ~ log_biom_1 + watershed*trt  + block + (1|startyear), data= data_i, family= "binomial"), silent= TRUE) # then fit a simpler model
      if (class(pf[[i]])[1] == "try-error" || lme4::isSingular(pf[[i]])) {  # if that model is singular or DN converge
        pf[[i]]  <- try(glm(prf_2 ~ log_biom_1 + watershed*trt  + block + startyear, data= data_i, family= "binomial"), silent=TRUE)
      
        if (class(pf[[i]])[1] == "try-error") {print(paste("species", i, "pf DN converge"))}}}
    
# number of fruits
    cf[[i]]<-try(lmer(log_cf_2 ~ log_biom_1 + watershed*trt + (1|block)+ (1|startyear), data= data_i)  , silent=TRUE) 
    if (class(cf[[i]])[1]== "try-error" || lme4::isSingular(cf[[i]])) { # if the most-complex model is a try-error, is singular or DN converge
      cf[[i]] <- try(lmer(log_cf_2 ~ log_biom_1 + watershed*trt  + block + (1|startyear), data= data_i), silent= TRUE) # then fit a simpler model
      if (class(cf[[i]])[1] == "try-error" || lme4::isSingular(cf[[i]])) {  # if that model is singular or DN converge
        cf[[i]] <- try(lm(log_cf_2 ~ log_biom_1 + watershed*trt  + block + startyear, data= data_i), silent=TRUE)
        if (class(cf[[i]])[1] == "try-error") {print(paste("species", i, "cf DN converge"))}}}
    
    
    } else { # if there is data in only one watershed, remove the watershed effect
      sur[[i]] <- try(
        glmer(sur_1_2 ~ log_biom_1 + trt  + (1|block)+ (1|startyear), data= data_i, family= "binomial")
        , silent=TRUE)
      if (class(sur[[i]])[1]== "try-error" || lme4::isSingular(sur[[i]])) { # if the most-complex model is a try-error, is singular or DN converge
        sur[[i]]  <- try(glmer(sur_1_2 ~ log_biom_1 + trt  + block + (1|startyear), data= data_i, family= "binomial"), silent= TRUE) # then fit a simpler model
        if (class(sur[[i]])[1] == "try-error" || lme4::isSingular(sur[[i]])) {  # if that model is singular or DN converge
          sur[[i]]  <- try(glm(sur_1_2 ~ log_biom_1 + trt  + block + startyear, data= data_i, family= "binomial"), silent=TRUE)
          if (class(sur[[i]])[1] == "try-error") {print(paste("species", i, "sur DN converge- trt only"))}}}
      
      # growth (mean only, no variatnce)
      gr[[i]] <-try(
        lmer(log_biom_2 ~ log_biom_1  + trt + (1|block)+ (1|startyear), data= data_i), 
        silent=TRUE)
      if (class(gr[[i]])[1]== "try-error" || lme4::isSingular(gr[[i]])) { # if the most-complex model is a try-error, is singular or DN converge
        gr[[i]] <- try(lmer(log_biom_2 ~ log_biom_1 + trt  + block + (1|startyear), data= data_i), silent= TRUE) # then fit a simpler model
        if (class(gr[[i]])[1] == "try-error" || lme4::isSingular(gr[[i]])) {  # if that model is singular or DN converge
          gr[[i]] <- try(glm(log_biom_2 ~ log_biom_1 + trt  + block + startyear, data= data_i), silent=TRUE)
          if (class(gr[[i]])[1] == "try-error") {print(paste("species", i, "gr DN converge- trt only"))}}}
      
      # prob of fruiting    
      pf[[i]] <-try(glmer(prf_2 ~ log_biom_1  + trt + (1|block)+ (1|startyear), data= data_i, family= "binomial"), silent=TRUE)
      if (class(pf[[i]])[1]== "try-error" || lme4::isSingular(pf[[i]])) { # if the most-complex model is a try-error, is singular or DN converge
        pf[[i]] <- try(glmer(prf_2 ~ log_biom_1 + trt  + block + (1|startyear), data= data_i, family= "binomial"), silent= TRUE) # then fit a simpler model
        if (class(pf[[i]])[1] == "try-error" || lme4::isSingular(pf[[i]])) {  # if that model is singular or DN converge
          pf[[i]] <- try(glm(prf_2 ~ log_biom_1 + trt  + block + startyear, data= data_i, family= "binomial"), silent=TRUE)
          if (class(pf[[i]])[1] == "try-error") {print(paste("species", i, "pf DN converge- trt only"))}}}
      
      # number of fruits
      cf[[i]] <- try(lmer(log_cf_2 ~ log_biom_1 + trt + (1|block)+ (1|startyear), data= data_i)  , silent=TRUE) 
      if (class(cf[[i]])[1]== "try-error" || lme4::isSingular(cf[[i]])) { # if the most-complex model is a try-error, is singular or DN converge
        cf[[i]] <- try(lmer(log_cf_2 ~ log_biom_1 + trt  + block + (1|startyear), data= data_i), silent= TRUE) # then fit a simpler model
        if (class(cf[[i]])[1] == "try-error" || lme4::isSingular(cf[[i]])) {  # if that model is singular or DN converge
          cf[[i]] <- try(glm(log_cf_2 ~ log_biom_1 + trt  + block + startyear, data= data_i), silent=TRUE)
          if (class(cf[[i]])[1] == "try-error") {print(paste("species", i, "cf DN converge- trt only"))}}}    }
  
  sur[[i]]$xlevels
print(i)}


# species 2 pf & cf DN converge; fix by hand-----
{data_i <- ecan_data
data_i <- left_join(data_i, trt)
data_i$watershed <- factor(data_i$watershed) %>% droplevels()
data_i$bison <- factor(data_i$bison) %>% droplevels()
data_i$small_mammal <- factor(data_i$small_mammal) %>% droplevels()
data_i$invertebrates <- factor( data_i$invertebrates) %>% droplevels()
data_i$trt <- factor(data_i$trt, levels= c("BSI","BSX", "XSI",  "XSX", "XXI", "XXX")) %>% droplevels()
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


pf[[2]] <- glm(prf_2 ~ log_biom_1 +  trt  + watershed + startyear, data= data_i, family= "binomial") # block will not converge due to ecan's patchy distribution & limited occurrence in N1A
cf[[2]] <- glm(log_cf_2 ~ log_biom_1 + trt  + block + startyear, data= data_i) # ditto here-- watershed will not converge


mods <- list("Survival"= sur[[1]] , 
             "Growth"= gr[[1]] , 
             "Prob. rep."= pf[[1]] , 
             "Amt. rep."= cf[[1]] , 
             "Survival"= sur[[2]] ,
             "Growth"= gr[[2]] , 
             "Prob. rep."= pf[[2]] , 
             "Amt. rep."= cf[[2]] , 
             "Survival"= sur[[3]] , 
             "Growth"= gr[[3]] , 
             "Prob. rep."= pf[[3]] , 
             "Amt. rep."= cf[[3]] )

texreg::wordreg(mods, digits=2, stars=0.05 , 
 custom.coef.map = list(
  "(Intercept)" = "Intercept", 
  "log_biom_1" = "Initial size",
  "trtXXX"= "Treatment (-B-S-I)", 
  "trtXXI"= "Treatment (-B-S+I)",
  "trtXSX"= "Treatment (-B+S-I)",
  "trtBSX"= "Treatment (+B+S-I)",
  "trtXSI"= "Treatment (-B+S+I)",
  "trtBSI"= "Treatment (+B+S+I)",
  "watershedN4B" = "4-year burn",
  "watershedN4B:trtXXS"= "4-year burn x treatment (-B-S-I)",
  "watershedN4B:trtXXI"= "4-year burn x treatment (-B-S+I)",
  "watershedN4B:trtXSX"= "4-year burn x treatment (-B+S-I)",
  "watershedN4B:trtBSX"= "4-year burn x treatment (+B+S-I)",
  "watershedN4B:trtXSI"= "4-year burn x treatment (-B+S+I)",
  "watershedN4B:trtBSI"= "4-year burn x treatment (+B+S+I)"
 ),
  custom.gof.rows = list("No. obs."= as.numeric(unlist(lapply(mods, nobs)))
  ), # add row for AICc weight & no obs,
  custom.header= list("A. canescens" = 1:4,"E. angustifolia" = 5:8,"B. eupatorioides" = 9:12 ), # this custom header isn't woring
  file = "derived_data/06_all_demo_rates.docx")

# now, make a figure-----
library(car)
mods <- list("Survival"= sur[[1]] , # amca
             "Growth"= gr[[1]] , 
             "Prob. rep."= pf[[1]] , 
             "Amt. rep."= cf[[1]] , 
             "Survival"= sur[[2]] , #ecan
             "Growth"= gr[[2]] , 
             "Prob. rep."= pf[[2]] , 
             "Amt. rep."= cf[[2]] , 
             "Survival"= sur[[3]] , #kueu
             "Growth"= gr[[3]] , 
             "Prob. rep."= pf[[3]] , 
             "Amt. rep."= cf[[3]] )
sp_ref <- c("amca", "amca", "amca", "amca", 
            "ecan", "ecan", "ecan", "ecan", 
            "kueu", "kueu", "kueu", "kueu" )
rate_ref <- rep(c("Survival", "Growth", "Prob. rep.", "Amt. rep."), 3)
all_coefs <- expand.grid(Species= c("amca", "ecan","kueu"), 
                         Rate= c("Survival", "Growth", "Prob. rep.", "Amt. rep."), 
                         Coefficient= c("BSX", "XXI", "XSI", "XXX","BSX_N4B", "XXI_N4B", "XSI_N4B"), Estimate= NA, 'Std. Error'= NA, 'pval' = NA)

for (i in 1: length(mods)){
  
my_mod <- mods[[i]]

BSX_coeffs_i <- XXI_coeffs_i <- XSI_coeffs_i <- XXX_coeffs_i <- BSX_coeffs_i_N4B <- XXI_coeffs_i_N4B <- XSI_coeffs_i_N4B <- c(NA, NA, NA)

if ("trtBSX" %in% rownames(summary(my_mod)$coefficients) ){
BSX_coeffs_i <- c(summary(my_mod)$coefficients["trtBSX", "Estimate"], 
                  summary(my_mod)$coefficients["trtBSX", "Std. Error"], 
                  summary(my_mod)$coefficients["trtBSX", which(grepl("Pr", dimnames(summary(my_mod)$coefficients)[2][[1]]))])}
if ("trtXXI" %in% rownames(summary(my_mod)$coefficients) ){
  XXI_coeffs_i <- c(summary(my_mod)$coefficients["trtXXI", "Estimate"], 
                    summary(my_mod)$coefficients["trtXXI", "Std. Error"], 
                    summary(my_mod)$coefficients["trtXXI", which(grepl("Pr", dimnames(summary(my_mod)$coefficients)[2][[1]]))])}
if ("trtXSI" %in% rownames(summary(my_mod)$coefficients) ){
  XSI_coeffs_i <- c(summary(my_mod)$coefficients["trtXSI", "Estimate"], 
                    summary(my_mod)$coefficients["trtXSI", "Std. Error"], 
                    summary(my_mod)$coefficients["trtXSI", which(grepl("Pr", dimnames(summary(my_mod)$coefficients)[2][[1]]))])}
if ("trtXXX" %in% rownames(summary(my_mod)$coefficients) ){
  XXX_coeffs_i <- c(summary(my_mod)$coefficients["trtXXX", "Estimate"], 
                    summary(my_mod)$coefficients["trtXXX", "Std. Error"], 
                    summary(my_mod)$coefficients["trtXXX", which(grepl("Pr", dimnames(summary(my_mod)$coefficients)[2][[1]]))])}

all_coefs[which(all_coefs$Species == sp_ref[i] & all_coefs$Rate == rate_ref[i] & all_coefs$Coefficient == "BSX"
                  ), c("Estimate", "Std. Error", "pval")] <- BSX_coeffs_i
all_coefs[which(all_coefs$Species == sp_ref[i] & all_coefs$Rate == rate_ref[i] & all_coefs$Coefficient == "XXI"
                  ), c("Estimate", "Std. Error", "pval")] <- XXI_coeffs_i
all_coefs[which(all_coefs$Species == sp_ref[i] & all_coefs$Rate == rate_ref[i] & all_coefs$Coefficient == "XSI"
                  ), c("Estimate", "Std. Error", "pval")] <- XSI_coeffs_i
all_coefs[which(all_coefs$Species == sp_ref[i] & all_coefs$Rate == rate_ref[i] & all_coefs$Coefficient == "XXX"
), c("Estimate", "Std. Error", "pval")] <- XXX_coeffs_i
}

library(ggplot2)
library(dplyr)

df_plot <- all_coefs %>%
  filter(!(Species == "amca")) %>%   # drop amca
  filter(!(Species == "ecan" & Rate== "Survival" & Coefficient== "XSI")) %>%   # drop amca
  filter(Coefficient %in% c("BSX", "XSI", "XXI", "XXX")) %>%
  mutate(
    Coefficient = factor(
      Coefficient,
      levels = c("BSX", "XSI", "XXI", "XXX"),
      labels = c("Insects",
                 "Bison",
                 "Bison +\nsm. mam.", 
                 "All")
    ),
    Rate = factor(Rate),
    sig = ifelse(pval <= 0.05, "*", "")   # mark significant results
  )

# base plot
p <- ggplot(df_plot, aes(x = Coefficient, y = Estimate, color = Rate, group = Rate)) +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(
    aes(ymin = Estimate - `Std. Error`, ymax = Estimate + `Std. Error`),
    width = 0.2,
    position = position_dodge(width = 0.6)
  ) +
  # add asterisks above error bars
  geom_text(
    aes(
      label = sig,
      y = Estimate + `Std. Error` + 0.05
    ),
    position = position_dodge(width = 0.6),
    vjust = 0,
    size = 8,
    show.legend = FALSE
  ) +
  facet_wrap(~ Species, nrow = 1, scales = "fixed") +
  coord_cartesian(ylim = c(-2, 5)) +
  theme_classic(base_size = 24) +
  theme(
    strip.text = element_blank(),
    legend.position = c(0.2, 0.85),        # inside panel A
    legend.background = element_rect(fill = alpha("white", 0.7))
  ) +
  labs(
    x = "Herbivore group removed",
    y = "Estimated effect",
    color = "Rate"
  )

# add manual A/B tags inside each panel
p + geom_text(
  data = data.frame(
    Species = c("ecan", "kueu"),
    x = 0.1,   # left-most x position
    y = 5,     # just above y-limit
    label = c("A", "B")
  ),
  aes(x = x, y = y, label = label),
  inherit.aes = FALSE,
  hjust = -0.5, vjust = 1,
  size = 10
)
