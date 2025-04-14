##  06vital_rate_analyses.R: tests for herbivore effects on vital rates of 3 forbs
##
##  Author: Allison Louthan
##  Date created: April 14, 2025
##  Modified : 
################################################################################

rm(list = ls())

# functions and libraries----
library("lme4") #NB this is DIFFERENT than what Kim used 

# load in data----

amca_data <- read.csv("derived_data/05_amca_clean_demo_data.RData")
ecan_data <- read.csv("derived_data/05_ecan_clean_demo_data.RData")
kueu_data <- read.csv("derived_data/05_kueu_clean_demo_data.RData")

trt <- read.csv('data/conSME_treatments.csv')

amca_data <- left_join(amca_data, trt)
ecan_data <- left_join(ecan_data, trt)
kueu_data <- left_join(kueu_data, trt)

# fitting vital rates

for (i in 1:3){
  if (i == 1) data_i <- amca_data
  if (i == 2) data_i <- ecan_data
  if (i == 3) data_i <- kueu_data
  
# survival
    glmer(sur_1_2 ~ log_biom_1 + I + B + S + WS + (1|block), data= data_i, family= "binomial")
# growth (mean only, no variatnce)
    glmer(log_biom_2 ~ log_biom_1  + I + B + S + WS + (1|block), data= data_i)
# prob of fruiting    
    glmer(prf_2 ~ log_biom_1  + I + B + S + WS + (1|block), data= data_i, family= "binomial")
# number of fruits
    glmer(log(cf_2) ~ log_biom_1 + I + B + S + WS + (1|block), data= data_i)
}