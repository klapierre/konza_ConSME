################################################################################
##  05demography.R: tests for herbivore effects on 3 forbs
##
##  Author: Allison Louthan
##  Date created: January 16, 2025
##  Modified : Feb 26, 2025
# this code is adapted from: /Users/allisonlouthan/Dropbox/konza/watershed demography/paper 1/codes/paper codes/01_demodatacleanup.R
################################################################################

rm(list = ls())
# functions and libraries-----
{
  library(readxl)   
  library(stringr)
 
  re_assign <- function(df, drop) {
    df <- df [, ! names(df) %in% drop, drop = FALSE]
    df
  }
  
  Numextract <- function(string){ # this function propogates any NA's it finds
    answer <- as.numeric(unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string))))
    if(is.na(string)){answer <- NA}
    if(grepl("NA", string)) {answer <- NA}
    if(grepl("<", string)) {answer <- 0.25} # if diameter is <0.5, replace with 0.25
    
    return(answer)
  }
}

# load in data----

amca_data_20 <- readxl::read_xlsx('/Users/allisonlouthan/Dropbox/konza/watershed demography/annual data/2020/consolidated all field data files/ConSME-amca-2020-all.xlsx')
ecan_data_20 <- readxl::read_xlsx('/Users/allisonlouthan/Dropbox/konza/watershed demography/annual data/2020/consolidated all field data files/ConSME-ecan-2020-all.xlsx')
ecan_data_20 <- ecan_data_20[,which(names(ecan_data_20) != "meas")]
kueu_data_20 <- readxl::read_xlsx('/Users/allisonlouthan/Dropbox/konza/watershed demography/annual data/2020/consolidated all field data files/ConSME-kueu-2020-all.xlsx')
kueu_data_20 <- kueu_data_20[,which(names(kueu_data_20) != "meas")]

amca_data_21 <- readxl::read_xlsx('/Users/allisonlouthan/Dropbox/konza/watershed demography/annual data/2021/consolidated all field data files/ConSME-amca-21-all.xlsx')
ecan_data_21 <- readxl::read_xlsx('/Users/allisonlouthan/Dropbox/konza/watershed demography/annual data/2021/consolidated all field data files/ConSME-ecan-21-all.xlsx')
kueu_data_21 <- readxl::read_xlsx('/Users/allisonlouthan/Dropbox/konza/watershed demography/annual data/2021/consolidated all field data files/ConSME-kueu-21-all.xlsx')

amca_data_22 <- readxl::read_xlsx('/Users/allisonlouthan/Dropbox/konza/watershed demography/annual data/2022/consolidated all field data files/ConSME-amca-22-all.xlsx')
ecan_data_22 <- readxl::read_xlsx('/Users/allisonlouthan/Dropbox/konza/watershed demography/annual data/2022/consolidated all field data files/ConSME-ecan-22-all.xlsx')
kueu_data_22 <- readxl::read_xlsx('/Users/allisonlouthan/Dropbox/konza/watershed demography/annual data/2022/consolidated all field data files/ConSME-kueu-22-all.xlsx')

amca_data_23 <- readxl::read_xlsx('/Users/allisonlouthan/Dropbox/konza/watershed demography/annual data/2023/consolidated all field data files/ConSME-amca-23-ds.xlsx')
ecan_data_23 <- readxl::read_xlsx('/Users/allisonlouthan/Dropbox/konza/watershed demography/annual data/2023/consolidated all field data files/ConSME-ecan-23-ds.xlsx')
kueu_data_23 <- readxl::read_xlsx('/Users/allisonlouthan/Dropbox/konza/watershed demography/annual data/2023/consolidated all field data files/ConSME-kueu-23-ds.xlsx')

# remove plants that notes specify should be removed----

amca_data_20$note[grepl("disc", amca_data_20$note)]
amca_data_21$note[grepl("disc", amca_data_21$note)]
amca_data_22$note[grepl("disc", amca_data_22$note)]
amca_data_23$note[grepl("disc", amca_data_23$note)]
amca_data_20$comm[grepl("disc", amca_data_20$comm)]
amca_data_21$comm[grepl("disc", amca_data_21$comm)] # this comm is noted in recruitment file!
amca_data_22$comm[grepl("disc", amca_data_22$comm)]
amca_data_23$comm[grepl("disc", amca_data_23$comm)]

kueu_data_20$note[grepl("disc", kueu_data_20$note)]
kueu_data_21$note[grepl("disc", kueu_data_21$note)]
kueu_data_22$note[grepl("disc", kueu_data_22$note)] # 
kueu_data_23$note[grepl("disc", kueu_data_23$note)] # 
kueu_data_20$comm[grepl("disc", kueu_data_20$comm)]
kueu_data_21$comm[grepl("disc", kueu_data_21$comm)]
kueu_data_22$comm[grepl("disc", kueu_data_22$comm)] # 
kueu_data_23$comm[grepl("disc", kueu_data_23$comm)]

ecan_data_20$note[grepl("disc", ecan_data_20$note)]
ecan_data_21$note[grepl("disc", ecan_data_21$note)]
ecan_data_22$note[grepl("disc", ecan_data_22$note)]
ecan_data_23$note[grepl("disc", ecan_data_23$note)]
ecan_data_20$comm[grepl("disc", ecan_data_20$comm)]
ecan_data_21$comm[grepl("disc", ecan_data_21$comm)]
ecan_data_22$comm[grepl("disc", ecan_data_22$comm)]
ecan_data_23$comm[grepl("disc", ecan_data_23$comm)]



# combining the three years' data into one file ----
sp_name <- c("amca", "kueu", "ecan")
sp_measurements <- list (c(), 
                         c("n", "h", "f"), 
                         c("n_ros","l_ros", "f", "n_f","l_f" ))
for (i in 1:3){
  i_data_20 <- get(paste(sp_name[i]), "_data_20", sep= "")
  i_data_21 <- get(paste(sp_name[i]), "_data_21", sep= "")
  i_data_22 <- get(paste(sp_name[i]), "_data_22", sep= "")
  i_data_23 <- get(paste(sp_name[i]), "_data_23", sep= "")
  
  # make a fully unique plant id
  i_data_20$plant_id <- paste(i_data_20$watershed, i_data_20$plant_id)
  i_data_21$plant_id <- paste(i_data_21$watershed, i_data_21$plant_id)
  i_data_22$plant_id <- paste(i_data_22$watershed, i_data_22$plant_id)
  i_data_23$plant_id <- paste(i_data_23$watershed, i_data_23$plant_id)
  

# testing to see if there are plants that came back to life & assign survival

unique_plant_ids <- unique(c(i_data_20$plant_id, i_data_21$plant_id, i_data_22$plant_id,i_data_23$plant_id ))

for ( ii in 1: length(unique_plant_ids)){ # cycles through all the individual plants of that particular species
  alive_sequence <- c ( # what is the sequence of alive values for that individual plant
    ifelse( length(i_data_20$alive_startyear[which(i_data_20$plant_id== unique_plant_ids[ii] )])>0, i_data_20$alive_startyear[which(i_data_20$plant_id== unique_plant_ids[ii] )], NA), 
    ifelse( length(i_data_21$alive_startyear[which(i_data_21$plant_id== unique_plant_ids[ii] )])>0, i_data_21$alive_startyear[which(i_data_21$plant_id== unique_plant_ids[ii] )], NA), 
    ifelse( length(i_data_22$alive_startyear[which(i_data_22$plant_id== unique_plant_ids[ii] )])>0, i_data_22$alive_startyear[which(i_data_22$plant_id== unique_plant_ids[ii] )], NA),
    ifelse( length(i_data_23$alive_startyear[which(i_data_23$plant_id== unique_plant_ids[ii] )])>0, i_data_23$alive_startyear[which(i_data_23$plant_id== unique_plant_ids[ii] )], NA)
  )
  if (sum(as.numeric(alive_sequence), na.rm= TRUE) != 0){ # if a plant was only measured once, as dead, this will stop the code below and we need not do any replacing
    zeros <- which(alive_sequence== 0 & !is.na(alive_sequence))
    NAs <- which(is.na(alive_sequence))
    ones <- which(alive_sequence== 1 & !is.na(alive_sequence))
    # two rules-- these will work even if you have many years of data:
    
    # 1. for any NA's or 0's with indices greater than the idex of the first 1, but less than the index of the last 1, change the NA or zero to 1. 
    replace_with_ones <- NAs[data.table::between(NAs, min(ones), max(ones))]
    replace_with_ones <- c(replace_with_ones, zeros[data.table::between(zeros, min(ones), max(ones))])
    
    # 2. for the first zero index, if that is less than any 1 index, change all intervening indices to 1 (including the zero) 
    replace_with_ones <- c(replace_with_ones, zeros[which(zeros< min(ones))])
    
    # then, replace the survival values in the dataset, adding rows if necessary...
    if (length(replace_with_ones)>0){
      if((1 %in% replace_with_ones) & # if you are supposed to replace the first entry (i.e. year 2020)
         (length(i_data_20$alive_startyear[which(i_data_20$plant_id== unique_plant_ids[i] )]) >0) # and the plantID was present in that year
      ){i_data_20$alive_startyear[which(i_data_20$plant_id== unique_plant_ids[i] )] <- 1} # then replace it
      if((1 %in% replace_with_ones) & # if you are supposed to replace the first entry (i.e. year 2020)
         (length(i_data_20$alive_startyear[which(i_data_20$plant_id== unique_plant_ids[i] )]) ==0)){# but the plantID was not present in that year
        i_data_20[nrow(i_data_20) +1,] <- NA # then add a row that contains an "alive" = 1 entry for that plantID x year
        i_data_20[nrow(i_data_20) +1,"plant_id"] <- unique_plant_ids[i]
        i_data_20[nrow(i_data_20) +1,"alive_startyear"] <- 1 
        i_data_20[nrow(i_data_20) +1,"startyear"] <- 2020
      }
      if((2 %in% replace_with_ones) & (length(i_data_21$alive_startyear[which(i_data_21$plant_id== unique_plant_ids[i] )]) >0)){i_data_21$alive_startyear[which(i_data_21$plant_id== unique_plant_ids[i] )] <- 1}
      if((2 %in% replace_with_ones) & (length(i_data_21$alive_startyear[which(i_data_21$plant_id== unique_plant_ids[i] )]) ==0)){
        i_data_21[nrow(i_data_21) +1,] <- NA
        i_data_21[nrow(i_data_21) +1,"plant_id"] <- unique_plant_ids[i]
        i_data_21[nrow(i_data_21) +1,"alive_startyear"] <- 1 
        i_data_21[nrow(i_data_21) +1,"startyear"] <- 2021
      }
      if((3 %in% replace_with_ones) & (length(i_data_22$alive_startyear[which(i_data_22$plant_id== unique_plant_ids[i] )]) >0)){i_data_22$alive_startyear[which(i_data_22$plant_id== unique_plant_ids[i] )] <- 1}
      if((3 %in% replace_with_ones) & (length(i_data_22$alive_startyear[which(i_data_22$plant_id== unique_plant_ids[i] )]) ==0)){
        i_data_22[nrow(i_data_22) +1,] <- NA
        i_data_22[nrow(i_data_22) +1,"plant_id"] <- unique_plant_ids[i]
        i_data_22[nrow(i_data_22) +1,"alive_startyear"] <- 1 
        i_data_22[nrow(i_data_22) +1,"startyear"] <- 2022
      }
      if((4 %in% replace_with_ones) & (length(i_data_23$alive_startyear[which(i_data_23$plant_id== unique_plant_ids[i] )]) >0)){i_data_23$alive_startyear[which(i_data_23$plant_id== unique_plant_ids[i] )] <- 1}
      if((4 %in% replace_with_ones) & (length(i_data_23$alive_startyear[which(i_data_23$plant_id== unique_plant_ids[i] )]) ==0)){
        i_data_23[nrow(i_data_23) +1,] <- NA
        i_data_23[nrow(i_data_23) +1,"plant_id"] <- unique_plant_ids[i]
        i_data_23[nrow(i_data_23) +1,"alive_startyear"] <- 1 
        i_data_23[nrow(i_data_23) +1,"startyear"] <- 2023
      }
    }}
}
# need to make the final year of data, which is currently unfinished. including this last year becomes important when calculating recruitment later
my_colnames <- c("plant_id", "watershed",sp_measurements[[i]], "note", "cov", "alive_startyear")
endingyear <- data.frame(matrix(nrow=0, ncol=length(my_colnames)))
colnames(endingyear) <- my_colnames
endingyear$plant_id <- as.character(endingyear$plant_id)
endingyear$watershed <- as.character(endingyear$watershed)

i_data_transition1 <- dplyr::full_join(i_data_20,
                                          i_data_21[,c("plant_id", "watershed",sp_measurements[[i]], "note", "cov", "alive_startyear")], by= c("plant_id", "watershed")) 
i_data_transition2 <- dplyr::full_join(i_data_21, 
                                          i_data_22[,c("plant_id", "watershed" ,sp_measurements[[i]], "note", "cov", "alive_startyear")], by= c("plant_id", "watershed")) 
i_data_transition3 <- dplyr::full_join(i_data_22,
                                          i_data_23[,c("plant_id", "watershed" ,sp_measurements[[i]], "note", "cov", "alive_startyear")], by= c("plant_id", "watershed")) 
i_data_transition4 <- dplyr::full_join(i_data_23,
                                          endingyear, by= c("plant_id", "watershed")) 

assign(paste(sp_name[i], "_data", sep= ""), # name of dataset; in the form of kueu_data
       rbind(i_data_transition1[which(!is.na(i_data_transition1$startyear)),order(names(i_data_transition1))], 
                   i_data_transition2[which(!is.na(i_data_transition2$startyear)),order(names(i_data_transition2))], 
                   i_data_transition3[which(!is.na(i_data_transition3$startyear)),order(names(i_data_transition3))],
                   i_data_transition4[which(!is.na(i_data_transition4$startyear)),order(names(i_data_transition4))]))
} # end loop that combines data acrtoss years


# calculating biomass for each species---- 
kueu_data$h.y[which(kueu_data$h.y == "40.8, 2.9, 0.3, 0.5, 0.4,")] <- kueu_data$h.x[which(kueu_data$h.x == "40.8, 2.9, 0.3, 0.5, 0.4,")] <-  "40.8" # fixing some typos
kueu_data$h.y[which(kueu_data$h.y == "39.6,")] <- kueu_data$h.x[which(kueu_data$h.x == "39.6,")] <- "39.6"
kueu_data$h.y[which(kueu_data$h.y == "48.1,")] <- kueu_data$h.x[which(kueu_data$h.x == "48.1,")] <- "48.1"
kueu_data$h.x <- as.numeric(kueu_data$h.x)
kueu_data$h.y <- as.numeric(kueu_data$h.y)

kueu_data$biom_1 <- kueu_data$h.x *kueu_data$n.x
kueu_data$biom_2 <- kueu_data$h.y *kueu_data$n.y
kueu_data$f_1 <-  NA
kueu_data$f_2 <-  NA
for (i in 1:dim(kueu_data)[1]){
  kueu_data$f_1[i] <-  sum(Numextract(kueu_data$f.x[i]))
  kueu_data$f_2[i] <-  sum(Numextract(kueu_data$f.y[i]))
}
kueu_data$new_startyear <- NA
kueu_data$new_startyear[which(grepl("new", kueu_data$comm))] <- "new" # denotes those plants that have a "new" in the comments for startyear
kueu_data <- kueu_data[, c( "watershed",  "transect","x","y", "plant_id" , "startyear" ,  "biom_1"   ,
                            "biom_2" ,   "f_1"   ,    "f_2"   ,"cov.x" ,    "cov.y" , "alive_startyear.x", "alive_startyear.y" , "new_startyear")]
names(kueu_data) <- c( "watershed", "transect","x","y",  "plant_id" , "startyear" ,  "biom_1"   ,
                       "biom_2" ,   "f_1"   ,    "f_2"   ,"cov_1" ,    "cov_2", "alive_1", "alive_2", "new_startyear" )
kueu_data$bison <- substr(kueu_data$watershed, 1, 1)
kueu_data$fri <- substr(kueu_data$watershed, 2, 2)




ecan_data$biom_1 <- ecan_data$biom_2 <- ecan_data$f_1 <- ecan_data$f_2 <- NA

# fixing some errors illuminated by the loop below:
ecan_data$l_f.y[which(ecan_data$plant_id== "K1A 197"  & ecan_data$startyear== 2020)] <- "17.5, 12.6"
ecan_data[which(ecan_data$plant_id== "K2A 118" & ecan_data$startyear== 2020), c("n_ros.x", "l_ros.x")] <- NA
ecan_data[which(ecan_data$plant_id== "K2A 119" & ecan_data$startyear== 2020), c("n_f.x", "l_f.x")] <- NA
ecan_data[which(ecan_data$plant_id== "K2A 134" & ecan_data$startyear== 2020), c("n_ros.y", "l_ros.y")] <- NA
ecan_data$l_f.x[which(ecan_data$plant_id== "K1A 197" & ecan_data$startyear== 2021)] <- "17.5, 12.6"
ecan_data[which(ecan_data$plant_id== "K2A 109" & ecan_data$startyear== 2021), c("n_ros.y", "l_ros.y")] <- NA
ecan_data[which(ecan_data$plant_id== "K2A 109" & ecan_data$startyear== 2022), c("n_ros.x", "l_ros.x")] <- NA
ecan_data[which(ecan_data$plant_id== "K2A 134" & ecan_data$startyear== 2021), c("l_ros.x", "n_ros.x")] <- NA
ecan_data[which(ecan_data$plant_id== "K1A 189" & ecan_data$startyear== 2022), c("l_ros.y", "n_ros.y")] <- ecan_data[which(ecan_data$plant_id== "K1A 189" & ecan_data$startyear== 2023), c("l_ros.x", "n_ros.x")] <- NA
ecan_data[which(ecan_data$plant_id== "K1A 189" & ecan_data$startyear== 2022), c("l_ros.y", "n_ros.y")] <- ecan_data[which(ecan_data$plant_id== "K1A 189" & ecan_data$startyear== 2023), c("l_ros.x", "n_ros.x")] <- NA
ecan_data[which(ecan_data$plant_id == "N1A 62" & ecan_data$startyear== 2022), c("l_ros.y", "n_ros.y")] <- ecan_data[which(ecan_data$plant_id == "N1A 62" & ecan_data$startyear== 2023), c("l_ros.x", "n_ros.x")] <- NA
ecan_data[which(ecan_data$plant_id == "K2A 109" & ecan_data$startyear== 2022), c("l_ros.y", "n_ros.y")] <- ecan_data[which(ecan_data$plant_id == "K2A 109" & ecan_data$startyear== 2023), c("l_ros.x", "n_ros.x")] <- NA
ecan_data$l_ros.y[which(ecan_data$plant_id == "K2A 195" & ecan_data$startyear== 2022)] <- ecan_data$l_ros.x[which(ecan_data$plant_id == "K2A 195" & ecan_data$startyear== 2023)] <- "23.5, 23.3"
ecan_data$l_ros.y[which(ecan_data$plant_id== "N1A 101" & ecan_data$startyear== 2022)] <- ecan_data$l_ros.x[which(ecan_data$plant_id== "N1A 101" & ecan_data$startyear== 2023)] <- "19.5, 20, 21"


for (i in 1:dim(ecan_data)[1]){
  biom_ros1 <- sum(Numextract(ecan_data$n_ros.x[i]) *Numextract(ecan_data$l_ros.x[i]))
  biom_f1 <- sum(Numextract(ecan_data$n_f.x[i]) *Numextract(ecan_data$l_f.x[i]))
  ecan_data$biom_1[i] <- ifelse(is.na(biom_f1) & is.na(biom_ros1), NA, 
                                sum(biom_ros1, biom_f1, na.rm= TRUE))
  biom_ros2 <- sum(Numextract(ecan_data$n_ros.y[i]) *Numextract(ecan_data$l_ros.y[i]))
  biom_f2 <- sum(Numextract(ecan_data$n_f.y[i]) *Numextract(ecan_data$l_f.y[i]))
  ecan_data$biom_2[i] <-  ifelse(is.na(biom_f2) & is.na(biom_ros2), NA, 
                                 sum(biom_ros2, biom_f2, na.rm= TRUE))
  ecan_data$f_1[i] <- sum(Numextract(ecan_data$f.x[i]))
  ecan_data$f_2[i] <- sum(Numextract(ecan_data$f.y[i]))
  rm(biom_ros1, biom_f1, biom_ros2, biom_f2)
}
ecan_data$new_startyear <- NA
ecan_data$new_startyear[which(grepl("new", ecan_data$comm))] <- "new" # denotes those plants that have a "new" in the comments for startyear
ecan_data <- ecan_data[, c( "watershed", "transect","x","y",   "plant_id" , "startyear" ,  "biom_1"   ,
                            "biom_2" ,   "f_1"   ,    "f_2"   ,"cov.x" ,    "cov.y" , "alive_startyear.x" , "alive_startyear.y", "new_startyear")]
names(ecan_data) <- c( "watershed", "transect","x","y",  "plant_id" , "startyear" ,  "biom_1"   ,
                       "biom_2" ,   "f_1"   ,    "f_2"   ,"cov_1" ,    "cov_2", "alive_1", "alive_2","new_startyear"    )
ecan_data$bison <- substr(ecan_data$watershed, 1, 1)
ecan_data$fri <- substr(ecan_data$watershed, 2, 2)


# the below code is not yet edited, just straight from the cleanup file-----
# combining the three years' data for amca-----
amca_data_20 <- amca_data_20[,which(names(amca_data_20) != "meas")]

# testing to see if there are plants that came back alive & assign survival
amca_data_20$plant_id <- paste(amca_data_20$watershed, amca_data_20$plant_id)
amca_data_21$plant_id <- paste(amca_data_21$watershed, amca_data_21$plant_id)
amca_data_22$plant_id <- paste(amca_data_22$watershed, amca_data_22$plant_id)
amca_data_23$plant_id <- paste(amca_data_23$watershed, amca_data_23$plant_id)

unique_plant_ids <- unique(c(amca_data_20$plant_id, amca_data_21$plant_id, amca_data_22$plant_id,amca_data_23$plant_id))

for ( i in 1: length(unique_plant_ids)){ # cycles through all the individual plants of that particular species
  alive_sequence <- c (
    ifelse( length(amca_data_20$alive_startyear[which(amca_data_20$plant_id== unique_plant_ids[i] )])>0, amca_data_20$alive_startyear[which(amca_data_20$plant_id== unique_plant_ids[i] )], NA), 
    ifelse( length(amca_data_21$alive_startyear[which(amca_data_21$plant_id== unique_plant_ids[i] )])>0, amca_data_21$alive_startyear[which(amca_data_21$plant_id== unique_plant_ids[i] )], NA), 
    ifelse( length(amca_data_22$alive_startyear[which(amca_data_22$plant_id== unique_plant_ids[i] )])>0, amca_data_22$alive_startyear[which(amca_data_22$plant_id== unique_plant_ids[i] )], NA),
    ifelse( length(amca_data_23$alive_startyear[which(amca_data_23$plant_id== unique_plant_ids[i] )])>0, amca_data_23$alive_startyear[which(amca_data_23$plant_id== unique_plant_ids[i] )], NA)
    
  )
  if (sum(as.numeric(alive_sequence), na.rm= TRUE) != 0){ # if a plant was only measured once, as dead, this wil lbreak the code below and we need not do any replacing
    zeros <- which(alive_sequence== 0 & !is.na(alive_sequence))
    NAs <- which(is.na(alive_sequence))
    ones <- which(alive_sequence== 1 & !is.na(alive_sequence))
    # two rules-- these will work even if you have many years of data:
    
    # for any NA's or 0's with indices greater than the idex of the first 1, but less than the index of the last 1, chagne to 1. 
    replace_with_ones <- NAs[data.table::between(NAs, min(ones), max(ones))]
    replace_with_ones <- c(replace_with_ones, zeros[data.table::between(zeros, min(ones), max(ones))])
    
    # for the first zero index, if that is less than any 1 index, change all intervening indices to 1 (including the zero) 
    replace_with_ones <- c(replace_with_ones, zeros[which(zeros< min(ones))])
    
    # then, making the actual chagnes in the dataset, adding rows if necessary...
    if (length(replace_with_ones)>0){
      if((1 %in% replace_with_ones) & (length(amca_data_20$alive_startyear[which(amca_data_20$plant_id== unique_plant_ids[i] )]) >0)){amca_data_20$alive_startyear[which(amca_data_20$plant_id== unique_plant_ids[i] )] <- 1}
      if((1 %in% replace_with_ones) & (length(amca_data_20$alive_startyear[which(amca_data_20$plant_id== unique_plant_ids[i] )]) ==0)){
        amca_data_20[nrow(amca_data_20) +1,] <- NA
        amca_data_20[nrow(amca_data_20) +1,"plant_id"] <- unique_plant_ids[i]
        amca_data_20[nrow(amca_data_20) +1,"alive_startyear"] <- 1 
        amca_data_20[nrow(amca_data_20) +1,"startyear"] <- 2020
      }
      if((2 %in% replace_with_ones) & (length(amca_data_21$alive_startyear[which(amca_data_21$plant_id== unique_plant_ids[i] )]) >0)){amca_data_21$alive_startyear[which(amca_data_21$plant_id== unique_plant_ids[i] )] <- 1}
      if((2 %in% replace_with_ones) & (length(amca_data_21$alive_startyear[which(amca_data_21$plant_id== unique_plant_ids[i] )]) ==0)){
        amca_data_21[nrow(amca_data_21) +1,] <- NA
        amca_data_21[nrow(amca_data_21) +1,"plant_id"] <- unique_plant_ids[i]
        amca_data_21[nrow(amca_data_21) +1,"alive_startyear"] <- 1 
        amca_data_21[nrow(amca_data_21) +1,"startyear"] <- 2021
      }
      if((3 %in% replace_with_ones) & (length(amca_data_22$alive_startyear[which(amca_data_22$plant_id== unique_plant_ids[i] )]) >0)){amca_data_22$alive_startyear[which(amca_data_22$plant_id== unique_plant_ids[i] )] <- 1}
      if((3 %in% replace_with_ones) & (length(amca_data_22$alive_startyear[which(amca_data_22$plant_id== unique_plant_ids[i] )]) ==0)){
        amca_data_22[nrow(amca_data_22) +1,] <- NA
        amca_data_22[nrow(amca_data_22) +1,"plant_id"] <- unique_plant_ids[i]
        amca_data_22[nrow(amca_data_22) +1,"alive_startyear"] <- 1 
        amca_data_22[nrow(amca_data_22) +1,"startyear"] <- 2022
      }
      if((4 %in% replace_with_ones) & (length(amca_data_23$alive_startyear[which(amca_data_23$plant_id== unique_plant_ids[i] )]) >0)){amca_data_23$alive_startyear[which(amca_data_23$plant_id== unique_plant_ids[i] )] <- 1}
      if((4 %in% replace_with_ones) & (length(amca_data_23$alive_startyear[which(amca_data_23$plant_id== unique_plant_ids[i] )]) ==0)){
        amca_data_23[nrow(amca_data_23) +1,] <- NA
        amca_data_23[nrow(amca_data_23) +1,"plant_id"] <- unique_plant_ids[i]
        amca_data_23[nrow(amca_data_23) +1,"alive_startyear"] <- 1 
        amca_data_23[nrow(amca_data_23) +1,"startyear"] <- 2023
      }
    }}
}
# need to make the final year of data, which is currently unfinished. including this last year becomes important when calculating recruitment later
endingyear <- data.frame(matrix(nrow=0, ncol=9))
colnames(endingyear) <- c("plant_id", "watershed", "H","N","D", "L_F", "note", "cov", "alive_startyear")
endingyear$plant_id <- as.character(endingyear$plant_id)
endingyear$watershed <- as.character(endingyear$watershed)

amca_data_transition1 <- dplyr::full_join(amca_data_20,
                                          amca_data_21[,c("plant_id", "watershed","H", "N","D", "L_F", "note", "cov", "alive_startyear")], by= c("plant_id", "watershed")) 
amca_data_transition1$'H.x' <- NA
names(amca_data_transition1)[which(names(amca_data_transition1) == "H")] <- "H.y"
amca_data_21 <- amca_data_21[!(amca_data_21$plant_id== "N1A 284" & amca_data_21$mrk== "t0017"), ]
amca_data_22 <- amca_data_22[!(amca_data_22$plant_id== "N1A 284" & amca_data_22$mrk== "t0017"), ]
amca_data_transition2 <- dplyr::full_join(amca_data_21, 
                                          amca_data_22[,c("plant_id", "watershed", "H","N","D", "L_F", "note", "cov", "alive_startyear")], by= c("plant_id", "watershed")) 
amca_data_transition3 <- dplyr::full_join(amca_data_22, 
                                          amca_data_23[,c("plant_id", "watershed", "H","N","D", "L_F", "note", "cov", "alive_startyear")], by= c("plant_id", "watershed"))
amca_data_transition4 <- dplyr::full_join(amca_data_23, 
                                          endingyear, by= c("plant_id", "watershed"))
amca_data <- rbind(amca_data_transition1[which(!is.na(amca_data_transition1$startyear)),order(names(amca_data_transition1))], # this is.na term takes out those plants that were newly marked in the second year
                   amca_data_transition2[which(!is.na(amca_data_transition2$startyear)),order(names(amca_data_transition2))], 
                   amca_data_transition3[which(!is.na(amca_data_transition3$startyear)),order(names(amca_data_transition3))], 
                   amca_data_transition4[which(!is.na(amca_data_transition4$startyear)),order(names(amca_data_transition4))])

# best measurement of biomass is h*s (r^2= 0.65), but you didn't measure h the first year
# the first year you only measured d*s (r^2= 0.61)
# not 100% clear what to do here. I'll just use d*s

amca_data$biom_2 <- amca_data$biom_1 <- amca_data$f_1 <- amca_data$f_2 <-NA
# fixing some errors found in the code:
amca_data[which(amca_data$plant_id=="K2A 98" & amca_data$startyear== 2020 ), "D.y"] <- 7.75
amca_data[which(amca_data$plant_id=="K2A 98" & amca_data$startyear== 2021 ), "D.x"] <- 7.75

for (i in 1:dim(amca_data)[1]){
  amca_data$biom_1[i] <- sum(Numextract(amca_data$D.x[i])*Numextract(amca_data$N.x[i]))
  amca_data$biom_2[i] <- sum(Numextract(amca_data$D.y[i])*Numextract(amca_data$N.y[i]))
  amca_data$f_1[i] <- sum(Numextract(amca_data$L_F.x[i]))
  amca_data$f_2[i] <- sum(Numextract(amca_data$L_F.y[i]))
}
amca_data$new_startyear <- NA
amca_data$new_startyear[which(grepl("new", amca_data$comm))] <- "new" # denotes those plants that have a "new" in the comments for startyear
amca_data <- amca_data[, c( "watershed", "transect","x","y",   "plant_id" , "startyear" ,  "biom_1"   ,
                            "biom_2" ,   "f_1"   ,    "f_2"  ,"cov.x" ,    "cov.y" , "alive_startyear.x" , "alive_startyear.y", "new_startyear")]
names(amca_data) <- c( "watershed", "transect","x","y",   "plant_id" , "startyear" ,  "biom_1"   ,
                       "biom_2" ,   "f_1"   ,    "f_2"   ,"cov_1" ,    "cov_2" ,"alive_1", "alive_2", "new_startyear" )
amca_data$bison <- substr(amca_data$watershed, 1, 1)
amca_data$fri <- substr(amca_data$watershed, 2, 2)

# doing some final touches here
amca_data$sur_1_2 <- NA
amca_data$sur_1_2[which(amca_data$alive_1 ==1 & amca_data$alive_2==1)] <- 1
amca_data$sur_1_2[which(amca_data$alive_1 ==1 & amca_data$alive_2==0)] <- 0
if(length(which(amca_data$alive_1 ==0 & amca_data$alive_2==1))>0){stop("amca alive_sequence has a 0, 1")}

kueu_data$sur_1_2 <- NA
kueu_data$sur_1_2[which(kueu_data$alive_1 ==1 & kueu_data$alive_2==1)] <- 1
kueu_data$sur_1_2[which(kueu_data$alive_1 ==1 & kueu_data$alive_2==0)] <- 0
if(length(which(kueu_data$alive_1 ==0 & kueu_data$alive_2==1))>0){stop("kueu alive_sequence has a 0, 1")}

ecan_data$sur_1_2 <- NA
ecan_data$sur_1_2[which(ecan_data$alive_1 ==1 & ecan_data$alive_2==1)] <- 1
ecan_data$sur_1_2[which(ecan_data$alive_1 ==1 & ecan_data$alive_2==0)] <- 0
if(length(which(ecan_data$alive_1 ==0 & ecan_data$alive_2==1))>0){stop("ecan alive_sequence has a 0, 1")}

asob_data$sur_1_2 <- NA
asob_data$sur_1_2[which(asob_data$alive_1 ==1 & asob_data$alive_2==1)] <- 1
asob_data$sur_1_2[which(asob_data$alive_1 ==1 & asob_data$alive_2==0)] <- 0
if(length(which(asob_data$alive_1 ==0 & asob_data$alive_2==1))>0){stop("asob alive_sequence has a 0, 1")}

amca_data$f_1[which(!is.na(amca_data$biom_1) & is.na(amca_data$f_1))] <- 0 # replacing flower # of any plants with mreasurements but no flower number with a zero
amca_data$f_2[which(!is.na(amca_data$biom_2) & is.na(amca_data$f_2))] <- 0 
kueu_data$f_1[which(!is.na(kueu_data$biom_1) & is.na(kueu_data$f_1))] <- 0 
kueu_data$f_2[which(!is.na(kueu_data$biom_2) & is.na(kueu_data$f_2))] <- 0 
ecan_data$f_1[which(!is.na(ecan_data$biom_1) & is.na(ecan_data$f_1))] <- 0 
ecan_data$f_2[which(!is.na(ecan_data$biom_2) & is.na(ecan_data$f_2))] <- 0 
asob_data$f_1[which(!is.na(asob_data$biom_1) & is.na(asob_data$f_1))] <- 0 
asob_data$f_2[which(!is.na(asob_data$biom_2) & is.na(asob_data$f_2))] <- 0 

ecan_data$cov_1 <- as.numeric(ecan_data$ cov_1)
ecan_data$cov_2 <- as.numeric(ecan_data$ cov_2)
ecan_data$fri <- factor(ecan_data$fri)
ecan_data$bison <- factor(ecan_data$bison)

kueu_data$cov_1 <- as.numeric(kueu_data$ cov_1)
kueu_data$cov_2 <- as.numeric(kueu_data$ cov_2)
kueu_data$fri <- factor(kueu_data$fri)
kueu_data$bison <- factor(kueu_data$bison)

amca_data$cov_1 <- as.numeric(amca_data$ cov_1)
amca_data$cov_2 <- as.numeric(amca_data$ cov_2)
amca_data$fri <- factor(amca_data$fri)
amca_data$bison <- factor(amca_data$bison)

asob_data$cov_1 <- as.numeric(asob_data$ cov_1)
asob_data$cov_2 <- as.numeric(asob_data$ cov_2)
asob_data$fri <- factor(asob_data$fri)
asob_data$bison <- factor(asob_data$bison)

save(amca_data, asob_data, ecan_data, kueu_data, file="/Users/allisonlouthan/Dropbox/konza/watershed demography/paper 1/code outputs/01_clean_demo_data.RData")




