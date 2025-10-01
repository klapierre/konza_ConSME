################################################################################
##  05demography_data_cleanup.R: cleans up data to test for consume effects on three forbs
##
##  Author: Allison Louthan
##  Date created: January 16, 2025
##  Modified : April 14, 2025
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

amca_data_20 <- readxl::read_xlsx('data/demo data/ConSME-amca-2020-all.xlsx')
amca_data_20 <- amca_data_20[,which(names(amca_data_20) != "meas")]
amca_data_20$H <- NA # incorrect size measurement in 2020
ecan_data_20 <- readxl::read_xlsx('data/demo data/ConSME-ecan-2020-all.xlsx')
ecan_data_20 <- ecan_data_20[,which(names(ecan_data_20) != "meas")]
kueu_data_20 <- readxl::read_xlsx('data/demo data/ConSME-kueu-2020-all.xlsx')
kueu_data_20 <- kueu_data_20[,which(names(kueu_data_20) != "meas")]

amca_data_21 <- readxl::read_xlsx('data/demo data/ConSME-amca-21-all.xlsx')
ecan_data_21 <- readxl::read_xlsx('data/demo data/ConSME-ecan-21-all.xlsx')
kueu_data_21 <- readxl::read_xlsx('data/demo data/ConSME-kueu-21-all.xlsx')

amca_data_22 <- readxl::read_xlsx('data/demo data/ConSME-amca-22-all.xlsx')
ecan_data_22 <- readxl::read_xlsx('data/demo data/ConSME-ecan-22-all.xlsx')
kueu_data_22 <- readxl::read_xlsx('data/demo data/ConSME-kueu-22-all.xlsx')

amca_data_23 <- readxl::read_xlsx('data/demo data/ConSME-amca-23-ds.xlsx')
ecan_data_23 <- readxl::read_xlsx('data/demo data/ConSME-ecan-23-ds.xlsx')
kueu_data_23 <- readxl::read_xlsx('data/demo data/ConSME-kueu-23-ds.xlsx')

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
kueu_data_21 <- dplyr::rename(kueu_data_21, alive21= alive20)


# combining the three years' data into one file ----

sp_name <- c("amca", "kueu", "ecan")
sp_measurements <- list (c("h", "n", "d", "l_f"), 
                         c("n", "h", "f"), 
                         c("n_ros","l_ros", "f", "n_f","l_f" ))
for (i in 1:3){
  i_data_20 <- get(paste(sp_name[i], "_data_20", sep= ""))
  i_data_21 <- get(paste(sp_name[i], "_data_21", sep= ""))
  i_data_22 <- get(paste(sp_name[i], "_data_22", sep= ""))
  i_data_23 <- get(paste(sp_name[i], "_data_23", sep= ""))
  names(i_data_20) <- tolower(names(i_data_20))
  names(i_data_21) <- tolower(names(i_data_21))
  names(i_data_22) <- tolower(names(i_data_22))
  names(i_data_23) <- tolower(names(i_data_23))
  i_data_20$block  <- toupper(i_data_20$block)
  i_data_21$block  <- toupper(i_data_21$block)
  i_data_22$block  <- toupper(i_data_22$block)
  i_data_23$block  <- toupper(i_data_23$block)
  i_data_20 <- dplyr::rename(i_data_20, alive_startyear= alive20)
  i_data_21 <- dplyr::rename(i_data_21, alive_startyear= alive21)
  i_data_22 <- dplyr::rename(i_data_22, alive_startyear= alive22)
  i_data_23 <- dplyr::rename(i_data_23, alive_startyear= alive23)
  
  
  i_data_20$startyear  <- 2020
  i_data_21$startyear  <- 2021
  i_data_22$startyear  <- 2022
  i_data_23$startyear  <- 2023
  
  # make a fully unique plant id
  i_data_20$plant_id <- paste(i_data_20$block, i_data_20$plot, i_data_20$x, i_data_20$y, i_data_20$mrk)
  i_data_21$plant_id <- paste(i_data_21$block, i_data_21$plot, i_data_21$x, i_data_21$y, i_data_21$mrk)
  i_data_22$plant_id <- paste(i_data_22$block, i_data_22$plot, i_data_22$x, i_data_22$y, i_data_22$mrk)
  i_data_23$plant_id <- paste(i_data_23$block, i_data_23$plot, i_data_23$x, i_data_23$y, i_data_23$mrk)  
  
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
         (length(i_data_20$alive_startyear[which(i_data_20$plant_id== unique_plant_ids[ii] )]) >0) # and the plantID was present in that year
      ){i_data_20$alive_startyear[which(i_data_20$plant_id== unique_plant_ids[ii] )] <- 1} # then replace it
      if((1 %in% replace_with_ones) & # if you are supposed to replace the first entry (i.e. year 2020)
         (length(i_data_20$alive_startyear[which(i_data_20$plant_id== unique_plant_ids[ii] )]) ==0)){# but the plantID was not present in that year
        i_data_20[nrow(i_data_20) +1,] <- NA # then add a row that contains an "alive" = 1 entry for that plantID x year
        i_data_20[nrow(i_data_20) +1,"plant_id"] <- unique_plant_ids[ii]
        i_data_20[nrow(i_data_20) +1,"alive_startyear"] <- 1 
        i_data_20[nrow(i_data_20) +1,"startyear"] <- 2020
      }
      if((2 %in% replace_with_ones) & (length(i_data_21$alive_startyear[which(i_data_21$plant_id== unique_plant_ids[ii] )]) >0)){i_data_21$alive_startyear[which(i_data_21$plant_id== unique_plant_ids[ii] )] <- 1}
      if((2 %in% replace_with_ones) & (length(i_data_21$alive_startyear[which(i_data_21$plant_id== unique_plant_ids[ii] )]) ==0)){
        i_data_21[nrow(i_data_21) +1,] <- NA
        i_data_21[nrow(i_data_21) +1,"plant_id"] <- unique_plant_ids[ii]
        i_data_21[nrow(i_data_21) +1,"alive_startyear"] <- 1 
        i_data_21[nrow(i_data_21) +1,"startyear"] <- 2021
      }
      if((3 %in% replace_with_ones) & (length(i_data_22$alive_startyear[which(i_data_22$plant_id== unique_plant_ids[ii] )]) >0)){i_data_22$alive_startyear[which(i_data_22$plant_id== unique_plant_ids[ii] )] <- 1}
      if((3 %in% replace_with_ones) & (length(i_data_22$alive_startyear[which(i_data_22$plant_id== unique_plant_ids[ii] )]) ==0)){
        i_data_22[nrow(i_data_22) +1,] <- NA
        i_data_22[nrow(i_data_22) +1,"plant_id"] <- unique_plant_ids[ii]
        i_data_22[nrow(i_data_22) +1,"alive_startyear"] <- 1 
        i_data_22[nrow(i_data_22) +1,"startyear"] <- 2022
      }
      if((4 %in% replace_with_ones) & (length(i_data_23$alive_startyear[which(i_data_23$plant_id== unique_plant_ids[ii] )]) >0)){i_data_23$alive_startyear[which(i_data_23$plant_id== unique_plant_ids[ii] )] <- 1}
      if((4 %in% replace_with_ones) & (length(i_data_23$alive_startyear[which(i_data_23$plant_id== unique_plant_ids[ii] )]) ==0)){
        i_data_23[nrow(i_data_23) +1,] <- NA
        i_data_23[nrow(i_data_23) +1,"plant_id"] <- unique_plant_ids[ii]
        i_data_23[nrow(i_data_23) +1,"alive_startyear"] <- 1 
        i_data_23[nrow(i_data_23) +1,"startyear"] <- 2023
      }
    }}
}
# need to make the final year of data, which is currently unfinished. including this last year becomes important when calculating recruitment later
my_colnames <- c("plant_id",sp_measurements[[i]], "comm", "note", "cov", "alive_startyear")
endingyear <- data.frame(matrix(nrow=0, ncol=length(my_colnames)))
colnames(endingyear) <- my_colnames
endingyear$plant_id <- as.character(endingyear$plant_id)
if (i == 1) {i_data_20$cov <- NA}
i_data_transition1 <- dplyr::full_join(i_data_20,
                                          i_data_21[,c("plant_id", sp_measurements[[i]], "comm", "note", "cov", "alive_startyear")], by= c("plant_id")) 
i_data_transition2 <- dplyr::full_join(i_data_21, 
                                          i_data_22[,c("plant_id", sp_measurements[[i]], "comm", "note", "cov", "alive_startyear")], by= c("plant_id")) 
i_data_transition3 <- dplyr::full_join(i_data_22,
                                          i_data_23[,c("plant_id", sp_measurements[[i]], "comm", "note", "cov", "alive_startyear")], by= c("plant_id")) 
i_data_transition4 <- dplyr::full_join(i_data_23,
                                          endingyear, by= c("plant_id")) 

assign(paste(sp_name[i], "_data", sep= ""), # name of dataset; in the form of kueu_data
       rbind(i_data_transition1[which(!is.na(i_data_transition1$startyear)),sort(names(i_data_transition1))], 
                   i_data_transition2[which(!is.na(i_data_transition2$startyear)),sort(names(i_data_transition2))], 
                   i_data_transition3[which(!is.na(i_data_transition3$startyear)),sort(names(i_data_transition3))],
                   i_data_transition4[which(!is.na(i_data_transition4$startyear)),sort(names(i_data_transition4))]))
} # end loop that combines data acrtoss years


# calculating biomass for each species---- 
kueu_data$h.x[which(kueu_data$h.x=="17, 5" )] <- NA
kueu_data$h.y[which(kueu_data$h.y=="17, 5" )] <- NA

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
kueu_data <- kueu_data[, c( "block", "plot","x","y", "plant_id" , "startyear" ,  "biom_1"   ,
                            "biom_2" , "comm.y",   "f_1"   ,    "f_2"   ,"cov.x" ,    "cov.y" , "alive_startyear.x", "alive_startyear.y" )]
names(kueu_data) <- c( "block", "plot","x","y",  "plant_id" , "startyear" ,  "biom_1"   ,
                       "biom_2" ,  "comm_2",  "f_1"   ,    "f_2"   ,"cov_1" ,    "cov_2", "alive_1", "alive_2")

ecan_data$biom_1 <- ecan_data$biom_2 <- ecan_data$f_1 <- ecan_data$f_2 <- NA

# fixing some errors illuminated by the loop below:
ecan_data$l_ros.y[which(ecan_data$l_ros.y== "25.2, 25.7, 22,2")] <- "25.2, 25.7, 22.2"
ecan_data$l_ros.x[which(ecan_data$l_ros.x== "25.2, 25.7, 22,2")] <- "25.2, 25.7, 22.2"

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
ecan_data <- ecan_data[, c( "block", "plot","x","y",   "plant_id" , "startyear" ,  "biom_1"   ,
                            "biom_2" , "comm.y",   "f_1"   ,    "f_2"   ,"cov.x" ,    "cov.y" , "alive_startyear.x" , "alive_startyear.y")]
names(ecan_data) <- c( "block", "plot","x","y",  "plant_id" , "startyear" ,  "biom_1"   ,
                       "biom_2" , "comm_2",   "f_1"   ,    "f_2"   ,"cov_1" ,    "cov_2", "alive_1", "alive_2"  )

#amca
# best measurement of biomass is h*s (r^2= 0.65), but you didn't measure h the first year
# the first year you only measured d*s (r^2= 0.61)
# not 100% clear what to do here. I'll just use d*s

amca_data$biom_2 <- amca_data$biom_1 <- amca_data$f_1 <- amca_data$f_2 <-NA

for (i in 1:dim(amca_data)[1]){
  amca_data$biom_1[i] <- sum(Numextract(amca_data$d.x[i])*Numextract(amca_data$n.x[i]))
  amca_data$biom_2[i] <- sum(Numextract(amca_data$d.y[i])*Numextract(amca_data$n.y[i]))
  amca_data$f_1[i] <- sum(Numextract(amca_data$l_f.x[i]))
  amca_data$f_2[i] <- sum(Numextract(amca_data$l_f.y[i]))
}
amca_data <- amca_data[, c( "block", "plot","x","y",   "plant_id" , "startyear" ,  "biom_1"   ,
                            "biom_2" ,"comm.y",    "f_1"   ,    "f_2"  ,"cov.x" ,    "cov.y" , "alive_startyear.x" , "alive_startyear.y")]
names(amca_data) <- c( "block", "plot","x","y",   "plant_id" , "startyear" ,  "biom_1"   ,
                       "biom_2" , "comm_2",   "f_1"   ,    "f_2"   ,"cov_1" ,    "cov_2" ,"alive_1", "alive_2")

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

amca_data$f_1[which(!is.na(amca_data$biom_1) & is.na(amca_data$f_1))] <- 0 # replacing flower # of any plants with mreasurements but no flower number with a zero
amca_data$f_2[which(!is.na(amca_data$biom_2) & is.na(amca_data$f_2))] <- 0 
kueu_data$f_1[which(!is.na(kueu_data$biom_1) & is.na(kueu_data$f_1))] <- 0 
kueu_data$f_2[which(!is.na(kueu_data$biom_2) & is.na(kueu_data$f_2))] <- 0 
ecan_data$f_1[which(!is.na(ecan_data$biom_1) & is.na(ecan_data$f_1))] <- 0 
ecan_data$f_2[which(!is.na(ecan_data$biom_2) & is.na(ecan_data$f_2))] <- 0 

ecan_data$cov_1 <- as.numeric(ecan_data$ cov_1) # can ignore NAs introduced by coercion
ecan_data$cov_2 <- as.numeric(ecan_data$ cov_2)
ecan_data$block <- factor(ecan_data$block)
ecan_data$plot <- factor(ecan_data$plot)

kueu_data$cov_1 <- as.numeric(kueu_data$ cov_1)
kueu_data$cov_2 <- as.numeric(kueu_data$ cov_2)
kueu_data$block <- factor(kueu_data$block)
kueu_data$plot <- factor(kueu_data$plot)

amca_data$cov_1 <- as.numeric(amca_data$ cov_1)
amca_data$cov_2 <- as.numeric(amca_data$ cov_2)
amca_data$block <- factor(amca_data$block)
amca_data$plot <- factor(amca_data$plot)

write.csv(amca_data,file= "derived_data/05_amca_clean_demo_data.RData")
write.csv(ecan_data,file= "derived_data/05_ecan_clean_demo_data.RData")
write.csv(kueu_data,file= "derived_data/05_kueu_clean_demo_data.RData")

# recruitment----
# you want the code to search comm_2 for "qns" or "qs"
# if qns, return NA for recruitment counts
# if qs: 
  # sum all the reproductive effort between x = 100-200 and y= 100-200 for f_1
  # count all the plants that contain "new" in comm_2 but do not contain "no new" 
for (j in 1:3){
  
j_data <- get(paste(sp_name[j], "_data", sep= ""))
#if QNS: 
j_rec_data <- j_data[which(grepl("qns",j_data$comm_2 )),c("block", "plot", "startyear")]
j_rec_data$f_1 <- j_rec_data$rec_2 <- NA

# if QS: 
qses <- j_data[which(grepl("qs",j_data$comm_2 )),c("block", "plot", "startyear")] # here, startyear is year 1-- the year the fruiting occurred
qses$f_1 <- qses$rec_2 <- NA
for (i in 1:dim(qses)[1]){
  
  block_i <- qses$block[i]
  plot_i <- qses$plot[i]
  startyear_i <- qses$startyear[i]
  
  qses$f_1[i] <- sum(j_data$f_1[which(
    j_data$block == block_i & j_data$plot == plot_i & j_data$startyear == startyear_i & #same block, plot, year
      j_data$x >= 100 & j_data$x <= 200 & # x bounds
      j_data$x >= 100 & j_data$x <= 200)], # y bounds
    na.rm=TRUE)
   qses$rec_2[i]  <- dim(j_data[
     which( grepl("new", j_data$comm_2) & !grepl("no new", j_data$comm_2) & # counts number of new plants
     j_data$block == block_i & j_data$plot == plot_i & j_data$startyear == startyear_i & #same block, plot, year
       j_data$x >= 100 & j_data$x <= 200 & # x bounds
       j_data$x >= 100 & j_data$x <= 200), # y bounds
     ])[1] #gives you the count
                   
 
}
assign(paste(sp_name[j], "_rec", sep= ""),
       rbind.data.frame(j_rec_data, qses))
}

write.csv(amca_rec,file= "derived_data/05_amca_clean_rec_data.RData")
write.csv(ecan_rec,file= "derived_data/05_ecan_clean_rec_data.RData")
write.csv(kueu_rec,file= "derived_data/05_kueu_clean_rec_data.RData")

