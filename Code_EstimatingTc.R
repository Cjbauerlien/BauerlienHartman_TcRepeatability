
##### Setting up the environment and data --------------------------------------

library(readxl)
library(reshape2)
library(dplyr)
library(lubridate)
library(doParallel)
library(ggplot2)
library(ggpubr)
library(gcplyr)


getwd()

# Calculation of Tc takes ~3 minutes when run in parallel with 7 cores. Rather than calculating Tc
# each time, the final data with calculated values were stored and are read in here (given the file
# exists in the directory). 

# Some of the data points needed to be appended because two individuals were missed during the first
# round of testing. Additionally, some individuals needed to be censored due to injury. Therefore, 
# the file with "Analyzed_Final" represents the complete and corrected data following initial calculations 
# of Tc and CTM.


TestID <- read_xlsx('Data_ThermalToleranceRepeatability_Trial.xlsx', sheet = 'TestID', na = 'NA') %>% filter(!is.na(TestNumber))
Trial  <- read_xlsx('Data_ThermalToleranceRepeatability_Trial.xlsx', sheet = 'Trial',  na = 'NA')

DAT <- merge(TestID, Trial,                         # Merges data frames 
             by = 'TestID',                         # By column Test ID
             suffixes = c('.TestID','.Trial'),      # Uses suffix to label repeated column names
             all.x = T)                             # Only keeping data for valid tests
  
dat <- DAT %>%
  mutate(Acc_t = difftime(StartTime, TimeAtHoldTank, units = 'min'),  # Calculating acclimation time
         Rec_t = t_Rec - t_Trans,                                     # Calculating recovery time
         Electrofished = as.factor(Electrofished),                    # Sets variables as factors
         SourceTank = as.factor(SourceTank),
         Short_PIT = as.factor(Short_PIT)) %>%
  filter(is.na(Mortality))                          # Removes all instances of test-related mortalities

##### Estimating tolerance metrics for each fish in parallel ---------------------

  # Set how many cores to use
  n.cores <- parallel::detectCores() - 1
  
  # Create the cluster
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  
  # Register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  
  start <- Sys.time()
  RES <- foreach(fish = 1:nrow(dat),
                 .packages = c('dplyr', 'reshape2', 'readxl', 'lubridate', 'gcplyr'),
                 .combine = 'rbind',
                 .inorder = T) %dopar% {
                   if(between(dat$Date[fish], as.Date('2023-03-12'), as.Date('2023-11-05')) | between(dat$Date[fish], as.Date('2024-03-10'), as.Date('2024-11-03'))){
                     Heat_Tank <- read_excel(paste('ThermalToleranceTemperatureData/', dat$HT_File[fish], ' EDT (Data EDT).xlsx', sep = ''))[,-1]; colnames(Heat_Tank) <- c('datetime', 'Temp')
                     Reco_Tank <- read_excel(paste('ThermalToleranceTemperatureData/', dat$RT_File[fish], ' EDT (Data EDT).xlsx', sep = ''))[,-1]; colnames(Reco_Tank) <- c('datetime', 'Temp')
                   } else {
                     Heat_Tank <- read_excel(paste('ThermalToleranceTemperatureData/', dat$HT_File[fish], ' EST (Data EST).xlsx', sep = ''))[,-1]; colnames(Heat_Tank) <- c('datetime', 'Temp')
                     Reco_Tank <- read_excel(paste('ThermalToleranceTemperatureData/', dat$RT_File[fish], ' EST (Data EST).xlsx', sep = ''))[,-1]; colnames(Reco_Tank) <- c('datetime', 'Temp')
                   }
                   
                   # Creates data frame of temperatures from both heating and recovery tanks
                   Temp_dat <- merge(Heat_Tank, Reco_Tank, 
                                     by = 'datetime', 
                                     suffixes = c('Heat', 'Reco'), 
                                     all = T) 
                   
                   # Filters combined temperature data to only include time range during individual
                   # fish trials, then assigns a temperature (i.e., heating tank or recovery tank)
                   # based on the time of transfer. Finally, temperature is log-transformed and 
                   # time over the test is used as an index for next step.
                   if(!is.na(dat$t_Rec[fish])){
                     
                     Trial_dat <- Temp_dat %>%
                       filter(datetime >= dat$StartTime[fish] & datetime <= dat$t_Rec[fish]) %>%
                       rowwise() %>%
                       mutate(Temp = ifelse(datetime < dat$t_Trans[fish], yes = TempHeat, no  = TempReco),
                              log.T = log(Temp),
                              t = difftime(datetime, dat$StartTime[fish], units = 'sec'))
                   } else {
                     Trial_dat <- Temp_dat %>%
                       filter(datetime >= dat$StartTime[fish] & datetime <= dat$EndTime[fish]) %>%
                       rowwise() %>%
                       mutate(Temp = ifelse(datetime < dat$t_Trans[fish], yes = TempHeat, no  = TempReco),
                              log.T = log(Temp),
                              t = difftime(datetime, dat$StartTime[fish], units = 'sec'))
                   }
                   
                   # Find index when fish exhibited LOE (used to find CTM and censor recovery temperatures
                   # during iterative minimization procedure).
                   if(is.na(dat$t_Trans[fish]) | is.na(dat$t_Rec[fish])){
                     ind_LOE <- NA} else {
                       ind_LOE <- min(which(Trial_dat$datetime == dat$t_Trans[fish]), na.rm = T)
                     }
                   
                   # Creates blank matrix to fill with iterative estimates of area between curves 
                   # based on different starting times.
                   mat <- matrix(NA, nrow = nrow(Trial_dat), ncol = 3)
                   
                   # Fills blank matrix with estimates of area between two curves (i.e., curve 1 = 
                   # log(Temperature at time t) and the curve 2 = estimate of log(Tc))
                   for(i in 1:nrow(Trial_dat)){
                     mat[i,1] <- Trial_dat$log.T[i] - sum(Trial_dat$log.T[i:nrow(Trial_dat)]) / (nrow(Trial_dat) - i)
                   }
                   
                   mat[,2] <- abs(mat[,1])
                   
                   # Computes the moving slope of the area function across a 20 unit window of indices
                   # to verify that the area function reached a true minimum (i.e., where the sign of the derivative
                   # switches from negative to positive).
                   for(i in 11:(nrow(Trial_dat)-11)){
                     mat[i,3] <- sum((mat[(i-10):(i+10),2] - mean(mat[c((i-10):(i+10)),2])) * (c((i - 10):(i + 10)) - mean(c((i - 10):(i + 10))))) / sum((c((i - 10):(i + 10)) - mean(c((i - 10):(i + 10))))^2)
                       
                   }
                   
                   # Identifies the index of the temperature function where the area under the curve
                   # function is minimized over the range of time indices (i.e., when the temperature
                   # function reaches Tc) without influence from the endpoints or when the fish was transferred
                   # to the recovery tank.
                   if(is.na(dat$t_Trans[fish]) | is.na(dat$t_Rec[fish])){
                     ind_Tc <- NA} else {
                       extrema <- find_local_extrema(mat[1:ind_LOE,2], window_width_n = 301, return_endpoints = F, return_minima = T, return_maxima = F)
                       ind_Tc <- extrema[which(mat[extrema,2] == min(mat[extrema, 2], na.rm = T))] # Absolute minimum
                       
                     }
                   
                   # Appends the original data frame with estimates of back-transformed Tc (i.e.,
                   # when the starting temperature at which the area between the two curves equals
                   # zero), as well as estimates of CTM (time at LOE).
                   res <- matrix(NA, nrow = 1, ncol = 4)
                   
                   # Estimate of Tc
                   if(is.na(dat$t_Trans[fish]) | is.na(dat$t_Rec[fish])){
                     res[, 1] <- NA } else {
                       res[, 1] <- exp(Trial_dat$log.T[ind_Tc])
                     }
                   
                   # Estimate of CTmax
                   if(is.na(dat$t_Trans[fish])){
                     res[, 2] <- NA } else {
                       res[, 2] <- max(Trial_dat$Temp)
                     }
                   
                   # Area under the curve to make sure estimation procedure worked for each fish
                   if(is.na(dat$t_Trans[fish]) | is.na(dat$t_Rec[fish])){
                     res[, 3] <- NA } else {
                       res[, 3] <- mat[ind_Tc, 2]
                     }
                  
                   # Storing indices to evaluate problem ranges during minimization
                   if(is.na(dat$t_Trans[fish]) | is.na(dat$t_Rec[fish])){
                     res[, 4] <- NA } else {
                       res[, 4] <- ind_Tc
                     }
                   
                   return(res)
                 }
  
  end <- Sys.time()
  
  # Stop Cluster
  parallel::stopCluster(cl = my.cluster)
  
  # Time Elapsed
  end - start
  
  RES <- as.data.frame(RES)
  colnames(RES) <- c('Tc','CTM','AreaUnderCurve','Ind')
  dat <- cbind(dat, RES)
  
  dat$sc_Tc <-   scale(dat$Tc)
  dat$sc_CTM <-  scale(dat$CTM)


##### Writing out final data frame ---------------------------------------------
stop('Make sure you are ready to write the final data')

# library(writexl)
# write_xlsx(dat, 'Data_ThermalToleranceRepeatability_Analyzed.xlsx')
