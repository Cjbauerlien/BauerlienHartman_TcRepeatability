
##### Setting up the environment and data --------------------------------------

library(readxl)
library(reshape2)
library(dplyr)
library(lubridate)
library(doParallel)
library(lme4)
library(ggplot2)
library(ggpubr)
library(rptR)


getwd()

##### Estimating Tc ------------------------------------------------------------

# Calculation of Tc takes ~3 minutes when run in parallel with 7 cores. Rather than calculating Tc
# each time the code is run, the final data with calculated values were stored and are read in here
# (given the file exists in the directory). Additionally, some of the data points needed to be appended
# because two individuals (highlighted in Excel sheet) were missed during the first round of testing 
# (see Table 1). Therefore, the file labeled "Analyzed_Final" represents the complete and corrected 
# data following initial calculations of Tc and CTMax (referred to as CTM throughout code).


if(file.exists('Data_ThermalToleranceRepeatability_Analyzed_Final.xlsx')){
  dat <- read_excel('Data_ThermalToleranceRepeatability_Analyzed_Final.xlsx')
} else {
  
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
  
  ##### Estimating tolerance metrics for each fish in parallel -----------------
  
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
                 .packages = c('dplyr', 'reshape2', 'readxl', 'lubridate'),
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
                   
                   # Creates blank matrix to fill with iterative estimates of area between curves 
                   # based on different starting times (i.e., tc).
                   mat <- matrix(NA, nrow = nrow(Trial_dat), ncol = 1)
                   
                   # Fills blank matrix with estimates of area between two curves (i.e., curve 1 = 
                   # log(Temperature at time t) and the curve 2 = estimate of log(Tc))
                   for(i in 1:nrow(Trial_dat)){
                     mat[i,1] <- sum((Trial_dat$log.T[i:nrow(Trial_dat)] - Trial_dat$log.T[i]))
                   }
                   
                   # Appends the original data frame with estimates of back-transformed Tc (i.e.,
                   # when the starting temperature at which the area between the two curves equals
                   # zero), as well as estimates of CTM (temperature at LOE).
                   res <- matrix(NA, nrow = 1, ncol = 2)
                   if(is.na(dat$t_Trans[fish]) | is.na(dat$t_Rec[fish])){
                     res[, 1] <- NA } else {
                       res[, 1] <- exp(Trial_dat$log.T[min(which(mat[,1] <= 0))])  
                     }
                   
                   if(is.na(dat$t_Trans[fish])){
                     res[, 2] <- NA } else {
                       res[, 2] <- max(Trial_dat$Temp)
                     }
                   
                   return(res)
                 }
  
  end <- Sys.time()
  
  # Stop Cluster
  parallel::stopCluster(cl = my.cluster)
  
  # Time Elapsed
  end - start
  
  RES <- as.data.frame(RES)
  colnames(RES) <- c('Tc','CTM')
  dat <- cbind(dat, RES)
  
  dat$sc_Tc <-   scale(dat$Tc)
  dat$sc_CTM <-  scale(dat$CTM)
  
}

##### Fitting all combinations of linear mixed effects models for Tc -----------

### Fitting models for Tc 

# Additive
lmer_Tc_01 <- lmer(Tc ~ TestNumber + T0 + Electrofished + Wt + (1|Short_PIT), dat)
lmer_Tc_02 <- lmer(Tc ~ TestNumber +      Electrofished + Wt + (1|Short_PIT), dat)
lmer_Tc_03 <- lmer(Tc ~ TestNumber + T0 +                 Wt + (1|Short_PIT), dat)
lmer_Tc_04 <- lmer(Tc ~ TestNumber + T0 + Electrofished      + (1|Short_PIT), dat)
lmer_Tc_05 <- lmer(Tc ~ TestNumber +                      Wt + (1|Short_PIT), dat)
lmer_Tc_06 <- lmer(Tc ~ TestNumber +      Electrofished      + (1|Short_PIT), dat)
lmer_Tc_07 <- lmer(Tc ~ TestNumber + T0                      + (1|Short_PIT), dat)
lmer_Tc_08 <- lmer(Tc ~ TestNumber                           + (1|Short_PIT), dat)
lmer_Tc_09 <- lmer(Tc ~              T0 + Electrofished + Wt + (1|Short_PIT), dat)
lmer_Tc_10 <- lmer(Tc ~                   Electrofished + Wt + (1|Short_PIT), dat)
lmer_Tc_11 <- lmer(Tc ~              T0 +                 Wt + (1|Short_PIT), dat)
lmer_Tc_12 <- lmer(Tc ~              T0 + Electrofished      + (1|Short_PIT), dat)
lmer_Tc_13 <- lmer(Tc ~                                   Wt + (1|Short_PIT), dat)
lmer_Tc_14 <- lmer(Tc ~                   Electrofished      + (1|Short_PIT), dat)
lmer_Tc_15 <- lmer(Tc ~              T0                      + (1|Short_PIT), dat)
lmer_Tc_16 <- lmer(Tc ~ 1                                    + (1|Short_PIT), dat)

# Interactions
lmer_Tc_17 <- lmer(Tc ~ TestNumber + T0 + Electrofished + Wt + TestNumber:T0            + (1|Short_PIT), dat)
lmer_Tc_18 <- lmer(Tc ~ TestNumber + T0 + Electrofished + Wt + TestNumber:Electrofished + (1|Short_PIT), dat)
lmer_Tc_19 <- lmer(Tc ~ TestNumber + T0 + Electrofished + Wt + TestNumber:Wt            + (1|Short_PIT), dat)
lmer_Tc_20 <- lmer(Tc ~ TestNumber + T0 + Electrofished + Wt + T0:Electrofished         + (1|Short_PIT), dat)
lmer_Tc_21 <- lmer(Tc ~ TestNumber + T0 + Electrofished + Wt + T0:Wt                    + (1|Short_PIT), dat)
lmer_Tc_22 <- lmer(Tc ~ TestNumber + T0 + Electrofished + Wt + Electrofished:Wt         + (1|Short_PIT), dat)
lmer_Tc_23 <- lmer(Tc ~ TestNumber +      Electrofished + Wt + TestNumber:Electrofished + (1|Short_PIT), dat)
lmer_Tc_24 <- lmer(Tc ~ TestNumber +      Electrofished + Wt + TestNumber:Wt            + (1|Short_PIT), dat)
lmer_Tc_25 <- lmer(Tc ~ TestNumber +      Electrofished + Wt + Electrofished:Wt         + (1|Short_PIT), dat)
lmer_Tc_26 <- lmer(Tc ~ TestNumber + T0 +                 Wt + TestNumber:T0            + (1|Short_PIT), dat)
lmer_Tc_27 <- lmer(Tc ~ TestNumber + T0 +                 Wt + TestNumber:Wt            + (1|Short_PIT), dat)
lmer_Tc_28 <- lmer(Tc ~ TestNumber + T0 +                 Wt + T0:Wt                    + (1|Short_PIT), dat)
lmer_Tc_29 <- lmer(Tc ~ TestNumber + T0 + Electrofished +      TestNumber:T0            + (1|Short_PIT), dat)
lmer_Tc_30 <- lmer(Tc ~ TestNumber + T0 + Electrofished +      TestNumber:Electrofished + (1|Short_PIT), dat)
lmer_Tc_31 <- lmer(Tc ~ TestNumber + T0 + Electrofished +      T0:Electrofished         + (1|Short_PIT), dat)
lmer_Tc_32 <- lmer(Tc ~              T0 + Electrofished + Wt + T0:Electrofished         + (1|Short_PIT), dat)
lmer_Tc_33 <- lmer(Tc ~              T0 + Electrofished + Wt + T0:Wt                    + (1|Short_PIT), dat)
lmer_Tc_34 <- lmer(Tc ~              T0 + Electrofished + Wt + Electrofished:Wt         + (1|Short_PIT), dat)
lmer_Tc_35 <- lmer(Tc ~ TestNumber + T0 +                      TestNumber:T0            + (1|Short_PIT), dat)
lmer_Tc_36 <- lmer(Tc ~ TestNumber +      Electrofished +      TestNumber:Electrofished + (1|Short_PIT), dat)
lmer_Tc_37 <- lmer(Tc ~ TestNumber +                      Wt + TestNumber:Wt            + (1|Short_PIT), dat)
lmer_Tc_38 <- lmer(Tc ~              T0 + Electrofished +      T0:Electrofished         + (1|Short_PIT), dat)
lmer_Tc_39 <- lmer(Tc ~              T0 +                 Wt + T0:Wt                    + (1|Short_PIT), dat)
lmer_Tc_40 <- lmer(Tc ~                   Electrofished + Wt + Electrofished:Wt         + (1|Short_PIT), dat)


### Fitting models for CTM

# Additive
lmer_CTM_01 <- lmer(CTM ~ TestNumber + T0 + Electrofished + Wt + (1|Short_PIT), dat)
lmer_CTM_02 <- lmer(CTM ~ TestNumber +      Electrofished + Wt + (1|Short_PIT), dat)
lmer_CTM_03 <- lmer(CTM ~ TestNumber + T0 +                 Wt + (1|Short_PIT), dat)
lmer_CTM_04 <- lmer(CTM ~ TestNumber + T0 + Electrofished      + (1|Short_PIT), dat)
lmer_CTM_05 <- lmer(CTM ~ TestNumber +                      Wt + (1|Short_PIT), dat)
lmer_CTM_06 <- lmer(CTM ~ TestNumber +      Electrofished      + (1|Short_PIT), dat)
lmer_CTM_07 <- lmer(CTM ~ TestNumber + T0                      + (1|Short_PIT), dat)
lmer_CTM_08 <- lmer(CTM ~ TestNumber                           + (1|Short_PIT), dat)
lmer_CTM_09 <- lmer(CTM ~              T0 + Electrofished + Wt + (1|Short_PIT), dat)
lmer_CTM_10 <- lmer(CTM ~                   Electrofished + Wt + (1|Short_PIT), dat)
lmer_CTM_11 <- lmer(CTM ~              T0 +                 Wt + (1|Short_PIT), dat)
lmer_CTM_12 <- lmer(CTM ~              T0 + Electrofished      + (1|Short_PIT), dat)
lmer_CTM_13 <- lmer(CTM ~                                   Wt + (1|Short_PIT), dat)
lmer_CTM_14 <- lmer(CTM ~                   Electrofished      + (1|Short_PIT), dat)
lmer_CTM_15 <- lmer(CTM ~              T0                      + (1|Short_PIT), dat)
lmer_CTM_16 <- lmer(CTM ~ 1                                    + (1|Short_PIT), dat)

# Interactions
lmer_CTM_17 <- lmer(CTM ~ TestNumber + T0 + Electrofished + Wt + TestNumber:T0            + (1|Short_PIT), dat)
lmer_CTM_18 <- lmer(CTM ~ TestNumber + T0 + Electrofished + Wt + TestNumber:Electrofished + (1|Short_PIT), dat)
lmer_CTM_19 <- lmer(CTM ~ TestNumber + T0 + Electrofished + Wt + TestNumber:Wt            + (1|Short_PIT), dat)
lmer_CTM_20 <- lmer(CTM ~ TestNumber + T0 + Electrofished + Wt + T0:Electrofished         + (1|Short_PIT), dat)
lmer_CTM_21 <- lmer(CTM ~ TestNumber + T0 + Electrofished + Wt + T0:Wt                    + (1|Short_PIT), dat)
lmer_CTM_22 <- lmer(CTM ~ TestNumber + T0 + Electrofished + Wt + Electrofished:Wt         + (1|Short_PIT), dat)
lmer_CTM_23 <- lmer(CTM ~ TestNumber +      Electrofished + Wt + TestNumber:Electrofished + (1|Short_PIT), dat)
lmer_CTM_24 <- lmer(CTM ~ TestNumber +      Electrofished + Wt + TestNumber:Wt            + (1|Short_PIT), dat)
lmer_CTM_25 <- lmer(CTM ~ TestNumber +      Electrofished + Wt + Electrofished:Wt         + (1|Short_PIT), dat)
lmer_CTM_26 <- lmer(CTM ~ TestNumber + T0 +                 Wt + TestNumber:T0            + (1|Short_PIT), dat)
lmer_CTM_27 <- lmer(CTM ~ TestNumber + T0 +                 Wt + TestNumber:Wt            + (1|Short_PIT), dat)
lmer_CTM_28 <- lmer(CTM ~ TestNumber + T0 +                 Wt + T0:Wt                    + (1|Short_PIT), dat)
lmer_CTM_39 <- lmer(CTM ~ TestNumber + T0 + Electrofished +      TestNumber:T0            + (1|Short_PIT), dat)
lmer_CTM_30 <- lmer(CTM ~ TestNumber + T0 + Electrofished +      TestNumber:Electrofished + (1|Short_PIT), dat)
lmer_CTM_31 <- lmer(CTM ~ TestNumber + T0 + Electrofished +      T0:Electrofished         + (1|Short_PIT), dat)
lmer_CTM_32 <- lmer(CTM ~              T0 + Electrofished + Wt + T0:Electrofished         + (1|Short_PIT), dat)
lmer_CTM_33 <- lmer(CTM ~              T0 + Electrofished + Wt + T0:Wt                    + (1|Short_PIT), dat)
lmer_CTM_34 <- lmer(CTM ~              T0 + Electrofished + Wt + Electrofished:Wt         + (1|Short_PIT), dat)
lmer_CTM_35 <- lmer(CTM ~ TestNumber + T0 +                      TestNumber:T0            + (1|Short_PIT), dat)
lmer_CTM_36 <- lmer(CTM ~ TestNumber +      Electrofished +      TestNumber:Electrofished + (1|Short_PIT), dat)
lmer_CTM_37 <- lmer(CTM ~ TestNumber +                      Wt + TestNumber:Wt            + (1|Short_PIT), dat)
lmer_CTM_38 <- lmer(CTM ~              T0 + Electrofished +      T0:Electrofished         + (1|Short_PIT), dat)
lmer_CTM_39 <- lmer(CTM ~              T0 +                 Wt + T0:Wt                    + (1|Short_PIT), dat)
lmer_CTM_40 <- lmer(CTM ~                   Electrofished + Wt + Electrofished:Wt         + (1|Short_PIT), dat)



# Creates list of all models in the environment that start with "lmer_"
lmer_Tc.model.list   <- mget(ls()[startsWith(ls(),'lmer_Tc_')])
lmer_CTM.model.list  <- mget(ls()[startsWith(ls(),'lmer_CTM_')])

# Creates blank data frame to be filled with AIC values for each lmer model
MOD_COMP_Tc   <- matrix(NA, nrow = length(lmer_Tc.model.list), ncol = 4)
MOD_COMP_CTM  <- matrix(NA, nrow = length(lmer_CTM.model.list), ncol = 4)

# Fills the blank matrix with AIC, delAIC, and df values, and names each model as a row
MOD_COMP_Tc[,1] <- as.numeric(unlist(lapply(lmer_Tc.model.list, FUN = AIC)))
MOD_COMP_Tc[,2] <- MOD_COMP_Tc[,1] - min(MOD_COMP_Tc[,1])
MOD_COMP_Tc[,3] <- as.numeric(unlist(lapply(lmer_Tc.model.list, FUN = df.residual)))
MOD_COMP_Tc[,3] <- nrow(dat[!is.na(dat$Tc),]) - MOD_COMP_Tc[,3] # workaround to calculate model df from residual df
rownames(MOD_COMP_Tc) <- names(lmer_Tc.model.list)
MOD_COMP_Tc <- as.data.frame(MOD_COMP_Tc)
for(i in 1:nrow(MOD_COMP_Tc)){MOD_COMP_Tc[i,4] <- paste(formula(get(rownames(MOD_COMP_Tc)[i]))[c(2,1,3)],collapse = '')}

MOD_COMP_CTM[,1] <- as.numeric(unlist(lapply(lmer_CTM.model.list, FUN = AIC)))
MOD_COMP_CTM[,2] <- MOD_COMP_CTM[,1] - min(MOD_COMP_CTM[,1])
MOD_COMP_CTM[,3] <- as.numeric(unlist(lapply(lmer_CTM.model.list, FUN = df.residual)))
MOD_COMP_CTM[,3] <- nrow(dat[!is.na(dat$CTM),]) - MOD_COMP_CTM[,3] # workaround to calculate model df from residual df
rownames(MOD_COMP_CTM) <- names(lmer_CTM.model.list)
MOD_COMP_CTM <- as.data.frame(MOD_COMP_CTM)
for(i in 1:nrow(MOD_COMP_CTM)){MOD_COMP_CTM[i,4] <- paste(formula(get(rownames(MOD_COMP_CTM)[i]))[c(2,1,3)],collapse = '')}

# Arranges matrix so that lowest AIC scores are at the top
MOD_COMP_Tc   <- MOD_COMP_Tc[order(MOD_COMP_Tc[,1]),]
MOD_COMP_CTM  <- MOD_COMP_CTM[order(MOD_COMP_CTM[,1]),]

# Selects the model with the lowest AIC score and stores it as "Top_Model"
Top_Model_Tc   <- get(rownames(MOD_COMP_Tc)[1])
Top_Model_CTM  <- get(rownames(MOD_COMP_CTM)[1])

summary(Top_Model_Tc)
summary(Top_Model_CTM)

##### Estimating repeatability -------------------------------------------------

# NOTE: 'rptR' calculates the intra-class correlation coefficient, which is estimated
# as (variance of interest) / (variance of interest + unwanted variance), where
# the "variance of interest" is the inter-individual variation (i.e., ability of
# an individual to maintain its rank in the analysis) and "unwanted variance"
# is the intra-individual variation between tests.

set.seed(1634)
rpt.Tc   <- rpt(formula = formula(Top_Model_Tc),   data = dat, parallel = T, datatype = 'Gaussian', grname = 'Short_PIT')
rpt.CTM  <- rpt(formula = formula(Top_Model_CTM),  data = dat, parallel = T, datatype = 'Gaussian', grname = 'Short_PIT')

summary(rpt.Tc)
summary(rpt.CTM)

##### Writing out final data frame ---------------------------------------------

# library(writexl)
# write_xlsx(dat, 'Data_ThermalToleranceRepeatability_Analyzed.xlsx')
