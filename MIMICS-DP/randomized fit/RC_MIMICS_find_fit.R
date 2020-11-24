
library(ggplot2)
library(dplyr)

# set working drive
setwd("C:/GitHub/MIMICS-RCrk/MIMICS-DP/")

# get MIMICS ftn
source("MIMICS_ftn.R")

# get litter characteristics
litter <- read.csv("Litter_characteristics/LIDET_LitterCharacteristics.csv")
litter <- litter[1:6,]   #subset foliar litter only, DP: last three are roots? 

# get site data
data <- read.csv("Site_data/RC_forMIMICS_all.csv", as.is=T) #site level forcing variables

### MODIFY DATA ###
# Remove Johnston Draw
data <- data %>% filter(L1 != "Johnston Draw") 

# Increase north facing precip
adj_North_MAP <- 2
data <- data %>%
  mutate(MAP = ifelse(ASPECT_CLASS == "North" |
                        ASPECT_CLASS == "Northwest" |
                        ASPECT_CLASS == "Northeast",
                      MAP*adj_North_MAP, MAP))

# Decrease north facing temp
adj_North_MAT <- 0.7
data <- data %>%
  mutate(MAT = ifelse(ASPECT_CLASS == "North" |
                        ASPECT_CLASS == "Northwest" |
                        ASPECT_CLASS == "Northeast",
                      MAT*adj_North_MAT, MAT))
#########

# Get parameters
parameters <- read.csv("MIMICS_parameters/parameters_LIDET-MIM-REV.csv")

#Parameters as list
parameter_list <- as.list(parameters)
parameter_list <- lapply(parameter_list, function(x) x[!is.na(x)])


#Create a randomized parameter list
#Starter parameters as list
parameter_list <- as.list(parameters)
parameter_list['rep'] <- 0
parameter_list <- lapply(parameter_list, function(x) x[!is.na(x)])


parm_ball <- list()
j=1

while(j<100) {
  
  add_parm <- parameter_list
  add_parm[['rep']] <- j
  add_parm[['tau_r']][1] <- runif(1,1,100)/100
  parm_ball[[j]] <- add_parm
  
  #iterate
  j=j+1
}




# Apply parameter lists to MIMICS and collect r^2 and fit
#df <- lapply(run_MIMICS(data, litter, parameter_list)

fit_MIMICS <- lapply(parm_ball, run_MIMICS, litter_in=litter, data_in = data)
fit_results <- do.call(rbind, fit_MIMICS)
fit_results
