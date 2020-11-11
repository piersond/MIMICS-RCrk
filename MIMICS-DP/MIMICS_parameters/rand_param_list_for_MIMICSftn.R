### MAKING RANDOM PARAMETER LISTS ###
library(dplyr)

# set working drive
setwd("C:/GitHub/MIMICS-RCrk/MIMICS-DP/")

# Get parameters
parameters <- read.csv("MIMICS_parameters/parameters_LIDET-MIM-REV.csv")

#Starter parameters as list
parameter_list <- as.list(parameters)
parameter_list <- lapply(parameter_list, function(x) x[!is.na(x)])


parm_ball <- list()
i=1

while(i<100) {
  
  new_parms <- parameter_list
  new_parms[['tau_r']][1] <- runif(1,1,100)/100000
  parm_ball[[i]] <- new_parms
  
  #iterate
  i=i+1
}



