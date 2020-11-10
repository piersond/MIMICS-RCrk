library(rootSolve)
library(boot)
library(ggplot2)


rm(list=ls())

# set working drive
setwd("C:/GitHub/MIMICS-RCrk/MIMICS-DP/")

# Gather req'd data, parameters, etc.
############################################################################################
# get RXEQ ftn
source("RXEQ_ftn.R")

# get MIMICS parameters
list.files("MIMICS_parameters")
para_file <- "parameters_LIDET-MIM-REV.csv"  # or use: list.files("Parameters")[#]
parameters <- read.csv(paste0("MIMICS_parameters/",para_file))
attach(parameters)
depth <- parameters$depth[1] 

# get litter characteristics
litter <- read.csv("Litter_characteristics/LIDET_LitterCharacteristics.csv")
litter <- litter[1:6,]   #subset foliar litter only, DP: last three are roots? 
attach(litter)

# get site data
data <- read.csv("Site_data/RC_forMIMICS_all.csv", as.is=T) #site level forcing variables

### MODIFY DATA ###

#Remove Johnston Draw
data <- data %>% filter(L1 != "Johnston Draw") 


#Increase north facing precip
adj_North_MAP <- 2
data <- data %>%
  mutate(MAP = ifelse(ASPECT_CLASS == "North" |
                      ASPECT_CLASS == "Northwest" |
                      ASPECT_CLASS == "Northeast",
                      MAP*adj_North_MAP, MAP))

#Decrease north facing temp
adj_North_MAT <- 0.7
data <- data %>%
  mutate(MAT = ifelse(ASPECT_CLASS == "North" |
                        ASPECT_CLASS == "Northwest" |
                        ASPECT_CLASS == "Northeast",
                      MAT*adj_North_MAT, MAT))

attach(data)

# Create/ set/ edit pools and parameters
############################################################################################

# revise model paramters from file
# Vslope[1:6] <- Vslope[1:6] 
# Vint[1:6] <- Vint[1:6]
# Kslope[1:6] <- Kslope[1:6]
# Kint[1:6] <- Kint[1:6]
# aK[1:6] <- aK[1:6]
# vMOD[1:6] <- vMOD[1:6]
# kMOD[1:6] <- kMOD[1:6]
# CUE[1:6] <- CUE[1:6]
#tau_r[1] <- tau_r[1] * 0.13#tau_r[1] is base value
#tau_r[2] <- tau_r[2] * 0.1  # tau_r[2] is e^ value

#tau_k[1] <- tau_k[1] * 0.5
#tau_k[2] <- tau_k[2] * 0.05

#Tau_MOD[1] <- Tau_MOD[1] * 0.5
#Tau_MOD[2] <- Tau_MOD[2] * 0.5
#Tau_MOD[3] <- Tau_MOD[3] * 0.5
#Tau_MOD[4] <- Tau_MOD[4] * 0.5

# fPHYS_r[1:6] <- fPHYS_r[1:6]
# fPHYS_K[1:6] <- fPHYS_K[1:6]
# fCHEM_r[1:6] <- fCHEM_r[1:6]
# fCHEM_K[1:6] <- fCHEM_K[1:6]
# fSOM_p[1:6] <- fSOM_p[1:6]
# PHYS_scalar[1:6] <- PHYS_scalar[1:6]
# FI[1:6] <- FI[1:6]
# fmet_p[1:6] <- fmet_p[1:6]

calcN    <- (1 / litCN) / 2.5 * 100  #DP: what's this calc N for?    
lit_fMET <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * litLIG / calcN) 

KO <- parameters$KO
ANPP  <- data$ANPP / 2           		# if needed convert to gC/m2/y from g/m2/y
clay  <- data$CLAY2/100  				    # if needed, convert from % clay to fraction
tsoi  <- MAT
nsites<- length(Site)

lig    <- LIG #/ 100
Nnew   <- N                                                  #N in litter additions
fMET1  <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * lig / Nnew)   #as partitioned in Daycent
MIMLIT <- rep(NA, nsites)           	                       #Vector for results
MIMMIC <- rep(NA, nsites)           
MIM_CO <- rep(NA, nsites)           
MIMSOC <- rep(NA, nsites)           

strSite <- as.character(data$Site)                           #convert site names to string
pools  <- c('site','LITm', 'LITs', 'MICr', 'MICK', 'SOMp', 'SOMc', 'SOMa')
table  <- array(NA, dim=c(nsites,8), dimnames=list(as.character(Site),pools)) 
table[,1]       <- as.character(Site)


POOLS  <- c("LIT","MIC","SOC")
npools <- length(POOLS)  
LITpool<- c('LITm', 'LITs') 
MICpool<- c('MICr', 'MICk') 
SOMpool<- c('SOMp', 'SOMa','SOMc') 

MIMLIT <- array(NA, dim=c(nsites))           	             # array for results
MIMMIC <- array(NA, dim=c(nsites))           
MIM_CO <- array(NA, dim=c(nsites))           
MIMSOC <- array(NA, dim=c(nsites))           

strSite <- as.character(data$Site)  #convert site names to string

#Make vectors to store model results that match obs temporal resolution
npts   <- 6*10*14   				#6 litter * 10 years * 14 sites
xyLIT  <- rep(NA, npts) 
xyTIME <- rep(NA, npts) 
xySITE <- rep(NA, npts) 
xyOBS  <- rep(NA, npts)
xyMIM  <- rep(NA, npts)
xyCount<- 1


# (B)       RXEQ for site using STODE function
# Starts Big loop over all sites (i) & all litter types (j)
#############################################################################
for (i in 1:nsites) {         #speeds up debugging     	
  # Read in site characteristics -----------
  print(paste("-------- starting ", data$Site[i], " --------") )
  
  ###debug
  #i <- 1
  
  fMET       <- mean(lit_fMET)           # uses mean litter fmet from LIDET
  # fMET       <- fMET1[i]                 # uses site estimate for fmet
  fCLAY      <- clay[i]
  TSOI       <- tsoi[i]
  EST_LIT_in <- ANPP[i] / (365*24)         # gC/m2/h (from gC/m2/y)
  h2y        <- 24*365
  MICROtoECO <- depth * 1e4 * 1e-3         # mgC/cm3 to g/m2
  EST_LIT    <- EST_LIT_in  * 1e3 / 1e4    # mgC/cm2/h(from gC/m2/h) 
  
  # ------------ caclulate parameters ---------------
  Vmax     <- exp(TSOI * Vslope + Vint) * aV
  Km       <- exp(TSOI * Kslope + Kint) * aK
  
  #ANPP strongly correlated with MAP
   #Adjusts tau to different set values based on ANPP
  Tau_MOD1 <- sqrt(ANPP[i]/Tau_MOD[1])          # basicaily standardize against NWT
  Tau_MOD2 <- Tau_MOD[4]                        # increased 3-fold for SS SOC pools
  Tau_MOD1[Tau_MOD1 < Tau_MOD[2]] <- Tau_MOD[2] # correction not used in LIDET resutls 
  Tau_MOD1[Tau_MOD1 > Tau_MOD[3]] <- Tau_MOD[3] 
  tau <- c(tau_r[1]*exp(tau_r[2]*fMET), 
           tau_K[1]*exp(tau_K[2]*fMET))   
  tau <- tau * Tau_MOD1 * Tau_MOD2
  
  fPHYS    <- c(fPHYS_r[1] * exp(fPHYS_r[2]*fCLAY), 
                fPHYS_K[1] * exp(fPHYS_K[2]*fCLAY)) 	            #fraction to SOMp
  fCHEM    <- c(fCHEM_r[1] * exp(fCHEM_r[2]*fMET) * fCHEM_r[3], 
                fCHEM_K[1] * exp(fCHEM_K[2]*fMET) * fCHEM_K[3]) 	#fraction to SOMc
  fAVAI    <- 1- (fPHYS + fCHEM)
  desorb   <- fSOM_p[1] * exp(fSOM_p[2]*(fCLAY))                  #CHANGED FOR GLOBAL RUN!!!     
  
  #desorb   <- desorb/10 # modified as in MIMdef from Zhang et al 2020
  #fPHYS    <- fPHYS/5  # to reduce allocation to physically protected pool 5x
  
  pSCALAR  <- PHYS_scalar[1] * exp(PHYS_scalar[2]*(sqrt(fCLAY)))  #Scalar for texture effects on SOMp
  v_MOD    <- vMOD  # to avoid writing over orig. parameters
  k_MOD    <- kMOD 
  k_MOD[3] <- k_MOD[3] * pSCALAR    
  k_MOD[6] <- k_MOD[6] * pSCALAR    
  
  VMAX     <- Vmax * v_MOD 
  KM       <- Km / k_MOD
  
  #initialize pools
  I       <- array(NA, dim=2)              #Litter inputs to MET/STR
  I[1]    <- (EST_LIT / depth) * fMET     
  I[2]    <- (EST_LIT / depth) * (1-fMET)
  lit     <- I   
  mic     <- I  
  som     <- rep(NA, 3) 
  som[1]  <- I[1]
  som[2]  <- I[2]
  som[3]  <- I[1] 
  LITmin  <- rep(NA, dim=4)
  MICtrn  <- rep(NA, dim=6)
  SOMmin  <- rep(NA, dim=2)
  DEsorb  <- rep(NA, dim=1)
  OXIDAT  <- rep(NA, dim=1)
  
  #Calculate RXEQ pools  
  Tpars <- c( I = I, VMAX = VMAX, KM = KM, CUE = CUE, 
              fPHYS = fPHYS, fCHEM = fCHEM, fAVAI = fAVAI, FI = FI, 
              tau = tau, LITmin = LITmin, SOMmin = SOMmin, MICtrn = MICtrn, 
              desorb = desorb, DEsorb = DEsorb, OXIDAT = OXIDAT, KO = KO)
  Ty    <- c( LIT_1 = lit[1], LIT_2 = lit[2], 
              MIC_1 = mic[1], MIC_2 = mic[2], 
              SOM_1 = som[1], SOM_2 = som[2], SOM_3 = som[3] )
  test  <- stode(y = Ty, time = 1e6, fun = RXEQ, parms = Tpars, positive = TRUE)
  if (i == 2) {print(test[1])}
  remove(lit, mic, som)
  
  table[i,2:8] <- as.numeric(test[[1]])
  MIMLIT[i]    <- (test[[1]][[1]]+test[[1]][[2]])  * depth *1e4 / 1e6#convert kgC/m2 from mgC/cm3 (0-30 cm) 
  MIMMIC[i]    <- (test[[1]][[3]]+test[[1]][[4]])  * depth *1e4 / 1e6
  MIM_CO[i]    <-  test[[1]][[3]]/test[[1]][[4]]
  MIMSOC[i]    <- sum(test[[1]])  * depth *1e4 / 1e6   
  
  remove(test, Ty, Tpars, LITmin, MICtrn, SOMmin)
  
}


# Plot results
###############################################################################
data$MIMSOC <- MIMSOC
data$MIMIC <- MIMMIC

# correlation
r_test <- cor.test(data$SOC, data$MIMSOC)
r_val <- round(as.numeric(unlist(r_test['estimate'])),3)

# linear fit
mdl <- lm(SOC~MIMSOC)
mdl_smry <- summary(mdl)
fit <- round(mdl_smry$coefficients[2, 1],2)

# measured vs predicted SOC stock
MIMSOC_plot <- ggplot(data, aes(x=MIMSOC, y=SOC, color=CURVATURE)) + geom_point(size=3) +
  scale_color_gradientn(colours = rainbow(5)) +
  #xlim(2,8) + #ylim(0,10) +
  geom_abline(mapping=aes(slope=1, intercept=0), linetype="dotted") +
  annotate("text", label = paste0("r^2 = ",r_val), x = 6, y = 8, size = 3, colour = "red") +
  annotate("text", label = paste0("fit = ",fit,"y"), x = 6, y = 8.5, size = 3, colour = "red") +
  annotate("text", label = paste0("1:1"), x = 8, y = 8.4, size = 3, colour = "black")
MIMSOC_plot

#Microbial pool vs SOC (aiming for microbes = 1-2% of SOC)
MIC_plot <- ggplot(data, aes(x=MIMSOC, y=MIMMIC, color=MAP)) + geom_point() +
  geom_abline(mapping=aes(slope=0.015, intercept=0), linetype="dotted")

# Aiming for ... 
LIT_plot <- ggplot(data, aes(x=MIMSOC, y=MIMLIT, color=MAP)) + geom_point()

# Misc plots
CLAY_plot <- ggplot(data, aes(x=CLAY, y=MIMSOC, color=MAP)) + geom_point()
SAND_plot <- ggplot(data, aes(x=SAND, y=MIMSOC, color=MAP)) + geom_point()
ANPP_plot <- ggplot(data, aes(x=ANPP, y=MIMSOC, color=MAP)) + geom_point()
MAP_plot <- ggplot(data, aes(x=MAP, y=MIMSOC, color=MAP)) + geom_point() 
MAT_plot <- ggplot(data, aes(x=MAT, y=MIMSOC, color=MAT)) + geom_point()

# print plots
LIT_plot
MIC_plot
MIMSOC_plot




