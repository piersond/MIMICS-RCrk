### MIMICS brute force best parameters

#Libraries
library(rootSolve)
library(boot)
library(ggplot2)
library(tidyverse)

#bring in RXEQ function
source("C:/github/MIMICS-RCrk/MIMICS-DP/MIMICS_brute_force/RXEQ_ftn.R")

########################################
# Set MIMICS parameters
########################################
Vslope  <- rep(0.063, 6)
Vint    <- rep(5.47, 6)
aV      <- rep(0.000008, 6)  
Kslope  <- rep(c(0.025, 0.035, 0.025),2)
Kint    <- rep(3.19, 6)
aK      <- rep(10, 6)
vMOD    <- c(10, 2, 10, 3, 3, 2)
kMOD    <- c(8, 2, 4, 2, 4, 6)
KO      <- c(6, 6)
CUE     <- c(0.55, 0.25, 0.75, 0.35)
tau_r   <- c(0.00052, 0.3)
tau_K   <- c(0.00024, 0.1)
Tau_MOD <- c(100, 0.8, 1.2, 2)
fPHYS_r <- c(0.3, 1.3)
fPHYS_K <- c(0.2, 0.8)
fCHEM_r <- c(0.1, -3, 1)
fCHEM_K <- c(0.3, -3, 1)
fSOM_p  <- c(0.000015, -1.5)
PHYS_scalar <- c(2, -2)
FI      <- c(0.05, 0.05)
fmet_p <- c(1, 0.85, 0.013)
depth <- 30

########################################
# Set litter parameters
########################################
litCN <- c(133.3, #TRAEf
           92.7,  #PIREf
           83.1,  #THPLf
           61.8,  #ACSAf
           50.5,  #QUPRf
           24.2)  #DRGLf

litLIG <- c(16.2, #TRAEf
            19.2, #PIREf
            26.7, #THPLf
            15.9, #ACSAf
            23.5, #QUPRf
            10.9) #DRGLf

calcN    <- (1 / litCN) / 2.5 * 100    
lit_fMET <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * litLIG / calcN) 

########################################
# Set forcing data variables
########################################
data <- read.csv(paste0(wd,"Site_data/MIMICS_pts_pMAP_fW_112320.csv")) #site level forcing variables

ANPP   <- data$pGPP+400
clay   <- data$CLAY/100  				    
tsoi   <- data$MAT
nsites <- length(data$Site)

lig    <- data$LIG #/ 100
Nnew   <- data$N                                                 
fMET1  <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * lig / Nnew)   

#Constants pulled out from MIMICS loop
fMET       <- mean(lit_fMET)  
h2y        <- 24*365
MICROtoECO <- depth * 1e4 * 1e-3         # mgC/cm3 to g/m2

########################################
# Apply rate curve brute force multipliers
########################################
# Vslope = Vslope * Vslope_mult
# Vint = Vint * Vint_mult
# Kslope = Kslope * Kslope_mult
# Kint = Kint * Kint_mult


###########################################
# MIMICS single point function
###########################################
MIMICS1 <- function(df){
  print(paste("-------- starting ", df$Site, " --------") )
  ANPP       <- df$pGPP+400
  fCLAY      <- df$CLAY/100
  TSOI       <- df$MAT
  FW         <- df$fW                #<-- DPierson added fW
  EST_LIT <- (ANPP / (365*24)) * 1e3 / 1e4         # gC/m2/h (from gC/m2/y) then mgC/cm2/h(from gC/m2/h) 

  # ------------ caclulate parameters ---------------
  Vmax     <- exp(TSOI * Vslope + Vint) * aV 
  Km       <- exp(TSOI * Kslope + Kint) * aK
  
  #ANPP strongly correlated with MAP
  Tau_MOD1 <- sqrt(df$ANPP/Tau_MOD[1])          # basicaily standardize against NWT
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
  desorb   <- fSOM_p[1] * exp(fSOM_p[2]*(fCLAY))                 #CHANGED FOR GLOBAL RUN!!!     
  
  #desorb   <- desorb/10 # modified as in MIMdef from Zhang et al 2020
  #fPHYS    <- fPHYS/5  # to reduce allocation to physically protected pool 5x
  
  pSCALAR  <- PHYS_scalar[1] * exp(PHYS_scalar[2]*(sqrt(fCLAY)))  #Scalar for texture effects on SOMp
  v_MOD    <- vMOD  # to avoid writing over orig. parameters
  k_MOD    <- kMOD 
  k_MOD[3] <- k_MOD[3] * pSCALAR    
  k_MOD[6] <- k_MOD[6] * pSCALAR    
  
  VMAX     <- Vmax * v_MOD #* FW    #<-- DPierson added fW
  KM       <- Km / k_MOD
  
  ### Apply tau brute force multiplier
  #tau <- tau * Tau_mult
  
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
  MICtrn  <- c(NA,NA,NA,NA,NA,NA) #rep(NA, dim=6)
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
              SOM_1 = som[1], SOM_2 = som[2], SOM_3 = som[3])

## Set global parameters to allow pass to stode function
  .GlobalEnv$VMAX <- VMAX
  .GlobalEnv$KM <- KM
  .GlobalEnv$fPHYS <- fPHYS
  .GlobalEnv$fCHEM <- fCHEM
  .GlobalEnv$fAVAI <- fAVAI
  .GlobalEnv$I <- I
  #.GlobalEnv$FI <- FI
  #.GlobalEnv$CUE <- CUE
  .GlobalEnv$tau <- tau
  .GlobalEnv$LITmin <- LITmin
  .GlobalEnv$SOMmin <- SOMmin
  .GlobalEnv$MICtrn <- MICtrn
  .GlobalEnv$desorb <- desorb
  .GlobalEnv$DEsorb <- DEsorb
  .GlobalEnv$OXIDAT <- OXIDAT
  

  test  <- stode(y = Ty, time = 1e6, fun = RXEQ, parms = Tpars, positive = TRUE)
  
###Return MIMICS output 
  #table[i,2:8] <- as.numeric(test[[1]])
  MIMLIT    <- (test[[1]][[1]]+test[[1]][[2]])  * depth *1e4 / 1e6 #convert kgC/m2 from mgC/cm3 (0-30 cm) 
  MIMMIC    <- (test[[1]][[3]]+test[[1]][[4]])  * depth *1e4 / 1e6
  MIM_CO    <-  test[[1]][[3]]/test[[1]][[4]]
  MIMSOC    <- sum(test[[1]])  * depth *1e4 / 1e6   
  
  MIMout <- data.frame(Site = df$Site,
                       fCLAY = fCLAY,
                       TSOI = TSOI,
                       ANPP = ANPP,
                       FW = FW,                
                       EST_LIT = EST_LIT,
                       SOC = df$SOC,
                       MIMSOC = MIMSOC,
                       MIMMIC = MIMMIC,
                       MIM_CO = MIM_CO
                       )
  
  #remove global variables set for stode ftn
  rm(I, VMAX, KM, fPHYS, fCHEM, fAVAI, tau, LITmin, SOMmin, MICtrn, desorb, DEsorb, OXIDAT)
  
  return(MIMout)
}

###############################
# Using the MIMICS1 ftn
###############################

#single point run
MIMout_single <- MIMICS1(data[1,])

#full run of forcing data csv
MIMrun <- data %>% split(1:nrow(data)) %>% map(MIMICS1) %>% bind_rows()

# test plot
ggplot(MIMrun, aes(x=MIMSOC, y=SOC)) + geom_point()



