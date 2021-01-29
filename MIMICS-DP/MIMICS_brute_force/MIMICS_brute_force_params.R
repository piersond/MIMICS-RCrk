### MIMICS brute force best parameters

#Libraries
library(rootSolve)
library(boot)
library(ggplot2)

#working drive address
wd <- "C:/github/MIMICS-RCrk/MIMICS-DP/"

#bring in RXEQ function
source(paste0(wd,"MIMICS_brute_force/RXEQ_ftn.R"))


# Create function to allow for MIMICS run loop
brute_MIMICS <- function(Tau_mult, Vslope_mult, Vint_mult, Kslope_mult, Kint_mult) {

  ########################################
  # Set MIMICS parameters
  ########################################
  Vslope  <- rep(0.063, 6)
  Vint    <- rep(5.47, 6)
  Av      <- rep(0.000008, 6)  
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
  MIMLIT <- rep(NA, nsites)           	                      
  MIMMIC <- rep(NA, nsites)           
  MIM_CO <- rep(NA, nsites)           
  MIMSOC <- rep(NA, nsites) 
  
  strSite <- as.character(data$Site)                          
  pools  <- c('site','LITm', 'LITs', 'MICr', 'MICK', 'SOMp', 'SOMc', 'SOMa')
  table  <- array(NA, dim=c(nsites,8), dimnames=list(as.character(data$Site),pools)) 
  table[,1]       <- as.character(data$Site)
  
  POOLS  <- c("LIT","MIC","SOC")
  npools <- length(POOLS)  
  LITpool<- c('LITm', 'LITs') 
  MICpool<- c('MICr', 'MICk') 
  SOMpool<- c('SOMp', 'SOMa','SOMc') 
  
  MIMLIT <- array(NA, dim=c(nsites))           	            
  MIMMIC <- array(NA, dim=c(nsites))           
  MIM_CO <- array(NA, dim=c(nsites))           
  MIMSOC <- array(NA, dim=c(nsites))           
  
  strSite <- as.character(data$Site) 
  
  #Make vectors to store model results that match obs temporal resolution
  npts   <- 6*10*14   				#6 litter * 10 years * 14 sites
  xyLIT  <- rep(NA, npts) 
  xyTIME <- rep(NA, npts) 
  xySITE <- rep(NA, npts) 
  xyOBS  <- rep(NA, npts)
  xyMIM  <- rep(NA, npts)
  xyCount<- 1
  
  ########################################
  # MIMICS LOOP
  ########################################

  for (i in 1:nsites) {         #speeds up debugging     	
    
  ### Debug single point run
    # i <- 1
    # Vslope_mult <- 1
    # Vint_mult <- 1
    # Kslope_mult <- 1
    # Kint_mult <- 1
    # Tau_mult <- 1
  ### 
     
    # Read in site characteristics -----------
    print(paste("-------- starting ", data$Site[i], " --------") )
    
    fMET       <- mean(lit_fMET)           # uses mean litter fmet from LIDET
    fCLAY      <- clay[i]
    TSOI       <- tsoi[i]
    FW         <- fW[i]                #<-- DPierson added fW
    EST_LIT_in <- (ANPP[i] / (365*24))        # gC/m2/h (from gC/m2/y)
    h2y        <- 24*365
    MICROtoECO <- depth * 1e4 * 1e-3         # mgC/cm3 to g/m2
    EST_LIT    <- EST_LIT_in  * 1e3 / 1e4    # mgC/cm2/h(from gC/m2/h) 
    
    # ------------ caclulate parameters ---------------
    
    # Apply rate curve brute force multipliers
    Vslope = Vslope * Vslope_mult
    Vint = Vint * Vint_mult
    Kslope = Kslope * Kslope_mult
    Kint = Kint * Kint_mult
  
    Vmax     <- exp(TSOI * Vslope + Vint) * aV 
    Km       <- exp(TSOI * Kslope + Kint) * aK
    
    #ANPP strongly correlated with MAP
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
    tau <- tau * Tau_mult
  
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
    

    #create RXEQ function for use in loop
    source(paste0(wd,"MIMICS_brute_force/RXEQ_ftn.R"))
    
    #Calculate RXEQ pools  
    Tpars <- c( I = I, VMAX = VMAX, KM = KM, CUE = CUE, 
                fPHYS = fPHYS, fCHEM = fCHEM, fAVAI = fAVAI, FI = FI, 
                tau = tau, LITmin = LITmin, SOMmin = SOMmin, MICtrn = MICtrn, 
                desorb = desorb, DEsorb = DEsorb, OXIDAT = OXIDAT, KO = KO)
    Ty    <- c( LIT_1 = lit[1], LIT_2 = lit[2], 
                MIC_1 = mic[1], MIC_2 = mic[2], 
                SOM_1 = som[1], SOM_2 = som[2], SOM_3 = som[3] )
    
    print("Tpars")
    print(Tpars)
    
    
    test  <- stode(y = Ty, time = 1e6, fun = RXEQ, parms = Tpars, positive = TRUE)
    if (i == 2) {print(test[1])}
    remove(lit, mic, som)
    
    table[i,2:8] <- as.numeric(test[[1]])
    MIMLIT[i]    <- (test[[1]][[1]]+test[[1]][[2]])  * depth *1e4 / 1e6#convert kgC/m2 from mgC/cm3 (0-30 cm) 
    MIMMIC[i]    <- (test[[1]][[3]]+test[[1]][[4]])  * depth *1e4 / 1e6
    MIM_CO[i]    <-  test[[1]][[3]]/test[[1]][[4]]
    MIMSOC[i]    <- sum(test[[1]])  * depth *1e4 / 1e6   
    
    remove(test, Ty, Tpars, LITmin, MICtrn, SOMmin)
    
  } # close i loop (sites)
  
  
  ########################################
  # MIMICS RUN DATA TO RETURN
  ########################################
  
  # Calculate correlation between field SOC and MIMSOC
  r2_test <- cor.test(data$SOC, MIMSOC)
  r_val <- round(as.numeric(unlist(r2_test ['estimate'])),3)
  
  # Calculate MIC pool size relative to SOC
  MICpoolRelSOC_mn <- mean(MIMMIC/MIMSOC)
  
  # Calculate LIT pool size relative to SOC
  LITpoolRelSOC_mn <- mean(MIMLIT/MIMSOC)
  
  # Calculate average residual value
  resid_avg <- mean(data$SOC - MIMSOC)
  resid_sd <- sd(data$SOC - MIMSOC)
  
  MIMOUT <- data.frame(
    r2 = r_val,
    resid_avg = resid_avg,
    resid_sd = resid_sd,
    MICpoolRelSOC_mn = MICpoolRelSOC_mn,
    LITpoolRelSOC_mn = LITpoolRelSOC_mn,
    Tau_mult = Tau_mult,
    Vslope_mult = Vslope_mult,
    Vint_mult = Vint_mult,
    Kslope_mult = Kslope_mult,
    Kint_mult = Kint_mult,
    Vslope = Vslope,
    Vint = Vint,
    Vmax = Vmax,
    Km = Km,
    CUE=CUE,
    tau_rc = tau_rc
  )
    
  return(MIMOUT)

}


### BRUTE FORCE MIMICS ###

##############
# Set brute force loop parameters
num_runs <- 100
Tau_mult_in <- runif(num_runs, 0.1, 2) 
Vslope_mult_in <- runif(num_runs, 1, 1) 
Vint_mult_in <- runif(num_runs, 1, 1) 
Kslope_mult_in <- runif(num_runs, 1, 1) 
Kint_mult_in <- runif(num_runs, 1, 1) 

test <- brute_MIMICS(Tau_mult = Tau_mult_in,
                      Vslope_mult = Vslope_mult_in,
                      Vint_mult = Vint_mult_in,
                      Kslope_mult = Kslope_mult_in,
                      Kint_mult = Kint_mult_in)
