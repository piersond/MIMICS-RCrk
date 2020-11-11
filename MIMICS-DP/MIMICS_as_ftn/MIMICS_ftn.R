library(rootSolve)
library(boot)

#run_MIMICS <- function(data_in, litter_in, params) {
  
  #debug
  setwd("C:/GitHub/MIMICS-RCrk/MIMICS-DP/")
  data_in <- read.csv("Site_data/LTER_SITE_1.csv", as.is=T)
  litter_in <- litter <- read.csv("Litter_characteristics/LIDET_LitterCharacteristics.csv")
  litter_in <- litter_in[1:6,] 
  params <- read.csv("MIMICS_parameters/parameters_LIDET-MIM-REV.csv")#parm_ball[[20]]
  
  
  #set parameter values
  Vslope <- params$Vslope
  Vint<- params$Vint
  aV <- params$aV
  Kslope <- params$Kslope
  Kint <- params$Kint
  aK <- params$aK
  vMOD <- params$vMOD
  kMOD <- params$kMOD
  KO <- params$KO[1:2]
  CUE <- params$CUE[1:4]
  tau_r <- params$tau_r[1:2]
  tau_K <- params$tau_K[1:2]
  Tau_MOD <- params$Tau_MOD[1:4]

  fPHYS_r <- params$fPHYS_r[1:2]
  fPHYS_K <- params$fPHYS_K[1:2]
  fCHEM_r <- params$fCHEM_r[1:3]
  fCHEM_K <- params$fCHEM_K[1:3]
  fSOM_p <- params$fSOM_p[1:2]
  PHYS_scalar <- params$PHYS_scalar[1:2]
  FI <- params$FI[1:2]
  fmet_p <- params$fmet_p[1:3]
  depth <- params$depth[1]

  calcN    <- (1 / litter_in$litCN) / 2.5 * 100  #DP: what's this calc N for?    
  lit_fMET <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * litter_in$litLIG / calcN) 
  
  ANPP  <- data_in$ANPP / 2           		# if needed convert to gC/m2/y from g/m2/y
  clay  <- data_in$CLAY2/100  				    # if needed, convert from % clay to fraction
  tsoi  <- data_in$MAT
  nsites<- length(data_in$Site)
  
  lig    <- data_in$LIG #/ 100
  Nnew   <- data_in$N                                                  #N in litter additions
  fMET1  <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * lig / Nnew)   #as partitioned in Daycent
  MIMLIT <- rep(NA, nsites)           	                       #Vector for results
  MIMMIC <- rep(NA, nsites)           
  MIM_CO <- rep(NA, nsites)           
  MIMSOC <- rep(NA, nsites)           
  
  strSite <- as.character(data_in$Site)                           #convert site names to string
  pools  <- c('site','LITm', 'LITs', 'MICr', 'MICK', 'SOMp', 'SOMc', 'SOMa')
  table  <- array(NA, dim=c(nsites,8), dimnames=list(as.character(data_in$Site),pools)) 
  table[,1] <- as.character(data_in$Site)
  
  
  POOLS  <- c("LIT","MIC","SOC")
  npools <- length(POOLS)  
  LITpool<- c('LITm', 'LITs') 
  MICpool<- c('MICr', 'MICk') 
  SOMpool<- c('SOMp', 'SOMa','SOMc') 
  
  MIMLIT <- array(NA, dim=c(nsites))           	             # array for results
  MIMMIC <- array(NA, dim=c(nsites))           
  MIM_CO <- array(NA, dim=c(nsites))           
  MIMSOC <- array(NA, dim=c(nsites))           
  
  strSite <- as.character(data_in$Site)  #convert site names to string
  
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
    
    #DEBUG
    i <- 1
    
    
    # Read in site characteristics -----------
    print(paste("-------- starting ", data_in$Site[i], " --------") )
    
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


    ### STODE ftn ######
    #test <- with (as.list(c(Ty, Tpars)),{
      
      #Flows to and from MIC_1
      LITmin[1] = Ty['MIC_1'] * Tpars['VMAX1'] * Ty['LIT_1'] / (Tpars['KM1'] + Ty['MIC_1'])   #MIC_1 decomp of MET lit
      LITmin[2] = Ty['MIC_1'] * Tpars['VMAX2'] * Ty['LIT_2'] / (Tpars['KM2'] + Ty['MIC_1'])   #MIC_1 decomp of STRUC lit
      MICtrn[1] = Ty['MIC_1'] * Tpars['tau1']  * Tpars['fPHYS1']                  #MIC_1 turnover to PHYSICAL SOM 
      MICtrn[2] = Ty['MIC_1'] * Tpars['tau1']  * Tpars['fCHEM1']                  #MIC_1 turnover to CHEMICAL SOM  
      MICtrn[3] = Ty['MIC_1'] * Tpars['tau1']  * Tpars['fAVAI1']                 #MIC_1 turnover to AVAILABLE SOM  
      SOMmin[1] = Ty['MIC_1'] * Tpars['VMAX3'] * Ty['SOM_3'] / (Tpars['KM3'] + Ty['MIC_1'])   #decomp of SOMa by MIC_1
      
      #Flows to and from MIC_2
      LITmin[3] = Ty['MIC_2'] * Tpars['VMAX4'] * Ty['LIT_1'] / (Tpars['KM4'] + Ty['MIC_2'])   #decomp of MET litter
      LITmin[4] = Ty['MIC_2'] * Tpars['VMAX5'] * Ty['LIT_2'] / (Tpars['KM5'] + Ty['MIC_2'])   #decomp of SRUCTURAL litter
      MICtrn[4] = Ty['MIC_2'] * Tpars['tau2']  * Tpars['fPHYS2']                  #MIC_2 turnover to PHYSICAL  SOM 
      MICtrn[5] = Ty['MIC_2'] * Tpars['tau2']  * Tpars['fCHEM2']                  #MIC_2 turnover to CHEMICAL  SOM  
      MICtrn[6] = Ty['MIC_2'] * Tpars['tau2']  * Tpars['fAVAI2']                 #MIC_2 turnover to AVAILABLE SOM  
      SOMmin[2] = Ty['MIC_2'] * Tpars['VMAX6'] * Ty['SOM_3'] / (Tpars['KM6'] + Ty['MIC_2'])   #decomp of SOMa by MIC_2
      
      DEsorb    = Ty['SOM_1'] * Tpars['desorb']  #* (MIC_1 + MIC_2)  	#desorbtion of PHYS to AVAIL (function of fCLAY)
      OXIDAT    = ((Ty['MIC_1'] * Tpars['VMAX2'] * Ty['SOM_2'] / (Tpars['KO1']*Tpars['KM2'] + Ty['MIC_1'])) +
                     (Ty['MIC_2'] * Tpars['VMAX5'] * Ty['SOM_2'] / (Tpars['KO2']*Tpars['KM5'] + Ty['MIC_2'])) )  #oxidation of C to A
      
      dLIT_1 = Tpars['I1']*(1-Tpars['FI1']) - LITmin[1] - LITmin[3]
      dMIC_1 = Tpars['CUE1']*(LITmin[1]+ SOMmin[1]) + Tpars['CUE2']*(LITmin[2]) - sum(MICtrn[1:3])
      dSOM_1 = Tpars['I1']*Tpars['FI1'] + MICtrn[1] + MICtrn[4]- DEsorb 
      
      dLIT_2 = Tpars['I2'] * (1-Tpars['FI2']) - LITmin[2] - LITmin[4]
      dMIC_2 = Tpars['CUE3']*(LITmin[3]+ SOMmin[2]) + Tpars['CUE4']*(LITmin[4]) - sum(MICtrn[4:6])  
      dSOM_2 = Tpars['I2']*Tpars['FI2'] + MICtrn[2] + MICtrn[5] - OXIDAT
      
      dSOM_3 = MICtrn[3] + MICtrn[6] + DEsorb + OXIDAT - SOMmin[1] - SOMmin[2]
      
      test <- list(c(dLIT_1, dLIT_2, dMIC_1, dMIC_2, dSOM_1, dSOM_2, dSOM_3))
    
    
    print(test[[1]])
    
    #if (i == 2) {print(test[1])}
    remove(lit, mic, som)
    
    table[i,2:8] <- as.numeric(test[[1]])
    MIMLIT[i]    <- (test[[1]][[1]]+test[[1]][[2]])  * depth *1e4 / 1e6#convert kgC/m2 from mgC/cm3 (0-30 cm) 
    MIMMIC[i]    <- (test[[1]][[3]]+test[[1]][[4]])  * depth *1e4 / 1e6
    MIM_CO[i]    <-  test[[1]][[3]]/test[[1]][[4]]
    MIMSOC[i]    <- sum(test[[1]])  * depth *1e4 / 1e6   
    
    #remove(test, Ty, Tpars, LITmin, MICtrn, SOMmin)
  } 
  
  mean(MIMSOC) #3.24184e-05
  Ty
  Tpars
  # correlation
  r_test <- NA
  r_test <- cor.test(data_in$SOC, MIMSOC)
  r_val <- round(as.numeric(unlist(r_test['estimate'])),3)
  
  # linear fit
  mdl <- lm(data_in$SOC~MIMSOC)
  mdl_smry <- summary(mdl)
  fit <- round(mdl_smry$coefficients[2, 1],2)
  
  output <- list(rep = params$rep,
                 tau = tau[1], 
                 SOC=data_in$SOC[1], 
                 MIMSOC=MIMSOC[12])
  return(output)
  
  #return(c(params$rep, params$tau_r[1], r_val, paste0(as.character(fit),"y")))
} 