library(tidyverse)
library(ggplot2)

# Calculate mean annual soil mositure (MASM) from ARS bimonthly data file
calc_MASM <- function(ars_sm_filename) {
  soil_mst <- read.csv(ars_sm_filename, as.is=T)
  soil_mst$datetime <- as.POSIXct(soil_mst$datetime)
  loc <- str_split(ars_sm_filename,"-")[[1]][3]
  
  #Remove rows with missing data
  soil_mst[soil_mst == -999] <- NA
  
  #Create columns for month and year
  soil_mst$Month <- as.numeric(format(soil_mst$datetime,'%m'))
  soil_mst$Year <- as.numeric(format(soil_mst$datetime,'%Y'))
  
  #set numeric columns to remove any non-numerics
  suppressWarnings(
    soil_mst[2:ncol(soil_mst)] <- lapply(soil_mst[2:ncol(soil_mst)], as.numeric)
  )
  
  #Create column for water year (October-Sept)
    #Also try by growing season (March-July)
  soil_mst$Year_water <- ifelse(soil_mst$Month > 9, soil_mst$Year+1, soil_mst$Year)
  
  #Create monthly average dataset
  soil_mst_monthly <- soil_mst %>% group_by(Year_water, Month, probe) %>% summarise_all(mean) 
  
  #Use only probe 564
  mst_monthly_564 <- soil_mst_monthly %>% filter(probe == 564) 
  
  #Use monthly average to create annual mean, max, min
    #Doesn't calc for years missing data
  mst_annual_mean <- mst_monthly_564 %>% group_by(Year_water) %>% summarize_all(mean)
  #mst_annual_max <- mst_monthly_564 %>% group_by(Year_water) %>% summarize_all(max)
  #mst_annual_min <- mst_monthly_564 %>% group_by(Year_water) %>% summarize_all(min)
  
  #store number of years in annual dataset
  years_n <- nrow(na.omit(mst_annual_mean))
  
  #mean, max, min moisture values
  mean_annual_mst <- mst_annual_mean %>% select(-probe, -Month, -datetime, -Year) %>% na.omit() %>% summarize_all(mean) %>% 
                      select(-Year_water) %>% add_column(met = loc, .before = "mst015") %>% add_column(yrs_n = years_n , .after = "met")
  
  return(mean_annual_mst)
}


#set path to soil moisture data
soil_water_data_path <- "C:/Users/drkpi/Google Drive/RCrk/RC_SoDaH/RC_data/ARS_climate/neutron probe soil water/headers_removed/"

#get ARS siol moisture data files from path
soil_mst_files <- list.files(path=soil_water_data_path)
soil_mst_files <- soil_mst_files[!soil_mst_files %in% c('desktop.ini')]
path_soil_mst_files <- paste0(soil_water_data_path, soil_mst_files)

#Compile from all ARS files
comp_df <- sapply(path_soil_mst_files,calc_MASM)
MASM <- do.call(bind_rows, comp_df)

#write MASM
write.csv(MASM, "C:/Users/drkpi/Google Drive/RCrk/RC_SoDaH/RC_data/ARS_climate/neutron probe soil water/MASM_pierson_fromR.csv")


# #plots
# ggplot(soil_mst, aes(x=datetime, y=mst015, color=as.character(probe))) + geom_point(alpha=0.5) + ylim(0,0.4) + geom_line()
# 
# ggplot(soil_mst_monthly, aes(x=datetime, y=mst015, color=as.character(probe))) + geom_point(alpha=0.5) + ylim(0,0.4) + geom_line()
# 
# ggplot(mst_monthly_564, aes(x=datetime, y=mst015, color=as.character(Year_water))) + geom_point(alpha=0.5) + ylim(0,0.4) + geom_line()
# 
# ggplot(mst_annual_mean, aes(x=Year_water, y=mst015, color=as.character(Year_water))) + geom_point(alpha=0.5) + ylim(0,0.4) + geom_line()
# 
