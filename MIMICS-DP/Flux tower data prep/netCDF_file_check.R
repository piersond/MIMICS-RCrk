require(ncdf4)
require(ncdump)

nc <- nc_open("C:/Users/Derek/Google Drive/RCrk/RC Flux Tower data/RC-Rls/gap-filled_monthly_netCDF/2015-02.nc")
print(nc)
time <- ncvar_get(nc, "time")


nc2 <- nc_open("B:/Downloads/2016-02.nc")
print(nc2)
time2 <- ncvar_get(nc2, "time")



nc <- nc_open("C:/Users/Derek/Google Drive/RCrk/RC Flux Tower data/RC-Rls/gap-filled_monthly_netCDF/2016-04.nc")
nc_time <- ncvar_get(nc, "time")
nc_fsds <- ncvar_get(nc, "FSDS")
nc_df <- data.frame(time=nc_time, FSDS = nc_fsds)

nc_df$time <- nc_df$time+31

nc2 <- nc_open("C:/Users/Derek/Google Drive/RCrk/RC Flux Tower data/RC-Rls/gap-filled_monthly_netCDF/2016-03.nc")
nc2_time <- ncvar_get(nc2, "time")
nc2_fsds <- ncvar_get(nc2, "FSDS")
nc_df2 <- data.frame(time=nc2_time, FSDS = nc2_fsds)



df <- rbind(nc_df2, nc_df)
nc_test_plt <- ggplot(df, aes(x=time, y=FSDS)) + geom_point() + geom_line() + ggtitle("netCDF data: March-April 2016") 
nc_test_plt

nc_test_plt <- ggplot(df, aes(x=time, y=FSDS)) + geom_point() + geom_line() + ggtitle("netCDF data: March-April 2016") + xlim(28,34)
nc_test_plt



