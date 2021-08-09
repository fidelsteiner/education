################################################################################
# Retrieveal of AWS data for snow analysis
# 
# AWS_snow_trainingCode.R
#
# ReadMe: 
# 
# Code for the ICIMOD Training Snow data analysis in the Hindu Kush Himalaya region using R and Google Earth Engine
# Used to learn how to retrieve, clean and analyse station data in R
#
# 1) We use raw AWS data from Yala BC (ICIMOD station in Nepal), which is processed for all variables, plotted and then saved for further analysis
# 2) We check the specific time period where snow is interesting

# Created:          2017/02/05
# Latest Revision:  2021/07/17
#
# Jakob F Steiner| ICIMOD | x-hydrolab.org
################################################################################

# clear entire workspace (excl. packages)
rm(list = ls())
gc()

# define &-sign for pasting string-elements
'&' <- function(...) UseMethod('&')
'&.default' <- .Primitive('&')
'&.character' <- function(...) paste(...,sep='')

# install necessary packages if not available yet via install.packages()
require('lubridate')

# Other code necessary
source("D:\\Work\\Code\\BasicCode\\timeseriesAnalysis\\TSaggregate.R")    # Download from: https://github.com/fidelsteiner/BasicCode/blob/master/timeseriesAnalysis/TSaggregate.R

pathAWS <- 'D:\\Work\\Education\\ICIMODtrainings\\snowcourse2021\\Data'   # path where raw data files are located locally
pathOutput <- 'D:\\Work\\Education\\ICIMODtrainings\\snowcourse2021\\output' # Specify a folder for all output from your code


finOutput <- 'finalYalaBC2019.csv'                                        # name of our final output file

raw10minFile <- '2019_YalaBCAWS_10min_RAW.csv'                            # raw 10 min data file
raw60minFile <- '2019_YalaBCAWS_60min_RAW.csv'                            # raw 60 min data file

datRaw <- read.table(pathAWS&'\\'&raw10minFile,header = T, sep = ",", dec = ".")  # read the actual data from the file

timestr <- as.POSIXct(datRaw$TIMESTAMP, format="%m/%d/%Y  %H:%M")         # create a readable date string


# Plot raw data in original format, as an example donw with the battery status and air temperature

plot(datRaw$BattV)            # plot the whole time series
plot(datRaw$BattV[1:144])     # plot the first day

plot(datRaw$AirTC_Avg)            # plot the whole time series
plot(datRaw$AirTC_Avg[1:144])     # plot the first day

# Produce hourly values
BVOLhourly <- TSaggregate(datRaw$BattV,timestr,60,0,'mean')               # BVOL, Battery status [V]
BCONhourly <- TSaggregate(datRaw$Bucket_RT,timestr,60,0,'mean')           # BCON, Bucket content [mm]
PVOLhourly <- TSaggregate(datRaw$Accumulated_RT_NRT_Tot,timestr,60,0,'sum') # PVOL, instantaneous precipitation, mm

TAIRhourly <- TSaggregate(datRaw$AirTC_Avg,timestr,60,0,'mean')           # TAIR, mean air temperature, [degC]

# Check quality of RH
datRaw$RH_Avg[datRaw$RH_Avg>100] <- 100
datRaw$RH_Avg[datRaw$RH_Avg<=0] <- NA
RHhourly <- TSaggregate(datRaw$RH_Avg,timestr,60,0,'mean')                # RH, mean air rel humidity, [%]
TCNR4hourly <- TSaggregate(datRaw$CNR4TC_Avg,timestr,60,0,'mean')         # TCNR4, mean CNR4 temperature, sensor for radiation correction [degC]

# Check radiation data
datRaw$SUp_Avg[abs(datRaw$SUp_Avg)>2000] <- NA # to catch NA values saved as -7999
datRaw$SDn_Avg[abs(datRaw$SDn_Avg)>2000] <- NA # to catch NA values saved as -7999
datRaw$SUp_Avg[datRaw$SUp_Avg<0]<-0
datRaw$SDn_Avg[datRaw$SDn_Avg<0]<-0
datRaw$SDn_Avg[which(datRaw$SDn_Avg>datRaw$SUp_Avg)] <- datRaw$SUp_Avg[which(datRaw$SDn_Avg>datRaw$SUp_Avg)]

KINChourly <- TSaggregate(datRaw$SUp_Avg,timestr,60,0,'mean')             # KINC, incoming SW radiation, [W m-2]
KUPWhourly <- TSaggregate(datRaw$SDn_Avg,timestr,60,0,'mean')             # KUPW, outgoing SW radiation, [W m-2]
LINChourly <- TSaggregate(datRaw$LUpCo_Avg,timestr,60,0,'mean')           # LINC, incoming LW radiation, [W m-2]
LUPWhourly <- TSaggregate(datRaw$LDnCo_Avg,timestr,60,0,'mean')           # LOUT, outgoing LW radiation, [W m-2]

PREShourly <- TSaggregate(datRaw$BP_mbar,timestr,60,0,'mean')             # PRES, air pressure, [mbar]

# Wind

# function for averages of wind directions
winddirmean <- function(ws,wd,timestr,timStep,timShift){
  u =-ws*sin(wd*pi/180)
  v =-ws*cos(wd*pi/180)
  
  umean <- TSaggregate(u,timestr,timStep,timShift,'mean')
  vmean <- TSaggregate(v,timestr,timStep,timShift,'mean')
  
  windvec <- umean
  windvec[,2] <- windvec[,2]*NA
  
  windvec[umean[,2]>0,2] <- (90-180/pi*atan(umean[umean[,2]>0,2]/vmean[umean[,2]>0,2])+180)
  windvec[umean[,2]<=0,2] <- (90-180/pi*atan(umean[umean[,2]<=0,2]/vmean[umean[,2]<=0,2]))
  windvec[,3] <- windvec[,3]*NA
  return(windvec)
}

datRaw$WS_ms_Avg[datRaw$WS_ms_Avg<0] <- 0
datRaw$WS_ms_Avg[datRaw$WS_ms_Avg>50] <- NA
WSPDhourly <- TSaggregate(datRaw$WS_ms_Avg,timestr,60,0,'mean')           # WSPD, wind speed, [m s-1]
WSPDmaxhourly <- TSaggregate(datRaw$WS_ms_Avg,timestr,60,0,'max')         # WSPDmax, max wind speed, [m s-1]

WINDDIRhourly <- winddirmean(datRaw$WS_ms_Avg,datRaw$WindDir_D1_WVT, timestr,60,0) # WINDDIR, hourly wind direction, [deg]

# SR50 Data for snow depth
datRaw_SR50 <- read.table(pathAWS&'\\'&raw60minFile,header = T, sep = ",", dec = ".")
timestrSR50 <- as.POSIXct(datRaw_SR50$TIMESTAMP, format="%m/%d/%Y  %H:%M")

#Quality check (Q needs to be between 152 and 200)
datRaw_SR50$SR50_HAS_cor.1.[which(datRaw_SR50$SR50_Q.1.<152|datRaw_SR50$SR50_Q.1.>210)] <- NA
datRaw_SR50$SR50_HAS_cor.1.[which(datRaw_SR50$SR50_HAS_cor.1.==0)] <- NA
datRaw_SR50$SR50_HAS_cor.2.[which(datRaw_SR50$SR50_Q.2.<152|datRaw_SR50$SR50_Q.2.>210)] <- NA
datRaw_SR50$SR50_HAS_cor.2.[which(datRaw_SR50$SR50_HAS_cor.2.==0)] <- NA
datRaw_SR50$SR50_HAS_cor.3.[which(datRaw_SR50$SR50_Q.3.<152|datRaw_SR50$SR50_Q.3.>210)] <- NA
datRaw_SR50$SR50_HAS_cor.3.[which(datRaw_SR50$SR50_HAS_cor.3.==0)] <- NA
datRaw_SR50$SR50_HAS_cor.4.[which(datRaw_SR50$SR50_Q.4.<152|datRaw_SR50$SR50_Q.4.>210)] <- NA
datRaw_SR50$SR50_HAS_cor.4.[which(datRaw_SR50$SR50_HAS_cor.4.==0)] <- NA
datRaw_SR50$SR50_HAS_cor.5.[which(datRaw_SR50$SR50_Q.5.<152|datRaw_SR50$SR50_Q.5.>210)] <- NA
datRaw_SR50$SR50_HAS_cor.5.[which(datRaw_SR50$SR50_HAS_cor.5.==0)] <- NA
datRaw_SR50$SR50_HAS_cor.6.[which(datRaw_SR50$SR50_Q.6.<152|datRaw_SR50$SR50_Q.6.>210)] <- NA
datRaw_SR50$SR50_HAS_cor.6.[which(datRaw_SR50$SR50_HAS_cor.6.==0)] <- NA
datRaw_SR50$SR50_HAS_cor.7.[which(datRaw_SR50$SR50_Q.7.<152|datRaw_SR50$SR50_Q.7.>210)] <- NA
datRaw_SR50$SR50_HAS_cor.7.[which(datRaw_SR50$SR50_HAS_cor.7.==0)] <- NA
datRaw_SR50$SR50_HAS_cor.8.[which(datRaw_SR50$SR50_Q.8.<152|datRaw_SR50$SR50_Q.8.>210)] <- NA
datRaw_SR50$SR50_HAS_cor.8.[which(datRaw_SR50$SR50_HAS_cor.8.==0)] <- NA
datRaw_SR50$SR50_HAS_cor.9.[which(datRaw_SR50$SR50_Q.9.<152|datRaw_SR50$SR50_Q.9.>210)] <- NA
datRaw_SR50$SR50_HAS_cor.9.[which(datRaw_SR50$SR50_HAS_cor.9.==0)] <- NA
datRaw_SR50$SR50_HAS_cor.10.[which(datRaw_SR50$SR50_Q.10.<152|datRaw_SR50$SR50_Q.10.>210)] <- NA
datRaw_SR50$SR50_HAS_cor.10.[which(datRaw_SR50$SR50_HAS_cor.10.==0)] <- NA

allSR50 <- cbind(datRaw_SR50$SR50_HAS_cor.1.,datRaw_SR50$SR50_HAS_cor.2.,datRaw_SR50$SR50_HAS_cor.3.,datRaw_SR50$SR50_HAS_cor.4.,
                 datRaw_SR50$SR50_HAS_cor.5.,datRaw_SR50$SR50_HAS_cor.6.,datRaw_SR50$SR50_HAS_cor.7.,datRaw_SR50$SR50_HAS_cor.8.,
                 datRaw_SR50$SR50_HAS_cor.9.,datRaw_SR50$SR50_HAS_cor.10.)


SR50_Final <- rowMeans(allSR50,na.rm=T)
SR50_Final[SR50_Final>2.15] <- NA
SR50_Final[is.nan(SR50_Final)]<-NA

SR50hourly <- TSaggregate(SR50_Final,timestrSR50,60,0,'mean')

# Match the time steps of the hourly climate with the snow depth data
matchAWSdata <- match(TAIRhourly[,1],SR50hourly[,1])

# Combine all Data
Date <- as.Date(as.POSIXct(TAIRhourly[,1],origin='1970-01-01'))
Time <- strftime(as.POSIXct(TAIRhourly[,1],origin='1970-01-01'),format="%H:%M:%S")

expData <- data.frame(matrix(ncol = 0, nrow = length(Date)))
expData$DATE <- Date
expData$TIME <- Time
expData$BVOL <- BVOLhourly[,2]
expData$TAIR <- TAIRhourly[,2]
expData$RH <- RHhourly[,2]
expData$PVOL <- PVOLhourly[,2]
expData$BCON <- BCONhourly[,2]
expData$TCNR4 <- TCNR4hourly[,2]
expData$KINC <- KINChourly[,2]
expData$KUPW <- KUPWhourly[,2]
expData$LINC <- LINChourly[,2]
expData$LUPW <- LUPWhourly[,2]
expData$TSOIL <- TAIRhourly[,2] * NA
expData$LSD <- TAIRhourly[,2] * NA
expData$PRES <- PREShourly[,2]
expData$WSPD <- WSPDhourly[,2] 
expData$WSPDmax <- WSPDmaxhourly[,2]
expData$WINDDIR <- WINDDIRhourly[,2]
expData$SR50 <- SR50hourly[matchAWSdata,2]

write.csv(expData,file=pathOutput&'//'&finOutput, row.names=FALSE)

# Plot your climate data


png(file=pathOutput&'//AWS_YalaBC.png', res = 300,width=3600,height=1800)

par(mar=c(3,2,1,1),mai = c(0.5, 1, 0.1, 0.2),cex.lab=1.5,cex.axis=1.5)
layout(matrix(c(1,2,3), nrow = 3, ncol = 1, byrow = T))

datAx1 <- as.POSIXct(TAIRhourly[,1],origin='1970-01-01')


plot(datAx1,expData$TAIR,type='l',lwd=2,ylim = c(-30,10),xlab ='', ylab = expression('T'['air']~~'[°C]'),xaxt='n')
abline(h = seq(-30,40,10),v = seq(datAx1[1],datAx1[length(datAx1)], by = "month"), col="gray", lty=3)

plot(datAx1,expData$KINC - expData$KUPW,type='l',lwd=2,ylim = c(0,1000),xlab ='', ylab = expression('SW'['net']~~'[W '~ m^{-2}~ ']'),xaxt='n')
abline(h = seq(0,1000,100),v = seq(datAx1[1],datAx1[length(datAx1)], by = "month"), col="gray", lty=3)


plot(datAx1,expData$BCON,type='l',xlab ='',ylab = expression('accum precipitation [mm]'),ylim=c(0,1500),xaxt='n')
abline(h = seq(0,1500,250),v = seq(datAx1[1],datAx1[length(datAx1)], by = "month"), col="gray", lty=3)


axis.POSIXct(1, at = seq(datAx1[1],datAx1[length(datAx1)], by = "month"), format = "%b-%d")


dev.off()

# Plot your snow depth data; the sensor is installed at 2.15 m height

plot(expData$SR50)                              # the straightfoward SR50 data

plot(2.15 - expData$SR50)                       # actual "depth"

actualSD <- 2.15 - expData$SR50


par(mar=c(3,2,1,1),mai = c(0.5, 1, 0.1, 0.2),cex.lab=1.5,cex.axis=1.5)
layout(matrix(c(1,2,3), nrow = 3, ncol = 1, byrow = T))

plot(datAx1,expData$BCON,type='l',xlab ='',ylab = expression('accum precipitation [mm]'),ylim=c(0,1500),xaxt='n')
abline(h = seq(0,1500,250),v = seq(datAx1[1],datAx1[length(datAx1)], by = "month"), col="gray", lty=3)

plot(datAx1,expData$TAIR,type='l',lwd=2,ylim = c(-30,10),xlab ='', ylab = expression('T'['air']~~'[°C]'),xaxt='n')
abline(h = seq(-30,40,10),v = seq(datAx1[1],datAx1[length(datAx1)], by = "month"), col="gray", lty=3)

plot(datAx1,actualSD,ylim = c(0,0.8),xlab ='', ylab = expression('snow depth [m]'),xaxt='n')
abline(h = seq(0,1,0.2),v = seq(datAx1[1],datAx1[length(datAx1)], by = "month"), col="gray", lty=3)
axis.POSIXct(1, at = seq(datAx1[1],datAx1[length(datAx1)], by = "month"), format = "%b-%d")

# remove impossible snow values due to season for example

mon_NoSnow <- c(6,7,8,9)
MonthID <- month(datAx1)
actualSD[which(!is.na(match(MonthID,mon_NoSnow)))] <- 0

# Focus only on period where it is interesting

which(actualSD>0)         # check the IDs where snow occurs

focustime <- seq(4300,7600,1)

png(file=pathOutput&'//AWS_YalaBC_SNOW.png', res = 300,width=3600,height=1800)


par(mar=c(3,2,1,1),mai = c(0.5, 1, 0.1, 0.2),cex.lab=1.5,cex.axis=1.5)
layout(matrix(c(1,2,3,4,5), nrow = 5, ncol = 1, byrow = T))

plot(datAx1[focustime],expData$RH[focustime],type='l',xlab ='',ylab = expression('RH [%]'),ylim=c(0,100),xaxt='n')
abline(h = seq(0,100,20),v = seq(datAx1[1],datAx1[length(datAx1)], by = "month"), col="gray", lty=3)

plot(datAx1[focustime],expData$TAIR[focustime],type='l',lwd=2,ylim = c(-30,10),xlab ='', ylab = expression('T'['air']~~'[°C]'),xaxt='n')
abline(h = seq(-30,40,10),v = seq(datAx1[1],datAx1[length(datAx1)], by = "month"), col="gray", lty=3)

netRad <- expData$KINC - expData$KUPW + expData$LINC - expData$LUPW
plot(datAx1[focustime],netRad[focustime],type='l',lwd=2,ylim = c(-100,800),xlab ='', ylab = expression('R'['net']~~'[W '~ m^{-2}~ ']'),xaxt='n')
abline(h = seq(-100,800,100),v = seq(datAx1[1],datAx1[length(datAx1)], by = "month"), col="gray", lty=3)

plot(datAx1[focustime],expData$WSPD[focustime],type='l',lwd=2,ylim = c(0,15),xlab ='', ylab = expression('u [m/s]'),xaxt='n')
abline(h = seq(0,15,5),v = seq(datAx1[1],datAx1[length(datAx1)], by = "month"), col="gray", lty=3)

plot(datAx1[focustime],actualSD[focustime],ylim = c(0,0.8),xlab ='', ylab = expression('snow depth [m]'),xaxt='n')
abline(h = seq(0,1,0.2),v = seq(datAx1[1],datAx1[length(datAx1)], by = "month"), col="gray", lty=3)
axis.POSIXct(1, at = seq(datAx1[1],datAx1[length(datAx1)], by = "month"), format = "%b-%d")

dev.off()

# Produce daily values

SD_daily <- TSaggregate(actualSD,datAx1,24,0,'mean')
SD_daily_max <- TSaggregate(actualSD,datAx1,24,0,'max')[,2]
SD_daily_min <- TSaggregate(actualSD,datAx1,24,0,'min')[,2]