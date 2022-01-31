################################################################################
# Running a simple energy balance model for snow melt simulations at the point scale
# 
# EB_snow_trainingCode.R
#
# ReadMe: 
# 
# Code for the ICIMOD Training Snow data analysis in the Hindu Kush Himalaya region using R and Google Earth Engine
# Used to learn how to build an energy balance model
#


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

raw10minFile <- '2019_YalaBCAWS_10min_RAW.csv'                            # raw 10 min data file
raw60minFile <- '2019_YalaBCAWS_60min_RAW.csv'                            # raw 60 min data file

datRaw <- read.table(pathAWS&'\\'&raw10minFile,header = T, sep = ",", dec = ".")  # read the actual data from the file

timestr <- as.POSIXct(datRaw$TIMESTAMP, format="%m/%d/%Y  %H:%M")         # create a readable date string

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

WINDDIRhourly <- winddirmean(datRaw$WS_ms_Avg,datRaw$WindDir_D1_WVT, timestr,60,5) # WINDDIR, hourly wind direction, [deg]

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

##### Prepare Parameters bulk method-- 
g = 9.81            # Gravity constant (ms-2)
k = 0.41            # Vonkarmann constant (-)
LL1 = 2830000       # Latent heat of sublimation (J kg-1)
LL2 = 2500000       # Latent heat of evaporation (J kg-1)
Zv = 2.5            # Height of wind speed measurements (m)
Zt = 1.75           # Height of temperature measurements (m)
Zrh = 1.75          # Height of relative humidity measurements (m)
Zq =1.75            # Height of relative humidity measurements (m)
z0 =0.012           # Momentum roughness length (m)
z0t =0.0012         # Thermal roughness length (m)
z0q =0.0012         # Humidity roughness length (m)

sigma <- 5.67 * 10^-8     # Stefan Boltzmann constant
#emiss_surf <- 0.8         # surface emissivity; changes with surface but for the simple purpose here we assume to be constant

#### Produce relevant time series for the energy balance

expData$TSURF <- (expData$LUPW / sigma )^(1/4) - 273.15 + 25


LL <- vector()      # initialize vector for latent heat (J kg-1)

LL[expData$TSURF>=0.0] <- LL2                   # If surface temperature >= 0.0, then Latent heat flux equals latent heat of evaporation
LL[expData$TSURF<0.0] <- LL1                    # If surface temperature < 0.0, then Latent heat flux equals the latent heat of sublimation


rho = 1.29*expData$PRES/1013                                      # Air density (Arrhenius substituded with rho at sea level)
e = (expData$RH/100)*0.6108*exp((21.875*expData$TAIR)/(expData$TAIR+265.5))*10      # Saturation vapour pressure [kPa](Tetens (1930) equation at level of humidity measurements) and converted to mbar
es = 0.6108*exp((21.875*expData$TSURF)/(expData$TSURF+265.5))*10           # Actual vapour pressure [kPa](at level of humidity measurements) and converted to mbar
q = e*0.622/expData$PRES                                          # Specific humidity of air [kg/kg]
qs = es*0.622/expData$PRES                                        # Specific humidity of surface [kg/kg]
Cp = 1005*(1+0.84*q)                                              # Specific heat of air
Ta_v = (expData$TAIR + 273.15) * (1 + 0.61 * q)                   # Virtual air temperature [K]
Ts_v = (expData$TSURF + 273.15) * (1 + 0.61 * qs)                 # Virtual surface temperature [K]

# Stability correction

# Richardson number
Ri_b <- (expData$TAIR - expData$TSURF) / (Zt - z0t) / ((expData$TAIR + expData$TSURF + 273.15*2)/2) / (expData$WSPD/(Zv-z0))^2*g
# stability factor
stab_f <- vector()
stab_f[which(Ri_b>0)] <- (1- 5 * Ri_b[which(Ri_b>0)])^2
stab_f[which(Ri_b<=0)] <- (1- 16 * Ri_b[which(Ri_b<=0)])^0.75
stab_f[which(abs(Ri_b)>0.2)] <- 1

#stability correction for the final equation of fluxes
c_bt <- k^2 / log(Zv/z0) / log(Zt/z0t) * stab_f

#### Calculate Turbulent fluxes

H = rho * Cp * c_bt * expData$WSPD * (expData$TAIR - expData$TSURF)      # sensible heat flux (W m-2)
LE = rho * LL * c_bt * expData$WSPD * (q - qs)                           # latent heat flux (W m-2)

expData$TSURF[expData$TSURF>=0] <- 0

#### Net Radiation

Rad_net <- expData$KINC - expData$KUPW + expData$LINC - expData$LUPW      # net radiative energy (W m-2)

Energy_n <- Rad_net + H + LE                                              # total net energy (W m-2)
#### Create snow depth vector
actualSD <- 2.15 - expData$SR50
actualSD[1:4200] <- 0
actualSD[8000:length(actualSD)] <- 0


################## Energy Balance calculations

# initialize vectors for all variables
melt <- rep(0,length(actualSD)+1)     # vector for melt [mm]
coldC <- melt                         # vector for actual cold Content [mm]
subl <- melt                          # vector for sublimation [mm]
evap <- melt                          # vector for evaporation [mm]
refr <- melt                          # vector for potential refreezing [mm]
actRefr <- melt                       # vector for actual refreezing [mm]
actStor <- melt                       # vector for actual storage - comparable to SWE [mm]
potStor <- melt                       # vector for potential storage [mm]


# initialize vectors for all cumulative variables
cumM <- melt                          # vector for cumulative melt [mm]
cumC <- melt                          # vector for cumulative cold content [mm]
cumS <- melt                          # vector for cumulative sublimation [mm]
cumE <- melt                          # vector for cumulative evaporation [mm]
cumR <- melt                          # vector for cumulative refreezing [mm]
cumD <- melt                          # vector for cumulative deposition [mm]

# add an empty first time step for the input data (this could be improved if you do a spin-up of the model)
actualSD[is.na(actualSD)] <- 0
actualSD <- c(0,actualSD) * 1000      # convert the SD from m to mm
TS <- c(0,expData$TSURF)
Energy_n <- c(0,Energy_n)
loss <- c(0,LE/LL)                    # Loss or gain due to sublimation and evaporation
# fill active Storage (SWE) with an assumed value to make sure that some initial water is present (could be replaced with field measurement)
actStor[which(TS>=0&actualSD>0)] <- 0.1

for(i in 2:length(actualSD)){

  if(actualSD[i]>=10){
  
    if(Energy_n[i]>0.0){                                   # If net energy is positive (adds energy to the surface)
      if(TS[i]>=0.0){                                   # If surface temperature is at melting point (positive net energy is used to melt snow)
        melt[i] = (Energy_n[i]/(0.334*10^6))*3600       # Calculate snow melt from positive net energy [mm h-1] 
        potStor[i] = 0.3*actualSD[i]                       # Potential storage of melt water is set to 30% of snowpack
        actStor[i] = min(potStor[i],(actStor[i]+melt[i]))  # Actual storage, limited to potential storage
        coldC[i] <- 0.0                                 # Change in cold content is 0.0 (as positive energy is used for snowmelt) [W m-2]
        
        cumM[i] = cumM[i-1] + melt[i]                   # Cumulative melt expressed in [mm h-1] or [W m-2] 
        cumC[i] = cumC[i-1] + cumC[i]                             
#browser()
      }else{                                            # If surface temperature is below melting point (positive net energy is affecting cold content and used to warm up the snowpack)
        coldC[i] = Energy_n[i]*(-1)                     # Cold content equals net energy (conditions: positive net energy and surface temperature is negative) expressed in [Wm-2]
        melt[i] <- 0.0                                  # Melt is 0.0 (because surface temperature <0.0)
        actRefr[i] <- 0.0                               # Actual refreezing 0.0 [W m-2]

        cumM[i] = cumM[i-1] + melt[i]
        cumC[i] = cumC[i-1] + coldC[i]                  # Decrease of cumulative cold content of the snowpack with positive net energy
        cumR[i] <- cumR[i-1] + actRefr[i]}
    }else{                                              # If net energy is negative or equals 0.0
      if(actStor[i]>0.0){                               # If liquid water is stored within the snowpack
        coldC[i] <- 0.0                                 # Change in cold content is 0.0
        cumC[i] = cumC[i-1] + coldC[i]
        refr[i] = (Energy_n[i]/(0.334*10^6))*3600       # Potential refreezing [mm h-1] calculated based on negative energy 
        actRefr[i] = min(refr[i],actStor[i])            # Actual refreezing, limited to available stored liquid water [mm]
        actStor[i] = max((actStor[i]-actRefr[i]),0.0)   # Update storage of liquid water with refreezing
        cumR[i] <- cumR[i-1] + actRefr[i]               # Cumulative refreezing expressed in [mm h-1] or [W m-2]
      }else{                                            # If no melt water is stored within snowpack
        coldC[i] <- Energy_n[i]*-1                      # Change in cold content equals negative net energy [W m-2]
        cumC[i] = cumC[i-1] + coldC[i]
        actRefr[i] <- 0.0                               # Actual refreezing is 0.0 (as no liquid water is available for refreezing)
        cumR[i] <- cumR[i-1] + actRefr[i]
        melt[i] = 0.0                                   # Melt is 0.0 as the net energy is negative
        cumM[i] = cumM[i-1] + melt[i]}
}
    if(TS[i]<0.0){                                      # If surface temperature smaller than 0.0, then ...

      subl[i] <- loss[i]*3600                           # Sublimation/deposition [mm]
      evap[i] <- 0.0                                    # Evaporation is 0.0 (as surface temperature is below 0.0)
      cumE[i] = cumE[i-1] + evap[i]
      if(subl[i]<0.0){                                  # If latent heat is negative, then ...
        cumS[i] = cumS[i-1] + subl[i]                   # Cumulative sublimation [mm]
      }else{                                            # If latent heat is positive, then ...
        cumD[i] = cumD[i-1] + subl[i]}                   # Cumulative deposition [mm]

      }else{                                            # If surface temperature 0.0 or higher, then
      evap[i] = loss[i]*3600                            # Evaporation/condensation [mm]
      subl[i] = 0.0}                                    # Sublimation is 0.0 (as surface temperature equals 0.0 or higher than 0.0)
    if(evap[i]<0.0){                                    # If latent heat is negative [W m-2]
      cumE[i] = cumE[i-1] + evap[i]                     # Cumulative evaporation expressed in [W m-2] or [mm]
    }else{                                              # If latent heat is positive [W m-2]
      cumD[i] = cumD[i-1] + evap[i]}                     # Cumulative condensation []
  }else{                   # If there are missing values, then ...
  subl[i] <- 0.0          # Sublimation is 0.0
  evap[i] <- 0.0          # Evaporation is 0.0
  melt[i] <- 0.0          # Melt is 0.0
  actRefr[i] <- 0.0       # Refreezing is 0.0
  coldC[i] <- 0.0}          # Change in cold content is 0.0
}

# Visualize your results

datAx1 <- as.POSIXct(TAIRhourly[,1],origin='1970-01-01')
# Zoom in on the period where snow is present

focustime <- seq(4500,8200,1)

# Do the results make sense? Let's check:

par(mar=c(3,2,1,1),mai = c(0.5, 1, 0.1, 0.2),cex.lab=1.5,cex.axis=1.5)
layout(matrix(c(1,2,3), nrow = 3, ncol = 1, byrow = T))

plot(datAx1[focustime],melt[focustime],type='l',xlab ='',ylab = expression('melt [mm]'),ylim=c(0,10),xaxt='n')
abline(h = seq(0,10,2),v = seq(datAx1[1],datAx1[length(datAx1)], by = "month"), col="gray", lty=3)

plot(datAx1[focustime],actualSD[focustime],ylim = c(0,600),xlab ='', ylab = expression('snow depth [mm]'),xaxt='n')
abline(h = seq(0,600,200),v = seq(datAx1[1],datAx1[length(datAx1)], by = "month"), col="gray", lty=3)

plot(datAx1[focustime],cumsum(melt[focustime]),type='l',xlab ='',ylab = expression('cum melt [mm]'),ylim=c(0,600),xaxt='n')
abline(h = seq(0,600,200),v = seq(datAx1[1],datAx1[length(datAx1)], by = "month"), col="gray", lty=3)

axis.POSIXct(1, at = seq(datAx1[1],datAx1[length(datAx1)], by = "month"), format = "%b-%d")





png(file=pathOutput&'//AWS_YalaBC_SNOWEB.png', res = 300,width=3600,height=1800)


par(mar=c(3,2,1,1),mai = c(0.5, 1, 0.1, 0.2),cex.lab=1.5,cex.axis=1.5)
layout(matrix(c(1,2,3,4,5), nrow = 5, ncol = 1, byrow = T))

plot(datAx1[focustime],melt[focustime],type='l',xlab ='',ylab = expression('melt [mm]'),ylim=c(0,10),xaxt='n')
abline(h = seq(0,10,2),v = seq(datAx1[1],datAx1[length(datAx1)], by = "month"), col="gray", lty=3)

plot(datAx1[focustime],expData$TAIR[focustime],type='l',lwd=2,ylim = c(-30,10),xlab ='', ylab = expression('T'['air']~~'[°C]'),xaxt='n')
abline(h = seq(-30,40,10),v = seq(datAx1[1],datAx1[length(datAx1)], by = "month"), col="gray", lty=3)

plot(datAx1[focustime],Energy_n[focustime],type='l',lwd=2,ylim = c(-300,1000),xlab ='', ylab = expression('E'['net']~~'[W '~ m^{-2}~ ']'),xaxt='n')
abline(h = seq(-300,1000,200),v = seq(datAx1[1],datAx1[length(datAx1)], by = "month"), col="gray", lty=3)

plot(datAx1[focustime],subl[focustime]*(-1),type='l',lwd=2,ylim = c(0,1),xlab ='', ylab = expression('sublimation [mm]'),xaxt='n')
abline(h = seq(0,1,0.2),v = seq(datAx1[1],datAx1[length(datAx1)], by = "month"), col="gray", lty=3)

plot(datAx1[focustime],actualSD[focustime]/1000,ylim = c(0,0.8),xlab ='', ylab = expression('snow depth [m]'),xaxt='n')
abline(h = seq(0,1,0.2),v = seq(datAx1[1],datAx1[length(datAx1)], by = "month"), col="gray", lty=3)
axis.POSIXct(1, at = seq(datAx1[1],datAx1[length(datAx1)], by = "month"), format = "%b-%d")

dev.off()
