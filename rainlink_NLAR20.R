#Rscript

# 0. Load R libraries, parameter values, and other settings.
rm(list = ls())
source("Config.R") 
source("functions.R")

Radius <- 20	# km

# 1. PreprocessingMinMaxRSL

# Load data:
load("data/Linkdata_vodafone2016v2.RData")
summary(Linkdata)

# Add column with polarization if this column is not supplied in the link data:
if ("Polarization" %in% names(Linkdata)==FALSE)
{
  Linkdata$Polarization <- rep(NA,nrow(Linkdata))
}


StartTime <- proc.time()

DataPreprocessed <- PreprocessingMinMaxRSL(Data=Linkdata,
                                           MaxFrequency=MaxFrequency,
                                           MinFrequency=MinFrequency,
                                           verbose=TRUE)

cat(sprintf("Preprocessing completed. (%.1f seconds)\n",round((proc.time()-StartTime)[3],digits=1))) # ~ 360 seconds



# 2. WetDryNearbyLinkApMinMaxRSL
StartTime <- proc.time()

WetDry <- WetDryNearbyLinkApMinMaxRSL(Data=DataPreprocessed,
                                      CoorSystemInputData=NULL, 
                                      MinHoursPmin=MinHoursPmin,
                                      PeriodHoursPmin=PeriodHoursPmin,
                                      Radius=Radius,
                                      Step8=Step8, 
                                      ThresholdMedian=ThresholdMedian,
                                      ThresholdMedianL=ThresholdMedianL,
                                      ThresholdNumberLinks=ThresholdNumberLinks, 
                                      ThresholdWetDry=ThresholdWetDry)

cat(sprintf("Classification completed. (%.1f seconds)\n",round((proc.time()-StartTime)[3],digits=1)))  # ~ 3100 seconds

# 3. RefLevelMinMaxRSL
StartTime <- proc.time()

Pref <- RefLevelMinMaxRSL(Data=DataPreprocessed,
                          Dry=WetDry$Dry,
                          HoursRefLevel=HoursRefLevel,
                          PeriodHoursRefLevel=PeriodHoursRefLevel)

cat(sprintf("Reference level found. (%.1f seconds)\n",round((proc.time()-StartTime)[3],digits=1))) # ~ 5610 seconds


# 4. OutlierFilterMinMax
DataOutlierFiltered <- OutlierFilterMinMaxRSL(Data=DataPreprocessed,
                                              F=WetDry$F,
                                              FilterThreshold=FilterThreshold)

# 5. CorrectMinMaxRSL
Pcor <- CorrectMinMaxRSL(Data=DataOutlierFiltered,
                         Dry=WetDry$Dry,
                         Pref=Pref)

# 6. RainRetrievalMinMaxRSL
kRPowerLawDataH <- read.table(FileRainRetrHorizontal)
colnames(kRPowerLawDataH) <- c("f", "a", "b")

kRPowerLawDataV <- read.table(FileRainRetrVertical)
colnames(kRPowerLawDataV) <- c("f", "a", "b")


StartTime <- proc.time()

Rmean <- RainRetrievalMinMaxRSL(Aa=Aa,
                                alpha=alpha,
                                Data=DataOutlierFiltered,
                                kRPowerLawDataH=kRPowerLawDataH,
                                kRPowerLawDataV=kRPowerLawDataV,
                                PmaxCor=Pcor$PmaxCor,
                                PminCor=Pcor$PminCor,
                                Pref=Pref)

cat(sprintf("Rain retrieval completed. (%.1f seconds)\n",round((proc.time()-StartTime)[3],digits=1))) # ~ 20 seconds

# Write path-average rainfall data to files:
ID <- unique(DataPreprocessed$ID)
t <- sort(unique(DataPreprocessed$DateTime))
t_sec <- as.numeric(as.POSIXct(as.character(t), format = "%Y%m%d%H%M"))
dt <- min(diff(t_sec))
save(list=ls(), file = "Cmldata_NLAR20_ER2016.RData")

## merge in a single dataset for analyses
CmlRainfall                   <- DataPreprocessed
CmlRainfall$Pref              <- Pref
CmlRainfall$PminCor           <- Pcor$PminCor
CmlRainfall$PmaxCor           <- Pcor$PmaxCor
CmlRainfall$DryClass          <- WetDry$Dry
CmlRainfall$RainfallMeanInt   <- Rmean
CmlRainfall$RainfallDepthPath <- Rmean * dt / 3600

save(CmlRainfall, file = "CmlRainfall_NLAR20_ER2016.RData")
# save(CmlRainfall, file = "CmlRainfall_NLAR20_ER2016v2.RData", version = 2)

# write.csv(x = CmlRainfall, file = "CmlRainfall_ER2016.csv")

## slow write-to-file, use tidyverse::write_delim() instead
ToFile = F
if (ToFile)
{	
  # Location of output link data:
  FolderRainEstimates <- paste("LinkPathRainDepths",TIMESTEP,"min",sep="")
  
  # Create directory for output files:
  if(!dir.exists(FolderRainEstimates)){ dir.create(FolderRainEstimates) }
  
  # Write output to file
  
  for (i in 1 : length(t))
  {
    ind <- which(DataPreprocessed$DateTime == t[i])
    int_data <- data.frame(ID = DataPreprocessed$ID[ind], 
                           RainfallDepthPath = Rmean[ind] * dt / 3600, 
                           PathLength = DataPreprocessed$PathLength[ind], 
                           XStart = DataPreprocessed$XStart[ind], 
                           YStart = DataPreprocessed$YStart[ind], 
                           XEnd = DataPreprocessed$XEnd[ind], 
                           YEnd = DataPreprocessed$YEnd[ind], 
                           IntervalNumber = rep(i, length(ind)), 
                           Frequency = DataPreprocessed$Frequency[ind])
    
    Filename <- paste(FolderRainEstimates, "/linkdata_", t[i], ".dat", sep="")
    write.table(int_data, Filename, row.names = FALSE, col.names = TRUE, append = FALSE, quote = FALSE)
  }
}

# 7. Interpolation

# load("CmlRainfall_ER2016.RData")

# Compute hourly accumulated rainfall as sum of the 15min rainfall depths
CmlHourlyData <- fast50x_accu1hr(CmlRainfall = CmlRainfall)

save(CmlHourlyData, file = "HourlyRainfall_NLAR20_ER2016.RData")

# load("HourlyRainfall_NLAR20_ER2016.RData")

# Read grid onto which data are interpolated
RainGrid <- read.table(FileGrid, header = TRUE, sep=",")

# # Location of output link data:
FolderRainMaps <- "HourlyRainMaps"

# Run R function:
StartTime <- proc.time()

RainFields <- Interpolation(Data = CmlHourlyData,
                            CoorSystemInputData = NULL,
                            idp = idp,
                            IntpMethod = IntpMethod,
                            nmax = nmax,
                            NUGGET = NUGGET,
                            RANGE = RANGE,
                            SILL = SILL,
                            Variogram = Variogram,
                            RainGrid = RainGrid,
                            Rmean = CmlHourlyData$HourlyRainfallDepth,
                            OutputDir = NULL)  # FolderRainMaps

save(RainFields, file = "IntpRainFields_NLAR20_ER2016.RData")

cat(sprintf("Interpolation finished. (%.1f seconds)\n",round((proc.time()-StartTime)[3],digits=1)))

# load("HourlyRainfall_NLAR20_ER2016.RData")
# load("IntpRainFields_NLAR20_ER2016.RData")
# RainGrid <- read.table(FileGrid, header = TRUE, sep=",")

# dimensional checks
stopifnot(dim(RainFields)[1] == length(unique(CmlHourlyData$DateTime)))
stopifnot(dim(RainFields)[2] == dim(RainGrid)[1])

# timestrings
row.names(RainFields) <- sort(unique(CmlHourlyData$DateTime))

# rain maps
mapXYZ <- RainGrid

RainMaps <- raster()
pb <- txtProgressBar(min = 0, max = nrow(RainFields), style = 3)

for(i in 1:nrow(RainFields)){
  mapXYZ$Z   <- RainFields[i,]
  rast.cml   <- rasterFromXYZ(mapXYZ, digits = 2)
  projection(rast.cml) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  # plot(rast.cml)
  RainMaps <- addLayer(RainMaps, rast.cml)
  setTxtProgressBar(pb, i)
}
names(RainMaps) <- row.names(RainFields)
close(pb)

save(RainMaps, file = "IntpRainMaps_NLAR20_ER2016.RData")
