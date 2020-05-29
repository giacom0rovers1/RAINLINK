#' ---
#' title: "RAINLINK Notebook"
#' output: html_notebook
#' ---
#' The RAINLINK package. Retrieval algorithm for rainfall mapping from microwave links
#' in a cellular communication network.
#' 
#' Version 1.14
#' Copyright (C) 2019 Aart Overeem
#' 
#' This program is free software: you can redistribute it and/or modify
#' it under the terms of the GNU General Public License as published by
#' the Free Software Foundation, either version 3 of the License, or
#' (at your option) any later version.
#' 
#' This program is distributed in the hope that it will be useful,
#' but WITHOUT ANY WARRANTY; without even the implied warranty of
#' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#' GNU General Public License for more details.
#' 
#' You should have received a copy of the GNU General Public License
#' along with this program. If not, see <http://www.gnu.org/licenses/>.
#' 
#' Note that it is not necessarily a problem if a function argument is not supplied to the function. If the
#' function argument is not used, then there is no problem. Only be aware that you should use e.g.
#' MaxFrequency=MaxFrequency. I.e. if you only supply MaxFrequency and the function argument before
#' MaxFrequency is missing, than the function will not execute properly.
#' 
#' 
#' # 0. Load R libraries, parameter values, and other settings.
#' This also loads the RAINLINK package.           
## ----Setup, include=FALSE-----------------------------------------------------
rm(list = ls())
source("Config.R") 
source("functions.R")

#' 
#' # 1. PreprocessingMinMaxRSL
#' 
## ----Data loading-------------------------------------------------------------

# Load data:
load("data/Linkdata_vodafone2016MJ.RData")
# summary(Linkdata)

# Add column with polarization if this column is not supplied in the link data:
if ("Polarization" %in% names(Linkdata)==FALSE)
{
  Linkdata$Polarization <- rep(NA,nrow(Linkdata))
}

#' 
#' When no information on polarization is provided, the above code creates a column of NA for Polarization. In the function "RainRetrievalMinMaxRSL.R" links with
#' NA values for polarization are processed with a & b values determined for vertically polarized signals.
#' If information on polarization of links is available, use H for horizontally polarized & V for vertically polarized in "Linkdata Polarization".
#' H, V & NA may occur in the same Linkdata file.
#' 
#' 
## ----Preprocessing------------------------------------------------------------
# Run R function:
StartTime <- proc.time()

DataPreprocessed <- PreprocessingMinMaxRSL(Data=Linkdata,
                                           MaxFrequency=MaxFrequency,
                                           MinFrequency=MinFrequency,
                                           verbose=TRUE)

cat(sprintf("Finished. (%.1f seconds)\n",round((proc.time()-StartTime)[3],digits=1))) # ~ 360 s

summary(DataPreprocessed)


#' 
#' 
#' # 2. WetDryNearbyLinkApMinMaxRSL
#' 
## ----Classification-----------------------------------------------------------
# Run R function:	
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

cat(sprintf("Finished. (%.1f seconds)\n",round((proc.time()-StartTime)[3],digits=1)))  # ~ 3100 s

summary(WetDry)

#' 
#' # 3. RefLevelMinMaxRSL
#' 
## ----Reference level----------------------------------------------------------
# Run R function:
StartTime <- proc.time()

Pref <- RefLevelMinMaxRSL(Data=DataPreprocessed,
                          Dry=WetDry$Dry,
                          HoursRefLevel=HoursRefLevel,
                          PeriodHoursRefLevel=PeriodHoursRefLevel)

cat(sprintf("Finished. (%.1f seconds)\n",round((proc.time()-StartTime)[3],digits=1))) # ~ 5610 s

summary(Pref)

#' 
#' # 4. OutlierFilterMinMax
## ----Outliers filter----------------------------------------------------------
# Run R function:
DataOutlierFiltered <- OutlierFilterMinMaxRSL(Data=DataPreprocessed,
                                              F=WetDry$F,
                                              FilterThreshold=FilterThreshold)

summary(DataOutlierFiltered)

#' 
#' # 5. CorrectMinMaxRSL
#' 
## ----Corrected powers---------------------------------------------------------
# Run R function:
Pcor <- CorrectMinMaxRSL(Data=DataOutlierFiltered,
                         Dry=WetDry$Dry,
                         Pref=Pref)

summary(Pcor)

#' 
#' 
#' # 6. RainRetrievalMinMaxRSL
#' 
## ----Rain retrival------------------------------------------------------------
kRPowerLawDataH <- read.table(FileRainRetrHorizontal)
colnames(kRPowerLawDataH) <- c("f", "a", "b")

kRPowerLawDataV <- read.table(FileRainRetrVertical)
colnames(kRPowerLawDataV) <- c("f", "a", "b")


# Run R function:
StartTime <- proc.time()

Rmean <- RainRetrievalMinMaxRSL(Aa=Aa,
                                alpha=alpha,
                                Data=DataOutlierFiltered,
                                kRPowerLawDataH=kRPowerLawDataH,
                                kRPowerLawDataV=kRPowerLawDataV,
                                PmaxCor=Pcor$PmaxCor,
                                PminCor=Pcor$PminCor,
                                Pref=Pref)

cat(sprintf("Finished. (%.1f seconds)\n",round((proc.time()-StartTime)[3],digits=1))) # ~ 20 s

summary(Rmean)
hist(log(Rmean))

#' 
#' 
#' # Write path-average rainfall data to files:
#' 
#' 
## ----Save to RData------------------------------------------------------------
ID <- unique(DataPreprocessed$ID)
t <- sort(unique(DataPreprocessed$DateTime))
t_sec <- as.numeric(as.POSIXct(as.character(t), format = "%Y%m%d%H%M"))
dt <- min(diff(t_sec))
save(list=ls(), file = "Cmldata_ER2016MJ.RData")

## merge in a single dataset for analyses
CmlRainfall                   <- DataPreprocessed
CmlRainfall$Pref              <- Pref
CmlRainfall$PminCor           <- Pcor$PminCor
CmlRainfall$PmaxCor           <- Pcor$PmaxCor
CmlRainfall$DryClass          <- WetDry$Dry
CmlRainfall$RainfallMeanInt   <- Rmean
CmlRainfall$RainfallDepthPath <- Rmean * dt / 3600

save(CmlRainfall, file = "CmlRainfall_ER2016MJ.RData")
save(CmlRainfall, file = "CmlRainfall_ER2016MJv2.RData", version = 2)

# write.csv(x = CmlRainfall, file = "CmlRainfall_ER2016.csv")

#' 
#' 
## ----Save to file-------------------------------------------------------------
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

#' Note that the output files contain rainfall depths (mm). If these data are to be used for the interpolation, they must first be read ("Interpolation.R" does not read these files).
#' Using the data for "Interpolation.R" requires a conversion from rainfall depth (mm) to rainfall intensity (mm/h).
#' 
#' 
#' # 7. Interpolation
#' Interpolation will be performed for hourly accumulated rainfall, so cumulative sums have to be performed
## -----------------------------------------------------------------------------
# load("CmlRainfall_ER2016.RData")

# Compute hourly accumulated rainfall as sum of the 15min rainfall depths
CmlHourlyData <- fast50x_accu1hr(CmlRainfall = CmlRainfall)

save(CmlHourlyData, file = "HourlyRainfall_ER2016MJ.RData")
summary(CmlHourlyData)

# plot(s2p(CmlRainfall$DateTime[220000:250000]), 
#      CmlRainfall$RainfallMeanInt[220000:250000], 
#      col = "red") +
#   points(s2p(CmlHourlyData$DateTime), 
#          CmlHourlyData$HourlyRainfallDepth, 
#          pch = "+")


#' 
#' 
#' Interpolation over the grid
## ---- include=FALSE-----------------------------------------------------------
# load("HourlyRainfall_ER2016.RData")

# Read grid onto which data are interpolated
RainGrid <- read.table(FileGrid, header = TRUE, sep=",")
# PolyGrid <- PolyGridGen(IntpGrid = RainGrid, SaveToFile = TRUE, FileName = FilePolygonsGrid)

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

save(RainFields, file = "IntpRainFields_ER2016MJ.RData")  # ~ 240 s

cat(sprintf("Finished. (%.1f seconds)\n",round((proc.time()-StartTime)[3],digits=1)))



#' 
#' 
#' Rasters
## ---- message=FALSE, warning=FALSE--------------------------------------------
# load("HourlyRainfall_ER2016.RData")
# load("IntpRainFields_ER2016.RData")
# RainGrid <- read.table(FileGrid, header = TRUE, sep=",")

# dimensional checks
stopifnot(dim(RainFields)[1] == length(unique(CmlHourlyData$DateTime)))
stopifnot(dim(RainFields)[2] == dim(RainGrid)[1])

# timestrings
row.names(RainFields) <- sort(unique(CmlHourlyData$DateTime))

# provinces shapefiles
borders   <- readOGR("PARMA.kml") + readOGR("BOLOGNA.kml")


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

save(RainMaps, file = "IntpRainMaps_ER2016MJ.RData")


# plot(rowSums(RainFields))
# plot(RainMaps$X201605112300)
# plot(RainMaps, "X201605112300")
plot(mask(RainMaps,borders), "X201605112300", zlim = c(0,16))


#' 
