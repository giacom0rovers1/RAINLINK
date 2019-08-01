# Giacomo Roversi 2019

# functions to manipulate time strings
s2p <- function(str, dtformat="%Y%m%d%H%M"){return(as.POSIXct(str, format=dtformat, tz='GMT'))}
p2s <- function(datetime, dtformat="%Y%m%d%H%M"){return(as.character(datetime, format=dtformat))}
n2s <- function(unixtime, dtformat="%Y%m%d%H%M"){return(format(as.POSIXct(unixtime, origin = "1970-01-01", tz = 'GMT'), dtformat))}

# hourly accumulation on CMLs
accu1hr <- function(CmlRainfall, TIMESTEP){
  frac <- 60/TIMESTEP
  ids  <- unique(CmlRainfall$ID)
  dts  <- sort(unique(CmlRainfall$DateTime))
  
  totdts <- length(dts)
  totids <- length(ids)
  
  CmlHourlyData <- data.frame()
  
  
  startTime <- proc.time()
  
  for(selid in ids){                       # Ciclo sugli ID
    bool1 <- CmlRainfall$ID == selid
    numids <- which(ids == selid)
    
    elapsed_start <- round((proc.time()-startTime)[3],digits=1)
    remaining_start <- elapsed_start/numids*(totids-numids)
    
    lapTime <- proc.time()
    hourlydepth   <- 0
    count         <- 0
    
    for(seldt in dts){                         # Ciclo sui TIMESTEP
      bool2 <- CmlRainfall$DateTime == seldt
      numdts <- which(dts == seldt)
      # print(selid, quote = F)             # [debug]
      
      j <- which(bool1 & bool2)
      stopifnot(length(j) <= 1)              # Verifico univocità
      
      depth <- CmlRainfall$RainfallDepthPath[j]
      if(length(j) == 0 ){
        depth <- NA
      }
      
      hourlydepth <- hourlydepth + depth
      count       <- count + 1
      
      if(substr(seldt,11,12)=='00'          # - Cerco le ore intere
         & numdts > frac - 1                # - Attendo che sia passata almeno la prima ora
         & count == frac                    # - Ho tutte le 4 misure dell'ora
         & length(j) != 0){                 # - L'ora esiste per il link selezionato
        # print(seldt, quote = F)             # [debug]
        CmlHourlyData_row <- cbind(CmlRainfall[j,
                                               c("Frequency", "PathLength",  "XStart", 
                                                 "YStart", "XEnd", "YEnd", "Label", 
                                                 "Polarization", "Direction",  "ID")],
                                   DateTime = seldt,
                                   HourlyRainfallDepth = hourlydepth)
        
        CmlHourlyData <- rbind(CmlHourlyData, CmlHourlyData_row)
        
        hourlydepth <- 0
        count       <- 0
      }
      
      if(count > frac){        # in caso di errore (ad esempio per ora piena mancante)
        hourlydepth <- depth   # riparto con la somma delle accumulate
        count       <- 1       # come se fosse il primo intervallo
      } 
      # Remaining time estimation
      elapsed_lap <- round((proc.time()-lapTime)[3],digits=1)
      remaining_lap <- elapsed_lap/numdts*(totdts-numdts)
      cat(sprintf("\r%d:%d\t[%.1f]\t%d:%d\t[%.1f]", 
                  numdts, 
                  totdts, 
                  remaining_lap,
                  numids, 
                  totids, 
                  remaining_start))
    }
  }
  cat(sprintf("\nElapsed %.1f s", elapsed_start))
  return(CmlHourlyData)
}# VERY SLOW



# hourly accumulation on 15 min CML data
fast50x_accu1hr <- function(CmlRainfall){
  require(zoo)
  startTime <- proc.time()
  
  CmlRainfall <- CmlRainfall[order(CmlRainfall$ID, CmlRainfall$DateTime),]
  
  log1 <- substr(CmlRainfall$DateTime,11,12)=='00'
  sign <- c(T,F,F,F,T)
  
  wholehour <- rollapply(log1, width = 5, by = 1,
                         FUN = function(x) identical(x,y=sign),
                         align="right", 
                         fill = NA)
  
  samelink <- rollapply(CmlRainfall$ID, width = 4, by = 1,
                        FUN = function(x) length(unique(x)) == 1,
                        align="right", 
                        fill = NA)

  sel <- wholehour & samelink
  
  
  fouravg <- rollapplyr(CmlRainfall$RainfallMeanInt, width=4, by=1,
                        FUN = mean,
                        align="right", 
                        fill = NA)
  
  CmlHourlyData <- cbind(CmlRainfall[sel,
                                     c("DateTime","Frequency", "PathLength",  "XStart", 
                                       "YStart", "XEnd", "YEnd", "Label", 
                                       "Polarization", "Direction",  "ID")],
                         HourlyRainfallDepth = fouravg[sel])
  
  elapsed_start <- round((proc.time()-startTime)[3],digits=1)
  cat(sprintf("\nElapsed %.1f s", elapsed_start))
  
  
  filter <- which(is.na(CmlHourlyData$Frequency) 
                  & is.na(CmlHourlyData$DateTime)
                  & is.na(CmlHourlyData$PathLength)
                  & is.na(CmlHourlyData$XStart)
                  & is.na(CmlHourlyData$YStart)
                  & is.na(CmlHourlyData$XEnd)
                  & is.na(CmlHourlyData$YEnd)
                  & is.na(CmlHourlyData$Label)
                  & is.na(CmlHourlyData$Direction)
                  & is.na(CmlHourlyData$ID))
  
  require(txtplot)
  txtboxplot(CmlHourlyData$HourlyRainfallDepth[CmlHourlyData$HourlyRainfallDepth > 0.1])
  
  return(CmlHourlyData[-filter,])
}


################################################################################################

## POLYGONS GRIDS GENERATOR FUNCTION
# Creates polygons for printing function from the interpolation grid.
# Giacomo Roversi, 13 sept 2016 
# Corrected 29 oct 2017 
# Revised 31 jul 2019

PolyGridGen <- function(IntpGrid, SaveToFile){
  cat(sprintf('Loading...\n'))
  require(sp)
  require(rgdal)
  
  ok   <- 0
  jump <- 0
  
  while (ok != 1){
    width <- 1
    print(jump)
    
    # look for X-width of the grid
    while(IntpGrid[(jump+width+1),1] != IntpGrid[(jump+1),1]){
      width <- width+1
    }
    
    # search for a width larger than 10, otherwise skip to the next line
    if (width < 5){
      jump <- jump + width
      next()
    }else{
      ok <- 1
    }
    
    cat(sprintf('First acceptable line width: %d points\n', width)) 
  }
  step <- abs(IntpGrid[(jump+2),] - IntpGrid[(jump+width+1),]) # vector of the grid-side increments
  halfstep <- step/2						                               # compute the semi-diagonal
  polygons <- data.frame(X=numeric(), Y=numeric())
  totalpoints <- length(IntpGrid[,1])
  #cat(sprintf('IntpGrid dimensions: \n lon %d ° \n lat %d °\n', step[1], step[2]))
  print(step)
  
  for(i in 1:totalpoints){
    center <- IntpGrid[i,]
    r      <- 6*(i-1)
    
    polygons[r+1,] <- center + halfstep*c(-1,1)
    polygons[r+2,] <- center + halfstep*c(1,1)
    polygons[r+3,] <- center + halfstep*c(1,-1)
    polygons[r+4,] <- center + halfstep*c(-1,-1)
    polygons[r+5,] <- center + halfstep*c(-1,1)
    polygons[r+6,] <- NA
    
    cat(sprintf('Generating %d of %d polygons\r', i, totalpoints))
  }
  cat(sprintf('Generating %d of %d polygons\n', i, totalpoints))
  stopifnot(length(polygons[,1])==(6*totalpoints))
  
  if(SaveToFile){
    write.table(polygons, file='PolyGrid5x5_BO+PR.dat', row.names=F)
  }
  
  return(polygons)
}
