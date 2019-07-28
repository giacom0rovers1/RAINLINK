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
}
