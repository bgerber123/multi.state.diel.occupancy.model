diel.occ.fun = function(photo.data, 
                    detection.matrix, 
                    time.zone, 
                    latitude, 
                    longitude, 
                    start.survey, 
                    end.survey, 
                    station.ID,
                    occasion.length){
  # Organize dates -------------------------------------------------------------
  # Here we combine and format the date-time in our defined time zone
  photo.data <- photo.data %>% 
    mutate(date = ymd_hms(paste(photo.data$Date,photo.data$Time), tz=time.zone)) %>% 
    arrange(date)
  
  # Transform time into 'Day' and 'Night' --------------------------------------
  # This function gives us the sunrise and sunset based on our study site's lat/long 
  sun_t = getSunlightTimes(date =(photo.data$Date), lat = latitude, 
                           lon = longitude, keep = c("sunrise", "sunset"),
                           tz=time.zone) %>% 
    select(sunrise, sunset)
  
  # Add sunrise and sunset to the metadata 
  photo.data <- tibble(cbind(photo.data, sun_t))
  
  # Create a new column, 'type,' which classifies time as 'Day' and'Night'
  photo.data <- photo.data %>%
    mutate(type=ifelse(date > sunrise & date < sunset, "Day", "Night"))
  
  # Add occasion number --------------------------------------------------------
  # create a vector with all the dates 
  start.survey = ymd(start.survey)
  end.survey = ymd(end.survey)
  
  start.survey = ymd(start.survey)
  end.survey = ymd(end.survey)
  full_date <- seq(start.survey, end.survey, by="day") 
  
  true_occ_length <- (length(full_date)/occasion.length) # true # of occasions
  occ_length <- ceiling(length(full_date)/occasion.length) # rounded occasions
  date.out=sort(rep(1:occ_length,times=occasion.length)) #index of all occasions 
  cut.this=length(date.out)-length(full_date) #this is used if true_occ_length is not an integer (see below)
  
  if (true_occ_length%%1 == 0) {
    occ <- sort(rep(1:occ_length,length.out=length(full_date)))} else
    {occ <- date.out[-c((length(date.out)-cut.this+1):length(date.out))]}
  
  full_date <- tibble(Date = full_date, occ) # combine days and occasion numbers
  photo.data <- left_join(photo.data, full_date) # add the occasion to the metadata
  
  # Build new diel detection matrix ----------------------------------------------------
  # add a column for stations (1ABJ01 to 1ABJ..)
  # add 0's to numbers from 1 to 9
  station.num <- 1:nrow(detection.matrix)
  if (length(station.num) <= 10){
    station.num <- paste("0",station.num)
  } else {
    station.num <- c(paste("0",station.num[1:9],sep=""),station.num[10:length(station.num)])
  }
  detection.matrix$Station <- paste(station.ID, station.num, sep="") 
  
  
  # longer format for the matrix: station, occasion, occasion type (0, 1, NA)
  detection.matrix <- detection.matrix %>% 
    pivot_longer(cols=matches("[0-9]"), values_to="old_occ_type", names_to="occ") %>% 
    mutate(occ = as.numeric(occ))
  
  # get if both or only 1 time is present for each station and occasion
  occasion <- photo.data %>% 
    group_by(Station,occ,type) %>% 
    summarize(rows = n()) %>% 
    arrange(Station)
  
  # assign new occasion types (2-4)
  occasion <- occasion %>% 
    group_by(Station,occ) %>% 
    do({
      p <- .
      # day + night is coded as 4
      if (nrow(p)==2){
        occ_type <- as.data.frame(4)
      }
      # only day is coded as 2
      else if (nrow(p)==1 & p$type=="Day"){
        occ_type <- as.data.frame(2)
      }
      # only night is coded as 3
      else if (nrow(p)==1 & p$type=="Night"){
        occ_type <- as.data.frame(3)
      }
    }) %>% 
    ungroup() %>% 
    pivot_longer(3:ncol(.), values_to="occ_type") %>% 
    select(-name) %>% 
    drop_na()
  
  # Join and add non detections (0 -> 1)
  diel.detection.matrix <- detection.matrix %>% 
    left_join(occasion) %>% 
    mutate(occ_type = replace(occ_type,is.na(occ_type) & !(is.na(old_occ_type)),1)) %>%
    select(-old_occ_type) %>% 
    pivot_wider(values_from=occ_type, names_from=occ)
  
  # Give your column names
  colnames(diel.detection.matrix) <- c("Station", paste("Occ_", 1:ncol(diel.detection.matrix), sep=""))
  
  #Output new diel detection matrix--------------------------------------------------
  out = list(diel.detection.matrix)
  return(out)
}
