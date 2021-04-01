#########################################
#
#
# Combining the detection data with datetime data to determine
#  what time of day coyote were detected
#
# Written by M. Fidino
#  
#########################################

library(lubridate)
library(tidyr)
library(dplyr)
library(suncalc)
library(stringr)

# read in the coyote datetime data
ctime <- read.csv("./Chicago coyote/data/coyote_datetime.csv")

# And the info on whether a camera was active on a given day
#  NA = not active
#   0 = active
cdet <- read.csv("./Chicago coyote/data/camera_active.csv")

# reduce ctime to only sites in cdet, which has already been
#  cleaned of sites that are not to be used for occupancy analyses.
#  (some of the sites are too close to eachother).

ctime <- ctime[which(ctime$Site %in% cdet$Site),]

# we need to make the detection data into long format so taht we can
#  connect it with ctime

cdet_long <- tidyr::pivot_longer(
  cdet,
  cols= tidyr::starts_with("Day"),
  values_to = "history"
)

# Make a new column that denotes the sampling day numerically so that we can
#  generate the sampling day from the Start column.
cdet_long$dayInt <- as.numeric(
  gsub(
    "Day_",
    "",
    cdet_long$name
  )
)

#convert Start to a date object

cdet_long$Start <- lubridate::ymd(
  cdet_long$Start,
  tz = "America/Chicago"
)


cdet_long$Date <- cdet_long$Start + lubridate::days(cdet_long$dayInt - 1)

# At time point cdet is ready to be joined to ctime.

ctime$Datetime <- lubridate::ymd_hms(ctime$Datetime, tz = "America/Chicago")

# now make an individual date column
ctime$Date <- as.Date(ctime$Datetime, tz = "America/Chicago")

# To join these two the date column in both data.frames must
#  become a character.

cdet_long$Date <- as.character(cdet_long$Date)
ctime$Date <- as.character(ctime$Date)

#join them together

coyote <- dplyr::left_join(
  cdet_long,
  ctime,
  by = c("Site", "Season", "Date")
) %>% 
  dplyr::select(
    "Site", "Lat", "Long", "Crs", "Season", "Start", 
    "name", "history", "Datetime",
  "Date") %>% 
  as.data.frame(.)


# collect the datetime info from the images

coyote$sunrise <- lubridate::as_datetime(NA, tz = "America/Chicago")
coyote$sunset <- lubridate::as_datetime(NA, tz = "America/Chicago")

# collect the sunrise and sunset time for each row in coyote
for(i in 1:nrow(coyote)){
  
  if(is.na(coyote$Datetime[i])){
    next
  }
  tmp_date <- as.Date(coyote$Datetime[i], tz = "America/Chicago")
  
  tmp <- suncalc::getSunlightTimes(
    tmp_date,
    coyote$Lat[i],
    coyote$Long[i], 
    keep = c("sunrise", "sunset"),
    tz = "America/Chicago"
  )
  coyote$sunrise[i] <- tmp$sunrise
  coyote$sunset[i] <- tmp$sunset 

  rm(tmp)
  rm(tmp_date)
}


# come up with the observation state for each row
coyote$day <- 0
coyote$night <- 0
coyote$day[is.na(coyote$history)] <- NA
coyote$night[is.na(coyote$history)] <- NA
for(i in 1:nrow(coyote)){
  if(is.na(coyote$Datetime[i])){
    next
  }
  if(
    dplyr::between(coyote$Datetime[i], coyote$sunrise[i], coyote$sunset[i])
  ){
    coyote$day[i] <- 1
  } 
  if(!dplyr::between(coyote$Datetime[i], coyote$sunrise[i], coyote$sunset[i])){
    coyote$night[i] <- 1
  }
}
# In the event that there is an NA where there should be a 0 (i.e.,
#  we detected a coyote but the pre-compiled camera activity matrix
#  has an NA in the associated cell).
coyote$history[which(coyote$night ==1 | coyote$day == 1)] <- 0


# summarise down to day. Will return -Inf if there are only NA values.
coyote <- coyote %>% 
  dplyr::group_by(Site, Season, name) %>% 
  dplyr::summarise(day = max(day, na.rm = TRUE),
                   night = max(night, na.rm = TRUE),
                   history = unique(history),
                   Date = unique(Date),
                   lat = unique(Lat),
                   Long = unique(Long),
                   Start = min(Start)) %>% 
  as.data.frame(.)

# reorder correctly, to do this I needed to pad some zeros
#  onto the numbers of the name vector
coyote$name <- paste0(
  "Day_", stringr::str_pad(
  gsub("Day_", "", coyote$name),
  2,
  pad = "0"
  )
)

# A function to order the seasons
order_seasons <- function(x){
  unq_seasons <- unique(x)
  month <- substr(unq_seasons,1,2)
  
  tmp <- sapply(month, function(y) switch(y, "JA" = 0.1,
                                          "AP" = 0.2,
                                          "JU" = 0.3,
                                          "OC" = 0.4)
  )
  
  year <- as.numeric(paste0("20", substr(unq_seasons, 3,4)))
  
  to_order <- year + tmp
  
  unq_seasons <- unq_seasons[order(to_order)]
  
  return(unq_seasons)
  
}

coyote$Season <- factor(
  coyote$Season,
  order_seasons(coyote$Season)
)

# order the data by season, site, and then sampling day
coyote <- coyote[order(coyote$Season, coyote$Site, coyote$name),]

# making the sampling category
# NA = camera not active
#  1 = coyote not detected
#  2 = coyote detected during day
#  3 = coyote detected at night
#  4 = coytoe detected during day and night
coyote$cat <- NA

# put a 0 if 0 if the -Inf was returned. This will happen when we only
#  picked up coyote during either the day or night. Following this
#  put NA values in day or night if sampling did not occur
coyote$day[is.infinite(coyote$day)] <- 0
coyote$night[is.infinite(coyote$night)] <- 0
coyote$day[is.na(coyote$history)] <- NA
coyote$night[is.na(coyote$history)] <- NA


# make the varying categories
tmp <- paste0(
  coyote$day,
  "-",
  coyote$night
)
# convert those categories to a numeric value
my_switch <- function(x){
  switch(x,
         "0-0" = 1,
         "1-0" = 2,
         "0-1" = 3,
         "1-1" = 4,
         "NA-NA" = NA
    
  )
}

coyote$cat <- as.numeric(
  sapply(
    tmp,
    my_switch
  )
)


# convert back to wide format, day night detections
dn_det <- tidyr::pivot_wider(
  coyote[,-which(colnames(coyote) %in% c("day", "night", "history","Date"))],
  names_from = name,
  values_from = cat
)

# convert down to sampling weeks to speed up the model, pull off just the
#  day columns that have the detection data
day_cols <- dn_det[,grep("^Day_",colnames(dn_det))]
# split them into seven day groups. It doesn't matter which day of the week
#  things really start on so I'm splitting everything into equal 7 day
#  samples (the last week is only 3 days so I'm removing it)
n_weeks <- floor(ncol(day_cols)/7)
week_groups <- rep(1:n_weeks, each = 7)

# a function to combine the days. If 2 & 3 happen in the same week the 
#  detection acutally becomes a 4, so there are a lot of if statements
#  in here to account for the different possible unique detections per week.
combine_days <- function(y, groups){
  ans <- rep(NA, max(groups))
  for(i in 1:length(groups)){
    tmp <- unique(as.numeric(y[groups == i]))
    if(all(is.na(tmp))){
      next
    } else {
      tmp <- tmp[!is.na(tmp)]
      if(length(tmp) == 1){
        ans[i] <- tmp
        next
      }
      if(4 %in% tmp){
        ans[i] <- 4
        next
      }
      if(2 %in% tmp & 3 %in% tmp){
        ans[i] <- 4
        next
      }
      # drop 1 now if it is in there
      if(1 %in% tmp){
        tmp <- tmp[-which(tmp == 1)]
        if(length(tmp)>1){
          stop()
        } else {
          ans[i] <- tmp
        }
      }
      
    }
  }
  return(ans)
}

# apply that function to each row of day_cols
week_summary <- t(
  apply(
    day_cols, 
    1, 
    combine_days,
    groups = week_groups
  )
)
# give logical names
colnames(week_summary) <- paste0(
  "Week_",
  stringr::str_pad(
    as.character(1:n_weeks),
    2,
    pad = "0"
  )
)
# drop the day data and tack on the week data
dn_det <- dn_det[,-grep("^Day_", colnames(dn_det))]
dn_det <- cbind(dn_det, week_summary)

# save the csv
write.csv(
  dn_det,
  "./Chicago coyote/data/day_night_detections.csv",
  row.names = FALSE
)
  