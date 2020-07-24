# Goals -----------------------------------------------------------------

# Converting capture data (1's and 0's) into diel capture data (1,2,3,4)
# 1 = 0 (not occupied), 2 = Only day, 3 = Only night, 4 = Day and Night use


# soon to come ->function
# df with columns Date, Time, Station 
# occ_number <- function(df,lat=-15.206893, lon=49.603431, start, end, 
#                        tz="Etc/GMT-3", occasion_step=6)

# Set-up -----------------------------------------------------------------

rm(list=ls())
gc() #garbage collection

library(tidyverse) # data manipulation 
library(lubridate) # dates 
library(suncalc) # sunlight and sunset 

# R project location
setwd("C:/Users/Kim Rivera/Documents/Summer Project/TOD data/Fosa")

# load metadata 
# Data Org = Date, Time (24 hr clock), Station (Camera) Name
fosa=read_csv(file="AJB_2008_TOD.csv")
fosa

# Identify and delete/ edit incomplete data  

# Organize Data Headers (if using read_csv)
# colnames(fosa)=c("Date","Time","Station")


# Organize dates ----------------------------------------------------------

# ATTENTION: dates are in Eastern African Timezone 
# if cameras malfunction and time is missing BUT
# photos are available, use light to distinguish if photos are 'Night' or 'Day'
# Here, raw data recorded 'Night' as PM and 'Day' as AM (see below)

# mutate(date_transit = case_when(is.na(date_transit) & Time=="AM" ~ ymd_hms(Date,"12:00:00"),
#                         is.na(date_transit) & Time=="PM" ~ ymd_hms(Date,"00:00:00"),
#                         !is.na(date_transit) ~ date_transit))

# combine date and time and set timezone
fosa <- fosa %>% 
  mutate(date = ymd_hms(paste(fosa$Date,fosa$Time), tz="Etc/GMT-3")) %>% 
  arrange(date)


# Transform time into 'Day' and 'Night' ---------------------------------------

# This studies Lat/ Long (Anjanaharibe): -15.206893, 49.603431
# get sunset and sunrise in EAT
sun_t = getSunlightTimes(date =(fosa$Date), lat = -15.206893, 
                                lon = 49.603431, keep = c("sunrise", "sunset"),
                                tz="Etc/GMT-3") %>% 
  select(sunrise, sunset)

# add sunrise and sunset to the metadata
fosa <- tibble(cbind(fosa, sun_t))

# Create a new column, 'type,' which classifies time as 'Day' and'Night'
fosa <- fosa %>%
  mutate(type=ifelse(date > sunrise & date < sunset, "Day", "Night"))


# Add occasion number --------------------------------------------------

# create a vector with all the dates 
# ATTENTION: specify start and end date 
start = ymd("2008-09-02")
end= ymd("2008-11-13")

full_date <- seq(start, end, by="day") 
occ_nb <- ceiling((length(full_date)/6)) # number of occasions (6 = occasion length)
occ <- sort(rep(1:occ_nb,length.out=length(full_date))) # indexes of occasions

full_date <- tibble(Date = full_date, occ) # combine days and occasion numbers
fosa <- left_join(fosa, full_date) # add the occasion to the metadata


# Update occupancy matrix ------------------------------------------------

# load capture history data 
setwd("C:/Users/Kim Rivera/Documents/Summer Project/Capture Data")
fosa_2008=read_csv(file="2008_occ.csv", na = ".") #fix missing data to NA
fosa_2008

# add a column for stations (1ABJ01 to 1ABJ..)
# add a 0 to numbers from 1 to 9
st_nb <- 1:nrow(fosa_2008)
if (length(st_nb) <= 10){
  st_nb <- paste("0",st_nb)
} else {
  st_nb <- c(paste("0",st_nb[1:9],sep=""),st_nb[10:length(st_nb)])
}
fosa_2008$Station <- paste("AJB", st_nb, sep="")


# longer format for the matrix: station, occasion, occasion type (1 or 0)
fosa_2008 <- fosa_2008 %>% 
  pivot_longer(cols=matches("[0-9]"), values_to="old_occ_type", names_to="occ") %>% 
  mutate(occ = as.numeric(occ))

# get if both or only 1 time is present for each station and occasion 
occasion <- fosa %>% 
  group_by(Station,occ,type) %>% 
  summarize(rows = n()) %>% 
  arrange(Station)

# assign new occasion types (2-4)
# 1 = 0 (not occupied), 2 = Only day, 3 = Only night, 4 = Day and Night use
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
  pivot_longer(3:5, values_to="occ_type") %>% 
  select(-name) %>% 
  drop_na()

# Join and add non detections (0 -> 1)
fosa_2008_new <- fosa_2008 %>% 
  left_join(occasion) %>% 
  mutate(occ_type = replace(occ_type,is.na(occ_type) & !(is.na(old_occ_type)),1)) %>%
  select(-old_occ_type) %>% 
  pivot_wider(values_from=occ_type, names_from=occ)

# change colnames (Occ1 to Occ..) 
colnames(fosa_2008_new)[2:length(fosa_2008_new)]<-paste("Occ_", 1:length(fosa_2008_new), sep="")


# Calculating Survey Day Difference ---------------------------------------

## 2008 end date = 
ED_2008 = as.Date("2008-11-13")
## 2010 start date = 
SD_2010 = as.Date("2010-09-16")
length(seq((ED_2008), SD_2010, by=1))-2
#this is the # of days between 11/13/2008 and 09/16/2010 (so 11/12 to 09/15)

## 2010 end date = 
ED_2010 = as.Date("2010-11-17")

## 2011 start date = 
SD_2011 = as.Date("2011-08-20")
length(seq((ED_2010), SD_2011, by=1))-2

## 2011 end date = 
ED_2011 = as.Date("2011-10-24")

## 2012 start date = 
SD_2012 = as.Date("2012-07-31")
length(seq((ED_2011), SD_2012, by=1))-2

## 2012 end date = 
ED_2012 = as.Date("2012-10-14")

## 2013 start date = 
SD_2013 = as.Date("2013-09-07")
length(seq((ED_2012), SD_2013, by=1))-2

## 2013 end date = 
ED_2013 = as.Date("2013-10-30")

## 2015 start date = 
SD_2015 = as.Date("2015-09-11")
length(seq((ED_2013), SD_2015, by=1))-2

## 2015 end date = 
ED_2015 = as.Date("2015-11-09")

