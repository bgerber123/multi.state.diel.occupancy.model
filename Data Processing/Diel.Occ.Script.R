# Function Purpose--------------------------------------------------------------
# The diel_occ function allows single-state occupancy detection matrices to 
# be converted into a 4 state matrix in order to better understand wildlife
# diel activity. Rather than 0's indicating no detection and 1's indicating
# an animal was detected, diel states are defined as the following:
# 1 = not detected
# 2 = detected during the day only (after sunrise, before sunset)
# 3 = detected during the night only (after sunset, before sunrise)
# 4 = detected during the day AND night 

# Lets get started!

# Clean Up ---------------------------------------------------------------------
rm(list=ls())
gc() #garbage collection

# Required Libraries -----------------------------------------------------------
library(tidyverse) # data manipulation 
library(lubridate) # dates 
library(suncalc) # sunlight and sunset 

# Set your variables -----------------------------------------------------------
time.zone = "Etc/GMT-3" # Time zone of your study site
latitude = -15.206893 # Latitude of study site
longitude = 49.603431 # Longitude of study site
start.survey = "2008-09-02" # Starting date of survey (YYYY-MM-DD)
end.survey = "2008-11-13" #Ending date of survey (YYYY-MM-DD)
occasion.length = 6 # Number of days per occasion 
station.ID = "AJB" # Name of Study Area

#Note- The detection matrix can be 1 (detected), 0 (not detected), and "." (not surveyed).
# The first row of the detection matrix and photo.data are column names       
detection.matrix=read_csv(file="Data Processing/Example_Detection_Matrix.csv") # Original Detection matrix

#Note- The photo.date dates need to be as YYYY-MM-DD
photo.data=read_csv(file="Data Processing/Example_Photo_Data.csv", na = ".") # CSV of 'Date', 'Time', and 'Station', of defections

# Run the function -------------------------------------------------------------
source("Data Processing/diel.occ.fun.R")

new.det.matrix=diel.occ.fun(photo.data,
                        detection.matrix, 
                        time.zone, 
                        latitude, 
                        longitude, 
                        start.survey, 
                        end.survey, 
                        station.ID, 
                        occasion.length)
 
new.det.matrix
