rm(list=ls())
setwd("C:/Users/Kim Rivera/Documents/Summer Project/Capture Data")

#############################################################################################
################################ 2008 #######################################################

fosa_2008=read.csv(file="2008_occ.csv", head=TRUE)
head(fosa_2008)

colnames(fosa_2008)<-paste("Occ", 1:ncol(fosa_2008), sep=" ")
rownames(fosa_2008)<-paste("1AJB", 1:nrow(fosa_2008), sep=" ")

library(dplyr)
fosa_2008 = na_if(fosa_2008[,], ".")

fosa_mat_2008=as.matrix(fosa_2008)
fosa_mat_2008

## Site Temporal Sub ##

#  1 = 0 (not occupied), 2 = Only day, 3 = Only night, 4 = Day and Night use
# first chage 1's to 2/3/4, then change all remaining 0's to 1's

fosa_mat_2008[12,1] = gsub("1", "3", fosa_mat_2008[12,1]) # Occ 1
fosa_mat_2008[12,2] = gsub("1", "3", fosa_mat_2008[12,2]) # Occ 2
fosa_mat_2008[12,4] = gsub("1", "3", fosa_mat_2008[12,4]) # Occ 4
fosa_mat_2008[12,5] = gsub("1", "2", fosa_mat_2008[12,5]) # Occ 5
fosa_mat_2008[12,6] = gsub("1", "4", fosa_mat_2008[12,6]) # Occ 6
fosa_mat_2008[12,7] = gsub("1", "3", fosa_mat_2008[12,7]) # Occ 7
fosa_mat_2008[12,12] = gsub("1", "3", fosa_mat_2008[12,12]) # Occ 8

fosa_mat_2008[1,6] = gsub("1", "3", fosa_mat_2008[1,6]) # Occ 6

fosa_mat_2008[2,11] = gsub("1", "3", fosa_mat_2008[2,11]) # Occ 11

fosa_mat_2008[3,8] = gsub("1", "3", fosa_mat_2008[3,8]) # Occ 8

fosa_mat_2008[8,4] = gsub("1", "4", fosa_mat_2008[8,4]) # Occ 4
fosa_mat_2008[8,8] = gsub("1", "3", fosa_mat_2008[8,8]) # Occ 8

fosa_mat_2008[9,7] = gsub("1", "2", fosa_mat_2008[9,7]) # Occ 7

fosa_mat_2008[11,6] = gsub("1", "2", fosa_mat_2008[11,6]) # Occ 10

fosa_mat_2008[13,4] = gsub("1", "4", fosa_mat_2008[13,4]) # Occ 4
fosa_mat_2008[13,6] = gsub("1", "4", fosa_mat_2008[13,6]) # Occ 6
fosa_mat_2008[13,8] = gsub("1", "2", fosa_mat_2008[13,8]) # Occ 8
fosa_mat_2008[13,10] = gsub("1", "4", fosa_mat_2008[13,10]) # Occ 10
fosa_mat_2008[13,12] = gsub("1", "2", fosa_mat_2008[13,12]) # Occ 12

fosa_mat_2008[14,4] = gsub("1", "2", fosa_mat_2008[14,4]) # Occ 4
fosa_mat_2008[14,8] = gsub("1", "3", fosa_mat_2008[14,8]) # Occ 8

fosa_mat_2008[16,13] = gsub("1", "3", fosa_mat_2008[16,13]) # Occ 13

fosa_mat_2008[18,2] = gsub("1", "3", fosa_mat_2008[18,2]) # Occ 2
fosa_mat_2008[18,3] = gsub("1", "3", fosa_mat_2008[18,3]) # Occ 3
fosa_mat_2008[18,4] = gsub("1", "3", fosa_mat_2008[18,4]) # Occ 4
fosa_mat_2008[18,6] = gsub("1", "3", fosa_mat_2008[18,6]) # Occ 6
fosa_mat_2008[18,11] = gsub("1", "2", fosa_mat_2008[18,11]) # Occ 11

fosa_mat_2008[19,3] = gsub("1", "2", fosa_mat_2008[19,3]) # Occ 3

fosa_mat_2008[20,4] = gsub("1", "3", fosa_mat_2008[20,4]) # Occ 4
fosa_mat_2008[20,8] = gsub("1", "3", fosa_mat_2008[20,8]) # Occ 8

as.matrix(fosa_mat_2008)
fosa_mat_2008[fosa_mat_2008 == 0] <-1

write.csv(fosa_mat_2008, '/Users/Kim Rivera/Documents/Summer Project/Fosa Diel Occ Excel/Fosa_Diel_Occ_2008.csv')
read.csv('/Users/Kim Rivera/Documents/Summer Project/Fosa Diel Occ Excel/Fosa_Diel_Occ_2008.csv')
#############################################################################################
################################ 2010 #######################################################

fosa_2010=read.csv(file="2010_occ.csv", head=TRUE)
head(fosa_2010)
dim(fosa_2010)

colnames(fosa_2010)<-paste("Occ", 1:ncol(fosa_2010), sep=" ")
rownames(fosa_2010)<-paste("2AJB", 1:nrow(fosa_2010), sep=" ")

library(dplyr)
fosa_2010 = na_if(fosa_2010[,], ".")

fosa_mat_2010=as.matrix(fosa_2010)
fosa_mat_2010

## Site Temporal Sub ##

#  1 = 0 (not occupied), 2 = Only day, 3 = Only night, 4 = Day and Night use
# first chage 1's to 2/3/4, then change all remaining 0's to 1's

fosa_mat_2010[3,6] = gsub("1", "2", fosa_mat_2010[3,6]) # Occ 6
fosa_mat_2010[3,7] = gsub("1", "2", fosa_mat_2010[3,7]) # Occ 7
fosa_mat_2010[3,9] = gsub("1", "3", fosa_mat_2010[3,9]) # Occ 9

fosa_mat_2010[6,3] = gsub("1", "3", fosa_mat_2010[6,3]) # Occ 3

fosa_mat_2010[11,11] = gsub("1", "2", fosa_mat_2010[11,11]) # Occ 11

fosa_mat_2010[12,6] = gsub("1", "3", fosa_mat_2010[12,6]) # Occ 6
fosa_mat_2010[12,7] = gsub("1", "4", fosa_mat_2010[12,7]) # Occ 7

fosa_mat_2010[14,7] = gsub("1", "3", fosa_mat_2010[14,7]) # Occ 7
fosa_mat_2010[14,11] = gsub("1", "0", fosa_mat_2010[14,11]) # Occ 11- MISSING OCC 11 IN TOD DATA (confirmed missing with Zach(deer cam))

fosa_mat_2010[15,1] = gsub("1", "2", fosa_mat_2010[15,1]) # Occ 1

fosa_mat_2010[16,2] = gsub("1", "3", fosa_mat_2010[16,2]) # Occ 2
fosa_mat_2010[16,7] = gsub("1", "4", fosa_mat_2010[16,7]) # Occ 7

fosa_mat_2010[18,1] = gsub("1", "2", fosa_mat_2010[18,1]) # Occ 1
fosa_mat_2010[18,8] = gsub("1", "3", fosa_mat_2010[18,8]) # Occ 8

fosa_mat_2010[19,4] = gsub("1", "3", fosa_mat_2010[19,4]) # Occ 4

fosa_mat_2010[21,2] = gsub("1", "4", fosa_mat_2010[21,2]) # Occ 2
fosa_mat_2010[21,10] = gsub("1", "2", fosa_mat_2010[21,10]) # Occ 10
fosa_mat_2010[21,11] = gsub("1", "2", fosa_mat_2010[21,11]) # Occ 11 

fosa_mat_2010[22,7] = gsub("1", "3", fosa_mat_2010[21,7]) # Occ 7
fosa_mat_2010[22,11] = gsub("1", "2", fosa_mat_2010[21,11]) # Occ 11

fosa_mat_2010[23,4] = gsub("1", "3", fosa_mat_2010[23,4]) # Occ 11 


as.matrix(fosa_mat_2010)
fosa_mat_2010[fosa_mat_2010 == 0] <-1

write.csv(fosa_mat_2010, '/Users/Kim Rivera/Documents/Summer Project/Fosa Diel Occ Excel/Fosa_Diel_Occ_2010.csv')
#############################################################################################
################################ 2011 #######################################################

fosa_2011=read.csv(file="2011_occ.csv", head=TRUE)
head(fosa_2011)
dim(fosa_2011)

colnames(fosa_2011)<-paste("Occ", 1:ncol(fosa_2011), sep=" ")
rownames(fosa_2011)<-paste("3AJB", 1:nrow(fosa_2011), sep=" ")

library(dplyr)
fosa_2011 = na_if(fosa_2011[,], ".")

fosa_mat_2011=as.matrix(fosa_2011)
fosa_mat_2011

## Site Temporal Sub ##

#  1 = 0 (not occupied), 2 = Only day, 3 = Only night, 4 = Day and Night use
# first chage 1's to 2/3/4, then change all remaining 0's to 1's

fosa_mat_2011[2,8] = gsub("1", "2", fosa_mat_2011[2,8]) # Occ 8

fosa_mat_2011[8,4] = gsub("1", "3", fosa_mat_2011[8,4]) # Occ 4

fosa_mat_2011[13,2] = gsub("1", "2", fosa_mat_2011[13,2]) # Occ 2
fosa_mat_2011[13,7] = gsub("1", "3", fosa_mat_2011[13,7]) # Occ 7

fosa_mat_2011[14,3] = gsub("1", "2", fosa_mat_2011[14,3]) # Occ 3
fosa_mat_2011[14,4] = gsub("1", "3", fosa_mat_2011[14,4]) # Occ 4
fosa_mat_2011[14,8] = gsub("1", "3", fosa_mat_2011[14,8]) # Occ 8

fosa_mat_2011[15,10] = gsub("1", "2", fosa_mat_2011[15,10]) # Occ 10

fosa_mat_2011[16,8] = gsub("1", "3", fosa_mat_2011[16,8]) # Occ 8
fosa_mat_2011[16,9] = gsub("1", "4", fosa_mat_2011[16,9]) # Occ 9
fosa_mat_2011[16,11] = gsub("1", "3", fosa_mat_2011[16,11]) # Occ 11

fosa_mat_2011[22,1] = gsub("1", "3", fosa_mat_2011[22,1]) # Occ 1
fosa_mat_2011[22,7] = gsub("1", "3", fosa_mat_2011[22,7]) # Occ 7
fosa_mat_2011[22,8] = gsub("1", "3", fosa_mat_2011[22,8]) # Occ 8
fosa_mat_2011[22,10] = gsub("1", "3", fosa_mat_2011[22,10]) # Occ 10

fosa_mat_2011[24,5] = gsub("1", "3", fosa_mat_2011[24,5]) # Occ 5

as.matrix(fosa_mat_2011)
fosa_mat_2011[fosa_mat_2011 == 0] <-1

write.csv(fosa_mat_2011, '/Users/Kim Rivera/Documents/Summer Project/Fosa Diel Occ Excel/Fosa_Diel_Occ_2011.csv')
#############################################################################################
################################ 2012 #######################################################

fosa_2012=read.csv(file="2012_occ.csv", head=TRUE)
head(fosa_2012)
dim(fosa_2012)

colnames(fosa_2012)<-paste("Occ", 1:ncol(fosa_2012), sep=" ")
rownames(fosa_2012)<-paste("4AJB", 1:nrow(fosa_2012), sep=" ")

library(dplyr)
fosa_2012 = na_if(fosa_2012[,], ".")

fosa_mat_2012=as.matrix(fosa_2012)
fosa_mat_2012

## Site Temporal Sub ##

#  1 = 0 (not occupied), 2 = Only day, 3 = Only night, 4 = Day and Night use
# first chage 1's to 2/3/4, then change all remaining 0's to 1's

fosa_mat_2012[1,12] = gsub("1", "3", fosa_mat_2012[1,12]) # Occ 12

fosa_mat_2012[3,12] = gsub("1", "3", fosa_mat_2012[3,12]) # Occ 12

fosa_mat_2012[6,3] = gsub("1", "2", fosa_mat_2012[6,3]) # Occ 3

fosa_mat_2012[7,3] = gsub("1", "2", fosa_mat_2012[7,3]) # Occ 3
fosa_mat_2012[7,4] = gsub("1", "3", fosa_mat_2012[7,4]) # Occ 4
fosa_mat_2012[7,5] = gsub("1", "4", fosa_mat_2012[7,5]) # Occ 5
fosa_mat_2012[7,10] = gsub("1", "3", fosa_mat_2012[7,10]) # Occ 10

fosa_mat_2012[8,1] = gsub("1", "3", fosa_mat_2012[8,1]) # Occ 1

fosa_mat_2012[9,5] = gsub("1", "3", fosa_mat_2012[9,5]) # Occ 5

fosa_mat_2012[11,4] = gsub("1", "2", fosa_mat_2012[11,4]) # Occ 4

fosa_mat_2012[13,10] = gsub("1", "3", fosa_mat_2012[13,10]) # Occ 10
fosa_mat_2012[13,13] = gsub("1", "2", fosa_mat_2012[13,13]) # Occ 13

fosa_mat_2012[14,4] = gsub("1", "3", fosa_mat_2012[14,4]) # Occ 4

fosa_mat_2012[15,1] = gsub("1", "3", fosa_mat_2012[15,1]) # Occ 1
fosa_mat_2012[15,5] = gsub("1", "3", fosa_mat_2012[15,5]) # Occ 5
fosa_mat_2012[15,6] = gsub("1", "3", fosa_mat_2012[15,6]) # Occ 6
fosa_mat_2012[15,10] = gsub("1", "3", fosa_mat_2012[15,10]) # Occ 10

fosa_mat_2012[17,4] = gsub("1", "3", fosa_mat_2012[17,4]) # Occ 4

fosa_mat_2012[18,4] = gsub("1", "3", fosa_mat_2012[18,4]) # Occ 4

fosa_mat_2012[20,2] = gsub("1", "2", fosa_mat_2012[20,2]) # Occ 2

fosa_mat_2012[21,5] = gsub("1", "3", fosa_mat_2012[21,5]) # Occ 2

as.matrix(fosa_mat_2012)
fosa_mat_2012[fosa_mat_2012 == 0] <-1

write.csv(fosa_mat_2012, '/Users/Kim Rivera/Documents/Summer Project/Fosa Diel Occ Excel/Fosa_Diel_Occ_2012.csv')
#############################################################################################
################################ 2013 #######################################################

fosa_2013=read.csv(file="2013_occ.csv", head=TRUE)
head(fosa_2013)
dim(fosa_2013)

colnames(fosa_2013)<-paste("Occ", 1:ncol(fosa_2013), sep=" ")
rownames(fosa_2013)<-paste("5AJB", 1:nrow(fosa_2013), sep=" ")

library(dplyr)
fosa_2013 = na_if(fosa_2013[,], ".")

fosa_mat_2013=as.matrix(fosa_2013)
fosa_mat_2013
fosa_mat_2013[13,9] = gsub("1", "0", fosa_mat_2013[13,9]) # Occ 9
#the above is due to lack of photographic evidence there was an occasion at this time

fosa_mat_2013[1,2] = gsub("1", "3", fosa_mat_2013[1,2]) # Occ 2
fosa_mat_2013[1,3] = gsub("1", "3", fosa_mat_2013[1,3]) # Occ 3
fosa_mat_2013[1,5] = gsub("1", "3", fosa_mat_2013[1,5]) # Occ 5
fosa_mat_2013[1,6] = gsub("1", "3", fosa_mat_2013[1,6]) # Occ 6
fosa_mat_2013[1,8] = gsub("1", "3", fosa_mat_2013[1,8]) # Occ 8

fosa_mat_2013[2,5] = gsub("1", "3", fosa_mat_2013[2,5]) # Occ 5
fosa_mat_2013[2,6] = gsub("1", "2", fosa_mat_2013[2,6]) # Occ 6
fosa_mat_2013[2,8] = gsub("1", "3", fosa_mat_2013[2,8]) # Occ 8

fosa_mat_2013[3,5] = gsub("1", "4", fosa_mat_2013[3,5]) # Occ 5

fosa_mat_2013[6,8] = gsub("1", "3", fosa_mat_2013[6,8]) # Occ 8

fosa_mat_2013[9,7] = gsub("1", "2", fosa_mat_2013[9,7]) # Occ 7
fosa_mat_2013[9,8] = gsub("1", "3", fosa_mat_2013[9,8]) # Occ 8

fosa_mat_2013[11,7] = gsub("1", "3", fosa_mat_2013[11,7]) # Occ 7

fosa_mat_2013[12,8] = gsub("1", "3", fosa_mat_2013[12,8]) # Occ 8
fosa_mat_2013[12,9] = gsub("1", "3", fosa_mat_2013[12,9]) # Occ 9

fosa_mat_2013[13,3] = gsub("1", "3", fosa_mat_2013[13,3]) # Occ 3
fosa_mat_2013[13,4] = gsub("1", "3", fosa_mat_2013[13,4]) # Occ 4
fosa_mat_2013[13,8] = gsub("1", "3", fosa_mat_2013[13,8]) # Occ 8
#Extra  in site 13 occ 9--removed after not finding photo evidence

fosa_mat_2013[14,4] = gsub("1", "2", fosa_mat_2013[14,4]) # Occ 4
fosa_mat_2013[14,6] = gsub("1", "4", fosa_mat_2013[14,6]) # Occ 6
fosa_mat_2013[14,8] = gsub("1", "3", fosa_mat_2013[14,8]) # Occ 8
fosa_mat_2013[14,9] = gsub("1", "3", fosa_mat_2013[14,9]) # Occ 9

fosa_mat_2013[15,3] = gsub("1", "3", fosa_mat_2013[15,3]) # Occ 3
fosa_mat_2013[15,7] = gsub("1", "3", fosa_mat_2013[15,7]) # Occ 7
fosa_mat_2013[15,9] = gsub("1", "3", fosa_mat_2013[15,9]) # Occ 9

fosa_mat_2013[16,2] = gsub("1", "3", fosa_mat_2013[16,2]) # Occ 2

fosa_mat_2013[18,7] = gsub("1", "3", fosa_mat_2013[18,7]) # Occ 7
fosa_mat_2013[18,9] = gsub("1", "2", fosa_mat_2013[18,9]) # Occ 9

fosa_mat_2013[20,1] = gsub("1", "3", fosa_mat_2013[20,1]) # Occ 20

fosa_mat_2013[21,2] = gsub("1", "3", fosa_mat_2013[21,2]) # Occ 2
fosa_mat_2013[21,4] = gsub("1", "2", fosa_mat_2013[21,4]) # Occ 4
fosa_mat_2013[21,5] = gsub("1", "3", fosa_mat_2013[21,5]) # Occ 5

fosa_mat_2013[23,2] = gsub("1", "2", fosa_mat_2013[23,2]) # Occ 2
fosa_mat_2013[23,5] = gsub("1", "4", fosa_mat_2013[23,5]) # Occ 5
fosa_mat_2013[23,7] = gsub("1", "3", fosa_mat_2013[23,7]) # Occ 7

fosa_mat_2013[24,7] = gsub("1", "2", fosa_mat_2013[24,7]) # Occ 7

as.matrix(fosa_mat_2013)
fosa_mat_2013[fosa_mat_2013 == 0] <-1

write.csv(fosa_mat_2013, '/Users/Kim Rivera/Documents/Summer Project/Fosa Diel Occ Excel/Fosa_Diel_Occ_2013.csv')
#############################################################################################
################################ 2015 #######################################################

fosa_2015=read.csv(file="2015_occ.csv", head=TRUE)
head(fosa_2015)
dim(fosa_2015)

colnames(fosa_2015)<-paste("Occ", 1:ncol(fosa_2015), sep=" ")
rownames(fosa_2015)<-paste("6AJB", 1:nrow(fosa_2015), sep=" ")

library(dplyr)
fosa_2015 = na_if(fosa_2015[,], ".")

fosa_mat_2015=as.matrix(fosa_2015)
fosa_mat_2015

fosa_mat_2015[2,1] = gsub("1", "3", fosa_mat_2015[2,1]) # Occ 1
fosa_mat_2015[2,2] = gsub("1", "3", fosa_mat_2015[2,2]) # Occ 2
fosa_mat_2015[2,3] = gsub("1", "3", fosa_mat_2015[2,3]) # Occ 3

fosa_mat_2015[4,6] = gsub("1", "3", fosa_mat_2015[4,6]) # Occ 6
fosa_mat_2015[4,7] = gsub("1", "3", fosa_mat_2015[4,7]) # Occ 7

fosa_mat_2015[6,1] = gsub("1", "3", fosa_mat_2015[6,1]) # Occ 1

fosa_mat_2015[7,2] = gsub("1", "3", fosa_mat_2015[7,2]) # Occ 2

fosa_mat_2015[9,1] = gsub("1", "2", fosa_mat_2015[9,1]) # Occ 1

fosa_mat_2015[12,4] = gsub("1", "2", fosa_mat_2015[12,4]) # Occ 4

fosa_mat_2015[13,2] = gsub("1", "3", fosa_mat_2015[13,2]) # Occ 2
fosa_mat_2015[13,3] = gsub("1", "3", fosa_mat_2015[13,3]) # Occ 3

fosa_mat_2015[15,7] = gsub("1", "3", fosa_mat_2015[15,7]) # Occ 7

fosa_mat_2015[17,2] = gsub("1", "3", fosa_mat_2015[17,2]) # Occ 2
fosa_mat_2015[17,4] = gsub("1", "3", fosa_mat_2015[17,4]) # Occ 4

fosa_mat_2015[21,2] = gsub("1", "4", fosa_mat_2015[21,2]) # Occ 2
fosa_mat_2015[21,3] = gsub("1", "3", fosa_mat_2015[21,3]) # Occ 3
fosa_mat_2015[21,7] = gsub("1", "2", fosa_mat_2015[21,7]) # Occ 7

fosa_mat_2015[23,2] = gsub("1", "2", fosa_mat_2015[23,2]) # Occ 2
fosa_mat_2015[23,5] = gsub("1", "3", fosa_mat_2015[23,5]) # Occ 5
fosa_mat_2015[23,6] = gsub("1", "4", fosa_mat_2015[23,6]) # Occ 6
fosa_mat_2015[23,7] = gsub("1", "3", fosa_mat_2015[23,7]) # Occ 7
fosa_mat_2015[23,8] = gsub("1", "3", fosa_mat_2015[23,8]) # Occ 8
fosa_mat_2015[23,9] = gsub("1", "2", fosa_mat_2015[23,9]) # Occ 9

fosa_mat_2015[24,1] = gsub("1", "2", fosa_mat_2015[24,1]) # Occ 1

as.matrix(fosa_mat_2015)
fosa_mat_2015[fosa_mat_2015 == 0] <-1

write.csv(fosa_mat_2015, '/Users/Kim Rivera/Documents/Summer Project/Fosa Diel Occ Excel/Fosa_Diel_Occ_2015.csv')
#############################################################################################
################################ DAYS BTW YEARS #############################################

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

## List of Matrices ##

as.list(fosa_mat_2008)
as.list(fosa_mat_2010)
as.list(fosa_mat_2011)
as.list(fosa_mat_2012)
as.list(fosa_mat_2013)
as.list(fosa_mat_2015)
