###################################################
#This script makes use of fosa detection non-detection data from the site AJB
#in northern Madagascar. These data come from Z.J. Farris.
#
#This script will process data to be used to fit multi-state temporal occupancy models to,
#taking a stacked occupancy approach for the primary survey periods.


# 4 States:
# 1: No use
# 2: Day only
# 3: Night Only
# 4: Night and Day


#Day was assigned as the period between sunrise and sunset
#Nigt was assigned as the period between sunset and sunrise

rm(list=ls())
library(dplyr)
logit=function(x){log(x/(1-x))}
expit=function(x){exp(x)/(exp(x)+1)}

#load csv data files
f2008=read.csv("AJB Fosa/data/Diel.Occ.AJB.csv",row.names = 1,header = TRUE)
f2010=read.csv("AJB Fosa/data/Diel.Occ.2AJB.csv",row.names = 1,header = TRUE)
f2011=read.csv("AJB Fosa/data/Diel.Occ.3AJB.csv",row.names = 1,header = TRUE)
f2012=read.csv("AJB Fosa/data/Diel.Occ.4AJB.csv",row.names = 1,header = TRUE)
f2013=read.csv("AJB Fosa/data/Diel.Occ.5AJB.csv",row.names = 1,header = TRUE)
f2015=read.csv("AJB Fosa/data/Diel.Occ.6AJB.csv",row.names = 1,header = TRUE)

#test of combination to ensure order is the same of sites
all(f2008$Site==f2015$Site)

dim(f2008)
dim(f2010)
dim(f2011)
dim(f2012)
dim(f2013)
dim(f2015)


f2015=cbind(f2015,matrix(NA,nrow=dim(f2015)[1],ncol=(14-dim(f2015)[2])))
colnames(f2015)=colnames(f2008)
dim(f2015)

fosa.ajb=bind_rows(f2008,f2010,f2011,f2012,f2013,f2015)
dim(fosa.ajb)
#For a stacked occupancy design, we need to stack primary sessions
rownames(fosa.ajb)=paste(sort(rep(1:6,24)),rep(rownames(f2008),6))

head(fosa.ajb)

#Also output data as an array
Nsurvey=6
fosa.array <- array(NA, dim=c(24,14,Nsurvey))
fosa.array[,,1]=matrix(unlist(f2008),nrow=24)
fosa.array[,,2]=matrix(unlist(f2010),nrow=24)
fosa.array[,,3]=matrix(unlist(f2011),nrow=24)
fosa.array[,,4]=matrix(unlist(f2012),nrow=24)
fosa.array[,,5]=matrix(unlist(f2013),nrow=24)
fosa.array[,,6]=matrix(unlist(f2015),nrow=24)
########################################
#Load the covariate data
cov2008=read.csv("AJB Fosa/data/AJB_Hum.csv",row.names = 1,header = TRUE)
cov2010=read.csv("AJB Fosa/data/2AJB_Hum.csv",row.names = 1,header = TRUE)
cov2011=read.csv("AJB Fosa/data/3AJB_Hum.csv",row.names = 1,header = TRUE)
cov2012=read.csv("AJB Fosa/data/4AJB_Hum.csv",row.names = 1,header = TRUE)
cov2013=read.csv("AJB Fosa/data/5AJB_Hum.csv",row.names = 1,header = TRUE)
cov2015=read.csv("AJB Fosa/data/6AJB_Hum.csv",row.names = 1,header = TRUE)

#HE= Human Events
#TN = Trap Night
#TS = Trap success
#day.TS = Daytime only Trap Success
#night.TS = Nighttime only Trap Success

head(cov2008)

#THEY ARE in the correct order, just different names
#test of combination to ensure order is the same of sites
#all(rownames(cov2008)==rownames(cov2011))



dim(cov2008)
dim(cov2010)
dim(cov2011)
dim(cov2012)
dim(cov2013)
dim(cov2015)

#combine all the covariates in the same order as the data by stacking
cov.ajb=bind_rows(cov2008,cov2010,cov2011,cov2012,cov2013,cov2015)
rownames(cov.ajb)=paste(sort(rep(1:6,24)),rep(rownames(f2008),6))


#drop the Human Events and Trap nights columns
cov.ajb=cov.ajb[,-c(1,2)]


################################
#Also put covariates in an array
Nsurvey=6
cov.array <- array(NA, dim=c(24,5,Nsurvey))
cov.array[,,1]=matrix(unlist(cov2008),nrow=24)
cov.array[,,2]=matrix(unlist(cov2010),nrow=24)
cov.array[,,3]=matrix(unlist(cov2011),nrow=24)
cov.array[,,4]=matrix(unlist(cov2012),nrow=24)
cov.array[,,5]=matrix(unlist(cov2013),nrow=24)
cov.array[,,6]=matrix(unlist(cov2015),nrow=24)

##################################

#define the survey covariate
survey.cov=as.factor(sort(rep(c("2008","2010","2011","2012","2013","2015"),24)))

survey.dummy.matrix=model.matrix(~survey.cov)

#This is a contrast matrix of surveys- estimated intercept as the grand mean
survey.constrast.matrix=model.matrix(~survey.cov,contrasts = list(survey.cov = "contr.sum"))
head(survey.constrast.matrix)

#******need to add the human covarite combined across sites to the list
AJB.data=list(data=fosa.ajb,cov=cov.ajb,conmat=survey.constrast.matrix,dummat=survey.dummy.matrix,
              cov.array=cov.array,fosa.array=fosa.array)
save(AJB.data,file="AJB Fosa/AJB.data")
