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
logit=function(x){log(x/(1-x))}
expit=function(x){exp(x)/(exp(x)+1)}

#load csv data files
f2008=read.csv("AJB Fosa/data/Fosa_Diel_Occ_2008.csv")
f2010=read.csv("AJB Fosa/data/Fosa_Diel_Occ_2010.csv")
f2011=read.csv("AJB Fosa/data/Fosa_Diel_Occ_2011.csv")
f2012=read.csv("AJB Fosa/data/Fosa_Diel_Occ_2012.csv")
f2013=read.csv("AJB Fosa/data/Fosa_Diel_Occ_2013.csv")
f2015=read.csv("AJB Fosa/data/Fosa_Diel_Occ_2015.csv")

#test of combination to ensure order is the same of sites
all(f2008$Site==f2015$Site)

dim(f2008)
dim(f2010)
dim(f2011)
dim(f2012)
dim(f2013)
dim(f2015)

#For a stacked occupancy design, we need to stack primary sessions

fosa.ajb=rbind(f2008,f2010,f2011,f2012,f2013,f2015)[,-1]
rownames(fosa.ajb)=paste(sort(rep(1:6,24)),rep(f2008[,1],6))

head(fosa.ajb)
########################################
#Load the covariate data
cov2008=read.csv("AJB Fosa/data/Human_COV_2008.csv")
cov2010=read.csv("AJB Fosa/data/Human_COV_2010.csv")
cov2011=read.csv("AJB Fosa/data/Human_COV_2011.csv")
cov2012=read.csv("AJB Fosa/data/Human_COV_2012.csv")
cov2013=read.csv("AJB Fosa/data/Human_COV_2013.csv")
cov2015=read.csv("AJB Fosa/data/Human_COV_2015.csv")

#HE= Human Events
#TN = Trap Night
#TS = Trap success
#day.TS = Daytime only Trap Success
#night.TS = Nighttime only Trap Success
#mean.TS = mean 

head(cov2008)

#test of combination to ensure order is the same of sites
all(cov2008$Site==cov2011$Site)


dim(cov2008)
dim(cov2010)
dim(cov2011)
dim(cov2012)
dim(cov2013)
dim(cov2015)

#combine all the covariates in the same order as the data by stacking
cov.ajb=rbind(cov2008,cov2010,cov2011,cov2012,cov2013,cov2015)[,-1]
rownames(cov.ajb)=paste(sort(rep(1:6,24)),rep(f2008[,1],6))

#drop the Human Events and Trap nights columns
cov.ajb=cov.ajb[,-c(1,2)]


#define the survey covariate
survey.cov=as.factor(sort(rep(c("2008","2010","2011","2012","2013","2015"),24)))

survey.dummy.matrix=model.matrix(~survey.cov)

#This is a contrast matrix of surveys- estimated intercept as the grand mean
survey.constrast.matrix=model.matrix(~survey.cov,contrasts = list(survey.cov = "contr.sum"))
head(survey.constrast.matrix)

#******need to add the human covarite combined across sites to the list
AJB.data=list(data=fosa.ajb,cov=cov.ajb,conmat=survey.constrast.matrix,dummat=survey.dummy.matrix)
save(AJB.data,file="AJB Fosa/AJB.data")
