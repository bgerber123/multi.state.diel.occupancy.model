###################################################
#This script makes use of fosa detection non-detection data from the Ranomafana
#in Madagascar. Surveys were part of the TEAM initiative
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

load("TEAM Fosa/TEAM.fosa.stacked")

dim(fosa.stacked[[1]])
dim(fosa.stacked[[2]])
dim(fosa.stacked[[3]])
dim(fosa.stacked[[4]])
dim(fosa.stacked[[5]])
dim(fosa.stacked[[6]])

#For a stacked occupancy design, we need to stack primary sessions

fosa.TEAM=rbind(fosa.stacked[[1]],fosa.stacked[[2]],
                fosa.stacked[[3]],fosa.stacked[[4]],
                fosa.stacked[[5]],fosa.stacked[[6]])

head(fosa.TEAM)

rownames(fosa.TEAM)=paste(rownames(fosa.TEAM),c(rep("2010",60),rep("2011",60),rep("2012",60),
                            rep("2013",60),rep("2014",60),rep("2015",60)),sep="-")

#Find cameras that were not used in a given survey period and remove them
index=which(apply(fosa.TEAM,1,FUN=function(x){(all(is.na(x)))}))

fosa.TEAM=fosa.TEAM[-index,]
dim(fosa.TEAM)

#Calculate max observed state detection for all surveys
obs.max.state=rep(1,dim(fosa.TEAM)[1])
for(i in 1:dim(fosa.TEAM)[1]){
  obs.max.state[i]=max(fosa.TEAM[i,],na.rm = TRUE)
  if(length(which(fosa.TEAM[i,]==2))>0 &
  length(which(fosa.TEAM[i,]==3))>0){obs.max.state[i]=4}
}

#Count of sites (all years) in observed states
table(obs.max.state)

############################
#define the survey covariate
survey.cov=as.factor(sort(rep(c("2010","2011","2012","2013","2014","2015"),60)))

survey.dummy.matrix=model.matrix(~survey.cov)

#This is a contrast matrix of surveys- estimated intercept as the grand mean
survey.constrast.matrix=model.matrix(~survey.cov,contrasts = list(survey.cov = "contr.sum"))
head(survey.constrast.matrix)

#******need to add the human covarite combined across sites to the list
TEAM.data=list(data=fosa.TEAM,conmat=survey.constrast.matrix,dummat=survey.dummy.matrix)
save(TEAM.data,file="TEAM Fosa/TEAM.data")
