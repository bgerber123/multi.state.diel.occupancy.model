#This data script will prepare detection data in a format useable for model fitting
#We will compile 6 survey-area detection histories into two formats for modeling.
#The first format is useful when there are no random effects/population-level inference
#across all surveys. The second is useful for population-level inference.


#Arrange occupancy data - formatting approach 1 - for no Random Effect
#load each survey areas occurence and covariate data
FRK.occ=read.csv("Makira Fosa/FRK/Diel.Occ.FRK.csv")
FRK.cov=read.csv("Makira Fosa/FRK/FRK_Hum.csv")

LKT.occ=read.csv("Makira Fosa/LKT/Diel.Occ.LKT.csv")
LKT.cov=read.csv("Makira Fosa/LKT/LKT_Hum.csv")

MGB.occ=read.csv("Makira Fosa/MGB/Diel.Occ.MGB.csv")
MGB.cov=read.csv("Makira Fosa/MGB/MGB_Hum.csv")

SLJ.occ=read.csv("Makira Fosa/SLJ/Diel.Occ.SLJ.csv")
SLJ.cov=read.csv("Makira Fosa/SLJ/SLJ_Hum.csv")

SOA.occ=read.csv("Makira Fosa/SOA/Diel.Occ.SOA.csv")
SOA.cov=read.csv("Makira Fosa/SOA/SOA_Hum.csv")

VIN.occ=read.csv("Makira Fosa/VIN/Diel.Occ.VIN.csv")
VIN.cov=read.csv("Makira Fosa/VIN/VIN_Hum.csv")


#Look at all survey data and compile it

dim(FRK.occ)
dim(LKT.occ)
dim(MGB.occ)
dim(SLJ.occ)
dim(SOA.occ)
dim(VIN.occ)

Makira.occ=rbind(FRK.occ,LKT.occ,MGB.occ,SLJ.occ,SOA.occ,VIN.occ)
rownames(Makira.occ)=Makira.occ[,1]
Makira.occ=Makira.occ[,-1]


dim(Makira.occ)
head(Makira.occ)

Makira.occ=data.matrix(Makira.occ)

is.matrix(Makira.occ)

#Compile to covaraite data
dim(FRK.cov)
head(FRK.cov)
dim(LKT.cov)
head(LKT.cov)

Makira.cov=rbind(FRK.cov,LKT.cov,MGB.cov,SLJ.cov,SOA.cov,VIN.cov)
dim(Makira.cov)
head(Makira.cov)


rownames(Makira.cov)=Makira.cov[,1]
Makira.cov=Makira.cov[,-1]

Makira.cov=data.matrix(Makira.cov)
dim(Makira.cov)
head(Makira.cov)

#look at a histograph of day and night TS
hist(Makira.cov[,4])
hist(Makira.cov[,5],add=TRUE,col=2)

#output the data
Makira.data=list()
Makira.data[[1]]=Makira.occ
Makira.data[[2]]=Makira.cov

save(Makira.data,file="Makira Fosa/Makira.data")


################################
#Make data in an array form, used for RE models
rownames(FRK.occ)=FRK.occ[,1]
rownames(LKT.occ)=LKT.occ[,1]
rownames(MGB.occ)=MGB.occ[,1]
rownames(SLJ.occ)=SLJ.occ[,1]
rownames(SOA.occ)=SOA.occ[,1]
rownames(VIN.occ)=VIN.occ[,1]

#drop the site names here from the matrix
FRK.occ=FRK.occ[,-1]
LKT.occ=LKT.occ[,-1]
MGB.occ=MGB.occ[,-1]
SLJ.occ=SLJ.occ[,-1]
SOA.occ=SOA.occ[,-1]
VIN.occ=VIN.occ[,-1]


FRK.occ=data.matrix(FRK.occ)
LKT.occ=data.matrix(LKT.occ)
MGB.occ=data.matrix(MGB.occ)
SLJ.occ=data.matrix(SLJ.occ)
SOA.occ=data.matrix(SOA.occ)
VIN.occ=data.matrix(VIN.occ)

#Which has the max number of sites
max.sites=max(nrow(FRK.occ),
nrow(LKT.occ),
nrow(MGB.occ),
nrow(SLJ.occ),
nrow(SOA.occ),
nrow(VIN.occ))

n.occs=ncol(FRK.occ)

#make all the survey data have the same number of sites
frk=max.sites-dim(FRK.occ)[1]
lkt=max.sites-dim(LKT.occ)[1]
mgb=max.sites-dim(MGB.occ)[1]
slj=max.sites-dim(SLJ.occ)[1]
soa=max.sites-dim(SOA.occ)[1]
vin=max.sites-dim(VIN.occ)[1]

FRK.occ2=rbind(FRK.occ,matrix(NA,nrow=frk,ncol=n.occs))
LKT.occ2=rbind(LKT.occ,matrix(NA,nrow=lkt,ncol=n.occs))
#MGB.occ2=rbind(MGB.occ,matrix(NA,nrow=mgb,ncol=n.occs)) #these are the max sites already
#SLJ.occ2=rbind(MGB.occ,matrix(NA,nrow=slj,ncol=n.occs)) #these are the max sites already
SOA.occ2=rbind(SOA.occ,matrix(NA,nrow=soa,ncol=n.occs))
VIN.occ2=rbind(VIN.occ,matrix(NA,nrow=vin,ncol=n.occs))

#put occupancy surveys in an array that is consistant in sizes
Makira.occ=array(NA,dim=c(24,14,6))
Makira.occ[,,1]=FRK.occ2
Makira.occ[,,2]=LKT.occ2
Makira.occ[,,3]=MGB.occ
Makira.occ[,,4]=SLJ.occ
Makira.occ[,,5]=SOA.occ2
Makira.occ[,,6]=VIN.occ2

dim(Makira.occ)
Makira.occ[,,1]
Makira.occ[,,3]

#Do the same for the covariates
FRK.cov2=FRK.cov[,-1]
names(FRK.cov2)=NULL
LKT.cov2=LKT.cov[,-1]
names(LKT.cov2)=NULL
MGB.cov2=MGB.cov[,-1]
names(MGB.cov2)=NULL
SLJ.cov2=SLJ.cov[,-1]
names(SLJ.cov2)=NULL
SOA.cov2=SOA.cov[,-1]
names(SOA.cov2)=NULL
VIN.cov2=VIN.cov[,-1]
names(VIN.cov2)=NULL

FRK.cov2=data.matrix(FRK.cov2)
LKT.cov2=data.matrix(LKT.cov2)
MGB.cov2=data.matrix(MGB.cov2)
SLJ.cov2=data.matrix(SLJ.cov2)
SOA.cov2=data.matrix(SOA.cov2)
VIN.cov2=data.matrix(VIN.cov2)


FRK.cov2=rbind(FRK.cov2,matrix(NA,nrow=frk,ncol=5))
LKT.cov2=rbind(LKT.cov2,matrix(NA,nrow=lkt,ncol=5))
#MGB.cov2=rbind(MGB.cov2,matrix(NA,nrow=mgb,ncol=5))
#SLJ.cov2=rbind(SLJ.cov2,matrix(NA,nrow=slj,ncol=5))
SOA.cov2=rbind(SOA.cov2,matrix(NA,nrow=soa,ncol=5))
VIN.cov2=rbind(VIN.cov2,matrix(NA,nrow=vin,ncol=5))

colnames(FRK.cov2)=colnames(FRK.cov)[-1]
colnames(LKT.cov2)=colnames(LKT.cov)[-1]
colnames(MGB.cov2)=colnames(MGB.cov)[-1]
colnames(SLJ.cov2)=colnames(SLJ.cov)[-1]
colnames(SOA.cov2)=colnames(SOA.cov)[-1]
colnames(VIN.cov2)=colnames(VIN.cov)[-1]



Makira.cov=array(NA,dim=c(24,5,6))
Makira.cov[,,1]=FRK.cov2
Makira.cov[,,2]=LKT.cov2
Makira.cov[,,3]=MGB.cov2
Makira.cov[,,4]=SLJ.cov2
Makira.cov[,,5]=SOA.cov2
Makira.cov[,,6]=VIN.cov2

#output the data
Makira.data2=list()
Makira.data2[[1]]=Makira.occ
Makira.data2[[2]]=Makira.cov

save(Makira.data2,file="Makira Fosa/Makira.data2")
