#This data script will prepare detection data in a format useable for model fitting
#We will compile 7 sites, in which multiple surveys within a site will be stacked.


#load each survey areas occurence and covariate data
FRK.occ=read.csv("./Makira Fosa2/FRK/Diel.Occ.FRK.csv",row.names = 1,header = TRUE)
FRK.cov=read.csv("./Makira Fosa2/FRK/FRK_Hum.csv",row.names = 1,header = TRUE)

LKT.occ=read.csv("./Makira Fosa2/LKT/Diel.Occ.LKT.csv",row.names = 1,header = TRUE)
LKT.cov=read.csv("./Makira Fosa2/LKT/LKT_Hum.csv",row.names = 1,header = TRUE)

MGB.occ=read.csv("./Makira Fosa2/MGB/Diel.Occ.MGB.csv",row.names = 1,header = TRUE)
MGB.cov=read.csv("./Makira Fosa2/MGB/MGB_Hum.csv",row.names = 1,header = TRUE)

MGB2.occ=read.csv("./Makira Fosa2/MGB/Diel.Occ.2MGB.csv",row.names = 1,header = TRUE)
MGB2.cov=read.csv("./Makira Fosa2/MGB/2MGB_Hum.csv",row.names = 1,header = TRUE)

SLJ.occ=read.csv("./Makira Fosa2/SLJ/Diel.Occ.SLJ.csv",row.names = 1,header = TRUE)
SLJ.cov=read.csv("./Makira Fosa2/SLJ/SLJ_Hum.csv",row.names = 1,header = TRUE)

SOA.occ=read.csv("./Makira Fosa2/SOA/Diel.Occ.SOA.csv",row.names = 1,header = TRUE)
SOA.cov=read.csv("./Makira Fosa2/SOA/SOA_Hum.csv",row.names = 1,header = TRUE)

VIN.occ=read.csv("./Makira Fosa2/VIN/Diel.Occ.VIN.csv",row.names = 1,header = TRUE)
VIN.cov=read.csv("./Makira Fosa2/VIN/VIN_Hum.csv",row.names = 1,header = TRUE)

#AJB has 6 surveys
AJB.occ1=read.csv("./Makira Fosa2/AJB/Diel.Occ.AJB.csv",row.names = 1,header = TRUE)
AJB.occ2=read.csv("./Makira Fosa2/AJB/Diel.Occ.2AJB.csv",row.names = 1,header = TRUE)
AJB.occ3=read.csv("./Makira Fosa2/AJB/Diel.Occ.3AJB.csv",row.names = 1,header = TRUE)
AJB.occ4=read.csv("./Makira Fosa2/AJB/Diel.Occ.4AJB.csv",row.names = 1,header = TRUE)
AJB.occ5=read.csv("./Makira Fosa2/AJB/Diel.Occ.5AJB.csv",row.names = 1,header = TRUE)
AJB.occ6=read.csv("./Makira Fosa2/AJB/Diel.Occ.6AJB.csv",row.names = 1,header = TRUE)

AJB.cov1=read.csv("./Makira Fosa2/AJB/AJB_Hum.csv",row.names = 1,header = TRUE)
AJB.cov2=read.csv("./Makira Fosa2/AJB/2AJB_Hum.csv",row.names = 1,header = TRUE)
AJB.cov3=read.csv("./Makira Fosa2/AJB/3AJB_Hum.csv",row.names = 1,header = TRUE)
AJB.cov4=read.csv("./Makira Fosa2/AJB/4AJB_Hum.csv",row.names = 1,header = TRUE)
AJB.cov5=read.csv("./Makira Fosa2/AJB/5AJB_Hum.csv",row.names = 1,header = TRUE)
AJB.cov6=read.csv("./Makira Fosa2/AJB/6AJB_Hum.csv",row.names = 1,header = TRUE)



#First is to stack the multi-year within site surveys
#Stack MGB- 2 surveys
dim(MGB.occ)
dim(MGB2.occ)
MGB.occ=rbind(MGB.occ,MGB2.occ)
MGB.cov=rbind(MGB.cov,MGB2.cov)


#Stack AJB- 6 surveys
AJB.occ=rbind(AJB.occ1,AJB.occ2,AJB.occ3,AJB.occ4,AJB.occ5,AJB.occ6)
AJB.cov=rbind(AJB.cov1,AJB.cov2,AJB.cov3,AJB.cov4,AJB.cov5,AJB.cov6)


#Look at all survey data and compile it

dim(FRK.occ)
dim(LKT.occ)
dim(MGB.occ)
dim(SLJ.occ)
dim(SOA.occ)
dim(VIN.occ)
dim(AJB.occ)

################################
#Put each site's data into an array format
FRK.occ=data.matrix(FRK.occ)
LKT.occ=data.matrix(LKT.occ)
MGB.occ=data.matrix(MGB.occ)
SLJ.occ=data.matrix(SLJ.occ)
SOA.occ=data.matrix(SOA.occ)
VIN.occ=data.matrix(VIN.occ)
AJB.occ=data.matrix(AJB.occ)

#Which has the max number of sites
max.sites=max(nrow(AJB.occ),nrow(FRK.occ),
nrow(LKT.occ),
nrow(MGB.occ),
nrow(SLJ.occ),
nrow(SOA.occ),
nrow(VIN.occ)
)

#AJB has the most number of sites

#They all have 14 occs
n.occs=ncol(FRK.occ)

#make all the survey data have the same number of sites
frk=max.sites-dim(FRK.occ)[1]
lkt=max.sites-dim(LKT.occ)[1]
mgb=max.sites-dim(MGB.occ)[1]
slj=max.sites-dim(SLJ.occ)[1]
soa=max.sites-dim(SOA.occ)[1]
vin=max.sites-dim(VIN.occ)[1]

#How many sites surveyed for each?
n.sites.per.survey=c(dim(AJB.occ)[1],dim(FRK.occ)[1],dim(LKT.occ)[1],dim(MGB.occ)[1],
                     dim(SLJ.occ)[1],dim(SOA.occ)[1],dim(VIN.occ)[1])

FRK.occ2=rbind(FRK.occ,matrix(NA,nrow=frk,ncol=n.occs))
LKT.occ2=rbind(LKT.occ,matrix(NA,nrow=lkt,ncol=n.occs))
MGB.occ2=rbind(MGB.occ,matrix(NA,nrow=mgb,ncol=n.occs)) 
SLJ.occ2=rbind(SLJ.occ,matrix(NA,nrow=slj,ncol=n.occs)) 
SOA.occ2=rbind(SOA.occ,matrix(NA,nrow=soa,ncol=n.occs))
VIN.occ2=rbind(VIN.occ,matrix(NA,nrow=vin,ncol=n.occs))

#put occupancy surveys in an array that is consistant in sizes
Makira.occ=array(NA,dim=c(max.sites,n.occs,7))
Makira.occ[,,1]=AJB.occ
Makira.occ[,,2]=FRK.occ2
Makira.occ[,,3]=LKT.occ2
Makira.occ[,,4]=MGB.occ2
Makira.occ[,,5]=SLJ.occ2
Makira.occ[,,6]=SOA.occ2
Makira.occ[,,7]=VIN.occ2

dim(Makira.occ)
Makira.occ[,,1]
Makira.occ[,,3]

#####################################
#now organize the covariates to match the detection data.
#Do the same for the covariates
AJB.cov2=data.matrix(AJB.cov)
FRK.cov2=data.matrix(FRK.cov)
LKT.cov2=data.matrix(LKT.cov)
MGB.cov2=data.matrix(MGB.cov)
SLJ.cov2=data.matrix(SLJ.cov)
SOA.cov2=data.matrix(SOA.cov)
VIN.cov2=data.matrix(VIN.cov)


FRK.cov2=rbind(FRK.cov2,matrix(NA,nrow=frk,ncol=5))
LKT.cov2=rbind(LKT.cov2,matrix(NA,nrow=lkt,ncol=5))
MGB.cov2=rbind(MGB.cov2,matrix(NA,nrow=mgb,ncol=5))
SLJ.cov2=rbind(SLJ.cov2,matrix(NA,nrow=slj,ncol=5))
SOA.cov2=rbind(SOA.cov2,matrix(NA,nrow=soa,ncol=5))
VIN.cov2=rbind(VIN.cov2,matrix(NA,nrow=vin,ncol=5))


Makira.cov=array(NA,dim=c(max.sites,5,7))
Makira.cov[,,1]=AJB.cov2
Makira.cov[,,2]=FRK.cov2
Makira.cov[,,3]=LKT.cov2
Makira.cov[,,4]=MGB.cov2
Makira.cov[,,5]=SLJ.cov2
Makira.cov[,,6]=SOA.cov2
Makira.cov[,,7]=VIN.cov2

#output the data
Makira.data2=list()
Makira.data2[[1]]=n.sites.per.survey
Makira.data2[[2]]=Makira.occ
Makira.data2[[3]]=Makira.cov

save(Makira.data2,file="./Makira Fosa2/Makira.data2")
