######RDA var partioning #################
###################################################################
#author: @mitramenon
#date: 8/4/2018
#last modiied: 11/9/2019 
#NOTE: If you use your own files and directories, make sure to change paths.
        #some of the input datasets are available at 10.6084/m9.figshare.c.5130104

###################################################################

library(data.table)
library(vegan)

#Provide list of pops to remove##
removePops<-c("CHI1H","CHI1L","CHI3","COLO2L","DAV1H","DAV1L","DAV2H","DAV2L","HUA2H","HUA2L","SANRIT1L","COLO3H")


########read data in##################################
AF<-fread("AF_98pops.txt",sep="\t",data.table=F)
dim(AF)

clim<-read.table("periphery_climate1981-2010.txt",header=T,sep="\t")
clim98<-clim[!(clim$Site%in%removePops), ]
dim(clim98)

Qscore<-read.table("Qscore_clean.txt",sep=" ",header=T)

bayenvMat<-read.table("eigenVect_out251500_run1Bay.txt",header=T,sep="\t")
bayenvMat<-cbind(Site=rownames(bayenvMat),bayenvMat)
dim(bayenvMat)

PopInd<-read.table("PopInd.txt",header=T,sep="\t")
coords<-read.table("coords_chpt2.txt",header=T,sep="\t")
colnames(coords)[2]<-"Pop"
PopInd<-merge(PopInd,coords,by="PopInd")

soil<-read.table("SoilGrid1km.txt",sep="\t",header=T)
soil98<-soil[!(soil$Site%in%removePops), ]
dim(soil98)

###########data transformations and making sure pops are in the same order across all DFs############################
TargetSeq<-AF$Site

GeoData<-clim[ ,c(1,5:7)]
XYtrans<-GeoData[match(TargetSeq,GeoData$Site), ]
head(XYtrans)
XYtrans<-XYtrans[ ,-1]

XYtrans<-cbind(XYtrans, XYtrans$Latitude^2, XYtrans$Longitude^2, XYtrans$Elevation^2, XYtrans$Latitude*XYtrans$Longitude, 
XYtrans$Latitude*XYtrans$Elevation, XYtrans$Longitude*XYtrans$Elevation, XYtrans$Latitude*XYtrans$Longitude*XYtrans$Elevation)
colnames(XYtrans)<-c("Latitude","Longitude","Elevation","lat2","long2","elev2","latlong","latelv","longelv","latlongelv")
dim(XYtrans)
XYtrans<-scale(XYtrans,scale = T,center = T)

AF<-merge(GeoData,AF,by="Site")
AF<-AF[match(TargetSeq,AF$Site), ]
AF<-as.matrix(AF[ ,-c(1:4)])
AF[is.na(AF)]<-0

#determine how many SNPs are monomorphic and remove them
num<-colSums(AF == 0) 
remove<-num[num==nrow(AF)]
length(remove)
AF<-AF[ ,!(colnames(AF)%in%names(remove))]
Gdata_hellinger<-decostand(AF, method="hellinger") #standardization for count/AF data


Env98<-merge(clim98,soil98,by="Site")
Env98<-Env98[match(TargetSeq,Env98$Site), ]
Env98RM<-c("Site","Latitude.x" ,"Longitude.x", "Elevation" , "Longitude.y", "Latitude.y")
Env98_PC<-prcomp(Env98[ ,!names(Env98) %in%Env98RM],scale=T,center=T)
summary(Env98_PC)
Env98_scores<-Env98_PC$x[ ,1:7]
write.table(Env98_PC$rotation,file="EnvPCLoadings.txt",row.names=T,quote=F,sep="\t")


QscorePer<-Qscore[Qscore$Species=="periphery", ]
Qscore_hybrids<-QscorePer[!(QscorePer$Pop%in%removePops), ]

QPop<-split(Qscore_hybrids,Qscore_hybrids$Pop,drop=TRUE)
length(QPop)
estMean<-as.data.frame(unlist(lapply(QPop,function(x) return(mean(x[ ,"PSanc"],na.rm=TRUE)))))
head(estMean)
K2<-cbind(rownames(estMean),estMean)
colnames(K2)<-c("Pop","score")
Unq<-Qscore_hybrids[!(duplicated(Qscore_hybrids$Pop)), ]
K2<-merge(K2,Unq[ ,c("Pop","lat")],by="Pop") 
K2<-K2[match(TargetSeq,K2$Pop), ]
head(K2)

bayenvMat<-merge(GeoData,bayenvMat,by="Site")
bayenvMat<-bayenvMat[match(TargetSeq,bayenvMat$Site), ]
 

######################################Function for R2##########################################
RsquareAdj2.rda <- function (x, type = 'semipartial', ...)
{
m <- x$CCA$qrank
n <- nrow(x$CCA$u)
if (is.null(x$pCCA)) {
R2 <- x$CCA$tot.chi/x$tot.chi
radj <- RsquareAdj(R2, n, m)
}
else if (type == 'semipartial') {
R2 <- x$CCA$tot.chi/x$tot.chi
R2p <- x$pCCA$tot.chi/x$tot.chi
p <- x$pCCA$rank
radj <- RsquareAdj(R2 + R2p, n, m + p) - RsquareAdj(R2p, n, p)
} else if (type == 'partial'){
R2 <- x$CCA$tot.chi/(x$tot.chi - x$pCCA$tot.chi)
p <- x$pCCA$rank
radj <- 1 - (1 - R2)*(n - p - 1)/(n - m - p - 1)
if (any(na <- m >= n - 1)) radj[na] <- NA
}
list(r.squared = R2, adj.r.squared = radj, type = type)
}

##########END OF FUNCTION########################################

###Model building for variance partioning#############
#These take a fair bit of time to run. i ran it on the commputing cluster
######################################################

RDAfull<-rda(Gdata_hellinger ~ Env98_scores + bayenvMat$x + XYtrans + K2$score)
anova(RDAfull,permutations=9999)
RsquareAdj2.rda(RDAfull,type='semipartial')


RDA1<-rda(Gdata_hellinger ~ Env98_scores)
anova(RDA1,permutations=9999)
RsquareAdj2.rda(RDA1,type='semipartial')

RDA2<-rda(Gdata_hellinger ~ bayenvMat$x )
anova(RDA2,permutations=9999)
RsquareAdj2.rda(RDA2,type='semipartial')

RDA3<-rda(Gdata_hellinger ~  XYtrans)
anova(RDA3,permutations=9999)
RsquareAdj2.rda(RDA3,type='semipartial')

RDA4<-rda(Gdata_hellinger ~ K2$score)
anova(RDA4,permutations=9999)
RsquareAdj2.rda(RDA4,type='semipartial')


pRDA14<-rda(Gdata_hellinger ~ Env98_scores + K2$score+ Condition(bayenvMat$x+ XYtrans ))
RsquareAdj2.rda(pRDA14,type="semipartial")

pRDA23<-rda(Gdata_hellinger ~ bayenvMat$x + XYtrans +Condition(Env98_scores+ K2$score ))
anova(pRDA23,permutations=9999)
RsquareAdj2.rda(pRDA23,type="semipartial")


pRDA24<-rda(Gdata_hellinger ~ bayenvMat$x + K2$score+ Condition(Env98_scores+ XYtrans ))
anova(pRDA24,permutations=9999)
RsquareAdj2.rda(pRDA24,type="semipartial")


pRDA34<-rda(Gdata_hellinger ~XYtrans+K2$score+ Condition(Env98_scores+bayenvMat$x))
anova(pRDA34,permutations=9999)
RsquareAdj2.rda(pRDA34,type="semipartial")


pRDA12<-rda(Gdata_hellinger ~Env98_scores+ bayenvMat$x +Condition(XYtrans+ K2$score ))
anova(pRDA12,permutations=9999)
RsquareAdj2.rda(pRDA12,type="semipartial")

pRDA13<-rda(Gdata_hellinger ~ Env98_scores+ XYtrans +Condition( bayenvMat$x+ K2$score ))
anova(pRDA13,permutations=9999)
RsquareAdj2.rda(pRDA13,type="semipartial")


pRDA_env<-rda(Gdata_hellinger ~ Env98_scores +Condition(bayenvMat$x + K2$score + XYtrans))
anova(pRDA_env,permutations=9999)
RsquareAdj2.rda(pRDA_env,type='semipartial')


pRDA_admix<-rda(Gdata_hellinger ~K2$score +Condition(Env98_scores + bayenvMat$x + XYtrans))
anova(pRDA_admix,permutations=9999)
RsquareAdj2.rda(pRDA_admix,type='semipartial')


pRDA_str<-rda(Gdata_hellinger ~ bayenvMat$x +Condition( Env98_scores+ K2$score + XYtrans))
anova(pRDA_str,permutations=9999)
RsquareAdj2.rda(pRDA_str,type='semipartial')

pRDA_geo<-rda(Gdata_hellinger ~ XYtrans +Condition( Env98_scores+ K2$score +  bayenvMat$x  ))
anova(pRDA_geo,permutations=9999)
RsquareAdj2.rda(pRDA_geo,type='semipartial')

pRDA123<-rda(Gdata_hellinger ~ Env98_scores+bayenvMat$x + XYtrans+Condition(K2$score))
anova(pRDA123,permutations=9999)
RsquareAdj2.rda(pRDA123,type='semipartial')


pRDA134<-rda(Gdata_hellinger ~ Env98_scores+ XYtrans+ K2$score+Condition(bayenvMat$x))
anova(pRDA134,permutations=9999)
RsquareAdj2.rda(pRDA134,type='semipartial')


pRDA124<-rda(Gdata_hellinger ~ Env98_scores+bayenvMat$x+ K2$score  +Condition(XYtrans))
anova(pRDA124,permutations=9999)
RsquareAdj2.rda(pRDA124,type='semipartial')


pRDA234<-rda(Gdata_hellinger ~ bayenvMat$x + XYtrans + K2$score+Condition(Env98_scores))
anova(pRDA234,permutations=9999)
RsquareAdj2.rda(pRDA234,type='semipartial')




############MODEL COMPARISION#########
RDAInter<-rda(Gdata_hellinger ~ Env98_scores * K2$score + Condition(bayenvMat$x + XYtrans))
RsquareAdj2.rda(RDAInter,type="semipartial")


anova(RDAInter,pRDA14,permutations=9999)


