library(ggplot2); theme_set(theme_bw())
library(ggpubr)
library(ggcorrplot)

library(plyr)
library(dplyr)
library(reshape)

library(FLCore)
library(FLBRP)
library(FLasher)
library(FLife)

library(ggplotFL)

library(statcomp)

dirOM  =file.path(getwd(),"data/om")
dirRes =file.path(getwd(),"data/results")
dirRuns=file.path(getwd(),"data/runs")

load(file.path(dirOM,"om.RData"))

mDet=iter(m(om4),1)

#### Uncertainty ###############################################################
devs=cbind(iter=seq(100),expand.grid(sd=seq(0.05,0.5,length.out=10),
                                     b =seq(0,   1.0,length.out=10))) 

## M Cohort
cDev=mdply(devs, 
      function(iter,sd,b) mdply(data.frame(cohort=dimnames(FLCohort(m(om4)))$cohort), function(cohort) 
        data.frame(cohort=cohort,expand.grid(dimnames(m(om4))[c(4,1)]),data=FLife:::noiseFn(prod(dim(m(om4))[c(1,4)]),sd=sd,b=sd))))
cDev=exp(as.FLQuant(transmute(cDev,year=an(ac(cohort))+an(ac(age)),age=age,iter=iter,season=season,data=data))[,dimnames(m(om4))$year])

## M Trend
tDev=as.FLQuant(mdply(devs,
               function(iter,sd,b) as.data.frame(rlnoise(1,m(om4)[1,,,1,,1]%=%0,sd=sd,b=b),drop=T))[,c("year","iter","data")])

## Recruitment
tm=as.FLQuant(mdply(devs,
                    function(iter,sd,b) as.data.frame(rlnoise(1,m(om4)[1,,,1,,1]%=%0,sd=sd,b=b),drop=T))[,c("year","iter","data")])
rDev=rec(om4)
rDev[]=tm

################################################################################

#### M trend ###################################################################
m(om4)=mDet%*%tDev

control=as(FLQuants("f"=fbar(om4)[,-1,,,,50]%=%0),"fwdControl")  
om4.annual.0=fwd(om4,control=control,sr=sr4)

control=as(FLQuants("f"=fbar(om4)[,-1,,,,50]%=%refpts(eq4)["msy","harvest"]/4),"fwdControl")  
om4.annual.1=fwd(om4,control=control,sr=sr4)

control=as(FLQuants("f"=fbar(om4)[,-1,,,,50]%=%refpts(eq4)["msy","harvest"]*0.125),"fwdControl")  
om4.annual.half=fwd(om4,control=control,sr=sr4)

plot(FLStocks("0"=om4.annual.0,"0.5"=om4.annual.half,"1"=om4.annual.1))
################################################################################

#### M Cohort ##################################################################
m(om4)=mDet%*%cDev

control=as(FLQuants("f"=fbar(om4)[,-1,,,,50]%=%refpts(eq4)["msy","harvest"]*0),"fwdControl")  
om4.cohort.0=fwd(om4,control=control,sr=sr4)

control=as(FLQuants("f"=fbar(om4)[,-1,,,,50]%=%refpts(eq4)["msy","harvest"]/4),"fwdControl")  
om4.cohort.1=fwd(om4,control=control,sr=sr4)

control=as(FLQuants("f"=fbar(om4)[,-1,,,,50]%=%refpts(eq4)["msy","harvest"]/8),"fwdControl")  
om4.cohort.half=fwd(om4,control=control,sr=sr4)

plot(FLStocks("0"=om4.cohort.0,"0.5"=om4.cohort.half,"1"=om4.cohort.1))
#################################################################################

#### Recruits ##################################################################
m(om4)=mDet

control=as(FLQuants("f"=fbar(om4)[,-1,,,,50]%=%0),"fwdControl")  
om4.rec.0=fwd(om4,control=control,sr=sr4,residuals=rDev)

control=as(FLQuants("f"=fbar(om4)[,-1,,,,50]%=%refpts(eq4)["msy","harvest"]/4),"fwdControl")  
om4.rec.1=fwd(om4,control=control,sr=sr4,residuals=rDev)

control=as(FLQuants("f"=fbar(om4)[,-1,,,,50]%=%refpts(eq4)["msy","harvest"]/8),"fwdControl")  
om4.rec.half=fwd(om4,control=control,sr=sr4,residuals=rDev)

plot(FLStocks("0"=om4.rec.0,"0.5"=om4.rec.half,"1"=om4.rec.1))
#################################################################################

save(om4,cDev,tDev,rDev,
     om4.annual.1    ,om4.annual.0,    om4.annual.half,
     om4.cohort.1    ,om4.cohort.0,    om4.cohort.half,
     om4.rec.1       ,om4.rec.0,       om4.rec.half,
     eq4,file=file.path(dirOM,"om4Devs.RData"))

#### All #######################################################################
load(file.path(dirOM,"om4.RData"))
om4 =iter(om4,-1)

rC<-sample(seq(100),100,T)
rT<-sample(seq(100),100,T)
rR<-sample(seq(100),100,T)

om4   =iter(om4,50)
m(om4)=mDet%*%tDev[,,,,,rT]%*%cDev[,,,,,rC]

control=as(FLQuants("f"=fbar(om4)[,-1]%=%0),"fwdControl")  
om4.all.0=fwd(om4,control=control,sr=sr4,residuals=rDev[,,,,,rR])

control=as(FLQuants("f"=fbar(om4)[,-1]%=%refpts(eq4)["msy","harvest"]/4),"fwdControl")  
om4.all.1=fwd(om4,control=control,sr=sr4,residuals=rDev[,,,,,rR])

control=as(FLQuants("f"=fbar(om4)[,-1]%=%refpts(eq4)["msy","harvest"]/8),"fwdControl")  
om4.all.half=fwd(om4,control=control,sr=sr4,residuals=rDev[,,,,,rR])

plot(FLStocks("0"=om4.all.0,"0.5"=om4.all.half,"1"=om4.all.1))
#################################################################################

save(rC,rT,rR,
     om4.all.1       ,om4.all.0,       om4.all.half,
     eq4,file=file.path(dirOM,"allDevs.RData"))

load(file.path(dirOM,"allDevs.RData"))
smry=rbind(cbind("Fmsy"=0.0,"Uncertainty"="all",omStock(om4.all.0)),
           cbind("Fmsy"=0.5,"Uncertainty"="all",omStock(om4.all.half)),
           cbind("Fmsy"=1.0,"Uncertainty"="all",omStock(om4.all.1)),
          
           cbind("Fmsy"=0.0,"Uncertainty"="Recruits",omStock(om4.rec.0)),
           cbind("Fmsy"=0.5,"Uncertainty"="Recruits",omStock(om4.rec.half)),
           cbind("Fmsy"=1.0,"Uncertainty"="Recruits",omStock(om4.rec.1)),
          
           cbind("Fmsy"=0.0,"Uncertainty"="Annual",omStock(om4.annual.0)),
           cbind("Fmsy"=0.5,"Uncertainty"="Annual",omStock(om4.annual.half)),
           cbind("Fmsy"=1.0,"Uncertainty"="Annual",omStock(om4.annual.1)),
          
           cbind("Fmsy"=0.0,"Uncertainty"="Cohort",omStock(om4.cohort.0)),
           cbind("Fmsy"=0.5,"Uncertainty"="Cohort",omStock(om4.cohort.half)),
           cbind("Fmsy"=1.0,"Uncertainty"="Cohort",omStock(om4.cohort.1)))

save(smry,file=file.path(dirOM,"results.RData"))


#### Fishing Trends ############################################################
load(file.path(dirOM,"om4.RData"))

burnin =propagate(window(iter(om4,51),start=1980,end=2040),100)
control=as(FLQuants("f"=fbar(burnin)[,-1]),"fwdControl") 
burnin =window(fwd(burnin,control=control,sr=sr4),#residuals=rDev),
               start=2001)

fnl=c(seq(.2,1,length.out=51)[-51],seq(1,2,length.out=50))
f  =as.FLQuant(transform(mdply(data.frame(final=fnl), function(final) data.frame(year=2021:2040,data=seq(1,final,length.out=20))),iter=factor(final))[,-1])
dimnames(f)$iter=1:100
f=fbar(burnin)[,ac(2021:2040)]%*%f

control=as(FLQuants("f"=f),"fwdControl")  
burnin.1=fwd(burnin,control=control,sr=sr4)#,residuals=rDev[,,,,,rR])

p=plot(window(burnin.1),metrics=list(
  Rec  =function(x) rec(x)[,,,3],
  SSB  =function(x) ssb(x)[,,,2],
  Catch=function(x) apply(catch(x),c(2,6), sum),
  F    =function(x) apply(fbar(x), c(2,6),mean)))+
  geom_flpar(data=FLPars(Rec  =FLPar("Rmsy"=refpts(eq4)[c("msy"),"rec",drop=T])*0.6,
                         SSB  =FLPar("Bmsy"=refpts(eq4)[c("msy"),"ssb",drop=T]),
                         F    =FLPar("Fmsy"=refpts(eq4)[c("msy"),"harvest",drop=T]),
                         Catch=FLPar("MSY" =refpts(eq4)[c("msy"),"yield",drop=T])),x=rep(c(2000),4))+
  theme_bw()+xlab("Year")+theme(legend.position="bottom")
p

#### Stochastic ################################################################
load(file.path(dirOM,"om4.RData"))

burnin =propagate(window(iter(om4,51),start=1980,end=2040),100)
control=as(FLQuants("f"=fbar(burnin)[,-1]),"fwdControl") 
burnin =window(fwd(burnin,control=control,sr=sr4,residuals=rDev),start=2001)

fnl=c(seq(.2,1,length.out=51)[-51],seq(1,2,length.out=50))
f  =as.FLQuant(transform(mdply(data.frame(final=fnl), function(final) data.frame(year=2021:2040,data=seq(1,final,length.out=20))),iter=factor(final))[,-1])
dimnames(f)$iter=1:100
f=fbar(burnin)[,ac(2021:2040)]%*%f

control=as(FLQuants("f"=f),"fwdControl")  
burnin.1=fwd(burnin,control=control,sr=sr4,residuals=rDev[,,,,,rR])

burnin.1=fwd(burnin,control=control,sr=sr4,residuals=rDev[,,,,,rR])

p=plot(window(burnin.1,start=2000),metrics=list(
  Rec  =function(x) rec(x)[,,,3],
  SSB  =function(x) ssb(x)[,,,2],
  Catch=function(x) apply(catch(x),c(2,6), sum),
  F    =function(x) apply(fbar(x), c(2,6),mean)))+
  geom_flpar(data=FLPars(Rec  =FLPar("Rmsy"=refpts(eq4)[c("msy"),"rec",1,drop=T])*0.6,
                         SSB  =FLPar("Bmsy"=refpts(eq4)[c("msy"),"ssb",1,drop=T]),
                         F    =FLPar("Fmsy"=refpts(eq4)[c("msy"),"harvest",1,drop=T]),
                         Catch=FLPar("MSY" =refpts(eq4)[c("msy"),"yield",1,drop=T])),x=rep(c(2000),4))+
  theme_bw()+xlab("Year")+theme(legend.position="bottom")
p

plot(FLStocks("Recs"=om4.rec.1[,,,2,,30],"M Annual"=om4.annual.1[,,,2,,30],"M Cohort"=om4.cohort.1[,,,2,,30]))+
  scale_x_continuous(limits=c(1981,2030))+
  theme(legend.position="bottom")

pe<-function(x,y)
  log(ebiomass(x)[,-1,,1]%/%ebiomass(x)[,-dim(x)[2],,1])

load(file.path(dirOM,"om4Devs.RData"))
     
pe.rec.0   =procerr(window(om4.rec.0,   start=1960),eq4)
pe.rec.half=procerr(window(om4.rec.half,start=1960),eq4)
pe.rec.1   =procerr(window(om4.rec.1,   start=1960),eq4)

pe.cohort.0   =procerr(window(om4.cohort.0,   start=1960),eq4)
pe.cohort.half=procerr(window(om4.cohort.half,start=1960),eq4)
pe.cohort.1   =procerr(window(om4.cohort.1,   start=1960),eq4)

pe.annual.0   =procerr(window(om4.annual.0,   start=1960),eq4)
pe.annual.half=procerr(window(om4.annual.half,start=1960),eq4)
pe.annual.1   =procerr(window(om4.annual.1,   start=1960),eq4)

pe.all.0   =procerr(window(om4.all.0,   start=1960),eq4)
pe.all.half=procerr(window(om4.all.half,start=1960),eq4)
pe.all.1   =procerr(window(om4.all.1,   start=1960),eq4)

dat=rbind(cbind(Deviates="Recruits","Fmsy"=0,  as.data.frame(apply(pe.rec.0,      c(4,6),var,na.rm=T)^0.5,drop=T)),
          cbind(Deviates="Recruits","Fmsy"=1,  as.data.frame(apply(pe.rec.1,      c(4,6),var,na.rm=T)^0.5,drop=T)),
          cbind(Deviates="Recruits","Fmsy"=0.5,as.data.frame(apply(pe.rec.half,   c(4,6),var,na.rm=T)^0.5,drop=T)),
          cbind(Deviates="Cohort",  "Fmsy"=0,  as.data.frame(apply(pe.cohort.0,   c(4,6),var,na.rm=T)^0.5,drop=T)),
          cbind(Deviates="Cohort",  "Fmsy"=1,  as.data.frame(apply(pe.cohort.1,   c(4,6),var,na.rm=T)^0.5,drop=T)),
          cbind(Deviates="Cohort",  "Fmsy"=0.5,as.data.frame(apply(pe.cohort.half,c(4,6),var,na.rm=T)^0.5,drop=T)),
          cbind(Deviates="Annual",  "Fmsy"=1,  as.data.frame(apply(pe.annual.0,   c(4,6),var,na.rm=T)^0.5,drop=T)),
          cbind(Deviates="Annual",  "Fmsy"=0,  as.data.frame(apply(pe.annual.1,   c(4,6),var,na.rm=T)^0.5,drop=T)),
          cbind(Deviates="Annual",  "Fmsy"=0.5,as.data.frame(apply(pe.annual.half,c(4,6),var,na.rm=T)^0.5,drop=T)))

dat=cbind(dat,expand.grid(sd=seq(0.05,0.5,length.out=10),
                          b =round(seq(0,   1.0,length.out=10),2))[dat$iter,])
dat=transform(subset(dat,!is.na(b)),Fmsy=paste("Fmsy", Fmsy),b=paste("AR",b))
ggplot(dat)+
  geom_point(aes(y=data,x=sd,col=Deviates))+facet_grid(Fmsy~b)+
  ylab("Process Error")+xlab("Deviates SD")+
  theme(legend.position="bottom")



om=FLStockR(burnin.1,refpts=refpts(eq4)[,,1])

p=plot(window(burnin.1,start=2000),metrics=list(
  Rec  =function(x) rec(x)[,,,3],
  SSB  =function(x) ssb(x)[,,,2],
  Catch=function(x) apply(catch(x),c(2,6), sum),
  F    =function(x) apply(fbar(x), c(2,6),mean)))+
  geom_flpar(data=FLPars(Rec  =FLPar("Rmsy"=refpts(eq4)[c("msy"),"rec",1,drop=T])*0.6,
                         SSB  =FLPar("Bmsy"=refpts(eq4)[c("msy"),"ssb",1,drop=T]),
                         F    =FLPar("Fmsy"=refpts(eq4)[c("msy"),"harvest",1,drop=T]),
                         Catch=FLPar("MSY" =refpts(eq4)[c("msy"),"yield",1,drop=T])),x=rep(c(2000),4))+
  theme_bw()+xlab("Year")+theme(legend.position="bottom")
p

om=FLStockR(window(burnin.1,start=2020),refpts=refpts(eq4)[,,1])
x =FLQuants(om, f=function(x) fbar(x)%/%refpts(x)["msy","harvest"])
y =indicators.len(om[,,,3], indicators=c('lbar', 'lmean'), params=pars[,1], metric='catch.n')

dat=mdply(data.frame(Year=seq(2022,2040)), function(Year) roc(x$f[,ac(Year),,3]<1,y$lmean[,ac(Year)]))
ggplot(dat)+
  geom_line(aes(FPR,TPR,col=Year,group=Year))


#  geom_point(aes(FPR,TPR),data=rtn[abs(rtn$ind-1)==min(abs(rtn$ind-1)),],col="red",size=3)

transform(rtn[abs(rtn$ind-1)==min(abs(rtn$ind-1)),],roc.TSS=TPR-FPR) 

ggplot(rtn)+
  geom_line(aes(ind,TSS))+
  geom_vline(aes(xintercept=1),col="red")


#Concerning faster-growing species, the conclusion was that trend-based management procedures (the ICES 2 over 3 rule, the rfb-rule or 
#any other combination of x over y rules with or without additional elements such as uncertainty caps or biomass safeguards) lead
#to poor management performance (high risks, low yields) for such species and should be avoided. The only way to comply with 
#precautionary principles for such rules and species is to apply very precautionary multipliers (very low catch advice). #
#Consequently, the recommendation would be to very cautious with trend-based rules for faster-growing species and consider
#abandoning them. Instead, alternative management procedures (e.g. harvest rate-based rules or escapement strategies) should be #
#explored for faster-growing species





  
