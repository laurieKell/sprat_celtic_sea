library(FLCore)
library(FLBRP)
library(FLasher)

library(ggplotFL)
library(plyr)

library(doParallel)
library(foreach)

runCHR<-function(htar,lag,yrs,om,sr,btrig,idxSeason,devRevs,idxDevs){
  
  ## HCR function ####
  hcr<-function(iYr,btrig,htar,idx){
    
    tac =FLQuant(0,dimnames=list(quant="all",year=iYr,unit=1,area=1,season=1:4,iter=seq(100)))
    
    abc=(idx*htar)%*%qmin(idx%/%btrig,1)
    
    #### Because fishing is done from 4th to 3rd Quarters, and fwd control is from 4,1,2,3
    tac[,,,1] =abc*0.8
    tac[,,,2] =abc*0.2
    
    return(tac)}
  
  for (iYr in yrs){
    if (lag==0)
      idx=apply(stock.n(om[,ac(iYr),,idxSeason-1])%*%exp(-z(om[,ac(iYr),,idxSeason-1]))%*%stock.wt(om[,ac(iYr),,idxSeason]),6,sum)
    else
      idx=biomass(om)[,ac(iYr-lag),,idxSeason]
    
    idx=idx%*%idxDevs[,ac(iYr-lag)]
    tac=hcr(iYr,btrig,htar,idx)
    
    cnl=fwdControl(year  =c(iYr, rep(iYr+1, 3)), 
                   season=c(4,1,2,3),
                   value=tac, quant="catch")
    save(cnl,om,sr4,devRecs,file="/home/laurie/Desktop/tmp/t.RData")
    #print(777)
    om =fwd(om,control=cnl,sr=sr4,residuals=devRecs,maxF=6.0)
  }
  
  om}


runXoY<-function(rtar,lag,ref,yrs,om,sr,btrig,idxSeason,devRevs,idxDevs){
  
  ## HCR function ####
  hcr2<-function(btrig,iYr,lag,ref,catch,idx){
    
    tac= FLQuant(0,dimnames=list(quant="all",year=iYr,unit=1,area=1,season=1:4,iter=seq(100)))
    abc =rtar*apply(idx[,ac(iYr-lag)],6,mean)%*%apply(catch[,ac(ref)],6,sum)%/%apply(idx[,ac(ref)],6,mean)
    flag=c(idx[,ac(iYr-lag[1])]>btrig)
    
    if (sum(flag>0)){
      #### Because fishing is done from 4th to 3rd Quarters, and fwd control is from 4,1,2,3
      tac[,,,1,,flag] =abc[,,,,,flag]*0.8
      tac[,,,2,,flag] =abc[,,,,,flag]*0.2
    }
    
    return(tac)}
  
  for (iYr in yrs){
    
    if (lag==0)
      idx=apply(stock.n(om[,ac(iYr-lag),,idxSeason-1])%*%exp(-z(om[,ac(iYr-lag),,idxSeason-1]))%*%stock.wt(om[,ac(iYr-lag),,idxSeason]), 6.,sum)
    else
      idx=biomass(om)[,ac(iYr-lag),,idxSeason]
    idx=idx%*%idxDevs
    
    tac=hcr2(btrig,iYr,lag=lag,ref=ref,catch(om),idx)
  
    cnl =fwdControl(year  =c(iYr, rep(iYr+1, 3)), 
                    season=c(4,1,2,3),
                    value=tac, quant="catch")
    om=fwd(om,control=cnl,sr=sr4,residuals=devRecs,maxF=6.0)
    
  }  
  om0=om}

load("~/pCloudDrive/papers/inPrep/sprat/data/om/omScenarios.RData")

devRecs  =rlnorm(dim(om)[6],iter(rec(om),1)%=%0,0.3)

#### Scenarios #################################################################
scen=expand.grid(
            targetF=c(0.1,0.25,0.5,1.0,1.5,2,4,6,8)*c(refpts(eq4)["msy","harvest"])/0.1, ## /0.1 cos burnin was for 0.1
            srr    =c(TRUE,FALSE)[1])  
scen$histF=scen$targetF

runIt<-function(targetF,histF,srr){
  ################################################################################
  load("~/pCloudDrive/papers/inPrep/sprat/data/om/omScenarios.RData")
  #load("~/pCloudDrive/papers/inPrep/sprat/data/om/dev-oms.RData")
  
  ## SRR #########################################################################
  if (!srr){
    model(eq4) =geomean()$model
    params(eq4)=FLPar("a"=refpts(eq4)["virgin","rec"])
    }
  
  sr.par=FLPar(NA, dimnames=list(params=dimnames(params(eq4))$params, season=1:4, iter=1))
  sr.par[,2]=params(eq4)[,1]
  sr4=predictModel(model=model(eq4), params=sr.par)
  
  ## Target F and burn in ########################################################
  F=fbar(burnin[[1]])*histF
  F[,ac(2021:2040)]=F[,ac(2021:2040)]/histF*targetF
  control=as(FLQuants("f"=F[,-1]),"fwdControl") 

  om=FLStocks(llply(burnin[4], function(x) window(fwd(x,control=control,sr=sr4,residuals=devRecs),start=1991))) 
  
  ## Index #######################################################################
  idxSeason=4
  idxDevs  =rlnoise(100,stock(burnin[[1]])[,,,idxSeason,,1]%=%0,0.3,0)%*%
       qmin(rlnoise(100,stock(burnin[[1]])[,,,idxSeason,,1]%=%0,0.5,0.6),1)
  
  ## Run #########################################################################
  for (i in 1){ #seq(length(burnin))){
    btrig=apply((biomass(om[[i]])[,ac(2001:2020),,4]%*%idxDevs[,ac(2001:2020)]),6,min)*1.4
    htar =mean(apply(catch(om[[i]])[,ac(2001:2020)],c(2,6),sum)%/%stock(om[[i]][,ac(2001:2020),,4]))
    
    HR0=runCHR(htar,              lag=0,         2021:2039,om[[i]],sr,btrig,idxSeason,recDevs,idxDevs)
    HR1=runCHR(htar,              lag=1,         2021:2039,om[[i]],sr,btrig,idxSeason,recDevs,idxDevs)
    XY0=runXoY(rtar=targetF/histF,lag=0,ref=2010,2021:2039,om[[i]],sr,btrig,idxSeason,recDevs,idxDevs)
    XY1=runXoY(rtar=targetF/histF,lag=1,ref=2010,2021:2039,om[[i]],sr,btrig,idxSeason,recDevs,idxDevs)
    
    save(i,HR0,HR1,XY0,XY1,om,eq4,sr4,
         file=paste("/home/laurie/pCloudDrive/papers/inPrep/sprat/data/runs/mse",round(histF,2),round(targetF,2),srr,i+3,"RData",sep="."))}}

#for (i in seq(dim(scen)[1]))
#  runIt(scen[i,"targetF"],scen[i,"histF"],scen[i,"srr"])

cl= makeCluster(2)
registerDoParallel(cl)

foreach(i=3:9, #seq(dim(scen)[1]),
        .combine="c",
        .packages=c("FLCore","FLasher","FLBRP","plyr"),
        .export  =c("scen","runCHR","runCHR","runIt")) %do%{ 
  runIt(scen[i,"targetF"],scen[i,"histF"],scen[i,"srr"])}

p<-function(){
  
  p1=plot(iter(window(FLStocks("base"=om[[i]],"CHR 0"=HR0,"CHR 1"=HR1,"XY 0"=XY0,"XY 1"=XY1),
                      start=2010,end=2030),1:100),metrics=list(
                        Rec   =function(x) rec(x)[,,,2], 
                        SSB   =function(x) ssb(x)[,,,3],
                        Catch =function(x) apply(catch(x),c(2,6), sum),
                        F     =function(x) apply(fbar(x), c(2,6),function(x) mean(x)),
                        Forage=function(x) apply((stock.wt(x)%*%stock.n(x)%*%(m(x)-0.025)%/%(z(x))%*%(1-exp(-z(x))))[-1],c(2,6),sum)),
          iter=c(11,17,21))+
    facet_grid(qname~stock,scale="free")+
    theme_bw()+xlab("Year")+theme(legend.position="bottom")+
    scale_color_manual(values=rep(c("grey50"),15))+
    scale_fill_manual(values=rep(c("red"),5))
  p1+
    geom_flpar(data=FLPars(Rec   =FLPar("Rmsy"=refpts(eq4)[c("msy"),"rec",1,drop=T])*0.6718,  #M at age 0 season 2 
                           SSB   =FLPar("Bmsy"=refpts(eq4)[c("msy"),"ssb",1,drop=T]),
                           F     =FLPar("Fmsy"=refpts(eq4)[c("msy"),"harvest",1,drop=T]),
                           Catch =FLPar("MSY" =refpts(eq4)[c("msy"),"yield",1,drop=T]),
                           Forage=FLPar("MSY" =refpts(eq4)[c("msy"),"yield",1,drop=T])),x=rep(ISOdate(2012,1,1),25))
}

load("/home/laurie/Desktop/inPrep/sprat/data/mse.0.68.0.68.TRUE.3.RData");om1=HR0 
load("/home/laurie/Desktop/inPrep/sprat/data/mse.1.7.1.7.TRUE.3.RData");om2=HR0 
load("/home/laurie/Desktop/inPrep/sprat/data/mse.3.39.3.39.TRUE.3.RData");om3=HR0 
load("/home/laurie/Desktop/inPrep/sprat/data/mse.6.78.6.78.TRUE.3.RData");om4=HR0 
load("/home/laurie/Desktop/inPrep/sprat/data/mse.13.57.13.57.TRUE.3.RData");om5=HR0 
load("/home/laurie/Desktop/inPrep/sprat/data/mse.27.13.27.13.TRUE.3.RData");om6=HR0 
load("/home/laurie/Desktop/inPrep/sprat/data/mse.40.7.40.7.TRUE.3.RData");om7=HR0 
load("/home/laurie/Desktop/inPrep/sprat/data/mse.67.83.67.83.TRUE.3.RData");om8=HR0 
load("/home/laurie/Desktop/inPrep/sprat/data/mse.101.74.101.74.TRUE.3.RData");om9=HR0 

p1=plot(iter(window(FLStocks("1"=om1,"2"=om2,"3"=om3,"4"=om4,"5"=om5,"6"=om6,"7"=om7,"8"=om8,"9"=om9),
                    start=2010,end=2030),1:100),metrics=list(
                      Rec   =function(x) rec(x)[,,,3], 
                      SSB   =function(x) ssb(x)[,,,2],
                      Catch =function(x) apply(catch(x),c(2,6), sum),
                      F     =function(x) apply(fbar(x), c(2,6),function(x) mean(x)), #pmin(sum(x),1.5)),
                      Forage=function(x) apply((stock.wt(x)%*%stock.n(x)%*%(m(x)-0.025)%/%(z(x))%*%(1-exp(-z(x))))[-1],c(2,6),sum)),
        iter=c(11,17,21))+
  facet_grid(qname~stock,scale="free")+
  theme_bw()+xlab("Year")+theme(legend.position="bottom")+
  scale_color_manual(values=rep(c("grey50"),15))+
  scale_fill_manual(values=rep(c("red"),9))
p1+
  geom_flpar(data=FLPars(Rec   =FLPar("Rmsy"=refpts(eq4)[c("msy"),"rec",1,drop=T])*0.6718,  #M at age 0 season 2 
                         SSB   =FLPar("Bmsy"=refpts(eq4)[c("msy"),"ssb",1,drop=T]),
                         F     =FLPar("Fmsy"=refpts(eq4)[c("msy"),"harvest",1,drop=T]),
                         Catch =FLPar("MSY" =refpts(eq4)[c("msy"),"yield",1,drop=T]),
                         Forage=FLPar("MSY" =refpts(eq4)[c("msy"),"yield",1,drop=T])),x=rep(ISOdate(2012,1,1),45))

cnl=as(FLQuants("catch"=catch(om[,ac(2020:2039)])*10),"fwdControl")
plot(fwd(om,control=cnl,sr=sr4,residuals=devRecs,maxF=6.0))

