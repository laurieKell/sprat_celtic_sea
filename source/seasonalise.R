#' @title Seasonalise, takes an annually structured `FLStock` and creates seasons
#' 
#' @description 
#'   While 'expand' adds seasons to an annually structured 'FLStock' it does not
#'   change the 'FLQuants' such as m, stock.wt, mat, stock.n, catch.n and harvest
#'   to take account of seasonal effects, such as growth and spawning. It also 
#'   does not account for changes in the population due to seasonal fishing. 
#'   Therefore 'seasonalise',  divides M into seasons, interpolates wts, and 
#'   projects a cohort across ages to estimate numbers and catch-at-age.
#'
#' @param object an \code{FLStock} object 
#' @param seasons a numeric with seasons
#' 
#' @aliases
#' 
#' @return \code{FLStock} object
#'
#' @seealso \code{\link{expand}}
#'
#' @export seasonalise
#' @docType methods
#' @rdname seasonalise
#'
#' 
#' @examples
#' \dontrun{
#' data(ple4)
#' ple4seasonal=seasonalise(ple4)
#' plot(FLStocks("Annual"=ple4,"Seasonal"=ple4seasonal)
#' }

seasonalise<-function(object, season=1:4){
  
  ## Stock and recruit                                   ###
  ## set expected to 1 and model variability as deviates ###
  sr=as.FLSR(object,model="geomean")
  params(sr)=FLPar(1,dimnames=list(params="a",iter=1))
  
  recs=FLCore:::expand(rec(object),season=season)
  
  ## Add seasons                                         ###
  object=FLCore:::expand(object,season=season)
  
  ## Divide up mortality by season                       ###
  m(      object)=wtInterp(m(object))/dim(object)[4]
  harvest(object)=harvest(object)/dim(object)[4]
  
  ## Seasonal growth                                     ###
  stock.wt(   object)=wtInterp(stock.wt(   object))
  catch.wt(   object)=wtInterp(catch.wt(   object))
  landings.wt(object)=wtInterp(landings.wt(object))
  discards.wt(object)=wtInterp(discards.wt(object))
  #m(          object)=wtInterp(          m(object))
  
  catch(object)=computeCatch(object,slot="all")

  object=adjust(object)
  
  ## Project for historic F                           ###
  #fbar=as(FLQuants("fbar"=fbar(object)[,-1]),"fwdControl")
  #object=fwd(object,control=fbar,sr=sr,residuals=recs)
  
  object}

adjust<-function (object) {
  dim = dim(object)
  un = units(catch.n(object))
  uwt = units(catch.wt(object))
  n = stock.n(object)
  m = m(object)
  f = harvest(object)
  pg = stock.n(object)[dim[1],,,dim[4]] * exp(-f[dim[1],
                                                    ,,dim[4]] - m[dim[1],,,dim[4]])
  for (i in seq(dim(object)[2] - 1)) for (j in seq(dim(object)[4])) {
    if (j != dim(object)[4]) 
      stock.n(object)[,i,,j + 1] = stock.n(object)[,
                                                      i,,j] * exp(-f[,i,,j] - m[,i,,j])
    else {
      stock.n(object)[-1,i + 1,,1] = stock.n(object)[-dim[1],
                                                        i,,j] * exp(-f[-dim[1],i,,j] - m[-dim[1],
                                                                                              i,,j])
      stock.n(object)[dim[1],i + 1,,1] = stock.n(object)[dim[1],
                                                            i + 1,,1] + pg[,i,,1]
    }
  }
  catch.n(object) = stock.n(object) * f/(m + f) * (1 - exp(-f-m))
  landings.n(object)[is.na(landings.n(object))|landings.n(object)<0] = 0
  discards.n(object)[is.na(discards.n(object))|discards.n(object)<0] = 0
  
  flag=discards.n(object)>0
  
  if (any(flag)){
     discards.n(object)[flag] = (catch.n(object) * discards.n(object)/(discards.n(object)+landings.n(object)))[flag]
     landings.n(object)[flag] = (catch.n(object) - discards.n(object))[flag]}
  else
    landings.n(object)[flag] = catch.n(object)
                                
  units(catch.n(object)) = un
  units(landings.n(object)) = un
  units(discards.n(object)) = un
  units(catch.wt(object)) = uwt
  units(landings.wt(object)) = uwt
  units(discards.wt(object)) = uwt
  catch(object) = computeCatch(object,"all")
  object}

## interpolates seasonal values based on years 
wtInterp<-function(wt){
  d4=dim(wt)[4]
  incmt=(-wt[-dim(wt)[1]]+wt[-1])%*%FLQuant(seq(0,1,length.out=d4+1)[-(d4+1)],dimnames=list(season=seq(d4)))
  
  wt[dimnames(incmt)$age]=wt[dimnames(incmt)$age]%+%incmt
  wt}

## used for extracting equilibrium values
qp<-function(stk,eql){
  
  dat=rbind.fill(
    cbind(What="Stock.wt",merge( model.frame(FLQuants("FLStock"=stock.wt(stk)),drop=T),
                                  transform(model.frame(FLQuants("FLBRP"  =stock.wt(eql)),drop=T),age=age%/%4,season=age-4*(age%/%4)+1))),
    cbind(What="Catch.wt",merge( model.frame(FLQuants("FLStock"=catch.wt(stk)),drop=T),
                                  transform(model.frame(FLQuants("FLBRP"  =catch.wt(eql)),drop=T),age=age%/%4,season=age-4*(age%/%4)+1))),
    cbind(What="Stock.n",merge( model.frame(FLQuants("FLStock"=stock.n(stk)),drop=T),
                                 transform(model.frame(FLQuants("FLBRP"  =stock.n(eql)),drop=T),age=age%/%4,season=age-4*(age%/%4)+1))),
    cbind(What="Catch.n",merge( model.frame(FLQuants("FLStock"=catch.n(stk)),drop=T),
                                 transform(model.frame(FLQuants("FLBRP"  =catch.n(eql)),drop=T),age=age%/%4,season=age-4*(age%/%4)+1))),
    cbind(What="M",      merge( model.frame(FLQuants("FLStock"=m(stk)),drop=T),
                                 transform(model.frame(FLQuants("FLBRP"  =m(eql)),drop=T),age=age%/%4,season=age-4*(age%/%4)+1))),
    cbind(What="Mat",    merge( model.frame(FLQuants("FLStock"=mat(stk)),drop=T),
                                 transform(model.frame(FLQuants("FLBRP"  =mat(eql)),drop=T),age=age%/%4,season=age-4*(age%/%4)+1))),
    cbind(What="Harvest",merge( model.frame(FLQuants("FLStock"=harvest(stk)),drop=T),
                                 transform(model.frame(FLQuants("FLBRP"  =harvest(eql)),drop=T),age=age%/%4,season=age-4*(age%/%4)+1))))
  
}


