require(sp)
require(raster)
require(rgeos)
require(gstat)
require(rgdal)
require(geoR)
require(maptools)
require(spsann)
require(mapview)
require(dplyr)
require(ggplot2)

### For mean spatial prediction error without weights :

objUSER_annealing_non_weight <- function (points,covars) {
  sm <- as.matrix(covars[points[, 1], ])
  xtx<-(t(sm)%*%sm)
  inv.xtx<-chol2inv(chol(as.matrix(xtx)))
  eng <- vector("numeric", nrow(covars))
  for(i in 1:nrow(covars)){eng[i]<-sum(as.matrix(covars[i,])%*%inv.xtx%*%t(as.matrix(covars[i,])))}
  return(meanreng<-mean(eng))
}

schedule <- scheduleSPSANN(initial.temperature = 20,initial.acceptance = 0.60,chains = 400)
newresUSER <- optimUSER(points = point.df, fun =objUSER_annealing_non_weight,covars=covars,candi = candi,schedule = schedule,boundary = boundary,plotit = T,track = T )
objSPSANN(newresUSER)
plot(newresUSER)



### For mean spatial prediction error with poputated housing area weight :

objUSER_annealing.pop <- function (points,covars,pop) {
  sm <- as.matrix(covars[points[, 1], ])
  xtx<-(t(sm)%*%sm)
  inv.xtx<-chol2inv(chol(as.matrix(xtx)))
  eng <- vector("numeric", nrow(covars))
  for(i in 1:nrow(covars)){eng[i]<-sum(as.matrix(pop[i])%*%(as.matrix(covars[i,])%*%inv.xtx%*%t(as.matrix(covars[i,]))))}
  return(meanreng<-mean(eng))
}


schedule <- scheduleSPSANN(initial.temperature = 20,initial.acceptance = 0.60,chains = 400)
newresUSER.pop <- optimUSER(points = point.df, fun =objUSER_annealing.pop,pop=buildings_1000df$pop,covars=covars,candi = candi,schedule = schedule,boundary = boundary,track=T)
objSPSANN(newresUSER.pop)
plot(newresUSER.pop)