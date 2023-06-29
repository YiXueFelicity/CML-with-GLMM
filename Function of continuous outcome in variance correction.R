
#####new pi ###
calc_pi = function(Y,Z,delta)
{
  #Z <- Z12
  int <-  cbind(1,Y,Z)
  pivalue <-  1/(1+exp(-as.vector(int%*%delta)))
  return(pivalue)
} 



############ new der_beta


der_betaM <- function(X,Z,beta){
  #Z=Z12
  #X=X12
  
  int <- cbind(1,X,Z)
  
  const <-  rep(1, nrow(int))
  
  Bvalue  <-int*const
  
  
  return(Bvalue) 
}



der_beta <- function(X,Z,beta){
  #Z=Z12
  #X=X12
  
  int <- cbind(1,X,Z)
  
  #eta <- as.vector(int%*%beta)
  
  int.list <- split(int,seq(nrow(int)))
  
  const <-  rep(1, nrow(int))
  
  Bvalue.l  <-Map ('*', int.list , const)
  
  
  return(Bvalue.l) 
}


############ new der_delta



derv_deltaM=function(Y,Z,delta){
  
  ZZ <- cbind(1,Y,Z)
  
  pi1 <- calc_pi(Y,Z,delta)
  
  Dvalue <- ZZ*pi1*(1-pi1)
  
  return(Dvalue)
  
}

derv_delta=function(Y,Z,delta){
  
  ZZ <- cbind(1,Y,Z)
  
  pi1 <- calc_pi(Y,Z,delta)
  
  const <- pi1*(1-pi1)
  
  zz.list <- split(ZZ,seq(nrow(ZZ)))
  
  Dvalue.l <- Map ('*', zz.list ,  const  )
  
  return(Dvalue.l)
  
}


derv_delta2 <- function(Y,Z,delta){
  #Z=Z1
  #
  
  ZZ <- cbind(1,Y,Z) ### create a big matrix (Nx3) 
  
  zz.list <- split(ZZ,seq(nrow(ZZ))) ##for each row of matrix, I make it as a list
  
  zz.listT <- lapply(zz.list, function(x)t(x))
  
  ZZ_L  <- Map ('%*%', zz.list, zz.listT )
  
  pi1 <- calc_pi(Y,Z,delta)
  
  const <- pi1*(1-pi1)*(1-2*pi1)
  
  value <- Map ('*', ZZ_L ,  const  )
  
  return(value)  ###The value which return the 300x3 matrix 
  
  
}



#################

der_beta2 <- function(X,Z,beta){
  #Z=Z1
  #X=X1
  #beta=beta
  int <- cbind(1,X,Z)
  
  #eta <- as.vector(int%*%beta)
  
  int.list <- split(int,seq(nrow(int)))
  
  int.listT <- lapply(int.list, function(x)t(x))
  
  XX_L  <- Map ('%*%', int.list , int.listT )
  
  const <-  rep(0, nrow(int))#mu(eta)*(1-mu(eta))*(1-2*mu(eta))
  
  value2  <-Map ('*', XX_L , const)
  
  return(value2)
}






#######   E(pi)###################

Epi  <- function(Z, delta, Y, hh, X, beta, sigma){
  #Z=Z12
  #X=X12
  #beta=beta0
  
  
  mar   <- 5*sd(Y)   ### find the marginal of the distribution (like the range of CI)
  yseq  <- seq(-mar, mar, hh)  ### chunck the margin into small width by setting hh 
  
  piExp<- NULL ### make an empty vector
  N<- length(Y) ### change based on N, will have N values
  
  id <- sort(rep(1: N, length(yseq)))  ###key point: create an index to sign each value. 
  Yseqi  <- rep(yseq, N)               #### Problem: whether they will be equal length for 
  
  Zi   <- Z[id, ]   #### Make a matrix based on the length of id 
  mui  <- cbind(1, X , Z )[id,]  %*%beta ## matrix of 
  fYi  <-  dnorm(Yseqi, mui, sigma ) # attention: must use Yseqi 
  piYi<-  calc_pi(Yseqi , Zi  ,delta)  
  piExp0  <- aggregate(piYi*fYi*hh~ id, FUN = sum)
  
  piExp0[,2]
}


#############E(pi*(Y-mu))############

EpiYmu <- function(Z, delta, Y, hh, X, beta, sigma){
  #Z=Z12
  #X=X12
  #beta=beta0  
  
  mar   <- 5*sd(Y)
  yseq  <- seq(-mar, mar, hh)
  
  piExp1<- NULL
  N<- length(Y)
  
  id <- sort(rep(1: N, length(yseq)))
  Yseqi  <- rep(yseq, N)
  
  Zi   <- Z[id, ]
  mui  <- cbind(1, X , Z )[id,]  %*%beta
  
  fYi  <-  dnorm(Yseqi, mui, sigma )
  piYi<-  calc_pi(Yseqi , Zi  ,delta)  
  piYiY <- piYi* (Yseqi- mui) 
  piExp1  <- aggregate(piYiY*fYi*hh~ id, FUN = sum)
  
  piExp1[,2]
}



EpiD  <- function(Z, delta, Y, hh, X, beta, sigma){
  #Z=Z12
  #X=X12
  #beta=beta0
  
  mar   <- 5*sd(Y)
  yseq  <- seq(-mar, mar, hh)
  
  ExpD<- NULL
  N<- length(Y)
  
  id <- sort(rep(1: N, length(yseq)))
  
  Zi   <- Z[id, ]
  mui  <- cbind(1, X , Z )[id,] %*%beta
  Xi   <- X[id, ]
  Yseqi  <- rep(yseq, N)
  fYi  <-  dnorm(Yseqi, mui, sigma )
  piYi<-  calc_pi(Yseqi, Zi  ,delta)  
  
  DY   <- der_betaM( Xi , Zi, beta)*as.vector(piYi)* fYi*hh   ### need to make the Y-mu as a vector
  
  colnames(DY  ) <- paste('D', 1:ncol(DY ), sep='')
  
  DYid  <-data.frame(DY , id=id)
  
  Dres <- as.data.frame(DYid %>% group_by(id) %>% summarize_all(sum))[,-1]
  
  Dres
}









EpiY2Ymu <- function(Z, delta, Y, hh, X, beta, sigma){
  #Z=Z12
  #X=X12
  #beta=beta0  
  
  mar   <- 5*sd(Y)
  yseq  <- seq(-mar, mar, hh)
  
  piExp1<- NULL
  N<- length(Y)
  
  id <- sort(rep(1: N, length(yseq)))
  Yseqi  <- rep(yseq, N)
  
  Zi   <- Z[id, ]
  mui  <- cbind(1, X , Z )[id,]  %*%beta
  
  fYi  <-  dnorm(Yseqi, mui, sigma )
  piYi<-  calc_pi(Yseqi , Zi  ,delta)  
  piYiY <- piYi* (Yseqi*(Yseqi- mui))/sigma 
  piExp1  <- aggregate(piYiY*fYi*hh~ id, FUN=sum)
  
  piExp1[,2]
}




############E(pi*(Y-mu)^2)###########


EpiYmu2  <- function(Z, delta, Y, hh, X, beta, sigma){
  mar   <- 5*sd(Y)
  yseq  <- seq(-mar, mar, hh)
  
  piExp2<- NULL
  N<- length(Y)
  
  id <- sort(rep(1: N, length(yseq)))
  Yseqi  <- rep(yseq, N)
  
  Zi   <- Z[id, ]
  mui  <- cbind(1, X , Z )[id,] %*%beta
  Yseqi  <- rep(yseq, N)
  fYi  <-  dnorm(Yseqi, mui, sigma )
  piYi<-  calc_pi(Yseqi, Zi  ,delta)  
  piYiY2 <- piYi* (Yseqi- mui)^2 #### E(pi*(Y-mu)^2)
  piExp2  <- aggregate(piYiY2*fYi*hh~ id, FUN=sum)
  
  piExp2[,2]
}


###########E(Dt)#######
library(dplyr)

EDt <- function(Z, delta, Y, hh, X, beta, sigma){
  #Z=Z12
  #X=X12
  #beta=beta0
  
  mar   <- 5*sd(Y)
  yseq  <- seq(-mar, mar, hh)
  
  ExpD<- NULL
  N<- length(Y)
  
  id <- sort(rep(1: N, length(yseq)))
  
  Zi   <- Z[id, ]
  mui  <- cbind(1, Y , Z )[id,]  %*%delta
  Yseqi  <- rep(yseq, N)
  fYi  <-  dnorm(Yseqi, mui, sigma ) 
  
  Dtilde <- derv_deltaM(Yseqi, Zi, delta)
  
  
  colnames(Dtilde) <- paste('D', 1:ncol(Dtilde), sep='')
  
  Dtilde  <-data.frame(Dtilde*fYi*hh, id=id)
  
  DtildeE <- as.data.frame(Dtilde %>% group_by(id) %>% summarize_all(sum))[,-1]
  
  DtildeE
}


#########E((y-mu)*Dtilde)###

EYmuDt <- function(Z, delta, Y, hh, X, beta, sigma){
  #Z=Z12
  #X=X12
  #beta=beta0
  
  mar   <- 5*sd(Y)
  yseq  <- seq(-mar, mar, hh)
  
  ExpD<- NULL
  N<- length(Y)
  
  id <- sort(rep(1: N, length(yseq)))
  
  Zi   <- Z[id, ]
  mui  <- cbind(1, Y , Z )[id,]  %*%delta
  Yseqi  <- rep(yseq, N)
  fYi  <-  dnorm(Yseqi, mui, sigma ) 
  
  DYmu <- derv_deltaM(Yseqi, Zi, delta)*as.vector(Yseqi - mui )  ### need to make the Y-mu as a vector
  
  colnames(DYmu ) <- paste('D', 1:ncol(DYmu), sep='')
  
  DYmu  <-data.frame(DYmu*fYi*hh, id=id)
  
  DtildeE2 <- as.data.frame(DYmu %>% group_by(id) %>% summarize_all(sum))[,-1]
  
  DtildeE2
}

#########E((y-mu)*DT)



EYmuDt <- function(Z, delta, Y, hh, X, beta, sigma){
  #Z=Z12
  #X=X12
  #beta=beta0
  
  mar   <- 5*sd(Y)
  yseq  <- seq(-mar, mar, hh)
  
  ExpD<- NULL
  N<- length(Y)
  
  id <- sort(rep(1: N, length(yseq)))
  
  Zi   <- Z[id, ]
  mui  <- cbind(1, Y , Z )[id,]  %*%delta
  Yseqi  <- rep(yseq, N)
  fYi  <-  dnorm(Yseqi, mui, sigma ) 
  
  DYmu <- derv_deltaM(Yseqi, Zi, delta)*as.vector(Yseqi - mui )  ### need to make the Y-mu as a vector
  
  colnames(DYmu ) <- paste('D', 1:ncol(DYmu), sep='')
  
  DYmu  <-data.frame(DYmu*fYi*hh, id=id)
  
  DtildeE2 <- as.data.frame(DYmu %>% group_by(id) %>% summarize_all(sum))[,-1]
  
  DtildeE2
}


########## E(piUci)#########
EpiUci <- function(R,Z, delta, Y, hh, X, beta, sigma){
  #Z=Z12
  #X=X12
  #beta=beta0
  
  mar   <- 5*sd(Y)
  yseq  <- seq(-mar, mar, hh)
  
  UciExp1<- NULL
  N<- length(Y)
  
  id <- sort(rep(1: N, length(yseq)))
  Yseqi  <- rep(yseq, N)
  
  Zi   <- Z[id, ]
  Xi   <- X[id, ]
  newR <- matrix(R,ncol = 1)
  Ri   <- newR[id]
  mui  <- cbind(1, X , Z )[id,]  %*%beta
  
  fYi  <-  dnorm(Yseqi, mui, sigma )
  piYi<-  calc_pi(Yseqi , Zi  ,delta)  
  
  expart1 <- EpiYmu(Z, delta, Y, hh, X, beta, sigma)
  expart2 <- Epi(Z,delta,Y,hh,X,beta,sigma)
  expi<-(expart1/ expart2)[id]    # need to add index for each fraction
  
  
  UciYi <- as.vector((Yseqi-mui)- expi)* der_betaM(Xi, Zi, beta)
  
  piUci <- piYi* (UciYi)
  
  
  colnames(piUci ) <- paste('W', 1:ncol(piUci), sep='')
  
  expUci  <-data.frame(piUci*fYi*hh, id=id)
  
  UciE1 <- as.data.frame(expUci %>% group_by(id) %>% summarize_all(sum))[,-1]
  
  #UciE2=UciE1%*%t(UciE1)
  
  UciE1
}


#############E2piUci
E2piUci <- function(R,Z, delta, Y, hh, X, beta, sigma){
  #Z=Z12
  #X=X12
  #beta=beta0
  
  mar   <- 5*sd(Y)
  yseq  <- seq(-mar, mar, hh)
  
  UciExp1<- NULL
  N<- length(Y)
  step <- length(yseq)
  
  #id <- sort(rep(1: N, length(yseq)))
  id <- rep(1:N,each=step)
  Yseqi  <- rep(yseq, N)
  
  Zi   <- Z[id, ]
  Xi   <- X[id, ]
  newR <- matrix(R,ncol = 1)
  Ri   <- newR[id]
  mui  <- cbind(1, X , Z )[id,]  %*%beta
  
  fYi  <-  dnorm(Yseqi, mui, sigma )
  piYi<-  calc_pi(Yseqi , Zi  ,delta)  
  
  expart1 <- EpiYmu(Z, delta, Y, hh, X, beta, sigma)
  expart2 <- Epi(Z,delta,Y,hh,X,beta,sigma)
  expi<-(expart1/ expart2)[id]    # need to add index for each fraction
  
  
  UciYi <- as.vector((Yseqi-mui)- expi)* der_betaM(Xi, Zi, beta)
  
  piUci <- piYi* (UciYi)
 
  lpiUci <- split(piUci,seq(nrow(piUci)))
  
  tpiUci <- lapply(lpiUci, function(x) t(x))
  
  E2 <- Map("%*%",lpiUci,tpiUci) 
  
  E2Uci <- Map("*",fYi*hh,E2)
  
  result.size <- length(E2Uci)/step
  
  split.value <- seq_along(id)
  
  splits.m <- split(split.value,id)
  
  E2piUci <- lapply(splits.m,function(x){
    Reduce("+",E2Uci[x[1]:x[length(x)]])
  })
  
  
  
  
  #UciE2=UciE1%*%t(UciE1)
  E2piUci
}


########## E(piUci^T Uci)#########
EpiUci2 <- function(R,Z, delta, Y, hh, X, beta, sigma){
  #Z=Z12
  #X=X12
  #beta=beta0
  
  mar   <- 5*sd(Y)
  yseq  <- seq(-mar, mar, hh)
  
  UciExp1<- NULL
  N<- length(Y)
  
  id <- sort(rep(1: N, length(yseq)))
  Yseqi  <- rep(yseq, N)
  
  Zi   <- Z[id, ]
  Xi   <- X[id, ]
  newR <- matrix(R,ncol = 1)
  Ri   <- newR[id,]
  mui  <- cbind(1, X , Z )[id,]  %*%beta
  
  fYi  <-  dnorm(Yseqi, mui, sigma )
  piYi<-  calc_pi(Yseqi , Zi  ,delta)  
  
  expart1 <- EpiYmu(Z, delta, Y, hh, X, beta, sigma)
  expart2 <- Epi(Z,delta,Y,hh,X,beta,sigma)
  piseq <- expart2[id]
  expi<-(expart1/ expart2)[id]
  
  
  UciYi <- as.vector((Yseqi-mui)- expi)* der_betaM(Xi, Zi, beta)
  
  piUci <- piYi* (UciYi)
  
  tUci <- split(UciYi, seq(nrow(UciYi))) ####where I changed from UciYi to piUci
  
  TUciYi <- lapply(tUci,function(x) t(x))
  
  Uci2 <- Map("%*%",tUci,TUciYi)
  
  piUci2<-Map("*",piYi*fYi*hh*Ri/piseq,Uci2)
  
  
  pop.exp <- Reduce("+",piUci2)
  
  return(pop.exp)
  
  
}



#####E(pi2Uci2)
Epi2Uci2 <- function(R,Z, delta, Y, hh, X, beta, sigma){
  #Z=Z12
  #X=X12
  #beta=beta0
  
  mar   <- 5*sd(Y)
  yseq  <- seq(-mar, mar, hh)
  step <- length(yseq)
  
  UciExp<- NULL
  N<- length(Y)
  
  #id <- sort(rep(1: N, length(yseq)))
  id <- rep(1:N, each=step)
  Yseqi  <- rep(yseq, N)
  
  Zi   <- Z[id, ]
  Xi   <- X[id, ]
  newR <- matrix(R,ncol = 1)
  Ri   <- newR[id]
  mui  <- cbind(1, X , Z )[id,]  %*%beta
  
  fYi  <-  dnorm(Yseqi, mui, sigma )
  piYi<-  calc_pi(Yseqi , Zi  ,delta)  
  
  expart1 <- EpiYmu(Z, delta, Y, hh, X, beta, sigma)
  expart2 <- Epi(Z,delta,Y,hh,X,beta,sigma)
  expi<-(expart1/ expart2)[id]    # need to add index for each fraction
  
  
  UciYi <- as.vector((Yseqi-mui)- expi)* der_betaM(Xi, Zi, beta)
  
  pi2Uci <- piYi* (UciYi)
  
  tpi2Uci <- split(pi2Uci ,seq(nrow(pi2Uci)))
  
  TUciYi <- lapply(tpi2Uci, function(x) t(x))
  
  Uci2 <- Map("%*%",tpi2Uci,TUciYi) ####got a matrix for repeated id
  
  expUci2 <- Map("*",fYi*hh,Uci2)
  
  result_size <- length(expUci2)/step
  
  split_value <- seq_along(id)
  
  splits <- split(split_value,id)
  
  Ep2Uci2 <- lapply(splits,function(x){
    Reduce("+",expUci2[x[1]:x[length(x)]])
  })
  
  
  # colnames(Uci2.d) <- paste("U",1:ncol(Uci2.d),sep="")
  # 
  # expUci2.d <- data.frame(Uci2.d*fYi*hh,id=id)
  # 
  # Epi2Uci2 <- as.data.frame(expUci2.d%>%group_by(id)%>%summarise_all(sum))[,-1]
  
  
  Ep2Uci2 
  
}



##################Another way to do 
# 
# EpiUci2 <- function(Z, delta, Y, hh, X, beta, sigma){
#   #Z=Z12
#   #X=X12
#   #beta=beta0
# 
#   mar   <- 5*sd(Y)
#   yseq  <- seq(-mar, mar, hh)
# 
#   UciExp1<- NULL
#   N<- length(Y)
# 
#   id <- sort(rep(1: N, length(yseq)))
#   Yseqi  <- rep(yseq, N)
# 
#   Zi   <- Z[id, ]
#   Xi   <- X[id, ]
#   newR <- matrix(R,ncol = 1)
#   Ri   <- newR[id,]
#   mui  <- cbind(1, X , Z )[id,]  %*%beta
# 
#   fYi  <-  dnorm(Yseqi, mui, sigma )
#   piYi<-  calc_pi(Yseqi , Zi  ,delta)
#   ciYi <- (Yseqi-mui)- EpiYmu(Z, delta, Y, hh, X, beta, sigma)/Epi(Z,delta,Y,hh,X,beta,sigma)
#   Uci <-  der_betaM(Xi, Zi, beta)*as.vector(ciYi)
# 
#   Tuci   <- t(Uci)
# 
# 
#   piUci2 <- (piYi*Tuci)%*%Uci
# 
#   colnames(piUci2 ) <- paste('W', 1:ncol(piUci2), sep='')
# 
#   expUci2  <-data.frame(piUci2*fYi*hh, id=id)
# 
#   UciE2 <- as.data.frame(expUci2 %>% group_by(id) %>% summarize_all(sum))[,-1]
# 
#   UciE2
# 
# }



##########









##################
##################
##################
##################
##################
##################
##################
##################
##################
##################
##################


normU  <- function(Y, X,Z, R, beta, delta,hh,sigma){
  # beta=beta0
  # delta=delta0
  # Z=Z12
  # X=X12
  # R=Ri
  
  int <- cbind(1,X,Z)  #int.list <- split(int,seq(nrow(int)))
  eta1 <- as.vector(int%*%beta)
  mu1 <- mu(eta1)
  
  
  e1 <- EpiYmu(Z,delta,Y,hh ,X,beta,sigma)
  e2 <- Epi(Z,delta,Y,hh ,X,beta,sigma)
  
  
  part.cont <- R*((Y-mu1)- (e1/e2))
  
  sc <- as.vector(part.cont)*der_betaM(X,Z,beta)
  
  scvalueM <- colSums(sc)
  
  
  
  #Sc <- Map("%*%", part.cont, der_beta(X,Z,beta))
  
  #scvalue <- Reduce("+",Sc)   ##### a vector contains 5 value
  
  return(scvalueM)
}




### For Uphi 

normUphi  <- function(Y, Z, R,  delta  ){
  #Z=Z1
  #R=Ri
  pi <- calc_pi(Y,Z,delta)
  
  ratio <- (R-pi)/ (pi*(1-pi))  ## a vector of constant
  
  der.valueM <- ratio*derv_deltaM(Y,Z,delta)
  
  phiM <- colSums(der.valueM)
  
  #der.value <- Map("*", ratio, derv_delta(Y,Z,delta))
  
  #phi <-  Reduce('+',der.value)
  
  return(phiM)
}


normS  <- function(Y, X,Z, R, beta, delta,hh,sigma ){
  #X <- X1
  #Z=Z1
  #R= Ri
  U1<-  normU(Y, X,Z, R, beta, delta,hh,sigma) 
  
  U2 <-   normUphi(Y, Z, R,  delta  )
  
  Svalue<-  c(U1 , U2)
  
  return(Svalue)
  
}



normdS_beta  <- function( Y, X,Z, R, beta, delta,hh,sigma){
  
  int <-  cbind(1,X,Z)
  #eta <- as.vector(int%*%beta)
  
  D <- der_betaM(X,Z,beta)
  
  
  exp1 <- Epi(Z,delta,Y,hh ,X,beta,sigma)
  exp2 <- EpiYmu(Z, delta, Y, hh , X, beta, sigma)
  exp3 <- EpiYmu2(Z, delta, Y, hh , X, beta, sigma)
  
  pii <- as.vector(1/exp1)
  pii2 <- as.vector(1/exp1^2)
  
  pR <- as.vector(R)
  
  part1 <- -(t(D*pR))%*%D
  
  part2 <- (t(D*exp1*R*pii))%*%( D)
  
  ## Note: part1 +part2 is zero for the identity link function, bt not for other links!
  
  
  pp3 <- as.vector(exp3*pii) - as.vector( (exp2*pii)^2)
  part3  <- -t(D*pp3) %*%D
  
  ## part4 is zero for the identity link function, but not otherwise. 
  part4<-0
  #Svalue <- part1+part2 
  
  Svalue <-   part1 +part2+part3 +part4
  
  ##   
  return(Svalue)
} 

normdS_delta <- function(Y, X,Z,R,beta,delta,hh,sigma){
  #Z=Z12
  #X=X12
  #R=Ri
  #beta=beta0
  #delta=delta0
  
  
  int <- cbind(1,X,Z)
  
  D <- der_betaM(X,Z,beta)
  
  exp1 <- Epi(Z,delta,Y,hh ,X,beta,sigma)
  exp2 <- EpiYmu(Z, delta, Y, hh , X, beta, sigma )
  exp3 <- EDt(Z, delta, Y, hh , X, beta, sigma )
  exp3M <- matrix(unlist(exp3),nrow(exp3),ncol(exp3) )
  exp4 <- EYmuDt(Z, delta, Y, hh , X, beta, sigma )
  exp4M <- matrix(unlist(exp4),nrow(exp4),ncol(exp4) )
  pii <- as.vector(1/exp1)
  pii2 <- as.vector(1/exp1^2)
  #exp4a <- matrix(unlist(exp4*R*pii), nrow(exp4), ncol(exp4)) ###check mode() should be numeric
  
  #part1 <- -t(D)%*%  exp4a 
  
  p1a <- as.vector(pii*R)
  part1a <- -t(D) %*% (exp4M*p1a)
  
  #exp3a <- matrix(unlist(exp3*pii^2*R*exp2), nrow(exp3), ncol(exp3))
  # part2 <- t(D)%*%exp3a
  p2b <- as.vector(pii2*R*exp2)
  part2b <- t(D)%*%(exp3M*p2b)
  
  
  #value<- part1+part2
  value <- part1a+part2b
  return(value)
}

dphi_delta <- function(R,Y,Z,delta){
  #Z=Z12
  #X=X12
  #R=Ri
  #beta=beta0
  #delta=delta0
  
  numerator <- R-calc_pi(Y,Z,delta) ## a vector 1:100
  
  denominator <- calc_pi(Y,Z,delta)*(1-calc_pi(Y,Z,delta)) ## a vector 1:100
  
  part1 <- numerator/denominator #a vector 1:100
  
  part2 <- (-1/denominator)- ((numerator*(1-2*calc_pi(Y,Z,delta)))/(denominator^2))#a vector 1:100
  
  v.derpi2 <- Map("*",part1,derv_delta2(Y,Z,delta)) # a 3x3 matrix for N
  
  derv_delta0 <-  derv_delta(Y,Z,delta) #list
  Tderv_delta0   <- lapply(derv_delta0 , function(x) t(x))
  Tmatrix  <- Map("%*%",derv_delta0, Tderv_delta0)
  v.derpi1 <- Map("*",part2,Tmatrix)
  
  
  #Tmatrix <- derv_delta0%*%t(derv_delta0)
  #v.derpi1 <- part2*Tmatrix
  
  
  Phivalue <- Reduce("+",v.derpi1) + Reduce("+",v.derpi2) ## matrix100x100
  
  return(Phivalue)
}





##################
################## Variance Component
##################
##################
##################
Uci <- function(Y, X,Z, R, beta, delta,hh,sigma){
  # beta=beta0
  # delta=delta0
  # Z=Z12
  # X=X12
  # R=Ri 
  
  int <- cbind(1,X,Z)  #int.list <- split(int,seq(nrow(int)))
  eta1 <- as.vector(int%*%beta)
  mu1 <- mu(eta1)
  
  
  e1 <- EpiYmu(Z,delta,Y,hh ,X,beta,sigma)
  e2 <- Epi(Z,delta,Y,hh ,X,beta,sigma)
  
  
  part.cont <- as.vector(Y-mu1- (e1/e2))
  
  uci <- Map("*",part.cont,der_beta(X,Z,beta))
  
  #uci <- part.cont*der_betaM(X,Z,beta)
  
  return(uci)   ###### list 5 col
} 

GUphi  <- function(Y, Z, R,  delta  ){
  # delta=delta0
  # Z=Z12
  # R=Ri 
  pi <- calc_pi(Y,Z,delta)
  
  ratio <- (R-pi)/ (pi*(1-pi))  ## a vector of constant
  
  uphi <- Map("*",ratio,derv_delta(Y,Z,delta))
  
  #uphi<- ratio*derv_deltaM(Y,Z,delta)
  
  return(uphi) ####list 4 col
} 





varbeta <- function(Y, X,Z, R, beta, delta,hh,sigma){
  
  coe1 <- as.vector(R/Epi(Z,delta,Y,hh,X,beta,sigma))
  
  Uc1M <- EpiUci2(R,Z, delta, Y, hh, X, beta, sigma)
  
  
  Uc2 <- EpiUci(R,Z, delta, Y, hh, X, beta, sigma)
  Uc2M <- matrix(unlist(Uc2),nrow(Uc2),ncol(Uc2)  )
  uc2l <- split(Uc2M,seq(nrow(Uc2M)))
  Tuc2l <- lapply(uc2l,function(x) t(x))
  
  part2a <- Map("*",coe1,uc2l)
  part2b <- Map("%*%",part2a,Tuc2l)
  part2 <- Reduce("+",part2b)
  
  
  varbeta <-  Uc1M-part2
  
  return(varbeta)
}




Covbeta <- function(Y, X,Z, R, beta, delta,hh,sigma){
  # beta=beta0
  # delta=delta0
  # Z=Z12
  # X=X12
  # R=Ri 
  
  cont <- R
  
  ucii <- Uci(Y, X,Z, R, beta, delta,hh,sigma)
  
  uphii <- GUphi(Y,Z,R,delta)
  
  
  part1 <- Map("*",cont,ucii)
  
  Tpart1 <- lapply(part1, function(x) t(x))  ####Uci^T
  
  final <- Map("%*%",uphii,Tpart1) #####Uphi*Uci^T
  
  covbeta <- Reduce("+",final)
  
  return(covbeta)
}


TCovbeta <- function(Y, X,Z, R, beta, delta,hh,sigma){
  # beta=beta0
  # delta=delta0
  # Z=Z12
  # X=X12
  # R=Ri 
  
  cont <- R
  
  ucii <- Uci(Y, X,Z, R, beta, delta,hh,sigma)
  
  uphii <- GUphi(Y,Z,R,delta)
  
  part1 <- Map("*",cont,uphii)
  
  Tpart1 <- lapply(part1, function(x) t(x)) ####Uphi^T
  
  final <- Map("%*%",ucii,Tpart1) ### #####Uci*Uphi^T
  
  Tcovbeta <- Reduce("+",final)
  
  
  return(Tcovbeta)
}


varphi <- function(Y, Z, R, X, delta, hh, beta, sigma  ){
  # delta=delta0
  # Z=Z12
  # R=Ri 
  
  cont <- R/Epi(Z,delta,Y,hh,X,beta,sigma)
  
  uphii <- GUphi(Y,Z,R,delta)
  
  tuphii <- lapply(uphii, function(x) t(x))
  
  part1 <- Map("*",cont,uphii)
  
  allmat <- Map("%*%",part1,tuphii)
  
  varphi2 <- Reduce("+",allmat)
  
  return(varphi2)
}







