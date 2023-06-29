# pi and link functions ---------------------------------------------------

calc_pi <- function(Y,Zdelta,delta){
  int <-  cbind(1,Y,Zdelta)
  #int=cbind(1,matrix(Y,ncol=1),matrix(Z, ncol = 2))
  pivalue <-  1/(1+exp(-int%*% t(t(delta))))
  
  return(pivalue)
} 



hlink <- function(b,X,Z,M,beta){
  
  fix.int <- cbind(1,X,Z)
  #fix.int = cbind(1, matrix(X, ncol = 2), matrix(Z, ncol = 2))
  link <- 1/(1+exp(- ((fix.int)%*%t(t(beta))) - M%*%t(t(b))))
  
  return(link)
}


calc_pi2 <- function(Y,Zdelta,delta){
  
  int=cbind(1,matrix(Y,ncol=1),matrix(Zdelta, ncol = 2))
  pivalue <-  1/(1+exp(-int%*% t(t(delta))))
  
  return(pivalue)
} 

#### for each scalar of bi, ONLY USE FOR EXPECTAION !!!!
hlink2 <- function(b,X,Z, beta){
  
  fix.int = cbind(1, matrix(X, ncol = 2), matrix(Z, ncol = 2))
  link <- 1/(1+exp(- ((fix.int)%*%t(t(beta))) -  b))
  
  return(link)
}


# second derivatives ------------------------------------------------------


sec.order.Q.numer <- function(b,X,Z,Y,beta,idks,sigma){
  
  #D <- -solve(diag(sigma, ncol=length(b),nrow=length(b)))
  
  nameK  <- unique(idks)
  
  frac.n2 <- NULL
  for (i in nameK) {
    #i=1
    K  <- idks== i
    
    fix <- cbind(1,X[K,],Z[K,])
    
    exp1 <- exp(  fix%*% beta + b[i])
    
    frac1 <- -(Y[K]*exp1)/((1+exp1)^2)
    
    frac2 <- ((Y[K]-1)*exp1)/((1+exp1)^2)
    
    
    frac.n2 [i]<- sum(frac1+frac2)
    
  }
  #sec.out <- D + frac
  
  sec.out <- as.vector(-sigma^(-2)+frac.n2)
  
  return(diag(sec.out))
  
}



# I total function --------------------------------------------------------------

Q.total.numer <- function(b,listM){
  
  
  #listM=list(p=p,d=d,thetaStart=c(theta,deltaFixed),idks=idks, X=X.S,Z=Z.S,Zdelta=Zdelta,M=M.S,R=R.S,Y=Y.S)
  
  listRun <- TRUE
  if(listRun){
    p <- listM$p
    d <- listM$d
    thetaStart <- listM$thetaStart
    sigmaStart <-  thetaStart[1] 
    betaStart <-thetaStart[2:(p+1)] 
    deltaStart<- thetaStart[(p+2):(1+p+d)]
    
    
    idks<-listM$idks
    X<- listM$X
    Z<- listM$Z
    Zdelta <- listM$Zdelta
    M<- listM$M
    K <- length(unique(idks))
    R<- listM$R
    Y<- listM$Y
    
  }
  
  # sigmaStart  <- sqrt(var(b))
  
  D.density <- mvtnorm::dmvnorm(b,mean = rep(0,length(b)),sigma = diag((sigmaStart^2), ncol=length(b),nrow=length(b)),log = TRUE)
  
  nameK  <- unique(idks)
  Q.numer  <- NULL
  
  
  for(i in nameK){
    
    condK  <- idks== i
    
    #condK <- idks == k  ## needed to be in the list
    
    #b[i] <- bk
    
    
    
    hk <-  hlink(b,X[condK,],Z[condK,],M[condK,],betaStart)
    
    pik <- calc_pi(Y[condK],Zdelta[condK,],deltaStart)
    
    
    part1 <-   log(as.vector(pik))
    
    
    logf.yi <-   Y[condK]*log(hk) + (1-Y[condK])*log(1-hk)
    
    
    Q.numer[i] <- sum(part1 +logf.yi)
    
    
    
  }
  
  
  Q.num.sum <- sum(Q.numer)+D.density
  
  
  return( Q.num.sum)
  
  
}




I.Total.numer3 <- function(logsigma,listC){ ##  k, b
  
  # listC= list(p=p,d=d,bnum=b,idks=idks,X=X.S,M=M.S,Z=Z.S,Zdelta=Zdelta,R=R.S,Y=Y.S, deltaFixed=deltaFixed,betaFixed=betaFixed)
  
  
  listRun <- TRUE
  if(listRun){
    p <- listC$p
    d <- listC$d 
    
    sigma <- exp(logsigma)
    beta <-listC$betaFixed 
    delta<- listC$deltaFixed
    
    
    bnum <- listC$bnum
    k <-  listC$k
    idks<- listC$idks 
    
    X<- listC$X
    Z<- listC$Z
    Zdelta <- listC$Zdelta
    M<- listC$M
    
    R<- listC$R
    Y<- listC$Y
    K <- length(unique(idks))
  }
  
  
  
  
  Q.bhat.k <-    Q.total.numer(b=bnum ,listM= list(p=p,d=d,thetaStart=c(sigma,beta,delta),idks=idks, X=X,Z=Z,Zdelta=Zdelta,M=M,R=R,Y=Y))
  Q.hess.n <-   sec.order.Q.numer(b=bnum,X=X,Z=Z,Y=Y,beta=beta,idks=idks,sigma=sigma)
  
  q <- nrow(Q.hess.n)
  #Ik <-  base::((2*pi)^(q/2))*((-det(Q.hess.n))^(-1/2))*exp(Q.bhat.k)) #### given q-dimension of b
  logIk <-   base::log(((2*pi)^(q/2))) +  (-1/2) *log(det(-Q.hess.n)) + Q.bhat.k#### given q=2)
  
  
  return(  logIk )
  
  #return(  logIk  )
  
}


I.Total.numer4 <- function(beta,listC){ ##  k, b
  
  # listC= list(p=p,d=d,bnum=b ,idks=idks,X=X.S,M=M.S,Z=Z.S,R=R.S,Y=Y.S, deltaFixed=deltaFixed,sigmaFixed=sigma)
  
  
  listRun <- TRUE
  if(listRun){
    p <- listC$p
    d <- listC$d 
    
    sigma <- listC$sigmaFixed
    #beta <-beta 
    delta<- listC$deltaFixed
    
    
    bnum <- listC$bnum
    k <-  listC$k
    idks<- listC$idks 
    
    X<- listC$X
    Z<- listC$Z
    Zdelta <- listC$Zdelta
    M<- listC$M
    
    R<- listC$R
    Y<- listC$Y
    K <- length(unique(idks))
  }
  
  
  
  
  Q.bhat.k <-    Q.total.numer(b=bnum ,listM= list(p=p,d=d,thetaStart=c(sigma,beta,delta),idks=idks, X=X,Z=Z,Zdelta=Zdelta,M=M,R=R,Y=Y))
  Q.hess.n <-   sec.order.Q.numer(b=bnum,X=X,Z=Z,Y=Y,beta=beta,idks=idks,sigma=sigma)
  
  q <- nrow(Q.hess.n)
  #Ik <-  base::((2*pi)^(q/2))*((-det(Q.hess.n))^(-1/2))*exp(Q.bhat.k)) #### given q-dimension of b
  logIk <-   base::log(((2*pi)^(q/2))) +  (-1/2) *log(det(-Q.hess.n)) + Q.bhat.k#### given q=2)
  
  
  return(  logIk )
  
  #return(  logIk  )
  
}



I.Total.numdenomEff2<- function(theta,listC){ ##  k,
  
  # listC= list(p=p,d=d,  idks=idks,X=X.S,M=M.S,Z=Z.S,Zdelta=Zdelta.S,R=R.S,Y=Y.S,deltaFixed=deltaFixed)
  
  
  listRun <- TRUE
  if(listRun){
    p <- listC$p
    d <- listC$d
    
    logsigma <- theta[1]
    beta <-theta[2:(p+1)]
    delta<- theta[(p+2):(1+p+d)]
    
    sigma <- exp(logsigma)
    
    
    k <-  listC$k
    idks<- listC$idks
    
    X<- listC$X
    Z<- listC$Z
    Zdelta <- listC$Zdelta
    M<- listC$M
    
    R<- listC$R
    Y<- listC$Y
    
  }
  
  sigmaAux  <-  2
  
  
  
  
  nameK  <- unique(idks)
  alldat   <-  data.frame(id=idks, X=X,Z=Z,Y=Y)
  groupdat <-  split(alldat,idks)
  I.denom<- NULL
  result.p0 <- list()
  result.p1 <- list()
  HK  <- list()
  HKn <- list()
  
  Exp <- NULL
  
  hh <- 0.5
  
  bran <- cbind(seq(-5*sigmaAux, 5*sigmaAux, hh))
  
  
  HK  <- apply(bran, 1,function(x) hlink2(x,as.matrix(alldat[,2:3]),as.matrix(alldat[,4:5]), beta))  
  
  fbk  <-  as.vector(dnorm(bran, 0, sd= sigma))
  
  
  HKk0  <- t(t(HK )*fbk*hh )
  
  Expk <-   (( rowSums(   HKk0 )))
  mu1 <- Expk
  mu0 <- 1-mu1
  pi1 <-  calc_pi2(Y*0+1,Zdelta,delta)
  pi0 <-  calc_pi2(Y*0,Zdelta,delta)
  
  
  piE <- mu1*pi1 +mu0*pi0
  
  # Ym  <- sapply(1:4000, function(x) rbinom(length(Expk), 1, Expk))
  
  # piY <- apply(Ym, 2,  function(x)  calc_pi2(x,Z,delta))
  
  
  #piE <-  rowMeans(piY)
  
  sum(log(piE))
  
  
  
}



I.Total.numer <- function(theta,listC){ ##  k, b
  
  # listC= list(p=p,d=d,b=bi,idks=idks,X=X.S,M=M.S,Z=Z.S,R=R.S,Y=Y.S,deltaFixed=deltaFixed)
  
  listRun <- TRUE
  if(listRun){
    p <- listC$p
    d <- listC$d 
    
    
    sigma <- theta[1]
    beta <-theta[2:(p+1)] 
    delta<- listC$deltaFixed
    
    
    b <- listC$b
    k <-  listC$k
    idks<- listC$idks 
    
    X<- listC$X
    Z<- listC$Z
    Zdelta <- listC$Zdelta
    M<- listC$M
    
    R<- listC$R
    Y<- listC$Y
    K <- length(unique(idks))
  }
  
  
  
  
  
  #  K.numer  <- optim(par=c(b=b ),fn = Q.total.numer , gr = NULL, 
  #                   listM= list(p=p,d=d,thetaStart=theta,idks=idks,b=b,X=X,Z=Z,M=M,R=R,Y=Y), 
  #                   method = "BFGS", hessian = TRUE,control= list(fnscale="-1"))
  
  
  
  Q.bhat.k <-    Q.total.numer(b ,listM= list(p=p,d=d,thetaStart=c(sigma,beta,delta),idks=idks, X=X,Z=Z,Zdelta=Zdelta,M=M,R=R,Y=Y))
  Q.hess.n <-   sec.order.Q.numer(b=b,X=X,Z=Z,Y=Y,beta=beta,idks=idks,sigma=sigma)
  
  q <- nrow(Q.hess.n)
  #Ik <-  base::((2*pi)^(q/2))*((-det(Q.hess.n))^(-1/2))*exp(Q.bhat.k)) #### given q-dimension of b
  logIk <-   base::log(((2*pi)^(q/2))) +  (-1/2) *log(det(-Q.hess.n)) + Q.bhat.k#### given q=2)
  
  
  return(  logIk )
  
  #return(  logIk  )
  
}




I_Total.theta <- function(theta,listC){ # Not Y
  
  
  listRun <- TRUE
  if(listRun){
    p <- listC$p
    d <- listC$d
    
    logsigma <- theta[1]
    sigma <- exp(logsigma)
    
    beta <-theta[-1]
    delta<- listC$deltaFixed
    
    
    bnum<- listC$b
    idks<- listC$idks
    K <- length(unique(idks))
    Y<- listC$Y
    X<- listC$X
    Z<- listC$Z
    Zdelta <- listC$Zdelta
    M<- listC$M
    
    R<- listC$R
    
    
    
  }
  
  
  
  #Aux.numer  <- I.Total.numer4(logsigma ,listC=list(p=p,d=d, betaFixed=beta, deltaFixed=delta, X=X, bnum=bnum,Z=Z  ,Y=Y  , R=R , M=M  ,idks=idks))
  Aux.numer  <- I.Total.numer(theta=c(sigma,beta) ,listC=list(p=p,d=d, deltaFixed=delta, X=X, bnum=bnum,Z=Z,Zdelta=Zdelta ,Y=Y,
                                                              R=R , M=M  ,idks=idks))
  
  
  Aux.numerD <- I.Total.numdenomEff2(theta=c(logsigma,beta, delta) ,listC=list(p=p,d=d,  X=X  , bnum=bnum ,Z=Z,Zdelta=Zdelta,Y=Y, R=R , M=M  ,idks=idks  )) 
  
  
  Aux  <-  Aux.numer  -  Aux.numerD 
  
  
  
  # final <- list( I.Total= Aux/K , bhat_n=bhat_n,  bhat_d=bhat_d) #### Aux0 contains 0
  final <-   Aux
  return(final)
}


# (need to check)New version of combined function  ------------------------



Q.total.numerbsigma <- function(bsigma,listM){
  
  
  #listM=list(p=p,d=d,thetaStart=c(1,betaFixed,deltaFixed),idks=idks, X=X.S,Z=Z.S,Zdelta=Zdelta,M=M.S,R=R.S,Y=Y.S)
  
  listRun <- TRUE
  if(listRun){
    b <- bsigma[-1]
    
    p <- listM$p
    d <- listM$d
    thetaStart <- listM$thetaStart
    sigmaStart <-  thetaStart[1] 
    betaStart <-thetaStart[2:(p+1)] 
    deltaStart<- thetaStart[(p+2):(1+p+d)]
    
    logsigma  <- bsigma[1]
    
    idks<-listM$idks
    X<- listM$X
    Z<- listM$Z
    Zdelta <- listM$Zdelta
    M<- listM$M
    K <- length(unique(idks))
    R<- listM$R
    Y<- listM$Y
    
  }
  
  
  sigmaStart  <- exp(logsigma)
  
  D.density <- mvtnorm::dmvnorm(b,mean = rep(0,length(b)),sigma = diag((sigmaStart^2), ncol=length(b),nrow=length(b)),log = TRUE)
  
  nameK  <- unique(idks)
  Q.numer  <- NULL
  
  
  for(i in nameK){
    
    condK  <- idks== i
    
    #condK <- idks == k  ## needed to be in the list
    
    #b[i] <- bk
    
    hk <-  hlink(b,X[condK,],Z[condK,],M[condK,],betaStart)
    
    pik <- calc_pi(Y[condK],Zdelta[condK,],deltaStart)
    
    
    part1 <-   log(as.vector(pik))
    
    
    logf.yi <-   Y[condK]*log(hk) + (1-Y[condK])*log(1-hk)
    
    
    Q.numer[i] <- sum(part1 +logf.yi)
    
    
    
  }
  
  
  Q.num.sum <- sum(Q.numer)+D.density
  
  
  return( Q.num.sum)
  
  
}


I_Total.beta_sig_nested_b <- function(beta,listC){ # Not Y
  
  #listC=list(p=length(betaOld),d=length(deltaFixed),bnum= binumold, idks=idks,X=X,Z=Z,Zdelta=Zdelta,M=M,R=R,Y=Y,deltaFixed=deltaFixed, logsigma=logsigma )
  listRun <- TRUE
  if(listRun){
    p <- listC$p
    d <- listC$d
    
    logsigma <- listC$logsigma
    delta<- listC$deltaFixed
    
    
    bnum <- listC$bnum
    idks<- listC$idks
    K <- length(unique(idks))
    Y<- listC$Y
    X<- listC$X
    Z<- listC$Z
    Zdelta <- listC$Zdelta
    M<- listC$M
    R<- listC$R
    
    
    sbnum<-    optim(par = c(logsigma, bnum  ),fn=Q.total.numerbsigma, listM=list(p=p,d=d,thetaStart=c(1,beta, delta),idks=idks,X=X,Z=Z,Zdelta=Zdelta,M=M,R=R,Y=Y),
                     control=list(fnscale=-1))$par
    
    
    logsigma <- sbnum[1]
    sigma <- exp(logsigma)
    bnum <- sbnum[-1]
    
    print(sigma)
    # print(beta)
    #logsigma <- log(sigma)
    
  }
  
  
  
  Aux.numer  <- I.Total.numer4(beta ,listC=list(p=p,d=d,sigmaFixed=sigma, deltaFixed=delta,
                                                bnum=bnum,X=X,Z=Z,Zdelta=Zdelta  ,Y=Y ,  R=R , M=M  ,idks=idks))
  
  Aux.numerD <-I.Total.numdenomEff2(c(logsigma,beta, delta) ,listC=list(p=p,d=d, bnum=bnum ,X=X  ,Z=Z,Zdelta=Zdelta,
                                                                        Y=Y ,R=R , M=M  ,idks=idks  ))
  
  
  Aux  <-  Aux.numer  -  Aux.numerD 
  
  # final <- list( I.Total= Aux/K , bhat_n=bhat_n,  bhat_d=bhat_d) #### Aux0 contains 0
  print(beta)
  final <-   Aux/K
  return(final)
  
}








# NR algorithm ------------------------------------------------------------
library(lme4)
library(numDeriv)
NR.full3<- function(thetaStart, bStart,R_F, Y_F,Zdelta_F, Z_F,Y, X, R,Z,Zdelta, M,p,d,idks )  ## does not inlcude beta or b
{
  # thetaStart=theta; bStart=b;R_F = R; Y_F=Y; Zdelta_F=Zdelta; Z_F=Z;Y=Y.S; X=X.S; Z=Z.S; Zdelta=Zdelta.S; M=M.S;R=R.S;p=p;d=d;idks=idks
  
  deltaFixed <- glm(R_F~ Y_F+Zdelta_F, family= "binomial")$coef
  
  auxmm  <-  glmer(Y~X+Z +(1|idks),  family="binomial" )
  KK     <-  1
  
  betaFixed  <- as.vector(fixef(auxmm)) 
  bFixed     <- as.vector(unlist(ranef(auxmm)))*KK
  sigmaFixed <- sqrt(unlist(VarCorr(auxmm)))
  
  
  sigmaOld <-sigmaFixed
  logsigma <- log(sigmaOld)
  betaOld<- betaFixed
  binumold <- b0 <- bFixed 
  
  
  ### only optim function is used
  
  obj <- optim(par = c(betaOld ),fn=I_Total.beta_sig_nested_b, listC=list(p=p,d=d,bnum= binumold, idks=idks,X=X,Z=Z,Zdelta=Zdelta,M=M,R=R,Y=Y,deltaFixed=deltaFixed, logsigma=logsigma ),
               control=list(fnscale=-1, reltol=0.001)  ) 
  # 
  betai <- obj$par
  
  sbnum<-    optim(par = c(logsigma, binumold  ),fn=Q.total.numerbsigma, listM=list(p=p,d=d,thetaStart=c(1,betai, deltaFixed),idks=idks,X=X,Z=Z,Zdelta=Zdelta,M=M,R=R,Y=Y),
                   control=list(fnscale=-1))$par
  
  
  logsigma.est <- sbnum[1]
  sigma.est <- exp(logsigma.est)
  bi <- sbnum[-1] 
  names(bi) <- unique(idks)
  
  
  
  print(var(bi))
  print(sigma.est^2)
  
  
  
  ####calculate the SD by emprical method
  thetahat <- c(logsigma.est,betai)
  listC= list(p=p,d=d,b=bi,idks=idks,X=X.S,M=M.S,Z=Z.S,Zdelta=Zdelta.S,R=R.S,Y=Y.S,deltaFixed=deltaFixed)
  Hess    <- numDeriv::hessian(I_Total.theta,thetahat,listC=listC)
  SD.beta <- sqrt(diag(-solve(Hess))) 
  
  print(SD.beta)
  #names(bi) <- 1:K
  
  return(list(thetaFinal=c(sigma.est,betai), bFinal=bi,SD.beta=SD.beta))
}



# Generate delta and beta automatically -----------------------------------

selectdelta0 <-function(   Z,Zdelta, X, M,  nTarget, delta, beta, sigma,K){
  seqi <- seq(-5, 5, 0.05)
  deltaM  <- matrix(delta, length(seqi) ,length(delta) , byrow=TRUE)
  deltaM[,1] <- deltaM[,1]  + seqi
  
  n_prop <- NULL
  
  bM <-  matrix((rnorm(K*50,0,sigma)), 50, K)
  fY<-  apply( bM ,1 , function(x){hlink(x, Z, X,M, beta )})
  YM <- matrix(rbinom(nrow(fY)*50, 1, fY), nrow(fY), 50)
  
  
  for(j in 1: nrow(deltaM)){
    pii <- rowMeans(apply(YM, 2,function(x) calc_pi(x,Zdelta,deltaM[j,])))
    n_prop[j] <- sum( pii)
  }
  
  cond <-  abs(n_prop - nTarget)
  alp0 <-   deltaM[cond==min(cond), ]
  
  list(alp0=alp0, n_prop= n_prop[cond==min(cond)])
}





selectABeta0 <-function(Z, X,M,   prevTarget, beta, sigma,K){
  K <- ncol(M)
  seqi <- seq(-5, 5, 0.05)
  
  betaM  <- matrix(beta, length(seqi) ,length(beta) , byrow=TRUE)
  betaM[,1] <- betaM[,1]  + seqi
  
  bM <-  matrix((rnorm(K*30,0,sigma)), 30, K)
  
  prop <- NULL
  for(j in 1: nrow(betaM)){
    
    probY<- rowMeans(apply( bM ,1 , function(x){hlink(x, Z, X,M, betaM[j,])}))
    prop[j] <- mean( probY)
  }
  
  cond <-  abs(prop - prevTarget)
  bet0 <-   betaM[cond==min(cond), ]
  
  list(bet0=bet0, prop= prop[cond==min(cond)])
}



#--------------iteration start----------------
IT <-2
N <- 5000



Kv      <- c(50)  ## multiples of N
SampleSizesv <- c(1000)
prev_v <- c(0.2 )
sigmauv      <- c(2)
Samplingv <-  c(-2,-1,-1,-1) 
link   <- 'logit'
family <- 'binomial'

beta   <-  c(-3, 1,  1,  1,  1)  
delta  <-  c(-3,1,1,1) 
p      <-  5
d      <-  4
#sigma  <-  0


CountSc <-0
CountErrors <- 0

Samplesize <- NULL
REAL.theta <- matrix(NA, nrow =  (IT), ncol =  (1+p+d))
rownames(REAL.theta) <- 1:IT

OUT.theta <- matrix(NA, nrow =  (IT), ncol =  (1+p))
rownames(OUT.theta) <- 1:IT
OUT.b <- list()
OUT.SD.beta <- matrix(NA, nrow =  (IT), ncol =  (1+p))

lme   <- matrix(NA, nrow =  (IT), ncol =  p)
lme.W1<- matrix(NA, nrow =  (IT), ncol =  p)
lme.W2<- matrix(NA, nrow =  (IT), ncol =  p)

lme.b <- lme.W1b <- lme.W2b <-list()
SD.lme <- SD.lme.W1 <- SD.lme.W2 <- list()
SE.lme <- SE.lme.W1 <- SE.lme.W2 <- matrix(NA, nrow=(IT), ncol= p)

library(tidyverse)
library(Matrix)

Z1  <- rbinom(N,1,0.4)   #Z1  <- rnorm(N,0,1) ####a population
Z2  <- rnorm(N, 0,1)
Z3  <- rnorm(N, 0,1)
vx1 <- -1
px1 <- 1/(1+exp(-(vx1)-Z1)) ;mean(px1)
vx2 <- -1
px2 <- 1/(1+exp(-(vx2)-Z2)) ;mean(px2)#px2 <- 1/(1+exp(-1-(vx2)-Z2)) ;mean(px2)
X1  <- rbinom(N, 1, px1) 
X2  <- rbinom(N, 1, px2)

Z    <- cbind(Z1, Z2)
Zdelta    <- cbind(Z1, Z3)
X   <- cbind(X1, X2)


idK <- sort(rep(1:Kv,N/Kv))
Model<- gl(Kv,N/Kv)
M <-as.matrix(t(as(Model, Class = "sparseMatrix")))

#SampleSizesvi=1
nTarget <- SampleSizesv
delta0  <- selectdelta0( Z=Z, Zdelta=Zdelta, X=X, M=M,  nTarget=nTarget, delta=delta, beta=beta, sigma=sigmauv ,K=Kv) 
alp0    <- delta0$alp0



##  ## save delta and n_prop
#prev_vi=1
prevTarget <- prev_v
beta0a <- selectABeta0(Z=Z, X=X,M=M,   prevTarget=prevTarget, beta=beta, sigma=sigmauv,K=Kv)
beta0 <- beta0a$bet0


i <- 0
while(i < IT){
  i <- i+1
  # generate Y with beta0
  #sigmauvi=1
  b <- as.vector(rnorm(Kv,0,sigmauv))
  names(b) <- 1:Kv
  
  eta   <- as.vector(cbind(1,X,Z)%*%beta0 + M%*%b)
  mu     <- function(eta)  1/(1+exp(-eta))  #logit link eta 
  
  
  Y <- rbinom(N,1, mu(eta))
  mean(Y)
  
  
  
  # save sample size- sum of pi s
  pi.ind <- calc_pi(Y=Y,Zdelta=Zdelta,delta=alp0)
  mean(pi.ind)
  
  R    <- rbinom(N, 1, pi.ind);
  samplesizer <- sum(R)
  #names(samplesizer) <- 1:Kv## real sample size
  
  idks <- idK[R==1]
  colnames<-unique(idks)
  ID=unique(idks)
  
  X.S <- cbind(idK,X)[idK%in%ID & R==1,][,-1]
  Z.S <- cbind(idK,Z)[idK%in%ID& R==1,][,-1]
  Zdelta.S <- cbind(idK,Zdelta)[idK%in%ID& R==1,][,-1]
  Y.S <- cbind(idK,Y)[idK%in%ID& R==1,][,-1]
  R.S <- cbind(idK,R)[idK%in%ID& R==1,][,-1]
  Pi.S <- cbind(idK,pi.ind)[idK%in%ID& R==1,][,-1]
  M.S <- M[idK%in% ID& R==1,ID]
  b.S <- b[ID]
  
  theta <- c(sigmauv,beta0,alp0)
  
  #-----------------Add weights for comparsion----------------
  dataF <- data.frame(idK=idK, pi.ind=pi.ind,  R=R)
  data_clusterF <- dataF %>%group_by(idK)%>%summarise(pi.k= 1-prod(1-pi.ind))
  data_clusterF <-left_join(dataF,data_clusterF , 'idK' )   %>%mutate(pi.ind.cond=pi.ind/pi.k, wki= 1/pi.ind.cond) 
  
  dataS <- data_clusterF[R==1,]
  
  data_cluster <- dataS %>%group_by(idK)%>%summarise(wk= sum(1/pi.ind.cond), wk2= sum(1/pi.ind.cond^2), nk= n())
  
  dataS <- left_join(dataS,data_cluster , 'idK' )  %>%mutate(wki_a=wki*nk/wk , wki_b=wki*wk/wk2, 
                                                             wki_aF=wki_a*(1/pi.k) , wki_bF=wki_b*(1/pi.k)  )
  
  #-------------------------------------------------------------------
  ### sample units in cluster k
  
  
  
  
  
  
  # apply method ------------------------------------------------------------
  
  theta <- c(sigmauv,beta0,alp0)
  #Newton[[i]] <- NR.mixed(thetaStart = theta,bStart = b,Z_F=Z,Y_F=Y,R_F = R,X=X.S,Z=Z.S,Y=Y.S,R=R.S,M=M.S,idks=idks,p=p,d=d)
  
  Newton.denom <- tryCatch(NR.full3(thetaStart = theta,bStart = b,R_F = R,Y_F = Y,Zdelta_F=Zdelta,Z_F = Z,Y=Y.S,X=X.S,R=R.S,Z=Z.S,Zdelta = Zdelta.S,M=M.S,p=p,d=d,idks = idks), error=function(e) 'Skip')
  #system.time(Newton.denom <-NR.full(thetaStart = theta,bStart = b,Z_F=Z,Y_F=Y,R_F = R,X=X.S,Z=Z.S,Y=Y.S,R=R.S,M=M.S,idks=idks,p=p,d=d))
  #system.time(NR.full2(thetaStart = theta,bStart = b,R_F = R,Y_F = Y,Z_F = Z,Y=Y.S,X=X.S,R=R.S,Z=Z.S,Zdelta = Zdelta.S,M=M.S,p=p,d=d,idks = idks))
  if(Newton.denom[1]=='Skip')
  {
    i=i-1
    CountErrors <- CountErrors+1
  }
  
  if (Newton.denom[1]!='Skip'){
    
    Samplesize[i] <- samplesizer
    REAL.theta[i,] <- theta
    
    OUT.theta[i,] <- Newton.denom$thetaFinal
    OUT.b[[i]]    <- Newton.denom$bFinal
    OUT.SD.beta[i,] <- Newton.denom$SD.beta
    
    #-----------glmer------------------------
    LMM.1 <-  glmer(Y.S~X.S+Z.S +(1|idks), family="binomial")
    LMM.2 <-  glmer(Y.S~X.S+Z.S +(1|idks), family="binomial",weights = dataS$wki_aF) 
    LMM.3 <-  glmer(Y.S~X.S+Z.S +(1|idks), family="binomial",weights =dataS$wki_bF )
    
    #------------standard error for beta in glmer------------
    SE.lme[i,]    <- coefficients(summary(LMM.1))[,2]
    SE.lme.W1[i,] <-coefficients(summary(LMM.2))[,2]
    SE.lme.W2[i,] <-coefficients(summary(LMM.3))[,2]
    
    #-----------Standard deviation for glmer------------------
    SD.lme[[i]]    <- VarCorr(LMM.1)
    SD.lme.W1[[i]] <- VarCorr(LMM.2)
    SD.lme.W2[[i]] <- VarCorr(LMM.3)
    
    #-------fixed effect-----------------
    lme[i,]    <- as.vector(fixef(LMM.1 ))###maybe an array will be much easier to compare in the end.
    lme.W1[i,] <- as.vector(fixef(LMM.2 ))
    lme.W2[i,] <- as.vector(fixef(LMM.3))
    
    #-------random effect
    lme.b[[i]]   <- unlist(ranef(LMM.1) )
    lme.W1b[[i]] <- unlist(ranef(LMM.2) )
    lme.W2b[[i]] <- unlist(ranef(LMM.3)) 
    
    
    CountSc  <- CountSc +1
    print(paste("CountSC=",CountSc))
  }
  
  #CountSc  <- CountSc +1
  #print(countSc)
  print(paste("CountErrors=",CountErrors))
  #return(OUT)
  print(paste("iteration=",i))
}

Results <- list(Sample.size=Samplesize, True.theta=REAL.theta,NR.theta=OUT.theta,NR.b=OUT.b,NR.SDbeta=OUT.SD.beta,lme.beta=lme,lme2.beta=lme.W1,lme3.beta=lme.W2,lme.b=lme.b,lme2.b=lme.W1b,lme3.b=lme.W2b, 
                SE.lme=SE.lme, SE.lme.W1=SE.lme.W1,SE.lme.W2=SE.lme.W2,SD.lme=SD.lme, SD.lme.W1=SD.lme.W1,SD.lme.W2=SD.lme.W2 )


saveRDS(Results, paste('~/Laplace_approximation/output','15-09-2022_K',Kv ,'_N',N,'_SampleSizesv',SampleSizesv,'_sigmauv',sigmauv,'_prevlance', prev_v, '_IT', IT , '_XY', '.rds', sep=""))

#saveRDS(Results, paste('/Volumes/GoogleDrive/My Drive/Approximation method/Laplacian output/all_','07-09-2022_K',Kv ,'_N',N,'_samplesize',SampleSizesv,'_sigma',sigmauv, '_prevlance', prev_v, '_steps', IT , '_XY',  '.rds', sep=""))
# 
# 
# out <- readRDS("/Volumes/GoogleDrive/My Drive/Approximation method/Laplacian output/numerator_only_02-07-2022_K20_N20000_samplesize2100_prevlance0.2_steps60_XY.rds")
# 
# 
# 

# 
# 
# 
Out <- readRDS("/Users/xyi504/Desktop/output15-09-2022_K50_N5000_SampleSizesv600_sigmauv1_prevlance0.2_IT500_XY.rds")

Out <- readRDS("/Users/xyi504/Desktop/output16-09-2022_K50_N5000_SampleSizesv1000_sigmauv1_prevlance0.2_IT500_XY.rds")

# 
#   upper_NR <- colMeans(Out$NR.theta[,2:6])+1.96*sqrt(diag(var(Out$NR.theta[,2:6])))
#   lower_NR <- colMeans(Out$NR.theta[,2:6])-1.96*sqrt(diag(var(Out$NR.theta[,2:6])))
# 
#   upper_lm <- colMeans(Out$lme.beta)+1.96*sqrt(diag(var(Out$lme.beta)))
#   lower_lm <- colMeans(Out$lme.beta)-1.96*sqrt(diag(var(Out$lme.beta)))
# 


# Summary statistics (first draft; idea) ----------------------------------


mean(Out2$Sample.size)
# head(unlist(Out$SD.lme.W2)
blah <- rbind(colMeans(Out$NR.theta[,2:6] ),colMeans(Out$lme.beta),colMeans(Out$lme2.beta),colMeans(Out$lme3.beta))
# CI <- cbind(c(lower_W1SE[1],upper_W1SE[1]),c(lower_NR[4],upper_NR[4]))
# library(xtable)
# xtable(blah)
# Out$NR.b

#### Empirical SE
lme1 <- cbind(unlist(Out$SD.lme),Out$lme.beta )
lme2 <- cbind(unlist(Out$SD.lme.W1),Out$lme2.beta )
lme3 <- cbind(unlist(Out$SD.lme.W2),Out$lme3.beta )

beta.NR <- sqrt(diag(var(Out$NR.theta))) ### the standard error for beta, need to compare with glmer
beta.lm <- sqrt(diag(var(lme1))) ### the standard error for beta, need to compare with glmer
beta.W1 <- sqrt(diag(var(lme2))) ### the standard error for beta, need to compare with glmer
beta.W2 <- sqrt(diag(var(lme2))) ### the standard error for beta, need to compare with glmer

Ep.SE <- rbind(beta.NR, beta.lm, beta.W1, beta.W2)
xtable(Ep.SE)
### standard error which comes from the summary table of the glmer

se1 <- unlist(Out$SD.lme)/sqrt(Out$Sample.size)
se.lme <- cbind(se1, Out$SE.lme)

se2 <- unlist(Out$SD.lme.W1)/sqrt(Out$Sample.size)
se.lme2 <- cbind(se2, Out$SE.lme.W1)

se3 <- unlist(Out$SD.lme.W2)/sqrt(Out$Sample.size)
se.lme3 <- cbind(se3, Out$SE.lme.W2)

SE.NR <- colMeans(Out$NR.SDbeta)
SE.lme <- colMeans(se.lme)
SE.W1 <- colMeans(se.lme2)
SE.W2 <- colMeans(se.lme3)

Et.SE <- rbind(SE.NR, SE.lme, SE.W1, SE.W2)
xtable(Et.SE)

######
xtable(Ep.SE/Et.SE)

#### Coverage

Beta.tr <- Out$True.theta[,1:6]
Beta.es <- Out$NR.theta
Beta.ep <- Out$NR.SDbeta

Coverage.NR <- matrix(NA,nrow = 1000, ncol = 6)
Coverage.lme <- matrix(NA,nrow = 500, ncol = 5)
Coverage.W1 <- matrix(NA,nrow = 500, ncol = 5)
Coverage.W2 <- matrix(NA,nrow = 500, ncol = 5)
for (i in 1:1000) {
  for (j in 1:6){
    upper_NR <- Beta.es[i,j] +1.96*Beta.ep[i,j]
    lower_NR <- Beta.es[i,j] -1.96*Beta.ep[i,j]
    
    # upper_lme <- Out$lme.beta[i,j] + 1.96*Out$SE.lme[i,j]
    # lower_lme <- Out$lme.beta[i,j] - 1.96*Out$SE.lme[i,j]
    # 
    # upper_W1 <- Out$lme2.beta[i,j] + 1.96*Out$SE.lme.W1[i,j]
    # lower_W1 <- Out$lme2.beta[i,j] - 1.96*Out$SE.lme.W1[i,j]
    # 
    # upper_W2 <- Out$lme3.beta[i,j] + 1.96*Out$SE.lme.W2[i,j]
    # lower_W2 <- Out$lme3.beta[i,j] - 1.96*Out$SE.lme.W2[i,j]
    
    Coverage.NR[i,j] <- ifelse(Beta.tr[j]<= upper_NR &Beta.tr[j]>= lower_NR,1,0)
    # Coverage.lme[i,j] <- ifelse(Beta.tr[j]<= upper_lme & Beta.tr[j]>= lower_lme,1,0)
    # Coverage.W1[i,j] <- ifelse(Beta.tr[j]<= upper_W1 & Beta.tr[j]>= lower_W1,1,0)
    # Coverage.W2[i,j] <- ifelse(Beta.tr[j]<= upper_W2 & Beta.tr[j]>= lower_W2,1,0)
  }
}

coverage <- rbind(colMeans(Coverage.NR),colMeans(Coverage.lme),colMeans(Coverage.W1),colMeans(Coverage.W2))
xtable(coverage)

#### variance
head(unlist(Out$SD.lme))
mean(unlist(Out$SD.lme)^2)
mean(Out$NR.theta[,1]) 

var <- rbind(mean((Out$NR.theta[,1])),mean(unlist(Out$SD.lme)^2),mean(unlist(Out$SD.lme.W1)^2),mean(unlist(Out$SD.lme.W2)^2))
xtable(var)

length(unlist(Out$SD.lme)^2)



#### MSE

