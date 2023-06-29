
var.clusterClaudia <- function(N,idK,  hh,sigma,beta,delta, Y, Z, X, R,pi){
  IdKs <- unique(idK[R==1])
  piK <- aggregate( pi~idK, FUN = function(pi)  1- prod(1-pi)  )
  piK0 <- piK[,2]
  names(piK0) <- piK[ , 1]
  piKvector <- NULL
  piKvector <- piK0[as.character(idK)]
  
  
  Uc <- Uci(Y, X,Z, R, beta, delta,hh,sigma)
  U1 <- R*matrix(unlist(Uc),ncol=length(beta),byrow = T ) 
  
  Up <- GUphi(Y, Z, R,  delta  ) 
  U2 <- matrix(unlist(Up),ncol=length(delta),byrow = T  )
  Rk <- as.numeric(idK %in% unique(idK[R==1]))
  
  UI  <- cbind(U1,U2)
  
  
  
  
  V1 <- 0
  
  
  
  for (k in 1: length(IdKs)){
    #k=1
    idK.r <-  (idK == k)
    piKvectorKK<- unique(piKvector[idK.r ])
    
    if(sum(idK.r)==1){
      Uk<-  UI[idK.r,]  
    }
    
    if(sum(idK.r) >  1){
      Uk<-  colSums(UI[idK.r,] )
    }
    
    V1 <- V1 +   Uk%*%t(Uk) 
  }
  
  ## Second-phase
  V2 <- V1*0 
  
  V2[1:length(beta), 1:length(beta)] <- t(R*U1* as.vector((pi-pi^2)/pi))%*%(U1    )
  
  
  return(V1+V2)
  
}

#NRest <- tryCatch(NR(Y,X,Z,R,beta,delta,hh,sigmaYk, eps0, N,idK ,pi), error=function(e) 'Skip')

#sqrt(diag(NRest$Covt2))
#sqrt(diag(NRest$CovTheta))


#NRest <- tryCatch(NR(Y,X,Z,R,beta,delta,hh,sigmaYk, eps0, N,idK ,pi), error=function(e) 'Skip')
#sqrt(diag(NRest$CovTheta))
#sqrt(diag(NRest$Covt2))


# 
# 
# ff00 <- function(betadelta,...){ ## pbeta, pdelta, Y, X,Z, R, hh,sigma 
#   beta <- betadelta[1:pbeta ]
#   delta <- betadelta[ (pbeta+1): (pdelta +pbeta)]
#   ff00<- normS(Y, X,Z, R, beta , delta,hh,sigma )
#  ff00
# }
#   
#   library(rootSolve)
#  df00 <-  gradient(ff00, x= c(betan_1,deltan_1), pbeta=pbeta, pdelta= pdelta, Y=Y, X=X,
#                   Z=Z, R=R, hh=hh,sigma=sigma  )
#  
#  






NR <- function(Y,X,Z,R,beta,delta,hh,sigma, eps0, N,idK, pi ){
  
  theta_n1 <- theta_n  <- c(beta,delta)
  epsiV <- NULL
  pdelta <- length(delta)
  pbeta <- length(beta)
  epsi <-  1  +eps0
  counti <- 1
  f00  <- normS(Y, X,Z, R, beta, delta ,hh,sigma ) 
  while(epsi  > eps0  )
  {
    
    
    betan_1 <- theta_n1[1:pbeta]
    deltan_1 <- theta_n1[(pbeta+1):(pbeta+pdelta)]
    
    #S_c <- normS(Y, X,Z, R, betan_1 , deltan_1,hh,sigma ) #
    
    f  <-   normS(Y, X,Z, R, betan_1 , deltan_1,hh,sigma )  ## S_C, psi
    
    #Derivative of dS_epslon w.r.t. beta:
    d11 <- normdS_beta(Y,X,Z,R,betan_1,deltan_1,hh,sigma)  ## d S_C/ d beta
    d12 <-  (normdS_delta(Y,X,Z,R,betan_1 ,deltan_1,hh,sigma) ) ## d S_C/ d delta
    
    d21  <-  matrix(0, ncol = pbeta ,nrow = pdelta) ## d psi/ d beta, whether it is a matrix
    
    d22  <- dphi_delta(R,Y,Z,deltan_1 )  ## d psi/ d delta
    
    
    df   <-  (rbind( cbind(d11, d12), cbind( (d21), d22)    ))
    
    
    ###where the function 
    dfI  <- solve(t(df) %*% df)
    df2 <-   dfI %*%t(df)
    
    
    theta_n1  <- theta_n  -    ((df2))%*%  f
    
    epsi  <- max(abs(theta_n1 - theta_n ))
    
    theta_n <- theta_n1
    #Update epsilon:
    
    
    epsiV[counti] <- epsi
    
    print(counti)
    #print(epsi)
    #print(theta_n1)
    counti  <- counti +1
    
  }
  
  
  
  v1= varbeta(Y, X,Z, R, theta_n1[1:pbeta],  theta_n1[(pbeta+1):(pbeta+pdelta)] ,hh,sigma) ####Need to check
  v2= TCovbeta(Y, X,Z, R, theta_n1[1:pbeta],  theta_n1[(pbeta+1):(pbeta+pdelta)] ,hh,sigma)
  v3= Covbeta(Y, X,Z, R, theta_n1[1:pbeta],  theta_n1[(pbeta+1):(pbeta+pdelta)],hh,sigma)
  v4= varphi(Y, Z, R, X,  theta_n1[(pbeta+1):(pbeta+pdelta)] ,hh, theta_n1[1:pbeta], sigma  )
  
  CovSTheta <-  rbind(  cbind(v1,v2), cbind(v3,v4))
  
  #Jtheta    <-  df 
  
  
  
  CovTheta   <-  ((df2))%*% CovSTheta  %*%  t( (df2))    # If does not work, try t(df) first
  
  #sqrt(diag(CovTheta))
  
  CovTheta2 <- -solve(df)
  
  #Uc <-  normU(Y, X,Z, R, beta=theta_n1[1:pbeta], delta=theta_n1[(pbeta+1):(pbeta+pdelta)],hh,sigma)
  
  #Up <- normUphi(Y, Z, R,  delta=theta_n1[(pbeta+1):(pbeta+pdelta)]  )
  # 
  
  #V1  <- clustereff(N, idK,  hh,sigma,beta=theta_n1[1:pbeta], delta=theta_n1[(pbeta+1):(pbeta+pdelta)], Y, Z, X, R)
  
  #V2 <- var.check(N, idK,  hh,sigma,beta=theta_n1[1:pbeta], delta=theta_n1[(pbeta+1):(pbeta+pdelta)], Y, Z, X, R)
  
  V3 <- var.clusterClaudia(N, idK,  hh,sigma ,beta=theta_n1[1:pbeta], delta=theta_n1[(pbeta+1):(pbeta+pdelta)], Y, Z, X, R, pi)
  
  #Vbeta <- solve((df))%*%V1  %*% t(solve(df))
  
  #Vbeta2 <- solve((df))%*%V2  %*% t(solve(df))
  
  Vbeta3 <-  ((df2))%*%V3  %*% t( (df2))
  
  list(theta=theta_n1,beta=theta_n1[1:pbeta],delta=theta_n1[(pbeta+1):(pbeta+pdelta)],
       cov1 = v1 ,cov2 =v2 , cov3= v3 ,cov4= v4, CovTheta= CovTheta,Covt2=Vbeta3,CovTheta2 = CovTheta2, UI= (f00  )   )
  
  
}

#######Simulation for variance
library(MASS)

#### Function for creating 

YCluster <- function(idK, mu0, sigma, sigmau){
  Y <- NULL
  idKA <- unique(idK)
  for(j in idKA){
    idKk  <- idK == j
    Nk   <- sum(idKk)  
    muk  <- mu0[idKk]
    uk<- rnorm(1, 0 , sigmau)  
    Y[idKk] <-  rnorm(Nk,muk,sigma) +uk
    
  }
  
  Y
  
}


steps <-1# 1000

K  <- 50
Nk <- 50###increase the cluster number to compare
N <- K*Nk
sigmau <- 1 # change to 0.5

beta  <- c(-2,1,1,1,1 )  

delta <- c(-6,-1,-1,-1 )

lambda <- c(-1,1)

sigma  <- 1

Z1  <- rnorm(N, 0,1) ####a population
Z2  <- rnorm(N, 0,1)
px1 <- 0.25
px2 <- 0.5
X1  <- rnorm(N, 0, 1)
X2  <- rbinom (N, 1,  1/(1+exp(-cbind(1, Z2)%*%lambda)))
mean(X2) ###about 0.3

Z    <- cbind(Z1,Z2 )
X   <- cbind(X1,X2 )
XX    <- cbind(1, X, Z)



eta   <- XX%*%beta
mu     <- function(eta)  eta  #identity link  # 1/(1+exp(-eta)) 





bvalue <-length(beta)
dvalue <-length(delta)
BetaM <- matrix(NA, nrow =  (steps), ncol =  (bvalue))
DeltaM = matrix(NA, nrow = (steps), ncol =  (dvalue))
Score<-  matrix(NA, nrow =  (steps), ncol =  (bvalue + dvalue))
cov=list() #empty list
covK=list()
covT1=list()
covT2=list()
covdf=list()

####for coverage
Cove_ind <- matrix(NA, nrow =  (steps), ncol =  (bvalue))
Cove_K <- matrix(NA, nrow =  (steps), ncol =  (bvalue))


#######
hh <- 0.1
epsi <- 1
eps0  <- 0.01
#####can change the number of cluster

#cluster.s <- id %in% sample(which(Ri == "1"),ncluster)
#beta <-  glm(Y ~X +Z ,family = gaussian)$coefficients


i=1
CountErrors <- 0
#Trybeta <- Results1$BetaM[or ,]
#Trydelta <- Results1$DeltaM[or ,]


idK <-  rep(1:K, Nk)
uK <- rep( rnorm(K, 0, sigmau), Nk)
uK <- uK[order(idK )]
idK <- sort(idK)

sigmaYk <-   sqrt( sigma^2 + sigmau^2) 



while (i <=steps){ 
  Y   <-  YCluster(idK,  mu0 = mu(eta) ,sigma , sigmau) 
  
  pi <- calc_pi(Y, Z, delta)
  
  
  id <- 1:nrow(XX)
  
  R    <- rbinom(N, 1, pi)
  sum(R) #1: 357 3:393 5: 302
  
  
  NRest <- tryCatch(NR(Y,X,Z,R,beta,delta,hh,sigmaYk, eps0, N,idK ,pi), error=function(e) 'Skip')
  if(class(NRest)!='list')
  {
    i=i
    CountErrors <- CountErrors+1
  }

  
  if(class(NRest) == 'list')
  {
    BetaM[i,] <- NRest$beta
    Score[i,] <- NRest$UI
    DeltaM[i,] <- NRest$delta
    cov[[i]] <-  NRest$CovTheta
    covK[[i]] <- NRest$Covnew ### Vbeta1 
    covdf[[i]] <- NRest$CovTheta2 ### Vbeta1 
    covT2[[i]] <- NRest$Covt2 ###Vbeta3 method combine (Uci,Uphi)
    
    
    for(j in 1:bvalue){
      
      upper_ind <- BetaM[i,j]+1.96*sqrt(diag(cov[[i]])[j])
      lower_ind <- BetaM[i,j]-1.96*sqrt(diag(cov[[i]])[j])
      
      
      upper_K <- BetaM[i,j]+1.96*sqrt(diag(covT2[[i]])[j])
      lower_K <- BetaM[i,j]-1.96*sqrt(diag(covT2[[i]])[j])
      
      
      Cove_ind[i,j] <-ifelse(beta[j]<=upper_ind & beta[j]>=lower_ind,1,0) ### mean() :matrix as well
      Cove_K[i,j] <- ifelse(beta[j]<=upper_K & beta[j]>=lower_K,1,0)
    }   
    
    
    
    print(i)
    print(CountErrors)
    
    # if(i>2){
    #   
    #   print(apply( BetaM[1:i,],2, mad)) 
    #   
    print(sqrt(diag(var(BetaM[1:i,]))))
    print(sqrt( diag(cov[[i]])))
    print(sqrt( diag(covT2[[i]])))
    print(sqrt( diag(covdf[[i]])))
    print(Cove_ind[i,j])
    print(Cove_K[i,j])
    
    # }
    
    
    
    
    
    
    i <- i+1
    
    
    
    
  }
  
  
}


Results <-list(BetaM=BetaM,DeltaM=DeltaM, CovInd= cov , CovN=covK,covcluster=covT2, Cove_ind=Cove_ind, Cove_K=Cove_K)

saveRDS(Results, paste('Output/Coverage8_','05-09-2021_K', K,'_Nk',Nk, '_sigmau', sigmau,'_hh',hh, '_steps', i , '_XY',  '.rds', sep=""))


###For the coverage 
            

GG <- i-1
### sigmau=1

HH<- HH2<-array(NA, c(6,6,GG))
for(j in 1:GG){
  HH[,,j]  <- Results$CovInd[[j]]
  HH2[,,j] <- Results$covcluster[[j]]
  
}

sqrt(diag( apply(HH, c(1,2 ), FUN='mean'))) #estimate for the independent
sqrt(diag( apply(HH2, c(1,2 ), FUN='mean'))) #estimate for the cluster effect


GG=10
round(sqrt(diag(var(Results$BetaM[1:GG,]))),3)
round(sqrt(diag(Reduce('+',Results$CovInd[1:GG])/GG)),3)
round(sqrt(diag(Reduce('+',Results$covcluster[1:GG])/GG)),3)
sqrt(diag( apply(HH2, c(1,2 ), FUN='mean')))

apply(BetaM[1:GG,], 2, mad)


