##### All Black-Scholes functions #####
## Efficient simulation of stock price process from Black-Scholes model (Geometric Brownian Motion) - inspired by https://bit.ly/37IROB4 #
# nSteps determines the number of steps in simulated paths
# nPaths determines the number of simulated paths
simulatedBlackScholesSPP <- function(T,S0,r,sigma,nSteps,nPaths){
  r <- (1+r)^(1/(252*24*60))-1                                                                    # converting annualized risk-free rate to minutely risk-free rate
  dt <- T/nSteps
  logS0 <- log(S0)
  eps <- matrix(rnorm(nPaths*nSteps,mean=0,sd=1),nrow=nPaths,ncol=nSteps)
  dlogS <- (r-sigma^2/2)*dt+sigma*sqrt(dt)*eps
  logS <- logS0+t(apply(dlogS,1,cumsum))                                                          # 't(apply(matrix,1,cumsum))' calculates cumsum across columns
  S <- exp(logS)
  return(cbind(S0,S))
}

## Pricing European call option via Black Scholes
BSCallPrice <- function(T,K,S,r,sigma){
  r <- (1+r)^(1/(252*24*60))-1                                                                    # converting annualized risk-free rate to minutely risk-free rate
  d1 <- (log(S/K)+(r+sigma^2/2)*T)/(sigma*sqrt(T))
  d2 <- d1-sigma*sqrt(T)
  value <- S*pnorm(d1)-K*exp(-r*T)*pnorm(d2)
  return(value)
}

## Delta via Black-Scholes
BSDelta <- function(T,K,S,r,sigma){
  r <- (1+r)^(1/(252*24*60))-1                                                                    # converting annualized risk-free rate to minutely risk-free rate
  d1 <- (log(S/K)+(r+sigma^2/2)*T)/(sigma*sqrt(T))
  delta <- pnorm(d1)
  return(delta)
}

##### All normal inverse Gaussian exponential Lévy functions ####
## Simulating random variables that follow inverse Gaussian distribution
simulatingIGVariable <- function(n,mu,lambda){
  Y <- rnorm(n)^2
  X <- mu+(mu^2*Y)/(2*lambda)-mu/(2*lambda)*sqrt(4*mu*lambda*Y+mu^2*Y^2)
  U <- runif(n)
  igvar <- ifelse(U<=mu/(mu+X),X,mu^2/X)
  return(igvar)
}

## Simulating normal inverse Gaussian Lévy process
simulatedNIGLP <- function(T,n,kappa,sigma,eta){
  mu <- T/n
  lambda <- mu^2/kappa
  deltaS <- simulatingIGVariable(n,mu,lambda)
  N <- rnorm(n)
  deltaX <- sigma*N*sqrt(deltaS)+eta*deltaS
  return(c(0,cumsum(deltaX)))
}

## Simulating normal inverse Gaussian exponential Lévy process
# nSteps determines the number of steps in simulated paths
# nPaths determines the number of simulated paths
simulatedNIGLSPP <- function(T,S0,r,kappa,sigma,eta,nSteps,nPaths){
  r <- (1+r)^(1/(252*24*60))-1                                                                    # converting annualized risk-free rate to minutely risk-free rate
  dt <- T/nSteps
  times <- seq(0,T,dt)
  logS0 <- log(S0)
  X <- t(replicate(nPaths,simulatedNIGLP(T,nSteps,kappa,sigma,eta)))
  logS <- log(S0)+r*times+X
  S <- exp(logS)
  return(S)
}

## Price of European call option via normal inverse Gaussian exponential Lévy process using Monte Carlo integration
# M determines the number of Monte Carlo repetitions
NIGLSPPCallPrice <- function(T,K,S,r,kappa,sigma,eta,nSteps,M,exclude.std=T){
  paths <- simulatedNIGLSPP(T,S,r,kappa,sigma,eta,nSteps,M)                                       # Beware nPaths = M
  payoff <- pmax(paths[,ncol(paths)]-K,0)
  value <- exp(-r*T)*mean(payoff)
  std <- sd(payoff)/sqrt(M)
  if(exclude.std){return(value)}else{return(list('value'=value,'std'=std))}
}

##### Black-Scholes offline grid #####
# To run this part of the code the functions 'BSCallPrice()' and 'BSDelta()' are required

# Trange is a vector of start time and end time, e.g., Trange = c(0,1)
# Srange is a vector of start price and end price, e.g., Srange = c(0,300)
# GridPartition is a tuning parameter determining the partition of the grid, e.g., GridPartition = 1/1000
simulatedBSGrid <- function(Trange,Srange,GridPartitionT,GridPartitionS,K,r,sigma){
  times <- seq(Trange[1],Trange[2],GridPartitionT)
  prices <- seq(Srange[1],Srange[2],GridPartitionS)
  #
  CallPriceGrid <- matrix(nrow=length(prices),ncol=length(times))
  for(t in 1:length(times)){
    for(s in 1:length(prices)){
      CallPriceGrid[s,t] <- BSCallPrice((Trange[2]-times[t]),K,prices[s],r,sigma)
    }
  }
  #
  DeltaGrid <- matrix(nrow=nrow(CallPriceGrid),ncol=(ncol(CallPriceGrid)-1))
  for(t in 1:(ncol(CallPriceGrid)-1)){
    for(s in 2:nrow(CallPriceGrid)){
      DeltaGrid[s,t] <- BSDelta(T-times[t],K,prices[s],r,sigma)
      #DeltaGrid[s,t] <- (CallPriceGrid[s,t]-CallPriceGrid[s-1,t])/GridPartitionS 
    }
  }
  return(list(CallPriceGrid,DeltaGrid,par=list('Trange'=Trange,'Srange'=Srange,
                                               'GridPartitionT'=GridPartitionT,'GridPartitionS'=GridPartitionS,
                                               'K'=K,'r'=r,'sigma'=sigma)))
}
system.time(BlackScholesGrid <- simulatedBSGrid(Trange=c(0,1),Srange=c(0,5),GridPartitionT=1/1000,
                                                GridPartitionS=1/100,K=1,r=0.02,sigma=0.18))
# http://www.sthda.com/english/wiki/saving-data-into-r-data-format-rds-and-rdata
# saveRDS(BlackScholesGrid, file = "BlackScholesGrid.rds")


##### Normal inverse Gaussian exponential Lévy offline grid #####
# To run this part of the code the function 'NIGLSPPCallPrice(s)' is required, and
# this function in turn depends on the functions 'simulatingIGVariable()', 'simulatedNIGLP()', and 'simulatedNIGLSPP()'

# Trange is a vector of start time and end time, e.g., Trange = c(0,1)
# Srange is a vector of start price and end price, e.g., Srange = c(0,300)
# GridPartition is a tuning parameter determining the partition of the grid, e.g., GridPartition = 1/1000

# Note that the following function is similar to 'simulatedBSGrid()'
# Parallelization is done following: https://www.blasbenito.com/post/02_parallelizing_loops_with_r/
library(foreach);library(doParallel);library(ranger);library(tidyverse);library(kableExtra)
simulatedNIGLGrid <- function(Trange,Srange,GridPartitionT,GridPartitionS,K,r,kappa,sigma,eta,M){
  times <- seq(Trange[1],Trange[2],GridPartitionT)
  prices <- seq(Srange[1],Srange[2],GridPartitionS)
  #
  my.cluster <- parallel::makeForkCluster(nnodes=48)   # creates cluster
  doParallel::registerDoParallel(cl = my.cluster)      # register cluster to be used by %dopar%
  CallPriceGrid <- 
    foreach(t = 1:length(times),.combine='cbind') %:%
    foreach(s = 1:length(prices),.combine='c') %dopar%{
      NIGLSPPCallPrice((Trange[2]-times[t]),K,prices[s],r,kappa,sigma,eta,nSteps=1000,M)    #nSteps?
    }
  parallel::stopCluster(cl = my.cluster)
  #
  DeltaGrid <- matrix(nrow=nrow(CallPriceGrid),ncol=(ncol(CallPriceGrid)-1))
  for(t in 1:(ncol(CallPriceGrid)-1)){
    for(s in 2:nrow(CallPriceGrid)){
      DeltaGrid[s,t] <- (CallPriceGrid[s,t]-CallPriceGrid[s-1,t])/GridPartitionS 
    }
  }
  return(list(CallPriceGrid,DeltaGrid,par=list('Trange'=Trange,'Srange'=Srange,'GridPartitionT'=GridPartitionT,'GridPartitionS'=GridPartitionS,
                                               'K'=K,'r'=r,'kappa'=kappa,'sigma'=sigma,'eta'=eta,'M'=M)))
}
system.time(NIGLGridBackup <- simulatedNIGLGrid(Trange=c(0,1),Srange=c(0,5),GridPartitionT=1/1000,GridPartitionS=1/100,
                                                K=1,r=0.02,kappa=0.98,sigma=0.18,eta=-0.0047,M=10000))
# saveRDS(NIGLGridBackup, file = "NIGLGridBackup.rds")

##### Black-Scholes Simulation study #####
# To run this part of the code the function 'simulatedBlackScholesSPP()' is required.

# nHedge determines the number of hedges
# offline is a list containing [[1]]CallPriceGrid, [[2]]DeltaGrid, 
# and [[3]] list of parameters (Trange,Srange,GridPartition,K,r,Sigma) used to construct grids [[1]] and [[2]]
# Suggestion: instead of specifying Trange and Srange a priori, make them based on Trange and Srange of simulated paths?
simulatedDeltaHedgingBS <- function(T,K,S0,r,sigma,nHedge=nSteps,nSteps,nPaths,offline){
  #if(nHedge>nSteps){
  #  stop('Number of hedges should be smaller than or equal to number of steps in simulated paths.')
  #} else if(nSteps %% nHedge){
  #  stop('nSteps modulo nHedge needs to be zero.')
  #} else if(offline[[3]]$Trange[2]<T){
  #  stop('Upper bound of Trange in supplied offline grid needs to be greater than or equal to option maturity time T.')
  #} else if(offline[[3]]$K!=K){
  #  stop('Strike price needs to match strike price used in supplied offline grid')
  #} else if(offline[[3]]$r!=r){
  #  stop('Risk free rate needs to match risk free rate used in supplied offline grid')
  #} else if(offline[[3]]$sigma!=sigma){
  #  stop('Volatility needs to match volatility used in supplied offline grid')
  #}
  #
  BankBalance <- matrix(nrow=nPaths,ncol=(nHedge+1))
  PortfolioValue <- matrix(nrow=nPaths,ncol=(nHedge+1))
  #
  dt <- T/nHedge
  times <- seq.int(0,T,dt)
  hedgeIndex <- seq.int(1,nSteps+1,nSteps/nHedge)
  #
  CallPriceGrid <- offline[[1]]
  DeltaGrid <- offline[[2]]
  GridPartition <- offline$par$GridPartitionS
  #
  breakcond <- TRUE
  while(breakcond){
    paths <- simulatedBlackScholesSPP(T,S0,r,sigma,nSteps,nPaths)
    paths <- round(paths,digits=(sapply(GridPartition,nchar)-2))                                  # 'sapply(GridPartition, nchar)-2' calculates number of digits in GridPartition, e.g. for 1/1000 output is 3
    breakcond <- max(paths)>offline$par$Srange[2]                                                 # ensuring that no simulated stock price is larger than largest stock price in supplied offline grid
  }
  #
  Sindices <- paths/GridPartition+1; Sindices <- Sindices[,hedgeIndex]
  Tindices <- tcrossprod(rep(1,nPaths),times/GridPartition+1)#; Tindices <- Tindices[,hedgeIndex]
  #
  CallPrices <- array(CallPriceGrid[cbind(c(Sindices),c(Tindices))],dim(Tindices))
  Deltas <- array(DeltaGrid[cbind(c(Sindices[,-(nHedge+1)]),c(Tindices[,-(nHedge+1)]))],dim(Tindices[,-(nHedge+1)]))  # Removing last row as this corresponds to delta at maturity time
  #
  BankBalance[,1] <- CallPrices[,hedgeIndex[1]]-Deltas[,hedgeIndex[1]]*S0                         # Bank balance at time 0
  PortfolioValue[,1] <- -BankBalance[,1]                                                          # Portfolio value at time 0
  #
  r <- (1+r)^(1/(252*24*60))-1                                                                    # converting annualized risk-free rate to minutely risk-free rate
  for(t in 2:nHedge){                                                                             # Loop through intermediate times, that is, 0<t<T
    S <- paths[,hedgeIndex[t]]
    BankBalance[,t] <- exp(r*dt)*BankBalance[,t-1]-(Deltas[,t]-Deltas[,t-1])*S
    PortfolioValue[,t] <- -CallPrices[,t]+Deltas[,t]*S
  }
  #
  S <- paths[,nSteps+1]                                                                           # Stock price at option maturity
  BankBalance[,nHedge+1] <- exp(r*dt)*BankBalance[,nHedge]+Deltas[,nHedge]*S-pmax(S-K,0)          # Bank balance at option maturity
  PortfolioValue[,nHedge+1] <- -pmax(S-K,0)+Deltas[,nHedge]*S                                     # Portfolio value at option maturity
  return(list(paths,Sindices,Tindices,CallPrices,Deltas,BankBalance,PortfolioValue,'payoff'=pmax(S-K,0),'times'=times,
              'par'=list('T'=T,'K'=K,'S0'=S0,'r'=r,'sigma'=sigma,'nHedge'=nHedge,'nSteps'=nSteps,'nPaths'=nPaths)))
}
system.time(BlackScholesSimHedging1 <- simulatedDeltaHedgingBS(T=1,K=1,S0=1,r=0.02,sigma=0.18,
                                                               nSteps=1000,nPaths=100000,offline=BlackScholesGrid))
 saveRDS(BlackScholesSimHedging1, file = "BlackScholesSimHedging1.rds")
system.time(BlackScholesSimHedging2 <- simulatedDeltaHedgingBS(T=1,K=1,S0=1.15,r=0.02,sigma=0.18,
                                                               nSteps=1000,nPaths=100000,offline=BlackScholesGrid))
 saveRDS(BlackScholesSimHedging2, file = "BlackScholesSimHedging2.rds")
system.time(BlackScholesSimHedging3 <- simulatedDeltaHedgingBS(T=1,K=1,S0=0.85,r=0.02,sigma=0.18,
                                                               nSteps=1000,nPaths=100000,offline=BlackScholesGrid))
 saveRDS(BlackScholesSimHedging3, file = "BlackScholesSimHedging3.rds")

###
system.time(BlackScholesSimHedging10 <- simulatedDeltaHedgingBS(T=1,K=1,S0=1,r=0.02,sigma=0.18,nHedge=10,
                                                               nSteps=1000,nPaths=1000,offline=BlackScholesGrid))
 saveRDS(BlackScholesSimHedging10, file = "BlackScholesSimHedging10.rds")
system.time(BlackScholesSimHedging100 <- simulatedDeltaHedgingBS(T=1,K=1,S0=1,r=0.02,sigma=0.18,nHedge=100,
                                                               nSteps=1000,nPaths=1000,offline=BlackScholesGrid))
 saveRDS(BlackScholesSimHedging100, file = "BlackScholesSimHedging100.rds")
system.time(BlackScholesSimHedging1000 <- simulatedDeltaHedgingBS(T=1,K=1,S0=1,r=0.02,sigma=0.18,nHedge=1000,
                                                               nSteps=1000,nPaths=1000,offline=BlackScholesGrid))
 saveRDS(BlackScholesSimHedging1000, file = "BlackScholesSimHedging1000.rds")





##### Normal inverse Gaussian exponential Lévy simulation study #####
# To run this part of the code the function 'simulatedNIGLSPP()' is needed, and
# this function in turn depends on the functions 'simulatingIGVariable()' and 'simulatedNIGLP'

# Note that the following function is EXACTLY the same as 'simulatedDeltaHedgingBS()'
# only difference is that 'simulatedBlackScholesSPP()' is replaced with 'simulatedNIGLSPP()'
simulatedDeltaHedgingNIGExpLevy <- function(T,K,S0,r,kappa,sigma,eta,nHedge=nSteps,nSteps,nPaths,offline){
  if(nHedge>nSteps){
    stop('Number of hedges should be smaller than or equal to number of steps in simulated paths.')
  } else if(nSteps %% nHedge){
    stop('nSteps modulo nHedge needs to be zero.')
  } else if(offline[[3]]$Trange[2]<T){
    stop('Upper bound of Trange in supplied offline grid needs to be greater than or equal to option maturity time T.')
  } else if(offline[[3]]$K!=K){
    stop('Strike price needs to match strike price used in supplied offline grid')
  } else if(offline[[3]]$r!=r){
    stop('Risk free rate needs to match risk free rate used in supplied offline grid')
  } else if(offline[[3]]$sigma!=sigma){
    stop('Volatility needs to match volatility used in supplied offline grid')
  }
  #
  BankBalance <- matrix(nrow=nPaths,ncol=(nHedge+1))
  PortfolioValue <- matrix(nrow=nPaths,ncol=(nHedge+1))
  #
  dt <- T/nHedge
  times <- seq.int(0,T,dt)
  hedgeIndex <- seq.int(1,nSteps+1,nSteps/nHedge)
  #
  CallPriceGrid <- offline[[1]]
  DeltaGrid <- offline[[2]]
  GridPartition <- offline$par$GridPartitionS
  #
  breakcond <- TRUE
  while(breakcond){
    paths <- simulatedNIGLSPP(T,S0,r,kappa,sigma,eta,nSteps,nPaths)
    paths <- round(paths,digits=(sapply(GridPartition,nchar)-2))
    breakcond <- max(paths)>offline$par$Srange[2]
  }
  #
  Sindices <- paths/GridPartition+1; Sindices <- Sindices[,hedgeIndex]
  Tindices <- tcrossprod(rep(1,nPaths),times/GridPartition+1)#; Tindices <- Tindices[,hedgeIndex]
  #
  CallPrices <- array(CallPriceGrid[cbind(c(Sindices),c(Tindices))],dim(Tindices))
  Deltas <- array(DeltaGrid[cbind(c(Sindices[,-(nHedge+1)]),c(Tindices[,-(nHedge+1)]))],dim(Tindices[,-(nHedge+1)]))
  #
  BankBalance[,1] <- CallPrices[,hedgeIndex[1]]-Deltas[,hedgeIndex[1]]*S0
  PortfolioValue[,1] <- -BankBalance[,1]
  #
  r <- (1+r)^(1/(252*24*60))-1                                                                      # converting annualized risk-free rate to minutely risk-free rate
  for(t in 2:nHedge){
    S <- paths[,hedgeIndex[t]]
    BankBalance[,t] <- exp(r*dt)*BankBalance[,t-1]-(Deltas[,t]-Deltas[,t-1])*S
    PortfolioValue[,t] <- -CallPrices[,t]+Deltas[,t]*S
  }
  #
  S <- paths[,nSteps+1]
  BankBalance[,nHedge+1] <- exp(r*dt)*BankBalance[,nHedge]+Deltas[,nHedge]*S-pmax(S-K,0)
  PortfolioValue[,nHedge+1] <- -pmax(S-K,0)+Deltas[,nHedge]*S
  return(list(paths,Sindices,Tindices,CallPrices,Deltas,BankBalance,PortfolioValue,'payoff'=pmax(S-K,0),'times'=times,
              'par'=list('T'=T,'K'=K,'S0'=S0,'r'=r,'sigma'=sigma,'nHedge'=nHedge,'nSteps'=nSteps,'nPaths'=nPaths)))
}
#
system.time(NIGExpLevySimHedging1 <- simulatedDeltaHedgingNIGExpLevy(T=1,K=1,S0=1,r=0.02,kappa=0.98,sigma=0.18,
                                                                    eta=-0.0047,nSteps=1000,nPaths=1000,
                                                                    offline=NIGLevyGrid))
 saveRDS(NIGExpLevySimHedging1, file = "NIGExpLevySimHedging1.rds")
system.time(NIGExpLevySimHedging2 <- simulatedDeltaHedgingNIGExpLevy(T=1,K=1,S0=1.15,r=0.02,kappa=0.98,sigma=0.18,
                                                                    eta=-0.0047,nSteps=1000,nPaths=1000,
                                                                    offline=NIGLevyGrid))
 saveRDS(NIGExpLevySimHedging2, file = "NIGExpLevySimHedging2.rds")
system.time(NIGExpLevySimHedging3 <- simulatedDeltaHedgingNIGExpLevy(T=1,K=1,S0=0.85,r=0.02,kappa=0.98,sigma=0.18,
                                                                    eta=-0.0047,nSteps=1000,nPaths=1000,
                                                                    offline=NIGLevyGrid))
 saveRDS(NIGExpLevySimHedging3, file = "NIGExpLevySimHedging3.rds")

##
system.time(NIGExpLevySimHedging10 <- simulatedDeltaHedgingNIGExpLevy(T=1,K=1,S0=1,r=0.02,kappa=0.98,sigma=0.18,nHedge=10,
                                                                     eta=-0.0047,nSteps=1000,nPaths=1000,
                                                                     offline=NIGLevyGrid))
 saveRDS(NIGExpLevySimHedging10, file = "NIGExpLevySimHedging10.rds")
system.time(NIGExpLevySimHedging100 <- simulatedDeltaHedgingNIGExpLevy(T=1,K=1,S0=1.15,r=0.02,kappa=0.98,sigma=0.18,nHedge=100,
                                                                     eta=-0.0047,nSteps=1000,nPaths=1000,
                                                                     offline=NIGLevyGrid))
 saveRDS(NIGExpLevySimHedging100, file = "NIGExpLevySimHedging100.rds")
system.time(NIGExpLevySimHedging1000 <- simulatedDeltaHedgingNIGExpLevy(T=1,K=1,S0=0.85,r=0.02,kappa=0.98,sigma=0.18,nHedge=1000,
                                                                     eta=-0.0047,nSteps=1000,nPaths=1000,
                                                                     offline=NIGLevyGrid))
 saveRDS(NIGExpLevySimHedging1000, file = "NIGExpLevySimHedging1000.rds")




##### Plots/Figures #####
# To run this part of the code you need to get a life :-)
library(parallel); library(ggplot2); library(plotly); library(processx)
library(orca)

## (1) Determining the best number of MC repetitions
library(parallel); library(ggplot2)
MRepTests <- c(1,seq(5,1000,5),seq(1100,10000,100),seq(11000,100000,1000))
set.seed(123)
system.time(MApproxRes <- parallel::mclapply(1:length(MRepTests),function(x){
  NIGLSPPCallPrice(T=(1-0.88),K=1,S=1.21,r=0.02,kappa=0.98,sigma=0.18,eta=-0.0047,
                   nSteps=1000,M=MRepTests[x],exclude.std=F)
},mc.cores=48))
# saveRDS(MApproxRes, file = "MApproxRes.rds")
dfpoints <- numeric(length(unlist(lapply(MApproxRes,'[[',1))))
dfpoints[which(MApproxRes %in% c(5,10,100,1000,5000,10000,25000,50000,75000,100000))] <- 1
MApproxRes2 <- cbind(MRepTests,unlist(lapply(MApproxRes,'[[',1)),dfpoints)
MApproxRes2.1 <- cbind(MRepTests[241:381],unlist(lapply(MApproxRes,'[[',1))[241:381])
MApproxRes2.2 <- cbind(MRepTests[331:381],unlist(lapply(MApproxRes,'[[',1))[331:381])


ggplot(data=as.data.frame(MApproxRes2), aes(x=MRepTests,y=V2)) + 
  geom_line() +
  labs(x='M',y='Call Option Price')

ggplot(data=as.data.frame(MApproxRes2.1), aes(x=MRepTests[241:381],y=V2)) + 
  geom_line() +
  labs(x='M',y='Call Option Price')

ggplot(data=as.data.frame(MApproxRes2.2), aes(x=MRepTests[331:381],y=V2)) + 
  geom_line() +
  labs(x='M',y='Call Option Price')


MsdRes.df <- unlist(lapply(MApproxRes,'[[',2))[which(MRepTests %in% c(5,10,100,1000,5000,10000,
                                                                   25000,50000,75000,100000))]
MApproxRes.df <- unlist(lapply(MApproxRes,'[[',1))[which(MRepTests %in% c(5,10,100,1000,5000,10000,
                                                                       25000,50000,75000,100000))]
MResdf <- cbind(MsdRes.df,MApproxRes.df)

## (2) Implied volatility
ImpliedVolatilityBS <- function(T,K,S,r,CM){
  f <- function(sigma){return(BSCallPrice(T,K,S,r,sigma)-CM)}
  value <- uniroot(f,lower=-100000,upper=100000,maxiter=1000000000,tol=1/1000000000000000)$root
  return(value)
}
ImpliedVolatilityBS(T=1,K=100,S=100,r=0.02,CM=10.5)
BSCallPrice(T=1,K=100,S=100,r=0.02,sigma=0.2639443)
times <- seq(0,1,1/1000)
prices <- seq(0,5,1/100)
IVdata <- matrix(nrow=length(prices),ncol=length(times))
CMdata <- BlackScholesGrid[[1]]; CMdata[,ncol(CMdata)] <- pmax(seq(0,5,1/100)-1,0)
CMdataTimesIndex <- seq(1,ncol(CMdata),1)
CMdataPricesIndex <- seq(1,nrow(CMdata),1)
k <- 1
for(t in 1:(length(times)-1)){
  for(s in 1:(length(prices-1))){
    IVdata[s,t] <- ImpliedVolatilityBS(T=(1-times[t]),K=1,S=prices[s],r=0.02,CM=CMdata[CMdataPricesIndex[s],CMdataTimesIndex[t]])
print(k); k <- k+1
  }
}

ImpliedVolatilityNIGExpLevy <- function(T,K,S,r,kappa,eta,nSteps,M,CM){
  f <- function(sigma){return(NIGLSPPCallPrice(T,K,S,r,kappa,sigma,eta,nSteps,M)-CM)}
  value <- uniroot(f,lower=0.05,upper=0.5)$root
  return(value)
}
ImpliedVolatilityNIGExpLevy(T=1,K=100,S=100,r=0.02,kappa=0.98,eta=-0.0047,nSteps=1000,M=1000,CM=10.5)
NIGLSPPCallPrice(T=1,K=100,S=100,r=0.02,kappa=0.98,sigma=0.2525051,eta=-0.0047,nSteps=1000,M=1000)
system.time(NIGLSPPCallPrice(T=1,K=100,S=100,r=0.02,kappa=0.98,sigma=0.2114684,eta=-0.0047,nSteps=1000,M=10000))
ImpliedVolatilityNIGExpLevy(T=1,K=1,S=5,r=0.02,kappa=0.98,eta=-0.0047,nSteps=1000,M=1000,CM=3.980237)

my.cluster <- parallel::makeForkCluster(nnodes=48)   # creates cluster
doParallel::registerDoParallel(cl = my.cluster)      # register cluster to be used by %dopar%
CallPriceGrid <- 
  foreach(t = 1:length(times),.combine='cbind') %:%
  foreach(s = 1:length(prices),.combine='c') %dopar%{
    IVdata[s,t] <- ImpliedVolatilityNIGExpLevy(T=(1-times[t]),K=1,S=prices[s],r=0.02,kappa=0.98,
                                               eta=-0.0047,nSteps=1000,M=10000,CM=CMdata[CMdataPricesIndex[s],CMdataTimesIndex[t]])
  }
parallel::stopCluster(cl = my.cluster)



## (3)

## (4) S shape
t <- seq(0,1,0.001)
s <- seq(0,5,0.01)[1:251]
data1 <- BlackScholesGrid[[2]]
l1.1 <- data1[,151][1:251]
l1.2 <- data1[,551][1:251]
l1.3 <- data1[,951][1:251]
df1 <- data.frame(x=s,y1=l1.1,y2=l1.2,y3=l1.3)
ggplot(data=df1, aes(x=s)) + 
  geom_line(aes(y=y1,color='t=0.15')) +
  geom_line(aes(y=y2,color='t=0.55')) +
  geom_line(aes(y=y3,color='t=0.95')) +
  labs(x='Stock Price',y='Delta') +
  xlim(0.4,1.8) +
  theme(legend.position = c(0.15,0.85), legend.key = element_blank(), legend.text = element_text(size=15), legend.title = element_blank(), legend.background=element_rect(color="black", linetype="solid")) +
  theme(legend.spacing.y = unit(-.01, "cm"))
#
data2 <- NIGLevyGrid[[2]]
l2.1 <- data2[,151][1:251]
l2.2 <- data2[,551][1:251]
l2.3 <- data2[,951][1:251]
df2 <- data.frame(x=s,y1=l2.1,y2=l2.2,y3=l2.3)
ggplot(data=df2, aes(x=s)) + 
  geom_line(aes(y=y1,color='t=0.15')) +
  geom_line(aes(y=y2,color='t=0.55')) +
  geom_line(aes(y=y3,color='t=0.95')) +
  labs(x='Stock Price',y='Delta') +
  xlim(0.4,2.1) +
  ylim(-0.5,2.5) +
  theme(legend.position = c(0.15,0.85), legend.key = element_blank(), legend.text = element_text(size=15), legend.title = element_blank(), legend.background=element_rect(color="black", linetype="solid")) +
  theme(legend.spacing.y = unit(-.01, "cm"))
# Ensuring that noise is due to MC

# saveRDS(DeltaGrid, file = "SShapeTest0150.rds")
# saveRDS(DeltaGrid, file = "SShapeTest0550.rds")
# saveRDS(DeltaGrid, file = "SShapeTest0950.rds")
SshapeTesting <- cbind(SShapeTest0150,SShapeTest0550,SShapeTest0950)
l3.1 <- SshapeTesting[,2][1:251]
l3.2 <- SshapeTesting[,4][1:251]
l3.3 <- SshapeTesting[,6][1:251]
df3 <- data.frame(x=s,y1=l3.1,y2=l3.2,y3=l3.3)
ggplot(data=df3, aes(x=s)) + 
  geom_line(aes(y=y1,color='t=0.15')) +
  geom_line(aes(y=y2,color='t=0.55')) +
  geom_line(aes(y=y3,color='t=0.95')) +
  labs(x='Stock Price',y='Delta') +
  xlim(0.4,2.1) +
  ylim(-0.5,2.5) +
  theme(legend.position = c(0.15,0.85), legend.key = element_blank(), legend.text = element_text(size=15), legend.title = element_blank(), legend.background=element_rect(color="black", linetype="solid")) +
  theme(legend.spacing.y = unit(-.01, "cm"))

## (5)
#1=paths, 2=Sindices, 3=Tindices, 4=CallPrices, 5=Deltas, 6=BankBalance, 7=Portfoliovalue

# a)
par(mfrow=c(1,1))
hist(BlackScholesSimHedging1[[6]][,ncol(BlackScholesSimHedging1[[6]])],
     xlim=c(-.2,.2),main='',xlab='Final Bank Balance',breaks=10)
hist(BlackScholesSimHedging2[[6]][,ncol(BlackScholesSimHedging2[[6]])],
     xlim=c(-.2,.2),main='',xlab='Final Bank Balance',breaks=20)
hist(BlackScholesSimHedging3[[6]][,ncol(BlackScholesSimHedging3[[6]])],
     xlim=c(-.2,.2),main='',xlab='Final Bank Balance',breaks=10)

hist(NIGExpLevySimHedging1[[6]][,ncol(NIGExpLevySimHedging1[[6]])],
     xlim=c(-.2,.2),main='',xlab='Final Bank Balance',breaks=20)
hist(NIGExpLevySimHedging2[[6]][,ncol(NIGExpLevySimHedging2[[6]])],
     xlim=c(-.2,.2),main='',xlab='Final Bank Balance',breaks=20)
hist(NIGExpLevySimHedging3[[6]][,ncol(NIGExpLevySimHedging3[[6]])],
     xlim=c(-.2,.2),main='',xlab='Final Bank Balance',breaks=50)


gghistogramhomemade <- function(series){
  df_data <- data.frame("series"=series)
  ggplot(df_data, aes(series)) + geom_histogram(aes(y=..density..), bins=nclass.FD(series)) +
    stat_function(fun=dnorm, color="red", args=list(mean = mean(series), sd = sd(series))) +    #theoretical Gaussian pdf
    geom_density(color="blue") +
    labs(x='Bank Balance at Maturity',y='Frequency')
}

# AT, IN, OUT, Respectively
gghistogramhomemade(BlackScholesSimHedging1[[6]][,ncol(BlackScholesSimHedging1[[6]])])
gghistogramhomemade(BlackScholesSimHedging2[[6]][,ncol(BlackScholesSimHedging2[[6]])])
gghistogramhomemade(BlackScholesSimHedging3[[6]][,ncol(BlackScholesSimHedging3[[6]])])

gghistogramhomemade(NIGExpLevySimHedging1[[6]][,ncol(NIGExpLevySimHedging1[[6]])])
gghistogramhomemade(NIGExpLevySimHedging2[[6]][,ncol(NIGExpLevySimHedging2[[6]])])
gghistogramhomemade(NIGExpLevySimHedging3[[6]][,ncol(NIGExpLevySimHedging3[[6]])])


library(reshape)
df.at_in_out <- data.frame(var1=BlackScholesSimHedging1[[6]][,ncol(BlackScholesSimHedging1[[6]])],
                           var2=BlackScholesSimHedging2[[6]][,ncol(BlackScholesSimHedging2[[6]])],
                           var3=BlackScholesSimHedging3[[6]][,ncol(BlackScholesSimHedging3[[6]])])
data.bs_at_in_out <- melt(df.at_in_out)
ggplot(data.bs_at_in_out,aes(x=value,fill=variable)) + 
  geom_density(alpha=.25) + labs(x='Bank Balance at Maturity',y='Empirical Density') +
  theme(legend.position="none")

df.at_in_out2 <- data.frame(var1=NIGExpLevySimHedging1[[6]][,ncol(NIGExpLevySimHedging1[[6]])],
                            var2=NIGExpLevySimHedging2[[6]][,ncol(NIGExpLevySimHedging2[[6]])],
                            var3=NIGExpLevySimHedging3[[6]][,ncol(NIGExpLevySimHedging3[[6]])])
data.bs_at_in_out2 <- melt(df.at_in_out2)
ggplot(data.bs_at_in_out2,aes(x=value,fill=variable)) + 
  geom_density(alpha=.25) + labs(x='Bank Balance at Maturity',y='Empirical Density') +
  theme(legend.position="none")



# 10, 100, 1000
gghistogramhomemade(BlackScholesSimHedging10[[6]][,ncol(BlackScholesSimHedging10[[6]])])
gghistogramhomemade(BlackScholesSimHedging100[[6]][,ncol(BlackScholesSimHedging100[[6]])])
gghistogramhomemade(BlackScholesSimHedging1000[[6]][,ncol(BlackScholesSimHedging1000[[6]])])

gghistogramhomemade(NIGExpLevySimHedging10[[6]][,ncol(NIGExpLevySimHedging10[[6]])])
gghistogramhomemade(NIGExpLevySimHedging100[[6]][,ncol(NIGExpLevySimHedging100[[6]])])
gghistogramhomemade(NIGExpLevySimHedging1000[[6]][,ncol(NIGExpLevySimHedging1000[[6]])])

df101001000bs <- data.frame(var1=BlackScholesSimHedging10[[6]][,ncol(BlackScholesSimHedging10[[6]])],
                             var2=BlackScholesSimHedging100[[6]][,ncol(BlackScholesSimHedging100[[6]])],
                             var3=BlackScholesSimHedging1000[[6]][,ncol(BlackScholesSimHedging1000[[6]])])
data101001000bs <- melt(df101001000bs)
ggplot(data101001000bs,aes(x=value,fill=variable)) + 
  geom_density(alpha=.25) + labs(x='Bank Balance at Maturity',y='Empirical Density') +
  theme(legend.position="none")

df101001000nig <- data.frame(var1=NIGExpLevySimHedging10[[6]][,ncol(NIGExpLevySimHedging10[[6]])],
                            var2=NIGExpLevySimHedging100[[6]][,ncol(NIGExpLevySimHedging100[[6]])],
                            var3=NIGExpLevySimHedging1000[[6]][,ncol(NIGExpLevySimHedging1000[[6]])])
data101001000nig <- melt(df101001000nig)
ggplot(data101001000nig,aes(x=value,fill=variable)) + 
  geom_density(alpha=.25) + labs(x='Bank Balance at Maturity',y='Empirical Density') +
  theme(legend.position="none")





# b)
par(mfrow=c(3,1))
plot(seq(0,1,1/1000),BlackScholesSimHedging1[[1]][333,],pch=20,
     xlab='Time',ylab='Asset Price'); abline(h=1)
plot(seq(0,1,1/1000),BlackScholesSimHedging1[[4]][333,],pch=20,
     xlab='Time',ylab='Call Option Price')
plot(seq(0,0.999,1/1000),BlackScholesSimHedging1[[5]][333,],pch=20,
     xlab='Time',ylab='Delta')
plot(BlackScholesSimHedging1[[6]][,ncol(BlackScholesSimHedging1[[6]])])

# c)
error <- BlackScholesSimHedging1[[6]][,ncol(BlackScholesSimHedging1[[6]])]
rmse <- sqrt(mean(error^2))
relativeerror <- rmse/1  #divided by initial price
# c.1)
h <- numeric(1000)
for(i in 1:length(h)){
  if(1000 %% i == 0){
    h[i] <- i
  }
}
h <- h[h!=0]
relativeerrorplotBS <- numeric(length(h))
for(i in 2:length(h)){
  res <- simulatedDeltaHedgingBS(T=1,K=1,S0=1,r=0.02,sigma=0.18,nHedge=,
                                 nSteps=1,nPaths=10000,offline=BlackScholesGrid)[[6]]
  error <- res[,ncol(res)]
  rmse <- sqrt(mean(error^2))
  relativeerrorplotBS[i] <- rmse/1
}
relativeerrorplotBS <- relativeerrorplotBS[relativeerrorplotBS!=0]
plot(log(h[-1]),log(relativeerrorplotBS))
# c.2)
reBS <- seq(10,1000,10)
for(i in 2:length(reBS)){
  res <- simulatedDeltaHedgingBS(T=1,K=1,S0=1,r=0.02,sigma=0.18,nHedge=reBS[i],nSteps=1000,nPaths=10000,offline=BlackScholesGrid)[[6]]
  error <- res[,ncol(res)]
  rmse <- sqrt(mean(error^2))
  reBS[i] <- rmse/1
  print(i)
}
plot(reBS[-1],xlab='Number of Hedges',ylab='Relative Error')
reBS

