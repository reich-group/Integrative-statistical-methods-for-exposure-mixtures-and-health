expit <- function(x){ifelse(x>10,1/(1+exp(-10)),1/(1+exp(-x)))}

link <- function(x){ifelse(x<10,exp(x),exp(10))}

log_like<-function(Y,r,eta,N=1){
  sum(dnbinom(Y,r,mu=N*link(eta),log=TRUE))
}
log_like<-function(Y,r,eta,N=1){
  sum(dpois(Y,N*link(eta),log=TRUE))
}

log_like_nosum<-function(Y,r,eta,N=1){
  dpois(Y,N*link(eta),log=TRUE)
}


rlog_like<-function(r,eta,N=1){
  rnbinom(length(eta),r,mu=N*link(eta))
}
rlog_like<-function(r,eta,N=1){
  rpois(length(eta),N*link(eta))
}

rtnorm<-function(n,m,s,l,u){
  lo <- pnorm(l,m,s)
  hi <- pnorm(u,m,s)
  Z  <- runif(n,lo,hi)
  Y  <- qnorm(Z,m,s)
return(Y)}

update_s <- function(Y,sig,scale=1/qt(0.99,df=1),MH=1){
    Y    <- as.vector(Y)
    A    <- sd(Y)
    B    <- 0.71*MH*A/sqrt(length(Y))
    can  <- abs(A+B*rt(1,df=2))
    R    <- sum(dnorm(Y,0,can,log=TRUE))-
            sum(dnorm(Y,0,sig,log=TRUE))+
            dt(can/scale,df=1,log=TRUE)-
            dt(sig/scale,df=1,log=TRUE)+
            dt((sig-A)/B,df=2,log=TRUE)-
            dt((can-A)/B,df=2,log=TRUE)
    if(!is.na(R)){if(can>0){if(can<Inf){
      sig  <- ifelse(log(runif(1))<R,can,sig)
    }}}
return(sig)}


newsd <- function(Y,sig,scale=1/qt(0.99,df=1),MH=1){
    Y    <- as.vector(Y)
    A    <- MH*length(Y)/2+0.1
    B    <- MH*sum(Y^2)/2+0.1
    can  <- 1/sqrt(rgamma(1,A,B))
    R    <- sum(dnorm(Y,0,can,log=TRUE))-
            sum(dnorm(Y,0,sig,log=TRUE))+
            dt(can/scale,df=1,log=TRUE)-
            dt(sig/scale,df=1,log=TRUE)+
            (dgamma(1/sig^2,A,B,log=TRUE)-3*log(sig))-
            (dgamma(1/can^2,A,B,log=TRUE)-3*log(can))
    if(!is.na(R)){
      sig  <- ifelse(log(runif(1))<R,can,sig)
    }
return(sig)}