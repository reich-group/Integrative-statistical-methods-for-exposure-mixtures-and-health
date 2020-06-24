#########################################################
#
# Code for the SSVS model with auxiliary data
#
#########################################################


#####################################
#
# Y[i]     ~ NB(r,N[i]*exp(eta[i]))
# eta      = a1 + X%*%beta
# beta     = bin[j]*cont[j]
# bin[j]   ~ Bern(expit(a1+a2*Z%*%b))
# cont[j]  ~ N(Z%*%b,sigc)
# b        ~ Normal(0,sigb)     
# 
# a1,a2,a3 ~ Normal(0,pri_sd_a^2)
# alpha    ~ Normal(0,pri_sd_a^2)
# log(r)   ~ HalfCauchy(1)
# sigma^2  ~ HalfCauchy(1)
# tau^2    ~ HalfCauchy(1)
#
#####################################

NB_SSVS<-function(Y,X,Z1,Z2,N=1,
                  common_b=FALSE,active_b1=TRUE,active_b2=TRUE,
                  iters=10000,burn=1000,update=1000000,thin=1,
                  can_sd=0.1,
                  pri_mn_r=0,pri_sd_r=10,
                  pri_sd_a=1){
  tick <- proc.time()[3]

  n    <- length(Y)
  p    <- ncol(X)
  q1   <- ncol(Z1)
  q2   <- ncol(Z2)
  q    <- q1+q2
  miss <- which(is.na(Y))
  nm   <- length(miss)

  if(common_b & ncol(Z1)!=ncol(Z2)){print("dim(Z1) needs to be dim(Z2)")}

  #Initial values:
  a1    <- log(mean((Y+.1)/(N+.01),na.rm=TRUE))
  a2    <- 0
  rho   <- 0
  bin   <- rep(0,p)
  cont  <- rep(0,p)
  b1    <- rep(0,q1)
  b2    <- rep(0,q2)
  r     <- 10
  sigb1 <- 0.1
  sigb2 <- 0.1
  sigc  <- 0.1

  # Keep track of the samples   

   keep_a        <- matrix(0,iters,2)
   keep_beta     <- matrix(0,iters,p)
   keep_b1       <- matrix(0,iters,q1)
   keep_b2       <- matrix(0,iters,q2)
   keep_sig      <- matrix(0,iters,3)
   keep_r        <- rep(0,iters)
   keep_rho      <- rep(0,iters)

   acc   <- att <- MH <- rep(can_sd,p+q+10)
   beta  <- bin*cont
   eta   <- a1 + X%*%beta
   Zb1   <- Z1%*%b1
   Zb2   <- Z2%*%b2

   Y_hat <- eta_hat <- 0

   for(iter in 1:iters){

     if(nm>0){
       Y[miss] <- rlog_like(r,eta[miss],N[miss])
     }
     curll <- log_like(Y,r,eta,N)

     for(thinthin in 1:thin){

      #Update the a's
      # a1
      att[p+q+1] <- att[p+q+1] + 1 
      cana     <- rnorm(1,a1,MH[p+q+1])
      caneta   <- eta - a1 + cana
      canll    <- log_like(Y,r,caneta,N)
      R        <- canll-curll+
                  dnorm(cana,0,pri_sd_a,log=TRUE)-
                  dnorm(a1,0,pri_sd_a,log=TRUE)
      if(!is.na(R)){if(log(runif(1)) < R){
       acc[p+q+1] <- acc[p+q+1] +1 
       a1         <- cana
       eta        <- caneta
       curll      <- canll
      }}

      # a2
       att[p+q+2] <- att[p+q+2] + 1 
       cana       <- rnorm(1,a2,MH[p+q+2])
       R          <- sum(dbinom(bin,1,expit(cana+rho*Zb1+Zb2),log=TRUE))-
                     sum(dbinom(bin,1,expit(a2+rho*Zb1+Zb2),log=TRUE))+
                     dnorm(cana,0,pri_sd_a,log=TRUE)-
                     dnorm(a2,0,pri_sd_a,log=TRUE)
       if(!is.na(R)){if(log(runif(1)) < R){
         acc[p+q+2] <- acc[p+q+2] + 1 
         a2         <- cana
       }}

      # Update beta using MH sampling:
      for(j in 1:p){
       inprob  <- expit(a2+rho*Zb1[j]+Zb2[j])
       eta     <- eta - X[,j]*beta[j]
       canll0  <- log_like(Y,r,eta,N)
       canll1  <- log_like(Y,r,eta+X[,j]*cont[j],N)
       R       <- expit(canll1-canll0+log(inprob)-log(1-inprob))
       bin[j]  <- rbinom(1,1,R)    

       if(bin[j]==0){cont[j]<-rnorm(1,Zb1[j],sigc)}
       if(bin[j]==1){
         att[j]   <- att[j] + 1 
         canc     <- rnorm(1,cont[j],MH[j])
         canll1   <- log_like(Y,r,eta+X[,j]*canc,N)
         canll2   <- log_like(Y,r,eta+X[,j]*cont[j],N)
         R        <- canll1-canll2+
                     dnorm(canc,Zb1[j],sigc,log=TRUE)-
                     dnorm(cont[j],Zb1[j],sigc,log=TRUE)
         if(!is.na(R)){if(log(runif(1)) < R){
           acc[j]  <- acc[j] + 1 
           cont[j] <- canc
         }}
       } 

       beta[j] <- bin[j]*cont[j]
       eta     <- eta + X[,j]*beta[j]
      }
      curll   <- log_like(Y,r,eta,N)


      #Update b1:
      if(!common_b & active_b1){for(j in 1:q1){
       att[p+j] <- att[p+j] + 1
       canb     <- rnorm(1,b1[j],MH[p+j])
       canZb    <- Zb1+Z1[,j]*(canb-b1[j])
       R        <- sum(dnorm(cont,canZb,sigc,log=TRUE))-
                   sum(dnorm(cont,Zb1,sigc,log=TRUE))+
                   sum(dbinom(bin,1,expit(a2+rho*canZb+Zb2),log=TRUE))-
                   sum(dbinom(bin,1,expit(a2+rho*Zb1+Zb2),log=TRUE))+
                   dnorm(canb,0,sigb1,log=TRUE)-
                   dnorm(b1[j],0,sigb1,log=TRUE)
       if(!is.na(R)){if(log(runif(1)) < R){
         acc[p+j] <- acc[p+j] + 1 
         b1[j]     <- canb
         Zb1       <- canZb
       }}
      }}
      

      #Update b2:
      if(!common_b & active_b2){for(j in 1:q2){
       att[p+q1+j] <- att[p+q1+j] + 1
       canb        <- rnorm(1,b2[j],MH[p+q1+j])
       canZb       <- Zb2+Z2[,j]*(canb-b2[j])
       R           <- sum(dbinom(bin,1,expit(a2+rho*Zb1+canZb),log=TRUE))-
                      sum(dbinom(bin,1,expit(a2+rho*Zb1+Zb2),log=TRUE))+
                      dnorm(canb,0,sigb2,log=TRUE)-
                      dnorm(b2[j],0,sigb2,log=TRUE)
       if(!is.na(R)){if(log(runif(1)) < R){
         acc[p+q1+j] <- acc[p+q1+j] + 1 
         b2[j]       <- canb
         Zb2         <- canZb
       }}
      }}
      
      if(common_b){
       b2  <- 0*b2
       Zb2 <- 0*Zb2
       for(j in 1:q1){
         att[p+j] <- att[p+j] + 1
         canb     <- rnorm(1,b1[j],MH[p+j])
         canZb    <- Zb1+Z1[,j]*(canb-b1[j])
         R        <- sum(dnorm(cont,canZb,sigc,log=TRUE))-
                     sum(dnorm(cont,Zb1,sigc,log=TRUE))+
                     sum(dbinom(bin,1,expit(a2+rho*canZb+Zb2),log=TRUE))-
                     sum(dbinom(bin,1,expit(a2+rho*Zb1+Zb2),log=TRUE))+
                     dnorm(canb,0,sigb1,log=TRUE)-
                     dnorm(b1[j],0,sigb1,log=TRUE)
         if(!is.na(R)){if(log(runif(1)) < R){
           acc[p+j] <- acc[p+j] + 1 
           b1[j]     <- canb
           Zb1       <- canZb
         }}
       }

      #Update rho:
       att[p+q+4] <- att[p+q+4] + 1
       can        <- rnorm(1,rho,MH[p+q+4])
       R          <- sum(dbinom(bin,1,expit(a2+can*Zb1+Zb2),log=TRUE))-
                     sum(dbinom(bin,1,expit(a2+rho*Zb1+Zb2),log=TRUE))+
                     dt(can,df=1,log=TRUE)-
                     dt(rho,df=1,log=TRUE)
       if(!is.na(R)){if(log(runif(1)) < R){
         acc[p+q+4] <- acc[p+q+4] + 1 
         rho        <- can
       }}

     }
      

      #Update r:
      att[p+q+3] <- att[p+q+3] + 1 
      canr     <- exp(rnorm(1,log(r),5*MH[p+q+3]))
      canll    <- log_like(Y,canr,eta,N)
      R        <- canll-curll+
                  dnorm(log(canr),pri_mn_r,pri_sd_r,log=TRUE)-
                  dnorm(log(r),pri_mn_r,pri_sd_r,log=TRUE)
      if(!is.na(R)){if(log(runif(1)) < R){
        acc[p+q+3] <- acc[p+q+3] + 1 
        r        <- canr
        curll    <- canll
      }}


      # Update sds  
      sigc  <- newsd(cont-Zb1,sigc)
      sigb1 <- newsd(b1,sigb1)
      sigb2 <- newsd(b2,sigb2)

     } # end thin

     keep_beta[iter,] <- beta
     keep_b1[iter,]   <- b1
     keep_b2[iter,]   <- b2
     keep_a[iter,]    <- c(a1,a2)
     keep_r[iter]     <- r
     keep_rho[iter]   <- rho
     keep_sig[iter,]  <- c(sigc,sigb1,sigb2)

     if(iter > burn){
          Y_hat <-   Y_hat +   Y/(iters-burn)
        eta_hat <- eta_hat + eta/(iters-burn)
     }

     if(iter%%update==0){
#       print(paste("Done with",iter,"of",iters,"iterations"))
       plot(Z[,4],beta,main=paste("Done with",iter,"of",iters,"iterations"))
     }

     if(iter < burn){for(j in which(att>50)){
        if(acc[j]/att[j]<0.2){MH[j] <- 0.8*MH[j]}
        if(acc[j]/att[j]>0.5){MH[j] <- 1.2*MH[j]}      
        acc[j] <- att[j] <- 0.5
     }}

   }
   tock <- proc.time()[3]

list(Y_hat=Y_hat,eta_hat=eta_hat,
     beta=keep_beta,b1=keep_b1,b2=keep_b2,
     sig=keep_sig,r=keep_r,a=keep_a,rho=keep_rho,
     min=(tock-tick)/60,acc_rate=acc/att)}
