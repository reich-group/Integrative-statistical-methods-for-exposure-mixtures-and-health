#########################################################
#
# Code for the multiple regression model with auxiliary data
#
#########################################################


#####################################
#
# Y[i]     ~ NB(r,N[i]*exp(eta[i]))
# eta      = alpha + X%*%beta
# beta     ~ Normal(Z%*%b,sigma^2)
# b        ~ Normal(0,tau^2D)     
# 
# alpha    ~ Normal(0,pri_sd_a^2)
# log(r)   ~ Normal(pre_mn_r,pri_sd_r^2)
# tau^2    ~ InvGamma(eps,eps)
# sigma^2  ~ InvGamma(eps,eps)
#
#####################################

NB_PCR<-function(Y,X1,X2,N,
                 iters=10000,burn=1000,update=1000,thin=1,
                 can_sd=0.1,
                 pri_mn_r=0,pri_sd_r=10,
                 pri_sd_a=10,eps=0.1){
  tick <- proc.time()[3]

  n    <- length(Y)
  p1   <- ncol(X1)
  p2   <- ncol(X2)
  p    <- p1+p2
  miss <- which(is.na(Y))
  nm   <- length(miss)

  #Initial values:
  alpha <- log(mean(Y,na.rm=TRUE)/mean(N,na.rm=TRUE))
  beta1 <- rep(0,p1)
  beta2 <- rep(0,p2)
  r     <- 10
  sig1  <- 0.1
  sig2  <- 0.1

  # Keep track of the samples   

   keep_alpha    <- rep(0,iters)
   keep_beta1    <- matrix(0,iters,p1)
   keep_beta2    <- matrix(0,iters,p2)
   keep_sig      <- matrix(0,iters,2)
   keep_r        <- rep(0,iters)

   acc   <- att <- MH <- rep(can_sd,p+2)
   eta   <- alpha + X1%*%beta1 + X2%*%beta2
   Y_hat  <- eta_hat <- rep(0,n)

   for(iter in 2:iters){

     if(nm>0){
       Y[miss] <- rlog_like(r,eta[miss],N[miss])
     }
     curll <- log_like(Y,r,eta,N)

     for(thinthin in 1:thin){

      #Update alpha using MH sampling:
      att[p+1] <- att[p+1] + 1 
      cana     <- rnorm(1,alpha,MH[p+1])
      caneta   <- eta - alpha + cana
      canll    <- log_like(Y,r,caneta,N)
      R        <- canll-curll+
                  dnorm(cana,0,pri_sd_a,log=TRUE)-
                  dnorm(alpha,0,pri_sd_a,log=TRUE)
      if(log(runif(1)) < R){
       acc[p+1] <- acc[p+1] +1 
       alpha    <- cana
       eta      <- caneta
       curll    <- canll
      }


      # Update beta using MH sampling:
      for(j in 1:p1){
       att[j]   <- att[j] + 1 
       canb     <- rnorm(1,beta1[j],MH[j])
       caneta   <- eta + X1[,j]*(canb-beta1[j])
       canll    <- log_like(Y,r,caneta,N)
       R        <- canll-curll+
                   dnorm(    canb,0,sig1,log=TRUE)-
                   dnorm(beta1[j],0,sig1,log=TRUE)
       if(log(runif(1)) < R){
         acc[j]   <- acc[j] + 1 
         beta1[j] <- canb
         eta      <- caneta
         curll    <- canll
       }
      }

      # Update beta using MH sampling:
      for(j in 1:p2){
       att[p1+j] <- att[p1+j] + 1 
       canb      <- rnorm(1,beta2[j],MH[p1+j])
       caneta    <- eta + X2[,j]*(canb-beta2[j])
       canll     <- log_like(Y,r,caneta,N)
       R         <- canll-curll+
                    dnorm(    canb,0,sig2,log=TRUE)-
                    dnorm(beta2[j],0,sig2,log=TRUE)
       if(log(runif(1)) < R){
         acc[p1+j]   <- acc[p1+j] + 1 
         beta2[j] <- canb
         eta      <- caneta
         curll    <- canll
       }
      }


      #Update r:
      att[p+2] <- att[p+2] + 1 
      canr     <- exp(rnorm(1,log(r),5*MH[p+2]))
      canll    <- log_like(Y,canr,eta,N)
      R        <- canll-curll+
                  dnorm(log(canr),pri_mn_r,pri_sd_r,log=TRUE)-
                  dnorm(log(r),pri_mn_r,pri_sd_r,log=TRUE)
      if(log(runif(1)) < R){
        acc[p+2] <- acc[p+2] + 1 
        r        <- canr
        curll    <- canll
      }

      sig1 <- newsd(beta1,sig1)
      sig2 <- newsd(beta2,sig2)


     } # end thin

     keep_beta1[iter,] <- beta1
     keep_beta2[iter,] <- beta2
     keep_alpha[iter]  <- alpha
     keep_r[iter]      <- r
     keep_sig[iter,]   <- c(sig1,sig2)

     if(iter > burn){
          Y_hat <-   Y_hat +   Y/(iters-burn)
        eta_hat <- eta_hat + eta/(iters-burn)
     }


     if(iter < burn){for(j in which(att>50)){
        if(acc[j]/att[j]<0.2){MH[j] <- 0.8*MH[j]}
        if(acc[j]/att[j]>0.5){MH[j] <- 1.2*MH[j]}      
        acc[j] <- att[j] <- 0.5
     }}

   }
   tock <- proc.time()[3]

list(Y_hat=Y_hat,eta_hat=eta_hat,
     beta1=keep_beta1,beta2=keep_beta2,
     sig=keep_sig,r=keep_r,alpha=keep_alpha,
     min=(tock-tick)/60,acc_rate=acc/att)}
