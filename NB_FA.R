#########################################################
#
# Code for the factor analysis model with auxiliary data
#
#########################################################

#####################################
#
# Y[i]     ~ NB(r,N[i]*exp(eta[i]))
# eta      = alpha + theta%*%gamma
# X[i]     ~ Normal(A%*%theta[i,],taux)
# A[i,j]   ~ Normal(b*G[i,j],taua)
#
# b        ~ Normal(0,pri_sd)      
# gamma[j] ~ Normal(0,pri_sd)      
# alpha    ~ Normal(0,pri_sd^2)
# log(r)   ~ Normal(pri_mn_r,pri_sd_r^2)
# taux     ~ InvGamma(eps,eps)
# taua     ~ InvGamma(eps,eps)
#
#####################################


NB_FA<-function(Y,X,G,N,
                iters=10000,burn=1000,update=1000,thin=1,
                can_sd=0.1,eps=0.1,pri_sd=100,pri_sd_b=100,pri_mn_r=0,pri_sd_r=10){
  tick <- proc.time()[3]

  # Bookeeping

  n    <- length(Y)
  p    <- ncol(X)
  d    <- ncol(G)
  miss <- which(is.na(Y))
  nm   <- length(miss)

  #Initial values:

  alpha <- log(mean(Y/N,na.rm=TRUE))
  theta <- matrix(0,n,d)
  gamma <- rep(0,d)
  r     <- 10
  b     <- 1
  taux  <- rep(10,p)
  taua  <- 1
  A     <- G

  # Keep track of the samples   

   keep_alpha     <- rep(0,iters)
   keep_beta      <- matrix(0,iters,p)
   keep_gamma     <- matrix(0,iters,d)
   keep_b         <- rep(0,iters)
   keep_sig       <- matrix(0,iters,p+1)
   keep_r         <- rep(0,iters)
   keep_A         <- array(0,c(iters,p,d))
   keep_gamma[1,] <- gamma
   keep_b[1]     <- b
   keep_alpha[1] <- alpha
   keep_r[1]     <- r
   keep_sig[1,]  <- 1/sqrt(c(taua,taux))
   keep_A[1,,]    <- A

   acc   <- att <- MH <- rep(can_sd,p+d+10)
   Pt    <- t(chol(solve(t(X)%*%X)))/sqrt(p)
   eta   <- alpha + theta%*%gamma

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
                  dnorm(cana,0,pri_sd,log=TRUE)-
                  dnorm(alpha,0,pri_sd,log=TRUE)
      if(log(runif(1)) < R){
       acc[p+1] <- acc[p+1] +1 
       alpha    <- cana
       eta      <- caneta
       curll    <- canll
      }


      # Update beta using MH sampling:
      for(j in 1:d){
       att[j]   <- att[j] + 1 
       cang     <- rnorm(1,gamma[j],MH[j])
       caneta   <- eta + theta[,j]*(cang-gamma[j])
       canll    <- log_like(Y,r,caneta,N)
       R        <- canll-curll+
                   dnorm(    cang,0,pri_sd,log=TRUE)-
                   dnorm(gamma[j],0,pri_sd,log=TRUE)
       if(log(runif(1)) < R){
         acc[j]   <- acc[j] + 1 
         gamma[j] <- cang
         eta      <- caneta
         curll    <- canll
       }
      }

     #Update theta:
     DDD    <- diag(taux)
     VVV    <- solve(t(A)%*%DDD%*%A + diag(d))
     PPP    <- t(chol(VVV))
     MMM    <- VVV%*%t(A)%*%DDD%*%t(X) 
     for(reps in 1:1){
        can    <- MMM + PPP%*%matrix(rnorm(n*d),d,n)
        can    <- t(can)
        caneta <- alpha + can%*%gamma
        ll1    <- log_like_nosum(Y,r,caneta,N)
        ll2    <- log_like_nosum(Y,r,eta,N)
        keep   <- log(runif(n)) < ll1-ll2
        if(sum(keep)>1){
         theta[keep,] <- can[keep,]
        }
      }
      eta          <- alpha + theta%*%gamma
      curll        <- log_like(Y,r,eta,N)


      #Update A:
      theta2 <- t(theta)%*%theta
      for(j in 1:p){
         VVV   <- solve(taux[j]*theta2 + taua*diag(d))
         PPP   <- t(chol(VVV))
         MMM   <- taux[j]*t(theta)%*%X[,j] + taua*b*G[j,]
         A[j,] <- VVV%*%MMM+PPP%*%rnorm(d)
      }

      #Update hyperparameters:
      VVV   <- taua*sum(G^2)+1/pri_sd_b^2
      MMM   <- taua*sum(G*A)
      b     <- rtnorm(1,MMM/VVV,1/sqrt(VVV),0,1000)

      taua <- 1/newsd(as.vector(A-b*G),1/sqrt(taua))^2
      for(j in 1:p){
         taux[j] <- 1/newsd(X[,j]-theta%*%A[j,],1/sqrt(taux[j]))^2
      }

     } # end thin

     AD                <- sweep(t(A),2,taux,"*")
     keep_beta[iter,]  <- t(AD)%*%solve(AD%*%A+diag(d))%*%gamma
     keep_gamma[iter,] <- gamma
     keep_A[iter,,]    <- A
     keep_b[iter]      <- b
     keep_alpha[iter]  <- alpha
     keep_r[iter]      <- r
     keep_sig[iter,]   <- 1/sqrt(c(taua,taux))

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
     beta=keep_beta,gamma=keep_gamma,
     b=keep_b,sig=keep_sig,r=keep_r,
     alpha=keep_alpha,A=keep_A,
     min=(tock-tick)/60,acc_rate=acc/att)}
