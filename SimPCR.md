## Supplemental materials for *``Integrative statistical methods for exposure mixtures and health''*

### Code for the PCR model

The code below simulated one dataset from the first simulation scenario and fits the PCR model with auxiliary data (PCR-XZ).

### Set-up


```r
dir <- "S:\\Documents\\AuxMixtures\\Demo\\"

source(paste0(dir,"functions.R"))
source(paste0(dir,"NB_PCR.R"))


iters  <- 6000 # Number of MCMC iterations
burn   <- 1000 # Length of burn-in
nsims  <- 1    # Number of simulated datasets
b_type <- 1    # True beta type (either 1 or 2)
rhoX   <- 0.5  # Autocorrelation of X
rhoZ   <- 0.5  # Autocorrelation of Z

n     <- 500   # Sample size
p     <- 40    # Number of exposure variables
q     <- 20    # Number of auxiliary variables
N     <- 10    # Poisson offset term
r     <- 2     # Overdispersion

true_beta <- function(X,Z,b_type){
    Sx <- eigen(cov(X))$vec[,1:3]
    Sz <- eigen(cov(t(Z)))$vec[,1:3]
    r  <- ifelse(b_type==1,0.1,0.9)
    b  <- 0.1*((1-r)*rowSums(Sx) + r*rowSums(Sz))
return(b)}
```

### Begin the simulation


```r
MSE            <- rep(0,nsims)

 for(sim in 1:nsims){    

    set.seed(0820*sim)

    # Generate data

    X <- matrix(rnorm(n*p),n,p)
    for(k in 0:4){
       RE    <- rnorm(n)
       for(j in which(1:p%%5==k)){
         X[,j] <- sqrt(rhoX)*RE + sqrt(1-rhoX)*X[,j]
       } 
     }
    Z <- matrix(rnorm(p*q),p,q)
    for(j in 2:p){
       Z[j,] <- rhoZ*Z[j-1,] + sqrt(1-rhoZ^2)*Z[j,]
    }

    beta <- true_beta(X,Z,b_type)
    Y    <- rlog_like(r,eta=X%*%beta,N=N)

    # Do PCA for the prior

    Ex <- eigen(cor(X))
    dx <- sum(cumsum(Ex$val)/sum(Ex$val)<0.9)+1
    dx <- ifelse(dx<2,2,dx)
    dx <- ifelse(dx>q,q,dx)
    Gx <- sweep(Ex$vec[,1:dx],2,sqrt(Ex$val[1:dx]),"*")

    Ez <- eigen(cor(t(Z)))
    dz <- sum(cumsum(Ez$val)/sum(Ez$val)<0.9)+1
    dz <- ifelse(dz<2,2,dz)
    dz <- ifelse(dz>q,q,dz)
    Gz <- sweep(Ez$vec[,1:dz],2,sqrt(Ez$val[1:dz]),"*")

    # Fit the model
    fit      <- NB_PCR(Y,X%*%Gx,X%*%Gz,N,iters=iters,burn=burn,thin=2,update=2*iters)
    fit$beta <- fit$beta1%*%t(Gx) + fit$beta2%*%t(Gz)

    MSE[sim] <- mean((colMeans(fit$beta[burn:iters,])-beta)^2)

}
```



### Summarize the results


```r
    fit$min
```

```
## elapsed 
##   0.763
```

```r
    boxplot(fit$beta,outline=FALSE)
    lines(beta,lwd=2,col=2)
```

![plot of chunk sum](figure/sum-1.png)

```r
    plot(fit$beta[,1],type="l")
```

![plot of chunk sum](figure/sum-2.png)


