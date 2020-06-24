## Supplemental materials for *``Integrative statistical methods for exposure mixtures and health''*

### Code for the SSVS model

The code below simulated one dataset from the first simulation scenario and fits the shared SSVS model.

### Set-up


```r
dir <- "S:\\Documents\\AuxMixtures\\Demo\\"
source(paste0(dir,"functions.R"))
source(paste0(dir,"NB_SSVS_Aux.R"))

# Offset
iters  <- 6000  # Number of MCMC iterations
burn   <- 1000  # Length of burn-in
nsims  <- 1     # Number of simulated datasets
b_type <- 1     # True beta type (either 1 or 2)
rhoX   <- 0.25  # Autocorrelation of X
rhoZ   <- 0     # Autocorrelation of Z

n     <- 500    # Sample size
p     <- 40     # Number of exposure variables
q     <- 5      # Number of auxiliary variables
N     <- 10     # Poisson offset term
r     <- 2      # Overdispersion

 pos_part <- function(x){ifelse(x>0,x,0)}
 true_beta <- function(Z,b_type){
    if(b_type==1){
       b <- 0.1*pos_part(Z[,4]+Z[,5])
    }
    if(b_type==2){
       b <- 0.05*pos_part(Z[,1:5]%*%(1-1:5/5)+1)
    }
 return(b)}
```

### Begin the simulation


```r
 MSE            <- rep(0,nsims)
 for(sim in 1:nsims){    

    set.seed(0820*sim)

    # Generate data

    X <- matrix(rnorm(n*p),n,p)
    for(j in 2:p){
       X[,j] <- rhoX*X[,j-1] + sqrt(1-rhoX^2)*X[,j]
    }
    Z <- matrix(rnorm(p*q),p,q)
    for(l in 2:q){
       Z[,l] <- rhoZ*Z[,l-1] + sqrt(1-rhoZ^2)*Z[,l]
    }
    beta <- true_beta(Z,b_type)
    Y    <- rlog_like(r,eta=X%*%beta,N=N)

    # Fit models
    fit <- NB_SSVS(Y,X,Z,Z,N,iters=iters,burn=burn,thin=2,common_b=TRUE)

    MSE[sim] <- mean((colMeans(fit$beta[burn:iters,])-beta)^2)
}
```



### Summarize the results


```r
    fit$min
```

```
##  elapsed 
## 2.736167
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

