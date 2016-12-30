#beta = 0 no discrepency

####################################################################################

######################### Sorbent Emulation MCMC ###################################

####################################################################################


###############################
##### Utility Functions #######
###############################


rwish <- function(v, S){
  S <- as.matrix(S)
  if (v < nrow(S)){
      stop(message = "v is less than the dimension of S in rwish().\n")
  }
  p <- nrow(S)
  CC <- chol(S)
  Z <- matrix(0, p, p)
  diag(Z) <- sqrt(rchisq(p, v:(v - p + 1)))
  if(p > 1){
    pseq <- 1:(p - 1)
    Z[rep(p*pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p*(p-1)/2)
  }
  return(crossprod(Z%*%CC))
}



riwish <- function(v, S){
  return(solve(rwish(v, solve(S))))
}


pos.pow.transform <- function(x,a=.4, knot=.001, c=.001){
  b <- (knot+c)^a
  const <- b*(1-log(b))
  ans <- (x+c)^a*(x>knot) + (a*b*log(x+c)+const)*(x<=knot)
  return(ans)
}



pos.pow.inv.transform <- function(y,a=.4, knot=.001, c=0.001){
  b <- (knot+c)^a
  const <- b*(1-log(b))
  x1 <- y^(1/a)-c
  x2 <- exp((y-const)/(a*b))-c
  ans <- ifelse(x2<knot,x2,x1)
  return(ans)
}



nextpow2 <- function(x) {

if (!(length(x) == 1))
stop("nextpow2: x must be a scalar or a vector")

t = length(x)
if (t > 1)
x = t

ceiling(log2(abs(x)))
}

var_of_mean <-function(x,mean_x) {

    # Georgios Karagiannis @ PNNL

    n = length(x)
    
    nn = 2^nextpow2(2*n)+1
    
    y = rep(0, nn)
    y[1:min(nn,n)] = x[1:min(nn,n)] - mean_x
    
    y = fft( y, inverse = FALSE ) 
    
    y = Re(y)^2+Im(y)^2 
    
    y = fft(y, inverse = TRUE)/nn
    
    y = Re(y[1:n]) 
    
    v = y[1]/n 
    
    y = y/y[1] 
    
    sumT = -1/3 
    i = 0 
    test = 0 
    while (i<=n & test==0)
        {
        i = i+1 
        sumT = sumT + y[i] - 1/6 
        if (sumT < 0.0) 
            { 
             test = 1
            }
        }
    var_mean = 2.0*(sumT+(i-1)/6)*v/n 
    
    var_mean
}



logit <- function(x, xlim=c(0, 1)){
  xt <- (x-xlim[1])/(xlim[2]-xlim[1])
  ans <- log(xt/(1-xt))
  return(ans)
}


inv.logit <- function(y, xlim=c(0,1)){
  xt <- exp(y)/(1+exp(y))
  x <- xt*(xlim[2]-xlim[1])+xlim[1]
  return(x)
}


index <- function(m,n){
  if(m<=n) return(m:n)
  else return(numeric(0))
}

which.equal <- function(x, y){

  n <- length(x)
  ans <- rep(0,n)
  for(i in 1:n){
    ans[i] <- any(approx.equal(y,x[i]))
  }
  return(as.logical(ans))
}

approx.equal <- function(x, y, tol=1E-9){

  return(abs(x - y) < tol)
}

robust.optimize <- function(fn, par.lim, npar=5, rel.tol=1E-2, ...){

 ## Conduct grid search on par
  par.min <- par.lim[1]
  par.max <- par.lim[2]

  inc <- (par.max - par.min)/(npar-1)
  par.vec <- seq(par.min, par.max, inc)
  par.old <- par.vec
  pen.best <- Inf
  pen.prev <- Inf

  repeat{
    npar.now <- length(par.vec)
    for(j in 1:npar.now){
      par <- par.vec[j]

#cat("\npar =",par)

      pen <- fn(par,...)

#cat(",   pen =",pen)


      if(pen < pen.best){
        pen.best <- pen
        par.best <- par
      }
    }

    if(inc/(abs(par.best)+1E-6) <= rel.tol)
      break
    else
      pen.prev <- pen.best

   ## create par.vec for next pass
    par.min <- par.best - floor(npar/2)/2*inc
    par.min <- max(par.min, par.lim[1])  
    par.max <- par.best + floor(npar/2)/2*inc
    par.max <- min(par.max, par.lim[2])  
    inc <- inc/2
    par.vec <- seq(par.min, par.max, inc) 
    ind <- which.equal(par.vec, par.old)
    par.vec <- par.vec[!ind]
    if(length(par.vec)==0)
      break
    par.old <- c(par.old, par.vec)
  }

  return(list(par=par.best, pen=pen.best))
}



## Generates a matrix whose columns are random samples from 1:n
get.obs.ind <- function(n, nfolds=5, seed=220){

  replace.seed <- T
  if(missing(seed))
    replace.seed <- F

  if(replace.seed){
   ## set seed to specified value
    if(!any(ls(name='.GlobalEnv', all.names=T)=='.Random.seed')){
      set.seed(1)
    }
    save.seed <- .Random.seed
    set.seed(seed)  
  }

  perm <- sample(1:n, n)
  n.cv <- rep(floor(n/nfolds),nfolds)
  rem <- n - n.cv[1]*nfolds
  n.cv[index(1,rem)] <- n.cv[index(1,rem)]+1
  obs.ind <- list()
  ind2 <- 0
  
  for(i in 1:nfolds){
    ind1 <- ind2+1
    ind2 <- ind2+n.cv[i]
    obs.ind[[i]] <- perm[ind1:ind2]
  }
  if(replace.seed){
   ## restore random seed to previous value
    .Random.seed <<- save.seed
  }
  return(obs.ind)
}



boscocov.mat <- function(x1, x2 = NULL){

    if(is.null(x2)){
      x2 <- x1
    }
    n1 <- length(x1)
    n2 <- length(x2)
    x1.ext <- rep(x1, times = n2)
    x2.ext <- rep(x2, each = n1)
    ans <- -1/24*B4(abs(x1.ext-x2.ext))
    Cov.mat <- matrix(ans, nrow = n1, ncol = n2)
    return(Cov.mat)
}


boscocov2.mat <- function(x1, x2 = NULL, sigma2){

    if(is.null(x2)){
      x2 <- x1
    }
    n1 <- length(x1)
    n2 <- length(x2)
    x1.ext <- rep(x1, times = n2)
    x2.ext <- rep(x2, each = n1)
    ans <- sigma2[1]+sigma2[2]*B1(x1.ext)*B1(x2.ext)+sigma2[3]*B2(x1.ext)*B2(x2.ext)-sigma2[4]/24*B4(abs(x1.ext-x2.ext))
    Cov.mat <- matrix(ans, nrow = n1, ncol = n2)
    return(Cov.mat)
}




ginv.gp <- function (X, eps = 1e-12){
    if(any(X==Inf))
      return(list(inv = diag(nrow(X)), log.det = Inf))
    eig.X <- eigen(X, symmetric = TRUE)
    P <- eig.X[[2]]
    lambda <- eig.X[[1]]
    ind <- lambda > eps
    lambda.inv <- lambda.sqrt <- lambda
    lambda.inv[ind] <- 1/lambda[ind]
    lambda.sqrt[ind] <- sqrt(lambda[ind])
    lambda.inv[!ind] <- lambda.sqrt[!ind] <- 0
    Q <- P%*%sqrt(diag(lambda.inv, nrow=length(lambda)))
    inv <- Q%*%t(Q)
    sqrt <- P %*% diag(lambda.sqrt, nrow=length(lambda)) %*% t(P)
    sqrt.inv <- P %*% diag(sqrt(lambda.inv), nrow=length(lambda)) %*% t(P)
    if(all(lambda>0))
      log.det <- sum(log(lambda))
    else
      log.det <- -Inf
    return(list(inv = inv, log.det = log.det, sqrt=sqrt, sqrt.inv=sqrt.inv))
}



rmvnorm <- function(mu=0, Sigma, S.sqrt=NA){

  if(is.na(S.sqrt[1])){
    ans <- ginv.gp(Sigma)
    S.sqrt <- ans$sqrt

#S.sqrt..<<-S.sqrt

  }
  if(length(mu)==1) 
    mu <- rep(mu,nrow(S.sqrt))
  x <- as.numeric(S.sqrt%*%rnorm(length(mu),0,1) + mu)
  return(x)
}


dmvnorm <- function(x, mu, Sigma, S.inv=NULL, log.det=NULL){

  if(is.null(S.inv[1]) || is.null(log.det)){
    ans <- ginv.gp(Sigma)
    S.inv <- ans$inv
    log.det <- ans$log.det
  }
  dens <- -length(x)/2*log(2*pi)-.5*log.det - 0.5*as.numeric(t(x-mu)%*%S.inv%*%(x-mu))
  return(dens)
}


get.comp2vars <- function(p, order=2, include2=NULL, include3=NULL, include4=matrix(0,0,4)){

  if(order==2)
    comp2vars <- matrix(NA,p+choose(p,2),2)
  if(order==3)
    comp2vars <- matrix(NA,p+choose(p,2)+choose(p,3),3)
  if(order==4)
    comp2vars <- matrix(NA,p+choose(p,2)+choose(p,3)+choose(p,4),4)
  for(j in 1:p)
    comp2vars[j,1] <- j
  comp.now <- p+1
  for(j in 1:(p-1)){
    for(k in (j+1):p){
      if(is.null(include2) || any(j==include2[,1] & k==include2[,2])){
        comp2vars[comp.now,1:2] <- c(j,k)
        comp.now <- comp.now+1
      }
    }
  }
  if(order>=3){
    for(j in 1:(p-2)){
      for(k in (j+1):(p-1)){
        for(l in (k+1):p){
          if(is.null(include3) || any(j==include3[,1] & k==include3[,2] & l==include3[,3])){
            comp2vars[comp.now,1:3] <- c(j,k,l)
            comp.now <- comp.now+1
          }
        }
      }
    }
  }
  if(order==4){
    for(j in 1:(p-3)){
      for(k in (j+1):(p-2)){
        for(l in (k+1):(p-1)){
          for(m in (l+1):p){
            if(is.null(include4) || any(j==include4[,1] & k==include4[,2] & l==include4[,3] & m==include4[,4])){
              comp2vars[comp.now,1:4] <- c(j,k,l,m)
              comp.now <- comp.now+1
            }
          }
        }
      }
    }
  }

  return(comp2vars[!is.na(comp2vars[,1]),])
}



### First Bernoulli Polynomial ###
B1 <- function(x)
  return(x-1/2)

### Second Bernoulli Polynomial ###
B2 <- function(x)
  return(x^2-x+1/6)

B3 <- function(x)
  return(x^3-3/2*x^2+x/2)

B4 <- function(x)
  return(x^4-2*x^3+x^2-1/30)

B5 <- function(x)
  return(x^5-5/2*x^4+5/3*x^3-1/6*x)

B6 <- function(x)
  return(x^6-3*x^5+5/2*x^4-1/2*x^2+1/42)


### Create discritized eigenfunctions ###
get.bss.basis <- function(t.grid){
  flag <- 0

 ##some logic to account for a singularity if 0 and 1 are both in t.grid
  if(min(t.grid)==0 && max(t.grid==1)){
    t.grid <- t.grid[-length(t.grid)]
    flag <- 1
  }
  N <- length(t.grid)
  Gamma <- boscocov.mat(t.grid)
  ans <- eigen(Gamma)
  lambda <- ans$values
  lambda[lambda<0] <- 0
  X <- ans$vectors

 ## Below I'm attaching the eigenvalue to the function from the get go,
 ##  so the coefficients (betas) can all have the same variance.
  Phi <- matrix(sqrt(lambda),N,N,byrow=TRUE)*X
  if(flag)
    Phi <- rbind(Phi,Phi[1,])
  return(Phi)
}




### Create discritized eigenfunctions ###
get.bss.basis2 <- function(t.grid){
  flag <- 0

  N <- length(t.grid)
  Gamma <- boscocov2.mat(t.grid,sigma2=c(0,1,1,1))
  ans <- eigen(Gamma)
  lambda <- ans$values
  lambda[lambda<0] <- 0
  X <- ans$vectors
  Phi <- matrix(sqrt(lambda),N,N,byrow=TRUE)*X
  return(Phi)
}



### Create discritized eigenfunctions ###
get.bss.basis3 <- function(t.grid){
  flag <- 0
 ##some logic to account for a singularity if 0 and 1 are both in t.grid
  if(min(t.grid)==0 && max(t.grid==1)){
    t.grid <- t.grid[-length(t.grid)]
    flag <- 1
  }

  N <- length(t.grid)
  Gamma <- boscocov2.mat(t.grid,sigma2=c(0,1,1,1))

#Gamma..<<-Gamma


  ans <- eigen(Gamma)
  lambda <- ans$values
  lambda[lambda<0] <- 0
  X <- ans$vectors
 ## Below I'm attaching the eigenvalue to the function from the get go,
 ##  so the coefficients (betas) can all have the same variance.
  Phi <- matrix(sqrt(lambda),N,N,byrow=TRUE)*X
  if(flag)
    Phi <- rbind(Phi,Phi[1,])
  return(Phi)
}



### "Evaluate" the eigen functions ###
get.main.X <- function(t,P,Phi,t.grid,mult=1,cat.var=0){

  n <- length(t)
  if(cat.var==0){
   ## for each value of t (in [0,1]), give the linear interpolation of 
   ##  the closest two rows of Phi 

    Phi <- Phi[,1:(P-2)]
    X <- matrix(1,n,P)
    X[,1] <- mult*B1(t)
    X[,2] <- mult*B2(t)
    for(i in 1:n){
      lo <- max(sum(t[i]>t.grid),1)
      diff.lo <- t[i]-t.grid[lo]
      diff.hi <- t.grid[lo+1]-t[i]
      w <- diff.hi/(diff.lo + diff.hi)
      X[i,3:P] <- w*Phi[lo,]+(1-w)*Phi[lo+1,]
    }
  }
  else{
    X <- matrix(0,n,cat.var)
    for(j in 1:cat.var)
      X[,j] <- (cat.var-1)/cat.var*(t==j) - 1/cat.var*(t!=j)
  }
  return(X)
}


### Produce 2-way interaction columns given the main effect columns for the pair ###
get.2way.X <- function(X1,X2,P,max.ord,var1,var2,cat.var1=0,cat.var2=0){

  n <- nrow(X1)
  P1 <- ifelse(cat.var1==0, P, cat.var1)
  P2 <- ifelse(cat.var2==0, P, cat.var2)

 ## Initialize X.12 (to be more columns than needed)
  maxcols <- P1*P2
  X.12 <- matrix(NA,nrow=n,ncol=maxcols)
  vars.basis <- matrix(NA,maxcols,4)
  col.now <- 1

 ## Take pairwise products of basis functions
 ## with the restriction that i+j only goes up to a maximum of max.ord
  for(i in 1:min(P1,max.ord)){
    for(j in index(1,min(P2,max.ord-i))){
      X.12[,col.now] <- X1[,i]*X2[,j]
      vars.basis[col.now,] <- c(var1,i,var2,j)
      col.now <- col.now+1
    }
  }
 ## Remove extra unused columns of X.12
  ind.keep <- which(!is.na(X.12[1,]))
  X.12 <- X.12[,ind.keep]
  X.12 <- as.matrix(X.12)
  if (dim(X.12)[2]==1){
  	X.12<-t(X.12)
  }
  vars.basis <- vars.basis[ind.keep,]
  return(list(X=as.matrix(X.12), vars.basis=vars.basis))
}



### Produce 3-way interaction columns given the main effect columns for the triplet ###
get.3way.X <- function(X1,X2,X3,P,max.ord,var1,var2,var3,cat.var1=0,cat.var2=0,cat.var3=0){

  n <- nrow(X1)
  P1 <- ifelse(cat.var1==0, P, cat.var1)
  P2 <- ifelse(cat.var2==0, P, cat.var2)
  P3 <- ifelse(cat.var3==0, P, cat.var3)
 ## Initialize X.123 (to be more columns than needed)
  maxcols <- P1*P2*P3
  X.123 <- matrix(NA,nrow=n,ncol=maxcols)
  vars.basis <- matrix(NA,maxcols,6)
  col.now <- 1

 ## Take pairwise products of basis functions
 ## with the restriction that i+j+k only goes up to a maximum of max.ord
  for(i in 1:min(P1,max.ord)){
    for(j in 1:min(P2,max.ord)){
      for(k in index(1,min(P3,max.ord-i-j))){
        X.123[,col.now] <- X1[,i]*X2[,j]*X3[,k]
        vars.basis[col.now,] <- c(var1,i,var2,j,var3,k)
        col.now <- col.now+1
      }
    }
  }
 ## Remove extra unused columns of X.123
  ind.keep <- which(!is.na(X.123[1,]))
  X.123 <- X.123[,ind.keep]
  X.123 <- as.matrix(X.123)
  if (dim(X.123)[2]==1){
  	X.123<-t(X.123)
  }
  vars.basis <- vars.basis[ind.keep,]
  return(list(X=as.matrix(X.123), vars.basis=vars.basis))
}




### Produce 4-way interaction columns given the main effect columns for the triplet ###
get.4way.X <- function(X1,X2,X3,X4,P,max.ord,var1,var2,var3,var4,cat.var1=0,cat.var2=0,cat.var3=0,cat.var4=0){

  n <- nrow(X1)
  P1 <- ifelse(cat.var1==0, P, cat.var1)
  P2 <- ifelse(cat.var2==0, P, cat.var2)
  P3 <- ifelse(cat.var3==0, P, cat.var3)
  P4 <- ifelse(cat.var4==0, P, cat.var4)
 ## Initialize X.1234 (to be more columns than needed)
  n.basis4 <- min(floor((max.ord+1)^4/24), P1*P2*P3*P4)
  X.1234 <- matrix(NA,nrow=n,ncol=n.basis4)
  vars.basis <- matrix(NA,n.basis4,8)
  col.now <- 1

 ## Take pairwise products of basis functions
 ## with the restriction that i+j+k only goes up to a maximum of max.ord
  for(i in 1:min(P1,max.ord)){
    for(j in 1:min(P2,max.ord)){
      for(k in 1:min(P3,max.ord)){
        for(l in index(1,min(P4,max.ord-i-j-k))){
          X.1234[,col.now] <- X1[,i]*X2[,j]*X3[,k]*X4[,l]
          vars.basis[col.now,] <- c(var1,i,var2,j,var3,k,var4,l)
          col.now <- col.now+1
        }
      }
    }
  }
 ## Remove extra unused columns of X.123
  ind.keep <- which(!is.na(X.1234[1,]))
  X.1234 <- X.1234[,ind.keep]
  X.1234 <- as.matrix(X.1234)
  if (dim(X.1234)[2]==1){
  	X.1234<-t(X.1234)
  }
  vars.basis <- vars.basis[ind.keep,]
  return(list(X=as.matrix(X.1234), vars.basis=vars.basis))
}



### Get one big design matrix ###
get.all.X <- function(t.mat, t.grid, Phi, P, max.ord, include2=NULL, include3=NULL, include4=matrix(0,0,4), cat.vars=rep(0,ncol(t.mat))){

 ## This function returns the entire design matrix of basis functions ...
 ## Inputs:
 ## t.mat - matrix of inputs (in our case columns are x,p,T,thetas)
 ## t.grid - the domain grid used for discritizing the basis functions on [0,1]
 ## Phi - the discretized eigenfunctions from K1 evaluated at each value of t.grid
 ## P - a 4-vector of how many eigeinfunctions to use for main effect, two-way, and
 ##     eventually 3-way interactions
 ## max.ord - the maximum total order (i.e., the sum of order of the baisis functions)
 ##           allowed in products of basis functions for interactions

  n <- nrow(t.mat)
  p <- ncol(t.mat)

 ## Make Main effect matrices
  X.list <- list()
  col.ind <- list()

 ## Store these length(t.grid) x P main effect matrices in a list (or 3D array) 

#t.mat..<<-t.mat
#Phi..<<-Phi
#t.grid..<<-t.grid
#cat.vars..<<-cat.vars

  for(j in 1:p){
    X.list[[j]] <- get.main.X(t.mat[,j],P[1],Phi,t.grid,mult=1,cat.var=cat.vars[j])
  }

 ## Make full design matrix with all basis vector columns included
  n.2way <- ifelse(is.null(include2),choose(p,2),nrow(include2))
  n.3way <- ifelse(is.null(include3),choose(p,3),nrow(include3))
  n.4way <- nrow(include4)
  n.basis2 <- min(floor((max.ord[2]+1)^2/2), P[2]^2)
  n.basis3 <- min(floor((max.ord[3]+1)^3/6), P[3]^3)
  n.basis4 <- min(floor((max.ord[4]+1)^4/24), P[4]^4)
  Q <- 1+p*P[1]+n.2way*n.basis2+n.3way*n.basis3+n.4way*n.basis4 
  X <- matrix(NA,n,Q)
  col2comp <- rep(NA,Q)
  vars.basis.ind <- matrix(NA,Q,8)
  vars.basis.ind[1,1:2] <- c(0,0)

 ## Intercept
  X[,1] <- 1
  col2comp[1] <- 0

 ## Main Effects
  ind.now <- 2
  for(j in 1:p){
    if(cat.vars[j]==0)
      cols.now <- ind.now:(P[1]+ind.now-1)
    else
      cols.now <- ind.now:(cat.vars[j]+ind.now-1)
   ## store the columns of X corresponding to each component
    P.now <- ifelse(cat.vars[j]==0, P[1], cat.vars[j])
    col.ind[[j]] <- cols.now
    vars.basis.ind[cols.now,1] <- j
    vars.basis.ind[cols.now,2] <- 1:P.now
   ## store the component corresponding to each column of X
    col2comp[cols.now] <- j
    X[,cols.now] <- X.list[[j]]
    ind.now <- ind.now + P.now
  }

 ## Two-Way Interactions
 ##  Loop through all p choose 2 input variable combos (if they are in include2)
  if(p>1){
    comp.now <- p+1 
    for(j in 1:(p-1)){
      for(k in (j+1):p){
        if(is.null(include2) || any(j==include2[,1] & k==include2[,2])){
          ans.2way <- get.2way.X(X.list[[j]],X.list[[k]],P=P[2],max.ord[2],j,k,cat.vars[j],cat.vars[k])
          X.now <- ans.2way$X
          cols.now <- ind.now:(ind.now+ncol(X.now)-1)
          X[,cols.now] <- X.now
          col.ind[[comp.now]] <- cols.now
          col2comp[cols.now] <- comp.now
          vars.basis.ind[cols.now,1:4] <- ans.2way$vars.basis
          comp.now <- comp.now+1
          ind.now <- ind.now+ncol(X.now)
        }
      }
    }
  }

 ## Three-Way Interactions  
 ##  Loop through all p choose 3 input variable combos (if they are in include3)
  if(p>2){
    for(j in 1:(p-2)){
      for(k in (j+1):(p-1)){
        for(l in (k+1):p){
          if(is.null(include3) || any(j==include3[,1] & k==include3[,2] & l==include3[,3])){
            ans.3way <- get.3way.X(X.list[[j]],X.list[[k]],X.list[[l]],P=P[3],max.ord[3],j,k,l,cat.vars[j],cat.vars[k],cat.vars[l])
            X.now <- ans.3way$X
            cols.now <- ind.now:(ind.now+ncol(X.now)-1)
            X[,cols.now] <- X.now
            col.ind[[comp.now]] <- cols.now
            col2comp[cols.now] <- comp.now
            vars.basis.ind[cols.now,1:6] <- ans.3way$vars.basis
            comp.now <- comp.now+1
            ind.now <- ind.now+ncol(X.now)
          }
        }
      }
    }
  }

 ## Four-Way Interactions  
 ##  Loop through all p choose 4 input variable combos (if they are in include4)
  if(p>3){
    for(j in 1:(p-3)){
      for(k in (j+1):(p-2)){
        for(l in (k+1):(p-1)){
          for(m in (l+1):p){
            if(any(j==include4[,1] & k==include4[,2] & l==include4[,3] & l==include4[,4])){
              ans.4way <- get.4way.X(X.list[[j]],X.list[[k]],X.list[[l]],X.list[[m]],P=P[4],max.ord[4],j,k,l,m,cat.vars[j],cat.vars[k],cat.vars[l],cat.vars[m])
              X.now <- ans.4way$X
              cols.now <- ind.now:(ind.now+ncol(X.now)-1)
              X[,cols.now] <- X.now
              col.ind[[comp.now]] <- cols.now
              col2comp[cols.now] <- comp.now
              vars.basis.ind[cols.now,1:8] <- ans.4way$vars.basis
              comp.now <- comp.now+1
              ind.now <- ind.now+ncol(X.now)
            }
          }
        }
      }
    }
  }
  ind.keep <- which(!is.na(X[1,]))
  X <- X[,ind.keep] 
  vars.basis.ind <- vars.basis.ind[ind.keep,]
  col2comp <- col2comp[!is.na(col2comp)]
  return(list(X=X,col.ind=col.ind,col2comp=col2comp, vars.basis.ind=vars.basis.ind))
}




###################################################
####### Conditional Deviates and Densities ########
###################################################



## observational error variance

rSigma.cond <- function(res.obs, P.Sigma, v.Sigma){

  v.post <- v.Sigma + nrow(res.obs)
  P.post <- P.Sigma + t(res.obs)%*%res.obs
  Sigma <- riwish(v.post, P.post)
  return(Sigma)
}



## simulator "solving" error variance

rUps.cond <- function(res.sim, P.Ups, v.Ups){

  v.post <- v.Ups + nrow(res.sim)
  P.post <- P.Ups + t(res.sim)%*%res.sim
  Ups <- riwish(v.post, P.post)
  return(Ups)
}



## discrepancy coefficient variance

rGamma.cond <- function(comp, col.ind, beta, P.Gamma, v.Gamma){

  if(comp==0)
    beta.comp <- matrix(beta[1,],nrow=1,ncol=ncol(beta))
  else
    beta.comp <- as.matrix(beta[col.ind[[comp]],])
  P.post <- P.Gamma + t(beta.comp)%*%beta.comp
  v.post <- v.Gamma + nrow(beta.comp)
  Gamma.comp <- riwish(v.post, P.post)
  return(Gamma.comp)
}



## emulator coefficient variance

rLambda.cond <- function(comp, col.ind, alpha, P.Lambda, v.Lambda){

  if(comp==0)
    alpha.comp <- matrix(alpha[1,],nrow=1,ncol=ncol(alpha))
  else
    alpha.comp <- as.matrix(alpha[col.ind[[comp]],])
  P.post <- P.Lambda + t(alpha.comp)%*%alpha.comp
  v.post <- v.Lambda + nrow(alpha.comp)
  Lambda.comp <- riwish(v.post, P.post)
  return(Lambda.comp)
}



## discrepancy intercept

rbeta0.cond <- function(res.obs, beta0.prev, Gamma.0, Sigma){

  I <- nrow(res.obs)
  pj <- 1
  M <- ncol(res.obs)
  y.rem <- as.numeric(res.obs + matrix(beta0.prev,I,M,byrow=TRUE))

  Sigma.inv <- ginv.gp(Sigma)$inv
  Gamma.inv <- ginv.gp(Gamma.0)$inv
  Xdis.0 <- as.matrix(rep(1,I))
  XtX.0 <- I

  sig.0.inv <- Sigma.inv %x% XtX.0 + Gamma.inv %x% diag(1,pj)
  ans <- ginv.gp(sig.0.inv)
  mu.0 <- ans$inv%*%(Sigma.inv %x% t(Xdis.0))%*%y.rem
  sig.0.sqrt <- ans$sqrt.inv
  log.det <- -ans$log.det

  #beta.0 <- matrix(rmvnorm(mu.0, S.sqrt=sig.0.sqrt),pj,M)
  beta.0 <- matrix(0,pj,M)
  #cat("\n beta.0",beta.0,"\n")
  dens.0 <- dmvnorm(as.numeric(beta.0), mu.0, S.inv=sig.0.inv, log.det=log.det)
  dens.0 <- 0*dens.0
  #cat("\n dens.0:", dens.0, "\n")
  #dens.0 <- matrix(0,nrow(dens.0),ncol(dens.0))
  dhat.0 <- rep(1,I)%*%beta.0
  res.obs <- matrix(y.rem,I,M) - dhat.0

  return(list(beta.0=as.numeric(beta.0), res.obs=res.obs, dens.0=dens.0))
}


dbeta0.cond <- function(res.obs, beta0.now, beta0.p, Gamma.0, Sigma){

  I <- nrow(res.obs)
  pj <- 1
  M <- ncol(res.obs)
  y.rem <- as.numeric(res.obs + matrix(beta0.now,I,M,byrow=TRUE))

  Sigma.inv <- ginv.gp(Sigma)$inv
  Gamma.inv <- ginv.gp(Gamma.0)$inv
  Xdis.0 <- as.matrix(rep(1,I))
  XtX.0 <- I

  sig.0.inv <- Sigma.inv %x% XtX.0 + Gamma.inv %x% diag(1,pj)
  ans <- ginv.gp(sig.0.inv)
  mu.0 <- ans$inv%*%(Sigma.inv %x% t(Xdis.0))%*%y.rem
  sig.0.sqrt <- ans$sqrt.inv
  log.det <- -ans$log.det

  #beta.0 <- matrix(beta0.p, pj, M)
  beta.0 <- matrix(0,pj,M)
  dens.0 <- dmvnorm(as.numeric(beta.0), mu.0, S.inv=sig.0.inv, log.det=log.det)
  ########
  dens.0 <- 0*dens.0
  #########
  dhat.0 <- rep(1,I)%*%beta.0
  res.obs <- matrix(y.rem,I,M) - dhat.0

  return(list(beta.0=as.numeric(beta.0), res.obs=res.obs, dens.0=dens.0))
}



## discrepancy component coefficients

rbeta.cond <- function(res.obs, Xdis.j, dhat.j, Gamma.j, Sigma){

  I <- nrow(res.obs)
  pj <- ncol(Xdis.j)
  M <- ncol(res.obs)
  y.rem <- as.numeric(res.obs + dhat.j)

  Sigma.inv <- ginv.gp(Sigma)$inv
  Gamma.inv <- ginv.gp(Gamma.j)$inv
  XtX.j <- t(Xdis.j)%*%Xdis.j

  sig.j.inv <- Sigma.inv %x% XtX.j + Gamma.inv %x% diag(1,pj)
  ans <- ginv.gp(sig.j.inv)
  mu.j <- ans$inv%*%(Sigma.inv %x% t(Xdis.j))%*%y.rem
  sig.j.sqrt <- ans$sqrt.inv
  log.det <- -ans$log.det

#res.obs..<<-res.obs
#x.j..<<-x.j
#y.rem..<<-y.rem
#beta.j.prev..<<-beta.j.prev
#Gamma.j..<<-Gamma.j
#Sigma..<<-Sigma
#mu.j..<<-mu.j
#sig.j.sqrt..<<-sig.j.sqrt


  #beta.j <- matrix(rmvnorm(mu.j, S.sqrt=sig.j.sqrt),pj,M)
  beta.j <- matrix(0,pj,M)
  dens.j <- dmvnorm(as.numeric(beta.j), mu.j, S.inv=sig.j.inv, log.det=log.det)
  #########
  dens.j <- 0*dens.j
  ########
  dhat.j <- Xdis.j%*%beta.j
  res.obs <- matrix(y.rem,I,M) - dhat.j

  return(list(beta.j=beta.j, res.obs=res.obs, dhat.j=dhat.j, dens.j=dens.j))
}


dbeta.cond <- function(res.obs, Xdis.j, betaj.p, dhat.j, Gamma.j, Sigma){

  I <- nrow(res.obs)
  pj <- ncol(Xdis.j)
  M <- ncol(res.obs)
  y.rem <- as.numeric(res.obs + dhat.j)

  Sigma.inv <- ginv.gp(Sigma)$inv
  Gamma.inv <- ginv.gp(Gamma.j)$inv
  XtX.j <- t(Xdis.j)%*%Xdis.j

  sig.j.inv <- Sigma.inv %x% XtX.j + Gamma.inv %x% diag(1,pj)
  ans <- ginv.gp(sig.j.inv)
  mu.j <- ans$inv%*%(Sigma.inv %x% t(Xdis.j))%*%y.rem
  sig.j.sqrt <- ans$sqrt.inv
  log.det <- -ans$log.det

  #beta.j <- matrix(betaj.p,pj,M)
  beta.j <- matrix(0,pj,M)
  dens.j <- dmvnorm(as.numeric(beta.j), mu.j, S.inv=sig.j.inv, log.det=log.det)
  #########
  dens.j <- 0*dens.j
  ########
  dhat.j <- Xdis.j%*%beta.j
  res.obs <- matrix(y.rem,I,M) - dhat.j

  return(list(beta.j=beta.j, res.obs=res.obs, dhat.j=dhat.j, dens.j=dens.j))
}



## emulator intercept

ralpha0.cond <- function(res.sim, res.obs, alpha0.prev, Lambda.0, Ups, Sigma){
  
  I <- nrow(res.obs)
  L <- nrow(res.sim)
  qj <- 1
  M <- ncol(res.obs)
  yobs.rem <- as.numeric(res.obs + matrix(alpha0.prev,I,M,byrow=TRUE))
  ysim.rem <- as.numeric(res.sim + matrix(alpha0.prev,L,M,byrow=TRUE))
  rem.star <- c(yobs.rem, ysim.rem)
  Sigma.inv <- ginv.gp(Sigma)$inv
  Ups.inv <- ginv.gp(Ups)$inv
  Lambda.inv <- ginv.gp(Lambda.0)$inv
  Xemd.0 <- as.matrix(rep(1,I))
  Xem.0 <- as.matrix(rep(1,L))
  XtX.obs <- I
  XtX.sim <- L

  blah <- Sigma.inv %x% XtX.obs + Ups.inv %x% XtX.sim + Lambda.inv %x% diag(1,qj)
  ans <- ginv.gp(blah)
  mu.0 <- ans$inv%*%((Sigma.inv %x% t(Xemd.0))%*%yobs.rem + (Ups.inv %x% t(Xem.0))%*%ysim.rem)
  sig.0.sqrt <- ans$sqrt.inv

  alpha.0 <- matrix(rmvnorm(mu.0, S.sqrt=sig.0.sqrt),qj,M)
  yhat.0 <- rep(1,L)%*%alpha.0
  res.sim <- matrix(ysim.rem,L,M) - yhat.0
  yhatd.0 <- rep(1,I)%*%alpha.0
  res.obs <- matrix(yobs.rem,I,M) - yhatd.0

  return(list(alpha.0=as.numeric(alpha.0), res.sim=res.sim, res.obs=res.obs))
}



## emulator component coefficients

ralpha.cond <- function(res.sim, res.obs, Xem.j, Xemd.j, yhat.j, yhatd.j, Lambda.j, Ups, Sigma){
  
  I <- nrow(res.obs)
  L <- nrow(res.sim)
  qj <- ncol(Xem.j)
  M <- ncol(res.obs)
  yobs.rem <- as.numeric(res.obs + yhatd.j)
  ysim.rem <- as.numeric(res.sim + yhat.j)

  Sigma.inv <- ginv.gp(Sigma)$inv
  Ups.inv <- ginv.gp(Ups)$inv
  Lambda.inv <- ginv.gp(Lambda.j)$inv
  XtX.obs <- t(Xemd.j)%*%Xemd.j
  XtX.sim <- t(Xem.j)%*%Xem.j

#res.sim..<<-res.sim
#yhat.j..<<-yhat.j
#Xem.j..<<-Xem.j
#Xemd.j..<<-Xemd.j
#ysim.rem..<<-ysim.rem
#yobs.rem..<<-yobs.rem
#Sigma..<<-Sigma
#Ups..<<-Ups
#Lambda.j..<<-Lambda.j

  blah <- Sigma.inv %x% XtX.obs + Ups.inv %x% XtX.sim + Lambda.inv %x% diag(1,qj)
  ans <- ginv.gp(blah)
  mu.j <- ans$inv%*%((Sigma.inv %x% t(Xemd.j))%*%yobs.rem + (Ups.inv %x% t(Xem.j))%*%ysim.rem)
  sig.j.sqrt <- ans$sqrt.inv

  alpha.j <- matrix(rmvnorm(mu.j, S.sqrt=sig.j.sqrt),qj,M)
  yhat.j <- Xem.j%*%alpha.j
  res.sim <- matrix(ysim.rem,L,M) - yhat.j
  yhatd.j <- Xemd.j%*%alpha.j
  res.obs <- matrix(yobs.rem,I,M) - yhatd.j

  return(list(alpha.j=alpha.j, res.sim=res.sim, yhat.j=yhat.j, res.obs=res.obs, yhatd.j=yhatd.j))
}



## Likelihood for theta given data and rest

get.like.theta <- function(res.obs, Sigma){

  ans <- ginv.gp(Sigma)
  S.inv <- ans$inv
  log.det <- ans$log.det
  I <- nrow(res.obs)
  like <- 0
  for(i in 1:I)
    like <- like+dmvnorm(res.obs[i,], 0, S.inv=S.inv, log.det=log.det)
  return(like)
}




## Theta proposal function

r.theta.proposal <- function(prev, scale, cat.theta){
  if(cat.theta==0){
    prev.x <- logit(prev)
    prop.x <- rnorm(1, prev.x, (abs(prev.x)+1)*scale)
    return(inv.logit(prop.x))
  }
  else{
#    vec <- (1:cat.theta)[-prev]
    vec <- (1:cat.theta)
    return(sample(vec, 1))
  }
}


dlogitnorm <- function(x, mu, sigma){
  ans <- dnorm(logit(x),mu,sigma, log=TRUE) - log(x) - log(1-x)
  return(ans)
}

d.theta.proposal <- function(prop, prev, scale, cat.theta){
  if(cat.theta==0){
    prev.x <- logit(prev)
    ans <- dlogitnorm(prop, prev.x, (abs(prev.x)+1)*scale)
  }
  else
#    ans <- -log(cat.theta-1)
    ans <- -log(cat.theta)
  return(ans)
}




## Update missing residuals
update.resid <- function(res.obs.now, Sigma.now, ind.na.obs){

#res.obs.now..<<-res.obs.now
#Sigma.now..<<-Sigma.now
#ind.na.obs..<<-ind.na.obs

  I <- nrow(res.obs.now)
  for(i in 1:I){
    if(any(ind.na.obs[i,])){
      ind.na.i <- ind.na.obs[i,]
      S11 <- as.matrix(Sigma.now[ind.na.i, ind.na.i])
      S22 <- as.matrix(Sigma.now[!ind.na.i, !ind.na.i])
      S12 <- matrix(Sigma.now[ind.na.i, !ind.na.i], nrow=sum(ind.na.i), ncol=sum(!ind.na.i))
      S22.inv <- ginv.gp(S22)$inv
      mu.miss <- S12%*%S22.inv%*%res.obs.now[i,!ind.na.i]
      S.miss <- S11-S12%*%S22.inv%*%t(S12)
      eps.miss <- rmvnorm(mu=mu.miss, Sigma=S.miss)
      res.obs.now[i,ind.na.i] <- eps.miss
    }
  }
  return(res.obs.now)
}



## Update alphas
update.alphas <- function(res.sim.now, res.obs.now, alpha.now, X.em, X.emd, col.ind.em, ncomp.em, yhat.mat, yhatd.mat, Lambda.now, Ups.now, Sigma.now){

  alpha0.ans <- ralpha0.cond(res.sim.now, res.obs.now, alpha.now[1,], as.matrix(Lambda.now[1,,]), Ups.now, Sigma.now)
  alpha.now[1,] <- alpha0.ans$alpha.0
  res.obs.now <- alpha0.ans$res.obs
  res.sim.now <- alpha0.ans$res.sim

  for(j in 1:ncomp.em){
    cols.j <- col.ind.em[[j]]
    Xem.j <- X.em[,cols.j]
    Xemd.j <- X.emd[,cols.j]

    alpha.ans <- ralpha.cond(res.sim.now, res.obs.now, Xem.j, Xemd.j, as.matrix(yhat.mat[,j,]), as.matrix(yhatd.mat[,j,]), as.matrix(Lambda.now[j+1,,]), Ups.now, Sigma.now)
    alpha.now[cols.j,] <- alpha.ans$alpha.j
    yhat.mat[,j,] <- alpha.ans$yhat.j
    yhatd.mat[,j,] <- alpha.ans$yhatd.j
    res.obs.now <- alpha.ans$res.obs
    res.sim.now <- alpha.ans$res.sim
  }
  return(list(alpha.now=alpha.now, res.obs.now=res.obs.now, res.sim.now=res.sim.now, yhat.mat=yhat.mat, yhatd.mat=yhatd.mat))
}



## Update betas
update.betas <- function(res.obs.now, beta.now, X.dis, col.ind.dis, ncomp.dis, dhat.mat, Gamma.now, Sigma.now){

  beta0.ans <- rbeta0.cond(res.obs.now, beta.now[1,], as.matrix(Gamma.now[1,,]), Sigma.now)
  beta.now[1,] <- beta0.ans$beta.0
  res.obs.now <- beta0.ans$res.obs
  dens <- beta0.ans$dens.0

  for(j in 1:ncomp.dis){
    cols.j <- col.ind.dis[[j]]
    Xdis.j <- X.dis[,cols.j]
    beta.ans <- rbeta.cond(res.obs.now, Xdis.j, as.matrix(dhat.mat[,j,]), as.matrix(Gamma.now[j+1,,]), Sigma.now)
    beta.now[cols.j,] <- beta.ans$beta.j
    dhat.mat[,j,] <- beta.ans$dhat.j
    res.obs.now <- beta.ans$res.obs
    dens <- dens + beta.ans$dens.j
  }
  return(list(beta.now=beta.now, res.obs.now=res.obs.now, dhat.mat=dhat.mat, dens=dens))
}



## proposal density for betas
dprop.betas <- function(res.obs.now, beta.now, beta.p, X.dis, col.ind.dis, ncomp.dis, dhat.mat.now, Gamma.now, Sigma.now){

  ans.dens <- dbeta0.cond(res.obs.now, beta.now[1,], beta.p[1,], as.matrix(Gamma.now[1,,]), Sigma.now)
  beta.now[1,] <- ans.dens$beta.0
  res.obs.now <- ans.dens$res.obs
  dens <- ans.dens$dens

  for(j in 1:ncomp.dis){
    cols.j <- col.ind.dis[[j]]
    Xdis.j <- X.dis[,cols.j]
    betaj.p <- beta.p[cols.j,]
    dens.ans <- dbeta.cond(res.obs.now, Xdis.j, betaj.p, as.matrix(dhat.mat.now[,j,]), as.matrix(Gamma.now[j+1,,]), Sigma.now)
    beta.now[cols.j,] <- dens.ans$beta.j
    dhat.mat.now[,j,] <- dens.ans$dhat.j
    res.obs.now <- dens.ans$res.obs
    dens <- dens + dens.ans$dens.j
  }
  return(dens)
}



pi.beta <- function(beta, Gamma, col.ind.dis){
  ans <- dmvnorm(beta[1,], 0, Gamma[1,,])
  for(j in 1:length(col.ind.dis)){
    col.ind.j <- col.ind.dis[[j]]
    for(l in col.ind.j){
      ans <- ans + dmvnorm(beta[l,], 0, Gamma[j+1,,])
    }
  }
  return(ans)
}



####################################################################
############# Calibration MCMC Estimation Function #################
####################################################################

Calibrate <- function(sim.dat, exp.dat, P, Q, M, N.mcmc, cat.locations=NULL, X.lim=NULL, T.lim=NULL, transform=NULL, pi.theta=NULL, mu.Sigma=NULL, v.Sigma=NULL, mu.Ups=NULL, v.Ups=NULL, theta.prop.sd=NULL, sim.out1=NULL, sim.out2=NULL, sim.out3=NULL, sim.out4=NULL, sim.out5=NULL, sim.out6=NULL, sim.out7=NULL, sim.out8=NULL, sim.out9=NULL, sim.out10=NULL, exp.out1=NULL, exp.out2=NULL, exp.out3=NULL, exp.out4=NULL, exp.out5=NULL, exp.out6=NULL, exp.out7=NULL, exp.out8=NULL, exp.out9=NULL, exp.out10=NULL, y.names=NULL, theta.now=NULL, nplot=100, Sigma.eq.Ups=FALSE){

  x.locs <- 1:P
  y.locs <- (P+1):(P+M)
  t.locs <- (P+1):(P+Q)
  ys.locs <- (P+Q+1):(P+Q+M)
  X <- as.matrix(exp.dat[,x.locs])
  Xs <- as.matrix(sim.dat[,x.locs])
  Ts <- sim.dat[,t.locs]
  I <- nrow(X)
  L <- nrow(Xs)

  trans.fun <- invtrans.fun <- list()
  pow <- rep(1,M)
  if(is.null(transform))
    transform <- rep(1,M)
  for(m in 1:M){
    if(transform[m]=="log"){
      trans.fun[[m]] <- function(y,a) log(y+1E-4)
      invtrans.fun[[m]] <- function(z,a) exp(z)-1E-4
    }
    else if(transform[m]=="logit"){
      trans.fun[[m]] <- function(y,a) logit(y+1E-4)
      invtrans.fun[[m]] <- function(z,a) inv.logit(z)-1E-4
    }
    else if(transform[m]==1){
      trans.fun[[m]] <- function(y,a) y
      invtrans.fun[[m]] <- function(z,a) z
    }  
    else{
      pow[m] <- as.numeric(transform[m])
      trans.fun[[m]] <- function(x, a){
        c <- .001
        knot <- .001
        b <- (knot+c)^a
        const <- b*(1-log(b))
        ans <- (x+c)^a*(x>knot) + (a*b*log(x+c)+const)*(x<=knot)
        return(ans)
      }
      invtrans.fun[[m]] <- function(y, a){
        c <- .001
        knot <- .001
        b <- (knot+c)^a
        const <- b*(1-log(b))
        x1 <- y^(1/a)-c
        x2 <- exp((y-const)/(a*b))-c
        ans <- ifelse(x2<knot,x2,x1)
        return(ans)
      }
    }
  }


  if(is.null(sim.out1)){
    y <- as.matrix(exp.dat[,y.locs])
    ys <- as.matrix(sim.dat[,ys.locs])
    y.names <- names(sim.dat)[ys.locs]
  }
  else{ ## convert time series to mean and SE (use SE to create prior for error var)
    if(M>10)
      stop("Currently only supports M<=10 outputs when using entire time series")
    y <- se2.y <- matrix(NA,I,M)
    ys <- se2.ys <- matrix(NA,L,M)
    for(m in 1:M){
      sim.m <- get(paste("sim.out",m,sep=""))
      exp.m <- get(paste("exp.out",m,sep=""))
      for(i in 1:I){
        y[i,m] <- mean(as.numeric(exp.m[i,]))
        se2.y[i,m] <- var_of_mean(as.numeric(exp.m[i,]), y[i,m])
      }
      for(l in 1:L){
        ys[l,m] <- mean(as.numeric(sim.m[l,]))
        se2.ys[l,m] <- var_of_mean(as.numeric(sim.m[l,]), ys[l,m])
      }
    }
    mu.sigma <- apply(se2.y,2,mean)
    var.sigma <- apply(se2.y,2,var)/I
    mu.ups <- apply(se2.ys,2,mean)
    var.ups <- apply(se2.ys,2,var)/L

    mu.Sigma <- diag(mu.sigma,M,M)
    mu.Ups <- diag(mu.ups,M,M)

    v.vec.sigma <- 2*mu.sigma^2/var.sigma+M+3
    v.Sigma <- min(v.vec.sigma)
    v.vec.ups <- 2*mu.ups^2/var.ups+M+3
    v.Ups <- min(v.vec.ups)

    if(is.null(y.names))
      y.names <- paste("Y_",1:M,sep="")

  }

#P.Sigma <- mu.Sigma*(v.Sigma-M-1)
#blah <- array(0, c(10000,1,1))
#for(i in 1:10000) blah[i,,] <- riwish(v.Sigma, P.Sigma)
#var(blah)


  x.names <- names(sim.dat)[x.locs]
  theta.names <- names(sim.dat)[t.locs]

  cat.theta <- rep(0,Q)
  cat.names <- list()
  if(is.null(cat.locations))
    cat.locations <- rep(0,0)
  if(length(cat.locations)>0){
    for(k in cat.locations){
      cat.theta[k] <- length(unique(Ts[,k]))
      cat.names[[k]] <- unique(Ts[,k])
    }
  }


  if(is.null(T.lim)){
    T.lim <- matrix(0,nrow=K,ncol=2)
    for(k in 1:K){
      if(cat.theta[k]==0){
        min.k <- min(Ts[,j])
        max.k <- max(Ts[,j])
        R.k <- max.k-min.k
        T.lim[k,] <- c(min.k-0.1*R.k, max.k+0.1*R.k)
      }
      else{
        T.lim[k,] <- c(0,1)
        labels <- unique(Ts[,k])
        Ts[,k] <- match(Ts[,k],labels)
      }
    }
  }
  T.lim[is.na(T.lim[,1]),1] <- 0
  T.lim[is.na(T.lim[,2]),2] <- 1


  if(is.null(pi.theta)){
    pi.theta <- function(theta){
      indc <- (cat.theta==0)
      ans.cont <- dunif(theta[indc], T.lim[indc,1], T.lim[indc,2], log=TRUE)
      ans.cat <- sum(-log(cat.theta[!indc]))
      return(ans.cont+ans.cat)
    }
  }
  
  if(Sigma.eq.Ups){
    mu.Sigma <- mu.Ups <- (mu.Sigma+mu.Ups)/2
    v.Sigma <- v.Ups <- min(v.Sigma, v.Ups)
  }

y..<<-y
X..<<-X
ys..<<-ys
Xs..<<-Xs
Ts..<<-Ts
cat.theta..<<-cat.theta
pi.theta..<<-pi.theta
X.lim..<<-X.lim
T.lim..<<-T.lim
transform..<<-trans.fun
inv.transform..<<-invtrans.fun
N.mcmc..<<-N.mcmc
mu.Sigma..<<-mu.Sigma
v.Sigma..<<-v.Sigma
mu.Ups..<<-mu.Ups
v.Ups..<<-v.Ups
theta.prop.sd..<<-theta.prop.sd
theta.names..<<-theta.names
cat.names..<<-cat.names
x.names..<<-x.names
y.names..<<-y.names

  ans.cal <- BSSANOVA.calibration(y, X, ys, Xs, Ts, cat.theta=cat.theta, pi.theta=pi.theta, X.lim=X.lim, T.lim=T.lim, transform=trans.fun, inv.transform=invtrans.fun, N.mcmc=N.mcmc, nplot=nplot, nback=5000, mu.Sigma=mu.Sigma, v.Sigma=v.Sigma, mu.Ups=mu.Ups, v.Ups=v.Ups, theta.prop.sd=theta.prop.sd, theta.names=theta.names, cat.names=cat.names, x.names=x.names, y.names=y.names, pow=pow, theta.now=theta.now)
  return(ans.cal)
}

  










BSSANOVA.calibration <- function(y, X, ys, Xs, Ts, cat.theta=rep(0,ncol(Ts)), pi.theta=NULL, X.lim=NULL, T.lim=NULL, transform=NULL, inv.transform=NULL, int.order=NULL, P=c(25,8,6,4), max.ord=P+1, P.Sigma=NULL, v.Sigma=NULL, P.Ups=NULL, v.Ups=NULL, P.Gamma=NULL, v.Gamma=NULL, P.Lambda=NULL, v.Lambda=NULL, N.grid=500, N.mcmc=6000, nplot=100, nback=5000, theta.now=NULL, alpha.now=NULL, beta.now=NULL, Lambda.now=NULL, Gamma.now=NULL, Ups.now=NULL, Sigma.now=NULL, theta.prop.sd=1, include2="default", include3=NULL, include4=NULL, mu.Sigma=NULL, mu.Ups=NULL, mu.Gamma=NULL, mu.Lambda=NULL, N.init=0, theta.names=NULL, cat.names=NULL, x.names=NULL, y.names=NULL, pow=rep(1,ncol(y))){

##INPUTS
##  y - matrix of experimental obs (rows are observed outputs at the rows of X)
##  X - matrix of experimental input values (e.g, angle and fluid velocity)
##  ys - matrix of model outputs (rows are outputs evaluated at the rows of Xs)
##  Xs - matrix of model input values (e.g., angle and fluid velocity)
##  Ts - matrix of model parameters values
##  X.lim - 2 column matrix each row has limits for the inputs
##  T.lim - 2 column matrix each row has limits for the theta params
##  tpow - a vector of length ncol(y) of transformation powers to use on y
##  cat.theta - vector of number of categories for thetas (0 indicates continuous theta)
##  pi.theta - function evaluating prior density for theta
##  P.Lambda, v.Lambda - Lambda (emulator magnitude) invgamma prior params
##  P.Gamma, v.Gamma - Gamma (discrepancy magnitude) invgamma prior params
##  P.Ups, v.Ups - Ups (simulator "jitter" error variance) inv-Wishart prior params, P is a positive definite matrix, n is an integer.
##  P.Sigma, v.Sigma - Sigma (observational error variance) inv-Wishart prior params
##  N.mcmc - Total number of MCMC iterations


##OUTPUTS
##  X.dis - matrix of basis functions for discrepancy
##  X.em <- matrix of basis functions for emulator
##  beta - posterior sample array of discrepancy coefs dim [N.mcmc, ncol(X), ncol(y)]
##  alpha - posterior sample array of emulator coefs dim [N.mcmc, ncol(X), ncol(y)]
##  Gamma - posterior sample array of dim [N.mcmc, ncomp.dis, ncol(y)]
##  Lambda - posterior sample array of dim [N.mcmc, ncomp.em, ncol(y)]
##  Sigma - posterior sample array of dim [N.mcmc, ncol(y), ncol(y)]
##  Ups - posterior sample array of dim [N.mcmc, ncol(y), ncol(y)]
##  theta - posterior sample matrix of dim [N.mcmc, q]

  y.o <- y
  ys.o <- ys
  X.o <- X
  Xs.o <- Xs
  Ts.o <- Ts

  y <- as.matrix(y)
  ys <- as.matrix(ys)
  I <- nrow(y)
  L <- nrow(ys)
  J <- ncol(X)
  K <- ncol(Ts)
  M <- ncol(y)

 ## Transform output (to stabilize variance, etc...)
  if(is.null(transform)){
    transform <- list()
    for(m in 1:M){
      transform[[m]] <- function(x, a){ return(x)}
      inv.transform[[m]] <- function(x, a){ return(x)}
    }
  }

  for(m in 1:M){
    y[,m] <- as.numeric(y[,m])
    ys[,m] <- as.numeric(ys[,m])
    y[,m] <- transform[[m]](y[,m], pow[m])
    ys[,m] <- transform[[m]](ys[,m], pow[m])
  }

 ## Transform inputs to [0,1]
  if(is.null(X.lim)){
    X.lim <- matrix(0,nrow=J,ncol=2)
    for(j in 1:J){
      min.j <- min(c(X[,j], Xs[,j]))
      max.j <- max(c(X[,j], Xs[,j]))
      R.j <- max.j-min.j
      X.lim[j,] <- c(min.j-0.1*R.j, max.j+0.1*R.j)
    }
  }
  for(j in 1:J){
    X[,j] <- (X[,j]-X.lim[j,1])/(X.lim[j,2]-X.lim[j,1])
    Xs[,j] <- (Xs[,j]-X.lim[j,1])/(X.lim[j,2]-X.lim[j,1])
  }

  if(is.null(y.names))
    y.names <- paste("Y_",1:M, sep="")
  if(is.null(x.names))
    x.names <- paste("X_",1:J, sep="")
  if(is.null(theta.names))
    theta.names <- paste("theta_",1:K, sep="")
  if(is.null(cat.names)){
    cat.names <- list()
    cat.flag <- 1
  }

 ## Transform continuous model params to [0,1]
  if(is.null(T.lim)){
    T.lim <- matrix(0,nrow=K,ncol=2)
    for(k in 1:K){
      if(cat.theta[k]==0){
        min.k <- min(Ts[,k])
        max.k <- max(Ts[,k])
        R.k <- max.k-min.k
        T.lim[k,] <- c(min.k-0.1*R.k, max.k+0.1*R.k)
      }
      else{
        T.lim[k,] <- c(0,1)
        labels <- unique(Ts[,k])
        if(cat.flag)
          cat.names[[k]] <- labels
        Ts[,k] <- match(Ts[,k],labels)
      }
    }
  }
  T.lim[is.na(T.lim[,1]),1] <- 0
  T.lim[is.na(T.lim[,2]),2] <- 1

  for(j in 1:K){
    if(cat.theta[j]==0)
      Ts[,j] <- (Ts[,j]-T.lim[j,1])/(T.lim[j,2]-T.lim[j,1])
    else
      Ts[,j] <- as.numeric(Ts[,j])
  }

 ## Convert prior pi.theta to act on [0,1]
  pi.theta.t <- function(theta){
    theta.o <- theta*(T.lim[,2]-T.lim[,1])+T.lim[,1]
    return(pi.theta(theta.o))
  }

#pi.theta.t..<<-pi.theta.t

 ## Set order of interactions
  if(is.null(int.order) && !is.null(include2) && include2=="default"){
    include2 <- NULL
    include4 <- matrix(0,0,4)
    if(J==1)
      include3 <- matrix(0,0,3)
    else if(J==2)
      include3 <- cbind(1,2,3:(J+K))
    else if(J==3)
      include3 <- rbind(cbind(1,2,4:(J+K)), cbind(1,3,4:(J+K)), cbind(2,3,4:(J+K)))
    else if(J==4)
      include3 <- rbind(cbind(1,2,5:(J+K)), cbind(1,3,5:(J+K)), cbind(1,4,5:(J+K)), cbind(2,3,5:(J+K)), cbind(2,4,5:(J+K)), cbind(3,4,5:(J+K)))
    else
      include <- matrix(0,0,3)
  }
  if(!is.null(int.order) && !is.null(include2) && include2=="default"){
    include2 <- include3 <- include4 <- NULL
    if(int.order<2)  include2 <- matrix(0,0,2)
    if(int.order<3)  include3 <- matrix(0,0,3)
    if(int.order<4)  include4 <- matrix(0,0,4)
  }

  SSTo <- SSTo.s <- rep(0,M)
  for(m in 1:M){
    SSTo[m] <- sum((y[,m]-mean(y[,m], na.rm=TRUE))^2, na.rm=TRUE)
    SSTo.s[m] <- sum((ys[,m]-mean(ys[,m], na.rm=TRUE))^2, na.rm=TRUE)
  }

 ## Initialize basis and design mat for discrepancy
  t.grid <- seq(0,1,length=N.grid)
  Phi <- get.bss.basis(t.grid)
  ans <- get.all.X(X, t.grid, Phi, P=P, max.ord=max.ord, include2, include3, include4, cat.vars=rep(0,ncol(X)))
  X.dis <- ans$X
  col.ind.dis <- ans$col.ind
  col2comp.dis <- ans$col2comp
  vars.basis.ind.dis <- ans$vars.basis.ind
  ncomp.dis <- length(col.ind.dis)

 ## Initialize basis and design mat for emulator at model runs
  ans <- get.all.X(cbind(Xs,Ts), t.grid, Phi, P=P, max.ord=max.ord, include2, include3, include4, cat.vars=c(rep(0,ncol(Xs)), cat.theta))
  X.em <- ans$X
  col.ind.em <- ans$col.ind
  col2comp.em <- ans$col2comp
  vars.basis.ind.em <- ans$vars.basis.ind
  ncomp.em <- length(col.ind.em)

#X.em..<<-X.em
#col.ind.em..<<-col.ind.em
#col.ind.dis..<<-col.ind.dis
#col2comp.em..<<-col2comp.em
#vars.basis.ind.em..<<-vars.basis.ind.em
#X..<<-X
#col2comp..<<-col2comp
#col.ind..<<-col.ind 
#vars.basis.ind.dis..<<-vars.basis.ind.dis
#vars.basis.ind.em..<<-vars.basis.ind.em
#Phi..<<-Phi

 ## Initialize posterior objects
  alpha <- array(0, c(N.mcmc, ncol(X.em), M))
  alpha[,1,] <- matrix(colMeans(ys),nrow=N.mcmc,ncol=M)
  beta <- array(0, c(N.mcmc, ncol(X.dis), M))
  Sigma <- array(1, c(N.mcmc, M, M))
  Ups <- array(1, c(N.mcmc, M, M))
  Gamma <- array(1,c(N.mcmc,ncomp.dis+1,M,M))
  Lambda <- array(1,c(N.mcmc,ncomp.em+1,M,M))
  theta <- matrix(.5, N.mcmc, K)
  theta[,cat.theta!=0] <- 1
  accept.count <- matrix(0, N.mcmc, K)
  res.obs <- array(0,c(N.mcmc, I, M))
  res.obs.sim <- array(0,c(N.mcmc, I, M))
  res.sim <- array(0,c(N.mcmc, L, M))
  y.complete <- array(0,c(N.mcmc, I, M))

  if(length(theta.prop.sd)==1)
    theta.prop.sd <- rep(theta.prop.sd,K)


  if(is.null(v.Sigma))
    v.Sigma <- I
  if(is.null(v.Ups))
    v.Ups <- I
  if(is.null(v.Gamma))
    v.Gamma <- M+3
  if(is.null(v.Lambda))
    v.Lambda <- M+3

  var.ys <- apply(ys, 2, var)
#####
cat("\nkkkkkkk\n ys is ",ys)
cat("\n\n var.ys is ", var.ys)
cat("\n\n mean.ys is ", apply(ys,2,mean))
cat("\n\n size ys is ", nrow(ys))
####
  if(is.null(mu.Lambda)){
    if(is.null(P.Lambda))
      mu.Lambda <- diag(var.ys/4,M,M)
    else
      mu.Lambda <- P.Lambda/(v.Lambda-M-1)
  }
  if(is.null(mu.Sigma)){
    if(is.null(P.Sigma))
      mu.Sigma <- mu.Lambda/10
    else
      mu.Sigma <- P.Sigma/(v.Sigma-M-1)
  }
  if(is.null(mu.Ups)){
    if(is.null(P.Ups))
      mu.Ups <- mu.Sigma
    else
      mu.Ups <- P.Ups/(v.Ups-M-1)
  }
  if(is.null(mu.Gamma)){
    if(is.null(P.Gamma))
      mu.Gamma <- mu.Lambda
    else
      mu.Gamma <- P.Gamma/(v.Gamma-M-1)
  }

  if(is.null(P.Sigma))
    P.Sigma <- mu.Sigma*(v.Sigma-M-1)
  if(is.null(P.Ups))
    P.Ups <- mu.Ups*(v.Ups-M-1)
  if(is.null(P.Gamma))  
    P.Gamma <-  mu.Gamma*(v.Gamma-M-1)
  if(is.null(P.Lambda))
    P.Lambda <- mu.Lambda*(v.Lambda-M-1)


 ## Initialize initial values if Null
  if(is.null(alpha.now))
    alpha.now <- as.matrix(alpha[1,,])
  if(is.null(beta.now))
    beta.now <- as.matrix(beta[1,,])
#    beta.now <- matrix(0,nrow(beta.now),ncol(beta.now))
  if(is.null(Sigma.now))
    Sigma.now <- mu.Sigma
  if(is.null(Ups.now))
    Ups.now <- mu.Ups
# set Gamma to be 0    
#  if(is.null(Gamma.now)){
    Gamma.now <- array(0,c(ncomp.dis+1,M,M))
#    for(j in 1:dim(Gamma.now)[1])
#      Gamma.now[j,,] <- mu.Gamma
#  }

  if(is.null(Lambda.now)){
    Lambda.now <- array(0,c(ncomp.em+1,M,M))
    for(j in 1:dim(Lambda.now)[1])
      Lambda.now[j,,] <- mu.Lambda
  }


  if(is.null(theta.now))
    theta.now <- theta[1,]
  else
    theta.now <- (theta.now-T.lim[,1])/apply(T.lim,1,diff)


 ## Initial values for Sigma and Ups prior params are tight for first N.init to speed up 
  
  v.Sigma.perm <- v.Sigma
  P.Sigma.perm <- P.Sigma
  v.Ups.perm <- v.Ups
  P.Ups.perm <- P.Ups
  mu.Sigma <- P.Sigma/(v.Sigma-M-1)
  v.Sigma <- 100
  P.Sigma <- mu.Sigma*(v.Sigma-M-1)
  mu.Ups <- P.Ups/(v.Ups-M-1)
  v.Ups <- 100
  P.Ups <- mu.Ups*(v.Ups-M-1)


 ## Initialize basis and design mat for emulator at experimental obs
  ans <- get.all.X(cbind(X,matrix(theta.now,nrow=nrow(X),ncol=K,byrow=TRUE)), t.grid, Phi, P=P, max.ord=max.ord, include2, include3, include4, cat.vars=c(rep(0,ncol(X)), cat.theta))
  X.emd <- ans$X

 ## Initialize residuals
  res.obs.now <- as.matrix(y-as.numeric(X.emd%*%alpha.now+X.dis%*%beta.now))
  ind.na.obs <- is.na(y)
  res.obs.now[ind.na.obs] <- 0
  res.sim.now <- as.matrix(ys-as.numeric(X.em%*%alpha.now))

 ## create arrays of component yhats 
#  XtX.list <- list()
  yhat.mat <- array(0,c(L,ncomp.em,M))
  yhatd.mat <- array(0,c(I,ncomp.em,M))
  dhat.mat <- array(0,c(I,ncomp.dis,M))
  for(j in 1:ncomp.em){
    cols.em <- col.ind.em[[j]]
    yhat.mat[,j,] <- X.em[,cols.em]%*%alpha.now[cols.em,]
    yhatd.mat[,j,] <- X.emd[,cols.em]%*%alpha.now[cols.em,]
  }
  for(j in 1:ncomp.dis){
    cols.dis <- col.ind.dis[[j]]
    dhat.mat[,j,] <- X.dis[,cols.dis]%*%beta.now[cols.dis,]
  }


 ## Begin MCMC ##
  cat("\n")
  for(it in 1:N.mcmc){

cat("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bIteration", it, "out of", N.mcmc)

    if(it > N.init){
      v.Sigma <- v.Sigma.perm
      P.Sigma <- P.Sigma.perm
      v.Ups <- v.Ups.perm
      P.Ups <- P.Ups.perm
    }


   ## Update missing residuals
    res.obs.now <- update.resid(res.obs.now, Sigma.now, ind.na.obs)


   ## Update alphas
    ans.alpha <- update.alphas(res.sim.now, res.obs.now, alpha.now, X.em, X.emd, col.ind.em, ncomp.em, yhat.mat, 		  yhatd.mat, Lambda.now, Ups.now, Sigma.now)
    alpha.now <- ans.alpha$alpha.now
    res.obs.now <- ans.alpha$res.obs.now
    res.sim.now <- ans.alpha$res.sim.now
    yhat.mat <- ans.alpha$yhat.mat
    yhatd.mat <- ans.alpha$yhatd.mat


   ## Update betas
# here we do not update the discrepancy mean, keep it as zero
#ans.beta <- update.betas(res.obs.now, beta.now, X.dis, col.ind.dis, ncomp.dis, dhat.mat, Gamma.now, Sigma.now)

#beta.now <- ans.beta$beta.now
#res.obs.now <- ans.beta$res.obs.now
#dhat.mat <- ans.beta$dhat.mat


   ## Update Lambdas
    for(j in 0:ncomp.em)
      Lambda.now[j+1,,] <- rLambda.cond(j, col.ind.em, alpha.now, P.Lambda, v.Lambda)


   ## Update Gammas

#here we set discrepancy into zero, the Gamma is set to zero!     
#for(j in 0:ncomp.dis)
#Gamma.now[j+1,,] <- rGamma.cond(j, col.ind.dis, beta.now, P.Gamma, v.Gamma)
#cat("\nGamma.now size ", ncol(Gamma.now));  
#cat("\ncomp.dis size ", ncomp.dis);  
#cat("Gamma.now is ", Gamma.now, "\n");

   ## Update Ups
#here we consider that there is no further discrepancy on simulation
#so we do not update  Ups, Ups keeps zero        
#   Ups.now <- rUps.cond(res.sim.now, P.Ups, v.Ups)


   ## Update Sigma
#here we fix Sigma as the default value        
#    Sigma.now <- rSigma.cond(res.obs.now, P.Sigma, v.Sigma)

   ## Update Thetas
    ytru.now <- as.numeric(X.emd%*%alpha.now) + as.numeric(X.dis%*%beta.now)
    yc.now <- ytru.now + res.obs.now
    for(k in 1:K){
      theta.p <- theta.now
      theta.p[k] <- r.theta.proposal(theta.now[k], theta.prop.sd[k], cat.theta[k])
      if(cat.theta[k]==0 && (theta.p[k]<0 || theta.p[k] > 1))
        next
     ## Recalculate basis and design mat for emulator at experimental obs w/new theta
      ans <- get.all.X(cbind(X,matrix(theta.p,nrow=nrow(X),ncol=length(theta.now),byrow=TRUE)), t.grid, Phi, P=P, max.ord=max.ord, include2, include3, include4, cat.vars=c(rep(0,ncol(X)), cat.theta))

      X.emd.p <- ans$X

      res.obs.prop <- yc.now - as.numeric(X.emd.p%*%alpha.now) - as.numeric(X.dis%*%beta.now)
      if(cat.theta[k]!=0 && it>N.init){    ### then update beta/Gamma in addition to residuals
        ans.beta <- update.betas(res.obs.prop, beta.now, X.dis, col.ind.dis, ncomp.dis, dhat.mat, Gamma.now, Sigma.now)
        beta.p <- ans.beta$beta.now
        res.obs.prop <- ans.beta$res.obs.now
        dhat.p <- ans.beta$dhat.mat
        dbeta.p.g.now <- ans.beta$dens
        dtheta.p.g.now <- d.theta.proposal(theta.p[k], theta.now[k], theta.prop.sd[k], cat.theta[k])
        dp.g.now <- dbeta.p.g.now + dtheta.p.g.now
        like.prop <- get.like.theta(res.obs.prop, Sigma.now)
        pi.prop <- pi.theta.t(theta.p)+pi.beta(beta.p,Gamma.now,col.ind.dis)
        

        res.obs.g.p <- yc.now - as.numeric(X.emd%*%alpha.now) - as.numeric(X.dis%*%beta.p)
        dbeta.now.g.p <- dprop.betas(res.obs.g.p, beta.p, beta.now, X.dis, col.ind.dis, ncomp.dis, dhat.p, Gamma.now, Sigma.now)
        dtheta.now.g.p <- d.theta.proposal(theta.now[k], theta.p[k], theta.prop.sd[k], cat.theta[k])
        dnow.g.p <- dbeta.now.g.p + dtheta.now.g.p
        like.now <- get.like.theta(res.obs.now, Sigma.now)
        pi.now <- pi.theta.t(theta.now)+pi.beta(beta.now,Gamma.now,col.ind.dis)

        MH.ratio <- exp(like.prop+pi.prop+dnow.g.p-like.now-pi.now-dp.g.now)
      }
      else{    ### just update residuals
        beta.p <- beta.now
        dhat.p <- dhat.mat
        like.now <- get.like.theta(res.obs.now, Sigma.now)
        like.prop <- get.like.theta(res.obs.prop, Sigma.now)
        pi.now <- pi.theta.t(theta.now)
        pi.prop <- pi.theta.t(theta.p)
        dnow.g.p <- d.theta.proposal(theta.now[k], theta.p[k], theta.prop.sd[k], cat.theta[k])
        dp.g.now <- d.theta.proposal(theta.p[k], theta.now[k], theta.prop.sd[k], cat.theta[k])
        MH.ratio <- exp(like.prop+pi.prop+dnow.g.p-like.now-pi.now-dp.g.now)
      }

      if(!is.na(MH.ratio) && runif(1) < MH.ratio){
        theta.now[k] <- theta.p[k]
        accept.count[it, k] <- 1
        res.obs.now <- res.obs.prop
        X.emd <- X.emd.p
        beta.now <- beta.p
        dhat.mat <- dhat.p
      }
    }

   ## create new matrix of component yhatd's for new thetas
    for(j in 1:ncomp.em){
      cols.em <- col.ind.em[[j]]
      yhatd.mat[,j,] <- as.numeric(X.emd[,cols.em]%*%alpha.now[cols.em,])
    }



   ## record params.now in posterior sample
    alpha[it,,] <- alpha.now
    beta[it,,] <- beta.now
    Lambda[it,,,] <- Lambda.now
    Gamma[it,,,] <- Gamma.now
    Ups[it,,] <- Ups.now
    Sigma[it,,] <- Sigma.now
    theta[it,] <- theta.now
    res.obs[it,,] <- res.obs.now
    res.sim[it,,] <- res.sim.now
    res.obs.sim[it,,] <- y-X.emd%*%alpha.now
    y.complete[it,,] <- yc.now


#X.dis..<<-X.dis
#alpha.now..<<-alpha.now
#beta.now..<<-beta.now
#X.emd..<<-X.emd
#Lambda.now..<<-Lambda.now
#Gamma.now..<<-Gamma.now
#Ups.now..<<-Ups.now
#Sigma.now..<<-Sigma.now
#theta.now..<<-theta.now

   ## Plot posterior
    if(it%%nplot==0){

      ys.lim <- apply(ys,2,range)
      N.plots <- M*((J+K+1)+(J+1)+5)+K
      cols <- min(8, ceiling(sqrt(N.plots)))
      rows <- min(8, ceiling(N.plots/cols))
      ind.now <- max(floor(it/2),it-nback+1):it
      par(mfrow=c(rows,cols), mar=c(2,2,2,1))
      ind.curves <- sample(ind.now,min(20,length(ind.now)))

      res.obs.now.act <- res.obs.now
      res.obs.now.act[ind.na.obs] <- NA

      Rsq <- 1-colSums(res.obs.now.act^2, na.rm=TRUE)/SSTo
      Rsq.sim.only <- 1-colSums(as.matrix(res.obs.sim[it,,]^2), na.rm=TRUE)/SSTo
      Rsq.s <- 1-colSums(res.sim.now^2)/SSTo.s

      res.obs.act <- array(res.obs[ind.now,,],c(length(ind.now),I, M))
      for(i in 1:length(ind.now)){
        blah <- res.obs.act[i,,]
        blah[ind.na.obs] <- NA
        res.obs.act[i,,] <- blah
      }
      Rsq.vec <- 1-apply(array(res.obs.act^2,c(length(ind.now),I,M)), c(1,3),sum,na.rm=TRUE)/matrix(SSTo,nrow=length(ind.now),ncol=M,byrow=TRUE)
      Rsq.sim.only.vec <-  1-apply(array(res.obs.sim[ind.now,,]^2,c(length(ind.now),I,M)),c(1,3),sum, na.rm=TRUE)/matrix(SSTo,nrow=length(ind.now),ncol=M,byrow=TRUE)
      Rsq.s.vec <- 1-apply(array(res.sim[ind.now,,]^2,c(length(ind.now),L,M)),c(1,3),sum)/matrix(SSTo.s,nrow=length(ind.now),ncol=M,byrow=TRUE)

      ys.hat <- yc.now - as.matrix(res.obs.sim[it,,])
      y.hat <- yc.now - res.obs.now.act
      Rel.err.sim <- colMeans(abs(1-abs(ys.hat)/(abs(y.hat)+1E-3)), na.rm=TRUE)
      Rel.max.sim <- apply(abs(1-abs(ys.hat)/(abs(y.hat)+1E-3)),2,max, na.rm=TRUE)



      cat("\n\nExp Rsq =",Rsq)
      cat("\nExp Rsq (recent avg) =",colMeans(Rsq.vec))
#      cat("\nExp Rsq (simulator only) =", Rsq.sim.only)
#      cat("\nExp Rsq (sim only recent avg) =",colMeans(Rsq.sim.only.vec))
      cat("\nEmulator Rsq =",Rsq.s)
      cat("\nEmulator Rsq (recent avg) =",colMeans(Rsq.s.vec))
      cat("\nSimulator Avg Rel Error =", Rel.err.sim)
      cat("\nSimulator Max Rel Error =", Rel.max.sim)
      cat("\naccept.pct = ", colMeans(accept.count[1:it,]),"\n\n")


 
     ## Plot Thetas
      for(k in 1:K){
        if(cat.theta[k]==0)
          plot(theta[ind.now,k]*diff(T.lim[k,])+T.lim[k,1],ylab="",main=paste("Theta.",k,sep=""),ylim=T.lim[k,],cex=.5)
        else
          plot(theta[ind.now,k],ylab="",main=paste("Theta.",k,sep=""),ylim=c(1,cat.theta[k]),cex=.5)
      }

      for(m in 1:M){
       ## Plot alpha0
        plot(alpha[ind.now,1,m],main=paste("Alpha.",0,".",m,sep=""),cex=.5)

       ## Plot Emulator Main Effect Components ...
        for(j in 1:(J+K)){
 
          if(j<=J){
            x.j <- Xs[,j]*(X.lim[j,2]-X.lim[j,1])+X.lim[j,1]
            plot(x.j, ys[,m], main=paste("eta.",j,".",m," by x.",j,sep=""),ylim=ys.lim[,m])  
          }
          else{
            x.j <- Ts[,j-J]*(T.lim[j-J,2]-T.lim[j-J,1])+T.lim[j-J,1]
            plot(x.j, ys[,m], main=paste("eta.",j,".",m," by theta.",j-J,sep=""),ylim=ys.lim[,m])  
          }
          ord.j <- order(x.j)
          X.j <- X.em[,c(1,col.ind.em[[j]])]
          for(i in ind.curves){
            alpha.j <- alpha[i,c(1,col.ind.em[[j]]),m]
            y.hat <- as.numeric(X.j%*%alpha.j)
            lines(x.j[ord.j],y.hat[ord.j],col=grey(.4))
          }
        }

       ## Plot beta0
        plot(beta[ind.now,1,m],main=paste("Beta.",0,".",m,sep=""),cex=.5)


       ## Plot Sim + Discrepancy Main Effect Components ...
        ind.na.m <- ind.na.obs[,m]
        for(j in 1:J){
          x.j <- X[,j]*(X.lim[j,2]-X.lim[j,1])+X.lim[j,1]
          yd <- yc.now[,m]
          R.yd <- diff(range(yd))
          yd.lim <- c(min(yd)-R.yd*.2, max(yd)+R.yd*.2)
          plot(x.j, yd, main=paste("Sim+Delta.",j,".",m," by x.",j,sep=""), ylim=yd.lim, pch='x', col=0)  
          ord.j <- order(x.j)
          for(i in ind.curves){
            ans <- get.all.X(cbind(X,matrix(theta[i,],nrow=nrow(X),ncol=K,byrow=TRUE)), t.grid, Phi, P=P, max.ord=max.ord, include2, include3, include4, cat.vars=c(rep(0,ncol(X)), cat.theta))
            Xem.i <- ans$X
#Xem.i..<<-Xem.i
            y.sim <- Xem.i[,c(1,col.ind.em[[j]])]%*%alpha[i,c(1,col.ind.em[[j]]),m]
            y.hat <- y.sim + X.dis[,c(1,col.ind.dis[[j]])]%*%beta[i,c(1,col.ind.dis[[j]]),m]
            lines(x.j[ord.j],y.hat[ord.j],col=2,cex=.5)
            lines(x.j[ord.j],y.sim[ord.j],col=4,cex=.5)
          }
          points(x.j[!ind.na.m], yd[!ind.na.m], pch='x', col=1, cex=1.5)  
          points(x.j[ind.na.m], yd[ind.na.m], pch='x', col=grey(.5), cex=1.5)  

        }

       ## Plot Ups
        plot(Ups[ind.now,m,m],main=paste("Ups.",m,sep=""),cex=.5)

       ## Plot Sigma
        plot(Sigma[ind.now,m,m],main=paste("Sigma.",m,sep=""),cex=.5)


       ## Plot Fitted Exp Data against observed
        y.hat <- rep(y[,m],times=length(ind.curves))-as.numeric(res.obs[ind.curves,,m])
        plot(rep(y[!ind.na.m,m],length(ind.curves)),y.hat[!ind.na.m],main=paste("Y.",m,"Fitted w/Discrep (Exp Data)",sep=""),cex=.5)
        abline(0,1)

       ## Plot Fitted (simulator Only) Exp Data against observed
        y.hat <- rep(y[,m],times=length(ind.curves))-as.numeric(res.obs.sim[ind.curves,,m])
        plot(rep(y[!ind.na.m,m],length(ind.curves)),y.hat[!ind.na.m],main=paste("Y.",m,"Fitted Simulator (Exp Data)",sep=""),cex=.5)
        abline(0,1)

       ## Plot Fitted Sim Data against observed runs
        ys.hat <- ys[,m]-res.sim.now[,m]
        lim <- range(ys[,m], ys.hat)
        plot(ys.hat,ys[,m],main=paste("Y.",m,"Fitted (Simulator)",sep=""),cex=.5,xlim=lim,ylim=lim)
        abline(0,1)
      }
    }

    
  }
  par(mfrow=c(1,1))

  return(list(alpha=alpha, beta=beta, Lambda=Lambda, Gamma=Gamma, Ups=Ups, Sigma=Sigma, theta=theta, y=y, ys=ys, X=X, Xs=Xs, Ts=Ts, Phi=Phi, t.grid=t.grid, P=P, max.ord=max.ord, include2=include2, include3=include3, include4=include4, P.Lambda=P.Lambda, v.Lambda=v.Lambda, P.Gamma=P.Gamma, v.Gamma=v.Gamma, P.Ups=P.Ups, v.Ups=v.Ups, P.Sigma=P.Sigma, v.Sigma=v.Sigma, vars.basis.ind.dis=vars.basis.ind.dis, vars.basis.ind.em=vars.basis.ind.em, X.em=X.em, X.emd=X.emd, X.dis=X.dis, res.obs=res.obs, res.sim=res.sim, res.obs.sim=res.obs.sim, pi.theta=pi.theta, X.lim=X.lim, T.lim=T.lim, cat.theta=cat.theta, transform=transform, inv.transform=inv.transform, col.ind.dis=col.ind.dis, col.ind.em=col.ind.em, y.complete=y.complete, theta.names=theta.names, cat.names=cat.names, x.names=x.names, y.names=y.names, pow=pow))
}









###########################################################################
###################### Posterior Plotting Functions #######################
###########################################################################

##### Plot maringal posterior of Theta #####

plot.thetas <- function(ans.calibration, post.ind=NULL, theta.names=NULL, cat.names=NULL, n=200, N=100, theta.true=NULL, plot.prior=TRUE, rows=NULL, cols=NULL, mar=NULL){

  lcol <- ifelse(plot.prior,4,0)
  n.MCMC <- nrow(ans.calibration$theta)
  T.lim <- ans.calibration$T.lim
  K <- ncol(ans.calibration$theta)
  cat.theta <- ans.calibration$cat.theta
  pi.theta <- ans.calibration$pi.theta

  if(is.null(rows)){
    cols <- min(8, ceiling(sqrt(K)))
    rows <- min(8, ceiling(K/cols))
  }

  if(is.null(post.ind))
    post.ind <- floor(.5*n.MCMC):n.MCMC
  if(is.null(theta.names))
    theta.names <- ans.calibration$theta.names
  if(is.null(cat.names))
    cat.names <- ans.calibration$cat.names
  theta.post <- ans.calibration$theta[post.ind,]
  for(k in 1:K)
    theta.post[,k] <- theta.post[,k]*diff(T.lim[k,])+T.lim[k,1]

  par(mfrow=c(rows,cols))
  if(!is.null(mar))
    par(mar=mar)

  for(k in 1:K){
    if(cat.theta[k]==0){
      x <- seq(T.lim[k,1],T.lim[k,2], length=n)
      if(plot.prior)
        p.x <- marginalize.pi(x, dim=k, pi.theta, cat.theta, T.lim, N=N)
      else 
        p.x <- rep(-1,length(x))
      p.hist <- hist(theta.post[,k],plot=FALSE)$density
      hist(theta.post[,k], ylab="density", main=theta.names[k], freq=FALSE, xlab="", xlim=T.lim[k,], ylim=c(0,max(p.x,p.hist)*1.1))
      lines(x,p.x, col=lcol)
      if(!is.null(theta.true)){
        par(col.axis=2)
        axis(side=1, at=theta.true[k], col=0, col.ticks=2, labels="",line=0, lwd.ticks=2)
        axis(side=1, at=theta.true[k], col=0, col.ticks=2, labels="",line=.5, lwd.ticks=2)
        axis(side=1, at=theta.true[k], col=0, col.ticks=2, labels="",line=1, lwd.ticks=2)
        axis(side=1, at=theta.true[k], col=0, col.ticks=2, labels=paste("theta_",k,sep=""),line=1.5, lwd.ticks=2)
        par(col.axis=1)

#        abline(v=theta.true[k], col=2)
      }
    }
    else{
      theta.k <- theta.post[,k]
      theta.k <- c(theta.k,1:cat.theta[k])
      props <- (table(theta.k)-1)/length(theta.k)
      x <- 1:cat.theta[k]
     if(plot.prior)
       p.x <- marginalize.pi(x, dim=k, pi.theta, cat.theta, T.lim, N=N)
      else 
        p.x <- rep(0,length(x))
      ylim <- c(0,max(p.x,props)*1.1)
      if(!is.null(theta.true)){
        names <- cat.names[[k]]
        names[theta.true[k]] <- ""
        barplot(props, names.arg=names, col=0, border=1, space=.4, main=theta.names[k], ylim=ylim)
        barplot(p.x, add=TRUE, col=0, border=lcol, space=.4, width=1, ylim=ylim)
        par(col.axis=2)
        names <- rep("", length(cat.names[[k]]))
        names[theta.true[k]] <- cat.names[[k]][theta.true[k]]
        barplot(p.x, add=TRUE, col=0, border=0, space=.4, names.arg=names, axes=F)
        par(col.axis=1)
      }
      else{
        barplot(props, names.arg=cat.names[[k]], col=0, border=1, space=.4, main=theta.names[k], ylim=ylim)
        barplot(p.x, add=TRUE, col=0, border=lcol, space=.4, width=1, ylim=ylim)
      }
    }
  }
}



marginalize.pi <- function(x, dim, pi.theta, cat.theta, T.lim, N){

  n <- length(x)
  K <- nrow(T.lim)
  X <- matrix(runif(K*N),N,K)
  ans <- matrix(0,n,N)
  pi.x <- rep(0,n)
  for(k in 1:K){
    if(cat.theta[k]==0)
      X[,k] <- X[,k]*diff(T.lim[k,])+T.lim[k,1]
    else
      X[,k] <- sample(1:cat.theta[k],N,replace=TRUE)
  }
  for(i in 1:length(x)){
    for(j in 1:N){
      theta.ij <- X[j,]
      theta.ij[dim] <- x[i]
      ans[i,j] <- exp(pi.theta(theta.ij))
      
    }
    pi.x[i] <- mean(ans[i,])
  }
  if(cat.theta[dim]==0)
    pi.x <- pi.x/(mean(pi.x)*diff(T.lim[dim,]))
  else
    pi.x <- pi.x/sum(pi.x)
  return(pi.x)
}

    


plot.thetas.bi <- function(ans.calibration, post.ind=NULL, theta.names=NULL, cat.names=NULL, theta.true=NULL){

  n.MCMC <- nrow(ans.calibration$theta)
  T.lim <- ans.calibration$T.lim
  K <- ncol(ans.calibration$theta)
  cat.theta <- ans.calibration$cat.theta

  if(is.null(post.ind))
    post.ind <- floor(.5*n.MCMC):n.MCMC
  if(is.null(theta.names))
    theta.names <- ans.calibration$theta.names
  if(is.null(cat.names))
    cat.names <- ans.calibration$cat.names
  theta.post <- ans.calibration$theta[post.ind,]
  for(k in 1:K)
    theta.post[,k] <- theta.post[,k]*diff(T.lim[k,])+T.lim[k,1]

  par(mfrow=c(K,K))
  for(k in 1:K){
    for(l in 1:K){
      theta.k <- theta.post[,k]
      theta.l <- theta.post[,l]
      title.kl <- paste(theta.names[k],"by",theta.names[l])

      if(k==l){
        if(cat.theta[k]==0)
          hist(theta.k,main=paste(theta.names[k]),cex.main=.8,prob=T)
        else{
          theta.k <- c(theta.k,1:cat.theta[k])
          props <- (table(theta.k)-1)/length(theta.k)
          barplot(props, names.arg=cat.names[[k]], col=0, border=1, space=.4, main=theta.names[k],cex.main=.8, cex.names=.8)
        }
      }
      else{
        if(cat.theta[k]==0 && cat.theta[l]==0)
          smoothScatter(theta.k,theta.l,main=title.kl, cex.main=.8)
        else
          plot(theta.k,theta.l,main=title.kl, cex.main=.8)
      }
    } 
  }
}





############ Plot Posterior Main Effect Simulator Fits & Discrepancy ##############




plot.fitted.2input <- function(ans.calibration, post.ind=NULL, input1=1, input2=2, nbin2=6, bin2=NULL, output=1, ncurves=500, nplot=ncurves, ylim=NULL, legend.loc="topright", nx.curve=50, transpose=FALSE,y.names=NULL,x.names=NULL, alpha.true=NULL, beta.true=NULL, theta.true=NULL, mar=NULL, digits=0, xlim=NULL){

  sim.col <- c(4,rgb(.7,.7,1),5)
  dis.col <- c(2,rgb(1,.7,.7),6)
  y <- ans.calibration$y
  X <- ans.calibration$X
  X.dis <- ans.calibration$X.dis
  X.lim <- ans.calibration$X.lim
  T.lim <- ans.calibration$T.lim
  theta <- ans.calibration$theta
  alpha <- ans.calibration$alpha
  beta <- ans.calibration$beta
  t.grid <- ans.calibration$t.grid
  Phi <- ans.calibration$Phi
  P <- ans.calibration$P
  max.ord <- ans.calibration$max.ord
  include2 <- ans.calibration$include2
  include3 <- ans.calibration$include3
  include4 <- ans.calibration$include4
  cat.theta <- ans.calibration$cat.theta
  col.ind.em <- ans.calibration$col.ind.em
  col.ind.dis <- ans.calibration$col.ind.dis
  transform <- ans.calibration$transform
  inv.transform <- ans.calibration$inv.transform
  res.obs <- ans.calibration$res.obs
  y.complete <- ans.calibration$y.complete
  Sigma <- ans.calibration$Sigma
  pow <- ans.calibration$pow

  ind.na.obs <- is.na(y)
  I <- nrow(y)
  M <- ncol(y)
  J <- ncol(X)
  K <- ncol(theta)
  n.MCMC <- nrow(ans.calibration$theta)

  x2.vals <- sort(unique(X[,input2]))
  if(length(x2.vals) > nbin2){
    if(is.null(bin2)){
      bin2 <- quantile(X[,input2], probs=c(seq(0,1,length=nbin2+1)))
    }
    x2.vals <- (bin2[-1]+bin2[-(nbin2+1)])/2
    n.x2 <- length(x2.vals)
    for(i in 1:I){
      ind.i <- max(sum(bin2<X[i,input2]),1)
      X[i,input2] <- x2.vals[ind.i]
    }
  }

  cols <- min(8, ceiling(sqrt(n.x2)))
  rows <- min(8, ceiling(n.x2/cols))
  if(is.null(post.ind))
    post.ind <- floor(.5*n.MCMC):n.MCMC
  ncurves <- min(ncurves,length(post.ind))
  ind.curves <- sample(post.ind,ncurves)
  if(is.null(y.names))
    y.names <- ans.calibration$y.names
  if(is.null(x.names))
    x.names <- ans.calibration$x.names
  if(is.null(xlim))
    xlim <- X.lim[input1,]

  if(transpose)
    par(mfcol=c(cols,rows))
  else
    par(mfrow=c(rows,cols))

  if(!is.null(mar))
    par(mar=mar)
  m <- output
  ## get obs error bounds
  ind.na <- ind.na.obs[,m]
  Sigma.hat <- apply(Sigma[post.ind,,], c(2,3), mean)
  y.conf <- inv.transform[[m]](rbind(y[,m]-2*sqrt(Sigma.hat[m,m]), y[,m]+2*sqrt(Sigma.hat[m,m])), pow[m])
  y.conf2 <-  inv.transform[[m]](apply(y.complete[post.ind,,m], 2, quantile, c(.025, .975)), pow[m])
#  y.conf2 <- (apply(y.complete[post.ind,,m], 2, quantile, c(.025, .975)))
#  y.conf2[1,] <- inv.transform[[m]](y.conf2[1,]-2*sqrt(Sigma.hat[m,m]))
#  y.conf2[2,] <- inv.transform[[m]](y.conf2[2,]+2*sqrt(Sigma.hat[m,m]))
  y.conf[is.na(y.conf)] <- 0
  y.conf2[is.na(y.conf2)] <- 0
  y.conf[,ind.na] <- y.conf2[,ind.na]
  yc <- apply(y.complete[post.ind,,m],2,mean)

#y.conf..<<-y.conf


  if(!is.null(theta.true))
    theta.true <- (theta.true-T.lim[,1])/apply(T.lim,1,diff)

  it <- 0
  for(x2 in x2.vals){
    it <- it+1
    x2t <- round(x2*(X.lim[input2,2]-X.lim[input2,1])+X.lim[input2,1],digits)
    ind.x2 <- (X[,input2]==x2)
    ind.na.x2 <- is.na(y[ind.x2,m])
    x.j <- X[ind.x2,input1]*(X.lim[input1,2]-X.lim[input1,1])+X.lim[input1,1]
    yd <- inv.transform[[m]](yc[ind.x2], pow[m])
    y.confd <- y.conf[,ind.x2]
    ord.j <- order(x.j)
    if(input1==1)
      X.plot <- matrix(c(seq(0,1, length=nx.curve), rep(x2,nx.curve)), nrow=nx.curve, ncol=J)
    if(input1==2)
      X.plot <- matrix(c(rep(x2,nx.curve),seq(0,1, length=nx.curve)), nrow=nx.curve, ncol=J)
    ans <- get.all.X(X.plot, t.grid, Phi, P=P, max.ord=max.ord, include2, include3, include4, cat.vars=rep(0,ncol(X.plot)))
    X.dis <- ans$X
    xj.plot <- X.plot[,input1]*(X.lim[input1,2]-X.lim[input1,1])+X.lim[input1,1]
    y.sim <- y.hat <- yt.sim <- yt.hat <- matrix(0, nrow(X.plot), ncurves)

    for(i in 1:length(ind.curves)){
      ans <- get.all.X(cbind(X.plot,matrix(theta[ind.curves[i],],nrow=nrow(X.plot),ncol=K,byrow=TRUE)), t.grid, Phi, P=P, max.ord=max.ord, include2, include3, include4, cat.vars=c(rep(0,ncol(X.plot)), cat.theta))
      Xem.i <- ans$X
      col.ind <- ans$col.ind
      yt.sim[,i] <- Xem.i%*%alpha[ind.curves[i],,m]
      yt.hat[,i] <- yt.sim[,i] + X.dis%*%beta[ind.curves[i],,m]
      y.hat[,i] <- inv.transform[[m]](yt.hat[,i], pow[m])
      y.sim[,i] <- inv.transform[[m]](yt.sim[,i], pow[m])
    }
    y.sim.mean <- inv.transform[[m]](apply(yt.sim,1,mean), pow[m])
    y.hat.mean <- inv.transform[[m]](apply(yt.hat,1,mean), pow[m])
    if(!is.null(alpha.true)){
      ans <- get.all.X(cbind(X.plot,matrix(theta.true,nrow=nrow(X.plot),ncol=K,byrow=TRUE)), t.grid, Phi, P=P, max.ord=max.ord, include2, include3, include4, cat.vars=c(rep(0,ncol(X.plot)), cat.theta))
      Xem.i <- ans$X
      y.sim.true <- Xem.i%*%alpha.true[,m]
      y.dis.true <- X.dis%*%beta.true[,m]
      y.hat.true <- inv.transform[[m]](y.sim.true + y.dis.true, pow[m])
      y.sim.true <- inv.transform[[m]](y.sim.true, pow[m])
    }

    if(is.null(ylim))
      yd.lim <- range(0,y.sim,y.hat,y.confd, na.rm=TRUE)
    else
      yd.lim <- ylim[it,]
    plot(x.j, yd, main=paste("Emulator+Discrep ",x.names[input2],"=",format(x2t,digits=1, nsmall=1),sep=""), ylim=yd.lim, pch='x', col=0, ylab="",xlab="",xlim=xlim)
    title(ylab=y.names[m], xlab=x.names[input1],line=2)  
    if(!is.null(alpha.true)){
      legend(legend.loc, legend=c("Simulator Truth","Emulator Mean", "Emulator Reals","Physical Reality", "Emu + Dis Mean", "Emu + Dis Reals", "Observed Data", "Missing Data"), lty=c(3,1,1,3,1,1,1,1), pch=c(NA,NA,NA,NA,NA,NA,'X','X'), lwd=c(2,2,.5,2,2,.5,1,1), col=c(sim.col[c(3,1:2)],dis.col[c(3,1:2)],1,grey(.5)), cex=.9)
    }
    else{
      legend(legend.loc, legend=c("Emulator Mean","Emulator Reals","Emu + Dis Mean", "Emu + Dis Reals", "Observed Data"), lty=c(1,1,1,1,1,1), pch=c(NA,NA,NA,NA,'X','X'), lwd=c(2,.5,2,.5,1,1), col=c(sim.col[1:2],dis.col[1:2],1,grey(.5)), cex=.9)
    }
    for(i in 1:nplot){
      lines(xj.plot,y.hat[,i],col=dis.col[2],lwd=.5)
      lines(xj.plot,y.sim[,i],col=sim.col[2],lwd=.5)
    }
    lines(xj.plot,y.hat.mean,col=dis.col[1],lwd=2)
    lines(xj.plot,y.sim.mean,col=sim.col[1],lwd=2)
    pt.col <-ifelse(ind.na.x2, grey(.5), 1)
    points(x.j, yd, pch='x', col=pt.col, cex=1.5)
    if(!is.null(alpha.true)){
      lines(xj.plot,y.hat.true,col=dis.col[3],lwd=2, lty=2)
      lines(xj.plot,y.sim.true,col=sim.col[3],lwd=2, lty=2)
    }
    for(i in 1:length(x.j)){
      lines(rep(x.j[i],2), y.confd[,i],col=pt.col[i],lwd=2)
     
    }
  }
}






plot.fitted.2input.cv <- function(ans.calibration, ans.cv, x.locs=NULL, post.ind=NULL, ncurves=500,transpose=FALSE,y.names=NULL,x.names=NULL, nx.curve=50, input1=1, input2=2, output=1, alpha.true=NULL, beta.true=NULL, theta.true=NULL, ylim=NULL, legend.loc="topright", conf=.90, df=15){

  sim.col <- c(4,rgb(.5,.5,1),5)
  dis.col <- c(2,rgb(1,.5,.5),6)
  y <- ans.calibration$y
  X <- ans.calibration$X
  X.dis <- ans.calibration$X.dis
  X.lim <- ans.calibration$X.lim
  T.lim <- ans.calibration$T.lim
  theta <- ans.calibration$theta
  t.grid <- ans.calibration$t.grid
  Phi <- ans.calibration$Phi
  P <- ans.calibration$P
  max.ord <- ans.calibration$max.ord
  include2 <- ans.calibration$include2
  include3 <- ans.calibration$include3
  include4 <- ans.calibration$include4
  cat.theta <- ans.calibration$cat.theta
  col.ind.em <- ans.calibration$col.ind.em
  col.ind.dis <- ans.calibration$col.ind.dis
  transform <- ans.calibration$transform
  inv.transform <- ans.calibration$inv.transform
  res.obs <- ans.calibration$res.obs
  y.complete <- ans.calibration$y.complete
  Sigma <- ans.calibration$Sigma
  pow <- ans.calibration$pow

  ind.na.obs <- is.na(y)
  I <- nrow(y)
  M <- ncol(y)
  J <- ncol(X)
  K <- ncol(theta)
  n.MCMC <- nrow(ans.cv[[1]]$theta)

  x2.vals <- sort(unique(X[,input2]))
  n.x2 <- length(x2.vals)

  x.locs <- (x.locs-X.lim[input2,1])/diff(X.lim[input2,])
  alpha <- beta <- list()
  for(f in 1:n.x2){
    ind.f <- which(x.locs==x2.vals[f])
    if(length(ind.f)==0){
print(ind.f)
      alpha[[f]] <- ans.calibration$alpha
      beta[[f]] <- ans.calibration$beta
    }
    else{
print(ind.f)
      alpha[[f]] <- ans.cv[[ind.f]]$alpha[1:n.MCMC,,]
      beta[[f]] <- ans.cv[[ind.f]]$beta[1:n.MCMC,,]
    }
  }

  cols <- min(8, ceiling(sqrt(n.x2)))
  rows <- min(8, ceiling(n.x2/cols))
  if(is.null(post.ind))
    post.ind <- floor(.5*n.MCMC):n.MCMC
  ncurves <- min(ncurves,length(post.ind))
  ind.curves <- sample(post.ind,ncurves)
  if(is.null(y.names))
    y.names <- ans.calibration$y.names
  if(is.null(x.names))
    x.names <- ans.calibration$x.names

  if(transpose)
    par(mfcol=c(cols,rows))
  else
    par(mfrow=c(rows,cols))

  m <- output
  ## get obs error bounds
  ind.na <- ind.na.obs[,m]
  Sigma.hat <- apply(Sigma[post.ind,,], c(2,3), mean)
  y.conf <- inv.transform[[m]](rbind(y[,m]-2*sqrt(Sigma.hat[m,m]), y[,m]+2*sqrt(Sigma.hat[m,m])), pow[m])
  y.conf2 <- (apply(y.complete[post.ind,,m], 2, quantile, c(.025, .975)))
  y.conf2[1,] <- inv.transform[[m]](y.conf2[1,]-2*sqrt(Sigma.hat[m,m]), pow[m])
  y.conf2[2,] <- inv.transform[[m]](y.conf2[2,]+2*sqrt(Sigma.hat[m,m]), pow[m])
  y.conf[is.na(y.conf)] <- 0
  y.conf2[is.na(y.conf2)] <- 0
  y.conf[,ind.na] <- y.conf2[,ind.na]
  yc <- apply(y.complete[post.ind,,m],2,mean)

  if(!is.null(theta.true))
    theta.true <- (theta.true-T.lim[,1])/apply(T.lim,1,diff)
  y.hat.all <- y.sim.all <- matrix(0,nx.curve,length(x2.vals))

  it <- 0
  for(x2 in x2.vals){
    it <- it+1
    x2t <- x2*(X.lim[input2,2]-X.lim[input2,1])+X.lim[input2,1]
    ind.x2 <- (X[,input2]==x2)
    ind.na.x2 <- is.na(y[ind.x2,m])
    x.j <- X[ind.x2,input1]*(X.lim[input1,2]-X.lim[input1,1])+X.lim[input1,1]
    yd <- inv.transform[[m]](yc[ind.x2], pow[m])
    y.confd <- y.conf[,ind.x2]
    ord.j <- order(x.j)
    if(input1==1)
      X.plot <- matrix(c(seq(0,1, length=nx.curve), rep(x2,nx.curve)), nrow=nx.curve, ncol=J)
    if(input1==2)
      X.plot <- matrix(c(rep(x2,nx.curve),seq(0,1, length=nx.curve)), nrow=nx.curve, ncol=J)
    ans <- get.all.X(X.plot, t.grid, Phi, P=P, max.ord=max.ord, include2, include3, include4, cat.vars=rep(0,ncol(X.plot)))
    X.dis <- ans$X
    xj.plot <- X.plot[,input1]*(X.lim[input1,2]-X.lim[input1,1])+X.lim[input1,1]
    y.sim <- y.hat <- yt.sim <- yt.hat <- matrix(0, nrow(X.plot), ncurves)

    for(i in 1:length(ind.curves)){
      ans <- get.all.X(cbind(X.plot,matrix(theta[ind.curves[i],],nrow=nrow(X.plot),ncol=K,byrow=TRUE)), t.grid, Phi, P=P, max.ord=max.ord, include2, include3, include4, cat.vars=c(rep(0,ncol(X.plot)), cat.theta))
      Xem.i <- ans$X
      col.ind <- ans$col.ind
      yt.sim[,i] <- Xem.i%*%alpha[[it]][ind.curves[i],,m]
      yt.hat[,i] <- yt.sim[,i] + X.dis%*%beta[[it]][ind.curves[i],,m]
      y.hat[,i] <- inv.transform[[m]](yt.hat[,i], pow[m])
      y.sim[,i] <- inv.transform[[m]](yt.sim[,i], pow[m])
    }
    sim.bands <- get.cred.bands(t(y.sim), alpha=1-conf)
    y.sim.conf <- sim.bands$band
    y.sim.conf[1,] <- smooth.spline(1:nx.curve,y.sim.conf[1,],df=df)$y
    y.sim.conf[2,] <- smooth.spline(1:nx.curve,y.sim.conf[2,],df=df)$y
    y.sim.mean <- sim.bands$mu
    hat.bands <- get.cred.bands(t(y.hat), alpha=1-conf)
    y.hat.conf <- hat.bands$band
    y.hat.conf[1,] <- smooth.spline(1:nx.curve,y.hat.conf[1,],df=df)$y
    y.hat.conf[2,] <- smooth.spline(1:nx.curve,y.hat.conf[2,],df=df)$y
    y.hat.mean <- hat.bands$mu

    y.hat.all[,it] <- y.hat.mean
    y.sim.all[,it] <- y.sim.mean


#    y.sim.mean <- inv.transform[[m]](apply(yt.sim,1,mean), pow[m])
#    y.hat.mean <- inv.transform[[m]](apply(yt.hat,1,mean), pow[m])
    if(!is.null(alpha.true)){
      ans <- get.all.X(cbind(X.plot,matrix(theta.true,nrow=nrow(X.plot),ncol=K,byrow=TRUE)), t.grid, Phi, P=P, max.ord=max.ord, include2, include3, include4, cat.vars=c(rep(0,ncol(X.plot)), cat.theta))
      Xem.i <- ans$X
      y.sim.true <- Xem.i%*%alpha.true[,m]
      y.dis.true <- X.dis%*%beta.true[,m]
      y.hat.true <- inv.transform[[m]](y.sim.true + y.dis.true, pow[m])
      y.sim.true <- inv.transform[[m]](y.sim.true, pow[m])
    }

    if(is.null(ylim))
      yd.lim <- range(y.sim.conf,y.hat.conf,y.confd, na.rm=TRUE)
    else
      yd.lim <- ylim[it,]
    plot(x.j, yd, main=paste("Simulator+Discrepancy, ",x.names[input2],"=",format(x2t,digits=1, nsmall=1),sep=""), ylim=yd.lim, pch='x', col=0, ylab=y.names[m], xlab=x.names[input1])  
    legend(legend.loc, legend=c("Emulator Mean",paste("Emulator ",conf*100,"% Bands",sep=""),"Sim + Dis Mean", paste("Sim + Dis ",conf*100,"% Bands",sep=""), "Data"), lty=c(1,2,1,2,1), pch=c(NA,NA,NA,NA,'X'), lwd=c(1.5,1.5,1.5,1.5,1), col=c(sim.col[c(1,1)],dis.col[c(1,1)],1), cex=.75)
#    for(i in 1:length(ind.curves)){
#      lines(xj.plot,y.hat[,i],col=dis.col[2],lwd=.5)
#      lines(xj.plot,y.sim[,i],col=sim.col[2],lwd=.5)
#    }
    lines(xj.plot, y.hat.conf[1,], col=dis.col, lty=2, lwd=1.5)
    lines(xj.plot, y.hat.conf[2,], col=dis.col, lty=2, lwd=1.5)
    lines(xj.plot, y.sim.conf[1,], col=sim.col, lty=2, lwd=1.5)
    lines(xj.plot, y.sim.conf[2,], col=sim.col, lty=2, lwd=1.5)
    lines(xj.plot,y.sim.mean,col=sim.col[1],lwd=2)
    lines(xj.plot,y.hat.mean,col=dis.col[1],lwd=2)
    pt.col <-ifelse(ind.na.x2, grey(.5), 1)
    points(x.j, yd, pch='x', col=pt.col, cex=1.5)
    if(!is.null(alpha.true)){
      lines(xj.plot,y.hat.true,col=dis.col[3],lwd=1.5, lty=2)
      lines(xj.plot,y.sim.true,col=sim.col[3],lwd=1.5, lty=2)
    }
    for(i in 1:length(x.j)){
      lines(rep(x.j[i],2), y.confd[,i],col=pt.col[i],lwd=2)
     
    }
  }
  return(list(y.hat=y.hat.all, y.sim=y.sim.all, xj.plot=xj.plot))
}











plot.fitted.all.output <- function(ans.calibration, post.ind=NULL, inputs=NULL, outputs=NULL, ncurves=500, nplot=ncurves, ylim=NULL,y.pred=y.pred ,legend.loc="topright", transpose=FALSE, y.names=NULL,x.names=NULL, nx.curve=50, mar=NULL, cex.leg=.8){

  sim.col <- c(4,rgb(.7,.7,1),5)
  dis.col <- c(2,rgb(1,.7,.7),6)
  y <- ans.calibration$y
  X <- ans.calibration$X
  X.dis <- ans.calibration$X.dis
  X.lim <- ans.calibration$X.lim
  T.lim <- ans.calibration$T.lim
  theta <- ans.calibration$theta
  alpha <- ans.calibration$alpha
  beta <- ans.calibration$beta
  t.grid <- ans.calibration$t.grid
  Phi <- ans.calibration$Phi
  P <- ans.calibration$P
  max.ord <- ans.calibration$max.ord
  include2 <- ans.calibration$include2
  include3 <- ans.calibration$include3
  include4 <- ans.calibration$include4
  cat.theta <- ans.calibration$cat.theta
  col.ind.em <- ans.calibration$col.ind.em
  col.ind.dis <- ans.calibration$col.ind.dis
  transform <- ans.calibration$transform
  inv.transform <- ans.calibration$inv.transform
  res.obs <- ans.calibration$res.obs
  y.complete <- ans.calibration$y.complete
  Sigma <- ans.calibration$Sigma
  pow <- ans.calibration$pow

  ind.na.obs <- is.na(y)
  I <- nrow(y)
  M <- ncol(y)
  J <- ncol(X)
  K <- ncol(theta)
  n.MCMC <- nrow(ans.calibration$theta)
  if(is.null(inputs))
    inputs <- 1:J
  if(is.null(outputs))
    outputs <- 1:M
  J.plot <- length(inputs)
  M.plot <- length(outputs)

  if(length(legend.loc)==1)
    legend.loc <- rep(legend.loc,M)

  if(J.plot >1){
    cols <- min(8, M.plot)
    rows <- min(8, J.plot)
  }
  else{
    cols <- min(8, ceiling(sqrt(M)))
    rows <- min(8, ceiling(M/cols))
  }
  if(is.null(post.ind))
    post.ind <- floor(.5*n.MCMC):n.MCMC
  ncurves <- min(ncurves,length(post.ind))
  ind.curves <- sample(post.ind,ncurves)
  if(is.null(y.names))
    y.names <- ans.calibration$y.names
  if(is.null(x.names))
    x.names <- ans.calibration$x.names

  if(transpose)
    par(mfcol=c(cols,rows))
  else
    par(mfrow=c(rows,cols))

  if(!is.null(mar))
    par(mar=mar)

for(j in 1:J.plot){
  X.plot <- matrix(.5, nrow=nx.curve, ncol=J)
  X.plot[,inputs[j]] <- seq(0,1, length=nx.curve)
  x.j <- X[,inputs[j]]*(X.lim[inputs[j],2]-X.lim[inputs[j],1])+X.lim[inputs[j],1]
  ord.j <- order(x.j)

  for(mp in 1:M.plot){
    m <- outputs[mp]
    ## get obs error bounds
    ind.na <- ind.na.obs[,m]
    Sigma.hat <- apply(array(Sigma[post.ind,,],c(length(post.ind),M,M)), c(2,3), mean)
    y.conf <- inv.transform[[m]](rbind(y[,m]-2*sqrt(Sigma.hat[m,m]), y[,m]+2*sqrt(Sigma.hat[m,m])), pow[m])
    y.conf2 <-  inv.transform[[m]](apply(y.complete[post.ind,,m], 2, quantile, c(.1, .9)), pow[m])
    y.conf[is.na(y.conf)] <- 0
    y.conf2[is.na(y.conf2)] <- 0
    y.conf[,ind.na] <- y.conf2[,ind.na]
    yc <- apply(y.complete[post.ind,,m],2,mean)
    yd <- inv.transform[[m]](yc, pow[m])
    ans <- get.all.X(X.plot, t.grid, Phi, P=P, max.ord=max.ord, include2, include3, include4, cat.vars=rep(0,ncol(X.plot)))
    X.dis <- ans$X
    xj.plot <- X.plot[,inputs[j]]*(X.lim[inputs[j],2]-X.lim[inputs[j],1])+X.lim[inputs[j],1]
    y.sim <- y.hat <- yt.sim <- yt.hat <- matrix(0, nrow(X.plot), ncurves)

    for(i in 1:length(ind.curves)){
      ans <- get.all.X(cbind(X.plot,matrix(theta[ind.curves[i],],nrow=nrow(X.plot),ncol=K,byrow=TRUE)), t.grid, Phi, P=P, max.ord=max.ord, include2, include3, include4, cat.vars=c(rep(0,ncol(X.plot)), cat.theta))
      Xem.i <- ans$X
      col.ind <- ans$col.ind
      yt.sim[,i] <- Xem.i%*%alpha[ind.curves[i],,m]
      yt.hat[,i] <- yt.sim[,i] + X.dis%*%beta[ind.curves[i],,m]
      y.hat[,i] <- inv.transform[[m]](yt.hat[,i], pow[m])
      y.sim[,i] <- inv.transform[[m]](yt.sim[,i], pow[m])
    }

    y.sim.mean <- inv.transform[[m]](apply(yt.sim,1,mean), pow[m])
    y.hat.mean <- inv.transform[[m]](apply(yt.hat,1,mean), pow[m])

    if(is.null(ylim))
      yd.lim <- range(0,y.sim,y.hat,y.conf, na.rm=TRUE)
    else
      yd.lim <- ylim[m,]
      ##
      cat("yd.lim",yd.lim,"\n")
      ##
    plot(x.j, yd, main=paste("Emulator, ",y.names[m],sep=""), ylim=yd.lim, pch='x', col=0, ylab="",xlab="")
    title(ylab=y.names[m], xlab=x.names[inputs[j]],line=2)  
    for(i in 1:nplot){
    #  lines(xj.plot,y.hat[,i],col=dis.col[2],lwd=.5)
      lines(xj.plot,y.sim[,i],col=sim.col[2],lwd=.5)
    }
    #lines(xj.plot,y.hat.mean,col=dis.col[1],lwd=2)
    lines(xj.plot,y.sim.mean,col=sim.col[1],lwd=2)
    pt.col <-ifelse(ind.na, grey(.5), 1)
    points(x.j, yd, pch='x', col=pt.col, cex=1.5)
    points(x.j,y.pred[m,],pch='o',col=dis.col[1],cex=1.5)
    ##
    #cat("yd",yd,"\n")
    ##
    for(i in 1:length(x.j)){
      lines(rep(x.j[i],2), y.conf[,i],col=pt.col[i],lwd=2)
    }
    legend(legend.loc[m], legend=c("Emulator Mean","Emulator Reals", "Observed Data","Predict Data"), lty=c(1,1,1,1), pch=c(NA,NA,'X','O'), lwd=c(2,.5,1,1), col=c(sim.col[1:2],1,dis.col[1]), cex=cex.leg)
  }
}
}









############ Plot Posterior Discrepancy ##############


plot.discrepancy <- function(ans.calibration, post.ind=NULL, nplot=50, ylim=NULL, legend.loc="topright", transpose=FALSE, y.names=NULL, x.names=NULL, Nit=NULL, nx.curve=25, alpha.true=NULL, beta.true=NULL, theta.true=NULL, df=12, rows=NULL, cols=NULL,cex.leg=.8){

  mu.col <- 4
  true.col <- 5
  cb.col <- 2
  dis.col <- grey(.7)
  y <- ans.calibration$y
  X <- ans.calibration$X
  X.dis <- ans.calibration$X.dis
  X.lim <- ans.calibration$X.lim
  T.lim <- ans.calibration$T.lim
  theta <- ans.calibration$theta
  alpha <- ans.calibration$alpha
  beta <- ans.calibration$beta
  t.grid <- ans.calibration$t.grid
  Phi <- ans.calibration$Phi
  P <- ans.calibration$P
  max.ord <- ans.calibration$max.ord
  include2 <- ans.calibration$include2
  include3 <- ans.calibration$include3
  include4 <- ans.calibration$include4
  cat.theta <- ans.calibration$cat.theta
  col.ind.dis <- ans.calibration$col.ind.dis
  transform <- ans.calibration$transform
  inv.transform <- ans.calibration$inv.transform
  res.obs <- ans.calibration$res.obs
  res.obs.sim <- ans.calibration$res.obs.sim
  pow <- ans.calibration$pow

  I <- nrow(y)
  M <- ncol(y)
  J <- ncol(X)
  K <- ncol(theta)
  if(is.null(rows)){
    rows <- M
    cols <- J
  }
  n.MCMC <- dim(beta)[1]
  if(is.null(Nit))
    Nit <- floor(n.MCMC/2)
  if(is.null(post.ind))
    post.ind <- sample(floor(.5*n.MCMC):n.MCMC, Nit)
  ind.curves <- sample(1:length(post.ind),nplot)
  if(is.null(y.names))
    y.names <- ans.calibration$y.names
  if(is.null(x.names))
    x.names <- ans.calibration$x.names
############
cols = 2
rows = 2
############
  if(transpose)
    par(mfcol=c(cols,rows))
  else
    par(mfrow=c(rows,cols))

  if(!is.null(theta.true))
    theta.true <- (theta.true-T.lim[,1])/apply(T.lim,1,diff)
  if(length(legend.loc)==1)
    legend.loc<-rep(legend.loc, J*M)

  X <- matrix(seq(0,1, length=nx.curve), nrow=nx.curve, ncol=J)
  ans <- get.all.X(X, t.grid, Phi, P=P, max.ord=max.ord, include2, include3, include4, cat.vars=rep(0,ncol(X)))
  X.dis <- ans$X

  d.post <- array(0,c(length(post.ind),nx.curve,J,M))
  if(!is.null(alpha.true))
    d.true <- array(0,c(nx.curve,J,M))
  for(i in 1:length(post.ind)){

    ans <- get.all.X(cbind(X,matrix(theta[post.ind[i],],nrow=nrow(X),ncol=K,byrow=TRUE)), t.grid, Phi, P=P, max.ord=max.ord, include2, include3, include4, cat.vars=c(rep(0,ncol(X)), cat.theta))
    Xem.i <- ans$X
    col.ind <- ans$col.ind

    for(m in 1:M){
      for(j in 1:J){
        y.sim <- Xem.i[,c(1,col.ind[[j]])]%*%alpha[post.ind[i],c(1,col.ind[[j]]),m]
        y.hat <- y.sim + X.dis[,c(1,col.ind.dis[[j]])]%*%beta[post.ind[i],c(1,col.ind.dis[[j]]),m]

        d.post[i,,j,m] <- inv.transform[[m]](y.hat, pow[m]) - inv.transform[[m]](y.sim, pow[m])
        if(any(is.na(d.post[i,,j,m])))
          d.post[i,,j,m] <- NA

        if(i==1 && !is.null(alpha.true)){
          ans <- get.all.X(cbind(X,matrix(theta.true,nrow=nrow(X),ncol=K,byrow=TRUE)), t.grid, Phi, P=P, max.ord=max.ord, include2, include3, include4, cat.vars=c(rep(0,ncol(X)), cat.theta))
          Xem.i <- ans$X
          y.sim.true <- Xem.i[,c(1,col.ind[[j]])]%*%alpha.true[c(1,col.ind[[j]]),m]
          y.hat.true <- y.sim.true + X.dis[,c(1,col.ind.dis[[j]])]%*%beta.true[c(1,col.ind.dis[[j]]),m]
          d.true[,j,m] <- inv.transform[[m]](y.hat.true, pow[m]) - inv.transform[[m]](y.sim.true, pow[m])
        }
      }
    }
  }

  it <- 0
  for(m in 1:M){
    for(j in 1:J){
      it <- it+1
      x.j <- X[,j]*(X.lim[j,2]-X.lim[j,1])+X.lim[j,1]
      d.post.jm <- d.post[!is.na(d.post[,1,j,m]),,j,m]
      ans.bands <- get.cred.bands(d.post.jm, alpha=.05)
      djm.conf <- ans.bands$band
      djm.conf[1,] <- smooth.spline(1:nx.curve,djm.conf[1,],df=df)$y
      djm.conf[2,] <- smooth.spline(1:nx.curve,djm.conf[2,],df=df)$y
      djm.mu <- ans.bands$mu
      if(is.null(ylim))
        ylim.jm <- range(0,d.post.jm)
      else{
        if(is.matrix(ylim))
          ylim.jm <- ylim[m,]
        else
          ylim.jm <- ylim
      }
#      plot(x.j, rep(0,length(x.j)), main=paste("Delta.",j,".",m," by x.",j,sep=""), ylim=ylim.jm, col=0, ylab=paste(y.names[m], "Discrepancy"), xlab=x.names[j])  
      plot(x.j, rep(0,length(x.j)), main="", ylim=ylim.jm, col=0, ylab=paste(y.names[m], "Discrepancy"), xlab=x.names[j])  
      if(!is.null(alpha.true)){
        legend(legend.loc[it], legend=c("True Discrepancy", "Posterior Mean", "95% Bands","Zero Line", "Posterior Realizations"), lty=c(3,1,1,2,1), lwd=c(2,2,2,2,.5), col=c(true.col,mu.col,cb.col,1,dis.col), cex=cex.leg)
      }
      else{
        legend(legend.loc[it], legend=c("Posterior Mean", "95% Bands","Zero Line", "Posterior Realizations"), lty=c(1,1,2,1), lwd=c(2,2,2,.5), col=c(mu.col,cb.col,1,dis.col), cex=cex.leg)
      }

      ord.j <- order(x.j)
      for(i in ind.curves){
        lines(x.j[ord.j],d.post[i,ord.j,j,m],col=dis.col,lwd=.5)
      }
      lines(x.j[ord.j], djm.mu[ord.j], col=mu.col, lwd=2)
      lines(x.j[ord.j], djm.conf[1,ord.j], col=cb.col, lwd=2)
      lines(x.j[ord.j], djm.conf[2,ord.j], col=cb.col, lwd=2)
      lines(x.j[ord.j], rep(0,length(x.j)), col=1, lty=2, lwd=2)
      if(!is.null(alpha.true)){
        lines(x.j[ord.j], d.true[ord.j,j,m], col=true.col, lwd=2, lty=2)
      }
    }
  }
}






get.cred.bands <- function(Y, alpha=.05, tol=1E-4){

  sim.conf <- function(gamma, Y){
    mu.Y <- colMeans(Y)
    sd.Y <- apply(Y,2,sd)
    band <- rbind(qnorm(gamma/2,mu.Y,sd.Y),qnorm(1-gamma/2,mu.Y,sd.Y))

    n.outside <- 0
    for(i in 1:nrow(Y))
      n.outside <- n.outside + any(Y[i,]<band[1,] | Y[i,]>band[2,])
    return((n.outside/nrow(Y)-alpha)^2)
  }

  gamma <- robust.optimize(fn=sim.conf,par.lim=c(0,alpha),npar=3,rel.tol=tol,Y=Y)$par
  band <- apply(Y,2,quantile,prob=c(gamma/2,1-gamma/2))
  mu <- colMeans(Y)
  return(list(mu=mu, band=band))
}








#################################################################################
####################### Predict new Simulator Runs ##############################
#################################################################################


predict.new <- function(x.new, ans.calibration, post.ind=NULL, Nit=NULL, trans.scale=TRUE){

  X.lim <- ans.calibration$X.lim
  theta <- ans.calibration$theta
  alpha <- ans.calibration$alpha
  beta <- ans.calibration$beta
  t.grid <- ans.calibration$t.grid
  Phi <- ans.calibration$Phi
  P <- ans.calibration$P
  max.ord <- ans.calibration$max.ord
  include2 <- ans.calibration$include2
  include3 <- ans.calibration$include3
  include4 <- ans.calibration$include4
  cat.theta <- ans.calibration$cat.theta
  col.ind.dis <- ans.calibration$col.ind.dis
  inv.transform <- ans.calibration$inv.transform
  pow <- ans.calibration$pow

  I <- nrow(x.new)
  M <- dim(beta)[3]
  J <- ncol(x.new)
  K <- ncol(theta)
  n.MCMC <- dim(beta)[1]
  if(is.null(Nit))
    Nit <- floor(n.MCMC/2)

  Sigma.hat <- apply(ans.calibration$Sigma[floor(.5*n.MCMC):n.MCMC,,],c(2,3),mean)
  if(is.null(post.ind))
    post.ind <- sample(floor(.5*n.MCMC):n.MCMC, Nit)

 ## Transform x.new to [0,1]
  X <- x.new
  for(j in 1:J)
    X[,j] <- (X[,j]-X.lim[j,1])/(X.lim[j,2]-X.lim[j,1])
  

 ## get design mat for discrepancy
  ans <- get.all.X(X, t.grid, Phi, P=P, max.ord=max.ord, include2, include3, include4, cat.vars=rep(0,ncol(X)))
  X.dis <- ans$X
  col.ind.dis <- ans$col.ind


  d.post <- array(0,c(length(post.ind),I,M))
  yem.post <- array(0,c(length(post.ind),I,M))
  for(i in 1:length(post.ind)){

    ans <- get.all.X(cbind(X,matrix(theta[post.ind[i],],nrow=nrow(X),ncol=K,byrow=TRUE)), t.grid, Phi, P=P, max.ord=max.ord, include2, include3, include4, cat.vars=c(rep(0,ncol(X)), cat.theta))
    Xem.i <- ans$X
    col.ind.em <- ans$col.ind
    yem.post[i,,] <- Xem.i%*%alpha[post.ind[i],,]
    d.post[i,,] <- X.dis%*%beta[post.ind[i],,]
  }
  y.hat <- (yem.post +  d.post)
  if(trans.scale){
    for(m in 1:M){
      y.hat[,,m] <- inv.transform[[m]](y.hat[,,m], pow[m])
      yem.post[,,m] <- inv.transform[[m]](yem.post[,,m], pow[m])
    }
  }
  return(list(y.hat=y.hat, y.sim=yem.post, theta=theta[post.ind,], Sigma.hat=Sigma.hat))
}






predict.sim <- function(x.new, t.new, ans.calibration, post.ind=NULL, Nit=1000, trans.scale=TRUE){

  X.lim <- ans.calibration$X.lim
  T.lim <- ans.calibration$T.lim
  alpha <- ans.calibration$alpha
  t.grid <- ans.calibration$t.grid
  Phi <- ans.calibration$Phi
  P <- ans.calibration$P
  max.ord <- ans.calibration$max.ord
  include2 <- ans.calibration$include2
  include3 <- ans.calibration$include3
  include4 <- ans.calibration$include4
  cat.theta <- ans.calibration$cat.theta
  inv.transform <- ans.calibration$inv.transform
  pow <- ans.calibration$pow

  I <- nrow(x.new)
  M <- dim(alpha)[3]
  J <- ncol(x.new)
  K <- ncol(t.new)
  n.MCMC <- dim(alpha)[1]
  Ups.hat <- apply(ans.calibration$Ups[floor(.5*n.MCMC):n.MCMC,,],c(2,3),mean)
  if(is.null(post.ind))
    post.ind <- sample(floor(.5*n.MCMC):n.MCMC, Nit)

 ## Transform x.new to [0,1]
  X <- x.new
  for(j in 1:J)
    X[,j] <- (X[,j]-X.lim[j,1])/(X.lim[j,2]-X.lim[j,1])
  
 ## Transform t.new to [0,1]
  T <- t.new
  for(k in 1:K)
    T[,k] <- (T[,k]-T.lim[k,1])/(T.lim[k,2]-T.lim[k,1])
  
 ## get design mat for emulator
  ans <- get.all.X(cbind(X,T), t.grid, Phi, P=P, max.ord=max.ord, include2, include3, include4, cat.vars=c(rep(0,ncol(X)), cat.theta))
  Xem <- ans$X
  col.ind.em <- ans$col.ind
  yem.post <- array(0,c(length(post.ind),I,M))
  for(i in 1:length(post.ind)){
    yem.post[i,,] <- Xem%*%alpha[post.ind[i],,]
  }
  if(trans.scale){
    for(m in 1:M){
      yem.post[,,m] <- inv.transform[[m]](yem.post[,,m], pow[m])
    }
  }
  return(list(y.sim=yem.post, Ups.hat=Ups.hat))
}






predict.at.mean <- function(x.new, t.new=NULL, ans.calibration, post.ind=NULL, trans.scale=TRUE){

  X.lim <- ans.calibration$X.lim
  theta <- ans.calibration$theta
  alpha <- ans.calibration$alpha
  beta <- ans.calibration$beta
  t.grid <- ans.calibration$t.grid
  Phi <- ans.calibration$Phi
  P <- ans.calibration$P
  max.ord <- ans.calibration$max.ord
  include2 <- ans.calibration$include2
  include3 <- ans.calibration$include3
  include4 <- ans.calibration$include4
  cat.theta <- ans.calibration$cat.theta
  col.ind.dis <- ans.calibration$col.ind.dis
  inv.transform <- ans.calibration$inv.transform

  I <- nrow(x.new)
  M <- dim(beta)[3]
  J <- ncol(x.new)
  K <- ncol(theta)
  n.MCMC <- dim(beta)[1]
  if(is.null(post.ind))
    post.ind <- floor(.5*n.MCMC):n.MCMC
  Sigma.hat <- apply(ans.calibration$Sigma[post.ind,,],c(2,3),mean)
  Ups.hat <- apply(ans.calibration$Ups[post.ind,,],c(2,3),mean)
  Lambda.hat <- apply(ans.calibration$Lambda[post.ind,,,],c(2,3,4),mean)
  Gamma.hat <- apply(ans.calibration$Gamma[post.ind,,,],c(2,3,4),mean)

 ## Transform x.new to [0,1]
  X <- x.new
  for(j in 1:J)
    X[,j] <- (X[,j]-X.lim[j,1])/(X.lim[j,2]-X.lim[j,1])

 ## Transform t.new to [0,1]
  if(!is.null(t.new)){
    T <- t.new
    for(k in 1:K)
      T[,k] <- (T[,k]-T.lim[k,1])/(T.lim[k,2]-T.lim[k,1])
  }
 ## posterior mean of theta
  theta.hat <- apply(theta[post.ind,], 2, mean)
  for(k in 1:K){
    if(cat.theta[k]!=0){
      blah <- c(theta[post.ind,k],1:cat.theta[k])
      theta.hat[k] <- order(-table(blah))[1]
    }
  }

 ## posterior mean alpha, beta
  alpha.hat <- apply(alpha[post.ind,,], c(2,3), mean)
  beta.hat <- apply(beta[post.ind,,], c(2,3), mean)

 ## get design mat for discrepancy
  ans <- get.all.X(X, t.grid, Phi, P=P, max.ord=max.ord, include2, include3, include4, cat.vars=rep(0,ncol(X)))
  X.dis <- ans$X

 ## get design mat for emulator
  if(is.null(t.new))
    ans <- get.all.X(cbind(X,matrix(theta.hat,nrow=nrow(X),ncol=K,byrow=TRUE)), t.grid, Phi, P=P, max.ord=max.ord, include2, include3, include4, cat.vars=c(rep(0,ncol(X)), cat.theta))
  else
    ans <- get.all.X(cbind(X,T), t.grid, Phi, P=P, max.ord=max.ord, include2, include3, include4, cat.vars=c(rep(0,ncol(X)), cat.theta))
    Xem <- ans$X
  

  d.hat <- X.dis%*%beta.hat
  yem.hat <- Xem%*%alpha.hat
  y.hat <- yem.hat +  d.hat
  if(trans.scale){
    for(m in 1:M){
      y.hat[,m] <- inv.transform[[m]](y.hat[,m], pow[m])
      yem.hat[,m] <- inv.transform[[m]](yem.hat[,m], pow[m])
      d.hat[,m] <- y.hat[,m] - yem.hat[,m]
    }
  }
  return(list(y.hat=y.hat, y.sim=yem.hat, d.hat=d.hat, theta.hat=theta.hat, alpha.hat=alpha.hat, beta.hat=beta.hat, Sigma.hat=Sigma.hat, Ups.hat=Ups.hat, Lambda.hat=Lambda.hat, Gamma.hat=Gamma.hat))
}









########################################################################
############## Bernoulli MCMC Example for SHort Course ##################
########################################################################


Bern.like <- function(y, theta){
  return(sum(dbinom(y,1,theta,log=TRUE)))
}


pi.theta.bern <- function(theta, a=.5, b=.5){
  return(dbeta(theta,a,b,log=TRUE))
}


Bern.MCMC <- function(y,a=.5,b=.5, n.MCMC=10000, prop.sd=.1, nplot=1000, nback=2000, delay=0, change=n.MCMC, tlim=NULL, plot.prop=TRUE){

  a.star <- a+sum(y)
  b.star <- b+n-sum(y)
  t.plot <- seq(0,1,length=500)
  d.plot <- dbeta(y.plot,a.star,b.star)
  tlim.h <- c(min(t.plot[d.plot>1E-3]),max(t.plot[d.plot>1E-3]))

  theta.now <- .5
  theta.r <- -1
  theta <- rep(.5,n.MCMC)
  theta.red <- rep(-1,n.MCMC)
  accept.count <- rep(0,n.MCMC)

  if(length(delay))
    delay <- rep(delay,n.MCMC)

  for(it in 1:n.MCMC){

Sys.sleep(delay[it])

    theta.p <- theta.now + rnorm(1,0,prop.sd)
    if(theta.p>=0 && theta.p<=1){
      pi.prop <- pi.theta.bern(theta.p)
      pi.now <- pi.theta.bern(theta.now)
      like.prop <- Bern.like(y,theta.p)
      like.now <- Bern.like(y,theta.now)
      MH.ratio <- exp(like.prop+pi.prop-like.now-pi.now)
    }
    else
      MH.ratio <- 0
    if(!is.na(MH.ratio) && runif(1) < MH.ratio){
      theta.now <- theta.p
      accept.count[it] <- 1
      theta.r <- -1
    }
    else
      theta.r <- theta.p

    theta.red[it] <- theta.r
    theta[it] <- theta.now

    if(it>1 && it%%nplot==0){
      cat("\niter =",it,"out of",n.MCMC)
      cat("\naccept.pct = ", mean(accept.count[1:it]))
      ind.now <- max(floor(it/2),it-nback+1):it
      par(mfrow=c(1,2))
#      plot(ind.now,theta[ind.now],xlab="MCMC iteration",ylab="Theta")
      plot(1:it,theta[1:it],xlab="MCMC iteration",ylab="Theta",ylim=tlim)
      if(plot.prop)
        points(theta.red[1:it],col=2)
      hist(theta[ind.now],prob=TRUE,main="Theta Posterior",xlim=tlim.h, breaks=seq(0,1,length=100),ylim=c(0,max(d.plot)*1.1),xlab="Theta")
      lines(t.plot,d.plot,col=4)

    }
  }
  return(theta)
}



######################################################################
############## Normal MCMC Example for SHort Course ##################
######################################################################


norm.like <- function(y, mu, sigma){
  return(sum(dnorm(y,mu,sigma,log=TRUE)))
}


pi.mu <- function(mu, mu0=0, sigma0=10){
  return(dnorm(mu,mu0,sigma0,log=TRUE))
}


norm.MCMC <- function(y,sigma=1, mu0=0,sigma0=100, n.MCMC=10000, prop.sd=.1, nplot=1000, nback=2000, delay=0, change=n.MCMC, tlim=NULL, plot.prop=TRUE){

  n <- length(y)
  sigma.star <- sqrt(1/(1/sigma0^2+n/sigma^2))
  mu.star <- (mu0/sigma0^2+sum(y)/sigma^2)*sigma.star^2
  mu.plot <- seq(-2,2,length=500)
  d.plot <- dnorm(mu.plot,mu.star,sigma.star)
  tlim.h <- c(min(mu.plot[d.plot>1E-3]),max(mu.plot[d.plot>1E-3]))

  mu.now <- 0
  mu.r <- -1000
  mu <- rep(0,n.MCMC)
  mu.red <- rep(-1000,n.MCMC)
  accept.count <- rep(0,n.MCMC)

  if(length(delay))
    delay <- rep(delay,n.MCMC)

  for(it in 1:n.MCMC){

Sys.sleep(delay[it])

    mu.p <- mu.now + rnorm(1,0,prop.sd)
    pi.prop <- pi.mu(mu.p)
    pi.now <- pi.mu(mu.now)
    like.prop <- norm.like(y,mu.p,sigma)
    like.now <- norm.like(y,mu.now,sigma)
    MH.ratio <- exp(like.prop+pi.prop-like.now-pi.now)
    if(!is.na(MH.ratio) && runif(1) < MH.ratio){
      mu.now <- mu.p
      accept.count[it] <- 1
      mu.r <- -1000
    }
    else
      mu.r <- mu.p

    mu.red[it] <- mu.r
    mu[it] <- mu.now

    if(it>0 && it%%nplot==0){
      cat("\niter =",it,"out of",n.MCMC)
      cat("\naccept.pct = ", mean(accept.count[1:it]))
      ind.now <- max(floor(it/2),it-nback+1):it
      par(mfrow=c(1,2))
#      plot(ind.now,mu[ind.now],xlab="MCMC iteration",ylab="Mu")
      plot(0:it,c(0,mu[1:it]),xlab="MCMC iteration",ylab="Mu",ylim=tlim)
      if(plot.prop)
        points(1:it,mu.red[1:it],col=2)
      hist(mu[ind.now],prob=TRUE,main="Mu Posterior",xlim=tlim.h, breaks=seq(-2,2,length=100),ylim=c(0,max(d.plot)*1.1),xlab="Mu")
      lines(mu.plot,d.plot,col=4)

    }
  }
  return(mu)
}


bootstrap.norm <- function(y,sigma,nboot=1000, conf=.95){

  n <- length(y)
  ybar <- mean(y)
  probs <- c((1-conf)/2, 1-(1-conf)/2)
  z.np <- z.p <- rep(0,nboot)
  for(i in 1:nboot){
    ynp <- sample(y,n,replace=TRUE)
    yp <- rnorm(n,ybar,sigma)
    z.np[i] <- mean(ynp)-ybar
    z.p[i] <- mean(yp)-ybar
  }
  xi.np <- quantile(z.np,prob=probs)
  xi.p <- quantile(z.p,prob=probs)
  CI.np <- c(ybar-xi.np[2],ybar-xi.np[1])
  CI.p <- c(ybar-xi.p[2],ybar-xi.p[1])
  return(list(CI.np=CI.np, CI.p=CI.p))

}




########################################################################
############## SA across Theta for specified input points ##############
########################################################################


get.SA.one.x <- function(ans.calibration, x.sa, n.mc.T=1000, n.post=100, n.burn=NULL, m=1){

  n.mcmc <- dim(ans.calibration$alpha)[1]
  if(is.null(n.burn))
    n.burn <- floor(n.mcmc/2)
  K <- ncol(ans.calibration$theta)
  J <- ncol(ans.calibration$X)

  post.ind <- sample((n.burn+1):n.mcmc, n.post, replace=TRUE)
  T.mat <- matrix(0,K,n.post)
  for(it in 1:n.post){
print(it)

#x.sa..<<-x.sa
#J..<<-J
#K..<<-K
#post.ind..<<-post.ind

    f.it <- function(X){
      ans <- predict.sim(matrix(x.sa,nrow(X),J,byrow=T), t.new=X, ans.calibration, post.ind=post.ind[it], trans.scale=TRUE)$y.sim
      return(ans[1,,m])
    }
    blah <- get.total.var(f.it,nx=K,N=n.mc.T, navg=1, distn="emp",X=Ts)$T.hat
    T.mat[,it] <- blah
  }
  return(T.mat)
}



get.SA <- function(ans.calibration, x.sa=NULL, n.mc.T=1000, n.post=100, n.burn=NULL, m=1){

  X.lim <- ans.calibration$X.lim
  if(is.null(x.sa)){
    x.sa <- unique(ans.calibration$Xs)
  }
  else{
    for(j in 1:ncol(x.sa))
      x.sa[,j] <- (x.sa[,j]-X.lim[j,1])/(X.lim[j,2]-X.lim[j,1])
  }
  J <- ncol(x.sa)
  if(J==1)
    x.sa <- x.sa[order(x.sa[,1]),]
  if(J==2)
    x.sa <- x.sa[order(x.sa[,1],x.sa[,2]),]
  if(J>=3)
    x.sa <- x.sa[order(x.sa[,1],x.sa[,2],x.sa[,3]),]

  x.sa2 <- x.sa
  for(j in 1:ncol(x.sa))
    x.sa2[,j] <- x.sa[,j]*(X.lim[j,2]-X.lim[j,1])+X.lim[j,1]
  

  nx <- nrow(x.sa)
  K <- ncol(ans.calibration$theta)
  T.mat <- matrix(0,nx,K)
  for(i in 1:nx){
cat("\ni =",i,"out of",nx,"\n")
cat("x =",x.sa[i,],"\n")
    ans.i <- get.SA.one.x(ans.calibration, x.sa[i,], n.mc.T=n.mc.T, n.post=n.post, n.burn=n.burn, m=1)
    T.mat[i,] <- rowMeans(ans.i)
  }
  return(list(T.mat=T.mat, x.sa=x.sa2))
}





  





















##################################################################################
################################## OLD CODE ######################################
##################################################################################





plot.fitted.old <- function(ans.calibration, post.ind=NULL, ncurves=25,transpose=FALSE,y.names=NULL,x.names=NULL, nx.curve=50){

  sim.col <- 5
  dis.col <- 2
  y <- ans.calibration$y
  X <- ans.calibration$X
  X.dis <- ans.calibration$X.dis
  X.lim <- ans.calibration$X.lim
  theta <- ans.calibration$theta
  alpha <- ans.calibration$alpha
  beta <- ans.calibration$beta
  t.grid <- ans.calibration$t.grid
  Phi <- ans.calibration$Phi
  P <- ans.calibration$P
  max.ord <- ans.calibration$max.ord
  include2 <- ans.calibration$include2
  include3 <- ans.calibration$include3
  include4 <- ans.calibration$include4
  cat.theta <- ans.calibration$cat.theta
  col.ind.em <- ans.calibration$col.ind.em
  col.ind.dis <- ans.calibration$col.ind.dis
  tpow <- ans.calibration$tpow
  res.obs <- ans.calibration$res.obs
  Sigma <- ans.calibration$Sigma

  ind.na.obs <- is.na(y)
  I <- nrow(y)
  M <- ncol(y)
  J <- ncol(X)
  K <- ncol(theta)
  rows <- M
  cols <- J
  n.MCMC <- nrow(theta)
  if(is.null(post.ind))
    post.ind <- floor(.5*n.MCMC):n.MCMC
  ind.curves <- sample(post.ind,ncurves)
  if(is.null(y.names))
    y.names <- paste("Y.",1:M)
  if(is.null(x.names))
    x.names <- paste("X.",1:J)

  if(transpose)
    par(mfcol=c(cols,rows))
  else
    par(mfrow=c(rows,cols))

  ## get obs error bounds
#   res.conf <-  apply(res.obs[post.ind,,], c(2,3), quantile, probs=c(.01, .99))
#   y.conf <- array(0,c(2,I,M))
#   y.conf[1,,] <- y - res.conf[2,,]
#   y.conf[2,,] <- y - res.conf[1,,]
  Sigma.hat <- apply(Sigma[post.ind,,], c(2,3), mean)
  y.conf <- array(0,c(2,I,M))
  for(m in 1:M)
    y.conf[,,m] <- rbind(y[,m]-2*sqrt(Sigma.hat[m,m]), y[,m]+2*sqrt(Sigma.hat[m,m]))^(1/tpow[m])


  for(m in 1:M){
    ind.na.m <- ind.na.obs[,m]
    for(j in 1:J){
      x.j <- X[,j]*(X.lim[j,2]-X.lim[j,1])+X.lim[j,1]
      yd <- y[!ind.na.m,m]^(1/tpow[m])
      R.yd <- diff(range(y.conf[,,m],na.rm=T))
      yd.lim <- c(min(y.conf[,,m], na.rm=T)-R.yd*.1, max(y.conf[,,m],na.rm=T)+R.yd*.1)
      plot(x.j[!ind.na.m], yd, main="Simulator+Discrepancy", ylim=yd.lim, pch='x', col=0, ylab=y.names[m], xlab=x.names[j])  
      legend("topright", legend=c("simulator", "sim + discrep"), lty=c(1,1), col=c(sim.col,dis.col))
      ord.j <- order(x.j)
      X.plot <- matrix(seq(0,1, length=nx.curve), nrow=nx.curve, ncol=J)
      ans <- get.all.X(X.plot, t.grid, Phi, P=P, max.ord=max.ord, include2, include3, include4, cat.vars=rep(0,ncol(X.plot)))
      X.dis <- ans$X
      xj.plot <- X.plot[,j]*(X.lim[j,2]-X.lim[j,1])+X.lim[j,1]
      for(i in ind.curves){
        ans <- get.all.X(cbind(X.plot,matrix(theta[i,],nrow=nrow(X.plot),ncol=K,byrow=TRUE)), t.grid, Phi, P=P, max.ord=max.ord, include2, include3, include4, cat.vars=c(rep(0,ncol(X.plot)), cat.theta))
        Xem.i <- ans$X
        col.ind <- ans$col.ind

        y.sim <- Xem.i[,c(1,col.ind[[j]])]%*%alpha[i,c(1,col.ind[[j]]),m]
        y.hat <- y.sim + X.dis[,c(1,col.ind[[j]])]%*%beta[i,c(1,col.ind[[j]]),m]
        lines(xj.plot,y.hat^(1/tpow[m]),col=dis.col,lwd=.5)
        lines(xj.plot,y.sim^(1/tpow[m]),col=sim.col,lwd=.5)
      }
      points(x.j[!ind.na.m], yd, pch='x', col=1, cex=1.5)  
      for(i in (1:I)[!ind.na.m]){
        lines(rep(x.j[i],2), y.conf[,i,m],col=1,lwd=2)
      }
    }
  }
}






plot.discrepancy.old <- function(ans.calibration, post.ind=NULL, ncurves=25,transpose=FALSE,y.names=NULL,x.names=NULL){

  mu.col <- 4
  cb.col <- 2
  dis.col <- grey(.4)
  y <- ans.calibration$y
  X <- ans.calibration$X
  X.dis <- ans.calibration$X.dis
  X.lim <- ans.calibration$X.lim
  beta <- ans.calibration$beta
  col.ind.dis <- ans.calibration$col.ind.dis
  tpow <- ans.calibration$tpow
  res.obs <- ans.calibration$res.obs
  res.obs.sim <- ans.calibration$res.obs.sim

  I <- nrow(y)
  M <- ncol(y)
  J <- ncol(X)
  rows <- M
  cols <- J
  n.MCMC <- dim(beta)[1]
  if(is.null(post.ind))
    post.ind <- floor(.5*n.MCMC):n.MCMC
  ind.curves <- sample(1:length(post.ind),ncurves)
  if(is.null(y.names))
    y.names <- paste("Y.",1:M)
  if(is.null(x.names))
    x.names <- paste("X.",1:J)

  if(transpose)
    par(mfcol=c(cols,rows))
  else
    par(mfrow=c(rows,cols))

  discrep <- res.obs.sim[ind.curves,,] - res.obs[ind.curves,,]
  for(m in 1:M){
    for(j in 1:J){
      x.j <- X[,j]*(X.lim[j,2]-X.lim[j,1])+X.lim[j,1]
      djm.post <- t(X.dis[,c(1,col.ind.dis[[j]])]%*%t(beta[post.ind,c(1,col.ind.dis[[j]]),m]))

      ans.bands <- get.cred.bands(djm.post, alpha=.05)
      djm.conf <- ans.bands$band
      djm.mu <- ans.bands$mu
#      djm.conf <- apply(djm.post,2,quantile,probs=c(.025, .975))
#      djm.mu <- apply(djm.post,2,mean)

      plot(x.j, rep(0,length(x.j)), main=paste("Delta.",j,".",m," by x.",j,sep=""), ylim=range(djm.post), col=0, ylab=y.names[m], xlab=x.names[j])  
      legend("topright", legend=c("Posterior Mean", "95% Bands","Zero Line"), lty=c(1,2,2), col=c(mu.col,cb.col,1))
      ord.j <- order(x.j)
      for(i in ind.curves){
        lines(x.j[ord.j],djm.post[i,ord.j],col=dis.col,lwd=.5)
      }
      lines(x.j[ord.j], djm.mu[ord.j], col=mu.col)
      lines(x.j[ord.j], djm.conf[1,ord.j], col=cb.col)
      lines(x.j[ord.j], djm.conf[2,ord.j], col=cb.col)
      lines(x.j[ord.j], rep(0,length(x.j)), col=1, lty=2)
    }
  }
}












