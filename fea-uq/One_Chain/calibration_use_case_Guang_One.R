source("BSSANOVA_calibration_One.R")
### Read in Data

sim.dat <- read.table("8-12_sim.csv", sep=",",header=TRUE)
exp.dat <- read.table("8-12_exp.csv", sep=",",header=TRUE)

head(exp.dat)
head(sim.dat)


### Assemble Data objects for input into calibration function
y <- exp.dat[,2:4]    
ys <- sim.dat[,5:7]  
X <- exp.dat[,1]
Xs <- sim.dat[,4]
Ts <- sim.dat[,2:3]   ### we have the constraint D1+D2 = 1.8, only have two calibration parameters D2, D3

X.lim <- rbind(c(0.16,0.44))
T.lim <- rbind(c(0,1.8), c(0,2))
cat.theta <- c(0,0) # all continuous density

## number of inputs variables
P <- 1

## number of model parameters
Q <- 2

## number of outputs
M <- 3

## Number of MCMC iterations
## it is set to 500 for illustration, should typically use 10000 or more
N.mcmc <- 50000
cat.locations = NULL

#transform <- c(.4, "logit")
transform <- c("log", "log","0.4","log")

### prior distribution for theta
# pi.theta <- function(theta){
  # R.T <- T.lim[,2]-T.lim[,1]
  # theta.t <- (theta - T.lim[,1])/R.T
  # ans <- -sum(log(R.T[1:2]))+dbeta(theta.t[1],1,3, log=TRUE) + dbeta(theta.t[2],3,1, log=TRUE)
  # return(ans)
# }

### uniform distribution
pi.theta <- function(theta){
 R.T <- T.lim[,2]-T.lim[,1]
 theta.t <- (theta - T.lim[,1])/R.T
 ans <- -sum(log(R.T[1:2]))+dunif(theta.t[1],0,1, log=TRUE) + dunif(theta.t[2],0,1, log=TRUE)
 return(ans)
}


## Prior mean for the measurement error variance Sigma
#mu.sigma <- rbind(c(.009,.001),c(.001,.009))
#mu.Sigma <- rbind(c(.0009,0,0,0),c(0,.00026,0,0),c(0,0,.022,0),c(0,0,0,.0004))
#mu.Sigma <- rbind(c(.09,0,0),c(0,.022,0),c(0,0,.2))
mu.Sigma <- rbind(c(.09,0,0),c(0,.005,0),c(0,0,.5))
v.Sigma <- 40

## Prior mean for the simulator output error variance Upsilon
mu.Ups <- rbind(c(.07,0,0),c(0,.07,0),c(0,0,.07))
v.Ups <- 20

theta.prop.sd <- c(.2,.2)

##############################
###### RUN Calibration #######sim
##############################

set.seed(11)

sim1.dat = cbind(Xs,Ts,ys)
exp1.dat = cbind(X,y)

ans.calibration <- Calibrate(sim.dat = sim1.dat, exp.dat = exp1.dat, P=P, Q=Q, M=M, N.mcmc=N.mcmc, X.lim=X.lim, T.lim=T.lim, transform = transform,pi.theta = pi.theta, mu.Sigma = mu.Sigma, v.Sigma=v.Sigma, mu.Ups = mu.Ups, v.Ups=v.Ups, theta.prop.sd=theta.prop.sd,nplot=100)

## To save results use
save(ans.calibration, file="results_Guang_C50000.Rdata")
## To come back into R later and use the results, use
