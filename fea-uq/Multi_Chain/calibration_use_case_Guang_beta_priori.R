source("BSSANOVA_calibration_Multi.R")
### Read in Data

sim.dat <- read.table("sim_Guang.csv", sep=",",header=FALSE)
exp.dat <- read.table("exp_Guang.csv", sep=",",header=FALSE)
### delete ID 13 in exp.data because of the N/A

head(exp.dat)
head(sim.dat)


### Assemble Data objects for input into calibration function

y <- exp.dat[,2:5]    ### kg mass transfer coefficient
ys <- sim.dat[,5:8]   ### kg mass transfer coefficient
X <- exp.dat[,1]
Xs <- sim.dat[,4]
Ts <- sim.dat[,1:3]

X.lim <- rbind(c(0.0,0.5))
T.lim <- rbind(c(-0.1,1.6), c(0,5.1), c(0,6.2))
cat.theta <- c(0,0,0) # all continuous density

## number of inputs variables
P <- 1

## number of model parameters
Q <- 3

## number of outputs
M <- 4

## Number of MCMC iterations
## it is set to 500 for illustration, should typically use 10000 or more
N.mcmc <- 40000
cat.locations = NULL

#transform <- c(.4, "logit")
#transform <- c(.4, 0.4, 0.4, 0.4)
transform <- c("log", "log", "log", "log")

### prior distribution for theta
### uniform distribution

pi.theta <- function(theta){
  R.T <- T.lim[,2]-T.lim[,1]
  theta.t <- (theta - T.lim[,1])/R.T
  ans <- -sum(log(R.T[1:3]))+dbeta(theta.t[1],1,1, log=TRUE) + dbeta(theta.t[2],2.4,1.0, log=TRUE) + dbeta(theta.t[3],1.5,1.0, log=TRUE) 
  return(ans)
}

#uniform priori is used!
#pi.theta <- function(theta){
#  R.T <- T.lim[,2]-T.lim[,1]
#  theta.t <- (theta - T.lim[,1])/R.T
#  ans <- -sum(log(R.T[1:3]))+dunif(theta.t[1],0,1, log=TRUE) + dunif(theta.t[2],0,1, log=TRUE) + dunif(theta.t[3],0,1, log=TRUE) 
#  return(ans)
#}


## Prior mean for the measurement error variance Sigma
#mu.Sigma <- c(0.09)
mu.Sigma <- rbind(c(.2,0,0,0),c(0,0.2,0,0),c(0,0,0.2,0),c(0,0,0,0.2));
#v.Sigma <- 50
v.Sigma <- 6

## Prior mean for the simulator output error variance Upsilon
mu.Ups <- mu.Sigma
#v.Ups <- 50
v.Ups <- 6

theta.prop.sd <- c(.2,.2,.2,.2,.2,.2,.2)

##############################
###### RUN Calibration #######sim
##############################

#set.seed(11)
set.seed(222)

sim1.dat = cbind(Xs,Ts,ys)
exp1.dat = cbind(X,y)

# we use default mu.Sigma, mu.Ups, which is comparable to the var of the simulation data
ans.calibration <- Calibrate(sim.dat = sim1.dat, exp.dat = exp1.dat, P=P, Q=Q, M=M, N.mcmc=N.mcmc, X.lim=X.lim, T.lim=T.lim, transform = transform, pi.theta = pi.theta, v.Sigma=v.Sigma, v.Ups=v.Ups, theta.prop.sd=theta.prop.sd)
#ans.calibration <- Calibrate(sim.dat = sim1.dat, exp.dat = exp1.dat, P=P, Q=Q, M=M, N.mcmc=N.mcmc, X.lim=X.lim, T.lim=T.lim, transform = transform, pi.theta = pi.theta,mu.Sigma=mu.Sigma, v.Sigma=v.Sigma, mu.Ups=mu.Ups, v.Ups=v.Ups, theta.prop.sd=theta.prop.sd)
#ans.calibration <- Calibrate(sim.dat = sim1.dat, exp.dat = exp1.dat, P=P, Q=Q, M=M, N.mcmc=N.mcmc, X.lim=X.lim, T.lim=T.lim, transform = transform, mu.Sigma=mu.Sigma, v.Sigma=v.Sigma, mu.Ups=mu.Ups, v.Ups=v.Ups, theta.prop.sd=theta.prop.sd)

## To save results use
save(ans.calibration, file="results.Rdata")

## To come back into R later and use the results, use
#load("results.Rdata")
