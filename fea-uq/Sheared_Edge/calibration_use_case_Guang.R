source("BSSANOVA_calibration_Sheared_Edge_Multi.R")
### Read in Data

sim.dat <- read.table("sim_Guang_full.csv", sep=",",header=FALSE)
exp.dat <- read.table("exp_Guang.csv", sep=",",header=FALSE)

head(exp.dat)
head(sim.dat)

### Assemble Data objects for input into calibration function

#y <- exp.dat[,2:5]    ### kg mass transfer coefficient
#yy <- exp.dat[,2:5]    ### kg mass transfer coefficient
#y <- cbind(exp.dat[,2], exp.dat[,5])    ### kg mass transfer coefficient
y <- cbind(exp.dat[,2], exp.dat[,4:5])    ### kg mass transfer coefficient
#y <- exp.dat[,5]    ### kg mass transfer coefficient
#ys <- sim.dat[,6:9]   ### kg mass transfer coefficient
ys <- cbind(sim.dat[,6], sim.dat[,8:9])   ### kg mass transfer coefficient
#ys <- cbind(sim.dat[,6], sim.dat[,9])   ### kg mass transfer coefficient
#ys <- sim.dat[,9]   ### kg mass transfer coefficient
X <- exp.dat[,1]
Xs <- sim.dat[,5]
##### Note that we only use two dimensional parameter t2, t3; t1+t2 is a constant!
Ts <- sim.dat[,3:4]

X.lim <- rbind(c(0.18,0.41))
T.lim <- rbind(c(0,1.8), c(0,2.1))
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
#transform <- c(.4, 0.4, 0.4, 0.4)
#transform <- c("log", "log", "log", "log")
transform <- c("log", "log", "log")

### prior distribution for theta
### uniform distribution

#pi.theta <- function(theta){
#  R.T <- T.lim[,2]-T.lim[,1]
#  theta.t <- (theta - T.lim[,1])/R.T
#  ans <- -sum(log(R.T[1:3]))+dbeta(theta.t[1],1,3, log=TRUE) + dbeta(theta.t[2],2.4,1.03, log=TRUE) + dbeta(theta.t[3],1.1,2.5, log=TRUE) 
#  return(ans)
#}

#uniform priori is used!
pi.theta <- function(theta){
  R.T <- T.lim[,2]-T.lim[,1]
  theta.t <- (theta - T.lim[,1])/R.T
  ans <- -sum(log(R.T[1:2]))+dunif(theta.t[1],0,1, log=TRUE) + dunif(theta.t[2],0,1, log=TRUE)  
  return(ans)
}


## Prior mean for the measurement error variance Sigma
#mu.Sigma <- c(0.09)
#mu.Sigma <- rbind(c(0.4,0),c(0,0.4));
mu.Sigma <- rbind(c(0.3,0,0),c(0,0.3,0),c(0,0,0.3));
#mu.Sigma <- 0.02;
#v.Sigma <- 50
v.Sigma <- 40

## Prior mean for the simulator output error variance Upsilon
#mu.Ups <- mu.Sigma
mu.Ups <- rbind(c(0.07,0,0),c(0,0.07,0),c(0,0,0.07))
#v.Ups <- 50
v.Ups <- 20

theta.prop.sd <- c(1.0,.7,.7,.7,.7,.7,.7)

##############################
###### RUN Calibration #######sim
##############################

#set.seed(11)
set.seed(8968)

sim1.dat = cbind(Xs,Ts,ys)

exp1.dat = cbind(X,y)
print(Xs);
print(Ts);
print(ys);
#cat("exp1.dat is ", exp1.dat);

# we use default mu.Sigma, mu.Ups, which is comparable to the var of the simulation data
#ans.calibration <- Calibrate(sim.dat = sim1.dat, exp.dat = exp1.dat, P=P, Q=Q, M=M, N.mcmc=N.mcmc, X.lim=X.lim, T.lim=T.lim,  pi.theta = pi.theta, mu.Sigma=mu.Sigma, v.Sigma=v.Sigma, mu.Ups=mu.Ups, v.Ups=v.Ups, theta.prop.sd=theta.prop.sd)
ans.calibration <- Calibrate(sim.dat = sim1.dat, exp.dat = exp1.dat, P=P, Q=Q, M=M, N.mcmc=N.mcmc, X.lim=X.lim, T.lim=T.lim, transform = transform, pi.theta = pi.theta, mu.Sigma=mu.Sigma, v.Sigma=v.Sigma, mu.Ups=mu.Ups, v.Ups=v.Ups, theta.prop.sd=theta.prop.sd)
#ans.calibration <- Calibrate(sim.dat = sim1.dat, exp.dat = exp1.dat, P=P, Q=Q, M=M, N.mcmc=N.mcmc, X.lim=X.lim, T.lim=T.lim, transform = transform, pi.theta = pi.theta,mu.Sigma=mu.Sigma, v.Sigma=v.Sigma, mu.Ups=mu.Ups, v.Ups=v.Ups, theta.prop.sd=theta.prop.sd)
#ans.calibration <- Calibrate(sim.dat = sim1.dat, exp.dat = exp1.dat, P=P, Q=Q, M=M, N.mcmc=N.mcmc, X.lim=X.lim, T.lim=T.lim, transform = transform, mu.Sigma=mu.Sigma, v.Sigma=v.Sigma, mu.Ups=mu.Ups, v.Ups=v.Ups, theta.prop.sd=theta.prop.sd)

## To save results use
save(ans.calibration, file="results.Rdata")

## To come back into R later and use the results, use
#load("results.Rdata")
