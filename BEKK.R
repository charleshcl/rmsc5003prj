#install.packages("MTS")
#install.packages("expm")
#install.packages("MASS")
library(MTS)
library(expm)
library(MASS)
library(zoo)
# Fit BEKK model (bivariate GARCH Model) use log return
#M2 <- BEKK11(R) # R is the 757 * 2 log-returns from independent GARCH

# Extract coefficients/parameters in BEKK model and store them in vectors/matrices
mu <- M2$estimates[1:2]
C <- matrix(c(M2$estimates[3:4],0,M2$estimates[5]), byrow=F, nrow=2)
A <- matrix(c(M2$estimates[6:9]), byrow=F, nrow=2)
B <- matrix(c(M2$estimates[10:13]), byrow=F, nrow=2)

# find the last day Residual vector and SIGMA matrix for simulation
Last.Residual <- R[length(R[,1]),] - mu
Last.SIGMA <- matrix(M2$Sigma.t[length(R[,1]),], nrow=2)


# Simulation based on BEKK Model
update_SIGMA <- function(SIGMA, Xi){
  update_SIGMA <- C %*% t(C) + A %*% (Xi %*% t(Xi)) %*% t(A) + B %*% SIGMA %*% t(B)
}

BEKK_Simu <- function(rept){
  Return_base <- matrix(rep((mu), T+1), nrow=2, byrow=F) # Constant mean vector
  SIGMA <- matrix(rep(Last.SIGMA,T+1), nrow=2, byrow=F) # For storing the conditional Variance Covariance Matrices
  Residual <- matrix(rep(Last.Residual,T+1), nrow=2, byrow=F) # For storing the residuals
  for (i in 2:(T+1)){
    SIGMA[,(2*i-1):(2*i)] <- update_SIGMA(SIGMA[,(2*i-3):(2*i-2)], Residual[,i-1])
    Residual[,i] <- mvrnorm(n = 1, mu = c(0,0), Sigma = SIGMA[,(2*i-1):(2*i)])
  }
  Sim_R <- Return_base + Residual
  Sim_R[,1] <- c(0,0)
  Sim_S <- exp(t(apply(Sim_R, 1, cumsum))) * c(AAPL[length(AAPL)], AMZN[length(AMZN)])
  return(Sim_S)
}
Output_BEKK_AAPL <- matrix(0, nrow = P, ncol=253)
Output_BEKK_AMZN <- Output_BEKK_AAPL
for (rept in 1:P){
  temp <- BEKK_Simu(rept)
  Output_BEKK_AAPL[rept,] <- temp[1,]
  Output_BEKK_AMZN[rept,] <- temp[2,]
}

plot.zoo(t(Output_BEKK_AAPL),plot.type="single")
plot.zoo(t(Output_BEKK_AMZN),plot.type="single")
