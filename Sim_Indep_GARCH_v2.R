#install.packages("quantmod")
#install.packages("fGarch")
library(quantmod)
library(fGarch)
############# Parameters for Adjustment ##############
P = 100000 # Number of Sample Path
T = 252  # Length of Simulated Series

############# Download stock data and calculate log-returns ##############
AAPL <- as.numeric(getSymbols("AAPL", from ="2015-11-10", to = "2018-11-12", auto.assign = FALSE)$AAPL.Close)
AMZN <- as.numeric(getSymbols("AMZN", from ="2015-11-10", to = "2018-11-12", auto.assign = FALSE)$AMZN.Close)
Data <- cbind(AAPL, AMZN)
R <- apply(log(Data),2,diff)
AAPL_Return <- R[,1]
AMZN_Return <- R[,2]
# Plot returns
par(mfrow=c(2,1))
plot(AAPL_Return ,type="l",xlab="Time",ylab = "AAPL", main="Log-returns for AAPL") # Plot log-return series to visualize heteroskedsticity
plot(AMZN_Return,type="l",xlab="Time",ylab = "AMZN",  main="Log-returns for AMZN")

# Plot acf & pacf
acf(AAPL_Return)
acf(AMZN_Return)
pacf(AAPL_Return)
pacf(AMZN_Return)
acf(AAPL_Return^2)
acf(AMZN_Return^2)
pacf(AAPL_Return^2)
pacf(AMZN_Return^2)

NormModel_1 <- vector("list", 6)
NormModel_2 <- vector("list", 6)

############# Fit GARCH Model with normal innovations ###################

NormModel_1[[1]] <- garchFit(formula ~ garch(1,0),data = AAPL_Return) 
NormModel_1[[2]] <- garchFit(formula ~ garch(1,1),data = AAPL_Return) 
NormModel_1[[3]] <- garchFit(formula ~ garch(1,2),data = AAPL_Return) 
NormModel_1[[4]] <- garchFit(formula ~ garch(2,0),data = AAPL_Return) 
NormModel_1[[5]] <- garchFit(formula ~ garch(2,1),data = AAPL_Return) 
NormModel_1[[6]] <- garchFit(formula ~ garch(2,2),data = AAPL_Return) 
NormModel_2[[1]] <- garchFit(formula ~ garch(1,0),data = AMZN_Return) 
NormModel_2[[2]] <- garchFit(formula ~ garch(1,1),data = AMZN_Return)
NormModel_2[[3]] <- garchFit(formula ~ garch(1,2),data = AMZN_Return)
NormModel_2[[4]] <- garchFit(formula ~ garch(2,0),data = AMZN_Return)
NormModel_2[[5]] <- garchFit(formula ~ garch(2,1),data = AMZN_Return)
NormModel_2[[6]] <- garchFit(formula ~ garch(2,2),data = AMZN_Return)
NormModel_2[[7]] <- garchFit(formula ~ garch(2,3),data = AMZN_Return)
NormModel_2[[8]] <- garchFit(formula ~ garch(3,2),data = AMZN_Return)

BIC1 <- as.numeric(sapply(1:6, function(i){NormModel_1[[i]]@fit$ics[2]}))
BIC2 <- as.numeric(sapply(1:8, function(i){NormModel_2[[i]]@fit$ics[2]}))
#tModel_1 <- garchFit(formula ~ garch(1,1),data = R[,1],cond.dist = "std",include.shape = F, shape  = 4)
#tModel_2 <- garchFit(formula ~ garch(1,1),data = R[,2],cond.dist = "std",include.shape = F, shape = 4)
NormModel_AAPL <- NormModel_1[[which.min(BIC1)]]
NormModel_AMZN <- NormModel_2[[2]]
mu_AAPL <- NormModel_AAPL@fit$coef[1]
Gcoef_AAPL <- NormModel_AAPL@fit$coef[-1]
mu_AMZN <- NormModel_AMZN@fit$coef[1]
Gcoef_AMZN <- NormModel_AMZN@fit$coef[-1]
########  Simulation Based on These Model #######

###  AAPL

Sim_AAPL <- function(rept){
  #set.seed(rept)
    Sim_AAPL_return_base <- as.numeric(rep(mu_AAPL, T+1))
    Sim_AAPL_residual <- rep(residuals(NormModel_AAPL)[length(AAPL_Return)], T+1)
    Sim_AAPL_sigma2 <- rep(volatility(NormModel_AAPL,type = "h")[length(AAPL_Return)], T+1)
    for (i in 2:(T+1)){
      #message("i=",i)
      temp <- c(1,Sim_AAPL_residual[i-1]^2,Sim_AAPL_sigma2[i-1])
      Sim_AAPL_sigma2[i] <- t(temp) %*% Gcoef_AAPL
      Sim_AAPL_residual[i] <- sqrt(Sim_AAPL_sigma2[i]) * rnorm(1)
    }
    Sim_AAPL_return <- Sim_AAPL_return_base + Sim_AAPL_residual
    Sim_AAPL_return[1] <- 0 ## This is artifically designed to make Sim_AAPL = the current price of AAPL
    Sim_AAPL <- as.numeric(exp(cumsum(Sim_AAPL_return)) * AAPL[length(AAPL)])
}
par(mfrow=c(1,1))
Output_AAPL <- t(sapply(1:P, Sim_AAPL))
plot.zoo(t(Output_AAPL),plot.type="single")

#### AMZN


Sim_AMZN <- function(rept){
  #set.seed(rept)
  Sim_AMZN_return <- as.numeric(rep(mu_AMZN, T+1))
  Sim_AMZN_residual <- rep(residuals(NormModel_AMZN)[length(AMZN_Return)], T+1)
  Sim_AMZN_residual[1] <- residuals(NormModel_AMZN)[length(AMZN_Return)-1]
  Sim_AMZN_sigma2 <- rep(volatility(NormModel_AMZN,type = "h")[length(AMZN_Return)], T+1)
  Sim_AMZN_sigma2[1] <- volatility(NormModel_AMZN,type = "h")[length(AMZN_Return)-1]
  for (i in 2:(T+1)){
    temp <- c(1,Sim_AMZN_residual[i-1]^2,Sim_AMZN_sigma2[i-1])
    Sim_AMZN_sigma2[i] <- t(temp) %*% Gcoef_AMZN
    Sim_AMZN_residual[i] <- sqrt(Sim_AMZN_sigma2[i]) * rnorm(1)
  }
  Sim_AMZN_return <- (Sim_AMZN_return + Sim_AMZN_residual)
  Sim_AMZN_return[1] <- 0
  Sim_AMZN <- as.numeric(exp(cumsum(Sim_AMZN_return)) * AMZN[length(AMZN)])
  Sim_AMZN
}
par(mfrow=c(1,1))

Output_AMZN <- t(sapply((1:P), Sim_AMZN))
plot.zoo(t(Output_AMZN),plot.type="single")
