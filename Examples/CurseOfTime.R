# CURSE OF TIME
# The curse of time appears also in dimension 1 meaning that lots of weight is only on one particle
# and the variance of the estimate grows exponentially with time! Particle filters do not show
# the same behaviour.

rm(list=ls())
require(SMC5)
require(parallel)
require(ggplot2)

######################
# SETTING EXPERIMENT #
######################

set.seed(1717)
n <- 50
ncores <- detectCores()
repetitions <- 10000
N <- 50
type_statistic <- "sum_squared"
comp_statistic <- NA #Used for type_statistic mean or mean_squared

##################
# GENERATE MODEL #
##################
dimension <- 1
A <- matrix(1) #If < 1 they go towards 0 -> bad approximation of relative variance
X0 <- 10 #Just so that it is far from 0 and we don't have problem when estimate relative variance
sigmaX <- matrix(1, nrow=1, ncol=1)
sigmaY <- matrix(1, nrow=1, ncol=1) #With a lower variance it should be easier for the algorithm to approximate the statistic
dataset <- generateData(n, A, X0, sigmaX, sigmaY)
X_data <- matrix(dataset$X_data, ncol=dimension)
fParams <- list(A=A, X0=X0, sigmaX=sigmaX)

##########################
# COMPUTE TRUE STATISTIC #
##########################
trueStatistic <- rep(NA, n)
for(time in 1:n) {
  trueStatistic[time] <- computeTrueStatisticFilter(X_data = X_data, n = time, type =  type_statistic,
                                                    comp = comp_statistic)
}

################
# RUN SIS/SIRS #
################
expRes <- mclapply(1:repetitions, function(rep) {
  #Generate Y always with the same X
  Y_data <- generateY(X_data, sigmaY)
  gParams <- list(Y=Y_data, sigmaY=sigmaY)
  return(list(SISRes = sequentialImportanceSampling(N=N, n=n, fParams = fParams,
                                                    gParams = gParams),
              SIRRes = bootstrapParticleFilter(N=N, n=n, fParams = fParams,
                                               gParams = gParams)))
}, mc.cores = detectCores())

#########################
# APPROXIMATE STATISTIC #
#########################

approxStatisticSIS <- matrix(NA, nrow=n, ncol=repetitions)
approxStatisticSIR <- matrix(NA, nrow=n, ncol=repetitions)
for(time in 1:n) {
  approxStatisticSIS[time,] <- sapply(expRes, function(res) {
    return(computeApproxStatisticFilter(particles = res$SISRes$filteringParticle,
                                        logWeights = res$SISRes$filteringLogWeights,
                                        n = time,
                                        type = type_statistic,
                                        comp = comp_statistic))

  })
  approxStatisticSIR[time,] <- sapply(expRes, function(res) {
    return(computeApproxStatisticFilter(particles = res$SIRRes$filteringParticle,
                                        logWeights = res$SIRRes$filteringLogWeights,
                                        n = time,
                                        type = type_statistic,
                                        comp = comp_statistic))

  })
}

########################
# COMPUTE BIAS/MSE/VAR #
########################

dfRes <- data.frame(Algorithm = character(),
                    Time = double(),
                    Type = character(),
                    Value = double(),
                    stringsAsFactors = FALSE)
for(time in 1:n) {
  estSIS <- list(Bias = mean(approxStatisticSIS[time,]-trueStatistic[time]),
              Var = mean((approxStatisticSIS[time,]-mean(approxStatisticSIS[time,]))^2),
              CoefVar = sd(approxStatisticSIS[time,])/mean(approxStatisticSIS[time,]),
              RelVar = var(approxStatisticSIS[time,])/abs(mean(approxStatisticSIS[time,])), #http://www.statisticshowto.com/relative-variance/
              RelVar2 = (sd(approxStatisticSIS[time,])/mean(approxStatisticSIS[time,]))^2,
              MSE = mean((approxStatisticSIS[time,]-trueStatistic[time])^2),
              RMSE = sqrt(mean((approxStatisticSIS[time,]-trueStatistic[time])^2)))
  estSIR <- list(Bias = mean(approxStatisticSIR[time,]-trueStatistic[time]),
                 Var = mean((approxStatisticSIR[time,]-mean(approxStatisticSIR[time,]))^2),
                 CoefVar = sd(approxStatisticSIR[time,])/mean(approxStatisticSIR[time,]),
                 RelVar = var(approxStatisticSIR[time,])/abs(mean(approxStatisticSIR[time,])),
                 RelVar2 = (sd(approxStatisticSIR[time,])/mean(approxStatisticSIR[time,]))^2,
                 MSE = mean((approxStatisticSIR[time,]-trueStatistic[time])^2),
                 RMSE = sqrt(mean((approxStatisticSIR[time,]-trueStatistic[time])^2)))
  for(type in c("Bias","Var", "CoefVar", "RelVar", "RelVar2", "MSE","RMSE")) {
    dfRes <- rbind(dfRes, c(Algorithm = "SIS",
                            Time = time,
                            Type = type,
                            Value = estSIS[[type]]),
                   stringsAsFactors = FALSE)
    dfRes <- rbind(dfRes, c(Algorithm = "SIR",
                            Time = time,
                            Type = type,
                            Value = estSIR[[type]]),
                   stringsAsFactors = FALSE)
  }
}
colnames(dfRes) <- c("Algorithm", "Time", "Type", "Value")
#Other relative variance (Definition 2) of http://www.statisticshowto.com/relative-variance/
for(time in 1:n) {
  variances <- as.numeric(dfRes[dfRes$Algorithm == "SIS" & dfRes$Type=="Var","Value"])
  dfRes <- rbind(dfRes, c(Algorithm = "SIS",
                          Time = time,
                          Type = "RelVar3",
                          Value = as.numeric(dfRes[dfRes$Algorithm == "SIS" &
                                                                    dfRes$Type=="Var" &
                                                                    dfRes$Time == time,"Value"])/sum(variances)),
                 stringsAsFactors = FALSE)
  dfRes <- rbind(dfRes, c(Algorithm = "SIR",
                          Time = time,
                          Type = "RelVar3",
                          Value = as.numeric(dfRes[dfRes$Algorithm == "SIR" &
                                                     dfRes$Type=="Var" &
                                                     dfRes$Time == time,"Value"])/sum(variances)),
                 stringsAsFactors = FALSE)
}

########
# PLOT #
########

ggplot(dfRes[dfRes$Type != "Var" & dfRes$Type != "MSE" &
               dfRes$Type != "Bias" & dfRes$Type != "RMSE",], aes(x = as.numeric(Time), y=as.numeric(Value), group=interaction(Algorithm,Type),
                    linetype = Algorithm, color = Type)) + geom_line()

