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
ncores <- 2
repetitions <- 10000
N <- 20
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
}, mc.cores = ncores)

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

dfRes <- computeDfBiasVar(approxStatisticPF = approxStatisticSIS, trueStatistics = trueStatistic,
                           algorithmName = "SIS", dependentVarName = "Time")
dfRes <- computeDfBiasVar(approxStatisticPF = approxStatisticSIR, trueStatistics = trueStatistic,
                           algorithmName = "SIR", dependentVarName = "Time", dfRes = dfRes)

########################
# MEAN MAXIMUM WEIGHTS #
########################

maxWeightSIS <- matrix(NA, nrow=n, ncol=repetitions)
maxWeightSIR <- matrix(NA, nrow=n, ncol=repetitions)
for(rep in 1:repetitions) {
  res <- expRes[[rep]]
  for(time in 1:n) {
    maxWeightSIS[time,rep] <- max(exp(res$SISRes$filteringLogWeights[time,1,]-logAdditionSum(res$SISRes$filteringLogWeights[time,1,])))
    maxWeightSIR[time,rep] <- max(exp(res$SIRRes$filteringLogWeights[time,1,]-logAdditionSum(res$SIRRes$filteringLogWeights[time,1,])))
  }
}

dfWeights <- data.frame(Algorithm = character(),
                        Time = numeric(),
                        Value = numeric(),
                        stringsAsFactors = FALSE)
for(time in 1:n) {
  dfWeights <- rbind(dfWeights,c(Algorithm = "SIS",
                                 Time = time,
                                 Value = mean(maxWeightSIS[time,])),
                     stringsAsFactors = FALSE)
  dfWeights <- rbind(dfWeights,c(Algorithm = "SIR",
                                 Time = time,
                                 Value = mean(maxWeightSIR[time,])),
                     stringsAsFactors = FALSE)
}
colnames(dfWeights) <- c("Algorithm","Time","Value")
dfWeights$Algorithm <- as.factor(dfWeights$Algorithm)

########
# PLOT #
########

#Variance
#Just to check the different type of relative variances:
ggplot(dfRes[dfRes$Type != "Var" & dfRes$Type != "MSE" &
               dfRes$Type != "Bias" & dfRes$Type != "RMSE" &
               dfRes$Type != "RelAbsBias",], aes(x = as.numeric(Time), y=as.numeric(Value), group=interaction(Algorithm,Type),
                    linetype = Algorithm, color = Type)) + geom_line()

#RelVar
ggplot(dfRes[dfRes$Type == "RelVar",], aes(x = as.numeric(Time), y=as.numeric(Value), color = Algorithm)) +
  geom_line(size = 1.2) +  scale_colour_brewer(palette = "Set1")


#Maximum Weights
ggplot(dfWeights, aes(x=as.numeric(Time), y=as.numeric(Value), color=Algorithm)) +
  geom_line(size = 1.2) + geom_hline(yintercept = 1, linetype="dashed") +
  scale_colour_brewer(palette = "Set1")
