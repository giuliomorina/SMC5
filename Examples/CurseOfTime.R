# CURSE OF TIME
# The curse of time appears also in dimension 1 meaning that lots of weight is only on one particle
# and the variance of the estimate grows exponentially with time! Particle filters do not show
# the same behaviour.

rm(list=ls())
require(SMC5)
require(parallel)
require(ggplot2)
require(latex2exp)

######################
# SETTING EXPERIMENT #
######################

set.seed(1717,"L'Ecuyer-CMRG")
n <- 50
ncores <- 2
repetitions <- 100
N <- 20
type_statistic_plot <- "sum_squared" #Which statistic to plot
comp_statistic_plot <- 2 #Which component of the statistic (must be 1 for sum and sum squared)
if(type_statistic_plot == "sum" || type_statistic_plot == "sum_squared") {
  comp_statistic_plot <- 1
}

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

dfResList <- lapply(c("mean","sum","mean_squared","sum_squared"), function(type_statistic){
  if(type_statistic == "sum" || type_statistic == "sum_squared") {
    possible_comp_statistic <- NA
  } else {
    possible_comp_statistic <- 1:dimension
  }
  dfResComp <- lapply(possible_comp_statistic, function(comp_statistic) {
    dataStatistics <- rep(NA, n)
    approxStatisticSIS <- matrix(NA, nrow=n, ncol=repetitions)
    approxStatisticSIR <- matrix(NA, nrow=n, ncol=repetitions)
    for(time in 1:n) {
      dataStatistics[time] <- computeTrueStatisticFilter(X_data = X_data, n = time, type =  type_statistic,
                                                         comp = comp_statistic)
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
    return(list(SIR=approxStatisticSIR,
                SIS=approxStatisticSIS,
                dataStatistics=dataStatistics))
  })
  return(dfResComp)
})
names(dfResList) <- c("mean","sum","mean_squared","sum_squared")

########################
# COMPUTE BIAS/MSE/VAR #
########################

dfResVarData <- computeDfVar(approxList = dfResList[[type_statistic_plot]][[comp_statistic_plot]][c("SIR","SIS")],
                             trueStatistics = dfResList[[type_statistic_plot]][[comp_statistic_plot]]$dataStatistics,
                             dependentVar = 1:n,
                             dependentVarName = "Time")

dfResBiasMSEData <- computeDfBiasMSE(approxList = dfResList[[type_statistic_plot]][[comp_statistic_plot]][c("SIR","SIS")],
                                     trueStatistics = dfResList[[type_statistic_plot]][[comp_statistic_plot]]$dataStatistics,
                                     dependentVar = 1:n,
                                     dependentVarName = "Time")


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
                        Repetition = numeric(),
                        Value = numeric(),
                        stringsAsFactors = FALSE)
for(rep in 1:repetitions) {
  for(time in 1:n) {
    dfWeights <- rbind(dfWeights,c(Algorithm = "SIS",
                                   Time = time,
                                   Repetition = rep,
                                   Value = maxWeightSIS[time,rep]),
                       stringsAsFactors = FALSE)
    dfWeights <- rbind(dfWeights,c(Algorithm = "SIR",
                                   Time = time,
                                   Repetition = rep,
                                   Value = maxWeightSIR[time,rep]),
                       stringsAsFactors = FALSE)
  }
}
colnames(dfWeights) <- c("Algorithm","Time","Repetition", "Value")
dfWeights$Algorithm <- as.factor(dfWeights$Algorithm)

########
# PLOT #
########

#Variance
#Just to check the different type of relative variances:
ggplot(dfResVarData[dfResVarData$Type != "Var",], aes(x = as.numeric(Time), y=as.numeric(Value), group=interaction(Algorithm,Type),
                                                      linetype = Algorithm, color = Type)) + geom_line()

#RelVar
ggplot(dfResVarData[dfResVarData$Type == "RelVar",], aes(x = Time, y=Value, color = Algorithm)) +
  geom_point(size = 3, shape = 4) +  scale_colour_brewer(palette = "Set1") +
  theme_grey(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Relative Variance") + xlab("n") +
  geom_smooth(method="lm",se=FALSE, size = 1.5) #Note that this is the relative variance, not the variance!


#Maximum Weights
ggplot(dfWeights, aes(x=as.numeric(Time), y=as.numeric(Value), color=Algorithm)) +
  geom_hline(yintercept = 1, linetype="dashed") +
  theme_grey(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Weight") + xlab("n") +
  stat_summary(geom="pointrange", fun.data=mean_cl_boot, fun.args=list(conf.int=1))+
  scale_colour_brewer(palette = "Set1")
