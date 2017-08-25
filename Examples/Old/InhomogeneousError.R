# INHOMOEGENOUS ERROR
# The blocked PF induces a bias that grows as the cardinality of the block decreases. Viceversa,
# the variance grows as the cardinality of the block increases.

rm(list=ls())
require(SMC5)
require(parallel)
require(ggplot2)

######################
# SETTING EXPERIMENT #
######################

set.seed(177, "L'Ecuyer-CMRG")
n <- 10  #Filter time. NOTE: The experiment is still done in 1 time step
N <- 100 #Number of particles
dimension <- 11
A <- generateA(c(0.1,0.3,0.2,0.1,0.05), dimension)
#A <- generateA(c(0.1,0.4,0.3), dimension)
A <- diag(dimension)

possible_cardinalities <- 1:dimension #Needs to start from 1 and be contigouos
repetitions <- 100
type_statistic <- "sum_squared"
comp_statistic <- 6 #Used for type_statistic mean or mean_squared
ncores <- 2

##################
# GENERATE MODEL #
##################
X0 <- rep(0,dimension) #Just so that it is far from 0 and we don't have problem when estimate relative variance
sigmaX <- diag(rep(1,dimension))
sigmaY <- diag(rep(0.00001,dimension)) #With a lower variance it should be easier for the algorithm to approximate the statistic
dataset <- generateData(n, A, X0, sigmaX, sigmaY)
X_data <- matrix(dataset$X_data, ncol=dimension)

##########################
# COMPUTE TRUE STATISTIC #
##########################
trueStatistic <-  computeTrueStatisticFilter(X_data = X_data, n = n, type =  type_statistic,
                                             comp = comp_statistic)

############################################
# RUN BLOCK PF FOR DIFFERENT CARDINALITIES #
############################################

expRes <- mclapply(1:repetitions, function(rep) {
  #Generate Y always with the same X
  Y_data <- generateY(X_data, sigmaY)
  fParams <- list(A=A, X0=X0, sigmaX=sigmaX)
  gParams <- list(Y=Y_data, sigmaY=sigmaY)
  #Just give the true value until time n-1
  particles <- array(NA, c(n-1,dimension,N))
  for(slice in 1:N) {
    particles[,,slice] <- X_data[(1:n-1),]
  }
  res <- lapply(possible_cardinalities, function(card) {
    blocks <- split(1:dimension,ceiling((1:dimension)/card))
    logWeights <- array(log(1/N), c(n-1,length(blocks),N))
    return(list(BlockRes = blockParticleFilterOnline(N=N, n=n,
                               particles = particles, logWeights = logWeights,
                               blocks = blocks,
                               fParams = fParams, gParams = gParams),
           blocks = blocks,
           card = card))
  })
  return(res)
}, mc.cores = ncores)

#########################
# APPROXIMATE STATISTIC #
#########################

approxStatisticBlock <- matrix(NA, nrow=length(possible_cardinalities), ncol=repetitions)
for(rep in 1:repetitions) {
  approxStatisticBlock[,rep] <- sapply(expRes[[rep]], function(res) {
    return(computeApproxStatisticFilter(particles = res$BlockRes$filteringParticle,
                                        logWeights = res$BlockRes$filteringLogWeights,
                                        n = n, blocks = res$blocks,
                                        type = type_statistic, comp = comp_statistic))
  })
}

########################
# COMPUTE BIAS/MSE/VAR #
########################

dfRes <- computeDfBiasVar(approxStatisticPF = approxStatisticBlock, trueStatistics = rep(trueStatistic, length(possible_cardinalities)),
                          algorithmName = "Block", dependentVarName = "Cardinality")

########
# PLOT #
########

#Variance
#Just to check the different type of relative variances:
ggplot(dfRes[dfRes$Type != "Var" & dfRes$Type != "MSE" &
               dfRes$Type != "Bias" & dfRes$Type != "RMSE" &
               dfRes$Type != "RelAbsBias",], aes(x = as.numeric(Cardinality), y=as.numeric(Value), color = Type)) + geom_line()

#Variance
ggplot(dfRes[dfRes$Type == "Var",], aes(x = as.numeric(Cardinality), y=as.numeric(Value), color = Type)) + geom_line()

#RMSE
ggplot(dfRes[dfRes$Type == "RMSE",], aes(x = as.numeric(Cardinality), y=abs(Value), color = Type)) + geom_line()

#Bias
ggplot(dfRes[dfRes$Type == "RelAbsBias" | dfRes$Type == "Bias",], aes(x = as.numeric(Cardinality), y=abs(Value), color = Type)) + geom_line()
