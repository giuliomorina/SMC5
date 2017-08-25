#This experiment shows that the bias is not homogeneous in the component
#in one block
rm(list=ls())
require(SMC5)
require(parallel)
require(ggplot2)

######################
# SETTING EXPERIMENT #
######################

set.seed(1771, "L'Ecuyer-CMRG")
n <- 10 #Filter time. Note: the experiment is still done in 1 time step
N <- 2 #Number of particles
dimension <- 12
#A <- generateA(c(0.1,0.3,0.2,0.1,0.05), dimension)
A <- generateA(c(0.1,0.4,0.3), dimension)
blocks <- list(c(1:3,11,12),4:10)
#blocks <- list(1:12)
repetitions <- 3000
type_statistic <- "mean"
ncores <- 2

##################
# GENERATE MODEL #
##################
X0 <- rep(10,dimension) #Just so that it is far from 0 and we don't have problem when estimate relative variance
sigmaX <- diag(rep(0.0001,dimension))
sigmaY <- diag(rep(1,dimension)) #With a lower variance it should be easier for the algorithm to approximate the statistic
dataset <- generateData(n, A, X0, sigmaX, sigmaY)
X_data <- matrix(dataset$X_data, ncol=dimension)

################
# RUN BLOCK PF #
################
expRes <- mclapply(1:repetitions, function(rep) {
  #Generate Y always with the same X
  Y_data <- generateY(X_data, sigmaY)
  fParams <- list(A=A, X0=X0, sigmaX=sigmaX)
  gParams <- list(Y=Y_data, sigmaY=sigmaY)
  #Kalman Filter
  kalmanFilterRes <- KalmanFilterCpp(m_0 = X0,
                                     C_0 = sigmaX,
                                     F_matrix = diag(dimension),
                                     G = A,
                                     V = sigmaY,
                                     W = sigmaX,
                                     y = Y_data)
  #Give particles sampled from the Kalman at the previous times 0,1,...,n-1
  particles <- array(NA, c(n-1,dimension,N))
  for(slice in 1:N) {
    particles[,,slice] <- mvrnormArma(1, mean = kalmanFilterRes$m[,n], sigma = kalmanFilterRes$C[,,n])
  }
  logWeights <- array(log(1/N), c(n-1,length(blocks),N))
  return(list(BlockRes = blockParticleFilterOnline(N=N, n=n,
                                                   particles = particles, logWeights = logWeights,
                                                   blocks = blocks,
                                                   fParams = fParams, gParams = gParams),
              blocks = blocks))
}, mc.cores = ncores)

#########################
# APPROXIMATE STATISTIC #
#########################

approxStatisticBlock <- matrix(NA, nrow=dimension, ncol=repetitions)
for(rep in 1:repetitions) {
  for(comp in 1:dimension) {
    approxStatisticBlock[comp,rep] <- computeApproxStatisticFilter(particles =
                                          expRes[[rep]]$BlockRes$filteringParticle,
                                          logWeights = expRes[[rep]]$BlockRes$filteringLogWeights,
                                          n = n, blocks = expRes[[rep]]$blocks,
                                          type = type_statistic, comp = comp)
  }
}

###################
# TRUE STATISTICS #
###################
trueStatistics <- rep(NA, dimension)
kalmanStatistics <- rep(NA, dimension)
for(comp in 1:dimension) {
  trueStatistics[comp] <- computeTrueStatisticFilter(X_data = X_data, n = n, type =  type_statistic,
                                                     comp = comp)
}

########################
# COMPUTE BIAS/MSE/VAR #
########################

dfRes <- computeDfBiasVar(approxStatisticPF = approxStatisticBlock, trueStatistics = trueStatistics,
                          algorithmName = "Block", dependentVarName = "Component")

########
# PLOT #
########
#Variance
#Just to check the different type of relative variances:
ggplot(dfRes[dfRes$Type != "Var" & dfRes$Type != "MSE" &
               dfRes$Type != "Bias" & dfRes$Type != "RMSE" &
               dfRes$Type != "RelAbsBias",], aes(x = as.numeric(Component), y=as.numeric(Value), color = Type)) + geom_line()

#Bias
ggplot(dfRes[dfRes$Type == "Bias",], aes(x = as.numeric(Component), y=abs(Value), color = Type)) + geom_line()
