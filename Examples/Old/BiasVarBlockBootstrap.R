#This experiment shows that the bias is not homogeneous in the component
#in one block
rm(list=ls())
require(SMC5)
require(parallel)
require(ggplot2)

#####################
# SET UP EXPERIMENT #
#####################

set.seed(8917, "L'Ecuyer-CMRG")
n <- 10 #Filter time. Note: the experiment is still done in 1 time step
N <- 500 #Number of particles
dimension <- 3
A <- generateA(c(0.5,0.2),dimension)
repetitions <- 3000
type_statistic <- "sum_squared"
comp_statistic <- 3
ncores <- 2

##################
# GENERATE MODEL #
##################
X0 <- rep(10,dimension) #Just so that it is far from 0 and we don't have problem when estimate relative variance
sigmaX <- diag(rep(1,dimension))
sigmaY <- diag(rep(1,dimension)) #With a lower variance it should be easier for the algorithm to approximate the statistic
dataset <- generateData(n, A, X0, sigmaX, sigmaY)
X_data <- matrix(dataset$X_data, ncol=dimension)

####################
# RUN BOOTSTRAP PF #
####################
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
  logWeights <- array(log(1/N), c(n-1,1,N))
  blocks <- split(1:dimension,ceiling((1:dimension)/1))
  logWeightsBlock <- array(log(1/N), c(n-1,length(blocks),N))
  return(list(SIRRes = bootstrapParticleFilterOnline(N=N, n=n,
                                                   particles = particles, logWeights = logWeights,
                                                   fParams = fParams, gParams = gParams),
              BlockRes = blockParticleFilterOnline(N=N, n=n, blocks = blocks,
                                             particles = particles, logWeights = logWeightsBlock,
                                             fParams = fParams, gParams = gParams),
              blocks = blocks
              ))
}, mc.cores = ncores)

####################
# APPROX STATISTIC #
####################

approxStatisticSIR <- matrix(NA, nrow=1, ncol=repetitions)
approxStatisticBlock <- matrix(NA, nrow=1, ncol=repetitions)
approxStatisticSIR[1,] <- sapply(expRes, function(res) {
  return(computeApproxStatisticFilter(particles = res$SIRRes$filteringParticle,
                                      logWeights = res$SIRRes$filteringLogWeights,
                                      n = n,
                                      type = type_statistic,
                                      comp = comp_statistic))
})
approxStatisticBlock[1,] <- sapply(expRes, function(res) {
  return(computeApproxStatisticFilter(particles = res$BlockRes$filteringParticle,
                                      logWeights = res$BlockRes$filteringLogWeights,
                                      n = n,
                                      blocks = res$blocks,
                                      type = type_statistic,
                                      comp = comp_statistic))
})

##################
# TRUE STATISTIC #
##################

trueStatistic <- computeTrueStatisticFilter(X_data = X_data,
                                            n = n, type = type_statistic, comp = comp_statistic)

print(var(as.vector(approxStatisticSIR)))
print(var(as.vector(approxStatisticBlock)))

print(abs(trueStatistic-rowMeans(approxStatisticSIR)))
print(abs(trueStatistic-rowMeans(approxStatisticBlock)))
