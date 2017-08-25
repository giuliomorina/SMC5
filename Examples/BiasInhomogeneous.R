# BLOCKED PF INDUCES BIAS
# This experiments show that the blocked PF induces ADDITIONAL bias (as the PF is not unbiased when
# approximating integral) but reduces variance. Moreover, the induces bias is not homogeneous in the
# components as it depends on the block structure.
# NOTE: in this experiment we keep both X and Y fix! See THE BOOTSTRAP PARTICLE FILTERING BIAS - Olsson

rm(list=ls())
require(SMC5)
require(parallel)
require(ggplot2)

######################
# SETTING EXPERIMENT #
######################
set.seed(17998, "L'Ecuyer-CMRG")
n <- 20 #Filter time. Note: the experiment is still done in 1 time step
N <- 1000 #Number of particles
dimension <- 6
A <- generateA(c(0.5,0.2,0.1),dimension)
card <- 3
blocks <- split(1:dimension,ceiling((1:dimension)/card))
repetitions <- 100
type_statistic <- "mean_squared"
ncores <- 2

##################
# GENERATE MODEL #
##################
X0 <- rep(10,dimension) #Just so that it is far from 0 and we don't have problem when estimate relative variance
sigmaX <- diag(rep(1,dimension))
sigmaY <- diag(rep(1,dimension)) #With a lower variance it should be easier for the algorithm to approximate the statistic
dataset <- generateData(n, A, X0, sigmaX, sigmaY)
X_data <- matrix(dataset$X_data, ncol=dimension)
Y_data <- matrix(dataset$Y_data, ncol=dimension)
fParams <- list(A=A, X0=X0, sigmaX=sigmaX)
gParams <- list(Y=Y_data, sigmaY=sigmaY)

######################
# PREVIOUS PARTICLES #
######################

#As the experiment is done in 1 time step, we pass as previous times particles
#taken from the Kalman filter distribution. Note: in each experiment different
#particles are passed
kalmanFilterRes <- KalmanFilterCpp(m_0 = X0,
                                   C_0 = sigmaX,
                                   F_matrix = diag(dimension),
                                   G = A,
                                   V = sigmaY,
                                   W = sigmaX,
                                   y = Y_data)

####################
# KALMAN STATISTIC #
####################

kalmanStatistics <- rep(NA, dimension)
dataStatistics <- rep(NA, dimension)
for(comp in 1:dimension) {
  kalmanStatistics[comp] <- computeKalmanStatisticFilter(fParams = fParams, gParams = gParams, n=n,
                                                         type = type_statistic, comp=comp)
  dataStatistics[comp] <- computeTrueStatisticFilter(X_data = X_data, n=n, type=type_statistic, comp=comp)
}

#######################
# RUN PARTICLE FILTER #
#######################

logWeightsBlocks <- array(log(1/N), c(n-1,length(blocks),N))
logWeights <- array(log(1/N), c(n-1,1,N))
expRes <- mclapply(1:repetitions, function(rep) {
  particles <- array(NA, c(n-1,dimension,N))
  for(slice in 1:N) {
    particles[,,slice] <- mvrnormArma(1, mean = kalmanFilterRes$m[,n], sigma = kalmanFilterRes$C[,,n])
  }
  return(list(SIRRes = bootstrapParticleFilter(N=N, n=n,
                                                     fParams = fParams, gParams = gParams),
              BlockRes = blockParticleFilter(N=N, n=n, blocks = blocks,
                                                   fParams = fParams, gParams = gParams)
  ))
}, mc.cores = ncores)

####################
# APPROX STATISTIC #
####################

approxStatisticSIR <- matrix(NA, nrow=dimension, ncol=repetitions)
approxStatisticBlock <- matrix(NA, nrow=dimension, ncol=repetitions)

for(comp in 1:dimension) {
  approxStatisticSIR[comp,] <- sapply(expRes, function(res) {
    return(computeApproxStatisticFilter(particles = res$SIRRes$filteringParticle,
                                        logWeights = res$SIRRes$filteringLogWeights,
                                        n = n,
                                        type = type_statistic,
                                        comp = comp))
  })
  approxStatisticBlock[comp,] <- sapply(expRes, function(res) {
    return(computeApproxStatisticFilter(particles = res$BlockRes$filteringParticle,
                                        logWeights = res$BlockRes$filteringLogWeights,
                                        n = n,
                                        blocks = blocks,
                                        type = type_statistic,
                                        comp = comp))
  })
}

print(apply(approxStatisticSIR,1,var))
print(apply(approxStatisticBlock,1,var))

print(abs(rowMeans(approxStatisticSIR)-kalmanStatistics))
print(abs(rowMeans(approxStatisticBlock)-kalmanStatistics))

print(apply(approxStatisticBlock,1,var) < apply(approxStatisticSIR,1,var))
print(abs(rowMeans(approxStatisticBlock)-kalmanStatistics) > abs(rowMeans(approxStatisticSIR)-kalmanStatistics))
