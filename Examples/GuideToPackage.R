# Guide to The Package
# A short guide to all the available functions in the package! (A vignette can be derived from this)
#

require(rbenchmark)
require(SMC5)

#####################
# UTILITY FUNCTIONS #
#####################

#Several C++ functions have been written to be called by other C++ or because
#faster than their R counterpart.

# logAddition
# Performs addition of x+y where x and y are in log scale. Returns the log.
logAddition(log(0.0002), log(0.000008)) - log(0.0002+0.000008) #Notice that it is more precise!

# logAdditionSum
# As above but with a vector as input
logAdditionSum(c(log(0.0002), log(0.000008))) - log(0.0002+0.000008)

# ProbSampleReplace
# Is an alternative to sample with replacement but in C++ (as it needs to be called by other C++)
sample(10,size=2,replace = TRUE, prob=1:10)
ProbSampleReplace(10, size = 2, prob = 1:10)
benchmark(sample(10000,size=10000,replace = TRUE, prob=1:10000),
          ProbSampleReplace(10000, size = 10000, prob = 1:10000)) #The R implementation is faster as it is probably made in FORTRAN.

#########################
# GAUSSIAN DISTRIBUTION #
#########################

# Computing density (also in log scale) or sampling from a Gaussian have been implemented
# in C++
require(MASS)
require(mvtnorm)

all.equal(dnorm(0.5,0.7,2,log=TRUE),dnrmArma(0.5,0.7,2,logd=TRUE))
all.equal(dmvnorm(1:10,rep(5,10),diag(rep(1,10))+0.5,log=TRUE),
          as.vector(dmvnrmArma(matrix(1:10,nrow=1),rep(5,10),diag(rep(1,10))+0.5,logd=TRUE))) #Requires a matrix as input

# Speed comparison in computing univariate gaussian
benchmark(dnorm(0.5,0.7,2,log=TRUE),dnrmArma(0.5,0.7,2,logd=TRUE), replications = 10000) #Negligible
benchmark(rnorm(100000,0.7,2),rnormArma(100000,0.7,2)) #C++ implementation is 1.5x faster

# Speed comparison in computing multivariate gaussian
benchmark(dmvnorm(1:10,rep(5,10),diag(rep(1,10))+0.5,log=TRUE),
          as.vector(dmvnrmArma(matrix(1:10,nrow=1),rep(5,10),diag(rep(1,10))+0.5,logd=TRUE)), replications = 1000) #C++ implementation is 10x faster

benchmark(mvrnorm(10000,rep(5,10),diag(rep(1,10))+0.5),
          mvrnormArma(10000,rep(5,10),diag(rep(1,10))+0.5),
          rmvnorm(10000,rep(5,10),diag(rep(1,10))+0.5)) #C++ implementation is 2.5x faster

# NOTE: If the covariance matrix is diagonal, it can be said to the C++ implementation to make it even faster
all.equal(dmvnorm(1:10,rep(5,10),diag(rep(2,10)),log=TRUE),
          as.vector(dmvnrmArma(matrix(1:10,nrow=1),rep(5,10),diag(rep(2,10)),logd=TRUE)),
          as.vector(dmvnrmArma(matrix(1:10,nrow=1),rep(5,10),diag(rep(2,10)),logd=TRUE, diag = TRUE)))

benchmark(dmvnorm(1:10,rep(5,10),diag(rep(2,10)),log=TRUE),
          as.vector(dmvnrmArma(matrix(1:10,nrow=1),rep(5,10),diag(rep(2,10)),logd=TRUE)),
          as.vector(dmvnrmArma(matrix(1:10,nrow=1),rep(5,10),diag(rep(2,10)),logd=TRUE, diag = TRUE)), replications = 1000)
benchmark(mvrnorm(10000,rep(5,10),diag(rep(2,10))),
          mvrnormArma(10000,rep(5,10),diag(rep(2,10))),
          mvrnormArma(10000,rep(5,10),diag(rep(2,10)), diag = TRUE),
          rmvnorm(10000,rep(5,10),diag(rep(2,10))))

##################
# MODEL BUILDING #
##################

# The model is of the following form:
# X0 ~ delta(x0) with x0 given
# Xt ~ N(AXt-1, sigmaX) where sigmaX must be diagonal
# Yt ~ N(Xt, sigmaY) where sigmaY must be diagonal

# The matrix A can theoretically be of any form, but it is common to be symmetric and
# with periodic condition on the boundary.
# It can be generated with the function generateA which takes as input the variane of the compoennt
# and the correlation with its neighbour
dimension <- 8
neighCorr <- c(0.8,0.4,0.1)
A <- generateA(neighCorr,dimension)
print(A)

# The data can be generated using the generateData function
# Note that X0 is NOT contained in X_data but it is in X0
n <- 10 #Time (Excluding 0)
X0 <- 1:dimension
dataset <- generateData(n, A, X0, diag(rep(1,dimension)), diag(rep(1,dimension)))
X_data <- dataset$X_data
Y_data <- dataset$Y_data

# NOTE: that with this choice of the matrix A and sigmaX, the process is pretty much deterministic
# and after few time steps also the Y are almost perfect.

###########################
# KALMAN FILTER/SMOOTHING #
###########################

# The package provides Kalman filter and smoothing implementation in C++

#Generate data
rm(list=ls())
dimension <- 5
n <- 10
A <- generateA(c(0.5,0.2), dimension)
X0 <- rep(0,dimension)
sigmaX <- diag(rep(1,dimension))
sigmaY <- diag(rep(1,dimension))
dataset <- generateData(n, A, X0, sigmaX, sigmaY)
X_data <- dataset$X_data
Y_data <- dataset$Y_data

#Kalman filter - notice that the first element is the one at time 0!
kalmanFilterRes <- KalmanFilterCpp(m_0 = X0,
                                   C_0 = sigmaX,
                                   F_matrix = diag(dimension),
                                   G = A,
                                   V = sigmaY,
                                   W = sigmaX,
                                   y = Y_data)

#Kalman smoothing
kalmanSmoothingRes <- KalmanSmoothingCpp(m_0 = X0,
                                         C_0 = sigmaX,
                                         F_matrix = diag(dimension),
                                         G = A,
                                         V = sigmaY,
                                         W = sigmaX,
                                         y = Y_data) #Note: inside it execute again the kalman filter

#This is a short example on how to estimate the smoothing
time_smoothing <- 5
reps <- 1000
res <- matrix(NA, nrow=reps, ncol = dimension)
sigmaY <- diag(rep(0.01,dimension)) #Reduce variance of Y to get better results with less repetitions (otherwise the number of repetitions needed might be huge!)
for(rep in 1:reps) {
  #Generate new Y with the same X
  Y_data_rep <- generateY(X_data,sigmaY)
  kalmanSmoothingRes_rep <- KalmanSmoothingCpp(m_0 = X0,
                                           C_0 = sigmaX,
                                           F_matrix = diag(dimension),
                                           G = A,
                                           V = sigmaY,
                                           W = sigmaX,
                                           y = Y_data_rep)
  res[rep,] <- KalmanSamplerCpp(1, kalmanSmoothingRes_rep, time_smoothing)
}

print(colMeans(res))
print(X_data[time_smoothing,])

##############################################################
# SEQUENTIAL IMPORTANCE SAMPLING / BOOTSTRAP PARTICLE FILTER #
##############################################################

# Generate model
rm(list=ls())
dimension <- 10
N <- 5 #Number of particles
n <- 4 #Time horizon (filter is returned for all time steps except time 0)
A <- generateA(c(0.5,0.2), dimension)
X0 <- rep(0,dimension)
sigmaX <- diag(rep(1,dimension))
sigmaY <- diag(rep(1,dimension))
dataset <- generateData(n, A, X0, sigmaX, sigmaY)

# The algorithms requires these two lists:
fParams <- list(A=A, X0=X0, sigmaX=sigmaX)
gParams <- list(Y=dataset$Y_data, sigmaY=sigmaY)

# Sequential Importance Sampling
seqImpSampRes <- sequentialImportanceSampling(N=N, n=n, fParams = fParams,
                                              gParams = gParams)
#Notice that the weights are not standardized

# Bootstrap Particle Filter
bootstrapPFRes <- bootstrapParticleFilter(N=N, n=n, fParams = fParams,
                                          gParams = gParams)

# ONLINE: There is also a function that allows for an online version of
# the algorithms

#Generate new additional data
n_new <- 2*n
dataset2 <- generateData(n_new-n, A, dataset$X_data[nrow(dataset$X_data),], sigmaX, sigmaY)

gParams2 <- list(Y=rbind(dataset$Y_data,dataset2$Y_data), sigmaY=sigmaY)
seqImpSampRes2 <- sequentialImportanceSamplingOnline(N=N, n=n_new,
                                                    particles = seqImpSampRes$filteringParticle,
                                                    logWeights = seqImpSampRes$filteringLogWeights,
                                                    fParams = fParams,
                                                    gParams = gParams2)
bootstrapPFRes2 <- bootstrapParticleFilterOnline(N=N, n=n_new,
                                                     particles = seqImpSampRes$filteringParticle,
                                                     logWeights = seqImpSampRes$filteringLogWeights,
                                                     fParams = fParams,
                                                     gParams = gParams2)

benchmark(bootstrapParticleFilterOnline(N=N, n=n_new,
                                        particles = seqImpSampRes$filteringParticle,
                                        logWeights = seqImpSampRes$filteringLogWeights,
                                        fParams = fParams,
                                        gParams = gParams2),
          bootstrapParticleFilter(N=N, n=n_new, fParams = fParams,
                                  gParams = gParams2), replications = 1000)

###########################
# BLOCKED PARTICLE FILTER #
###########################

#Blocked particle filter requires to define a list with all the blocks
#within eachother there has to be the components of the variable in the block
blocks <- list(1:3,c(4,8),5,c(6,7,9,10))
#There is no need for the  blocks to have the same number of compoennts or for them
#to be contiguous
blockPFRes <- blockParticleFilter(N=N, n=n, blocks = blocks,
                                  fParams = fParams, gParams = gParams)
#The logweights are for each block

#NOTE: The most convinent way to create contiguous block of equal cardinality is:
card <- 3
blocks2 <- split(1:dimension,ceiling((1:dimension)/card))
print(blocks2)

#########################
# GIBBS PARTICLE FILTER #
#########################

# Gibbs particle filter requires to define a number m of sweeps+
# as well as a radius
m <- 100
radius <- 1

gibbsPFRes <- gibbsParticleFilter(N=N, n=n, m=m, radius=radius,
                                  fParams = fParams, gParams=gParams)
#Note that this function returns a list with particle and weihgts to be consistent
#with all the other functions. Hoewver, the weights are the same for all the particles.

#An online version is also available
gibbsPFRes2 <- gibbsParticleFilterOnline(N=N, n=n_new, m=m, radius=radius,
                                         particles = gibbsPFRes$filteringParticle,
                                         fParams = fParams, gParams=gParams2)
benchmark(gibbsParticleFilter(N=N, n=n_new, m=m, radius=radius,
                              fParams = fParams, gParams=gParams2),
          gibbsParticleFilterOnline(N=N, n=n_new, m=m, radius=radius,
                                    particles = gibbsPFRes$filteringParticle,
                                    fParams = fParams, gParams=gParams2))

#There is an old implementation that still works but is slower and does not work
#in an on-line fashion.
gibbsPFResOld <- gibbsParticleFilterSMC4(N=N, n=n, m=m, radius=radius,
                                         fParams = fParams, gParams=gParams)
benchmark(gibbsParticleFilterSMC4(N=N, n=n, m=m, radius=radius,
                                  fParams = fParams, gParams=gParams),
          gibbsParticleFilter(N=N, n=n, m=m, radius=radius,
                                  fParams = fParams, gParams=gParams))


############################
# APPROXIMATING STATISTICS #
############################

# The function computeApproxStatisticFilter provides
# a way to compute 4 different statistics for the filter corresponding to the following type:
# 1. "mean" -> expected value of a specific component
# 2. "sum" -> sum over all the component of their expected value
# 3. "mean_squared" -> expected value of the square specific component
# 4. "sum_squared " -> sum over all the component of the expected value of their square

#True statistic is just computed on the dataset that generated the data
trueStatistics <- list(meanStatistic = computeTrueStatisticFilter(X_data = dataset$X_data, n = 2,
                                                                 type="mean", comp = 3),
                      sumStatistic = computeTrueStatisticFilter(X_data = dataset$X_data, n = 2,
                                                                 type="sum"),
                      meanSquaredStatistic = computeTrueStatisticFilter(X_data = dataset$X_data, n = 2,
                                                                 type="mean_squared", comp = 3),
                      SumSquaredStatistic = computeTrueStatisticFilter(X_data = dataset$X_data, n = 2,
                                                                 type="sum_squared"))
#This is the estimate given by the Kalman filter
approxStatisticKalman <- list(meanStatistic = computeKalmanStatisticFilter(fParams = fParams, gParams = gParams,
                                                                           n = 2, type = "mean", comp = 3),
                              sumStatistic = computeKalmanStatisticFilter(fParams = fParams, gParams = gParams,
                                                                           n = 2, type = "sum"),
                              meanSquaredStatistic = computeKalmanStatisticFilter(fParams = fParams, gParams = gParams,
                                                                           n = 2, type = "mean_squared", comp = 3),
                              sumSquaredStatistic = computeKalmanStatisticFilter(fParams = fParams, gParams = gParams,
                                                                           n = 2, type = "sum_squared"))

approxStatisticsBootstrap <- list(meanStatistic = computeApproxStatisticFilter(particles = bootstrapPFRes$filteringParticle,
                                                                      logWeights = bootstrapPFRes$filteringLogWeights,
                                                                      n = 2, type = "mean", comp = 3),
                                  sumStatistic = computeApproxStatisticFilter(particles = bootstrapPFRes$filteringParticle,
                                                                               logWeights = bootstrapPFRes$filteringLogWeights,
                                                                               n = 2, type = "sum"),
                                  meanSquaredStatistic = computeApproxStatisticFilter(particles = bootstrapPFRes$filteringParticle,
                                                                               logWeights = bootstrapPFRes$filteringLogWeights,
                                                                               n = 2, type = "mean_squared", comp = 3),
                                  sumSquaredStatistic = computeApproxStatisticFilter(particles = bootstrapPFRes$filteringParticle,
                                                                               logWeights = bootstrapPFRes$filteringLogWeights,
                                                                               n = 2, type = "sum_squared"))

approxStatisticBlock <-  list(meanStatistic = computeApproxStatisticFilter(particles = blockPFRes$filteringParticle,
                                                                           logWeights = blockPFRes$filteringLogWeights,
                                                                           n = 2, blocks = blocks, type = "mean", comp = 3),
                              sumStatistic = computeApproxStatisticFilter(particles = blockPFRes$filteringParticle,
                                                                           logWeights = blockPFRes$filteringLogWeights,
                                                                           n = 2, blocks = blocks, type = "sum"),
                              meanSquaredStatistic = computeApproxStatisticFilter(particles = blockPFRes$filteringParticle,
                                                                           logWeights = blockPFRes$filteringLogWeights,
                                                                           n = 2, blocks = blocks, type = "mean_squared", comp = 3),
                              sumSquaredStatistic = computeApproxStatisticFilter(particles = blockPFRes$filteringParticle,
                                                                           logWeights = blockPFRes$filteringLogWeights,
                                                                           n = 2, blocks = blocks, type = "sum_squared"))
approxStatisticGibbsPF <-  list(meanStatistic = computeApproxStatisticFilter(particles = gibbsPFRes$filteringParticle,
                                                                           logWeights = gibbsPFRes$filteringLogWeights,
                                                                           n = 2, type = "mean", comp = 3),
                              sumStatistic = computeApproxStatisticFilter(particles = gibbsPFRes$filteringParticle,
                                                                          logWeights = gibbsPFRes$filteringLogWeights,
                                                                          n = 2, type = "sum"),
                              meanSquaredStatistic = computeApproxStatisticFilter(particles = gibbsPFRes$filteringParticle,
                                                                                  logWeights = gibbsPFRes$filteringLogWeights,
                                                                                  n = 2, type = "mean_squared", comp = 3),
                              sumSquaredStatistic = computeApproxStatisticFilter(particles = gibbsPFRes$filteringParticle,
                                                                                 logWeights = gibbsPFRes$filteringLogWeights,
                                                                                 n = 2, type = "sum_squared"))
print(trueStatistics)
print(approxStatisticKalman)
print(approxStatisticsBootstrap)
print(approxStatisticBlock)
print(approxStatisticGibbsPF)

# The function computeDfBiasVar can also be used to quickly generate a dataset with bias/var/MSE and other
# information of the estimator obtained through several replications. See CurseOfTime.R or CurseOfDimensionality.R
# to see how it is used.
