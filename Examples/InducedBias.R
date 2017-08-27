# INDUCED BIAS
# This experiment shows that the block and the Gibbs PF algorithms induce an
# additional bias that is not present in the SIR algorithm.
# Note that the number of particles does not affect the induced bias of the block/Gibbs PF,
# but affects the "base" bias. A modest number of particles is then used.
# We still want to show the performance in one time step, but we can not do n=1
# as the dependency among all the components would not kick in already.
# So we give for time 1,...,n-1 particles sampled from the true filter
# and just use the algorithm for the nth step.

rm(list=ls())
require(SMC5)
require(parallel)
require(ggplot2)

######################
# SETTING EXPERIMENT #
######################
set.seed(878,"L'Ecuyer-CMRG")
n <- 10 #Filter time (NOTE: the experiment is still performed in one time step)
dimension <- 8
ncores <- 2
repetitions <- 10
N <- 50 #Number of particles
A <- generateA(c(0,0,0,0.6), dimension)
sigmaX <- diag(rep(1,dimension))
sigmaY <- diag(rep(1,dimension))
X0 <- rep(0, dimension)
possible_cardinalities <- 1:dimension
possible_radius <- (0:floor(dimension/2))
type_statistic <- "sum_squared"
comp_statistic <- 1 #Used for type_statistic mean or mean_squared
m <- 10 #Number of sweeps

###################
# GENERATE X DATA #
###################

dataset <- generateData(n, A, X0, sigmaX, sigmaY)
X_data <- matrix(dataset$X_data, ncol=dimension)
fParams <- list(A=A, X0=X0, sigmaX=sigmaX)

################################
# RUN BLOCK & APPROX STATISTIC #
################################

#We immediately approximate the statistics as otherwise it becomes too big.

resStats <- mclapply(1:repetitions, function(rep) {
  #Generate new Y data
  Y_data <- generateY(X_data, sigmaY)
  gParams <- list(Y=Y_data, sigmaY=sigmaY)
  #Generate Kalman Particles
  kalmanFilterRes <- KalmanFilterCpp(m_0 = X0,
                                     C_0 = sigmaX,
                                     F_matrix = diag(dimension),
                                     G = A,
                                     V = sigmaY,
                                     W = sigmaX,
                                     y = Y_data)
  kalmanParticles <- array(NA,c(n,dimension,N))
  for(aux in 1:n) {
    kalmanParticles[aux,,] <- t(mvrnormArma(N,mean = kalmanFilterRes$m[,aux+1], sigma = as.matrix(kalmanFilterRes$C[,,aux+1])))
  }
  #Approximate statistic for each cardinality
  res <- lapply(possible_cardinalities, function(card) {
    blocks <- split(1:dimension,ceiling((1:dimension)/card))
    #Run SIR
    BlockRes <- blockParticleFilterOnline(N=N, n=n, blocks = blocks,
                                          particles = kalmanParticles[1:(n-1),,],
                                          logWeights = array(log(1/N),c(n-1,length(blocks),N)),
                                   fParams = fParams, gParams = gParams)
    #Approx statistic
    BlockStat <- computeApproxStatisticFilter(particles = BlockRes$filteringParticle,
                                              logWeights = BlockRes$filteringLogWeights,
                                              n = n,
                                              blocks = blocks,
                                              type = type_statistic,
                                              comp = comp_statistic)

    return(list(BlockStat = BlockStat))
  })
  resGibbs <- lapply(possible_radius, function(radius) {
    #Rub Gibbs PF
    GibbsPFRes <- gibbsParticleFilter(N=N, n=n, m=m, radius=radius,
                                      fParams = fParams, gParams=gParams)
    #Approx statistic
    GibbsPFStat <- computeApproxStatisticFilter(particles = GibbsPFRes$filteringParticle,
                                              logWeights = GibbsPFRes$filteringLogWeights,
                                              n = n,
                                              type = type_statistic,
                                              comp = comp_statistic)

    return(list(GibbsPFStat = GibbsPFStat))
  })
  return(list(gParams=gParams,
              fParams = fParams,
              X_data = X_data,
         res=res,
         resGibbs = resGibbs))
}, mc.cores = ncores)

#########################
# APPROXIMATE STATISTIC #
#########################

approxStatisticBlock <- matrix(NA, nrow=length(possible_cardinalities), ncol=repetitions)
approxStatisticGibbsPF <- matrix(NA, nrow=length(possible_radius), ncol=repetitions)

for(rep in 1:repetitions) {
  approxStatisticBlock[,rep] <- sapply(resStats[[rep]]$res, function(res) {
    return(res$BlockStat)
  })
  for(radius in possible_radius) {
    approxStatisticGibbsPF[radius+1,rep] <- resStats[[rep]]$resGibbs[[radius+1]]$GibbsPFStat
  }
}
trueStatistics <- rep(computeTrueStatisticFilter(X_data = resStats[[1]]$X_data, n=n,
                                                 type = type_statistic, comp = comp_statistic), length(possible_cardinalities))

########################
# COMPUTE BIAS/MSE/VAR #
########################

dfRes <- computeDfBiasVar(approxStatisticPF = approxStatisticBlock, trueStatistics = trueStatistics,
                          algorithmName = "Block", dependentVarName = "Cardinality")
dfResGibbs <- computeDfBiasVar(approxStatisticPF = approxStatisticGibbsPF, trueStatistics = trueStatistics[1:nrow(approxStatisticGibbsPF)],
                          algorithmName = "GibbsPF", dependentVarName = "Cardinality")

########
# PLOT #
########

dfToPlot <- dfRes

#Bias
ggplot(dfToPlot[dfToPlot$Type == "Bias",], aes(x = Cardinality, y=Value, color = Algorithm)) +
  geom_point(size = 1.2) +  scale_colour_brewer(palette = "Set1") +
  geom_smooth(method="loess",se=TRUE) + geom_hline(yintercept = 0, linetype="dashed")

#Variance
ggplot(dfToPlot[dfToPlot$Type == "Var",], aes(x = Cardinality, y=Value, color = Algorithm)) +
  geom_point(size = 1.2) +  scale_colour_brewer(palette = "Set1") +
  geom_smooth(method="loess",se=TRUE)

