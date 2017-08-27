# CURSE OF DIMENSIONALITY
# The curse of dimensionality appears as the dimension of the space grows as, again, all the weight
# is put on few particles leading to a big variance of the estimator. This happens also at time step n=1.
# This happen also when the matrix A is diagonal, i.e. there is no actual correlation between the components
# Blocked PF and Gibbs PF do not suffer from this.

rm(list=ls())
require(SMC5)
require(parallel)
require(ggplot2)

######################
# SETTING EXPERIMENT #
######################
set.seed(88,"L'Ecuyer-CMRG")
n <- 1 #Filter time
possible_dimension <- 3:25 #Possible dimensions. Need to be contiguous and start from 1!
ncores <- 5
repetitions <- 100
N <- 20 #Number of particles
A_diag <- 1 #Elements on the diagonal of A
varX <- 1 #Variance of sigmaX
varY <- 1 #Variance of sigmaY
X0_point <- 10 #Initial point, just away from 0 to avoid problem when estimate relative values
card_block <- 3 #Cardinality of the block
radius <- 1 #Radius of Gibbs PF
m <- 10 #Number of sweeps

#######################
# RUN SIR/BLOCK/GIBBS #
#######################

expRes <- lapply(possible_dimension, function(dimension) {
  #GENERATE MODEL
  A <- diag(rep(A_diag,dimension))
  X0 <- rep(X0_point, dimension)
  sigmaX <- diag(rep(varX, dimension))
  sigmaY <- diag(rep(varY, dimension))
  dataset <- generateData(n, A, X0, sigmaX, sigmaY)
  X_data <- matrix(dataset$X_data, ncol=dimension)
  fParams <- list(A=A, X0=X0, sigmaX=sigmaX)
  #RUN EXPERIMENT
  res <- mclapply(1:repetitions, function(rep) {
    #GENERATE NEW Y, different for every repetition, always with the same X
    Y_data <- generateY(X_data, sigmaY)
    gParams <- list(Y=Y_data, sigmaY=sigmaY)
    blocks <- split(1:dimension,ceiling((1:dimension)/card_block))
    #Kalman
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
    return(list(SIRRes = bootstrapParticleFilter(N=N, n=n, fParams = fParams,
                                                 gParams = gParams),
                BlockRes = blockParticleFilter(N=N, n=n, blocks = blocks,
                                               fParams = fParams, gParams = gParams),
                GibbsRes = gibbsParticleFilter(N=N, n=n, m=m, radius=radius,
                                               fParams = fParams, gParams=gParams),
                KalmanRes = list(filteringParticle = kalmanParticles,
                                 filteringLogWeights = array(log(1/N),c(n,1,N))),
                gParams = gParams,
                Y_data = Y_data,
                blocks = blocks))

  }, mc.cores = ncores)
  return(list(X_data = X_data,
              fParams = fParams,
              res = res,
              dimension = dimension))
})

#########################
# APPROXIMATE STATISTIC #
#########################

if(Sys.info()["nodename"] == "greyplover.stats.ox.ac.uk" || Sys.info()["nodename"] == "greypartridge.stats.ox.ac.uk" ||
   Sys.info()["nodename"] == "greyheron.stats.ox.ac.uk" || Sys.info()["nodename"] == "greywagtail.stats.ox.ac.uk") {
  #If on the server, compute all the statistics!
  type_statistic <- c("mean", "sum", "mean_squared", "sum_squared")
}

dfResList <- lapply(c("mean","sum","mean_squared","sum_squared"), function(type_statistic){
  approxStatisticSIR <- matrix(NA, nrow=length(possible_dimension), ncol=repetitions)
  approxStatisticBlock <- matrix(NA, nrow=length(possible_dimension), ncol=repetitions)
  approxStatisticGibbsPF <- matrix(NA, nrow=length(possible_dimension), ncol=repetitions)
  approxStatisticKalman <- matrix(NA, nrow=length(possible_dimension), ncol=repetitions) #Baseline

  dataStatistics <- rep(NA, length(possible_dimension))
  trueStatistics <- matrix(NA, nrow=length(possible_dimension), ncol=repetitions) #Using Kalman E[theta|Y]
  counter <- 1
  for(dimension in possible_dimension) {
    # if(expRes[[dimension]]$dimension != dimension) {
    #   stop("Uh-oh something's wrong! Perhaps possible_dimension vector is not continuous and starting from one?")
    # }
    approxStatisticSIR[counter,] <- sapply(expRes[[counter]]$res, function(res) {
      return(computeApproxStatisticFilter(particles = res$SIRRes$filteringParticle,
                                          logWeights = res$SIRRes$filteringLogWeights,
                                          n = n,
                                          type = type_statistic,
                                          comp = comp_statistic))
    })
    approxStatisticBlock[counter,] <- sapply(expRes[[counter]]$res, function(res) {
      return(computeApproxStatisticFilter(particles = res$BlockRes$filteringParticle,
                                          logWeights = res$BlockRes$filteringLogWeights,
                                          n = n,
                                          blocks = res$blocks,
                                          type = type_statistic,
                                          comp = comp_statistic))
    })
    approxStatisticGibbsPF[counter,] <- sapply(expRes[[counter]]$res, function(res) {
      return(computeApproxStatisticFilter(particles = res$GibbsRes$filteringParticle,
                                          logWeights = res$GibbsRes$filteringLogWeights,
                                          n = n,
                                          type = type_statistic,
                                          comp = comp_statistic))
    })
    approxStatisticKalman[counter,] <- sapply(expRes[[counter]]$res, function(res) {
      return(computeApproxStatisticFilter(particles = res$KalmanRes$filteringParticle,
                                          logWeights = res$KalmanRes$filteringLogWeights,
                                          n = n,
                                          type = type_statistic,
                                          comp = comp_statistic))
    })
    trueStatistics[counter,] <- sapply(expRes[[counter]]$res, function(res) {
      return(computeKalmanStatisticFilter(fParams=expRes[[counter]]$fParams,
                                          gParams = res$gParams,
                                          n = n,
                                          type = type_statistic,
                                          comp = comp_statistic))
    })
    dataStatistics[counter] <- computeTrueStatisticFilter(X_data = expRes[[counter]]$X_data,
                                                          n = n, type = type_statistic, comp = comp_statistic)
    counter <- counter+1
  }

  return(list(approxStatisticSIR=approxStatisticSIR,
              approxStatisticBlock=approxStatisticBlock,
              approxStatisticGibbsPF=approxStatisticGibbsPF,
              approxStatisticKalman=approxStatisticKalman,
              dataStatistics=dataStatistics,
              trueStatistics=trueStatistics))
})


if(Sys.info()["nodename"] == "greyplover.stats.ox.ac.uk" || Sys.info()["nodename"] == "greypartridge.stats.ox.ac.uk" ||
   Sys.info()["nodename"] == "greyheron.stats.ox.ac.uk" || Sys.info()["nodename"] == "greywagtail.stats.ox.ac.uk") {
  saveRDS(dfResList, file = paste0("curse_of_dimensionality_res_",sample(1e7,size = 1),".RDS"))
}

########################
# COMPUTE BIAS/MSE/VAR #
########################

dfRes <- computeDfBiasVar(approxStatisticPF = approxStatisticSIR, trueStatistics = trueStatistics,
                          algorithmName = "SIR", dependentVarName = "Dimension")
dfRes <- computeDfBiasVar(approxStatisticPF = approxStatisticBlock, trueStatistics = trueStatistics,
                          algorithmName = "Block", dependentVarName = "Dimension", dfRes = dfRes)
dfRes <- computeDfBiasVar(approxStatisticPF = approxStatisticGibbsPF, trueStatistics = trueStatistics,
                          algorithmName = "GibbsPF", dependentVarName = "Dimension", dfRes = dfRes)
dfRes <- computeDfBiasVar(approxStatisticPF = approxStatisticKalman, trueStatistics = trueStatistics,
                          algorithmName = "Kalman", dependentVarName = "Dimension", dfRes = dfRes)
dfRes$Dimension <- dfRes$Dimension + (possible_dimension[1] - 1)
dfRes2 <- computeDfBiasVar(approxStatisticPF = approxStatisticSIR, trueStatistics = trueStatistics2,
                           algorithmName = "SIR", dependentVarName = "Dimension")
dfRes2 <- computeDfBiasVar(approxStatisticPF = approxStatisticBlock, trueStatistics = trueStatistics2,
                           algorithmName = "Block", dependentVarName = "Dimension", dfRes = dfRes2)
dfRes2 <- computeDfBiasVar(approxStatisticPF = approxStatisticGibbsPF, trueStatistics = trueStatistics2,
                           algorithmName = "GibbsPF", dependentVarName = "Dimension", dfRes = dfRes2)
dfRes2 <- computeDfBiasVar(approxStatisticPF = approxStatisticKalman, trueStatistics = trueStatistics2,
                           algorithmName = "Kalman", dependentVarName = "Dimension", dfRes = dfRes2)
dfRes2$Dimension <- dfRes2$Dimension + (possible_dimension[1] - 1)

########
# PLOT #
########

dfToPlot <- dfResList[[1]]$dfRes2

#Variance
#Just to check the different type of relative variances:
ggplot(dfToPlot[dfToPlot$Type != "Var" & dfToPlot$Type != "MSE" &
                  dfToPlot$Type != "Bias" & dfToPlot$Type != "RMSE" &
                  dfToPlot$Type != "RelAbsBias",], aes(x = as.numeric(Dimension), y=as.numeric(Value), group=interaction(Algorithm,Type),
                                                       linetype = Algorithm, color = Type)) + geom_line()
#RelVar
ggplot(dfToPlot[dfToPlot$Type == "RelVar",], aes(x = Dimension, y=Value, color = Algorithm)) +
  geom_line(size = 1.2) +  scale_colour_brewer(palette = "Set1")

ggplot(dfToPlot[dfToPlot$Type == "RelVar",], aes(x = Dimension, y=Value, color = Algorithm)) +
  geom_point(size = 1.2) +  scale_colour_brewer(palette = "Set1") +
  geom_smooth(method="loess",se=FALSE)

#Bias
ggplot(dfToPlot[dfToPlot$Type == "Bias",], aes(x = as.numeric(Dimension), y=abs(Value), color = Algorithm)) + geom_line() +
  scale_colour_brewer(palette = "Set1")

ggplot(dfToPlot[dfToPlot$Type == "Bias",], aes(x = Dimension, y=Value, color = Algorithm)) +
  geom_point(size = 1.2) +  scale_colour_brewer(palette = "Set1") +
  geom_smooth(method="loess",se=TRUE) + geom_hline(yintercept = 0, linetype="dashed")

PROVA <- data.frame(Repetition = numeric(),
                    Dimension = numeric(),
                    Bias = numeric(),
                    stringsAsFactors = FALSE)
for(i in 1:nrow(dfResList[[1]]$approxStatisticSIR)) {
  for(j in 1:ncol(dfResList[[1]]$approxStatisticSIR)) {
    PROVA <- rbind(PROVA,c(Repetition = j,
                           Dimension = possible_dimension[i],
                           Bias = dfResList[[1]]$approxStatisticSIR[i,j] - dfResList[[1]]$trueStatistics2[i,j]),
                   stringsAsFactors = FALSE)
  }
}
colnames(PROVA) <- c("Repetition","Dimension","Bias")
ggplot(PROVA, aes(x=Dimension, y=Bias)) +  scale_colour_brewer(palette = "Set1") +
  geom_smooth(method="lm" ,se=TRUE)

PROVA2 <- computeDfBiasMSE(approxList = dfResList[[1]][3:6], trueStatistics = dfResList[[1]]$trueStatistics2,
                           dependentVar = possible_dimension, dependentVarName = "Dimension")
ggplot(PROVA2[PROVA2$Type=="Bias",], aes(x=Dimension, y=Value, color = Algorithm)) +  scale_colour_brewer(palette = "Set1") +
  geom_smooth(method="lm" ,se=TRUE)
