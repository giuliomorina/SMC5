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
id <- sample(1e7,size = 1)
set.seed(878,"L'Ecuyer-CMRG")
n <- 10 #Filter time (NOTE: the experiment is still performed in one time step)
dimension <- 8
ncores <- 4
repetitions <- 1000
N <- 20 #Number of particles
A <- generateA(c(0,0,0,0.6), dimension)
sigmaX <- diag(rep(1,dimension))
sigmaY <- diag(rep(1,dimension))
X0 <- rep(10, dimension)
possible_cardinalities <- 1:dimension
possible_radius <- (0:floor(dimension/2))
m <- 10 #Number of sweeps
type_statistic_plot <- "sum_squared" #Which statistic to plot
comp_statistic_plot <- 6 #Which component of the statistic (must be 1 for sum and sum squared)
if(type_statistic_plot == "sum" || type_statistic_plot == "sum_squared") {
  comp_statistic_plot <- 1
}

###################
# GENERATE X DATA #
###################

dataset <- generateData(n, A, X0, sigmaX, sigmaY)
X_data <- matrix(dataset$X_data, ncol=dimension)
fParams <- list(A=A, X0=X0, sigmaX=sigmaX)

##############################
# RUN BLOCK/GIBBS/SIR/KALMAN #
##############################

expRes <- mclapply(1:repetitions, function(rep) {
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
  resBlock <- lapply(possible_cardinalities, function(card) {
    blocks <- split(1:dimension,ceiling((1:dimension)/card))
    #Run Block
    BlockRes <- blockParticleFilterOnline(N=N, n=n, blocks = blocks,
                                          particles = kalmanParticles[1:(n-1),,],
                                          logWeights = array(log(1/N),c(n-1,length(blocks),N)),
                                          fParams = fParams, gParams = gParams)


    return(list(BlockRes = BlockRes, blocks = blocks))
  })
  resGibbs <- lapply(possible_radius, function(radius) {
    #Rub Gibbs PF
    GibbsPFRes <- gibbsParticleFilterOnline(N=N, n=n, m=m, radius=radius,
                                            particles = kalmanParticles[1:(n-1),,],
                                            fParams = fParams, gParams=gParams)


    return(list(GibbsPFRes = GibbsPFRes))
  })
  resSIR <- bootstrapParticleFilterOnline(N=N, n=n, fParams = fParams,
                                               particles = kalmanParticles[1:(n-1),,],
                                               logWeights = array(log(1/N),c(n-1,1,N)),
                                         gParams = gParams)
  resGlobalGibbs <- gibbsParticleFilterOnline(N=N, n=n, m=m, radius=dimension,
                                              particles = kalmanParticles[1:(n-1),,],
                                              fParams = fParams, gParams=gParams)
  resKalman <- list(filteringParticle = kalmanParticles,
                    filteringLogWeights = array(log(1/N),c(n,1,N)))
  return(list(gParams=gParams,
              fParams = fParams,
              X_data = X_data,
              resBlock=resBlock,
              resGibbs = resGibbs,
              resKalman=resKalman,
              resGlobalGibbs=resGlobalGibbs,
              resSIR = resSIR))
}, mc.cores = ncores)

#########################
# APPROXIMATE STATISTIC #
#########################

dfResList <- mclapply(c("mean","sum","mean_squared","sum_squared"), function(type_statistic){
  if(type_statistic == "sum" || type_statistic == "sum_squared") {
    possible_comp_statistic <- NA
  } else {
    possible_comp_statistic <- 1:dimension
  }
  dfResComp <- lapply(possible_comp_statistic, function(comp_statistic) {
    approxStatisticBlock <- matrix(NA, nrow=length(possible_cardinalities), ncol=repetitions)
    approxStatisticGibbs <- matrix(NA, nrow=length(possible_radius), ncol=repetitions)
    approxStatisticGlobalGibbs <- matrix(NA, nrow=1, ncol=repetitions)
    approxStatisticGlobalGibbs[1,] <- sapply(expRes, function(res) {
      computeApproxStatisticFilter(particles = res$resGlobalGibbs$filteringParticle,
                                   logWeights = res$resGlobalGibbs$filteringLogWeights,
                                   n = n,
                                   type = type_statistic,
                                   comp = comp_statistic)
    })
    approxStatisticKalman <- matrix(NA, nrow=1, ncol=repetitions)
    approxStatisticKalman[1,] <- sapply(expRes, function(res) {
      computeApproxStatisticFilter(particles = res$resKalman$filteringParticle,
                                   logWeights = res$resKalman$filteringLogWeights,
                                   n = n,
                                   type = type_statistic,
                                   comp = comp_statistic)
    })
    approxStatisticSIR <- matrix(NA, nrow=1, ncol=repetitions)
    approxStatisticSIR[1,] <- sapply(expRes, function(res) {
      computeApproxStatisticFilter(particles = res$resSIR$filteringParticle,
                                   logWeights = res$resSIR$filteringLogWeights,
                                   n = n,
                                   type = type_statistic,
                                   comp = comp_statistic)
    })
    trueStatistics <- matrix(NA, nrow=1, ncol=repetitions)
    trueStatistics[1,] <- sapply(expRes, function(res) {
      return(computeKalmanStatisticFilter(fParams=res$fParams,
                                          gParams = res$gParams,
                                          n = n,
                                          type = type_statistic,
                                          comp = comp_statistic))
    })
    dataStatistics <- computeTrueStatisticFilter(X_data = expRes[[1]]$X_data,
                                                  n = n, type = type_statistic, comp = comp_statistic)
    dataStatistics <- matrix(dataStatistics, ncol=repetitions)
    counter <- 1
    for(rep in 1:repetitions) {
      approxStatisticBlock[,counter] <- sapply(expRes[[rep]]$resBlock, function(res) {
        return(computeApproxStatisticFilter(particles = res$BlockRes$filteringParticle,
                                            logWeights = res$BlockRes$filteringLogWeights,
                                            n = n,
                                            blocks = res$blocks,
                                            type = type_statistic,
                                            comp = comp_statistic))
      })
      approxStatisticGibbs[,counter] <- sapply(expRes[[rep]]$resGibbs, function(res) {
        return(computeApproxStatisticFilter(particles = res$GibbsPFRes$filteringParticle,
                                            logWeights = res$GibbsPFRes$filteringLogWeights,
                                            n = n,
                                            type = type_statistic,
                                            comp = comp_statistic))
      })
      counter <- counter + 1
    }
    return(list(Block=approxStatisticBlock,
                GibbsPF=approxStatisticGlobalGibbs,
                LocGibbsPF=approxStatisticGibbs,
                Kalman=approxStatisticKalman,
                SIR=approxStatisticSIR,
                trueStatistics=trueStatistics,
                dataStatistics=dataStatistics))
  })
  return(dfResComp)
}, mc.cores = ncores)
names(dfResList) <- c("mean","sum","mean_squared","sum_squared")

#dfResList[[TYPE]][[COMPONENT]]$Block[CARDINALITY,REPETITION]
#dfResList[[TYPE]][[COMPONENT]]$GibbsPF[RADIUS,REPETITION]
#dfResList[[TYPE]][[COMPONENT]]$Kalman
#dfResList[[TYPE]][[COMPONENT]]$trueStatistics[REPETITION]

if(Sys.info()["nodename"] == "greyplover.stats.ox.ac.uk" || Sys.info()["nodename"] == "greypartridge.stats.ox.ac.uk" ||
   Sys.info()["nodename"] == "greyheron.stats.ox.ac.uk" || Sys.info()["nodename"] == "greywagtail.stats.ox.ac.uk") {
  saveRDS(dfResList, file = paste0("induced_bias_res_",id,".RDS"))
}

########################
# COMPUTE BIAS/MSE/VAR #
########################

dfResBiasVar <- mclapply(c("mean","sum","mean_squared","sum_squared"), function(type_statistic){
  if(type_statistic == "sum" || type_statistic == "sum_squared") {
    possible_comp_statistic <- 1
  } else {
    possible_comp_statistic <- 1:dimension
  }
  dfResBiasVarComp <- lapply(possible_comp_statistic, function(comp_statistic) {
    dfResVar <- as.data.frame(rbindlist(list(computeDfVar(approxList = dfResList[[type_statistic]][[comp_statistic]][c("Block")],
                                                          trueStatistics = dfResList[[type_statistic]][[comp_statistic]]$trueStatistics[rep(1:nrow(dfResList[[type_statistic]][[comp_statistic]]$trueStatistics), times = length(possible_cardinalities)),],
                                                          dependentVar = possible_cardinalities,
                                                          dependentVarName = "CardinalityRadius"),
                                             computeDfVar(approxList = dfResList[[type_statistic]][[comp_statistic]][c("LocGibbsPF")],
                                                          trueStatistics = dfResList[[type_statistic]][[comp_statistic]]$trueStatistics[rep(1:nrow(dfResList[[type_statistic]][[comp_statistic]]$trueStatistics), times = length(possible_radius)),],
                                                          dependentVar = pmin((possible_radius*2+1),rep(dimension,length(possible_radius))), #Radius 0 corresponds to cardinality 1, radius 1 corresponds to cardinality 3 and so on
                                                          dependentVarName = "CardinalityRadius")
    )))
    dfResVar <- as.data.frame(rbindlist(lapply(c(0,possible_cardinalities), function(card) {
      if(card == 0) return(dfResVar)
      computeDfVar(approxList = dfResList[[type_statistic]][[comp_statistic]][c("Kalman","SIR","GibbsPF")],
                   trueStatistics = dfResList[[type_statistic]][[comp_statistic]]$trueStatistics,
                   dependentVar = card,
                   dependentVarName = "CardinalityRadius")
    })))

    dfResBiasMSETrue <- as.data.frame(rbindlist(list(computeDfBiasMSE(approxList = dfResList[[type_statistic]][[comp_statistic]][c("Block")],
                                                                      trueStatistics = dfResList[[type_statistic]][[comp_statistic]]$trueStatistics[rep(1:nrow(dfResList[[type_statistic]][[comp_statistic]]$trueStatistics), times = length(possible_cardinalities)),],
                                                                      dependentVar = possible_cardinalities,
                                                                      dependentVarName = "CardinalityRadius"),
                                                     computeDfBiasMSE(approxList = dfResList[[type_statistic]][[comp_statistic]][c("LocGibbsPF")],
                                                                      trueStatistics = dfResList[[type_statistic]][[comp_statistic]]$trueStatistics[rep(1:nrow(dfResList[[type_statistic]][[comp_statistic]]$trueStatistics), times = length(possible_radius)),],
                                                                      dependentVar = pmin((possible_radius*2+1),rep(dimension,length(possible_radius))),
                                                                      dependentVarName = "CardinalityRadius")
    )))
    dfResBiasMSETrue <- as.data.frame(rbindlist(lapply(c(0,possible_cardinalities), function(card) {
      if(card == 0) return(dfResBiasMSETrue)
      computeDfBiasMSE(approxList = dfResList[[type_statistic]][[comp_statistic]][c("Kalman","SIR","GibbsPF")],
                       trueStatistics = dfResList[[type_statistic]][[comp_statistic]]$trueStatistics,
                       dependentVar = card,
                       dependentVarName = "CardinalityRadius")
    })))

    dfResBiasMSEData <- as.data.frame(rbindlist(list(computeDfBiasMSE(approxList = dfResList[[type_statistic]][[comp_statistic]][c("Block")],
                                                                      trueStatistics = dfResList[[type_statistic]][[comp_statistic]]$dataStatistics[rep(1:nrow(dfResList[[type_statistic]][[comp_statistic]]$dataStatistics), times = length(possible_cardinalities)),],
                                                                      dependentVar = possible_cardinalities,
                                                                      dependentVarName = "CardinalityRadius"),
                                                     computeDfBiasMSE(approxList = dfResList[[type_statistic]][[comp_statistic]][c("LocGibbsPF")],
                                                                      trueStatistics = dfResList[[type_statistic]][[comp_statistic]]$dataStatistics[rep(1:nrow(dfResList[[type_statistic]][[comp_statistic]]$dataStatistics), times = length(possible_radius)),],
                                                                      dependentVar = pmin((possible_radius*2+1),rep(dimension,length(possible_radius))),
                                                                      dependentVarName = "CardinalityRadius")
    )))
    dfResBiasMSEData <- as.data.frame(rbindlist(lapply(c(0,possible_cardinalities), function(card) {
      if(card == 0) return(dfResBiasMSEData)
      computeDfBiasMSE(approxList = dfResList[[type_statistic]][[comp_statistic]][c("Kalman","SIR","GibbsPF")],
                       trueStatistics = dfResList[[type_statistic]][[comp_statistic]]$dataStatistics,
                       dependentVar = card,
                       dependentVarName = "CardinalityRadius")
    })))
    return(list(dfResVar=dfResVar,
                dfResBiasMSEData=dfResBiasMSEData,
                dfResBiasMSETrue=dfResBiasMSETrue))
  })
  return(dfResBiasVarComp)
}, mc.cores = ncores)
names(dfResBiasVar) <- c("mean","sum","mean_squared","sum_squared")

if(Sys.info()["nodename"] == "greyplover.stats.ox.ac.uk" || Sys.info()["nodename"] == "greypartridge.stats.ox.ac.uk" ||
   Sys.info()["nodename"] == "greyheron.stats.ox.ac.uk" || Sys.info()["nodename"] == "greywagtail.stats.ox.ac.uk") {
  saveRDS(dfResBiasVar, file = paste0("induced_bias_res_biasvar_",id,".RDS"))
}

#########
# PLOTS #
#########

dfVarToPlot <- dfResBiasVar[[type_statistic_plot]][[comp_statistic_plot]]$dfResVar
dfBiasMSEToPlot <- dfResBiasVar[[type_statistic_plot]][[comp_statistic_plot]]$dfResBiasMSEData

#Variance
ggplot(dfVarToPlot[dfVarToPlot$Type == "Var",], aes(x = CardinalityRadius, y=Value, color = Algorithm)) +
  geom_point(size = 3, shape = 4) +  scale_colour_brewer(palette = "Dark2") +
  geom_smooth(method="lm",se=FALSE, size = 1.5) +
  ylab("Variance")

#Bias
ggplot(dfBiasMSEToPlot[dfBiasMSEToPlot$Type == "Bias" & (dfBiasMSEToPlot$Algorithm == "Block" |
                         dfBiasMSEToPlot$Algorithm == "LocGibbsPF"),],
        aes(x = CardinalityRadius, y=Value, color = Algorithm)) +
  scale_color_manual(values=c("Block" = "#e62d89", "LocGibbsPF" = "#66a423")) +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) as.numeric(quantile(x, probs=0.25)),
               fun.ymax = function(x) as.numeric(quantile(x, probs=0.75)),
               geom = "pointrange") +
  #geom_smooth(method="loess",se=TRUE) + geom_hline(yintercept = 0, linetype="dashed") +
  stat_summary(fun.y = mean,
               geom = "line") +
  theme_bw() +
  #geom_smooth(method="lm",se=TRUE) + geom_hline(yintercept = 0, linetype="dashed")
  ylab("Bias") + xlab("Cardinality or 2*Radius+1")

####################################
# OTHER PLOTS (INHOMOGENEOUS BIAS) #
####################################

#For a fixed cardinality and a fixed radius, the bias is not homogeneous in the components
# -> need to show the relative bias

fixed_card <- 8
fixed_radius <- 2
fixed_type_statistic <- "mean_squared"
dfResBiasMSEComponentData <- data.frame(Repetition = numeric(),
                                        DependentVar = numeric(),
                                        Algorithm = character(),
                                        Type = character(),
                                        Value= numeric(),
                                        stringsAsFactors = FALSE)
colnames(dfResBiasMSEComponentData) <- c("Repetition", "Component", "Algorithm", "Type", "Value")
for(comp in 1:dimension) {
  dfResListFixed <- list(dfResList[[fixed_type_statistic]][[comp]][c("Block")][[1]][fixed_card,],
                    dfResList[[fixed_type_statistic]][[comp]][c("LocGibbsPF")][[1]][fixed_radius,])
  names(dfResListFixed) <- c("Block","LocGibbsPF")

  dfResBiasMSEComponentData <- as.data.frame(rbindlist(list(dfResBiasMSEComponentData,
                                                computeDfBiasMSE(approxList = dfResListFixed,
                                                trueStatistics = dfResList[[fixed_type_statistic]][[comp]]$dataStatistics,
                                                dependentVar = comp,
                                                dependentVarName = "Component"))))
}



dfBiasMSEComponentToPlot <- dfResBiasMSEComponentData
ggplot(dfBiasMSEComponentToPlot[dfBiasMSEComponentToPlot$Type == "RelAbsBias",], aes(x = Component, y=Value, color = Algorithm)) +
  scale_colour_brewer(palette = "Set1") +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x),
               fun.ymax = function(x) mean(x) + sd(x),
               geom = "pointrange") +
  #geom_smooth(method="loess",se=TRUE) + geom_hline(yintercept = 0, linetype="dashed") +
  stat_summary(fun.y = mean,
               geom = "line") +
  ylab("Rel. Abs. Bias")

