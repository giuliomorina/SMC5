# GIBBS PF
# This example wants to show that both Gibbs PF induce
# an additional bias that increases as the radius decreases.
# We run the experiment fixing n=10

rm(list=ls())
require(SMC5)
require(parallel)
require(ggplot2)

######################
# SETTING EXPERIMENT #
######################
id <- sample(1e7,size = 1)
set.seed(87817,"L'Ecuyer-CMRG")
n <- 10 #Filter time
dimension <- 60 #Dimension of the state space
N <- 200 #Number of particles
ncores <- 20 #Number of cores used
repetitions <- 1000 #Repetitions of the experiment to estimate bias/var
alpha <- -1 #Power for the decay of correlation. Should be negative. The bigger it is, less correlation there is!
A <- generateA_power(alpha, dimension)
sigmaX <- diag(rep(1,dimension))
sigmaY <- diag(rep(1,dimension))
m_0 <- rep(0, dimension)
X0 <- mvrnormArma(1, mean = m_0, sigma = sigmaX) #Initial distribution
possible_radius <- seq(0,floor(dimension/2)+1, by=10)
m <- 100 #Number of sweeps
type_statistic_plot <- "sum" #Which statistic to plot
comp_statistic_plot <- 1 #Which component of the statistic (must be 1 for sum and sum squared)
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
  #Generate iid sample from Kalman filter
  kalmanFilterRes <- KalmanFilterCpp(m_0 = m_0,
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
  resKalman <- list(filteringParticle = kalmanParticles,
                    filteringLogWeights = array(log(1/N),c(n,1,N)))
  #Run local Gibbs for each radius
  resGibbs <- lapply(possible_radius, function(radius) {
    #Rub Gibbs PF
    GibbsPFRes <- gibbsParticleFilter(N=N, n=n, m=m, radius=radius,
                                      fParams = fParams, gParams=gParams)
    return(list(GibbsPFRes = GibbsPFRes, radius = radius))
  })

  outlist <- list(gParams=gParams,
                  fParams = fParams,
                  resGibbs = resGibbs,
                  resKalman=resKalman)

  if(Sys.info()["nodename"] == "greyplover.stats.ox.ac.uk" || Sys.info()["nodename"] == "greypartridge.stats.ox.ac.uk" ||
     Sys.info()["nodename"] == "greyheron.stats.ox.ac.uk" || Sys.info()["nodename"] == "greywagtail.stats.ox.ac.uk") {
    saveRDS(dfResList, file = paste0("GibbsPF_rep_",rep,"_",id,".RDS"))
  }

  return(outlist)
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
    approxStatisticGibbs <- matrix(NA, nrow=length(possible_radius), ncol=repetitions)
    approxStatisticKalman <- matrix(NA, nrow=1, ncol=repetitions)
    approxStatisticKalman[1,] <- sapply(expRes, function(res) {
      computeApproxStatisticFilter(particles = res$resKalman$filteringParticle,
                                   logWeights = res$resKalman$filteringLogWeights,
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
    dataStatistics <- computeTrueStatisticFilter(X_data = X_data,
                                                 n = n, type = type_statistic, comp = comp_statistic)
    dataStatistics <- matrix(dataStatistics, ncol=repetitions)
    counter <- 1
    for(rep in 1:repetitions) {
      approxStatisticGibbs[,counter] <- sapply(expRes[[rep]]$resGibbs, function(res) {
        return(computeApproxStatisticFilter(particles = res$GibbsPFRes$filteringParticle,
                                            logWeights = res$GibbsPFRes$filteringLogWeights,
                                            n = n,
                                            type = type_statistic,
                                            comp = comp_statistic))
      })
      counter <- counter + 1
    }
    return(list(LocGibbsPF=approxStatisticGibbs,
                Kalman=approxStatisticKalman,
                trueStatistics=trueStatistics,
                dataStatistics=dataStatistics))
  })
  return(dfResComp)
}, mc.cores = ncores)
names(dfResList) <- c("mean","sum","mean_squared","sum_squared")

if(Sys.info()["nodename"] == "greyplover.stats.ox.ac.uk" || Sys.info()["nodename"] == "greypartridge.stats.ox.ac.uk" ||
   Sys.info()["nodename"] == "greyheron.stats.ox.ac.uk" || Sys.info()["nodename"] == "greywagtail.stats.ox.ac.uk") {
  saveRDS(dfResList, file = paste0("GibbsPF_res_id_",id,".RDS"))
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
    dfResVar <- as.data.frame(rbindlist(list(computeDfVar(approxList = dfResList[[type_statistic]][[comp_statistic]][c("LocGibbsPF")],
                                                          trueStatistics = dfResList[[type_statistic]][[comp_statistic]]$trueStatistics[rep(1:nrow(dfResList[[type_statistic]][[comp_statistic]]$trueStatistics), times = length(possible_radius)),],
                                                          dependentVar = pmin((possible_radius*2+1),rep(dimension,length(possible_radius))), #Radius 0 corresponds to cardinality 1, radius 1 corresponds to cardinality 3 and so on
                                                          dependentVarName = "CardinalityRadius")
    )))
    dfResVar <- as.data.frame(rbindlist(lapply(c(0,possible_radius), function(card) {
      if(card == 0) return(dfResVar)
      computeDfVar(approxList = dfResList[[type_statistic]][[comp_statistic]][c("Kalman")],
                   trueStatistics = dfResList[[type_statistic]][[comp_statistic]]$trueStatistics,
                   dependentVar = card,
                   dependentVarName = "CardinalityRadius")
    })))

    dfResBiasMSETrue <- as.data.frame(rbindlist(list(computeDfBiasMSE(approxList = dfResList[[type_statistic]][[comp_statistic]][c("LocGibbsPF")],
                                                                      trueStatistics = dfResList[[type_statistic]][[comp_statistic]]$trueStatistics[rep(1:nrow(dfResList[[type_statistic]][[comp_statistic]]$trueStatistics), times = length(possible_radius)),],
                                                                      dependentVar = pmin((possible_radius*2+1),rep(dimension,length(possible_radius))),
                                                                      dependentVarName = "CardinalityRadius")
    )))
    dfResBiasMSETrue <- as.data.frame(rbindlist(lapply(c(0,possible_radius), function(card) {
      if(card == 0) return(dfResBiasMSETrue)
      computeDfBiasMSE(approxList = dfResList[[type_statistic]][[comp_statistic]][c("Kalman")],
                       trueStatistics = dfResList[[type_statistic]][[comp_statistic]]$trueStatistics,
                       dependentVar = card,
                       dependentVarName = "CardinalityRadius")
    })))

    dfResBiasMSEData <- as.data.frame(rbindlist(list(computeDfBiasMSE(approxList = dfResList[[type_statistic]][[comp_statistic]][c("LocGibbsPF")],
                                                                      trueStatistics = dfResList[[type_statistic]][[comp_statistic]]$dataStatistics[rep(1:nrow(dfResList[[type_statistic]][[comp_statistic]]$dataStatistics), times = length(possible_radius)),],
                                                                      dependentVar = pmin((possible_radius*2+1),rep(dimension,length(possible_radius))),
                                                                      dependentVarName = "CardinalityRadius")
    )))
    dfResBiasMSEData <- as.data.frame(rbindlist(lapply(c(0,possible_radius), function(card) {
      if(card == 0) return(dfResBiasMSEData)
      computeDfBiasMSE(approxList = dfResList[[type_statistic]][[comp_statistic]][c("Kalman")],
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
  saveRDS(dfResBiasVar, file = paste0("GibbsPF_res_biasvar_id_",id,".RDS"))
}

#########
# PLOTS #
#########

dfVarToPlot <- dfResBiasVar[[type_statistic_plot]][[comp_statistic_plot]]$dfResVar
dfBiasMSEToPlot <- dfResBiasVar[[type_statistic_plot]][[comp_statistic_plot]]$dfResBiasMSETrue

#Variance
ggplot(dfVarToPlot[dfVarToPlot$Type == "Var",], aes(x = CardinalityRadius, y=Value, color = Algorithm)) +
  geom_point(size = 3, shape = 4) +  scale_colour_brewer(palette = "Dark2") +
  geom_smooth(method="lm",se=FALSE, size = 1.5) +
  ylab("Variance")

#Bias
ggplot(dfBiasMSEToPlot[dfBiasMSEToPlot$Type == "Bias" & dfBiasMSEToPlot$Algorithm != "Block",],
       aes(x = CardinalityRadius, y=Value, color = Algorithm)) +
  scale_color_manual(values=c("Block" = "#e62d89", "LocGibbsPF" = "#66a423", "Kalman" = "#000000")) +
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

fixed_card <- 5
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

