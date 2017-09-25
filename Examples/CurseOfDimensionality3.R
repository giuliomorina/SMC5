# CURSE OF DIMENSIONALITY 3
# In this example we want to show that the curse of dimensionality is basically
# intrinsic in the correction operator when it is applied to an empirical measure.
# To do so we start by sampling correctly from C_1P\pi_0 (Kalman Filter) and then we
# sample from the empirical measure.
# Sampling from the empirical measure is equivalent as sampling from a mixture of
# multivariate gaussian (all with equal probability) whose parameters are
# mean = 1/2*I*(A*X_n^{(i)}+Y_n)
# cov. = 1/2*I
# In turns this implies that the Gibbs PF suffers from the curse of dimensionality
# as it converges to this setting!

rm(list=ls())
require(SMC5)
require(parallel)
require(ggplot2)
require(MASS)
######################
# SETTING EXPERIMENT #
######################
id <- sample(1e7,size = 1)
set.seed(8877,"L'Ecuyer-CMRG")
n <- 1 #Filter time
possible_dimension <- seq(from=1, to=50, by=1)
N <- 1000 #Number of particles
ncores <- 4
repetitions <- 5
A_diag <- 1 #Elements on the diagonal of A
varX <- 1 #Variance of sigmaX
varY <- 1 #Variance of sigmaY
X0_point <- 10 #Initial point, just away from 0 to avoid problem when estimate relative values

type_statistic_plot <- "sum_squared" #Which statistic to plot
comp_statistic_plot <- 2 #Which component of the statistic (must be 1 for sum and sum squared)
if(type_statistic_plot == "sum" || type_statistic_plot == "sum_squared") {
  comp_statistic_plot <- 1
}

##################
# RUN EXPERIMENT #
##################

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
    #Optimal SIR
    optimalSIRParticles <- array(NA,c(n,dimension,N))
    #Fill the particles at time 1 with the particles from Kalman
    covMatrix <- 0.5*diag(dimension)
    optimalSIRParticles[1,,] <- t(mvrnormArma(N,mean = kalmanFilterRes$m[,aux+1], sigma = as.matrix(kalmanFilterRes$C[,,aux+1])))
    #Fill the remaining times
    if(n > 1) {
      for(k in 2:n) {
        for(i in 1:N) {
          #Choose a particle at random
          J <- sample(N, size=1)
          #Sample from multivariate gaussian
          optimalSIRParticles[k,,i] <- mvrnorm(n=1, mu = covMatrix%*%(A%*%optimalSIRParticles[k-1,,J]+Y_data[k]),
                                               Sigma = covMatrix)
        }
      }
    }

    return(list(OptimalSIR = list(filteringParticle = optimalSIRParticles,
                                  filteringLogWeights = array(log(1/N),c(n,1,N))),
                KalmanRes = list(filteringParticle = kalmanParticles,
                                 filteringLogWeights = array(log(1/N),c(n,1,N))),
                gParams = gParams,
                Y_data = Y_data))

  }, mc.cores = ncores)
  return(list(X_data = X_data,
              fParams = fParams,
              res = res,
              dimension = dimension))
})

#########################
# APPROXIMATE STATISTIC #
#########################

#dfResList <- lapply(c("mean","sum","mean_squared","sum_squared"), function(type_statistic){
dfResList <- lapply(c("sum","sum_squared"), function(type_statistic){
  if(type_statistic == "sum" || type_statistic == "sum_squared") {
    possible_comp_statistic <- NA
  } else {
    possible_comp_statistic <- 1:min(possible_dimension)
  }
  dfResComp <- lapply(possible_comp_statistic, function(comp_statistic) {
    approxStatisticOptimalSIR <- matrix(NA, nrow=length(possible_dimension), ncol=repetitions)
    approxStatisticKalman <- matrix(NA, nrow=length(possible_dimension), ncol=repetitions) #Baseline
    dataStatistics <- rep(NA, length(possible_dimension))
    trueStatistics <- matrix(NA, nrow=length(possible_dimension), ncol=repetitions) #Using Kalman E[theta|Y]

    counter <- 1
    for(dimension in possible_dimension) {
      approxStatisticOptimalSIR[counter,] <- sapply(expRes[[counter]]$res, function(res) {
        return(computeApproxStatisticFilter(particles = res$OptimalSIR$filteringParticle,
                                            logWeights = res$OptimalSIR$filteringLogWeights,
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

    return(list(OptimalSIR=approxStatisticOptimalSIR,
                Kalman=approxStatisticKalman,
                dataStatistics=dataStatistics,
                trueStatistics=trueStatistics))
  })
  return(dfResComp)
})
#names(dfResList) <- c("mean","sum","mean_squared","sum_squared")
names(dfResList) <- c("sum","sum_squared")

########################
# COMPUTE BIAS/MSE/VAR #
########################

#dfResBiasVar <- lapply(c("mean","sum","mean_squared","sum_squared"), function(type_statistic){
dfResBiasVar <- lapply(c("sum","sum_squared"), function(type_statistic){
  if(type_statistic == "sum" || type_statistic == "sum_squared") {
    possible_comp_statistic <- 1
  } else {
    possible_comp_statistic <- 1:min(possible_dimension)
  }
  dfResBiasVarComp <- lapply(possible_comp_statistic, function(comp_statistic) {

    dfResVar <- computeDfVar(approxList = dfResList[[type_statistic]][[comp_statistic]][c("OptimalSIR","Kalman")],
                             trueStatistics = dfResList[[type_statistic]][[comp_statistic]]$dataStatistics,
                             dependentVar = possible_dimension,
                             dependentVarName = "Dimension")

    dfResBiasMSEData <- computeDfBiasMSE(approxList = dfResList[[type_statistic]][[comp_statistic]][c("OptimalSIR","Kalman")],
                                         trueStatistics = dfResList[[type_statistic]][[comp_statistic]]$dataStatistics,
                                         dependentVar = possible_dimension,
                                         dependentVarName = "Dimension")
    dfResBiasMSETrue <- computeDfBiasMSE(approxList = dfResList[[type_statistic]][[comp_statistic]][c("OptimalSIR","Kalman")],
                                         trueStatistics = dfResList[[type_statistic]][[comp_statistic]]$trueStatistics,
                                         dependentVar = possible_dimension,
                                         dependentVarName = "Dimension")

    return(list(dfResVar=dfResVar,
                dfResBiasMSEData=dfResBiasMSEData,
                dfResBiasMSETrue=dfResBiasMSETrue))
  })
  return(dfResBiasVarComp)
})
names(dfResBiasVar) <- c("sum","sum_squared")

########
# PLOT #
########

dfResVarToPlot <-  dfResBiasVar[[type_statistic_plot]][[comp_statistic_plot]]$dfResVar
dfResVarToPlot$Algorithm <- factor(dfResVarToPlot$Algorithm, c("OptimalSIR","Kalman"))
ggplot(dfResVarToPlot[dfResVarToPlot$Type == "RelVar",], aes(x = Dimension, y=Value, color = Algorithm)) +
  geom_point(size = 3, shape = 4) +
  scale_color_manual(values=c("Kalman" = "#e5a900", "OptimalSIR" = "#7470b2")) +
  theme_grey(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Relative Variance") +
  theme_bw() +
  geom_smooth(method="lm",se=FALSE, size = 1.5)
