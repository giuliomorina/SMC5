# CURSE OF DIMENSIONALITY 2

rm(list=ls())
require(SMC5)
require(parallel)
require(ggplot2)

######################
# SETTING EXPERIMENT #
######################
id <- sample(1e7,size = 1)
set.seed(88,"L'Ecuyer-CMRG")
n <- 10 #Filter time - Still done in 1 time step
possible_dimension <- seq(from=1, to=50, by=1) #Possible dimensions. Need to be contiguous and start from 1!
ncores <- 4
repetitions <- 10
N <- 20 #Number of particles
A_diag <- 1 #Elements on the diagonal of A
varX <- 1 #Variance of sigmaX
varY <- 1 #Variance of sigmaY
X0_point <- 10 #Initial point, just away from 0 to avoid problem when estimate relative values
card_block <- 1 #Cardinality of the block
radius <- 0 #Radius of Gibbs PF
m <- 10 #Number of sweeps

type_statistic_plot <- "sum_squared" #Which statistic to plot
comp_statistic_plot <- 2 #Which component of the statistic (must be 1 for sum and sum squared)
if(type_statistic_plot == "sum" || type_statistic_plot == "sum_squared") {
  comp_statistic_plot <- 1
}

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
    return(list(SISRes = sequentialImportanceSamplingOnline(N=N, n=n, fParams = fParams,
                                                            particles = array(kalmanParticles[1:(n-1),,], c(n-1,dimension,N)),
                                                            logWeights = array(log(1/N),c(n-1,1,N)),
                                                            gParams = gParams),
                SIRRes = bootstrapParticleFilterOnline(N=N, n=n, fParams = fParams,
                                                 particles = array(kalmanParticles[1:(n-1),,], c(n-1,dimension,N)),
                                                 logWeights = array(log(1/N),c(n-1,1,N)),
                                                 gParams = gParams),
                BlockRes = blockParticleFilterOnline(N=N, n=n, blocks = blocks,
                                               particles = array(kalmanParticles[1:(n-1),,], c(n-1,dimension,N)),
                                               logWeights = array(log(1/N),c(n-1,length(blocks),N)),
                                               fParams = fParams, gParams = gParams),
                GibbsRes = gibbsParticleFilterOnline(N=N, n=n, m=m, radius=radius,
                                               particles = array(kalmanParticles[1:(n-1),,], c(n-1,dimension,N)),
                                               fParams = fParams, gParams=gParams),
                GlobalGibbsRes = gibbsParticleFilterOnline(N=N, n=n, m=m, radius=dimension,
                                                     particles = array(kalmanParticles[1:(n-1),,], c(n-1,dimension,N)),
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

dfResList <- lapply(c("mean","sum","mean_squared","sum_squared"), function(type_statistic){
  if(type_statistic == "sum" || type_statistic == "sum_squared") {
    possible_comp_statistic <- NA
  } else {
    possible_comp_statistic <- 1:min(possible_dimension)
  }
  dfResComp <- lapply(possible_comp_statistic, function(comp_statistic) {
    approxStatisticSIS <- matrix(NA, nrow=length(possible_dimension), ncol=repetitions)
    approxStatisticSIR <- matrix(NA, nrow=length(possible_dimension), ncol=repetitions)
    approxStatisticBlock <- matrix(NA, nrow=length(possible_dimension), ncol=repetitions)
    approxStatisticGibbsPF <- matrix(NA, nrow=length(possible_dimension), ncol=repetitions)
    approxStatisticGlobalGibbsPF <- matrix(NA, nrow=length(possible_dimension), ncol=repetitions)
    approxStatisticKalman <- matrix(NA, nrow=length(possible_dimension), ncol=repetitions) #Baseline
    dataStatistics <- rep(NA, length(possible_dimension))
    trueStatistics <- matrix(NA, nrow=length(possible_dimension), ncol=repetitions) #Using Kalman E[theta|Y]

    counter <- 1
    for(dimension in possible_dimension) {
      approxStatisticSIS[counter,] <- sapply(expRes[[counter]]$res, function(res) {
        return(computeApproxStatisticFilter(particles = res$SISRes$filteringParticle,
                                            logWeights = res$SISRes$filteringLogWeights,
                                            n = n,
                                            type = type_statistic,
                                            comp = comp_statistic))
      })
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
      approxStatisticGlobalGibbsPF[counter,] <- sapply(expRes[[counter]]$res, function(res) {
        return(computeApproxStatisticFilter(particles = res$GlobalGibbsRes$filteringParticle,
                                            logWeights = res$GlobalGibbsRes$filteringLogWeights,
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

    return(list(SIS=approxStatisticSIS,
                SIR=approxStatisticSIR,
                Block=approxStatisticBlock,
                GibbsPF=approxStatisticGlobalGibbsPF,
                LocGibbsPF=approxStatisticGibbsPF,
                Kalman=approxStatisticKalman,
                dataStatistics=dataStatistics,
                trueStatistics=trueStatistics))
  })
  return(dfResComp)
})

names(dfResList) <- c("mean","sum","mean_squared","sum_squared")

if(Sys.info()["nodename"] == "greyplover.stats.ox.ac.uk" || Sys.info()["nodename"] == "greypartridge.stats.ox.ac.uk" ||
   Sys.info()["nodename"] == "greyheron.stats.ox.ac.uk" || Sys.info()["nodename"] == "greywagtail.stats.ox.ac.uk") {
  saveRDS(dfResList, file = paste0("curse_of_dimensionality2_res_",id,".RDS"))
}

########################
# COMPUTE BIAS/MSE/VAR #
########################

dfResBiasVar <- lapply(c("mean","sum","mean_squared","sum_squared"), function(type_statistic){
  if(type_statistic == "sum" || type_statistic == "sum_squared") {
    possible_comp_statistic <- 1
  } else {
    possible_comp_statistic <- 1:min(possible_dimension)
  }
  dfResBiasVarComp <- lapply(possible_comp_statistic, function(comp_statistic) {

    dfResVar <- computeDfVar(approxList = dfResList[[type_statistic]][[comp_statistic]][c("SIR","Block","GibbsPF","Kalman","SIS","LocGibbsPF")],
                             trueStatistics = dfResList[[type_statistic]][[comp_statistic]]$dataStatistics,
                             dependentVar = possible_dimension,
                             dependentVarName = "Dimension")

    dfResBiasMSEData <- computeDfBiasMSE(approxList = dfResList[[type_statistic]][[comp_statistic]][c("SIR","Block","GibbsPF","Kalman","SIS","LocGibbsPF")],
                                         trueStatistics = dfResList[[type_statistic]][[comp_statistic]]$dataStatistics,
                                         dependentVar = possible_dimension,
                                         dependentVarName = "Dimension")
    dfResBiasMSETrue <- computeDfBiasMSE(approxList = dfResList[[type_statistic]][[comp_statistic]][c("SIR","Block","GibbsPF","Kalman","SIS","LocGibbsPF")],
                                         trueStatistics = dfResList[[type_statistic]][[comp_statistic]]$trueStatistics,
                                         dependentVar = possible_dimension,
                                         dependentVarName = "Dimension")

    return(list(dfResVar=dfResVar,
                dfResBiasMSEData=dfResBiasMSEData,
                dfResBiasMSETrue=dfResBiasMSETrue))
  })
  return(dfResBiasVarComp)
})
names(dfResBiasVar) <- c("mean","sum","mean_squared","sum_squared")

if(Sys.info()["nodename"] == "greyplover.stats.ox.ac.uk" || Sys.info()["nodename"] == "greypartridge.stats.ox.ac.uk" ||
   Sys.info()["nodename"] == "greyheron.stats.ox.ac.uk" || Sys.info()["nodename"] == "greywagtail.stats.ox.ac.uk") {
  saveRDS(dfResBiasVar, file = paste0("curse_of_dimensionality2_res_biasvar_",id,".RDS"))
}


########################
# MEAN MAXIMUM WEIGHTS #
########################

maxWeightSIR <- matrix(NA, nrow=length(possible_dimension), ncol=repetitions)
maxWeightSIS <- matrix(NA, nrow=length(possible_dimension), ncol=repetitions)
maxWeightBlock <- matrix(NA, nrow=length(possible_dimension), ncol=repetitions)

for(dim in possible_dimension) {
  for(rep in 1:repetitions) {
    maxWeightSIR[dim,] <- sapply(expRes[[dim]]$res, function(res_rep) {
      max(exp(res_rep$SIRRes$filteringLogWeights[n,1,]-logAdditionSum(res_rep$SIRRes$filteringLogWeights[n,1,])))
    })
    maxWeightSIS[dim,] <- sapply(expRes[[dim]]$res, function(res_rep) {
      max(exp(res_rep$SISRes$filteringLogWeights[n,1,]-logAdditionSum(res_rep$SISRes$filteringLogWeights[n,1,])))
    })
    maxWeightBlock[dim,] <- sapply(expRes[[dim]]$res, function(res_rep) {
      max(apply(matrix(res_rep$BlockRes$filteringLogWeights[n,,],nrow=card_block), 1, function(x) {max(exp(x-logAdditionSum(x)))}))
    })
  }
}

dfWeights <- data.frame(Algorithm = character(),
                        Time = numeric(),
                        Repetition = numeric(),
                        Value = numeric(),
                        stringsAsFactors = FALSE)
for(rep in 1:repetitions) {
  for(dim in possible_dimension) {
    dfWeights <- rbind(dfWeights,c(Algorithm = "SIR",
                                   Dimension = dim,
                                   Repetition = rep,
                                   Value = maxWeightSIR[dim,rep]),
                       stringsAsFactors = FALSE)
    dfWeights <- rbind(dfWeights,c(Algorithm = "SIS",
                                   Dimension = dim,
                                   Repetition = rep,
                                   Value = maxWeightSIS[dim,rep]),
                       stringsAsFactors = FALSE)
  }
}
colnames(dfWeights) <- c("Algorithm","Dimension","Repetition", "Value")
dfWeights$Algorithm <- as.factor(dfWeights$Algorithm)

if(Sys.info()["nodename"] == "greyplover.stats.ox.ac.uk" || Sys.info()["nodename"] == "greypartridge.stats.ox.ac.uk" ||
   Sys.info()["nodename"] == "greyheron.stats.ox.ac.uk" || Sys.info()["nodename"] == "greywagtail.stats.ox.ac.uk") {
  saveRDS(dfWeights, file = paste0("curse_of_dimensionality2_weight_",id,".RDS"))
}

########
# PLOT #
########

dfBiasMSEToPlot <- dfResBiasVar[[type_statistic_plot]][[comp_statistic_plot]]$dfResBiasMSEData
dfResVarToPlot <-  dfResBiasVar[[type_statistic_plot]][[comp_statistic_plot]]$dfResVar
dfResVarToPlot$Algorithm <- factor(dfResVarToPlot$Algorithm, c("SIR","SIS","GibbsPF","Block","LocGibbsPF","Kalman"))

#Variance
#Just to check the different type of relative variances:
# ggplot(dfVarToPlot[dfVarToPlot$Type != "Var",], aes(x = as.numeric(Dimension), y=as.numeric(Value), group=interaction(Algorithm,Type),
#                                                        linetype = Algorithm, color = Type)) + geom_line()
# #RelVar
# ggplot(dfVarToPlot[dfVarToPlot$Type == "RelVar",], aes(x = Dimension, y=Value, color = Algorithm)) +
#   geom_line(size = 1.2) +  scale_colour_brewer(palette = "Set1")

ggplot(dfResVarToPlot[dfResVarToPlot$Type == "RelVar",], aes(x = Dimension, y=Value, color = Algorithm)) +
  geom_point(size = 3, shape = 4) +  scale_colour_brewer(palette = "Dark2") +
  theme_grey(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Relative Variance") +
  theme_bw() +
  geom_smooth(method="lm",se=FALSE, size = 1.5)

#Bias
# ggplot(dfBiasMSEToPlot[dfBiasMSEToPlot$Type == "Bias",], aes(x = as.numeric(Dimension), y=abs(Value), color = Algorithm)) + geom_line() +
#   scale_colour_brewer(palette = "Set1")

ggplot(dfBiasMSEToPlot[dfBiasMSEToPlot$Type == "RelAbsBias",], aes(x = Dimension, y=Value, color = Algorithm)) +
  scale_colour_brewer(palette = "Set1") +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x),
               fun.ymax = function(x) mean(x) + sd(x),
               geom = "pointrange") +
  stat_summary(fun.y = mean,
               geom = "line") +
  #geom_smooth(method="lm",se=TRUE) + geom_hline(yintercept = 0, linetype="dashed") +
  ylab("Rel. Abs. Bias")

#Maximum Weights
ggplot(dfWeights, aes(x=as.numeric(Dimension), y=as.numeric(Value), color=Algorithm)) +
  geom_hline(yintercept = 1, linetype="dashed") +
  theme_grey(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Weight") + xlab("Dimension") +
  theme_bw() +
  stat_summary(geom="pointrange", fun.data=mean_cl_boot, fun.args=list(conf.int=1))+
  scale_colour_brewer(palette = "Dark2")
