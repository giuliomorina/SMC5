# SMOOTHING CURSE OF DIMENSIONALITY
# This is like in the paper of Finke and Singh, with the addition of the
# Gibbs Smoothing. This is done only when using the perfect filter (Kalman)

rm(list=ls())
require(SMC5)
require(parallel)
require(ggplot2)

if(Sys.info()["nodename"] == "greyplover.stats.ox.ac.uk" || Sys.info()["nodename"] == "greypartridge.stats.ox.ac.uk" ||
   Sys.info()["nodename"] == "greyheron.stats.ox.ac.uk" || Sys.info()["nodename"] == "greywagtail.stats.ox.ac.uk") {
  setwd("/homes/morina/Local_Project/Smoothing_res")
}

######################
# SETTING EXPERIMENT #
######################
id <- sample(1e7,size = 1)
#set.seed(872328,"L'Ecuyer-CMRG")
n <- 20
possible_dimension <- seq(from=20, to=300, by=40)
possible_cardinalities <- c(1,4,10)
possible_radius <- c(1,4)
ncores <- 20
repetitions <- 40
N <- 500 #Number of particles
A_generator <- c(0.5,0.2)
varX <- 1 #Variance of sigmaX
varY <- 1 #Variance of sigmaY
X0_point <- 0 #Initial point

#########################
# RUN FORWARD SMOOTHING #
#########################

expRes <- lapply(possible_dimension, function(dimension) {
  #GENERATE MODEL
  A <- generateA(A_generator, dimension)
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
    #GENERATE KALMAN PARTICLES
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
    kalmanLogWeights <- array(log(1/N),c(n,1,N))
    #STANDARD FORWARD SMOOTHING
    resStandard <- list(forwardSmoothing(n=n,filteringParticlesTemp = kalmanParticles,
                                         filteringLogWeightsTemp = kalmanLogWeights,
                                         fParams = fParams)$statistic)
    #BLOCKED FORWARD SMOOTHING
    resBlock <- lapply(possible_cardinalities, function(card) {
      blocks <- split(1:dimension,ceiling((1:dimension)/card))
      return(blockForwardSmoothing(n=n, filteringParticlesTemp = kalmanParticles,
                                   filteringLogWeightsTemp = kalmanLogWeights, blocks=blocks,
                                   fParams = fParams)$statistic)
    })
    #GIBBS FORWARD SMOOTHING
    resGibbs <- lapply(possible_radius, function(radius) {
      return(gibbsForwardSmoothing(n=n, filteringParticlesTemp = kalmanParticles,
                                   filteringLogWeightsTemp = kalmanLogWeights, radius=radius,
                                   fParams = fParams)$statistic)
    })

    res_out <- list(resStandard = resStandard,
                    resBlock = resBlock,
                    resGibbs = resGibbs,
                    dataStatistic = smoothingStatistic(X_data),
                    trueStatistic = smoothingStatisticKalman(fParams=fParams, gParams = gParams))

    if(Sys.info()["nodename"] == "greyplover.stats.ox.ac.uk" || Sys.info()["nodename"] == "greypartridge.stats.ox.ac.uk" ||
       Sys.info()["nodename"] == "greyheron.stats.ox.ac.uk" || Sys.info()["nodename"] == "greywagtail.stats.ox.ac.uk") {
      saveRDS(res_out, file = paste0("smoothing_curse_dimensionality_dim_",dimension,"_rep_",rep,"_id_",id,".RDS"))
    }

    return(res_out)
  }, mc.cores = ncores)
  if(Sys.info()["nodename"] == "greyplover.stats.ox.ac.uk" || Sys.info()["nodename"] == "greypartridge.stats.ox.ac.uk" ||
     Sys.info()["nodename"] == "greyheron.stats.ox.ac.uk" || Sys.info()["nodename"] == "greywagtail.stats.ox.ac.uk") {
    saveRDS(res, file = paste0("smoothing_curse_dimensionality_FINAL_dim_",dimension,"_id_",id,".RDS"))
  }
  return(list(res=res))
})

###########################
# RECOVER PREVIOUS RESULT #
###########################

if(!exists("expRes")) {
  id <- 7903543
  directory_path <- "/Users/Giulio/Dropbox/OxWaSP/2nd Mini Project/Server_result/Smoothing_res/"
  expRes <- lapply(possible_dimension, function(dimension) {
    res <- lapply(1:repetitions, function(rep) {
      res_out <- readRDS(paste0(directory_path,"smoothing_curse_dimensionality_dim_",dimension,"_rep_",rep,"_id_",id,".RDS"))
      return(res_out)
    })
    return(list(res=res,
                trueStatistic = res[[1]]$trueStatistic))
  })
} else {
  print("expRes already exists.")
}


#########################
# APPROXIMATE STATISTIC #
#########################

dfResList <- vector("list")
counter_list <- 1
#Standard
approxStatistic <- matrix(NA, nrow=length(possible_dimension), ncol=repetitions)
for(counter in 1:length(possible_dimension)) {
  approxStatistic[counter,] <- sapply(expRes[[counter]]$res, function(res) {
    return(res$resStandard[[1]]/possible_dimension[counter])
  })
}
dfResList[[counter_list]] <- approxStatistic
counter_list <- counter_list + 1
#Block
for(i in 1:length(possible_cardinalities)) {
  for(counter in 1:length(possible_dimension)) {
    approxStatistic[counter,] <- sapply(expRes[[counter]]$res, function(res) {
      return(res$resBlock[[i]]/possible_dimension[counter])
    })
  }
  dfResList[[counter_list]] <- approxStatistic
  counter_list <- counter_list + 1
}
#Gibbs
for(i in 1:length(possible_radius)) {
  for(counter in 1:length(possible_dimension)) {
    approxStatistic[counter,] <- sapply(expRes[[counter]]$res, function(res) {
      return(res$resGibbs[[i]]/possible_dimension[counter])
    })
  }
  dfResList[[counter_list]] <- approxStatistic
  counter_list <- counter_list + 1
}
#True Statistic
for(counter in 1:length(possible_dimension)) {
  approxStatistic[counter,] <- sapply(expRes[[counter]]$res, function(res) {
    return(res$trueStatistic/possible_dimension[counter])
  })
}
dfResList[[counter_list]] <- approxStatistic
counter_list <- counter_list + 1
#Data statistic
for(counter in 1:length(possible_dimension)) {
  approxStatistic[counter,] <- sapply(expRes[[counter]]$res, function(res) {
    return(res$dataStatistic/possible_dimension[counter])
  })
}
dfResList[[counter_list]] <- approxStatistic
counter_list <- counter_list + 1
names(dfResList)[1] <- "Standard"
names(dfResList)[2:(length(possible_cardinalities)+1)] <- paste0("BlockCard",possible_cardinalities)
names(dfResList)[(length(possible_cardinalities)+2):(length(dfResList)-2)] <- paste0("GibbsRadius",possible_radius)
names(dfResList)[length(dfResList)-1] <- "TrueStatistic"
names(dfResList)[length(dfResList)] <- "DataStatistic"


if(Sys.info()["nodename"] == "greyplover.stats.ox.ac.uk" || Sys.info()["nodename"] == "greypartridge.stats.ox.ac.uk" ||
   Sys.info()["nodename"] == "greyheron.stats.ox.ac.uk" || Sys.info()["nodename"] == "greywagtail.stats.ox.ac.uk") {
  saveRDS(dfResList, file = paste0("smoothing_curse_dimensionality_",id,".RDS"))
}

###############
# COMPUTE MSE #
###############

dfResMSE <- computeDfBiasMSE(approxList = dfResList[1:(length(dfResList)-2)],
                             trueStatistics =  dfResList[["DataStatistic"]],
                             dependentVar = possible_dimension,
                             dependentVarName = "Dimension")

dfResRMSEToPlot <- data.frame(Repetitions = numeric(),
                              Dimension = numeric(),
                              Algorithm = character(),
                              MSE = numeric(),
                              RMSE = numeric(),
                              stringsAsFactors = FALSE)
for(dimension in possible_dimension) {
  subset <- dfResMSE[dfResMSE$Dimension == dimension,]
  for (algorithm in unique(subset$Algorithm)) {
    subset2 <- subset[subset$Algorithm == algorithm & subset$Type == "MSE",]
    dfResRMSEToPlot <- rbind(dfResRMSEToPlot,c(Repetitions = nrow(subset2),
                                               Dimension = dimension,
                                               Algorithm = algorithm,
                                               MSE = mean(subset2$Value),
                                               RMSE = sqrt(mean(subset2$Value))),
                             stringsAsFactors = FALSE)
  }
}
names(dfResRMSEToPlot) <- c("Repetitions","Dimension","Algorithm","MSE","RMSE")
dfResRMSEToPlot$Repetitions <- as.numeric(dfResRMSEToPlot$Repetitions)
dfResRMSEToPlot$Dimension <- as.numeric(dfResRMSEToPlot$Dimension)
dfResRMSEToPlot$Algorithm <- as.factor(dfResRMSEToPlot$Algorithm)
dfResRMSEToPlot$MSE <- as.numeric(dfResRMSEToPlot$MSE)
dfResRMSEToPlot$RMSE <- as.numeric(dfResRMSEToPlot$RMSE)

########
# PLOT #
########
ggplot(dfResRMSEToPlot, aes(x = Dimension, y=RMSE, color = Algorithm)) +
  scale_colour_brewer(palette = "Set1") +
  geom_point() + geom_line() +
  ylab("RMSE")


