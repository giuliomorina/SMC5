computeTrueStatisticFilter <- function(X_data, n, type = c("sum","sum_squared","mean","mean_squared"), comp = NA) {
  type <- match.arg(type)
  if(type == "sum") {
    return(sum(X_data[n,]))
  } else if(type == "sum_squared") {
    return(sum(X_data[n,]^2))
  } else if(type == "mean" && !is.null(comp) && !is.na(comp)) {
    return(X_data[n,comp])
  } else if(type == "mean_squared" && !is.null(comp) && !is.na(comp)) {
    return(X_data[n,comp]^2)
  }
  else {
    stop("comp must be specified for this statistic.")
  }
}

computeKalmanStatisticFilter <- function(fParams, gParams, n, type = c("sum","sum_squared","mean","mean_squared"), comp = NA) {
  dimension <- ncol(fParams$A)
  #Kalman Filter
  kalmanFilterRes <- KalmanFilterCpp(m_0 = fParams$X0,
                                     C_0 = fParams$sigmaX,
                                     F_matrix = diag(dimension),
                                     G = fParams$A,
                                     V = gParams$sigmaY,
                                     W = fParams$sigmaX,
                                     y = gParams$Y)
  type <- match.arg(type)
  if(type == "sum") {
    #True statistic is the sum of the mean values over all the components
    return(sum(kalmanFilterRes$m[,n+1]))
  } else if(type == "mean") {
    #True statistic is the mean in the specified component
    if(is.null(comp) || is.na(comp)) { stop("comp must be specified for this statistic.")}
    return(kalmanFilterRes$m[comp,n+1])
  } else if(type == "mean_squared") {
    #True statistic is the mean of the specified component squared
    #which is the variance + the mean squared
    return(kalmanFilterRes$C[comp,comp,n+1] + kalmanFilterRes$m[comp,n+1]^2)
  } else if(type == "sum_squared") {
    #True statistic is the sum of the mean of each component squared
    return(sum(diag(kalmanFilterRes$C[,,n+1]) + kalmanFilterRes$m[,n+1]^2))
  }
}

computeApproxStatisticFilter <- function(particles, logWeights, n, blocks = NULL, type = c("sum","sum_squared","mean","mean_squared"), comp = NA) {
  type <- match.arg(type)
  dimension <- dim(particles)[2]
  if(is.null(blocks) && dim(logWeights)[2] == 1) {
    blocks <- list(1:dimension) #Corresponds to SIS or SIR
  } else if (is.null(blocks)) {
    stop("blocks must be specified.")
  }
  statistic_components <- rep(NA, dimension) #Compute the statistic in each component
  for(v in 1:dimension) {
    #Find in which block it is
    for(b in 1:length(blocks)) {
      if(v %in% blocks[[b]]) {
        #Compute the statistic x^2 of that component
        if(type == "sum" || type == "mean") {
          statistic_components[v] <- weighted.mean(particles[n,v,],exp(logWeights[n,b,] - max(logWeights[n,b,])))
        } else {
          statistic_components[v] <- weighted.mean(particles[n,v,]^2,exp(logWeights[n,b,] - max(logWeights[n,b,])))
        }
        break
      }
    }
  }

  if(anyNA(statistic_components)) {
    stop("Impossible computing the statistic (perhaps infinite weights?")
  }
  if(type == "sum" || type == "sum_squared") {
    return(sum(statistic_components))
  } else if(is.null(comp) || is.na(comp)) {
    stop("comp must be specified for this statistic.")
  }
  else {
    return(statistic_components[comp])
  }
}

computeDfBiasVar <- function(approxStatisticPF, trueStatistics, algorithmName, dependentVarName, dfRes = NULL) {
  if(is.vector(trueStatistics)) {
    trueStatistics <- matrix(trueStatistics, nrow = nrow(approxStatisticPF), ncol = ncol(approxStatisticPF), byrow = FALSE)
  }
  if(nrow(approxStatisticPF) != nrow(trueStatistics)) {
    stop("Incompatible dimensions.")
  }
  if(is.null(dfRes)) {
    dfRes <- data.frame(Algorithm = character(),
                        Time = double(),
                        Type = character(),
                        Value = double(),
                        stringsAsFactors = FALSE)
  }
  for(i in 1:nrow(trueStatistics)) {
    estAll <- list(Bias = mean(approxStatisticPF[i,]-trueStatistics[i,]),
                   RelAbsBias = mean(abs((approxStatisticPF[i,]-trueStatistics[i,])/trueStatistics[i,])),
                   Var = mean((approxStatisticPF[i,]-mean(approxStatisticPF[i,]))^2),
                   CoefVar = sd(approxStatisticPF[i,])/mean(approxStatisticPF[i,]),
                   RelVar = var(approxStatisticPF[i,])/abs(mean(approxStatisticPF[i,])), #http://www.statisticshowto.com/relative-variance/
                   RelVar2 = (sd(approxStatisticPF[i,])/mean(approxStatisticPF[i,]))^2,
                   RelVar4 = var(approxStatisticPF[i,])/abs(mean(trueStatistics[i])), #A little modification I made from RelVar
                   MSE = mean((approxStatisticPF[i,]-trueStatistics[i,])^2),
                   RMSE = sqrt(mean((approxStatisticPF[i,]-trueStatistics[i,])^2)))
    for(type in c("Bias", "RelAbsBias","Var", "CoefVar", "RelVar", "RelVar2", "RelVar4", "MSE","RMSE")) {
      dfRes <- rbind(dfRes, c(Algorithm = algorithmName,
                              Time = i,
                              Type = type,
                              Value = estAll[[type]]),
                     stringsAsFactors = FALSE)
    }
  }
  colnames(dfRes) <- c("Algorithm", "Time", "Type", "Value")
  #Other relative variance (Definition 2) of http://www.statisticshowto.com/relative-variance/
  for(i in 1:nrow(trueStatistics)) {
    variances <- as.numeric(dfRes[dfRes$Algorithm == algorithmName & dfRes$Type=="Var","Value"])
    dfRes <- rbind(dfRes, c(Algorithm = algorithmName,
                            Time = i,
                            Type = "RelVar3",
                            Value = as.numeric(dfRes[dfRes$Algorithm == algorithmName &
                                                       dfRes$Type=="Var" &
                                                       dfRes$Time == i,"Value"])/sum(variances)),
                   stringsAsFactors = FALSE)
  }
  dfRes$Time <- as.numeric(dfRes$Time)
  dfRes$Value <- as.numeric(dfRes$Value)
  dfRes$Type <- as.factor(dfRes$Type)
  colnames(dfRes) <- c("Algorithm", dependentVarName, "Type", "Value")
  return(dfRes)
}
