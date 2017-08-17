computeTrueStatisticFilter <- function(X_data, n, type = c("sum","sum_squared","mean","mean_squared"), comp = NA) {
  type <- match.arg(type)
  if(type == "sum") {
    return(sum(X_data[n,]))
  } else if(type == "sum_squared") {
    return(sum(X_data[n,]^2))
  } else if(type == "mean" && !is.null(comp) && !is.na(comp)) {
    return(sum(X_data[n,comp]))
  } else if(type == "mean_squared" && !is.null(comp) && !is.na(comp)) {
    return(sum(X_data[n,comp]^2))
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
