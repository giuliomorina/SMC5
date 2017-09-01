#Compute the smoothing statistic of Finke Singh with r=0
smoothingStatistic <- function(X_data) {
  return(sum(X_data[2:nrow(X_data),]*X_data[1:(nrow(X_data)-1),]))
}

smoothingStatisticKalman <- function(fParams, gParams){
  A <- fParams$A
  X0 <- fParams$X0
  sigmaX <- fParams$sigmaX
  Y_data <- gParams$Y
  sigmaY <- gParams$sigmaY
  dimension <- length(X0)
  #Kalman filter
  kalmanFilterRes <- KalmanFilterCpp(m_0 = X0,
                                     C_0 = sigmaX,
                                     F_matrix = diag(dimension),
                                     G = A,
                                     V = sigmaY,
                                     W = sigmaX,
                                     y = Y_data)
  #Kalman Smoothing
  kalmanSmoothingRes <- KalmanSmoothingCpp(m_0 = X0,
                                           C_0 = sigmaX,
                                           F_matrix = diag(dimension),
                                           G = A,
                                           V = sigmaY,
                                           W = sigmaX,
                                           y = Y_data)
  sigmaX_inv <- solve(sigmaX)
  aux <- t(A) %*% sigmaX_inv %*% A
  time <- ncol(kalmanSmoothingRes$m)-1
  result <- sapply(2:time, function(t){
    sigma <- solve(aux + solve(kalmanFilterRes$C[,,t]))
    part1 <- t(kalmanSmoothingRes$m[,t+1]) %*% sigma %*% solve(kalmanFilterRes$C[,,t]) %*% kalmanFilterRes$m[,t]
    B <- sigma %*% t(A) %*% sigmaX_inv
    corrected_cov <- kalmanSmoothingRes$C[,,t+1] + kalmanSmoothingRes$m[,t+1]%*%t(kalmanSmoothingRes$m[,t+1]) #E[XiXj] = cov(Xi,Xj) + EXiExj
    part2 <- sum(corrected_cov*B)
    return(part1+part2)
  })
  return(sum(result))
}
