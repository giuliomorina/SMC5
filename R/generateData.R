generateData <- function(n, A, X0, sigmaX, sigmaY) {
  if(!checkDiagonal(sigmaY)) stop("sigmaY must be diagonal.")
  if(!checkDiagonal(sigmaX)) stop("sigmaX must be diagonal.")
  dimension <- nrow(A)
  X_data_temp <- matrix(NA, ncol = dimension, nrow = n+1) #Cause it starts from n=0
  Y_data <- matrix(NA, ncol = dimension, nrow = n)
  X_data_temp[1,] <- X0
  for(i in 2:(n+1)) {
    X_data_temp[i,] <- rmvnorm(1, mean = A%*%X_data_temp[i-1,], sigma = sigmaX)
    Y_data[i-1,] <- rmvnorm(1, mean = X_data_temp[i,], sigma = sigmaY)
  }
  if(n == 1) {
    X_data <- matrix(X_data_temp[2,], nrow=1)
  } else {
    X_data <- X_data_temp[2:nrow(X_data_temp),]
  }
  return(list(X_data = X_data, Y_data = Y_data))
}

generateY <- function(X_data, sigmaY) {
  Y_data <- matrix(NA, ncol = ncol(X_data), nrow = nrow(X_data))
  for(i in 1:nrow(X_data)) {
    Y_data[i,] <- rmvnorm(1, mean = X_data[i,], sigma = sigmaY)
  }
  return(Y_data)
}

generateA <- function(neighbour_correlation, dimension) {
  if((length(neighbour_correlation)*2-2) > dimension) stop("Dimension must be bigger.")
  if(dimension == 1) return(matrix(neighbour_correlation, nrow=1, ncol=1))
  if(length(neighbour_correlation)*2-2 == dimension) {
    first_row <- neighbour_correlation
  } else {
    first_row <- c(neighbour_correlation,rep(0,dimension-length(neighbour_correlation)*2+1), rev(neighbour_correlation[-1]))
  }
  A <- matrix(first_row, nrow=length(first_row), ncol=length(first_row), byrow = TRUE)
  for(i in 2:length(first_row)) {
    A[i,] <- shift(first_row,i-1)
  }
  return(A)
}