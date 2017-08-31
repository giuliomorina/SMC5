#Compute the smoothing statistic of Finke Singh with r=0
smoothingStatistic <- function(X_data) {
  return(sum(X_data[2:nrow(X_data),]*X_data[1:(nrow(X_data)-1),]))
}
