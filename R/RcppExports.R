# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

KalmanFilterCpp <- function(m_0, C_0, F_matrix, G, V, W, y) {
    .Call('SMC5_KalmanFilterCpp', PACKAGE = 'SMC5', m_0, C_0, F_matrix, G, V, W, y)
}

KalmanSmoothingCpp <- function(m_0, C_0, F_matrix, G, V, W, y) {
    .Call('SMC5_KalmanSmoothingCpp', PACKAGE = 'SMC5', m_0, C_0, F_matrix, G, V, W, y)
}

KalmanSamplerCpp <- function(n, kalmanRes, time) {
    .Call('SMC5_KalmanSamplerCpp', PACKAGE = 'SMC5', n, kalmanRes, time)
}

#' Multivariate Normal Density
#'
#' This function computes the density of a multivariate normal.
#'
#' @param x a nxp matrix whose columns are the components of the multivariate normal.
#' @param mean a p-dimensional mean vector
#' @param sigma a pxp covariance matrix
#' @param logd if \code{TRUE} log densities are returned
#' @details This function is taken from \url{http://gallery.rcpp.org/articles/dmvnorm_arma/}
#' @return A n-dimensional column vector with the values of the density.
#' @author Nino Hardt, Dicko Ahmadou
#' @export
dmvnrmArma <- function(x, mean, sigma, logd = FALSE, diag = FALSE) {
    .Call('SMC5_dmvnrmArma', PACKAGE = 'SMC5', x, mean, sigma, logd, diag)
}

dnrmArma <- function(x, mean, sigma, logd = FALSE) {
    .Call('SMC5_dnrmArma', PACKAGE = 'SMC5', x, mean, sigma, logd)
}

#' Simulate from a Multivariate Normal Density
#'
#' This function produces samples from a Multivariate Normal distribution.
#'
#' @param n size of sample required
#' @param mean a p-dimensional mean vector
#' @param sigma a pxp covariance matrix
#' @details This function is taken from \url{http://gallery.rcpp.org/articles/simulate-multivariate-normal/}
#' @return A nxp matrix with the generated sample
#' @author Dicko Ahmadou
#' @export
mvrnormArma <- function(n, mean, sigma, diag = FALSE) {
    .Call('SMC5_mvrnormArma', PACKAGE = 'SMC5', n, mean, sigma, diag)
}

rnormArma <- function(n, mean, sigma) {
    .Call('SMC5_rnormArma', PACKAGE = 'SMC5', n, mean, sigma)
}

blockForwardSmoothingSMC4 <- function(filteringResults, filteringLogWeights, blocks, fParams) {
    .Call('SMC5_blockForwardSmoothingSMC4', PACKAGE = 'SMC5', filteringResults, filteringLogWeights, blocks, fParams)
}

blockForwardSmoothingOnline <- function(n, filteringParticlesTemp, filteringLogWeightsTemp, previous_alpha, blocks, fParams) {
    .Call('SMC5_blockForwardSmoothingOnline', PACKAGE = 'SMC5', n, filteringParticlesTemp, filteringLogWeightsTemp, previous_alpha, blocks, fParams)
}

gibbsForwardSmoothingOnline <- function(n, filteringParticlesTemp, filteringLogWeightsTemp, previous_log_alpha, radius, fParams) {
    .Call('SMC5_gibbsForwardSmoothingOnline', PACKAGE = 'SMC5', n, filteringParticlesTemp, filteringLogWeightsTemp, previous_log_alpha, radius, fParams)
}

forwardSmoothingOnline <- function(n, filteringParticlesTemp, filteringLogWeightsTemp, previous_log_alpha, fParams) {
    .Call('SMC5_forwardSmoothingOnline', PACKAGE = 'SMC5', n, filteringParticlesTemp, filteringLogWeightsTemp, previous_log_alpha, fParams)
}

blockForwardSmoothing <- function(n, filteringParticlesTemp, filteringLogWeightsTemp, blocks, fParams) {
    .Call('SMC5_blockForwardSmoothing', PACKAGE = 'SMC5', n, filteringParticlesTemp, filteringLogWeightsTemp, blocks, fParams)
}

gibbsForwardSmoothing <- function(n, filteringParticlesTemp, filteringLogWeightsTemp, radius, fParams) {
    .Call('SMC5_gibbsForwardSmoothing', PACKAGE = 'SMC5', n, filteringParticlesTemp, filteringLogWeightsTemp, radius, fParams)
}

forwardSmoothing <- function(n, filteringParticlesTemp, filteringLogWeightsTemp, fParams) {
    .Call('SMC5_forwardSmoothing', PACKAGE = 'SMC5', n, filteringParticlesTemp, filteringLogWeightsTemp, fParams)
}

blockParticleFilterOnline <- function(N, n, particles, logWeights, blocks, fParams, gParams, resampling = TRUE, init = FALSE) {
    .Call('SMC5_blockParticleFilterOnline', PACKAGE = 'SMC5', N, n, particles, logWeights, blocks, fParams, gParams, resampling, init)
}

blockParticleFilter <- function(N, n, blocks, fParams, gParams, resampling = TRUE) {
    .Call('SMC5_blockParticleFilter', PACKAGE = 'SMC5', N, n, blocks, fParams, gParams, resampling)
}

sequentialImportanceSampling <- function(N, n, fParams, gParams) {
    .Call('SMC5_sequentialImportanceSampling', PACKAGE = 'SMC5', N, n, fParams, gParams)
}

sequentialImportanceSamplingOnline <- function(N, n, particles, logWeights, fParams, gParams) {
    .Call('SMC5_sequentialImportanceSamplingOnline', PACKAGE = 'SMC5', N, n, particles, logWeights, fParams, gParams)
}

bootstrapParticleFilter <- function(N, n, fParams, gParams) {
    .Call('SMC5_bootstrapParticleFilter', PACKAGE = 'SMC5', N, n, fParams, gParams)
}

bootstrapParticleFilterOnline <- function(N, n, particles, logWeights, fParams, gParams) {
    .Call('SMC5_bootstrapParticleFilterOnline', PACKAGE = 'SMC5', N, n, particles, logWeights, fParams, gParams)
}

gibbsParticleFilterOnline <- function(N, n, m, radius, particles, fParams, gParams, init = FALSE) {
    .Call('SMC5_gibbsParticleFilterOnline', PACKAGE = 'SMC5', N, n, m, radius, particles, fParams, gParams, init)
}

gibbsParticleFilter <- function(N, n, m, radius, fParams, gParams) {
    .Call('SMC5_gibbsParticleFilter', PACKAGE = 'SMC5', N, n, m, radius, fParams, gParams)
}

neighbourSMC4 <- function(component, radius, dimension) {
    .Call('SMC5_neighbourSMC4', PACKAGE = 'SMC5', component, radius, dimension)
}

gibbsParticleFilterSMC4 <- function(N, n, m, radius, fParams, gParams) {
    .Call('SMC5_gibbsParticleFilterSMC4', PACKAGE = 'SMC5', N, n, m, radius, fParams, gParams)
}

logAddition <- function(x, y) {
    .Call('SMC5_logAddition', PACKAGE = 'SMC5', x, y)
}

logAdditionSum <- function(x) {
    .Call('SMC5_logAdditionSum', PACKAGE = 'SMC5', x)
}

ProbSampleReplace <- function(nOrig, size, prob) {
    .Call('SMC5_ProbSampleReplace', PACKAGE = 'SMC5', nOrig, size, prob)
}

checkDiagonal <- function(X) {
    .Call('SMC5_checkDiagonal', PACKAGE = 'SMC5', X)
}

checkSymmetric <- function(X) {
    .Call('SMC5_checkSymmetric', PACKAGE = 'SMC5', X)
}

