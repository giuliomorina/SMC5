// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// KalmanFilterCpp
List KalmanFilterCpp(arma::vec& m_0, arma::mat& C_0, arma::mat& F_matrix, arma::mat& G, arma::mat& V, arma::mat& W, arma::mat& y);
RcppExport SEXP SMC5_KalmanFilterCpp(SEXP m_0SEXP, SEXP C_0SEXP, SEXP F_matrixSEXP, SEXP GSEXP, SEXP VSEXP, SEXP WSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type m_0(m_0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type C_0(C_0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type F_matrix(F_matrixSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(KalmanFilterCpp(m_0, C_0, F_matrix, G, V, W, y));
    return rcpp_result_gen;
END_RCPP
}
// KalmanSmoothingCpp
List KalmanSmoothingCpp(arma::colvec& m_0, arma::mat& C_0, arma::mat& F_matrix, arma::mat& G, arma::mat& V, arma::mat& W, arma::mat& y);
RcppExport SEXP SMC5_KalmanSmoothingCpp(SEXP m_0SEXP, SEXP C_0SEXP, SEXP F_matrixSEXP, SEXP GSEXP, SEXP VSEXP, SEXP WSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec& >::type m_0(m_0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type C_0(C_0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type F_matrix(F_matrixSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(KalmanSmoothingCpp(m_0, C_0, F_matrix, G, V, W, y));
    return rcpp_result_gen;
END_RCPP
}
// KalmanSamplerCpp
arma::mat KalmanSamplerCpp(int n, List kalmanRes, int time);
RcppExport SEXP SMC5_KalmanSamplerCpp(SEXP nSEXP, SEXP kalmanResSEXP, SEXP timeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< List >::type kalmanRes(kalmanResSEXP);
    Rcpp::traits::input_parameter< int >::type time(timeSEXP);
    rcpp_result_gen = Rcpp::wrap(KalmanSamplerCpp(n, kalmanRes, time));
    return rcpp_result_gen;
END_RCPP
}
// dmvnrmArma
arma::vec dmvnrmArma(arma::mat& x, arma::rowvec& mean, arma::mat& sigma, bool logd, bool diag);
RcppExport SEXP SMC5_dmvnrmArma(SEXP xSEXP, SEXP meanSEXP, SEXP sigmaSEXP, SEXP logdSEXP, SEXP diagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::rowvec& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type logd(logdSEXP);
    Rcpp::traits::input_parameter< bool >::type diag(diagSEXP);
    rcpp_result_gen = Rcpp::wrap(dmvnrmArma(x, mean, sigma, logd, diag));
    return rcpp_result_gen;
END_RCPP
}
// dnrmArma
double dnrmArma(double x, double mean, double sigma, bool logd);
RcppExport SEXP SMC5_dnrmArma(SEXP xSEXP, SEXP meanSEXP, SEXP sigmaSEXP, SEXP logdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type logd(logdSEXP);
    rcpp_result_gen = Rcpp::wrap(dnrmArma(x, mean, sigma, logd));
    return rcpp_result_gen;
END_RCPP
}
// mvrnormArma
arma::mat mvrnormArma(int n, arma::rowvec& mean, arma::mat& sigma, bool diag);
RcppExport SEXP SMC5_mvrnormArma(SEXP nSEXP, SEXP meanSEXP, SEXP sigmaSEXP, SEXP diagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::rowvec& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type diag(diagSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnormArma(n, mean, sigma, diag));
    return rcpp_result_gen;
END_RCPP
}
// rnormArma
arma::colvec rnormArma(int n, double mean, double sigma);
RcppExport SEXP SMC5_rnormArma(SEXP nSEXP, SEXP meanSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(rnormArma(n, mean, sigma));
    return rcpp_result_gen;
END_RCPP
}
// blockForwardSmoothingSMC4
arma::cube blockForwardSmoothingSMC4(arma::cube filteringResults, arma::cube filteringLogWeights, List blocks, List fParams);
RcppExport SEXP SMC5_blockForwardSmoothingSMC4(SEXP filteringResultsSEXP, SEXP filteringLogWeightsSEXP, SEXP blocksSEXP, SEXP fParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type filteringResults(filteringResultsSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type filteringLogWeights(filteringLogWeightsSEXP);
    Rcpp::traits::input_parameter< List >::type blocks(blocksSEXP);
    Rcpp::traits::input_parameter< List >::type fParams(fParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(blockForwardSmoothingSMC4(filteringResults, filteringLogWeights, blocks, fParams));
    return rcpp_result_gen;
END_RCPP
}
// blockForwardSmoothingOnline
List blockForwardSmoothingOnline(int n, arma::cube filteringParticlesTemp, arma::cube filteringLogWeightsTemp, arma::cube previous_alpha, List blocks, List fParams);
RcppExport SEXP SMC5_blockForwardSmoothingOnline(SEXP nSEXP, SEXP filteringParticlesTempSEXP, SEXP filteringLogWeightsTempSEXP, SEXP previous_alphaSEXP, SEXP blocksSEXP, SEXP fParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type filteringParticlesTemp(filteringParticlesTempSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type filteringLogWeightsTemp(filteringLogWeightsTempSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type previous_alpha(previous_alphaSEXP);
    Rcpp::traits::input_parameter< List >::type blocks(blocksSEXP);
    Rcpp::traits::input_parameter< List >::type fParams(fParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(blockForwardSmoothingOnline(n, filteringParticlesTemp, filteringLogWeightsTemp, previous_alpha, blocks, fParams));
    return rcpp_result_gen;
END_RCPP
}
// gibbsForwardSmoothingOnline
List gibbsForwardSmoothingOnline(int n, arma::cube filteringParticlesTemp, arma::cube filteringLogWeightsTemp, arma::cube previous_log_alpha, int radius, List fParams);
RcppExport SEXP SMC5_gibbsForwardSmoothingOnline(SEXP nSEXP, SEXP filteringParticlesTempSEXP, SEXP filteringLogWeightsTempSEXP, SEXP previous_log_alphaSEXP, SEXP radiusSEXP, SEXP fParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type filteringParticlesTemp(filteringParticlesTempSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type filteringLogWeightsTemp(filteringLogWeightsTempSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type previous_log_alpha(previous_log_alphaSEXP);
    Rcpp::traits::input_parameter< int >::type radius(radiusSEXP);
    Rcpp::traits::input_parameter< List >::type fParams(fParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(gibbsForwardSmoothingOnline(n, filteringParticlesTemp, filteringLogWeightsTemp, previous_log_alpha, radius, fParams));
    return rcpp_result_gen;
END_RCPP
}
// forwardSmoothingOnline
List forwardSmoothingOnline(int n, arma::cube filteringParticlesTemp, arma::cube filteringLogWeightsTemp, arma::cube previous_log_alpha, List fParams);
RcppExport SEXP SMC5_forwardSmoothingOnline(SEXP nSEXP, SEXP filteringParticlesTempSEXP, SEXP filteringLogWeightsTempSEXP, SEXP previous_log_alphaSEXP, SEXP fParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type filteringParticlesTemp(filteringParticlesTempSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type filteringLogWeightsTemp(filteringLogWeightsTempSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type previous_log_alpha(previous_log_alphaSEXP);
    Rcpp::traits::input_parameter< List >::type fParams(fParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(forwardSmoothingOnline(n, filteringParticlesTemp, filteringLogWeightsTemp, previous_log_alpha, fParams));
    return rcpp_result_gen;
END_RCPP
}
// blockForwardSmoothing
List blockForwardSmoothing(int n, arma::cube filteringParticlesTemp, arma::cube filteringLogWeightsTemp, List blocks, List fParams);
RcppExport SEXP SMC5_blockForwardSmoothing(SEXP nSEXP, SEXP filteringParticlesTempSEXP, SEXP filteringLogWeightsTempSEXP, SEXP blocksSEXP, SEXP fParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type filteringParticlesTemp(filteringParticlesTempSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type filteringLogWeightsTemp(filteringLogWeightsTempSEXP);
    Rcpp::traits::input_parameter< List >::type blocks(blocksSEXP);
    Rcpp::traits::input_parameter< List >::type fParams(fParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(blockForwardSmoothing(n, filteringParticlesTemp, filteringLogWeightsTemp, blocks, fParams));
    return rcpp_result_gen;
END_RCPP
}
// gibbsForwardSmoothing
List gibbsForwardSmoothing(int n, arma::cube filteringParticlesTemp, arma::cube filteringLogWeightsTemp, int radius, List fParams);
RcppExport SEXP SMC5_gibbsForwardSmoothing(SEXP nSEXP, SEXP filteringParticlesTempSEXP, SEXP filteringLogWeightsTempSEXP, SEXP radiusSEXP, SEXP fParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type filteringParticlesTemp(filteringParticlesTempSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type filteringLogWeightsTemp(filteringLogWeightsTempSEXP);
    Rcpp::traits::input_parameter< int >::type radius(radiusSEXP);
    Rcpp::traits::input_parameter< List >::type fParams(fParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(gibbsForwardSmoothing(n, filteringParticlesTemp, filteringLogWeightsTemp, radius, fParams));
    return rcpp_result_gen;
END_RCPP
}
// forwardSmoothing
List forwardSmoothing(int n, arma::cube filteringParticlesTemp, arma::cube filteringLogWeightsTemp, List fParams);
RcppExport SEXP SMC5_forwardSmoothing(SEXP nSEXP, SEXP filteringParticlesTempSEXP, SEXP filteringLogWeightsTempSEXP, SEXP fParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type filteringParticlesTemp(filteringParticlesTempSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type filteringLogWeightsTemp(filteringLogWeightsTempSEXP);
    Rcpp::traits::input_parameter< List >::type fParams(fParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(forwardSmoothing(n, filteringParticlesTemp, filteringLogWeightsTemp, fParams));
    return rcpp_result_gen;
END_RCPP
}
// blockParticleFilterOnline
List blockParticleFilterOnline(int N, int n, arma::cube& particles, arma::cube& logWeights, List blocks, List fParams, List gParams, bool resampling, bool init);
RcppExport SEXP SMC5_blockParticleFilterOnline(SEXP NSEXP, SEXP nSEXP, SEXP particlesSEXP, SEXP logWeightsSEXP, SEXP blocksSEXP, SEXP fParamsSEXP, SEXP gParamsSEXP, SEXP resamplingSEXP, SEXP initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type particles(particlesSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type logWeights(logWeightsSEXP);
    Rcpp::traits::input_parameter< List >::type blocks(blocksSEXP);
    Rcpp::traits::input_parameter< List >::type fParams(fParamsSEXP);
    Rcpp::traits::input_parameter< List >::type gParams(gParamsSEXP);
    Rcpp::traits::input_parameter< bool >::type resampling(resamplingSEXP);
    Rcpp::traits::input_parameter< bool >::type init(initSEXP);
    rcpp_result_gen = Rcpp::wrap(blockParticleFilterOnline(N, n, particles, logWeights, blocks, fParams, gParams, resampling, init));
    return rcpp_result_gen;
END_RCPP
}
// blockParticleFilter
List blockParticleFilter(int N, int n, List blocks, List fParams, List gParams, bool resampling);
RcppExport SEXP SMC5_blockParticleFilter(SEXP NSEXP, SEXP nSEXP, SEXP blocksSEXP, SEXP fParamsSEXP, SEXP gParamsSEXP, SEXP resamplingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< List >::type blocks(blocksSEXP);
    Rcpp::traits::input_parameter< List >::type fParams(fParamsSEXP);
    Rcpp::traits::input_parameter< List >::type gParams(gParamsSEXP);
    Rcpp::traits::input_parameter< bool >::type resampling(resamplingSEXP);
    rcpp_result_gen = Rcpp::wrap(blockParticleFilter(N, n, blocks, fParams, gParams, resampling));
    return rcpp_result_gen;
END_RCPP
}
// sequentialImportanceSampling
List sequentialImportanceSampling(int N, int n, List fParams, List gParams);
RcppExport SEXP SMC5_sequentialImportanceSampling(SEXP NSEXP, SEXP nSEXP, SEXP fParamsSEXP, SEXP gParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< List >::type fParams(fParamsSEXP);
    Rcpp::traits::input_parameter< List >::type gParams(gParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(sequentialImportanceSampling(N, n, fParams, gParams));
    return rcpp_result_gen;
END_RCPP
}
// sequentialImportanceSamplingOnline
List sequentialImportanceSamplingOnline(int N, int n, arma::cube& particles, arma::cube& logWeights, List fParams, List gParams);
RcppExport SEXP SMC5_sequentialImportanceSamplingOnline(SEXP NSEXP, SEXP nSEXP, SEXP particlesSEXP, SEXP logWeightsSEXP, SEXP fParamsSEXP, SEXP gParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type particles(particlesSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type logWeights(logWeightsSEXP);
    Rcpp::traits::input_parameter< List >::type fParams(fParamsSEXP);
    Rcpp::traits::input_parameter< List >::type gParams(gParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(sequentialImportanceSamplingOnline(N, n, particles, logWeights, fParams, gParams));
    return rcpp_result_gen;
END_RCPP
}
// bootstrapParticleFilter
List bootstrapParticleFilter(int N, int n, List fParams, List gParams);
RcppExport SEXP SMC5_bootstrapParticleFilter(SEXP NSEXP, SEXP nSEXP, SEXP fParamsSEXP, SEXP gParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< List >::type fParams(fParamsSEXP);
    Rcpp::traits::input_parameter< List >::type gParams(gParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(bootstrapParticleFilter(N, n, fParams, gParams));
    return rcpp_result_gen;
END_RCPP
}
// bootstrapParticleFilterOnline
List bootstrapParticleFilterOnline(int N, int n, arma::cube& particles, arma::cube& logWeights, List fParams, List gParams);
RcppExport SEXP SMC5_bootstrapParticleFilterOnline(SEXP NSEXP, SEXP nSEXP, SEXP particlesSEXP, SEXP logWeightsSEXP, SEXP fParamsSEXP, SEXP gParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type particles(particlesSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type logWeights(logWeightsSEXP);
    Rcpp::traits::input_parameter< List >::type fParams(fParamsSEXP);
    Rcpp::traits::input_parameter< List >::type gParams(gParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(bootstrapParticleFilterOnline(N, n, particles, logWeights, fParams, gParams));
    return rcpp_result_gen;
END_RCPP
}
// gibbsParticleFilterOnline
List gibbsParticleFilterOnline(int N, int n, int m, int radius, arma::cube& particles, List fParams, List gParams, bool init);
RcppExport SEXP SMC5_gibbsParticleFilterOnline(SEXP NSEXP, SEXP nSEXP, SEXP mSEXP, SEXP radiusSEXP, SEXP particlesSEXP, SEXP fParamsSEXP, SEXP gParamsSEXP, SEXP initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type radius(radiusSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type particles(particlesSEXP);
    Rcpp::traits::input_parameter< List >::type fParams(fParamsSEXP);
    Rcpp::traits::input_parameter< List >::type gParams(gParamsSEXP);
    Rcpp::traits::input_parameter< bool >::type init(initSEXP);
    rcpp_result_gen = Rcpp::wrap(gibbsParticleFilterOnline(N, n, m, radius, particles, fParams, gParams, init));
    return rcpp_result_gen;
END_RCPP
}
// gibbsParticleFilter
List gibbsParticleFilter(int N, int n, int m, int radius, List fParams, List gParams);
RcppExport SEXP SMC5_gibbsParticleFilter(SEXP NSEXP, SEXP nSEXP, SEXP mSEXP, SEXP radiusSEXP, SEXP fParamsSEXP, SEXP gParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type radius(radiusSEXP);
    Rcpp::traits::input_parameter< List >::type fParams(fParamsSEXP);
    Rcpp::traits::input_parameter< List >::type gParams(gParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(gibbsParticleFilter(N, n, m, radius, fParams, gParams));
    return rcpp_result_gen;
END_RCPP
}
// neighbourSMC4
arma::uvec neighbourSMC4(int component, int radius, int dimension);
RcppExport SEXP SMC5_neighbourSMC4(SEXP componentSEXP, SEXP radiusSEXP, SEXP dimensionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type component(componentSEXP);
    Rcpp::traits::input_parameter< int >::type radius(radiusSEXP);
    Rcpp::traits::input_parameter< int >::type dimension(dimensionSEXP);
    rcpp_result_gen = Rcpp::wrap(neighbourSMC4(component, radius, dimension));
    return rcpp_result_gen;
END_RCPP
}
// gibbsParticleFilterSMC4
List gibbsParticleFilterSMC4(int N, int n, int m, int radius, List fParams, List gParams);
RcppExport SEXP SMC5_gibbsParticleFilterSMC4(SEXP NSEXP, SEXP nSEXP, SEXP mSEXP, SEXP radiusSEXP, SEXP fParamsSEXP, SEXP gParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type radius(radiusSEXP);
    Rcpp::traits::input_parameter< List >::type fParams(fParamsSEXP);
    Rcpp::traits::input_parameter< List >::type gParams(gParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(gibbsParticleFilterSMC4(N, n, m, radius, fParams, gParams));
    return rcpp_result_gen;
END_RCPP
}
// logAddition
double logAddition(double x, double y);
RcppExport SEXP SMC5_logAddition(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(logAddition(x, y));
    return rcpp_result_gen;
END_RCPP
}
// logAdditionSum
double logAdditionSum(arma::vec& x);
RcppExport SEXP SMC5_logAdditionSum(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logAdditionSum(x));
    return rcpp_result_gen;
END_RCPP
}
// ProbSampleReplace
arma::uvec ProbSampleReplace(int nOrig, int size, arma::vec prob);
RcppExport SEXP SMC5_ProbSampleReplace(SEXP nOrigSEXP, SEXP sizeSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nOrig(nOrigSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(ProbSampleReplace(nOrig, size, prob));
    return rcpp_result_gen;
END_RCPP
}
// checkDiagonal
bool checkDiagonal(arma::mat& X);
RcppExport SEXP SMC5_checkDiagonal(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(checkDiagonal(X));
    return rcpp_result_gen;
END_RCPP
}
// checkSymmetric
bool checkSymmetric(arma::mat& X);
RcppExport SEXP SMC5_checkSymmetric(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(checkSymmetric(X));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"SMC5_KalmanFilterCpp", (DL_FUNC) &SMC5_KalmanFilterCpp, 7},
    {"SMC5_KalmanSmoothingCpp", (DL_FUNC) &SMC5_KalmanSmoothingCpp, 7},
    {"SMC5_KalmanSamplerCpp", (DL_FUNC) &SMC5_KalmanSamplerCpp, 3},
    {"SMC5_dmvnrmArma", (DL_FUNC) &SMC5_dmvnrmArma, 5},
    {"SMC5_dnrmArma", (DL_FUNC) &SMC5_dnrmArma, 4},
    {"SMC5_mvrnormArma", (DL_FUNC) &SMC5_mvrnormArma, 4},
    {"SMC5_rnormArma", (DL_FUNC) &SMC5_rnormArma, 3},
    {"SMC5_blockForwardSmoothingSMC4", (DL_FUNC) &SMC5_blockForwardSmoothingSMC4, 4},
    {"SMC5_blockForwardSmoothingOnline", (DL_FUNC) &SMC5_blockForwardSmoothingOnline, 6},
    {"SMC5_gibbsForwardSmoothingOnline", (DL_FUNC) &SMC5_gibbsForwardSmoothingOnline, 6},
    {"SMC5_forwardSmoothingOnline", (DL_FUNC) &SMC5_forwardSmoothingOnline, 5},
    {"SMC5_blockForwardSmoothing", (DL_FUNC) &SMC5_blockForwardSmoothing, 5},
    {"SMC5_gibbsForwardSmoothing", (DL_FUNC) &SMC5_gibbsForwardSmoothing, 5},
    {"SMC5_forwardSmoothing", (DL_FUNC) &SMC5_forwardSmoothing, 4},
    {"SMC5_blockParticleFilterOnline", (DL_FUNC) &SMC5_blockParticleFilterOnline, 9},
    {"SMC5_blockParticleFilter", (DL_FUNC) &SMC5_blockParticleFilter, 6},
    {"SMC5_sequentialImportanceSampling", (DL_FUNC) &SMC5_sequentialImportanceSampling, 4},
    {"SMC5_sequentialImportanceSamplingOnline", (DL_FUNC) &SMC5_sequentialImportanceSamplingOnline, 6},
    {"SMC5_bootstrapParticleFilter", (DL_FUNC) &SMC5_bootstrapParticleFilter, 4},
    {"SMC5_bootstrapParticleFilterOnline", (DL_FUNC) &SMC5_bootstrapParticleFilterOnline, 6},
    {"SMC5_gibbsParticleFilterOnline", (DL_FUNC) &SMC5_gibbsParticleFilterOnline, 8},
    {"SMC5_gibbsParticleFilter", (DL_FUNC) &SMC5_gibbsParticleFilter, 6},
    {"SMC5_neighbourSMC4", (DL_FUNC) &SMC5_neighbourSMC4, 3},
    {"SMC5_gibbsParticleFilterSMC4", (DL_FUNC) &SMC5_gibbsParticleFilterSMC4, 6},
    {"SMC5_logAddition", (DL_FUNC) &SMC5_logAddition, 2},
    {"SMC5_logAdditionSum", (DL_FUNC) &SMC5_logAdditionSum, 1},
    {"SMC5_ProbSampleReplace", (DL_FUNC) &SMC5_ProbSampleReplace, 3},
    {"SMC5_checkDiagonal", (DL_FUNC) &SMC5_checkDiagonal, 1},
    {"SMC5_checkSymmetric", (DL_FUNC) &SMC5_checkSymmetric, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_SMC5(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
