// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double logAddition(double x, double y) {
  //Perform logarithmic addition, i.e. log(exp(x)+exp(y))
  if(x > y) {
    return(x+log(1+exp(y-x)));
  } else {
    return(y+log(1+exp(x-y)));
  }
}

// [[Rcpp::export]]
double logAdditionSum(arma::vec &x) {
  //Perform logarithmic addition of all the elements in the vector
  double res = x(0);
  for(int i=1; i<x.size(); i++) {
    res = logAddition(res,x(i));
  }
  return(res);
}

void standardizeLogWeights(arma::cube& logWeights) {
  //Standardize the weights for each time
  for(int i=0; i<logWeights.n_rows; i++) {
    for(int j=0; j<logWeights.n_cols; j++) {
      vec logWeightsCurrent = logWeights(span(i),span(j),span());
      logWeights(span(i),span(j),span()) = logWeights(span(i),span(j),span())- logAdditionSum(logWeightsCurrent);
    }
  }
}

void standardizeLogVector(arma::vec& vector){
  vector = vector - logAdditionSum(vector);
}

// Unequal probability sampling with replacement
// [[Rcpp::export]]
arma::uvec ProbSampleReplace(int nOrig, int size, arma::vec prob){
  //arma_rng::set_seed_random(); //Set a random seed
  arma::uvec index(size, fill::zeros);
  double rU;
  int ii, jj;
  int nOrig_1 = nOrig - 1;
  arma::uvec perm = arma::sort_index(prob, "descend"); //descending sort of index
  prob = prob/sum(prob);
  prob = arma::sort(prob, "descend");  // descending sort of prob
  // cumulative probabilities
  prob = arma::cumsum(prob);
  // compute the sample
  for (ii = 0; ii < size; ii++) {
    rU = unif_rand();
    for (jj = 0; jj < nOrig_1; jj++) {
      if (rU <= prob[jj])
        break;
    }
    index[ii] = perm[jj];
  }
  return(index);
}

// [[Rcpp::export]]
bool checkDiagonal(arma::mat& X) {
  for(int i = 0; i < X.n_rows; i++) {
    for(int j = 0; j < X.n_cols; j++) {
      if(i != j && X(i,j) != 0) return(false);
    }
  }
  return(true);
}

// [[Rcpp::export]]
bool checkSymmetric(arma::mat& X) {
  return(approx_equal(trimatu(X),trimatl(X).t(),"absdiff", 0.002));
}

int modulo(int a, int b) {
  const int result = a % b;
  return result >= 0 ? result : result + b;
}
