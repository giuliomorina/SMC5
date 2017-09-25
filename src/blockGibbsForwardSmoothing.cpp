// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
#include "smc.h"

double fFunctional(arma::cube& filteringParticles, int time, int component,
                   int particleIndCurrent, int particleIndNext) {
  //This is f_t^v(x_{t-1},x_t) = x_t^v*x_{t-1}^v
  return(filteringParticles(time-1,component,particleIndCurrent) *
         filteringParticles(time,component,particleIndNext));
}

// [[Rcpp::export]]
List blockForwardSmoothingOnline(int n, arma::cube filteringParticlesTemp,
                                 arma::cube filteringLogWeightsTemp,
                                 arma::cube previous_alpha,
                                 List blocks, List fParams) {
  //Unroll fParams
  arma::mat A = fParams["A"];
  arma::rowvec X0 = fParams["X0"];
  arma::mat sigmaX = fParams["sigmaX"];
  //Check
  if(!checkDiagonal(sigmaX)) stop("sigmaX must be diagonal");
  if(!checkSymmetric(A)) stop("A must be symmetric.");
  if(filteringLogWeightsTemp.n_cols > 1) stop("the logWeights have wrong dimension.");
  //Unroll blocks into an array
  int number_blocks = blocks.size();
  arma::uvec *coord_array;
  coord_array = new arma::uvec[number_blocks]; // coord_array similar to a list in Rccp but the element need to be of the same type (not dimension)
  for(int i=0; i < number_blocks; i++) {
    arma::uvec coord = blocks[i];
    coord_array[i] = coord - 1; //-1 cause components passed by R start from 1
  }
  //Set up variables
  int dimension = A.n_rows;
  int N = filteringParticlesTemp.n_slices;
  arma::cube alpha(n+1, dimension, N, fill::zeros); //n+1 cause it starts from 0
  //Fill up with the already provided results
  alpha.subcube(0,0,0,previous_alpha.n_rows-1,dimension-1,N-1) = previous_alpha;
  //Needs to add back the results at time 0
  arma::cube filteringParticles(n+1, dimension, N, fill::ones);
  for(int p = 0; p < N; p++) {
    filteringParticles.slice(p).row(0) = X0;
  }
  filteringParticles.subcube(1,0,0,n,dimension-1,N-1) = filteringParticlesTemp;
  arma::cube filteringLogWeights(n+1, 1, N, fill::ones);
  filteringLogWeights *= log(1.0/N);
  filteringLogWeights.subcube(1,0,0,n,0,N-1) = filteringLogWeightsTemp;
  standardizeLogWeights(filteringLogWeights); //Standardize them
  //Auxiliary variables
  arma::colvec temp_log_weights(N);
  arma::colvec temp_denominator(N);
  arma::colvec temp_alpha(N);

  for(int b=0; b < number_blocks; b++) {
    arma::uvec coord = coord_array[b]; //Coordinates of this block
    for(int t=previous_alpha.n_rows; t < n+1; t++) { //+1 cause the first one is all zeros
      for(int i=0; i<N; i++) {
        temp_log_weights = filteringLogWeights(span(t-1),span(0),span());
        //Compute weights
        for(int j=0; j<N; j++) {
          for(int v=0; v < coord.size(); v++) {
            temp_log_weights(j) += fDensityComponent(filteringParticles,j,i,t,coord(v),A,sigmaX);
          }
        }
        //Standardize
        standardizeLogVector(temp_log_weights);

        //Compute alpha in each component of the block
        for(int v=0; v < coord.size(); v++) {
          temp_alpha = exp(temp_log_weights);
          for(int j=0; j<N; j++) {
            temp_alpha(j) *= (alpha(t-1,coord(v),j)+fFunctional(filteringParticles, t, coord(v), j, i));
          }
          alpha(t,coord(v),i) = sum(temp_alpha);
        }

      }
    }
  }

  //Compute final result
  double res = 0.0;
  for(int i=0; i<N; i++) {
    for(int v=0; v<dimension; v++) {
      res += exp(filteringLogWeights(n,0,i))*alpha(n,v,i);
    }
  }

  arma::cube alpha_res = alpha.subcube(1,0,0,n,dimension-1,N-1); //Remove the first so that it is consistent with filter results
  return(List::create(Named("statistic") = res,
                      Named("alpha") = alpha_res));
}

// [[Rcpp::export]]
List gibbsForwardSmoothingOnline(int n, arma::cube filteringParticlesTemp,
                                 arma::cube filteringLogWeightsTemp,
                                 arma::cube previous_log_alpha,
                                 int radius,
                                 List fParams) {
  //The blocks are the neighbour!
  int dimension = filteringParticlesTemp.n_cols;
  Rcpp::List blocks(dimension);
  for(int v=0; v < dimension; v++) {
    blocks[v] = neighbour(v,radius,dimension)+1; //+1 to emulate R where the coordinates start from 1
  }

  return(blockForwardSmoothingOnline(n, filteringParticlesTemp,
                                     filteringLogWeightsTemp,
                                     previous_log_alpha,
                                     blocks, fParams));
}

// [[Rcpp::export]]
List forwardSmoothingOnline(int n, arma::cube filteringParticlesTemp,
                            arma::cube filteringLogWeightsTemp,
                            arma::cube previous_log_alpha,
                            List fParams) {
  int dimension = filteringParticlesTemp.n_cols;
  List blocks = List::create(linspace(1,dimension,dimension).t());
  return(blockForwardSmoothingOnline(n, filteringParticlesTemp,
                                     filteringLogWeightsTemp,
                                     previous_log_alpha,
                                     blocks,
                                     fParams));
}
// [[Rcpp::export]]
List blockForwardSmoothing(int n, arma::cube filteringParticlesTemp,
                           arma::cube filteringLogWeightsTemp,
                           List blocks, List fParams) {
  arma::cube previous_log_alpha(1, filteringParticlesTemp.n_cols, filteringParticlesTemp.n_slices, fill::zeros);
  return(blockForwardSmoothingOnline(n, filteringParticlesTemp,
                                     filteringLogWeightsTemp,
                                     previous_log_alpha,
                                     blocks, fParams));
}

// [[Rcpp::export]]
List gibbsForwardSmoothing(int n, arma::cube filteringParticlesTemp,
                           arma::cube filteringLogWeightsTemp,
                           int radius,
                           List fParams) {
  arma::cube previous_log_alpha(1, filteringParticlesTemp.n_cols, filteringParticlesTemp.n_slices, fill::zeros);
  return(gibbsForwardSmoothingOnline(n, filteringParticlesTemp,
                                     filteringLogWeightsTemp,
                                     previous_log_alpha,
                                     radius,
                                     fParams));
}

// [[Rcpp::export]]
List forwardSmoothing(int n, arma::cube filteringParticlesTemp,
                      arma::cube filteringLogWeightsTemp,
                      List fParams) {
  int dimension = filteringParticlesTemp.n_cols;
  List blocks = List::create(linspace(1,dimension,dimension).t());
  return(blockForwardSmoothing(n, filteringParticlesTemp,
                               filteringLogWeightsTemp,
                               blocks,
                               fParams));
}
