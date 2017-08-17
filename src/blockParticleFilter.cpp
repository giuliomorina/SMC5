// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
#include "smc.h"

arma::rowvec fSample(arma::rowvec currentX, arma::mat& A, arma::mat& sigmaX) {
  arma::rowvec meanXNew = currentX*A;
  return(mvrnormArma(1,meanXNew,sigmaX, true).row(0));
}

double gDensityLog(double currentX, double Y, double var) {
  //This is the g^v (only one component)
  return(dnrmArma(Y,currentX,sqrt(var),true));
}

// [[Rcpp::export]]
List blockParticleFilter(int N, int n, List blocks, List fParams, List gParams,
                         bool resampling = true) {
  //Unroll fParams and gParams
  arma::mat A = fParams["A"];
  arma::rowvec X0 = fParams["X0"];
  arma::mat sigmaX = fParams["sigmaX"];
  arma::mat Y = gParams["Y"];
  arma::mat sigmaY = gParams["sigmaY"];
  //Check
  if(!checkDiagonal(sigmaX)) stop("sigmaX must be diagonal");
  if(!checkDiagonal(sigmaY)) stop("sigmaY must be diagonal");
  if(!checkSymmetric(A)) stop("A must be symmetric.");
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
  arma::cube results_array(n+1, dimension, N, fill::zeros); //n+1 because we start from time zero
  arma::cube log_weights(n+1, number_blocks, N, fill::ones);
  log_weights *= log(1.0/N); //Fill all with equal weight
  //Fill the first time with X0 for all particles and same weights for everything
  arma::rowvec ones(dimension,fill::ones);
  for(int particle = 0; particle < N; particle++) {
    results_array.slice(particle).row(0) = X0;
  }
  //Auxiliary variable
  arma::mat resampled_particles(dimension,N);
  arma::colvec current_log_weights(N);
  arma::uvec which_sampled(N);

  //Main loop
  for(int t=1; t < n+1; t++) { //Iterating over time
    //RESAMPLING STEP
    //For each block resample the particle
    for(int b=0; b<number_blocks;b++) {
      //Get weight of this block
      arma::uvec coord = coord_array[b];
      if(resampling) {
        current_log_weights = log_weights.subcube(span(t-1), span(b), span());
        which_sampled = ProbSampleReplace(N, N, exp(current_log_weights-max(current_log_weights))); //Sample which particles given the probabilities
      } else {
        //No resampling scheme
        which_sampled = linspace<uvec>(0, N-1, N);
      }
      for(int v=0; v < coord.size(); v++) {
        for(int particle = 0; particle < N; particle++) {
          resampled_particles(coord(v),particle) = results_array(t-1,coord(v),which_sampled(particle));
        }
      }
    }

    for(int i=0; i<N; i++) {
      //PROPAGATE PARTICLE
      results_array.slice(i).row(t) = fSample(resampled_particles.col(i).t(), A, sigmaX);

      //CORRECT ITS WEIGHT
      for(int b=0; b<number_blocks; b++) {
        log_weights(t, b, i) = 0.0; // reset
        arma::uvec coord = coord_array[b];
        for(int v=0; v < coord.size(); v++) {
          log_weights(t, b, i) += gDensityLog(results_array(t,coord(v),i), Y(t-1,coord(v)), sigmaY(coord(v),coord(v))); //Y(t-1) because the time starts from 0
        }
      }
      if(!resampling) {
        //If there is no resampling scheme the weight is the product of the g times the weights at the previous step
        log_weights.slice(i).row(t) = log_weights.slice(i).row(t) % log_weights.slice(i).row(t-1); //% is element wise multiplication
      }
    }


  }

  arma::cube results_array_res = results_array.subcube(1,0,0,n,dimension-1,N-1); //Remove time 0
  arma::cube log_weights_res = log_weights.subcube(1,0,0,n,number_blocks-1,N-1); //Remove time 0

  return(List::create(Named("filteringParticle") = results_array_res,
                      Named("filteringLogWeights") = log_weights_res));
}

// [[Rcpp::export]]
List sequentialImportanceSampling(int N, int n, List fParams, List gParams) {
  arma::mat A = fParams["A"];
  int dimension = A.n_rows;
  List blocks = List::create(linspace(1,dimension,dimension).t());
  return(blockParticleFilter(N,n,blocks,fParams,gParams,false));
}

// [[Rcpp::export]]
List bootstrapParticleFilter(int N, int n, List fParams, List gParams) {
  arma::mat A = fParams["A"];
  int dimension = A.n_rows;
  List blocks = List::create(linspace(1,dimension,dimension).t());
  return(blockParticleFilter(N,n,blocks,fParams,gParams,true));
}
