// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
#include "smc.h"

double fDensityBlock(arma::rowvec xCurrent, arma::rowvec xNew, arma::uvec coord, arma::mat& A,
                     arma::mat& sigma){
  //WARNING: Tailored on a particular example
  //X_{t+1} = AX_t + W_t, with W_t ~ N(mean,sigma)
  //Meaning that X_{t+1} ~ N(AX_t+mean,sigma)
  arma::rowvec newPointBlock = (xNew.elem(coord-1)).t();
  arma::uvec indices =  linspace<uvec>(0, A.n_rows - 1, A.n_rows);
  arma::rowvec meanXBlock = xCurrent*A.submat(coord-1,indices).t();
  arma::mat sigmaReduced = sigma.submat(coord-1,coord-1);

  return(dmvnrmArma(newPointBlock, meanXBlock, sigmaReduced, true, true).at(0,0));
}

arma::colvec fFunctional(arma::cube& results_array, arma::uvec& coord, int& t ,int& n_alpha){
  int v;
  colvec result(results_array.n_slices, fill::zeros);
  double temp1;
  colvec temp2(results_array.n_slices);

  if(t==0){
    return(result);
  }

  for(int i = 0; i < coord.size(); i++){ //summing over dimensions belonging to a given block
    v = coord(i)-1;
    temp1 = results_array.subcube(span(t), span(v), span(n_alpha)).at(0,0,0);
    temp2 = results_array.subcube(span(t-1), span(v), span());
    result += temp1*temp2;
  }
  //Rcout << result << endl;
  return(result); //we return vector of length N
}

// [[Rcpp::export]]
arma::cube blockForwardSmoothingSMC4(arma::cube filteringResults,
                                 arma::cube filteringLogWeights,
                                 List blocks, List fParams) {
  //Prepare output
  int number_blocks = blocks.size();
  int n = filteringResults.n_rows;
  int dimension = filteringResults.n_cols;
  int N = filteringResults.n_slices;
  arma::cube alpha(n, number_blocks, N, fill::zeros);
  //Setup blocks for C++
  arma::uvec *coord_array;
  coord_array = new arma::uvec[number_blocks];
  for(int i=0; i < number_blocks; i++) {
    arma::uvec coord = blocks[i];
    coord_array[i] = coord;
  }
  //Unroll fParams
  arma::mat A = fParams["A"];
  arma::mat sigmaX = fParams["sigmaX"];

  //Main
  for(int block = 0; block < number_blocks; block++){ //Parallelizable
    arma::rowvec x_current(dimension, fill::zeros);
    arma::rowvec x_new(dimension, fill::zeros);
    arma::vec temp_previous_alpha(N, fill::zeros);
    arma::vec temp_values_alpha(N, fill::zeros);
    arma::vec temp_weights_alpha(N, fill::zeros);
    for(int t=1; t < n; t++){ //Skip the first time as alpha is equal to 0
      arma::uvec coord = coord_array[block];
      for(int n_alpha=0; n_alpha < N; n_alpha++){ //Compute alpha of each particle
        x_new = filteringResults.subcube(span(t), span(), span(n_alpha));
        for(int m_alpha=0; m_alpha < N; m_alpha++){
          x_current = filteringResults.subcube(span(t-1), span(), span(m_alpha));
          temp_weights_alpha(m_alpha) = fDensityBlock(x_current, x_new, coord,A, sigmaX) + filteringLogWeights(t-1, block, m_alpha);
        }
        standardizeLogVector(temp_weights_alpha);

        temp_previous_alpha = alpha.subcube(span(t-1),span(block), span());
        temp_values_alpha =  temp_previous_alpha + fFunctional(filteringResults, coord, t, n_alpha);
        temp_values_alpha =  exp(temp_weights_alpha) % temp_values_alpha;
        alpha(t, block, n_alpha) = sum(temp_values_alpha);
      }
    }
  }


  return(alpha);
}
