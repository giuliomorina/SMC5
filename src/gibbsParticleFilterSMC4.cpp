// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
using namespace arma;
#include "smc.h"

double fDensityComponent(arma::rowvec& xCurrent, arma::rowvec& xNew, int& component,  arma::mat& A,
                         arma::mat& sigma, arma::uvec& neighbour_vec) {
  //WARNING: Tailored on a particular example
  //X_{t+1} = AX_t + W_t, with W_t ~ N(mean,sigma)
  //Meaning that X_{t+1^^component ~ N((AX_t+mean)^component,sigma^component)

  arma::rowvec newMean = xCurrent*A.t();
  return(dnrmArma(xNew(component), newMean(component), sqrt(sigma(component,component)), true));
}

double fgSampleComponent(arma::rowvec& xCurrent, arma::rowvec Y, int& component, arma::mat& A,
                         arma::mat& sigmaX, arma::mat& sigmaY) {
  //WARNING: Tailored on a particular example
  //Returns p^v*g^v which is still a gaussian with
  //1/sigma^2 = 1/sigma^2_p + 1/sigma^2_g
  //mu = (mu_p/sigma^2_p + mu_g/sigma^2_g)*sigma^2

  //p_v is a Gaussian with
  //mu_p = (AX_t+mean)^component
  //sigma_p = sigmaX^component

  double fg_sigma = sigmaX(component,component)*sigmaY(component,component)/(sigmaX(component,component)+sigmaY(component,component));
  arma::rowvec aux = xCurrent*A.t();
  double fg_mean = fg_sigma * (aux(component)/sigmaX(component,component) + Y(component)/sigmaY(component,component));
  arma::colvec res = rnormArma(1, fg_mean, sqrt(fg_sigma));

  return(res.at(0));

}

double fgScaleFactorComponent(arma::rowvec& xCurrent, arma::rowvec Y, int& component, arma::mat& A,
                              arma::mat& sigmaX, arma::mat& sigmaY) {
  //See fgSampleComponent
  //The scale factor is given by
  //1/sqrt(2*pi*(sigma_f^2+sigma_g^2))*exp(-(mu_f-mu_g)^2/(2*(sigma_f^2+sigma_g^2)))
  //This function returns the log
  arma::rowvec aux = xCurrent*A.t();

  // double res = -0.5*log(2.0 * M_PI * (sigmaX(component,component) + sigmaY(component,component)));
  // res = res - pow(aux(component) - Y(component),2)/(2.0 * (sigmaX(component,component) + sigmaY(component,component)));

  return(dnrmArma(0,aux(component) - Y(component),sqrt(sigmaX(component,component) + sigmaY(component,component)),true));
}

// [[Rcpp::export]]
arma::uvec neighbourSMC4(int component, int radius, int dimension) {
  if(radius <= 0) stop("Radius must be bigger than 0.");
  //We're assuming that the structure is quite easy, that is a line
  //If A has periodic condition it is like a circle.
  int counter = 0;
  if(radius >= ceil(dimension/2.0)) {
    arma::uvec neighbour(dimension);
    for(int i = 0; i<dimension; i++) {
      neighbour(i) = i;
    }
    return(neighbour);
  }

  arma::uvec neighbour(2*radius+1);
  for(int i=1; i<=radius; i++) {
    neighbour(counter) = modulo(component + i,dimension);
    neighbour(counter+1) = modulo(component - i,dimension);
    counter += 2;
  }
  neighbour(counter) = component;
  return(neighbour);

  // if(radius >= ceil(dimension/2.0) && A(0,A.n_cols-1) != 0 && A(A.n_rows-1,0) != 0) {
  //   arma::uvec neighbour(dimension-1);
  //   for(int i = 0; i<dimension; i++) {
  //     if(i != component) {
  //       neighbour(counter) = i;
  //       counter++;
  //     }
  //   }
  //   return(neighbour);
  // }
  //
  // arma::vec neighbour(2*radius);
  // for(int i=1; i<=radius; i++) {
  //   neighbour(counter) = component + i;
  //   neighbour(counter+1) = component - i;
  //   counter += 2;
  // }
  //
  // if(A(0,A.n_cols-1) != 0 && A(A.n_rows-1,0) != 0) {
  //   //Circle line
  //   uvec neighbourRes(2*radius);
  //   for(int i=0; i < neighbour.size(); i++) {
  //     neighbourRes(i) = modulo((int)neighbour(i),dimension);
  //   }
  //   return(neighbourRes);
  // }
  //
  // //Straight line
  // warning("Not homogeneous errors!");
  // uvec neg = find(neighbour < 0);
  // uvec big = find(neighbour >= dimension);
  // int length = 2*radius - neg.size() - big.size();
  // uvec neighbourRes(length);
  // counter = 0;
  // for(int i=0; i < neighbour.size(); i++) {
  //   if(neighbour(i) >= 0 && neighbour(i) < dimension) {
  //     neighbourRes(counter) = (int)neighbour(i);
  //     counter++;
  //   }
  // }
  // return(neighbourRes);

}


// [[Rcpp::export]]
List gibbsParticleFilterSMC4(int N, int n, int m, int radius,
                               List fParams, List gParams) {

  //Unroll fParams and gParams
  arma::mat A = fParams["A"];
  int dimension = A.n_rows;
  arma::rowvec X0 = fParams["X0"];
  arma::mat sigmaX = fParams["sigmaX"];
  arma::mat Y = gParams["Y"];
  //arma::rowvec meanY = gParams["mean"]; //IGNORED
  arma::mat sigmaY = gParams["sigmaY"];
  //Check if sigmaX and sigmaY are diagonal
  bool diagonal = checkDiagonal(sigmaX) && checkDiagonal(sigmaY);
  if(!diagonal) stop("This implementation works only for diagonal sigmaX and sigmaY.");

  //Auxiliary variables
  arma::rowvec zero_vec(dimension, fill::zeros);
  arma::cube results_array(n+1, dimension, N, fill::zeros);
  arma::rowvec R(dimension);
  int random_number;
  srand (time(NULL)); // initialise random seed
  arma::colvec log_weights(N, fill::zeros);
  arma::vec prob;
  int J, w;
  arma::rowvec X_sm1_j(dimension);
  int length = (radius >= ceil(dimension/2.0)) ? dimension-1 : 2*radius;
  arma::uvec neighbour_vec(length);

  //Fill time zero
  for(int p=0; p<N; p++) {
    results_array.slice(p).row(0) = X0;
  }

  for(int t = 1; t<n+1; t++) { //Loop over time
    for(int i=0; i<N; i++) { //Loop over particles (parallelizable)
      //Gibbs Sampler
      //Sample initial points of the Gibbs Sampler
      random_number = rand() % N;
      R = (results_array.slice(random_number)).row(t-1);
      for(int l=0; l<m; l++) { //Loop over sweeps
        for(int k=0; k<dimension; k++) { //Loop over components
          //Compute weights
          log_weights = zeros<arma::colvec>(N); //Reset
          for(int j=0; j<N; j++) {
            //Compute the log product
            neighbour_vec = neighbourSMC4(k,radius,dimension);
            X_sm1_j = results_array.slice(j).row(t-1);
            for(int auxi = 0; auxi<neighbour_vec.size(); auxi++) {
              w = neighbour_vec(auxi);
              if(w != k) log_weights(j) = log_weights(j) + fDensityComponent(X_sm1_j, R, w,  A, sigmaX, neighbour_vec);
            }
            //Add the log scale factor
            log_weights(j) = log_weights(j) + fgScaleFactorComponent(X_sm1_j, Y.row(t-1), k, A, sigmaX, sigmaY); //Y(t-1) because t starts from 1
          }
          //Sample j between 1 and N with probability given by the weights
          prob =  exp(log_weights-max(log_weights));
          J = (int)ProbSampleReplace(N, 1, prob).at(0);

          //Sample R(k)
          X_sm1_j = results_array.slice(J).row(t-1);
          R(k) = fgSampleComponent(X_sm1_j, Y.row(t-1), k, A, sigmaX, sigmaY);
        }
      }
      //Save result from Gibbs
      results_array.slice(i).row(t) = R;
      Rcpp::checkUserInterrupt(); //Check user interruption
    }
  }

  arma::cube log_weights_final(n, 1, N, fill::ones);
  log_weights_final *= log(1.0/N);
  arma::cube results_array_res = results_array.subcube(span(1,n),span(),span());

  return(List::create(Named("filteringParticle") = results_array_res,
                      Named("filteringLogWeights") = log_weights_final));
}
