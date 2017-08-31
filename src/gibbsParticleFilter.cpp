// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
#include "smc.h"

double fDensityComponent(arma::cube& results_array,
                         int particleIndCurrent, int particleIndNew, int time,
                         int component,
                         arma::mat& A, arma::mat& sigma) {
  //WARNING: Tailored on a particular example
  //X_{t+1} = AX_t + W_t, with W_t ~ N(mean,sigma)
  //Meaning that X_{t+1^^component ~ N((AX_t+mean)^component,sigma^component)
  arma::rowvec newMean = results_array.slice(particleIndCurrent).row(time-1)*A.t();
  return(dnrmArma(results_array(time,component,particleIndNew), newMean(component), sqrt(sigma(component,component)), true));
}

double fgSampleComponent(arma::cube& results_array, arma::mat& Y,
                         int particleCurrent, int timeX, int timeY,
                         int component, arma::mat& A,
                         arma::mat& sigmaX, arma::mat& sigmaY) {
  //WARNING: Tailored on a particular example
  //Returns p^v*g^v which is still a gaussian with
  //1/sigma^2 = 1/sigma^2_p + 1/sigma^2_g
  //mu = (mu_p/sigma^2_p + mu_g/sigma^2_g)*sigma^2

  //p_v is a Gaussian with
  //mu_p = (AX_t+mean)^component
  //sigma_p = sigmaX^component

  double fg_sigma = sigmaX(component,component)*sigmaY(component,component)/(sigmaX(component,component)+sigmaY(component,component));
  arma::rowvec aux = results_array.slice(particleCurrent).row(timeX)*A.t();
  double fg_mean = fg_sigma * (aux(component)/sigmaX(component,component) + Y(timeY,component)/sigmaY(component,component));
  arma::colvec res = rnormArma(1, fg_mean, sqrt(fg_sigma));

  return(res.at(0));

}


double fgScaleFactorComponent(arma::cube& results_array, arma::mat& Y,
                              int particleInd, int timeX, int timeY,
                              int component, arma::mat& A,
                              arma::mat& sigmaX, arma::mat& sigmaY) {
  arma::rowvec aux = results_array.slice(particleInd).row(timeX)*A.t();
  return(dnrmArma(0,aux(component) - Y(timeY,component),sqrt(sigmaX(component,component) + sigmaY(component,component)),true));
}

arma::uvec neighbour(int component, int radius, int dimension) {
  if(radius < 0) stop("Radius must be bigger than 0.");
  //We're assuming that the structure is quite easy, that is a line
  //If A has periodic condition it is like a circle.
  if(radius == 0) {
    arma::uvec neighbour(1);
    neighbour(0) = component;
    return(neighbour);
  }
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

}


// [[Rcpp::export]]
List gibbsParticleFilterOnline(int N, int n, int m,  int radius,
                               arma::cube& particles,
                               List fParams, List gParams, bool init = false) {
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
  //Set up variables
  int dimension = A.n_rows;
  int time_matrix = n;
  if(init) {
    time_matrix = n+1; //It's starting from time zero so we need one more
  }
  arma::cube results_array(time_matrix, dimension, N, fill::zeros); //n+1 because we start from time zero
  //Fill up with the already provided results
  results_array.subcube(0,0,0,particles.n_rows-1,dimension-1,N-1) = particles;

  //Neighbours
  arma::uvec *neighbour_array; //Contain the neighbour of each component (exluding itself)
  neighbour_array = new arma::uvec[dimension];
  for(int k=0; k<dimension; k++) {
    neighbour_array[k] = neighbour(k,radius,dimension);
  }

  //Auxiliary variables
  arma::mat logScalingFactor(dimension,N);
  arma::uvec randomNumbers(N);
  arma::mat logWeights(dimension,N);
  arma::colvec logWeightsParticles(N);
  arma::vec prob;
  int w, J;

  // arma::rowvec R_ELIMINAMI(dimension);
  // arma::colvec log_weights_ELIMINAMI(N, fill::zeros);
  // arma::rowvec X_sm1_j_ELIMINAMI(dimension);
  // int length = (radius >= ceil(dimension/2.0)) ? dimension-1 : 2*radius;
  // arma::uvec neighbour_vec_ELIMINAMI(length);
  // arma::vec ELIMINAMI(N);
  // arma::vec ELIMINAMI2(N);


  for(int t=particles.n_rows; t < time_matrix; t++) { //Loop over time
    //Sample from initial distribution
    randomNumbers = ProbSampleReplace(N,N,ones<arma::colvec>(N));
    for(int j=0; j<N; j++) results_array.subcube(span(t),span(),span(j)) = results_array.subcube(span(t-1),span(),span(randomNumbers(j)));
    //Compute log scaling factor for the weights
    for(int j=0; j<N; j++) {
      for(int k=0; k<dimension; k++) {
        logScalingFactor(k,j) = init ? fgScaleFactorComponent(results_array, Y, j, t-1, t-1, k, A, sigmaX, sigmaY) : fgScaleFactorComponent(results_array, Y, j, t-1, t, k, A, sigmaX, sigmaY);
      }
    }

    for(int i=0; i<N; i++) { //Loop over particles (parallelizable)
      // R_ELIMINAMI = (results_array.slice(randomNumbers(i))).row(t-1);

      //Compute initial weights
      if(radius > 0) {
        for(int j=0; j<N; j++) {
          for(int w=0; w<dimension; w++) {
            logWeights(w,j) = fDensityComponent(results_array, j, i, t, w, A, sigmaX);
          }
        }
      }

      for(int l=0; l<m; l++) { //Loop over sweeps
        for(int k=0; k<dimension; k++) { //Loop over components
          //Compute weights
          logWeightsParticles = zeros<arma::colvec>(N); //Reset
          if(radius > 0) {
            for(int j=0; j<N; j++) {
              for(int auxi = 0; auxi < neighbour_array[k].size(); auxi++) {
                w = neighbour_array[k].at(auxi);
                if(w!=k) logWeightsParticles(j) += logWeights(w,j);
              }
            }
          }
          logWeightsParticles += logScalingFactor.row(k).t(); //Add log scale factor

          //ELIMINAMI Compute weights
          // log_weights_ELIMINAMI = zeros<arma::colvec>(N); //Reset
          // for(int j=0; j<N; j++) {
          //   //Compute the log product
          //   neighbour_vec_ELIMINAMI = neighbourSMC4(k,radius,dimension);
          //   X_sm1_j_ELIMINAMI = results_array.slice(j).row(t-1);
          //   for(int auxi = 0; auxi<neighbour_vec_ELIMINAMI.size(); auxi++) {
          //     w = neighbour_vec_ELIMINAMI(auxi);
          //     if(w != k) log_weights_ELIMINAMI(j) = log_weights_ELIMINAMI(j) + fDensityComponent(X_sm1_j_ELIMINAMI, R_ELIMINAMI, w,  A, sigmaX, neighbour_vec_ELIMINAMI);
          //     //if(fDensityComponent(X_sm1_j_ELIMINAMI, R_ELIMINAMI, w,  A, sigmaX, neighbour_vec_ELIMINAMI) != logWeights(w,j)) Rcout << "Uh-Oh! w = " << w << " - j = " << j << endl;
          //   }
          //   //Add the log scale factor
          //   log_weights_ELIMINAMI(j) = log_weights_ELIMINAMI(j) + fgScaleFactorComponent(X_sm1_j_ELIMINAMI, Y.row(t-1), k, A, sigmaX, sigmaY); //Y(t-1) because t starts from 1
          //   ELIMINAMI(j) = fgScaleFactorComponent(X_sm1_j_ELIMINAMI, Y.row(t-1), k, A, sigmaX, sigmaY);
          // }
          // Rcout << t << endl << i << endl << l << endl << k << endl;
          // Rcout << logWeightsParticles << endl << "----" << endl;
          // Rcout << log_weights_ELIMINAMI << endl;
          // Rcout << "------" << endl;
          // Rcout << logScalingFactor.row(k) << endl;
          // Rcout << ELIMINAMI << endl;
          // Rcout << "*****" << endl;
          // if(!approx_equal(logWeightsParticles, log_weights_ELIMINAMI, "absdiff", 0.002)) {
          //   Rcout << t << endl << i << endl << l << endl << k << endl;
          //   stop("Error");
          // }
          //if(k == 1) stop("Ciao");
          //ELIMINAMI FINO A QUA

          //Sample j between 1 and N with probability given by the weights
          prob =  exp(logWeightsParticles-max(logWeightsParticles));
          J = (int)ProbSampleReplace(N, 1, prob).at(0);

          //Sample new point
          results_array(t,k,i) = init ? fgSampleComponent(results_array, Y, J, t-1, t-1, k, A, sigmaX, sigmaY) : fgSampleComponent(results_array, Y, J, t-1, t, k, A, sigmaX, sigmaY);
          //ELIMINAMI
          // arma::rowvec X_sm1_j(dimension);
          // X_sm1_j = results_array.slice(J).row(t-1);
          // results_array(t,k,i) = fgSampleComponent(X_sm1_j, Y.row(t-1), k, A, sigmaX, sigmaY);
          // R_ELIMINAMI(k) = results_array(t,k,i);
          //FINO A QUI


          //Update logweights
          if(radius > 0) {
            for(int j=0; j<N; j++) {
              logWeights(k,j) = fDensityComponent(results_array, j, i, t, k, A, sigmaX);
            }
          }

          Rcpp::checkUserInterrupt(); //Check user interruption

        }
      }
    }
  }

  arma::cube log_weights(time_matrix, 1, N, fill::ones);
  log_weights *= log(1.0/N);

  if(init) {
    arma::cube results_array_res = results_array.subcube(1,0,0,n,dimension-1,N-1); //Remove time 0
    arma::cube log_weights_res = log_weights.subcube(1,0,0,n,0,N-1); //Remove time 0
    return(List::create(Named("filteringParticle") = results_array_res,
                        Named("filteringLogWeights") = log_weights_res));
  } else {
    return(List::create(Named("filteringParticle") = results_array,
                        Named("filteringLogWeights") = log_weights));
  }

}

// [[Rcpp::export]]
List gibbsParticleFilter(int N, int n, int m,  int radius,
                         List fParams, List gParams) {
  //Fill the first time with X0 for all particles
  arma::mat A = fParams["A"];
  arma::rowvec X0 = fParams["X0"];
  int dimension = A.n_rows;
  arma::rowvec ones(dimension,fill::ones);
  arma::cube particles(1, dimension, N, fill::zeros);
  for(int p = 0; p < N; p++) {
    particles.slice(p).row(0) = X0;
  }
  return(gibbsParticleFilterOnline(N,n,m,radius,particles,fParams,gParams,true));
}
