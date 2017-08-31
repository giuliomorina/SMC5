#ifndef SMC_HEADER
#define SMC_HEADER

arma::vec dmvnrmArma(arma::mat&, arma::rowvec&, arma::mat&, bool, bool);
arma::vec dmvnrmArma(arma::rowvec&, arma::rowvec&, arma::mat&, bool, bool);
double dnrmArma(double, double, double, bool);
arma::mat mvrnormArma(int, arma::rowvec&, arma::mat&, bool);
arma::colvec rnormArma(int, double, double);
double logAddition(double, double);
double logAdditionSum(arma::vec&);
void standardizeLogWeights(arma::cube&);
void standardizeLogVector(arma::vec& );
arma::uvec ProbSampleReplace(int, int, arma::vec);
bool checkDiagonal(arma::mat&);
bool checkSymmetric(arma::mat&);
int modulo(int,int);
double fDensityComponent(arma::cube&,int, int, int, int, arma::mat&, arma::mat&);
double fgSampleComponent(arma::rowvec& , arma::rowvec , int& , arma::mat& , arma::mat& , arma::mat& );
double fgScaleFactorComponent(arma::rowvec& , arma::rowvec , int& , arma::mat& , arma::mat& , arma::mat& );
arma::uvec neighbourSMC4(int, int, int);
arma::uvec neighbour(int, int, int);
#endif
