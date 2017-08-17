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

#endif
