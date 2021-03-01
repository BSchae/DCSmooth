#ifndef DCSMOOTH_H
#define DCSMOOTH_H

arma::mat weightMatrix(arma::colvec weights, arma::mat matrix);
arma::mat xMatrix(arma::colvec xVector, int polyOrder);
arma::mat npMatrix(SEXP kernFcnPtr, int p, int n);
arma::vec mWeights(arma::mat npMatrix, arma::vec u, int drv);
int factorialFunction(int value);

arma::vec kernFkt_MW200(arma::vec&, double);
arma::vec kernFkt_MW210(arma::vec&, double);
arma::vec kernFkt_MW220(arma::vec&, double);
arma::vec kernFkt_MW320(arma::vec&, double);
arma::vec kernFkt_MW420(arma::vec&, double);
arma::vec kernFkt_MW422(arma::vec&, double);


#endif