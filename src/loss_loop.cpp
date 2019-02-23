

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


//' helper function for the loss computation 
//' @useDynLib DIDdesign
// [[Rcpp::export]]
arma::mat loss_loop(
  const arma::mat  &X,
  const arma::vec  &y,
  const arma::vec  &par,
  const arma::uvec &id_subject,
  const arma::uvec &uid,
  const arma::uvec &is_na,
  const int &nobs,
  const int &p
) {
  
  arma::mat Xe(nobs, p);
  for (int i = 0; i < nobs; i++) {
    // subset id_subject 
    arma::uvec use_id  = arma::find(is_na == 0);
    arma::uvec use_id2 = arma::find(id_subject(use_id) == uid(i));
    arma::mat Xi = X.rows(use_id2);
    arma::vec yi = y(use_id2);
    
    Xe.row(i) = arma::trans(Xi.t() * (yi - Xi * par));
  }

  return Xe;
  
}
