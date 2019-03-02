

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


//' helper function for the loss computation
//' @useDynLib DIDdesign
//' @keywords internal
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




//' helper function for repeated cross-section data
//' @param outcome outcome vector
//' @param treatment treatment vector
//' @param time_index time id for each unit
//' @param time_unique a vector of unique time index (should be increasing order)
//' @param y1mean a vector of mean outcome for the treated
//' @param y0mean a vector of mean outcome for the control
//' @value a matrix of outcome filled with mean outcome for each group
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List rcs_mean_fill(
  const arma::vec &outcome,
  const arma::ivec &treatment,
  const arma::uvec &time_index,
  const arma::uvec &time_unique,
  const arma::rowvec &y1mean,
  const arma::rowvec &y0mean
) {

  int n_obs = outcome.n_elem;
  arma::mat y_demean = arma::zeros(n_obs, time_unique.n_elem);
  arma::uvec time_match(n_obs);
  for (int i = 0; i < n_obs; i++) {
    int di = treatment(i);
    if (di == 1) {
      y_demean.row(i) = y1mean;
    } else {
      y_demean.row(i) = y0mean;
    }

    // get time index
    arma::uvec tmp = arma::find(time_unique == time_index(i));
    time_match(i) = tmp(0);
    y_demean(i,time_match(i)) = outcome(i);
  }


  return Rcpp::List::create(y_demean, time_match);
}
