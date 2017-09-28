#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat subdistsmult(const Rcpp::List& subcentroids,
                       const Rcpp::List& query_chunks,
                       const int& M,
                       const int& k,
                       const int& cellnum){
  int m,j;
  arma::vec subqm;
  arma::mat dists(M, k);
  for (m = 0; m < M; m++) {
    arma::mat chunk = query_chunks[m];
    arma::vec subqm = chunk.col(cellnum);
    arma::mat subcm = subcentroids[m];
    for (j = 0; j < k; j++) {
      dists(m,j) = arma::dot(subqm, subcm.col(j));
      //dists(m,j) = arma::as_scalar(subqm.t() *subcm.col(j));
    }
  }
  return dists;
}

// [[Rcpp::export]]
Rcpp::List nearcs(const int& w,
                 const int& M,
                 const arma::Col<int> clus_ind, 
                 const arma::mat & subclusters,
                 const arma::mat& dists,
                 const int& len) {
  arma::vec nn(w); 
  Rcpp::NumericVector NN(w);
  int ind, r, m, min_index;
  double celldist;
  for (r = 0; r < len; r++) {
    celldist = 0;
    ind = clus_ind(r);
    
    for (m = 0; m < M; m++) {
      celldist += dists(m, subclusters(m,ind-1));
    }
    
    // double total_mag =0;
    // for (m = 0; m < (M+1); m++) {
    //   total_mag += mags(m,subclusters(m,ind-1));
    // }
    // celldist = celldist/pow(total_mag, 0.5);
    if (r < w) {
      NN(r) = celldist;
      nn(r) = ind;
    } else {
      //double mn = Rcpp::min(NN);
      min_index = Rcpp::which_min(NN);
      if (celldist > NN(min_index)) {
        nn(min_index) = ind;
        NN(min_index) = celldist;
      }
      
    }
  }
  return List::create(_["nn"] = nn,
                      _["distances"]   = NN);
}


// [[Rcpp::export]]
arma::mat normalise(const arma::mat& dat) {
  arma::mat res = arma::normalise(dat);
  return res;
}


// [[Rcpp::export]]
Rcpp::List NNfirst(const int& w,
                const int& k,
                const Rcpp::List& subcentroids,
                const arma::mat& subclusters,
                const Rcpp::List& query_chunks,
                const int& M,
                const arma::vec& SqNorm) {
  
  arma::mat dists;
  //arma::vec nn;
  arma::Col<double> query;
  int cols = subclusters.n_cols;
  int n;
  arma::mat queryc = query_chunks[0];
  int numcells = queryc.n_cols;
  arma::mat best_cells(w, numcells);
  Rcpp::NumericMatrix best_distances(w, numcells);
  arma::mat subqueries;
  arma::vec nn(w); 
  double celldist;
  Rcpp::NumericVector NN(w);
  int ind, r, m, min_index;
  arma::Col<int> clus_ind = arma::linspace<arma::Col<int> >(1, cols, cols);
  for (n = 0; n < numcells; n++) {
    dists = subdistsmult(subcentroids, query_chunks, M, k, n);
    for (r = 0; r < cols; r++) {
      celldist = 0;
      ind = clus_ind(r);
      
      for (m = 0; m < M; m++) {
        celldist += dists(m, subclusters(m,ind-1)-1);
      }
      // only divide through by the Eucl. norm of the query if it is non-zero
      if (SqNorm(n) > 0) {
        celldist = celldist/(sqrt(M)*sqrt(SqNorm(n)));
      }
      if (r < w) {
        NN(r) = celldist;
        nn(r) = ind;
      } else {
        min_index = Rcpp::which_min(NN);
        if (celldist > NN(min_index)) {
          nn(min_index) = ind;
          NN(min_index) = celldist;
        }
        
      }
    }
    best_cells.col(n) = nn;
    best_distances.column(n) = NN;
  }
  arma::mat dataset_inds = arma::ones(w, numcells);
  return List::create(_["cells"] = best_cells,
                      _["distances"]   = best_distances,
                      _["dataset_inds"] = dataset_inds);
}

// [[Rcpp::export]]
Rcpp::List NNmult(const int& w,
                  const int& k,
                  const Rcpp::List& subcentroids,
                  const arma::mat& subclusters,
                  const Rcpp::List& query_chunks,
                  const int& M,
                  const arma::vec& SqNorm,
                  const arma::mat& best_cells_so_far,
                  const Rcpp::NumericMatrix& best_distances_so_far,
                  arma::mat dataset_inds,
                  const int& dat_num) {
  
  arma::mat dists;
  //arma::vec nn;
  arma::Col<double> query;
  int cols = subclusters.n_cols;
  int n;
  arma::mat queryc = query_chunks[0];
  int numcells = queryc.n_cols;
  arma::mat best_cells(w, numcells);
  Rcpp::NumericMatrix best_distances(w, numcells);
  arma::mat subqueries;
  arma::vec nn(w); 
  double celldist;
  Rcpp::NumericVector NN(w);
  double best_poss = 0;
  int ind, r, m, min_index;
  arma::Col<int> clus_ind = arma::linspace<arma::Col<int> >(1, cols, cols);
  for (n = 0; n < numcells; n++) {
    NN = best_distances_so_far.column(n);
    nn = best_cells_so_far.col(n);
    dists = subdistsmult(subcentroids, query_chunks, M, k, n);
    // calculate total distance from best possible comb. of subcentroids
    for (m = 0; m < M; m++) {
      best_poss += arma::max(dists.row(m));
    }
    if (SqNorm(n) > 0) {
      best_poss = best_poss/(sqrt(M)*sqrt(SqNorm(n)));
    }
    // if the best possible isn't better than what we have already, then skip the search
    if (best_poss > min(NN)) {
      for (r = 0; r < cols; r++) {
        celldist = 0;
        ind = clus_ind(r);
        
        for (m = 0; m < M; m++) {
          celldist += dists(m, subclusters(m,ind-1));
        }
        // only divide through by the Eucl. norm of the query if it is non-zero
        if (SqNorm(n) > 0) {
          celldist = celldist/(sqrt(M)*sqrt(SqNorm(n)));
        }
        
        min_index = Rcpp::which_min(NN);
        if (celldist > NN(min_index)) {
          nn(min_index) = ind;
          NN(min_index) = celldist;
          dataset_inds(min_index, n) = dat_num;
        }
          
      }
    } 
    best_cells.col(n) = nn;
    best_distances.column(n) = NN;
  }
  return List::create(_["cells"] = best_cells,
                      _["distances"]   = best_distances,
                      _["dataset_inds"] = dataset_inds);
}

// [[Rcpp::export]]
Rcpp::NumericVector EuclSqNorm(const arma::mat& dat) {
  int cols = dat.n_cols;
  Rcpp::NumericVector res(cols);
  int c;
  for (c = 0; c < cols; c++) {
    res(c) = pow(norm(dat.col(c)),2.0);
  }
  return res;
}


