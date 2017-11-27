#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Computes the dot product between the subcentroids from the indexed reference and the subvectors of an element of the query dataset. 
//' Returns an M by k matrix. Used as an intermediate step (in NNfirst and NNmult) for calculating
//' an approximation of the cosine similarity between the query and the reference. 
//' 
//' @param subcentroids A list of matrices containing the subcentroids of each chunk.
//' @param query_chunks A list of matrices containing the chunks of the query dataset after it has been split
//' according to the product quantization method
//' @param M An integer specifying the number of chunks
//' @param k An integer specifying the number of subcentroids per chunk
//' @param cellnum An integer specifying the column of the query dataset we wish to consider
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
    }
  }
  return dists;
}

//' Normalises each column of a matrix
//' 
//' @param dat A numerical matrix
// [[Rcpp::export]]
arma::mat normalise(const arma::mat& dat) {
  arma::mat res = arma::normalise(dat);
  return res;
}

//' Main nearest neighbour calculation function. Used on the first reference dataset.
//' Returns a list of three objects:
//' 1) the cell indices of the w nearest neighbours
//' 2) the corresponding approx. cosine similarities
//' 3) the corresponding datasets they came from (in this case, the first dataset only)
//' 
//' @param w An integer specifying the number of nearest neighbours
//' @param k An integer specifying the number of subcentroids for each product quantization chunk
//' @param subcentroids A list of matrices containing the subcentroids of each chunk.
//' @param subclusters A matrix containing the subcentroid assignments of each reference cell. See scf_index.
//' @param query_chunks A list of matrices containing the chunks of the query dataset after it has been split
//' according to the product quantization method
//' @param M An integer specifying the number of chunks
//' @param SqNorm A numerical vector containing the Euclidean Squared Norm of each query cell.
// [[Rcpp::export]]
Rcpp::List NN(const int& w,
                const int& k,
                const Rcpp::List& subcentroids,
                const arma::mat& subclusters,
                const Rcpp::List& query_chunks,
                const int& M,
                const arma::vec& SqNorm) {
  
  // Initialization of variables
  arma::mat dists;
  arma::Col<double> query;
  int cols = subclusters.n_cols; //number of cells in the reference
  int n;
  arma::mat queryc = query_chunks[0];
  int numcells = queryc.n_cols; //number of cells in the query dataset
  arma::mat best_cells(w, numcells);
  Rcpp::NumericMatrix best_distances(w, numcells);
  arma::mat subqueries;
  arma::vec nn(w); 
  double celldist;
  Rcpp::NumericVector NN(w);
  int r, m, min_index;
  
  for (n = 0; n < numcells; n++) {
    dists = subdistsmult(subcentroids, query_chunks, M, k, n); // see description of subdistsmult above
    for (r = 0; r < cols; r++) {
      celldist = 0;
      
      // add up corresponding dot products from each chunk
      for (m = 0; m < M; m++) {
        celldist += dists(m, subclusters(m,r)-1);
      }
      
      // divide the sum by the Euclidean norms of both the query and the concatenated subcentroids (sqrt(M) due to normalization)
      // only divide through by the Eucl. norm of the query if it is non-zero
      if (SqNorm(n) > 0) {
        celldist = celldist/(sqrt(M)*sqrt(SqNorm(n)));
      }
      
      // start by assigning the first w reference cells as nearest neighbours
      // subsequently, update the nearest neighbours as we search through the rest of the reference
      if (r < w) {
        NN(r) = celldist;
        nn(r) = r+1; // "+1" required due to the difference in indexing in R and C++
      } else {
        min_index = Rcpp::which_min(NN);
        if (celldist > NN(min_index)) {
          nn(min_index) = r+1;
          NN(min_index) = celldist;
        }
        
      }
    }
    best_cells.col(n) = nn;
    best_distances.column(n) = NN;
  }
  
  return List::create(_["cells"] = best_cells,
                      _["distances"]   = best_distances);
}

//' The Euclidean Squared Norm of each column of a matrix is computed and the whole result is returned as a vector.
//' Used as part of the approx. calculations of the cosine similarity between the query and the reference.
//' 
//' @param dat A numerical matrix
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


