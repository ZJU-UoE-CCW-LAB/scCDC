#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
double VectorToEntropy(IntegerVector x) {
  
  IntegerVector col = x;
  IntegerVector col_ = x;

  // obtain the set of numbers which removes duplicates   
  col.sort();
  col.erase(unique(col.begin(), col.end()), col.end());

  // obtain how many times does any number appear in the vector  
  NumericVector out = rep(NumericVector::create(0), col.size());
  
  for (int i = 0; i < col_.size(); i++){
    for (int j = 0; j < col.size(); j++){
      if (col_[i] == col[j]){
        out[j]++;
        break;
      }
    }
  }

  // obtain the total size of the vector  
  int sum = 0;
  for (int i = 0; i < out.size(); i++) {
    sum += out[i];
  }

  // calculate the entropy for the vector  
  double entropy=0;
  for (int i = 0; i < out.size(); i++) {
    entropy = entropy - (out[i] / sum) * log2(out[i] / sum);
  }
  
  return entropy;
}
