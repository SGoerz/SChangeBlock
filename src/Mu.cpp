#include <Rcpp.h>
using namespace Rcpp;

//' Mu
//' 
//' This function returns a vector of the block means for a given random field X.
//' 
//' @param X Numeric vector or matrix.
//' @param l block length. Numeric vector of length 1 or 2, depending on the number of dimensions of X.
//' 
//' @return A numeric vector of length \code{floor(n[1] / l[1]) * floor(n[2] / l[2])}.
//' 
//' @examples 
//' X <- genField(c(50, 100), H = 100, type = 2)
//' M <- Mu(X, c(10, 20))
//' 
//' plot(X)
//' image(matrix(M, ncol = 5))
//' 
//' @export
// [[Rcpp::export]]
 NumericVector Mu(NumericMatrix X, IntegerVector l)//, Nullable<IntegerVector> e = R_NilValue)
 {
   int l_n, l_m; //ex, ey;
   
   l_n = l[0]; 
   if(l.length() == 1)
   {
     l_m = 1;
   } else
   {
     l_m = l[1];
   }
   
   // if(e.isNull())
   // {
   //   ex = ey = 0;
   // } else
   // {
   //   IntegerVector e_(e);
   //   ex = e_(0);
   //   
   //   if(e_.length() == 1)
   //   {
   //     ey = 0; 
   //   } else
   //   {
   //     ey = e_(1);
   //   }
   // }
   
   int b_n = floor(X.nrow() / l_n);
   int b_m = floor(X.ncol() / l_m);
   
   NumericVector MuVec(b_n * b_m);
   int i, j, i2, j2, blocklength;
   
   blocklength = l_n * l_m;
   
   // MuVec(0) = 0;
   // count = 0;
   // for(i2 = ex; i2 < l_n; i2++)
   // {
   //   for(j2 = ey; j2 < l_m; j2++)
   //   {
   //     MuVec(0) += X(i2, j2);
   //     count++;
   //   }
   // }
   // MuVec(0) /= count;
   
   for(i = 0; i < b_n; i++)
   {
     for(j = 0; j < b_m; j++)
     {
       MuVec(i * b_m + j) = 0;
       
       for(i2 = i * l_n; i2 < (i + 1) * l_n; i2++)
       {
         for(j2 = j * l_m; j2 < (j + 1) * l_m; j2++)
         {
           MuVec(i * b_m + j) += X(i2, j2);
         }
       }
       
       MuVec(i * b_m + j) /= blocklength;
     }
   }
   
   return MuVec;
 }
