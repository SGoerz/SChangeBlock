#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

double kTH(double x)
{
  if(fabs(x) <= 1) return (1 + cos(M_PI * x)) / 2;
  return 0;
}

//' @export
// [[Rcpp::export]]
double gamma(NumericMatrix X, int h1, int h2)
{
  int n = X.nrow();
  int m = X.ncol();
  
  if(h1 < 0 && h2 < 0) 
  {
    h1 = -h1;
    h2 = -h2;
  }
  
  int i, j;
  double res = 0; 
  
  for(i = std::max(0, -h1); i < std::min(n - h1, n); i++)
  {
    for(j = std::max(0, -h2); j < std::min(m - h2, m); j++)
    {
      res += (X(i, j) * X(i + h1, j + h2));
    }
  }
    
  return res / (n * m);
}

//' Long run variance
//' 
//' Estimates the long run variance of a 2-dimensional matrix \code{X} using kernel
//' density estimation and the Tukey-Hanning kernel function.
//' 
//' @param X numeric matrix,
//' @param b numeric vector containing exactly two values: the bandwidths for the row- reps. column-wise estimation.
//' 
//' @details ??????? Tukey-Hanning kernel ????????
//' 
//' @return A numeric value.
//' 
//' @examples
//' X1 <- genField(c(50, 50), Phi = genPhi(1, 0.4))
//' b <- lrvBandwidth(X1, 1/3, 2/3)
//' lrv(X1, b)
//' 
//' Phi <- matrix(c(0.08, 0.1, 0.08, 0.8, 1, 0.8, 0.08, 0.1, 0.08), ncol = 3)
//' X2 <- genField(c(50, 50), Phi = Phi)
//' b <- lrvBandwidth(X2, 1/3, 2/3)
//' lrv(X2, b)
//' 
//' @export
// [[Rcpp::export]]
double lrv(NumericMatrix X, NumericVector b = 1)
{
  NumericMatrix Y = clone(X);
  Y = X - mean(X);
  
  int b_n, b_m; 
  
  b_n = b[0]; 
  if(b.length() == 2)
  {
    b_m = b[1];
  } else
  {
    stop("l has to contain exactly two values!");
  }
  
  double rows = 0, cols = 0, both = 0;
  int h1, h2;
  
  if(b_n > 1 || b_m > 1)
  {
    for(h1 = 1; h1 < b_n; h1++)
    {
      rows += gamma(Y, h1, 0) * kTH(h1 / b_n);
      
      for(h2 = 1; h2 < b_m; h2++)
      {
        both += gamma(Y, h1, h2) * kTH(h1 / b_n) * kTH(h2 / b_n);
      }
    }
    for(h2 = 1; h2 < b_m; h2++)
    {
      cols += gamma(Y, 0, h2) * kTH(h2 / b_m);
    }
  }
  
  double v = gamma(X, 0, 0);
  double res = v + 2 * rows + 2 * cols + 4 * both;
  if(res <= 0) return v;
  return res;
}


//' Autocorrelation matrix
//' 
//' @param X numeric matrix,
//' @param b numeric vector containing exactly two values: the bandwidths for the row- reps. column-wise estimation. 
//' 
//' @export
// [[Rcpp::export]]
NumericMatrix autocov(NumericMatrix X, NumericVector b, int direction = 0)
{
 NumericMatrix Y = clone(X);
 Y = X - mean(X);

 int n = X.nrow();
 int m = X.ncol(); 
 
 if(direction == 1) n = 1; else if(direction == 2) m = 1; 
 int N = n * m; 
 
 int i, j; 
 int b_n, b_m; 
 
 double thresh = R_NegInf;
 
 b_n = b[0]; 
 if(b.length() == 2)
 {
   b_m = b[1];
 } else
 {
   stop("b has to contain exactly two values!");
 }
 
 if(b_n > n + 1 || b_m > m + 1 || (direction == 1 && b_n > 1) || (direction == 2 && b_m > 1))
 {
   stop("Bandwidths cannot be larger than X.");
 }
 
 NumericMatrix A(N, N); 
 double cov; 
 int h1, h2;
 
 // main diagonal: 
 double var = gamma(Y, 0, 0);
 for(i = 0; i < N; i++)
 {
   A(i, i) = var;
 }

 if(b_n > 1 || b_m > 1)
 {
   // x direction
   for(h1 = 1; h1 < b_n; h1++)
   {
     cov = gamma(Y, h1, 0); // * kTH(h1 / b_n);
     if(cov / var >= thresh) 
     {
       for(i = 0; i < N - h1; i++)
       {
         if(i % n < n - h1)
         {
           A(i + h1, i) = cov; 
           A(i, i + h1) = cov;
         }
       }
     }
   }
   
   // y direction
   for(h2 = 1; h2 < b_m; h2++)
   {
     cov = gamma(Y, 0, h2); // * kTH(h2 / b_m);
     
     if(cov / var >= thresh) 
     {
       for(i = 0; i < N - h2 * n; i++)
       {
         A(i + h2 * n, i) = cov; 
         A(i, i + h2 * n) = cov;
       }
     }
     // x and y; only x goes into the negative
     for(h1 = -b_n + 1; h1 < b_n; h1++)
     {
       if(h1 != 0)
       {
         cov = gamma(Y, h1, h2); // * kTH(h1 / b_n) * kTH(h2 / b_n);
         if(cov / var >= thresh) 
         {
           for(j = std::max(0, -h1); j < N - h2 * n - std::max(0, h1); j++)
           {
             if((j - std::max(0, -h1)) % n < n - abs(h1))
             {
               A(j + h1 + h2 * n, j) = cov;
               A(j, j + h1 + h2 * n) = cov;
             }
           }
         }
       }
     }
   }
 }
 
 return A;
}
