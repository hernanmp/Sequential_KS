



#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector empirical_cdf_cpp(NumericVector  x)
{
  int m = x.size();
  NumericVector ecd(m,0.0); 
  
  double N;
  double  sum_x  = 0.0;
  
  for(int j = 0; j < m; j++)
  {
    sum_x = sum_x + x[j];
  }
  N =  sum_x;
  
  for(int j = 0; j<m; j++)
  {
    if(j == 0)
    {
      ecd[0] = x[0]/N;
    }
    if(j > 0)
    {
      ecd[j] = ecd[j-1] + x[j]/N;
    }
  }
  
  return(ecd);
}


// [[Rcpp::export]]

double kolmogorov_statistic_cpp(NumericVector x, NumericVector F0)
{
   int m = x.size();
   NumericVector ecd(m); 
   ecd = empirical_cdf_cpp(x);
   
   double sum_x  = 0.0;
   double max_diff = -pow(10,9);
   double aux = 0.0;
   for(int j = 0; j <m;j++)
   {
     sum_x = sum_x + x[j];
     
     aux = F0[j] - ecd[j];
     if(aux <0)
     {
       aux = -aux;
     }
     if(aux >  max_diff)
     {
       max_diff = aux;
     }
     
   }
   
   
  return( sqrt(sum_x)*max_diff  );
} 



// [[Rcpp::export]]
NumericVector quantileCpp(NumericVector x, NumericVector probs) {
  Environment stats("package:stats");
  Function quantile = stats["quantile"];
  int npr = probs.size();
  NumericVector ans(npr);
  for(int i=0; i<npr; i++){
    ans[i] = as<double>(quantile(x, probs[i]));
  }
return ans;
}

