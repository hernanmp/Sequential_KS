library(Rcpp)
sourceCpp("functions.cpp")



##  For each location generates, with probability prop, a uniform draw from 0:M. With probability 1-prop the
## the rate in a given location is 1-prop
generate_sparse_uniform_data = function(M,prop)
{
  
  ind0  = rbinom(m, size=1,prob= prop)
  
  lambda0 =  M*runif(m)*ind0 + 10^-3
  theta0 = lambda0  / sum(lambda0)
  F0 = rep(0,m)  ## empirical known distribution
  
  F0 = cumsum(theta0)
  
  return( list(theta = theta0, lambda  = lambda0, F = F0 ))
}




############Generate data from DP prior
## a shape, b rate. M is a normalization constant
generate_data_from_DP_prior = function(a,b,M)
{

   lambda0 =  rgamma(m,a,b)/M 
   theta0 = lambda0  / sum(lambda0)
    F0 = rep(0,m)  ## empirical known distribution
  # 
   F0 = cumsum(theta0)

   return( list(theta = theta0, lambda  = lambda0, F = F0 ))
}


###########################################################3333333
###########  Generate data from normal distribution

generate_normal_data = function(mu,sigma,M,xlim,ylim)
{
  loc =  seq(xlim,ylim,length = m)
  lambda0 =  M*dnorm(loc,mu,sigma)  
  theta0 = lambda0  / sum(lambda0)
  F0 = rep(0,m)  ## empirical known distribution
  
  F0 = cumsum(theta0)
  
  
  
  return( list(theta = theta0, lambda  = lambda0, F = F0 ))
}

################################################################3

gridprep = function(y, gridsize=250) {
  n = length(y)
  
  # Set up the discrete grid
  ymin = min(y) - 1
  ymax = max(y) + 1
  mybreaks = seq(ymin, ymax, length=gridsize+1)
  k = length(mybreaks)
  mids = (mybreaks[2:k] + mybreaks[1:(k-1)])/2
  
  yhist = hist(y, mybreaks,plot=FALSE)
  ycounts = yhist$counts
  # G = outer(mids,mids, FUN=function(x,y){dnorm(x-y)*diff(mybreaks)})
  
  list(mids=mids, ycounts=ycounts, breaks=mybreaks)
}


#####################################################################

grid_counts = function( x, d  )
{ 
  if( log(length(x))/log(2)  !=  d )
  { 
    return (NULL);
  }
  
  n = sum(x)
  
  counts =  rep(0,2^{d+1}-1)
  
  counts[1] = n
  cc = 2
  
  for(j  in 2: (d+1))
  {
    num = 2^{j-1} 
    
    counts[cc] = sum(x[1: (length(x)/num)])
    last_ind = length(x)/num
    cc = cc+ 1
    
    for(i in 2: num)
    {
      counts[cc] = sum(x[(last_ind +1): (last_ind + length(x)/num) ])
      cc = cc+ 1
      last_ind =  last_ind + length(x)/num
    }
  }
  
  return( counts)
}
##################33##################33##################33##################33##################33
##################33##################33##################33##################33
##################33v##################33##################33##################33v


log_poyla_tree_prob = function(ycounts, alpha,d)
{
  bf = 0
#   ycounts2 = ycounts
#   ind = which(ycounts == 0)
#   ycounts2[ind] = 1 
#   bf = sum(log(1:ycounts2[1] + 10^-9))
#   
#   cc = 2^d 
#   for(j in 1:2^d)
#   {
#     bf = bf - sum(log( 1:ycounts2[cc]  + 10^-9))
#     cc = cc + 1
#   }
#   #####
  for(j in 1:(2^d-1))
  {
     bf =  bf  + lbeta( ycounts[2*j] + alpha[2*j] + 10^{-99}  ,ycounts[2*j+1] + alpha[2*j+1] + 10^{-99}  )
     bf =  bf -  lbeta(alpha[2*j]+ 10^{-99}  , alpha[2*j+1] + 10^{-99}  )
  }
  
  return(bf)
}

################################################

empirical_cdf = function(x)
{
  ecd = rep(0,length(x))
  N =  sum(x)
  
  return(cumsum(x/N))
}


#############333
kolmogorov_statistic = function(x,F0)
{
  ecd = empirical_cdf(x)
  
  return( sqrt(sum(x))*max(abs(F0 -ecd)) )
} 
##########################################################################################33
######################################################################################333
##
###
## function to genare 1 -  level confident bands for poyla tree with parameters alpha
## N is the number of samples used to build the confidence bands
##  d satisfies 2^d = m where m is the length of the count vectorsat each time point

generate_poyla_tree = function(alpha,N,m,d,level)
{
    samples_cdf = matrix(0, nrow = N, ncol =  m)
    alpha[1] = 1
    
    
    for(iter in 1:N)
    {
      p = alpha
      for(j in 1:(2^d-1))
      {
        temp = rbeta( 1,  alpha[2*j], alpha[2*j+1] )
        p[2*j] = temp*p[j]
        p[2*j+1] =  (1-temp)*p[j]
      }
      p = p[2^d:(2^{d+1}-1)]
      p = cumsum(p)
      samples_cdf[iter,] = p
    } 
     posterior_lower = rep(0,m)
     posterior_upper = rep(0,m)  
     posterior_mean = rep(0,m)
     
     
     for(j in 1:m)
     {
       posterior_lower[j] = quantile(samples_cdf[,j],prob= level/2)
       posterior_upper[j] = quantile(samples_cdf[,j],prob= 1-level/2)
       posterior_mean[j] = mean(samples_cdf[,j])
     }
     
     return( list(lower = posterior_lower,  mean = posterior_mean , upper = posterior_upper ))
}

################# update alpha

update_alpha = function(ycounts,alpha,d)
{
  
  for(i in 1:(length(alpha)/2 -1))
  {
    alpha[2*i]  =  alpha[2*i] + ycounts[2*i]
    alpha[2*i+1]  =  alpha[2*i+1] + ycounts[2*i+1]
  }
  return(alpha)
}
##########################################################################################33
######################################################################################333
##
###
## function to genare 1 -  level confident bands for poyla tree with parameters alpha
## N is the number of samples used to build the confidence bands
##  d satisfies 2^d = m where m is the length of the count vectorsat each time point

#  generate_poyla_tree_KSstatistic(alpha,kappa*alpha,theta0,N,m,d,.75,5000)

# alphaN = alpha
# alpha = kappa*alpha

generate_poyla_tree_KSstatistic = function(alphaN,alpha,theta0,N,m,d,level,max_rate)
{
 # max_rate = 4000
  
  alphaN[1] = 1
  alpha[1] = 1  
  
  posterior_mean = cumsum(theta0)

  Ks_statistic = rep(0,floor(N/2))
    
  for(iter in 1:floor(N/2))
  {
    p = alpha
    for(j in 1:(2^d-1))
    {
      temp = rbeta( 1,  alpha[2*j]+ 10^{-99}, alpha[2*j+1] + 10^{-99})
      p[2*j] = temp*p[j]
      p[2*j+1] =  (1-temp)*p[j]
    }
    p = p[2^d:(2^{d+1}-1)]
    
    #total_rate = max(runif(1)*max_rate, 2*m)
  #  x = rpois(length(p),p*max_rate)
      x = rmultinom(1, max_rate, p )
  #  x = rmultinom(1, max_rate, theta0 )
      #rpois(length(p),p*total_rate)
  
    Ks_statistic[iter] = kolmogorov_statistic(x,posterior_mean)    
  } 

  return(list(tresd_hold =  quantile(Ks_statistic,prob= 1-level)))
}
#####################################################################################################
####################################################################################################
##################################################

Training_period = function(ntraining = 500,m,d,lambda0) 
{
  temp = rep(0,m )
  alpha = rep(1,2^{d+1}-1)
  alpha[1] = 1
  
  training_counts = matrix(0, ntraining, m)
  total_rate = 0
  
  for(t in 1: ntraining )
  {
  #  k = rgamma(length(lambda0),lambda0,rep(1,length(lambda0))  )
   # y = rpois(length(lambda0),k)   
    y =rpois(length(lambda0),lambda0)
    #y = rmultinom(1, , prob = lambda0/sum( lambda0))
    temp = temp +  y/ ntraining
    training_counts[t,] = y
    ycounts =  grid_counts(y,d)
    total_rate = total_rate + sum(y)/ntraining
  
    for(i in 1:(length(alpha)/2 -1))
    {
      alpha[2*i]  =  alpha[2*i] + ycounts[2*i]
      alpha[2*i+1]  =  alpha[2*i+1] + ycounts[2*i+1]
    }
   }

   theta0 = temp / sum(temp)
   dir_par =  1+ ntraining*temp  
   F0 = rep(0,m) 
   F0 = cumsum(theta0)
  
  p = alpha
  for(j in 1:(2^d-1))
  {
    temp =  alpha[2*j]/( alpha[2*j] + alpha[2*j+1])# rbeta( 1,  alpha[2*j], alpha[2*j+1] )
    p[2*j] = temp*p[j]
    p[2*j+1] =  (1-temp)*p[j]
  }
  #p = p[2^d:(2^{d+1}-1)]
 
  
  return(list(  F0_hat = F0, alpha = alpha,p =p,theta0 = theta0,dir_par = dir_par,total_rate = total_rate))
}


########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################


sequential_KS_statistics = function(x,t,F0,L)
{
  delta = rep(0,t)
  
  if(t ==1)
  {
    x_prime = x[1,]
    delta[1] = kolmogorov_statistic(x_prime,F0) 
  }
  if(t > 1)
  {
    x_prime = colSums(x[max(1,t-L):t,])
    delta[1] = kolmogorov_statistic(x_prime,F0)      
    
    
    for(s in max(1,t-L):(t-1))
    {
      x_prime = x_prime - x[s,]
      delta[s+1] = kolmogorov_statistic(x_prime,F0)  
    }
    
   }## case t>1
  
  return(max(delta))
}
  
########################################################################################################
########################################################################################################

kolmogorv_distribution_cdf = function(x,N)
{
    k = 1:N
    
    a  = sqrt(2*3.1416)/x*sum(  exp(  -(3.1416*3.1416/(8*x*x))*(2*k-1)^2) )
    return(a)
}
########################################################################################################
########################################################################################################

log_poyla_tree_prob2 = function(ycounts, alpha,d,Num_simulations)
{
   aux = 0
   p = alpha
   
   for(j in 1:Num_simulations)
   {
     for(j in 1:(2^d-1))
     {
       temp =  rbeta( 1,  alpha[2*j], alpha[2*j+1] )
       p[2*j] = temp*p[j]
       p[2*j+1] =  (1-temp)*p[j]
     }
     p = p[2^d:(2^{d+1}-1)]
    
     aux = sum( ycounts*log(p + 10^-40) )     
   }
   return(aux)
}
########################################################################################################
########################################################################################################

log_prob_p_given_alpha  = function(p,alpha,d)
{
    temp = 0
   for(j in 1:(2^d-1))
   {
     temp = temp + dbeta( p[2*j],  alpha[2*j], alpha[2*j+1] ,log=TRUE)
  #   p[2*j] = temp*p[j]
  #   p[2*j+1] =  (1-temp)*p[j]
   }
  return(temp)
}  
  
## log_prob_p_given_alpha(pc,alpha,d) - log_prob_p_given_alpha(pc,kappa*alpha,d)  

#########################################################

poylaTreeBF  = function(x,t,l ,alpha,kappa)
{
  alpha2  = rep(1,length(alpha))
  s = max(t-l,1)
  kappa_alpha = kappa*alpha
  
  log_bayes_factor = rep(0,t-s)
  cumulative_counts = rep(0,length(x))
 # ycounts = rep(0# grid_counts(x[tt,],d)
  ind = 1
  for(tt in  t:s )
  {    
       cumulative_counts = cumulative_counts + x[tt,]
       ycounts =  gridCounts(cumulative_counts ,d) 
       log_bayes_factor[ind] = logPoylaTreeProb(ycounts, alpha,d) - logPoylaTreeProb(ycounts, kappa_alpha,d)
      ind = ind +1
  }
   
  return( min(log_bayes_factor[1:(ind-1)]) )
}

################################################################################
################################################################################
################################################################################

#### Dirichlet chossing kappa


confidence_level_dirichlet_prior = function(dir_par,F0,N,m,level,total_rate)
{ 

  Ks_statistic = rep(0,floor(N/2))
  
  for(iter in 1:floor(N/2))
  {
    p = rgamma(length(dir_par),dir_par,rep(1,length(dir_par)))
    p = p/sum(p)
    
  #  total_rate = max(runif(1)*max_rate, 2*m)
   # x = rpois(length(p),p*total_rate)
    x = rmultinom(1, floor(total_rate), p )
    #rmultinom(1, floor(total_rate), p )
    #rpois(length(p),p*total_rate)
    
    Ks_statistic[iter] = kolmogorov_statistic(x,F0)  
  }
  return(list(tresd_hold =  quantile(Ks_statistic,prob= 1-level)))
}
  

#dir_par = Dir_par
#max_rate = total_rate
#level = level_grid_dir

choosing_kappa_dirichlet_prior = function(dir_par,kappa_grid,F0,N,m,level,max_rate)
{
  
   range = matrix(0,length(kappa_grid),length(level))
    
   for(ind  in 1:length(kappa_grid))
   {
     kappa = kappa_grid[ind]
     range[ind,] = confidence_level_dirichlet_prior(kappa*dir_par,F0,N,m,level,max_rate)$tresd_hold     
   }
   max_range = 1.36
   
   kappa = rep(0,length(level))
   for(j in 1:length(level))
   {
     ind = which((range[,j] - max_range)>0)
     best_ind = which.min(range[ind,j] - max_range)  
     ind  = ind[best_ind]
     kappa[j] = kappa_grid[ind]
   }
   
    
   return(kappa)
}

log_prob_factor_dirichlet_prior = function(dir_par,x,m)
{
   aux = lgamma(sum(dir_par)) - lgamma(sum(x)+sum(dir_par))
   aux = aux + sum(lgamma(x + dir_par)) - sum(lgamma(dir_par))
   return(aux)
}

BF_dirichlet_prior = function(kappa,dir_par,x,m,l,t)
{
  s = max(t-l,1)
  
  log_bayes_factor = rep(0,t-s)
  # ycounts = rep(0# grid_counts(x[tt,],d)
  ind = 1
  for(tt in  t:s )
  {    
    if(tt< t)
    {
      ycounts =  colSums(x[tt:t,]) 
    }
    else{
      ycounts =  x[t,]
    }
    log_bayes_factor[ind] = log_prob_factor_dirichlet_prior(dir_par,ycounts,m) - log_prob_factor_dirichlet_prior(kappa*dir_par,ycounts,m)
    ind = ind +1
  }
  
  return( min(log_bayes_factor[1:(ind-1)]) )
}

##############

##########
cppFunction(' NumericVector gridCounts(NumericVector x, int d) {
		int m = pow(2,d);
    int n  = 0;
    int frac;  
    int cc;
    int num, last_ind;


    for(int i = 0; i <m; i++)
    {
      n = n+ x[i];
    }
   
    NumericVector  counts(  2*m - 1 , 0.0);

    counts[0] = n;
    cc = 2;


    for(int j =  2; j < (d+2); j++)
    {
        num = pow(2,j-1);
        frac = int(m/num);
       
        counts[cc-1] = 0.0;
         
        for(int k = 1; k < (frac+1); k++ ) 
        {
           counts[cc-1] = counts[cc-1] + x[k-1];
        }

        last_ind = int(m/num);
        cc = cc+ 1;

        for(int i  = 2; i< (num+1); i++)
        {
           counts[cc-1] = 0.0;

           for(int l = (last_ind +1); l <  (last_ind + frac +1); l++)
           {
                 counts[cc-1] = counts[cc-1] + x[l-1];
           }

           cc = cc+ 1;
           last_ind =  last_ind + frac;
        }
    }
    return counts;
}')

############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################

# NumericVector  ycounts2( pow(2,d+1)-1,1);
# 
# for(int j = 1; j < pow(2,d+1); j++)
# {
#   if( ycounts[j-1] !=  0)
#   {
#     ycounts2[j-1]  =  ycounts[j-1]; 
#   }
# }


# for(int j = 1; j < int(ycounts2[0]+1); j++)
# {
#   bf = bf + log(j +  pow(10,-9)  );
# }
# 
# int cc = int(pow(2,d)); 
# for(int j = 1;   j < int(pow(2,d)+1);j++)
# {
#   
#   for(int k = 1; k < int(ycounts2[cc-1]+1); k++)
#   {
#     bf = bf - log(k+ pow(10,-9));
#   }
#   
#   cc = cc + 1;
# }    



cppFunction(' double logPoylaTreeProb(NumericVector ycounts, NumericVector alpha,int d) {

  double bf = 0.0;
   

  double aux;
  for(int j = 1;  j <  int(pow(2,d)); j++)
  {
    aux = R::lbeta( ycounts[2*j-1] + alpha[2*j-1] + pow(10,-99) ,ycounts[2*j] + alpha[2*j] + pow(10,-99) );
    bf =  bf  +  aux;
    aux =   R::lbeta(alpha[2*j-1]+ pow(10,-99), alpha[2*j] + pow(10,-99));
    bf =  bf -  aux;
  } 

  return bf; 
}')


#system.time( log_poyla_tree_prob(ycounts, alpha,d))
#system.time(logPoylaTreeProb(ycounts, alpha,d))

########################################################################
##################################################################################################################




#poylaTreeBF  = function(x,t,l ,alpha,kappa)
# {
#   alpha2  = rep(1,length(alpha))
#   s = max(t-l,1)
#   
#   log_bayes_factor = rep(0,t-s)
#   # ycounts = rep(0# grid_counts(x[tt,],d)
#   ind = 1
#   for(tt in  t:s )
#   {    
#     if(tt< t)
#     {
#       ycounts =  grid_counts(floor(colSums(x[tt:t,]))  ,d)  
#     }
#     else{
#       ycounts =  grid_counts(x[t,]  ,d) 
#     }
#     log_bayes_factor[ind] = log_poyla_tree_prob(ycounts, alpha,d) - log_poyla_tree_prob(ycounts, kappa*alpha,d)
#     ind = ind +1
#   }
#   
#   return( min(log_bayes_factor[1:(ind-1)]) )
# }


# 
# 
# cppFunction('
# 
#  double logPoylaTreeProb2(NumericVector ycounts, NumericVector alpha,int d) {
# 
#   NumericVector  ycounts2( pow(2,d+1)-1,1);
#             
#             for(int j = 1; j < pow(2,d+1); j++)
#             {
#                  if( ycounts[j-1] !=  0)
#                  {
#                    ycounts2[j-1]  =  ycounts[j-1];  
#                  }
#             }
#             
#             double bf = 0.0;
#             
#             for(int j = 1; j < int(ycounts2[0]+1); j++)
#             {
#                  bf = bf + log(j +  pow(10,-9)  );
#             }
#             
#             int cc = int(pow(2,d)); 
#             for(int j = 1;   j < int(pow(2,d)+1);j++)
#             {
#             
#                  for(int k = 1; k < int(ycounts2[cc-1]+1); k++)
#                  {
#                      bf = bf - log(k+ pow(10,-9));
#                  }
#             
#                   cc = cc + 1;
#              }    
#             
#              double aux;
#              for(int j = 1;  j <  int(pow(2,d)); j++) 
#              {
#                   aux = R::lbeta( ycounts[2*j-1] + alpha[2*j-1] + pow(10,-9) ,ycounts[2*j] + alpha[2*j] + pow(10,-9) );
#                   bf =  bf  +  aux;
#                   aux =   R::lbeta(alpha[2*j-1]+ pow(10,-9), alpha[2*j] + pow(10,-9));
#                   bf =  bf -  aux;
#              } 
#             
#              return bf; 
#             }
# 
# 
# NumericVector gridCounts2(NumericVector x, int d) {
#     int m = pow(2,d);
#     int n  = 0;
#             int frac;  
#             int cc;
#             int num, last_ind;
#             
#             
#             for(int i = 0; i <m; i++)
#             {
#               n = n+ x[i];
#             }
#             
#             NumericVector  counts(  2*m - 1 , 0.0);
#             
#             counts[0] = n;
#             cc = 2;
#             
#             
#             for(int j =  2; j < (d+2); j++)
#             {
#                 num = pow(2,j-1);
#                 frac = int(m/num);
#             
#                 counts[cc-1] = 0.0;
#             
#                 for(int k = 1; k < (frac+1); k++ ) 
#                 {
#                    counts[cc-1] = counts[cc-1] + x[k-1];
#                 }
#             
#                  last_ind = int(m/num);
#                  cc = cc+ 1;
#             
#                  for(int i  = 2; i< (num+1); i++)
#                  {
#                      counts[cc-1] = 0.0;
#             
#                      for(int l = (last_ind +1); l <  (last_ind + frac +1); l++)
#                      {
#                        counts[cc-1] = counts[cc-1] + x[l-1];
#                      }
#             
#                     cc = cc+ 1;
#                     last_ind =  last_ind + frac;
#                   }
#             }
#             return counts;
#             }
# 
# double poylaTreeBF_rcpp(NumericMatrix x,int t,int l , NumericVector alpha, NumericVector kappa_alpha, int d){
# 
#   NumericVector  alpha2( pow(2,d+1)-1,1);
#   NumericVector   cumulative_counts(pow(2,d),0);
#   NumericVector ycounts( pow(2,d+1)-1,0);
#   int m = pow(2,d) ;
#   double min_value = 10^9;
#    
#     int s = t-l;
#     if(s< 1)
#     {
#        s = 1;
#     }
#   
#     NumericVector log_bayes_factor(t-s,0.0);
#      int ind = 1;
#      int tt;
# 
# 
#     ycounts =  gridCounts( cumulative_counts,d); 
#     return  min_value;
#      
#      for(tt = t;  tt> s-1; tt-- )
#      {    
# 
#        for(int j = 1; j<(m+1); j++)
#        {
#           cumulative_counts[j-1] = cumulative_counts[j-1] + x(tt-1,j-1);
#        }
#       
#        ycounts =  gridCounts2( cumulative_counts,d); 
#        log_bayes_factor[ind-1] = logPoylaTreeProb2(ycounts, alpha,d) - logPoylaTreeProb2(ycounts, kappa_alpha,d);
#        ind = ind +1;
#      }
# 
# 
#       
#       for(int j = 0; j < ind; j++)
#       {
#           if(log_bayes_factor[j] < min_value ) 
#           {
#              min_value = log_bayes_factor[j] ;
#           }
#       }
#   
#   
# 
# }')
# 

#poylaTreeBF_rcpp(x,t,l ,alpha, kappa*alpha, d) 
  

####################################################################################
####################################################################################
####################################################################################
####################################################################################


logRnf = function(t,beta,alpha,x,L,lambda0)
{
  m = length(lambda0)
  if(t == 1)
  {
    log_Rnf = alpha*m*log(beta  +  10^-99) - m*lgamma(alpha   +  10^-99) +  sum( lgamma(x[1,] + alpha  +  10^-99  ))  
    log_Rnf = log_Rnf  - (alpha + sum(x[1,]))*log(1+beta) + sum(lambda0)  - sum(x[1,]*log(lambda0 +  +  10^-99))
    
  }##### first case
  
  if(t > 1)
  {
    
    ind = max(t-L ,1):t
    log_Rnf = rep(0,length(ind))
    
    xprime = colSums(x[max(1,t-L):t,])
    
    ## first 
    log_Rnf[1] = alpha*m*log(beta  +  10^-99) - m*lgamma(alpha +  10^-99) +  sum( lgamma(xprime + alpha   +  10^-99 ))  
    log_Rnf[1]  = log_Rnf[1]   - (alpha + sum(xprime ))*log( t - max(1,t-L)  +  1+beta  +  10^-99) 
    log_Rnf[1]  =  log_Rnf[1] +  ( t - max(1,t-L)  +  1)*sum(lambda0)  - sum(xprime*log(lambda0 + 10^-99))
    
    cc = 2
    for(s in max(1,t-L):(t-1))
    {
      xprime = xprime - x[s,]
      
      log_Rnf[cc] = alpha*m*log(beta  +  10^-99) - m*lgamma(alpha +  10^-99) +  sum( lgamma(xprime + alpha   +  10^-99))  
      log_Rnf[cc]  = log_Rnf[cc]   - (alpha + sum(xprime ))*log( t-(s+1)  + 1+beta) + ( t-(s+1) + 1)*sum(lambda0)  - sum(xprime*log(lambda0 +  10^-99))
      cc = 1+ cc
    }  
    
    log_Rnf_max = max(log_Rnf)
    
    log_Rnf = log(sum(exp(log_Rnf -log_Rnf_max   )     )) +   log_Rnf_max
    
  }## closing general case
  
  return(log_Rnf)
}

###############################################################################################
##################################################################################################
###############################################################################################
##############################################################################################



poylaTreeBF_d  = function(x,t,l ,alpha,kappa)
{
  s = max(t-l,1)
  kappa_alpha = alpha
  
  for(j in 1:d)
  {
    kappa_alpha[(2^j):(2^{j+1}-1)    ]  = kappa[j]*alpha[(2^j):(2^{j+1}-1)] 
  }
  
  log_bayes_factor = rep(0,t-s)
  cumulative_counts = rep(0,length(x))
  # ycounts = rep(0# grid_counts(x[tt,],d)
  ind = 1
  for(tt in  t:s )
  {    
    cumulative_counts = cumulative_counts + x[tt,]
    ycounts =  gridCounts(cumulative_counts ,d) 
    log_bayes_factor[ind] = logPoylaTreeProb(ycounts, alpha,d) - logPoylaTreeProb(ycounts, kappa_alpha,d)
    ind = ind +1
  }
  
  return( min(log_bayes_factor[1:(ind-1)]) )
}
####################################################################################
############################################################################################
###########################################################################################33


generate_poyla_tree_KSstatistic_d = function(alpha2,theta0,N,m,d,level,max_rate,ind,current_kappa,posterior_mean)
{
  
  alpha2[2^{ind} : (2^{ind+1}-1)] = current_kappa*alpha2[2^{ind} : (2^{ind+1}-1)]
  alpha2[1] = 1
  
  if(ind == 1)
  {
    posterior_mean = rep(1,2)
    posterior_mean[1] = sum(theta0[1:(length(theta0)/2)])
    
    p = rep(0,2)
    Ks_statistic = rep(0,N)
    for(iter in 1:N)
    {
      temp =   rbeta(1, alpha2[2],alpha2[3])
      p[1] = temp
      p[2] = 1- temp
      
      x = rmultinom(1, max_rate, p )
      Ks_statistic[iter] = kolmogorov_statistic(x,posterior_mean) 
    }
  }
  if(ind > 1)
  {
    
    Ks_statistic = rep(0,N)
    
    for(iter in 1:N)
    {      
#       p = alpha2
#       
#       for(j in 1:(2^{d}-1))
#       {
#         temp =  rbeta(1, alpha2[2*j]+ 10^{-99} , alpha2[2*j+1]+ 10^{-99} )
#         p[2*j] = temp*p[j]
#         p[2*j+1] =  (1-temp)*p[j]
#         
#         if(j > 2^{ind+1})
#           break;
#       }
#       p = p[2^{ind}:(2^{ind+1}-1)]
#       
#       x = rmultinom(1, max_rate, p )
#       Ks_statistic[iter] = kolmogorov_statistic(x,posterior_mean) 
      #######################################################33333
      p = alpha2
      counts_p = alpha2
      counts_p[1] =  rpois(1,max_rate)#rbinom(1,2,.2)

      for(j in 1:(2^{d}-1))
      {
        temp =  rbeta(1, alpha2[2*j]+ 10^{-99} , alpha2[2*j+1]+ 10^{-99} )
        p[2*j] = temp*p[j]
        p[2*j+1] =  (1-temp)*p[j]
        counts_p[2*j] = rbinom(1, counts_p[j],p[2*j])
        counts_p[2*j+1] = counts_p[j] - counts_p[2*j]
 
        if(j > 2^{ind+1})
          break;
      }
      ###3
      x = counts_p[2^{ind}:(2^{ind+1}-1)]
      Ks_statistic[iter] = kolmogorov_statistic(x,posterior_mean) 
    }## close iter 
    
  }
  ################
  return(list(tresd_hold =  quantile(Ks_statistic,prob= 1-level), KS =       Ks_statistic))
}
#####################################################################################33
#######################################################################################

#finding_kappa_alpha(level_grid,theta0_hat,alpha,d,m,total_rate)
finding_kappa_alpha = function(level_grid,theta0,alpha,d,m,total_rate)
{
  kappa_d =  matrix(0,d,length(level_grid))  #a
  kappa_d_mean = kappa_d
  posterior_mean = cumsum(theta0)
  kappa_alpha = alpha
  
  theta0_dist = gridCounts(theta0,d)
  kappa_alpha = alpha
  
  Niter = 60
  for(iter in 1:Niter)
  {
    for(j in 1:length(level_grid))
    {
      kappa_alpha = alpha
      for(ind in 1:d)
      {
        posterior_mean_ind = cumsum(theta0_dist[2^{ind}:(2^{ind+1}-1)])
        kappa_grid = 10^seq(-5,-1,length= 20)
        Ah_range = rep(0,length(kappa_grid))
        max_range =  1.36
        
        for(i in 1:length(kappa_grid))
        {
          #temp = generate_poyla_tree_KSstatistic_d(kappa_alpha,theta0,N= 100,m,d,level_grid[j],total_rate,ind,kappa_grid[i] ,posterior_mean_ind)
           # 100
          KS = generate_poyla_tree_KSstatistic_d_cpp(kappa_alpha,theta0,N= 100,m,d,level_grid[j],total_rate,ind,kappa_grid[i] ,posterior_mean_ind)
          Ah_range[i] =  quantile(KS,prob= 1-level_grid[j])
          #Ah_range[i] =  temp$tresd_hold
         
        }###  iter = 1
        
        
        index= which.min(    abs(Ah_range - max_range)   )  
        index  = index[1]
        kappa_d[ind,j] = kappa_grid[index]
        # print(j)
  
#         print(kappa_grid[index])
#         temp = generate_poyla_tree_KSstatistic_d(kappa_alpha,theta0,N= 50,m,d,level_grid[j],total_rate,ind,kappa_grid[index] ,posterior_mean_ind)
#         hist(temp$KS ,20)
#         print( length(which(temp$KS < 1.36   ))/50)
#         
        kappa_alpha[2^{ind} : (2^{ind+1}-1)] = kappa_d[ind,j]*alpha[2^{ind} : (2^{ind+1}-1)]
        
      }###  close for ind 
    }### close  for level_grid
    kappa_d_mean =  kappa_d_mean + kappa_d/Niter
  }###  close for iter
  kappa = kappa_d_mean
  return( list(kappa_d = kappa_d_mean) )
}


########################################################################################3
########################################################################################3
########################################################################################3
########################################################################################3



poylaTreeBF_d2  = function(x,t,l ,alpha,kappa_alpha)
{
  s = max(t-l,1)
  
  log_bayes_factor = rep(0,t-s)
  cumulative_counts = rep(0,length(x))

  ind = 1
   for(tt in  t:s )
   {    
     cumulative_counts = cumulative_counts + x[tt,]
     ycounts =  gridCounts(cumulative_counts ,d) 
     log_bayes_factor[ind] = logPoylaTreeProb(ycounts, alpha,d) - logPoylaTreeProb(ycounts, kappa_alpha,d)
     ind = ind +1
   }
 
#     aa = -log_bayes_factor[1:(ind-1)]
#     
#     max_aa = max(aa)
#     cc =  sum(exp(aa - max_aa))
#      
#     cc = log(cc ) + max_aa
#    
#   print(-log_bayes_factor[1:(ind-1)])
    
#     aa = -log_bayes_factor[1:(ind-1)]
#     bb = aa
#    
#     plot(aa)
    
#   
#   aa =  rev(cumsum(rev(aa)))
#   
#    max_aa = max(aa)
#    cc =  sum(exp(aa - max_aa))
#    
#   cc = log(cc ) + max_aa

  return( min(log_bayes_factor[1:(ind-1)]) )
 # return(cc)
}

####################################################################
####################################################################
###################################################

# choose_tres_hold(N = 1000,theta0_hat,total_rate,d,alpha,kappa_alpha[,i], size_prob,L)
choose_tres_hold = function(N,theta0_hat,total_rate,d,alpha,kappa_alpha, size_prob,L)
{
  x = matrix(0,N,length(theta0_hat))
  
  log_bayes_factor =   rep(0,N)
  for(iter in 1:N)
  {
     x[iter,] = rpois( length(theta0_hat), total_rate*theta0_hat   )

     log_bayes_factor[iter] = poylaTreeBF_d2(x,iter,l = L,alpha,kappa_alpha)
       #logPoylaTreeProb(ycounts, alpha,d) - logPoylaTreeProb(ycounts, kappa_alpha,d)
  }
  
   #plot(-log_bayes_factor)
   return(quantile(log_bayes_factor,prob = 1-size_prob))
}

####################################################################
####################################################################
###################################################

choose_tres_hold_poyla_H = function(N,theta0_hat,total_rate,d,alpha,kappa, size_prob = .95,L)
{
  x = matrix(0,N,length(theta0_hat))
  
  log_bayes_factor =   rep(0,N)
  for(iter in 1:N)
  {
    x[iter,] = rpois( length(theta0_hat), total_rate*theta0_hat   )
    
    log_bayes_factor[iter] = poylaTreeBF(x,iter,l = L,alpha,kappa)
  }
  
 # plot(-log_bayes_factor)
  return(quantile(log_bayes_factor,prob = 1 - size_prob))
}
####################################################################
####################################################################
###################################################
choose_tres_hold_dir  = function(N,theta0_hat,total_rate,d,Dir_par,kappa_dir, size_prob = .95,L)
{
  x = matrix(0,N,length(theta0_hat))
  m = 2^d
  
  log_bayes_factor =   rep(0,N)
  for(iter in 1:N)
  {
    x[iter,] = rpois( length(theta0_hat), total_rate*theta0_hat   )
    
    log_bayes_factor[iter] = BF_dirichlet_prior(kappa_dir,Dir_par,x,m,l=L,iter)
  }
  
  #plot(log_bayes_factor)
  return(quantile(log_bayes_factor,prob = 1 - size_prob))
}

####################################################################
####################################################################
###################################################
choose_tres_hold_KS = function(N = 200,theta0_hat,total_rate,d,F0, size_prob = .95,L)
{
  x = matrix(0,N,length(theta0_hat))
  m = 2^d
  
  KS =   rep(0,N)
  for(iter in 1:N)
  {
    x[iter,] = rpois( length(theta0_hat), total_rate*theta0_hat)
    
    KS[iter] = sequential_KS_statistics(x,iter,F0,L)
  }
  
  #plot(KS)
  return(quantile(KS,prob = size_prob))
}

####################################################################
####################################################################
##################################################
choose_tres_hold_log_Rnf = function(N = 200,theta0_hat,total_rate,d, size_prob = .95,L)
{
  x = matrix(0,N,length(theta0_hat))
  m = 2^d
  
  log_Rnf_statistic  =   rep(0,N)
  for(iter in 1:N)
  {
    x[iter,] = rpois( length(theta0_hat), total_rate*theta0_hat   )
    
    log_Rnf_statistic[iter] =  logRnf(iter,1,1,x,L,total_rate*theta0_hat )
  }
  
 # plot(  log_Rnf_statistic)
  return(quantile(log_Rnf_statistic,prob = size_prob))
}
########################################################################33
###################################################################################
####################################################################################
#####################################################################################33


mle_detection = function(x,t,L,lambda0)
{
  m = length(lambda0)
  
  init = max(t-L,1)
#  x = x[init:t,]
  
  if(  t> init)
  {
    total_sums  = colSums(x[init:t,])   
  }

  if(t == init)
  {
    total_sums = x[t,]
  }

  
  test = rep(0,t-init +1)
  
  for(s in init:t)
  {
    lambda_hat = total_sums/(t-s +1)
    
    test[ s - init +1] = (t-s +1)*sum(lambda0 - lambda_hat)
    test[ s - init +1] =  test[ s - init +1] + sum(x[init:t,] %*%  (-log(lambda0+ 10^-99) + log(lambda_hat+ 10^-99)) )    
    
    if(s< t + 1)
    {
      total_sums = total_sums  - x[s,] 
    }
  }
  return(max(test))
} 

####################################################################
####################################################################
###################################################
choose_tres_hold_mle = function(N = 200,theta0_hat,lambda0_hat,d, size_prob = .95,L)
{
  x = matrix(0,N,length(theta0_hat))
  m = 2^d
  
  mle =   rep(0,N)
  for(iter in 1:N)
  {
    x[iter,] = rpois( length(theta0_hat), lambda0_hat)
    
    mle[iter] = mle_detection(x,iter,L,lambda0_hat)
  }
 # plot(mle)
  #plot(KS)
  return(quantile(mle,prob = size_prob))
}



####################################################################
####################################################################
###################################################
# 
# finding_kappa_alpha_withData = function(level_grid,theta0,alpha,d,m,x,F0)
# {
#   kappa_d =  matrix(0,d,length(level_grid))  #a
#   kappa_d_mean = kappa_d
#   posterior_mean = cumsum(theta0)
#   kappa_alpha = alpha
#   
#   theta0_dist = gridCounts(theta0,d)
#   kappa_alpha = alpha
#   
#   Niter = 4
#   for(iter in 1:Niter)
#   {
#     for(j in 1:length(level_grid))
#     {
#       kappa_alpha = alpha
#       for(ind in 1:d)
#       {
#         posterior_mean_ind = cumsum(theta0_dist[2^{ind}:(2^{ind+1}-1)])
#         kappa_grid = 10^seq(-7,-1,length= 20)
#         Ah_range = rep(0,length(kappa_grid))
#         max_range =  1.36
#         
#         for(i in 1:length(kappa_grid))
#         {
#            KS = rep(0, floor(dim(x)[1]/Niter))
#            
#            for(jj in 1:floor(dim(x)[1]/Niter))
#            {
#              KS[jj] = kolmogorov_statistic(x[ jj + (iter-1)* floor(dim(x)[1]/Niter), ],F0)
#            }
#            
#            Ah_range[i] =  quantile(KS,prob= 1-level_grid[j])
#         }###  iter = 1
#         
#         
#         index= which.min(    abs(Ah_range - max_range)   )  
#         index  = index[1]
#         kappa_d[ind,j] = kappa_grid[index]
#         # print(j)
#         
#         kappa_alpha[2^{ind} : (2^{ind+1}-1)] = kappa_d[ind,j]*alpha[2^{ind} : (2^{ind+1}-1)]
#         
#       }###  close for ind 
#     }### close  for level_grid
#     kappa_d_mean =  kappa_d_mean + kappa_d/Niter
#   }###  close for iter
#   kappa = kappa_d_mean
#   return( list(kappa_d = kappa_d_mean) )
# }
# 
choose_tres_hold_withData = function(theta0_hat,x_tilde,d,alpha,kappa_alpha, size_prob,L)
{
  
  N = dim(x_tilde)[1]
   
  log_bayes_factor =   rep(0,N)
   
  for(iter in 1:N)
  {
    log_bayes_factor[iter] = poylaTreeBF_d2(x_tilde,iter,l = L,alpha,kappa_alpha)
  }
  
  #plot(-log_bayes_factor)
  return(quantile(log_bayes_factor,prob = 1-size_prob))
}

####################################################################
####################################################################
###################################################
choose_tres_hold_dir_withData  = function(theta0_hat,x_tilde,d,Dir_par,kappa_dir, size_prob = .95,L)
{
  #x = matrix(0,N,length(theta0_hat))
  m = 2^d
  N = dim(x_tilde)[1]
  
  log_bayes_factor =   rep(0,N)
  for(iter in 1:N)
  {
    #x[iter,] = rpois( length(theta0_hat), total_rate*theta0_hat   )
    log_bayes_factor[iter] = BF_dirichlet_prior(kappa_dir,Dir_par,x_tilde,m,l=L,iter)
  }
  
  #plot(log_bayes_factor)
  return(quantile(log_bayes_factor,prob = 1 - size_prob))
}


####################################################################
####################################################################
###################################################
choose_tres_hold_KS_withData = function(theta0_hat,x_tilde,d,F0, size_prob = .95,L)
{
  #x = matrix(0,N,length(theta0_hat))

  N = dim(x_tilde)[1]
  
  KS =   rep(0,N)
  for(iter in 1:N)
  {
  #  x[iter,] = rpois( length(theta0_hat), total_rate*theta0_hat)
    KS[iter] = sequential_KS_statistics(x_tilde,iter,F0,L)
  }
  
  quantile(KS,prob = size_prob)
  
  #plot(KS)
  return(quantile(KS,prob = size_prob))
}

#choose_tres_hold_KS_withData(theta0_hat,x_tilde,d,F0, size_prob = .99,L)


####################################################################
####################################################################
##################################################
choose_tres_hold_log_Rnf_withData = function(theta0_hat,x_tilde,total_rate,d, size_prob = .95,L)
{
  
  m = 2^d
  N = dim(x_tilde)[1]
  
  log_Rnf_statistic  =   rep(0,N)
  for(iter in 1:N)
  {
    log_Rnf_statistic[iter] =  logRnf(iter,1,1,x_tilde,L,total_rate*theta0_hat )
  }
  
  # plot(  log_Rnf_statistic)
  return(quantile(log_Rnf_statistic,prob = size_prob))
}

####################################################################
####################################################################
###################################################
choose_tres_hold_mle_withData = function(x_tilde,lambda0,d, size_prob = .95,L)
{
  m = 2^d
  N = dim(x_tilde)[1]
  
  mle =   rep(0,N)
  for(iter in 1:N)
  {
    mle[iter] = mle_detection(x_tilde,iter,L,lambda0)
  }
  
  #plot(KS)
  return(quantile(mle,prob = size_prob))
}
####################################################################
####################################################################
###################################################

log_p_value  = function(x,t,l ,prop)
{
  s = max(t-l,1)
  
  log_bayes_factor = rep(0,t-s)
  cumulative_counts = rep(0,length(x))
  
  ind = 1
  for(tt in  t:s )
  {    
    cumulative_counts = cumulative_counts + x[tt,]
    ycounts =  gridCounts(cumulative_counts ,d) 
    
    ####################################
    
    #aa = rep(0,2^d-1)
    # for(j in 1:(2^d-1))
    #  {
    indices  =  1:(2^d-1) 
    aa = dbinom(ycounts[2*indices],ycounts[indices], prop[2*indices]  ,log = TRUE )
    # }
    
    
    log_bayes_factor[ind] = min(aa) 
    ind = ind +1
  }
  
  return( min(log_bayes_factor[1:(ind-1)]) )
  # return(cc)
}



########################################################################################3
########################################################################################3
########################################################################################3
########################################################################################3


choose_tres_p_value= function(N = 200,prop,theta0_hat,lambda0_hat,d, size_prob = .95,L)
{
  x = matrix(0,N,length(theta0_hat))
  m = 2^d
  
  pvalue =   rep(0,N)
  for(iter in 1:N)
  {
    x[iter,] = rpois( length(theta0_hat), lambda0_hat)
    
    pvalue[iter] = log_p_value(x,iter,L,prop)
  }
  # plot(mle)
  #plot(KS)
  return(quantile(pvalue,prob = 1-size_prob))
}



########################################################################################3
################################################################

####################################################################
####################################################################
###################################################
choose_tres_hold_mle_withData = function(x_tilde,lambda0,d, size_prob = .95,L)
{
  m = 2^d
  N = dim(x_tilde)[1]
  
  mle =   rep(0,N)
  for(iter in 1:N)
  {
    mle[iter] = mle_detection(x_tilde,iter,L,lambda0)
  }
  
  #plot(KS)
  return(quantile(mle,prob = size_prob))
}

####################################################################
####################################################################
#######################################################################################################################
####################################################################
###################################################
log_p_value2  = function(x,t,l ,prop)
{
  s = max(t-l,1)
  
  log_bayes_factor = rep(0,t-s)
  cumulative_counts = rep(0,length(x))
  
  ind = 1
  for(tt in  t:s )
  {    
    cumulative_counts = cumulative_counts + x[tt,]
    ycounts =  gridCounts(cumulative_counts ,d) 
    
    ####################################
    
    #aa = rep(0,2^d-1)
    # for(j in 1:(2^d-1))
    #  {
    indices  =  1:(2^d-1) 
    aa = dbinom(ycounts[2*indices],ycounts[indices], prop[2*indices]  ,log = TRUE )
    # }
    
    
    log_bayes_factor[ind] = sum(aa) 
    ind = ind +1
  }
  
  return( min(log_bayes_factor[1:(ind-1)]) )
  # return(cc)
}



########################################################################################3
########################################################################################3
########################################################################################3
########################################################################################3


choose_tres_p_value2= function(N = 200,prop,theta0_hat,lambda0_hat,d, size_prob = .95,L)
{
  x = matrix(0,N,length(theta0_hat))
  m = 2^d
  
  pvalue =   rep(0,N)
  for(iter in 1:N)
  {
    x[iter,] = rpois( length(theta0_hat), lambda0_hat)
    
    pvalue[iter] = log_p_value2(x,iter,L,prop)
  }
  # plot(mle)
  #plot(KS)
  return(quantile(pvalue,prob = 1-size_prob))
}



########################################################################################3
################################################################




# choose_tres_hold_KS_old = function(begining = 100,N = 200,theta0_hat,total_rate,d,F0, size_prob = .95,L)
# {
#   x = matrix(0,N,length(theta0_hat))
#   m = 2^d
#   
#   KS =   rep(0,N)
#   cum_x = rep(0,m)
#   
#   for(iter in 1:N)
#   {
#     x[iter,] = rpois( length(theta0_hat), total_rate*theta0_hat)
#     
#     cum_x = cum_x + x[iter,]
#     
#     KS[iter] = kolmogorov_statistic(cum_x,F0)
#   }
#   
#  #plot(KS)
#   return(quantile(KS[begining:N],prob = size_prob))
# }

choose_tres_hold_KS_old = function(begining = 100,N = 200,theta0_hat,total_rate,d,F0, size_prob = .95,L)
{
  x = matrix(0,N,length(theta0_hat))
  m = 2^d
  
  # print(m)
  KS =   rep(0,N)
  cum_x = rep(0,m)
  
  for(iter in 1:N)
  {
    x[iter,] = rpois( length(theta0_hat), total_rate*theta0_hat)
    
    cum_x = cum_x + x[iter,]
    
    KS[iter] = sqrt(sum(cum_x))*kolmogorov_statistic(cum_x,F0)
  }
  
  #plot(KS)
  return(quantile(KS[begining:N],prob = size_prob))
}
####################################################################
####################################################################
###################################################
new_choose_tres_hold_KS = function(N = 200,theta0_hat,total_rate,d,F0, size_prob = .95,L)
{
  x = matrix(0,N,length(theta0_hat))
  m = 2^d
  
  KS =   rep(0,N)
  for(iter in 1:N)
  {
    x[iter,] = rpois( length(theta0_hat), total_rate*theta0_hat)
    
    KS[iter] = sequential_KS_statistics(x,iter,F0,L)
  }
  
  #plot(KS)
  return(KS)
}

###########################################################################
##################################################################################
##################################################################################

new_choose_tres_hold_KS_old = function(begining = 100,N = 200,theta0_hat,total_rate,d,F0, size_prob = .95,L)
{
  x = matrix(0,N,length(theta0_hat))
  m = 2^d
  
  # print(m)
  KS =   rep(0,N)
  cum_x = rep(0,m)
  
  for(iter in 1:N)
  {
    x[iter,] = rpois( length(theta0_hat), total_rate*theta0_hat)
    
    cum_x = cum_x + x[iter,]
    
    KS[iter] = sqrt(sum(cum_x))*kolmogorov_statistic(cum_x,F0)
  }
  
  #plot(KS)
  return(KS)
}

####################################################################
####################################################################
##################################################
new_choose_tres_hold_log_Rnf = function(N = 200,theta0_hat,total_rate,d, size_prob = .95,L)
{
  x = matrix(0,N,length(theta0_hat))
  m = 2^d
  
  log_Rnf_statistic  =   rep(0,N)
  for(iter in 1:N)
  {
    x[iter,] = rpois( length(theta0_hat), total_rate*theta0_hat   )
    
    log_Rnf_statistic[iter] =  logRnf(iter,1,1,x,L,total_rate*theta0_hat )
  }
  
  # plot(  log_Rnf_statistic)
  return(log_Rnf_statistic)
}


####################################################################
####################################################################
###################################################
new_choose_tres_hold_mle = function(N = 200,theta0_hat,lambda0_hat,d, size_prob = .95,L)
{
  x = matrix(0,N,length(theta0_hat))
  m = 2^d
  
  mle =   rep(0,N)
  for(iter in 1:N)
  {
    x[iter,] = rpois( length(theta0_hat), lambda0_hat)
    
    mle[iter] = mle_detection(x,iter,L,lambda0_hat)
  }
  # plot(mle)
  #plot(KS)
  return(mle)
}
####################################################################
###################################################


scr_update = function(lambda_t,omega_t,lambda0_new,xx,lamb)
{
  Tt = diag(- lambda_t[1]/(lambda_t[2:length(lambda_t)] )  )
  Tt = cbind(rep(1,length(lambda0_new)-1),Tt)
  
  alpha_t = Tt%*% xx
  
  S =  Tt %*% omega_t %*% t(Tt)
  aux = solve( S,  alpha_t)
  D2 =   sum(alpha_t*drop(aux))
  
  
  lambda_t =  (1-lamb)*lambda_t  + lamb*xx
  
  omega_t =  diag(lambda_t )
  
  return(list(D2 = D2,lambda_t = lambda_t,omega_t = omega_t))
}

########################################################################

donwsample = function(x)
{
   x2 = matrix(0,dim(x),285)
   
   for(t  in 1:dim(x)[1])
   {
     for( i in 1:285)
     {
       for(j in 1:7)
       {
         x2[t,i] = x2[t,i] + x[t,i*7 - j + 1] 
       }
     }
     x2[t,285] = x2[t,285] + sum(x[t,1996:1999])  
     #x2[t,401]  = x[t,1997] +x[t,1998] +x[t,1999]  
   }
   
  
   return(x2)  
}


########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################

L_infinity_cdf =  function(x,F0) 
{
  ecd = empirical_cdf(x)
  
  return( max(abs(F0 -ecd)) )  
}



sequential_KS_statistics_overdispersion = function(x,t,F0,L,phi)
{
  delta = rep(0,t)
  
  if(t ==1)
  {
    x_prime = x[1,]
    delta[1] =    L_infinity_cdf(x_prime,F0)/sqrt(phi) 
  }
  if(t > 1)
  {
    x_prime = colSums(x[max(1,t-L):t,])
    delta[1] = L_infinity_cdf(x_prime,F0)* sqrt((t - max(1,t-L) + 1)/phi)     
    
    
    for(s in max(1,t-L):(t-1))
    {
      x_prime = x_prime - x[s,]
      delta[s+1] = L_infinity_cdf(x_prime,F0) * sqrt((t - (s+1) + 1)/phi)   
    }
    
  }## case t>1
  
  return(max(delta))
}



####################################################################
####################################################################
###################################################

kolmogorov_distribution_cdf =  function(x_array)
{
  
  aux = rep(0,length(x_array))
  for(ii in 1:length(x_array))
  {
    k = 1:200
    x = x_array[ii]
    temp = -(1/(8*x*x))*pi*pi*(2*k -1)^2 
    aux[ii] = sqrt(2*pi)/x *sum(exp(temp)) 
    
  }  
  return(   aux  ) 
}


####################################################################
####################################################################
###################################################
####################################################################
####################################################################
###################################################


  choose_treshold_np_cusum = function(begining = 100,N = 200,lambda0,d, size_prob = .95,L)
  {
      
    x = matrix(0,N,length(lambda0))
    m = 2^d
    
    k = 0
    # print(m)
    
    d_prob = matrix(0,N,length(lambda0)+1)  
    for(iter in 1:N)
    {
      x[iter,] = rpois( length(theta0_hat), lambda0)
      
      xtilde =     x[iter,] -  lambda0
      ind  =  which.max(c(xtilde,0))
      d_prob[iter,ind ] =  d_prob[iter,ind] + 1
      
    }
    d_prob  =   colMeans(d_prob) + 10^-4
  #  d_prob = d_prob/sum(d_prob)
   # non_zero_ind  =   which(d_prob>0)  
    
    S1 =  matrix(0,N,m+1) #length(non_zero_ind))
    S2 =  matrix(0,N,m+1) #length(non_zero_ind))
    y =  rep(0,N)

    for(iter in 1:N)
    {
      x[iter,] = rpois( length(theta0_hat), lambda0)
      
      xtilde =     x[iter,] -  lambda0
      ind  =  which.max(c(xtilde,0))#which.min(xtilde)
      eta_1n = rep(0,length(lambda0)+1)
      eta_1n[ind] = 1
       
       
      if(iter ==1)
      {
         c_n = sum((eta_1n-d_prob)^2  / d_prob)
           #drop((eta_1n-d_prob)%*%diag(1/d_prob)%*% (eta_1n-d_prob)) 
         
         if(c_n > k)
         {
           S1[iter,] = eta_1n*(c_n - k)/c_n 
           S2[iter,] = d_prob*(c_n - k)/c_n 
         }
             
        y[iter] = sum((S1[iter,]-S2[iter,])^2   / S2[iter,]  )
          #drop((S1[iter,]-S2[iter,])%*%diag(1/S2[iter,])%*%(S1[iter,]-S2[iter,])) 
      }##########
      
      if(iter > 1)
      {
          temp1 = S1[iter-1,] - S2[iter-1,] +  eta_1n  - d_prob
          temp2 =   S2[iter-1,] +   d_prob
          
          c_n = sum((temp1)^2 /temp2)
            #drop(temp1%*%diag(temp2)%*%temp1 ) 
           
          if(c_n >k)
          {
            S1[iter,] = (S1[iter-1,] + eta_1n)*(c_n - k)/c_n 
            S2[iter,] = (S2[iter-1,] + d_prob)*(c_n - k)/c_n
          }
          y[iter] = sum((S1[iter,]-S2[iter,])^2   / S2[iter,]  ) 
            #drop((S1[iter,]-S2[iter,])%*%diag(1/S2[iter,])%*%(S1[iter,]-S2[iter,])) 
#           if(is.nan(y[iter]))
#           {break;}
      }
    }### close for iter
    
    #plot(KS)
    return(list(y=y, d_prob = d_prob ))
  }
  
cusum_statistic = function(x_vec,d_prob,S1_prev,S2_prev,k=0)
{ 
 # non_zero_ind  =   which(d_prob>0)  
  
  ind  =  which.max(x_vec)
  eta_1n = rep(0,length(x_vec))
  eta_1n[ind] = 1
  
  temp1 = S1_prev - S2_prev +  eta_1n  - d_prob
  temp2 =   S2_prev +   d_prob
  
  c_n = sum((temp1)^2 /temp2)
  #drop(temp1%*%diag(temp2)%*%temp1 ) 
  
  S1_now =  rep(0,length(x_vec))
  S2_now =  rep(0,length(x_vec))
  
  if(c_n >k)
  {
    S1_now = (S1_prev + eta_1n)*(c_n - k)/c_n 
    S2_now = (S2_prev + d_prob)*(c_n - k)/c_n
  }
  y = sum((S1_now -S2_now)^2   /S2_now  )
  
  return(list(y =y,S1_now = S1_now,S2_now = S2_now))
  #drop((S1[iter,]-S2[iter
}