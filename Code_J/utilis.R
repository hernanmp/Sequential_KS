


choose_tres_hold_KS   = function(N = 1000,mu,sigma,total_rate, size_prob,L)
{
  
  x = matrix(0,N,total_rate)
  m = 2^d
  
  KS =   rep(0,N)
  for(iter in 1:N)
  {
    
    x[iter,] = rnorm(total_rate,mu,sigma)
    #sqrt(total_rate)*(pnorm(x[iter,]) -   (1:total_rate)/total_rate)
    if(iter ==1)
    {
      x[iter,] =  sort(x[iter,])
      KS[iter] = sqrt(total_rate)*max(pnorm(x[iter,],mu,sigma) -   (1:total_rate)/total_rate)
    }
    if(iter >1)
    {
      xx = sort(as.vector(x[max(1,iter-L):iter,])) 
      KS[iter] = sqrt(length(xx))*max(pnorm(xx,mu,sigma)  - (1:length(xx))/length(xx))     
    }

  }
  
  #plot(KS)
  return(KS)
  #(quantile(KS,prob = size_prob))
  
  
}

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
choose_tres_hold_log_Rnf = function(N = 1000,mu,sigma,total_rate,size_prob,L,tau)
{
  
  
  x = rep(0,N)
  m = 2^d
  
  sigma_new = sigma/sqrt(total_rate) 
  
  Rnf =   rep(0,N)
  for(iter in 1:N)
  {
    
    x[iter] = sum(rnorm(total_rate,mu,sigma))/total_rate
    
    if(iter ==1)
    {
       Rnf[iter] = exp(.5*((x[1]^2)/sigma_new^4)/(1/sigma_new^2 + 1/tau^2))*(1/sigma_new^2 + 1/tau^2  )^{-1/2}
    }
    if(iter >1)
    {
      init = max(1,iter-L)
      xx = x[init:iter]
      
      aa = rev(cumsum(rev(xx)))
      bb = exp( .5*( (aa^2)/sigma_new^4)/(rev(1:length(xx))/sigma_new^2 + 1/tau^2)  )
      Rnf[iter] = sum( bb*( rev(1:length(xx))/sigma_new^2 + 1/tau^2 )^{-1/2})
        #sum(  exp( ( aa^2/sigma_new^4)/(length(xx)/sigma_new^2 + 1/tau^2)*( rev(1:length(xx))/sigma_new^2 + 1/tau^2 )^{-1/2})  )
    }
    
  }
  
  #plot(Rnf)
  #return(quantile(Rnf,prob = size_prob))
  return(Rnf)
}




################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
choose_tres_hold_log_mle = function(N = 1000,mu,sigma,total_rate,size_prob,L)
{
  
  
  x = rep(0,N)
  m = 2^d
  
  sigma_new = sigma/sqrt(total_rate) 
  
  mle =   rep(0,N)
  for(iter in 1:N)
  {
    
    x[iter] = sum(rnorm(total_rate,mu,sigma))/total_rate
    
    if(iter ==1)
    {
      yhat = x[1]
      mle[iter] = -1/2*yhat^2  + x[1]*yhat  
    }
    if(iter >1)
    {
      init = max(1,iter-L)
      xx = x[init:iter]
      
      aa = rev(1:length(xx))
      bb =  rev(cumsum(rev(xx)))
      yhat =bb/aa
      mle[iter] = max(-.5*aa*(yhat)^2 + bb*yhat) 
    }
    
  }
  
#  plot(mle)
  #return(quantile(mle,prob = size_prob))

  return(mle)  
}


##############################33##############################33##############################33
##############################33##############################33##############################33
##############################33##############################33##############################33
##############################33##############################33
##############################33##############################33
##############################33##############################33##############################33


sequential_KS_statistics = function(x,t,mu,sigma,L)
{
  xx = sort(as.vector(x[max(1,t-L):t,])) 
  KS = sqrt(length(xx))*max(pnorm(xx,mu,sigma)  - (1:length(xx))/length(xx))     
  return(KS)
}

##############################33##############################33##############################33
##############################33##############################33##############################33
##############################33##############################33##############################33
##############################33##############################33
##############################33##############################33
##############################33##############################33##############################33

logRnf  = function(t,tau,x,L,mu,sigma,total_rate)
{
  sigma_new = sigma/sqrt(total_rate)
  
  if(t ==1)
  {
    aa = mean(x[t,])
    Rnf = exp(.5*((aa^2)/sigma_new^4)/(1/sigma_new^2 + 1/tau^2))*(1/sigma_new^2 + 1/tau^2  )^{-1/2}
  }
  if(t >1)
  {
    init = max(1,t-L)
    xx = rowMeans(x[init:t,])
    
    aa = rev(cumsum(rev(xx)))
    
    bb = exp( .5*( (aa^2)/sigma_new^4)/(rev(1:length(xx))/sigma_new^2 + 1/tau^2)  )
    Rnf = sum( bb*( rev(1:length(xx))/sigma_new^2 + 1/tau^2 )^{-1/2})
    
    #bb =  exp( .5*( (aa^2)/sigma_new^4)/(length(xx)/sigma_new^2 + 1/tau^2))
  #  Rnf = sum( *( rev(1:length(xx))/sigma_new^2 + 1/tau^2 )^{-1/2})  )
  }
  return(Rnf)
}

##############################33##############################33##############################33
##############################33##############################33##############################33
##############################33##############################33##############################33
##############################33##############################33##############################33
##############################33##############################33##############################33
##############################33##############################33##############################33
##############################33##############################33##############################33




logmle  = function(t,x,L,mu,sigma,total_rate)
{
  sigma_new = sigma/sqrt(total_rate)
  
  
  if(t ==1)
  {
    yhat = mean(x[1,])
    mle = -1/2*yhat^2  + yhat*yhat  
  }
  if(t >1)
  {
    init = max(1,t-L)
    xx = x[init:t,]
  #  print(dim(xx))
    xx = rowMeans(xx)
    
    aa = rev(1:length(xx))
    bb =  rev(cumsum(rev(xx)))
    yhat =bb/aa
    mle = max(-.5*aa*(yhat)^2 + bb*yhat) 
  }
  
  return(mle)
}


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

choose_tres_hold_PKS   = function(N = 1000,mu,sigma,total_rate, size_prob,L)
{
  
  x = matrix(0,N,total_rate)
  m = 2^d
  
  KS =   rep(0,N)
  for(iter in 1:N)
  {
    
    x[iter,] = rnorm(total_rate,mu,sigma)
    #sqrt(total_rate)*(pnorm(x[iter,]) -   (1:total_rate)/total_rate)
    if(iter ==1)
    {
      x[iter,] =  sort(x[iter,])
      KS[iter] = total_rate*max(pnorm(x[iter,],mu,sigma) -   (1:total_rate)/total_rate)
    }
    if(iter >1)
    {
      xx = sort(as.vector(x[1:iter,])) 
      KS[iter] = length(xx)*max(pnorm(xx,mu,sigma)  - (1:length(xx))/length(xx))     
    }
    
  }
  
  #plot(KS)
  return(KS)
  #(quantile(KS,prob = size_prob))
  
  
}

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
choose_tres_hold_log_mle2 = function(N = 1000,mu,sigma,total_rate,size_prob,L)
{
  
  N = N+1
  x = rep(0,N)
  m = 2^d
  
  sigma_new = sigma/sqrt(total_rate) 
  
  mle =   rep(0,N)
  for(iter in 1:N)
  {
    
    x[iter] = sum(rnorm(total_rate,mu,sigma))/total_rate
    
#     if(iter ==1)
#     {
#       yhat = x[1]
#       mle[iter] = -1/2*yhat^2  + x[1]*yhat  
#     }
    if(iter >1)
    {
      init = max(1,iter-L)
      xx = x[init:iter]
      
      aa = rev(1:length(xx))
      bb =  rev(cumsum(rev(xx)))
      yhat = bb[1:(length(bb)-1)]/aa[1:(length(bb)-1)]
      sigma_hat =  rep(0,length(xx)-1)
      
      for(jj in 1:(length(xx)-1))
      {
        sigma_hat[jj] =  sum((yhat[jj] - xx[jj: length(xx)])^2)/aa[jj]
      }
      
      
      log_lik = rep(0,(length(xx)-1))
      
      for(jj in 1:(length(xx)-1))
      {
        log_lik[jj] =  sum(dnorm(xx[jj:length(xx)],yhat[jj],sigma_hat[jj],log = TRUE) - dnorm(xx[jj:length(xx)],mu,sigma,log = TRUE)  )
      }
      
      mle[iter] =  max(log_lik)
    }
    
  }
  
  #  plot(mle)
  #return(quantile(mle,prob = size_prob))
  
  return(mle[2:N])  
}




