

#rm(list = ls())

#setwd("/home/hernanmp/Desktop/change_point")



setwd(paste(working_directory,"/Code",sep=""))

theta0 = background_train
thetac = sim_spectrum
d= 11

theta0_new =  rep(0,2^{d-3})
thetac_new =  rep(0,2^{d-3})

for( i in 1:(2^{d-3}))
{
  for(j in 1:8)
  {
    theta0_new[i] = theta0_new[i] + theta0[i*8 - j + 1] 
    thetac_new[i] =  thetac_new[i]+ thetac[i*8 - j + 1]
  }
}

theta0 = theta0_new
thetac =  thetac_new

matplot(cbind(theta0,thetac),type ="l")

d = 8
lamb = 0.002


matplot(cbind(theta0,thetac),type ="l")

source("Generate_data.R")


matplot(cbind(theta0,thetac),type="l")

####  setting parameters for simulations
T_grid   =   c(300)
d_grid = c(d)
Total_rate_factor_grid = c(1000/2^d,500/2^d,100/2^d) #c(2,3,4)
NMC = 100


L = 50   ## window length for KS
k_alpha = 3 ## KS


k_alpha_grid = c(1.36) ## this is changed later
scr_tresdhold  = 0


size_prob  = .999### 1 - size_prob  is the probability of a false alarm







#v = min(floor(runif(1)*T),T)
######################################################### ######################################   
ind_T  = 1
ind_d  = 1
ind_TRF = 1
T = T_grid[ind_T]
d = d_grid[ind_d]


average_false_alarms_scr = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid)))
average_delay_time_scr = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid)))




for(ind_T in 1:length(T_grid))
{
  T = T_grid[ind_T]
  ######################################  
  for(ind_d in 1:length(d_grid))
  {
    d = d_grid[ind_d]
    ######################################  
    for(ind_TRF in 1:length(Total_rate_factor_grid))
    {
      Total_rate_factor = Total_rate_factor_grid[ind_TRF ]
      
      ###  True normalized pre and post densities 
      ## pre change
      m = 2^d
      
      
      
      lambda0 = Total_rate_factor*m*theta0
      lambdac = Total_rate_factor*m*thetac
      
      
      ###  Training period
      temp = Training_period(ntraining = 1000,m,d,lambda0) 
      F0 = temp$F0
      alpha = temp$alpha
      p =  temp$p  
      Dir_par =  temp$dir_par
      total_rate = temp$total_rate
      matplot(cbind(lambda0,lambdac),type="l")
      theta0_hat = temp$theta0
      
      ################################################################      
      
      scr_tresdhold  = 0
      for(ii in 1:20)
      {
        D2 = rep(0,1000)
        lambda_t = lambda0
        omega_t =  diag(lambda0)  
        for(t in 1:1000)
        {
          xx = rpois(length(theta0), lambda0 )
          
          Tt = diag(- lambda_t[1]/(lambda_t[2:length(lambda_t)]  + 10^-99)  )
          Tt = cbind(rep(1,length(lambda0)-1),Tt)
          
          alpha_t = Tt%*% xx
          
          S =  Tt %*% omega_t %*% t(Tt)
          aux = solve( S,  alpha_t)
          D2[t] =   sum(alpha_t*drop(aux))
          
          
          lambda_t =  (1-lamb)*lambda_t  + lamb*xx
          
          omega_t =  diag(lambda_t)
          
        }
        aa = quantile(D2,prob =  size_prob)
        scr_tresdhold =  scr_tresdhold + aa/20
      }
      
      ################################################################        
      
      
      ######################################    
      
      false_alarms_scr = rep(0,NMC)
      delay_time_scr =   rep(0,NMC) 
      
      
      v_array = rep(0,NMC )
      
      for(iter in 1:NMC)
      { 
        print("iter")
        print(iter) 
        
        scr = rep(0,T) 
        
        
        x =  matrix(0,T,m)
        
        ## v = 100
        v = min(floor( runif(1)*((T-1-100)-100) + 100    ),T)
        v_array[iter] = v
        
        
        best_scr =  -10^15
        # best_KS =  -10^15
        
        lambda_t = lambda0
        omega_t =  diag(lambda0)
        
        ptm <- proc.time()
        for(t in 1:T)
        {
          lambda0_t = lambda0 #+ (.1)*rnorm(length(lambda0)) 
          #indices = which(lambda0_t<0)
          #lambda0_t[indices] = lambda0[indices]
          
          lambdac_t = lambdac #+ (.1)*rnorm(length(lambdac)) 
          #indices = which(lambdac_t<0)
          #lambdac_t[indices] = lambdac[indices]
          
          if(t <= v) 
          { 
            x[t, ] =rpois(length(theta0), lambda0_t )
          }  
          if(t > v)
          { 
            x[t,] = rpois(length(theta0),lambdac_t)
          } #rpois(length(theta0),lambdac_t)}
          #ycounts =  gridCounts(x[t,],d)
          
          ## csr
          if(t <= v)
          {
            xx = x[t,]
            
            Tt = diag(- lambda_t[1]/(lambda_t[2:length(lambda_t)]  + 10^-99)  )
            Tt = cbind(rep(1,length(lambda0)-1),Tt)
            
            alpha_t = Tt%*% xx
            S =  Tt %*% omega_t %*% t(Tt)
            aux = solve( S,  alpha_t)
            
            lambda_t =  (1-lamb)*lambda_t  + lamb*xx
            
            omega_t =  diag(lambda_t)
            scr[t] =   sum(alpha_t*drop(aux))
          }
          
          if(t > v)
          {
            if(best_scr <  scr_tresdhold  )
            {
              xx = x[t,]
              
              Tt = diag(- lambda_t[1]/(lambda_t[2:length(lambda_t)]  + 10^-99)  )
              Tt = cbind(rep(1,length(lambda0)-1),Tt)
              
              alpha_t = Tt%*% xx
              
              lambda_t =  (1-lamb)*lambda_t  + lamb*xx
              
              omega_t =   diag(lambda_t)  
              #lamb*outer(xx - lambda_t, xx - lambda_t, FUN = function(x,y){x*y} )  +  (1-lamb)*omega_t
              
              S =  Tt %*% omega_t %*% t(Tt)
              aux = solve( S,  alpha_t)
              scr[t] =   sum(alpha_t*drop(aux))
            }
            if( scr[t] > best_scr && scr[t] != 0)
            {
              best_scr = scr[t] 
            }
            
          }
          
          
          
        }## close for 1:T      
        proc.time() - ptm
        
        
        ## update delayed time            
        
        
        ind2 = which(scr[(1+v):T] >   scr_tresdhold )
        delay_time_scr[iter] = T - (v+1)
        
        if( length(ind2) > 0)
        {
          delay_time_scr[iter] = min(ind2) -1
        }
        
        
        
        
      }## close for NMC simulations, iter
      
      
      average_delay_time_scr[ind_T,ind_d,ind_TRF,] = mean(  delay_time_scr)
      print(average_delay_time_scr[ind_T,ind_d,ind_TRF,])
      
    }## close for Total rate
  }## close for d
}## close for T