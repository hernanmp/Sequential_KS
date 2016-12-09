
rm(list = ls())

mu_c = .3;
mu = 0

sigma= 6
sigma_c = 6

tau_grid = c(.1,1,5,10)


setwd("/Users/user/Desktop/JASA_code_Anomaly/Experiments/Code")
source("utilis.R")


####  setting parameters for simulations
T_grid   =   c(1100)
d_grid = c(11)
d=  11
Total_rate_factor_grid = c(100/2^d,500/2^d,1000/2^d) 
NMC = 100


L = 50   ## window length for KS
k_alpha = 3 ## KS

k_alpha_grid = c(1.36) ## this is changed later
k_alpha_grid_old = k_alpha_grid
k_alpha_grid_star  = k_alpha_grid
level_grid =  c(.9,.7,.5,.3) #c(.95,.99)
level_grid_dir = c(.9,.7,.5,.3) 
log_A_grid = rep(0,length(tau_grid)) ## this is changed later
log_A_grid_mle = c(0) 

size_prob  = .999### 1 - size_prob  is the probability of a false alarm


tres_hold_KS = rep(0,length(k_alpha_grid))
tres_hold_log_Rnf = rep(0,length(log_A_grid))


#v = min(floor(runif(1)*T),T)
######################################################### ######################################   
ind_T  = 1
ind_d  = 1
ind_TRF = 1
T = T_grid[ind_T]
d = d_grid[ind_d]


average_false_alarms_KS = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid)))
average_false_alarms_PKS = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid_old)))
average_false_alarms_Rnf = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(log_A_grid)))
average_false_alarms_mle = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(log_A_grid_mle)))
average_false_alarms_KS_star = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid_star)))



average_delay_time_KS = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid)))
average_delay_time_Rnf = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(log_A_grid)))
average_delay_time_mle = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(log_A_grid_mle)))
average_delay_time_PKS = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid_old)))
average_delay_time_KS_star = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid_star)))


NN  =   100
T2 = 1000  #

ks =  matrix(0,NN,T2)
rnf = array(0,c(length(tau_grid),NN,T2))
mle_tr = matrix(0,NN,T2)
pks = matrix(0,NN,T2)

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
      
      total_rate = m*Total_rate_factor
      ###  Training period
      
      
      ################################################################        
      ## KS^* threshold
      
      k_alpha_grid_star[1]   = sqrt(log(2*1000*L)/2)
      ################################################################      
      ### KS threshold
      #       
      #       
      k_alpha_grid[1]  = 0
      for(ii in 1:NN)
      {
        #bb = choose_tres_hold_KS(N = 1000,mu,sigma,total_rate, size_prob,L)
        
        ks[ii, ]  = choose_tres_hold_KS(N = 1000,mu,sigma,total_rate, size_prob,L)
        
        # k_alpha_grid[1] = k_alpha_grid[1] + bb/10
      }
      
      lower = min(ks)
      upper = max(ks)
      
      threshold_grid  =  seq(lower,upper,length = 1000)
      Expected_false_alarms = rep(0,length(threshold_grid ))
      
      for(jj in 1:length(threshold_grid))
      {
        aux = rep(0,NN)
        for(ii in 1:NN)
        {
          aux[ii] =    length(which(ks[ii,]> threshold_grid[jj]))
        }
        Expected_false_alarms[jj] =  mean(aux)     
      }
      
      jjj =  which.min(abs(Expected_false_alarms - 1))
      
      k_alpha_grid[1]  =  threshold_grid[jjj]
      
      ################################################################      
      ### PKS tres_hold
      #       
      k_alpha_grid_old[1]  = 0
      for(ii in 1:NN)
      {
        pks[ii, ]  = choose_tres_hold_PKS(N = 1000,mu,sigma,total_rate, size_prob,L)
      }
      
      lower = min(pks)
      upper = max(pks)
      
      threshold_grid  =  seq(lower,upper,length = 1000)
      Expected_false_alarms = rep(0,length(threshold_grid ))
      
      for(jj in 1:length(threshold_grid))
      {
        aux = rep(0,NN)
        for(ii in 1:NN)
        {
          aux[ii] =    length(which(pks[ii,]> threshold_grid[jj]))
        }
        Expected_false_alarms[jj] =  mean(aux)     
      }
      
      jjj =  which.min(abs(Expected_false_alarms - 1))
      
      k_alpha_grid_old[1]  =  threshold_grid[jjj]
      
      
      ################################################################      
      ### log_Rnf tresd_hold
      
      for(kk in 1:length(tau_grid))
      {
        for(ii in 1:NN)
        {
         rnf[kk,ii,] =  choose_tres_hold_log_Rnf(N = 1000,mu,sigma,total_rate,size_prob,L,tau_grid[kk])
        }
        
        lower = min(rnf)
        upper = max(rnf)
        
        threshold_grid  =  seq(lower,upper,length = 4000)
        Expected_false_alarms = rep(0,length(threshold_grid ))
        
        for(jj in 1:length(threshold_grid))
        {
          aux = rep(0,NN)
          for(ii in 1:NN)
          {
            aux[ii] =    length(which(rnf[kk,ii,]> threshold_grid[jj]))
          }
          Expected_false_alarms[jj] =  mean(aux)     
        }
        
        jjj =  which.min(abs(Expected_false_alarms - 1))
        
        log_A_grid[kk]  =  threshold_grid[jjj]
        
      }
      
      ################################################################      
      ### mle
      log_A_grid_mle[1]  =  0
      for(ii in 1:NN)
      {
          mle_tr[ii,] = choose_tres_hold_log_mle(N = 1000,mu,sigma,total_rate,size_prob,L)
      } 
      
      
      lower = min(mle_tr)
      upper = max(mle_tr)
      
      threshold_grid  =  seq(lower,upper,length = 1000)
      Expected_false_alarms = rep(0,length(threshold_grid ))
      
      for(jj in 1:length(threshold_grid))
      {
        aux = rep(0,NN)
        for(ii in 1:NN)
        {
          aux[ii] =    length(which(mle_tr[ii,]> threshold_grid[jj]))
        }
        Expected_false_alarms[jj] =  mean(aux)     
      }
      
      jjj =  which.min(abs(Expected_false_alarms - 1))
      
      log_A_grid_mle[1]  =  threshold_grid[jjj]
      
      
      ######################################    
      
      false_alarms_PKS = matrix(0,NMC,length(k_alpha_grid_old))
      false_alarms_KS = matrix(0,NMC,length(k_alpha_grid))
      false_alarms_Rnf = matrix(0,NMC,length(log_A_grid))
      false_alarms_mle = matrix(0,NMC,length(log_A_grid_mle))
      false_alarms_KS_star = matrix(0,NMC,length(k_alpha_grid_star))
      
      delay_time_KS =   matrix(0,NMC,length(k_alpha_grid)) 
      delay_time_PKS =   matrix(0,NMC,length(k_alpha_grid_old)) 
      delay_time_Rnf = matrix(0,NMC,length(log_A_grid))
      delay_time_mle = matrix(0,NMC,length(log_A_grid_mle))
      delay_time_KS_star =   matrix(0,NMC,length(k_alpha_grid_star)) 
      
      
      v_array = rep(0,NMC )
      
      for(iter in 1:NMC)
      { 
        if(iter %%20 ==0)
        {
          print("iter")
          print(iter) 
        }
        
 
        KS = rep(0,T) 
        PKS = rep(0,T) 
        log_Rnf = matrix(0,T,length(log_A_grid)) 
        log_mle = rep(0,T)
        KS_star = rep(0,T) 
        
        x =  matrix(0,T,total_rate)
        
        v = 1001
        v_array[iter] = v
        
        
        ptm <- proc.time()
        for(t in 1:T)
        {
          if(t <= v) 
          { 
            x[t, ] =  rnorm(total_rate,mu,sigma)
          }  
          if(t > v)
          { 
            x[t,] = rnorm(total_rate,mu_c,sigma)
          } #rpois(length(theta0),lambdac_t)}
          # ycounts =  gridCounts(x[t,],d)
          
          
          ## KS
          KS[t] = sequential_KS_statistics(x,t,mu,sigma,L)
          
          if(t==1)
          {
            PKS[t] = total_rate*max(pnorm(x[iter,],mu,sigma) -   (1:total_rate)/total_rate)
          }
          if(t>1)
          {
            xx = sort(as.vector(x[1:t,])) 
            PKS[t] = length(xx)*max(pnorm(xx,mu,sigma)  - (1:length(xx))/length(xx))
          }
          # log_rnf
          for(jj in 1:length(log_A_grid))
          {
            log_Rnf[t,jj] = logRnf(t,tau_grid[jj],x,L,mu,sigma,total_rate)  
          }
          
          ### mle
          log_mle[t]  =  logmle(t,x,L,mu,sigma,total_rate)
        }## close for 1:T      
        proc.time() - ptm
        
        ## update false alarms
        
        for(j in 1:length(k_alpha_grid))
        {
          false_alarms_KS[iter,j] = length(which( KS[1:v]> k_alpha_grid[j]))
        }
        
        for(j in 1:length(k_alpha_grid))
        {
          false_alarms_KS_star[iter,j] = length(which( KS[1:v]> k_alpha_grid_star[j]))
        }
        
        for(j in 1:length(k_alpha_grid_old))
        {
          false_alarms_PKS[iter,j] = length(which( PKS[1:v]> k_alpha_grid_old[j]))
        }
        
        for(j in 1:length(log_A_grid))
        {
          false_alarms_Rnf[iter,j] = length(which( log_Rnf[1:v,j]> log_A_grid[j]  ))
        } 
        
        for(j in 1:length(log_A_grid_mle))
        {
          false_alarms_mle[iter,j] = length(which( log_mle[1:v]> log_A_grid_mle[j]  ))
        } 
        ## update delayed time            
        
        
        
        for(j in 1:length(k_alpha_grid))
        {
          ind2 = which(KS[(1+v):T] > k_alpha_grid[j])
          delay_time_KS[iter ,j] = T - (v+1)
          
          if( length(ind2) > 0)
          {
            delay_time_KS[iter,j ] = min(ind2) -1
          }
        }
        for(j in 1:length(k_alpha_grid_star))
        {
          ind2 = which(KS[(1+v):T] > k_alpha_grid_star[j])
          delay_time_KS_star[iter ,j] = T - (v+1)
          
          if( length(ind2) > 0)
          {
            delay_time_KS_star[iter,j ] = min(ind2) -1
          }
        }
        
        for(j in 1:length(k_alpha_grid_old))
        {
          ind2 = which(PKS[(1+v):T] > k_alpha_grid_old[j])
          delay_time_PKS[iter ,j] = T - (v+1)
          
          if( length(ind2) > 0)
          {
            delay_time_PKS[iter,j ] = min(ind2) -1
          }
        }
        
        for(j in 1:length(log_A_grid))
        {
          ind2 = which(log_Rnf[(1+v):T,j] > log_A_grid[j])
          delay_time_Rnf[iter ,j] = T - (v+1)
          if( length(ind2) > 0)
          {
            delay_time_Rnf[iter,j ] = min(ind2) -1
          }
        }
        for(j in 1:length(log_A_grid_mle))
        {
          ind2 = which(log_mle[(1+v):T] > log_A_grid_mle[j])
          delay_time_mle[iter ,j] = T - (v+1)
          if( length(ind2) > 0)
          {
            delay_time_mle[iter,j ] = min(ind2) -1
          }
        }
        
        
      }## close for NMC simulations, iter
      
      for(j in 1:length(k_alpha_grid))
      {
        average_delay_time_KS[ind_T,ind_d,ind_TRF,j] = mean( delay_time_KS[,j])
        average_false_alarms_KS[ind_T,ind_d,ind_TRF,j] = mean(false_alarms_KS[,j])
      }
      
      for(j in 1:length(k_alpha_grid))
      {
        average_delay_time_KS_star[ind_T,ind_d,ind_TRF,j] = mean( delay_time_KS_star[,j])
        average_false_alarms_KS_star[ind_T,ind_d,ind_TRF,j] = mean(false_alarms_KS_star[,j])
      }
      
      for(j in 1:length(k_alpha_grid_old))
      {
        average_delay_time_PKS[ind_T,ind_d,ind_TRF,j] = mean( delay_time_PKS[,j])
        average_false_alarms_PKS[ind_T,ind_d,ind_TRF,j] = mean(false_alarms_PKS[,j])
      }
      for(j in 1:length(log_A_grid))
      {
        average_delay_time_Rnf[ind_T,ind_d,ind_TRF,j] = mean( delay_time_Rnf[,j])
        average_false_alarms_Rnf[ind_T,ind_d,ind_TRF,j] = mean(false_alarms_Rnf[,j])
      }
      
      for(j in 1:length(log_A_grid_mle))
      {
        average_delay_time_mle[ind_T,ind_d,ind_TRF,j] = mean( delay_time_mle[,j])
        average_false_alarms_mle[ind_T,ind_d,ind_TRF,j] = mean(false_alarms_mle[,j])
      }
      
      print("average v")
      print(mean(v_array))
      print("average_delay_time_KS")
      print(average_delay_time_KS[ind_T,ind_d,ind_TRF,])
      print("average_delay_time_KS_star")
      print(average_delay_time_KS_star[ind_T,ind_d,ind_TRF,])
      print("average_delay_time_PKS")
      print(average_delay_time_PKS[ind_T,ind_d,ind_TRF,])
      print("average_delay_time_Rnf")
      print(average_delay_time_Rnf[ind_T,ind_d,ind_TRF,] )
      print("average_delay_time_mle")
      print(average_delay_time_mle[ind_T,ind_d,ind_TRF,] )
     
      
       print(" average_false_alarms_KS")
      print( average_false_alarms_KS[ind_T,ind_d,ind_TRF,])
      print(" average_false_alarms_PKS")
      print( average_false_alarms_PKS[ind_T,ind_d,ind_TRF,])
      print("average_false_alarms_Rnf")
      print(average_false_alarms_Rnf[ind_T,ind_d,ind_TRF,] )
      print("average_false_alarms_mle")
      print(average_false_alarms_mle[ind_T,ind_d,ind_TRF,] )
      
    }## close for Total rate
  }## close for d
}## close for T


ind_TRF = 3
print("average_delay_time_KS")
print(average_delay_time_KS[ind_T,ind_d,ind_TRF,])
print("average_delay_time_Rnf")
print(average_delay_time_Rnf[ind_T,ind_d,ind_TRF,] )
print("average_delay_time_mle")
print(average_delay_time_mle[ind_T,ind_d,ind_TRF,] )
print(" average_false_alarms_KS")
print( average_false_alarms_KS[ind_T,ind_d,ind_TRF,])
print("average_false_alarms_Rnf")
print(average_false_alarms_Rnf[ind_T,ind_d,ind_TRF,] )
print("average_false_alarms_mle")
print(average_false_alarms_mle[ind_T,ind_d,ind_TRF,] )