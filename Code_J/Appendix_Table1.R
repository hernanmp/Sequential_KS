
###  code for obtianing results in Table 1 in the appendix
###  the average delay time  for method X is saved to array:
###  "average_delay_time_X".  From these arrays one can get the results in Table 1
###  The dimensions of average_delay_time_X are:
### number of possible time horizons (we only use  one) for testing methods
### number of possible lengths of count vectors 
### number of possible total rates of photon counts
### number of different choices of hyperparameters for method X


rm(list = ls())

working_directory = "/Users/user/Desktop/JASA_code_Anomaly/Experiments/Code"
setwd(working_directory )

## loading all the auxiliary functions
source("utilis.R")



##  mean after change point
mu_c = .3;
### mean before change point
mu = 0

##  standard deviation before change point
sigma= 6
##  standard deviation after change point
sigma_c = 6

tau_grid = c(.1,1,5,10)



####  setting parameters for simulations
##  Time horizon
T_grid   =   c(1100)
###  the number of bins is 2^d
d_grid = c(11)
d=  11
###  Average number of pothon measurements accros bins
Total_rate_factor_grid = c(100/2^d,500/2^d,1000/2^d) 

### Number of Monte Carlo simulations
NMC = 100


L = 50   ## window length for KS
k_alpha = 3 ## KS

k_alpha_grid = c(1.36) ## array for KS threshold
k_alpha_grid_old = k_alpha_grid ## array for precursor KS threshold
k_alpha_grid_star  = k_alpha_grid ## array for KS threshold 
level_grid =  c(.9,.7,.5,.3) #c(.95,.99)
level_grid_dir = c(.9,.7,.5,.3) 
log_A_grid = rep(0,length(tau_grid)) ## this is changed later
log_A_grid_mle = c(0) 

size_prob  = .999### 1 - size_prob  is the probability of a false alarm


#tres_hold_KS = rep(0,length(k_alpha_grid))
#tres_hold_log_Rnf = rep(0,length(log_A_grid))


#v = min(floor(runif(1)*T),T)
######################################################### ######################################   
ind_T  = 1
ind_d  = 1
ind_TRF = 1
T = T_grid[ind_T]
d = d_grid[ind_d]

###  Arrays to save average false alarms rates
average_false_alarms_KS = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid)))
average_false_alarms_PKS = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid_old)))
average_false_alarms_Rnf = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(log_A_grid)))
average_false_alarms_mle = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(log_A_grid_mle)))
average_false_alarms_KS_star = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid_star)))


### Arrays to save average delay times
average_delay_time_KS = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid)))
average_delay_time_Rnf = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(log_A_grid)))
average_delay_time_mle = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(log_A_grid_mle)))
average_delay_time_PKS = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid_old)))
average_delay_time_KS_star = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid_star)))

## Number of training data sets 
N_training =   100
## time horizon for training
T_training = 1000  #

### Arrays to save statistics for training different methods
ks =  matrix(0,N_training,T_training)
rnf = array(0,c(length(tau_grid),N_training,T_training))
mle_tr = matrix(0,N_training,T_training)
pks = matrix(0,N_training,T_training)

for(ind_T in 1:length(T_grid))
{
  ##  Time horizon for evaluating performance
  T = T_grid[ind_T]
  ######################################  
  for(ind_d in 1:length(d_grid))
  {
    ####  2^d  is the number of bins
    d = d_grid[ind_d]
    ######################################  
    for(ind_TRF in 1:length(Total_rate_factor_grid))
    {
      ### average number of photons per channel
      Total_rate_factor = Total_rate_factor_grid[ind_TRF ]
      
      ###  True normalized pre and post densities 
      ## pre change
      m = 2^d
      
      total_rate = m*Total_rate_factor
      
      ###  Training period
      ################################################################        
      ## KS^* threshold, based on Corrollary 2
      k_alpha_grid_star[1]   = sqrt(log(2*1000*L)/2)
      ################################################################      
      ### KS threshold, proposed approach,  callibrating with MC simulations
      k_alpha_grid[1]  = 0
      for(i in 1:N_training)
      {
        ## computing test statistic for trial i
        ks[i, ]  = choose_tres_hold_KS(N = T_training,mu,sigma,total_rate, size_prob,L)
      }
      ### now we compute the expected false alarm rate for different thresholds
      lower = min(ks)
      upper = max(ks)
      
      threshold_grid  =  seq(lower,upper,length = 1000)
      Expected_false_alarms = rep(0,length(threshold_grid ))
      
      for(j in 1:length(threshold_grid))
      {
        ### array to store expected false alarm rates, when
        ### using threshold_grid[i] as threshold
        aux = rep(0,N_training)
        for(i in 1:N_training)
        {
          aux[i] =    length(which(ks[i,]> threshold_grid[j]))
        }
        Expected_false_alarms[j] =  mean(aux)     
      }
      j_best =  which.min(abs(Expected_false_alarms - 1))
      k_alpha_grid[1]  =  threshold_grid[j_best]
      
      ################################################################      
      ### Precusor KS tres_hold
      k_alpha_grid_old[1]  = 0
      for(i in 1:N_training)
      {
        ### compute, for differnt MC, the statistic of the precusor ks
        pks[i, ]  = choose_tres_hold_PKS(N = T_training,mu,sigma,total_rate, size_prob,L)
      }
      ### now we compute the expected false alarm rate for different thresholds
      lower = min(pks)
      upper = max(pks)
      
      threshold_grid  =  seq(lower,upper,length = 1000)
      Expected_false_alarms = rep(0,length(threshold_grid ))
      
      for(j in 1:length(threshold_grid))
      {
        ### array to store expected false alarm rates, when
        ### using threshold_grid[i] as threshold
        aux = rep(0,N_training)
        for(i in 1:N_training)
        {
          aux[i] =    length(which(pks[i,]> threshold_grid[j]))
        }
        Expected_false_alarms[j] =  mean(aux)     
      }
      
      j_best =  which.min(abs(Expected_false_alarms - 1))
      
      k_alpha_grid_old[1]  =  threshold_grid[j_best]
      
      
      ################################################################      
      ### log_Rnf (method integrating prior)  tresd_hold 
      ## looping over different prior parameters
      for(k in 1:length(tau_grid))
      {
        for(i in 1:N_training)
        {
          ## compute statistic
          rnf[k,i,] =  choose_tres_hold_log_Rnf(N = T_training,mu,sigma,total_rate,size_prob,L,tau_grid[k])
        }
        
        ## compute expected false alarms rate for different thresholds
        
        lower = min(rnf[k,,])
        upper = max(rnf[k,,])
        
        threshold_grid  =  seq(lower,upper,length = 4000)
        Expected_false_alarms = rep(0,length(threshold_grid ))
        
        for(j in 1:length(threshold_grid))
        {
          aux = rep(0,N_training)
          for(i in 1:N_training)
          {
            aux[i] =    length(which(rnf[k,i,]> threshold_grid[j]))
          }
          Expected_false_alarms[j] =  mean(aux)     
        }
        
        j_best =  which.min(abs(Expected_false_alarms - 1))
        
        log_A_grid[k]  =  threshold_grid[j_best]
        
      }
      
      ################################################################      
      ### mle,   method based on generalized likelihood ration statistic
      log_A_grid_mle[1]  =  0
      for(i in 1:N_training)
      {
        ## compute statistic for different MC simulations
        mle_tr[i,] = choose_tres_hold_log_mle(N = T_training,mu,sigma,total_rate,size_prob,L)
      } 
      
      ##  Next choosing threshold
      
      lower = min(mle_tr)
      upper = max(mle_tr)
      
      threshold_grid  =  seq(lower,upper,length = 1000)
      Expected_false_alarms = rep(0,length(threshold_grid ))
      
      for(j in 1:length(threshold_grid))
      {
        ###  aray to store false alarm rates
        aux = rep(0,N_training)
        for(i in 1:N_training)
        {
          aux[i] =    length(which(mle_tr[i,]> threshold_grid[j]))
        }
        Expected_false_alarms[j] =  mean(aux)     
      }
      
      j_best =  which.min(abs(Expected_false_alarms - 1))
      
      log_A_grid_mle[1]  =  threshold_grid[j_best]
      
      
      ######################################    
      
      #### arrays to store false alarms rates for different MC tests
      false_alarms_PKS = matrix(0,NMC,length(k_alpha_grid_old))
      false_alarms_KS = matrix(0,NMC,length(k_alpha_grid))
      false_alarms_Rnf = matrix(0,NMC,length(log_A_grid))
      false_alarms_mle = matrix(0,NMC,length(log_A_grid_mle))
      false_alarms_KS_star = matrix(0,NMC,length(k_alpha_grid_star))
      
      #### arrays to store delay times for different MC tests
      delay_time_KS =   matrix(0,NMC,length(k_alpha_grid)) 
      delay_time_PKS =   matrix(0,NMC,length(k_alpha_grid_old)) 
      delay_time_Rnf = matrix(0,NMC,length(log_A_grid))
      delay_time_mle = matrix(0,NMC,length(log_A_grid_mle))
      delay_time_KS_star =   matrix(0,NMC,length(k_alpha_grid_star)) 
      
      ###vector to store change point locations
      v_array = rep(0,NMC )
      
      for(iter in 1:NMC)
      { 
        if(iter %%20 ==0)
        {
          print("iter")
          print(iter) 
        }
        
        ###  array to store statistics for current trial
        KS = rep(0,T) 
        PKS = rep(0,T) 
        log_Rnf = matrix(0,T,length(log_A_grid)) 
        log_mle = rep(0,T)
        KS_star = rep(0,T) 
        
        ### array to store photon measurements
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
          }
          
          ## KS
          KS[t] = sequential_KS_statistics(x,t,mu,sigma,L)
          
          if(t==1)
          {
            PKS[t] = total_rate*max(pnorm(x[1,],mu,sigma) -   (1:total_rate)/total_rate)
          }
          if(t>1)
          {
            xx = sort(as.vector(x[1:t,])) 
            PKS[t] = length(xx)*max(pnorm(xx,mu,sigma)  - (1:length(xx))/length(xx))
          }
          # log_rnf
          for(j in 1:length(log_A_grid))
          {
            log_Rnf[t,j] = logRnf(t,tau_grid[j],x,L,mu,sigma,total_rate)  
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
      
      #### averaging delay times and false alarm rates over different MC simulations
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