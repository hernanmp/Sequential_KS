###  Code that can be used obtianing results in Table 1 in the paper.
###  the average delay time  for method X is saved to array:
###  "average_delay_time_X".  From these arrays one can manually get the results in Table 1
###  The dimensions of average_delay_time_X are:
### number of possible time horizons (we use only one) for testing methods
### number of possible lengths of count vectors 
### number of possible total rates of photon counts
### number of different choices of hyperparameters for method X


rm(list = ls())
working_directory = "/Users/user/Desktop/JASA_code_Anomaly/Experiments"
setwd(paste(working_directory,"/Simulated_example",sep=""))



## read parameters of experiment
mu = read.table("example2_L50_mu.txt")
mu = mu[,1]
mu
tau2 =  read.table("example2_L50_tau2.txt")
tau2 = tau2[,1]
tau2
pi =  read.table("example2_L50_pi.txt")
pi = pi[,1]
pi
weight = read.table("example2_L50_weight.txt")
weight = weight[1,1]
weight
mu_c = read.table("example2_L50_mu_c.txt")
mu_c = mu_c[1,1]
mu_c
sigma_c = read.table("example2_L50_sigma_c.txt")
sigma_c = sigma_c[1,1] 
sigma_c
theta0 =  read.table("example2_L50_theta0.txt")
theta0 =  theta0[,1]
plot(theta0)
thetac = read.table("example2_L50_thetac.txt")
thetac = thetac[,1]
plot(thetac)



setwd(paste(working_directory,"/Code",sep=""))

### 2^d is the number of bins
d = 11


source("Generate_data.R")
matplot(cbind(theta0,thetac),type="l")

####  setting parameters for simulations
####  setting parameters for simulations
##  Time horizon
T_grid   =   c(700)
###  the number of bins is 2^d
d_grid = c(9,10,11)
###  Average number of pothon measurements accros bins
Total_rate_factor_grid = c(100/2^d,500/2^d,1000/2^d) 
### Number of Monte Carlo simulations
NMC = 100


L_grid = c(50)  ## window length for KS
k_alpha = 3 ## KS

k_alpha_grid = c(1.36) ## this is changed later
k_alpha_grid_star  = c(1.36)
level_grid_dir = c(.5) 
log_A_grid = c(0) ## this is changed later
k_alpha_grid_old = c(1.36)
tresd_hold_mle_grid = c(0)

size_prob  = .999 ### 1 - size_prob  is the probability of a false alarm


######################################################### ######################################   
ind_T  = 1
ind_d  = 1
ind_TRF = 1
T = T_grid[ind_T]
d = d_grid[ind_d]


###  Arrays to save average false alarms rates
##### Predecesor KS
average_false_alarms_KS_old = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid_old),length(L_grid)))
##### Proposed KS
average_false_alarms_KS = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid),length(L_grid)))
##### Priro based method
average_false_alarms_Rnf = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(log_A_grid),length(L_grid)))
##### Mle based method
average_false_alarms_mle = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(tresd_hold_mle_grid),length(L_grid)))
##### Proposed KS with Corollary 2 rule
average_false_alarms_KS_star = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid_star),length(L_grid)))
##### SCR method
average_false_alarms_scr = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid),length(L_grid)))


### Arrays to save average delay times
##### Mle based method
average_delay_time_mle = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(tresd_hold_mle_grid),length(L_grid)))
##### Predecesor KS
average_delay_time_KS = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid),length(L_grid)))
##### Predecesor KS
average_delay_time_KS_old = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid_old),length(L_grid)))
average_delay_time_Rnf = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(log_A_grid),length(L_grid)))
##### Proposed KS with Corollary 2 rule
average_delay_time_KS_star = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid_star),length(L_grid)))
##### SCR method
average_delay_time_scr = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid),length(L_grid)))

## Number of training data sets 
N_training =   100
## time horizon for training
T_training = 1000  #


for(ind_L in 1:length(L_grid))
{
  L = L_grid[ind_L]
  
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
        
        loc =  seq(-8,8,length = m)
        theta0 =  pi[1]*dnorm(loc,mu[1],sqrt(tau2[1])) + pi[2]*dnorm(loc,mu[2],sqrt(tau2[2])) + pi[3]*dnorm(loc,mu[3],sqrt(tau2[3]))+pi[4]*dnorm(loc,mu[4],sqrt(tau2[4]))
        thetac = weight*theta0 + (1-weight)*dnorm(loc,mu_c,sigma_c)
        matplot(cbind(theta0,thetac),type="l")
        theta0 = theta0/sum(theta0)
        thetac = thetac/sum(thetac)
        
        
        lambda0 = Total_rate_factor*m*theta0
        lambdac = Total_rate_factor*m*thetac
        
        
        ###  Training period
        temp = Training_period(ntraining = T_training,m,d,lambda0) 
        #  F0 = temp$F0
        alpha = temp$alpha
        
        total_rate = temp$total_rate
        
        theta0_hat = temp$theta0
        F0 =  cumsum(theta0)
        ###########################
        
        ################################################################  
        kappa_grid = 10^seq(-5,-1,length= 20)
        
        ################################################################        
        ## KS^* threshold, based on Corrollary 2
        k_alpha_grid_star[1]   = sqrt(log(2*1000*(L+1))/2)
        
        ################################################################      
        ### KS threshold, proposed approach,  callibrating with MC simulations
        
        k_alpha_grid[1]  = 0
        KS = matrix(0,N_training,T_training)
        for(i in 1:N_training)
        {
          ## computing test statistic for trial i
          KS[i,] =   new_choose_tres_hold_KS(N = T_training,theta0 ,total_rate,d,F0, size_prob,L)
        }
        ### now we compute the expected false alarm rate for different thresholds
        lower = min(KS)
        upper = max(KS)
        
        threshold_grid  =  seq(lower,upper,length = 1000)
        Expected_false_alarms = rep(0,length(threshold_grid ))
        
        for(j in 1:length(threshold_grid))
        {
          ### array to store expected false alarm rates, when
          ### using threshold_grid[i] as threshold
          aux = rep(0,N_training)
          for(i in 1:N_training)
          {
            aux[i] =    length(which(KS[i,]> threshold_grid[j]))
          }
          Expected_false_alarms[j] =  mean(aux)     
        }
        j_best =  which.min(abs(Expected_false_alarms - 1))
        
        k_alpha_grid[1]  =  threshold_grid[j_best]
        
        ################################################################      
        ## Ks old
        ### Precusor KS tres_hold
        k_alpha_grid_old[1]  = 0
        
        KS = matrix(0,N_training,T_training)
        for(i in 1:N_training)
        {
          ### compute, for differnt MC, the statistic of the precusor ks
          KS[i,] = new_choose_tres_hold_KS_old(1,N = T_training,theta0,total_rate,d,F0, size_prob,L)
        }
        ### now we compute the expected false alarm rate for different thresholds
        lower = min(KS)
        upper = max(KS)
        
        threshold_grid  =  seq(lower,upper,length = 1000)
        Expected_false_alarms = rep(0,length(threshold_grid ))
        
        for(j in 1:length(threshold_grid))
        {
          ### array to store expected false alarm rates, when
          ### using threshold_grid[i] as threshold
          aux = rep(0,N_training)
          for(i in 1:N_training)
          {
            aux[i] =    length(which(KS[i,]> threshold_grid[j]))
          }
          Expected_false_alarms[j] =  mean(aux)     
        }
        
        j_best =  which.min(abs(Expected_false_alarms - 1))
        
        k_alpha_grid_old[1]  =  threshold_grid[j_best]  
        
        
        ################################################################      
        ### log_Rnf (method integrating prior)  tresd_hold 
        ## looping over different prior parameters
        log_A_grid[1]  =  0
        log_Rnf  = matrix(0,N_training,T_training)
        for(i in 1:N_training)
        {
          ## compute statistic
          log_Rnf[i,] = new_choose_tres_hold_log_Rnf(N = T_training,theta0,sum(lambda0),d, size_prob,L)
        }    
        
        ## compute expected false alarms rate for different thresholds
        lower = min( log_Rnf)
        upper = max( log_Rnf)
        
        threshold_grid  =  seq(lower,upper,length = 1000)
        Expected_false_alarms = rep(0,length(threshold_grid ))
        
        for(j in 1:length(threshold_grid))
        {
          aux = rep(0,N_training)
          for(i in 1:N_training)
          {
            aux[i] =    length(which(log_Rnf[i,]> threshold_grid[j]))
          }
          Expected_false_alarms[j] =  mean(aux)     
        }
        
        j_best =  which.min(abs(Expected_false_alarms - 1))
        
        log_A_grid[1] =  threshold_grid[j_best]  
        
        
        ################################################################        
        ### mle,   method based on generalized likelihood ration statistic
        tresd_hold_mle_grid[1]   = 0
        mle  = matrix(0,N_training,T_training)
        for(i in 1:N_training)
        {
          ## compute statistic for different MC simulations
          mle[i,] =    new_choose_tres_hold_mle(N = T_training,theta0,lambda0,d, size_prob,L)
        }   
        ##  Next choosing threshold
        
        lower = min( mle)
        upper = max( mle)
        
        threshold_grid  =  seq(lower,upper,length = 1000)
        Expected_false_alarms = rep(0,length(threshold_grid ))
        
        for(j in 1:length(threshold_grid))
        {
          ###  aray to store false alarm rates
          aux = rep(0,N_training)
          for(i in 1:N_training)
          {
            aux[i] =    length(which(mle[i,]> threshold_grid[j]))
          }
          Expected_false_alarms[j] =  mean(aux)     
        }
        
        j_best =  which.min(abs(Expected_false_alarms - 1))
        tresd_hold_mle_grid[1] =  threshold_grid[j_best]  
        
        ######################################    
        #### arrays to store false alarms rates for different MC tests
        
        false_alarms_KS = matrix(0,NMC,length(k_alpha_grid))
        false_alarms_KS_old = matrix(0,NMC,length(k_alpha_grid_old))
        false_alarms_Rnf = matrix(0,NMC,length(log_A_grid))
        false_alarms_KS_star = matrix(0,NMC,length(k_alpha_grid_star))
        
        #### arrays to store delay times for different MC tests
        delay_time_KS =   matrix(0,NMC,length(k_alpha_grid)) 
        delay_time_KS_star =   matrix(0,NMC,length(k_alpha_grid_star)) 
        delay_time_KS_old =   matrix(0,NMC,length(k_alpha_grid_old)) 
        delay_time_Rnf = matrix(0,NMC,length(log_A_grid))
        false_alarms_mle = matrix(0,NMC,length(tresd_hold_mle_grid))
        delay_time_mle =   matrix(0,NMC,length(tresd_hold_mle_grid)) 
        
        ###vector to store change point locations
        v_array = rep(0,NMC )
        
        for(iter in 1:NMC)
        { 
          if(iter %%20 ==0)
          {
            print("iter")
            print(iter) 
          }
          ###  arrays to store statistics for current trial
          KS = rep(0,T) 
          KS_old = rep(0,T) 
          log_Rnf = rep(0,T) 
          mle = rep(0,T) 
          
          ### array to store photon measurements
          x =  matrix(0,T,m)
          
          v = min(floor( runif(1)*((T-1-100)-100) + 100    ),T)
          v_array[iter] = v
          
          cum_x = rep(0,length(lambda0))
          
          
          ptm <- proc.time()
          for(t in 1:T)
          {
            lambda0_t = lambda0 
            
            lambdac_t = lambdac 
            
            if(t <= v) 
            { 
              x[t, ] =rpois(length(theta0), lambda0_t )
            }  
            if(t > v)
            { 
              x[t,] = rpois(length(theta0),lambdac_t)
            } 
            
            cum_x = cum_x + x[t,]
            
            KS_old[t] =  sqrt(sum(cum_x))*kolmogorov_statistic(cum_x,F0)
            
            #
            log_Rnf[t] = logRnf(t,1,1,x,L,lambda0)
            
            mle[t] = mle_detection(x,t,L,lambda0)
            
            KS[t] = sequential_KS_statistics(x,t,F0,L)
            
            
            
          }## close for 1:T      
          proc.time() - ptm
          
          ## update false alarms
          
          for(j in 1:length(k_alpha_grid))
          {
            false_alarms_KS[iter,j] = length(which( KS[1:v]> k_alpha_grid[j]))
          }
          
          for(j in 1:length(k_alpha_grid_star))
          {
            false_alarms_KS_star[iter,j] = length(which( KS[1:v]> k_alpha_grid_star[j]))
          }
          
          for(j in 1:length(k_alpha_grid_old))
          {
            false_alarms_KS_old[iter,j] = length(which( KS_old[1:v]> k_alpha_grid_old[j]))
          }
          
          for(j in 1:length(log_A_grid))
          {
            false_alarms_Rnf[iter,j] = length(which( log_Rnf[1:v]> log_A_grid[j]  ))
          } 
          for(j in 1:length(tresd_hold_mle_grid))
          {
            false_alarms_mle[iter,j] = length(which( mle[1:v]> tresd_hold_mle_grid[j]))
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
            ind2 = which(KS_old[(1+v):T] > k_alpha_grid_old[j])
            delay_time_KS_old[iter ,j] = T - (v+1)
            
            if( length(ind2) > 0)
            {
              delay_time_KS_old[iter,j ] = min(ind2) -1
            }
          }
          for(j in 1:length(tresd_hold_mle_grid  ))
          {
            ind2 = which(mle[(1+v):T] > tresd_hold_mle_grid[j])
            delay_time_mle[iter ,j] = T - (v+1)
            
            if( length(ind2) > 0)
            {
              delay_time_mle[iter,j ] = min(ind2) -1
            }
          } 
          for(j in 1:length(log_A_grid))
          {
            ind2 = which(log_Rnf[(1+v):T] > log_A_grid[j])
            delay_time_Rnf[iter ,j] = T - (v+1)
            if( length(ind2) > 0)
            {
              delay_time_Rnf[iter,j ] = min(ind2) -1
            }
          }
          
        }## close for NMC simulations, iter
        #### averaging delay times and false alarm rates over different MC simulations
        
        
        for(j in 1:length(k_alpha_grid_old))
        {
          average_delay_time_KS_old[ind_T,ind_d,ind_TRF,j,ind_L] = mean( delay_time_KS_old[,j])
          average_false_alarms_KS_old[ind_T,ind_d,ind_TRF,j,ind_L] = mean(false_alarms_KS_old[,j])
        }
        for(j in 1:length(k_alpha_grid))
        {
          average_delay_time_KS[ind_T,ind_d,ind_TRF,j,ind_L] = mean( delay_time_KS[,j])
          average_false_alarms_KS[ind_T,ind_d,ind_TRF,j,ind_L] = mean(false_alarms_KS[,j])
        }
        for(j in 1:length(k_alpha_grid_star))
        {
          average_delay_time_KS_star[ind_T,ind_d,ind_TRF,j,ind_L] = mean( delay_time_KS_star[,j])
          average_false_alarms_KS_star[ind_T,ind_d,ind_TRF,j,ind_L] = mean(false_alarms_KS_star[,j])
        }
        for(j in 1:length( tresd_hold_mle_grid ))
        {
          average_delay_time_mle[ind_T,ind_d,ind_TRF,j,ind_L] = mean( delay_time_mle[,j])
          average_false_alarms_mle[ind_T,ind_d,ind_TRF,j,ind_L] = mean(false_alarms_mle[,j])
        }
        for(j in 1:length(log_A_grid))
        {
          average_delay_time_Rnf[ind_T,ind_d,ind_TRF,j,ind_L] = mean( delay_time_Rnf[,j])
          average_false_alarms_Rnf[ind_T,ind_d,ind_TRF,j,ind_L] = mean(false_alarms_Rnf[,j])
        }
        
        
        print("L")
        print(L)
        print("average v")
        print(mean(v_array))
        print("average_delay_time_KS")
        print(average_delay_time_KS[ind_T,ind_d,ind_TRF,,ind_L])
        print("average_delay_time_KS_star")
        print(average_delay_time_KS_star[ind_T,ind_d,ind_TRF,,ind_L])
        print("average_delay_time_KS_old")
        print(average_delay_time_KS_old[ind_T,ind_d,ind_TRF,,ind_L] )
        print("average_delay_time_Rnf")
        print(average_delay_time_Rnf[ind_T,ind_d,ind_TRF,,ind_L] )
        print("average_delay_time_mle")
        print(average_delay_time_mle[ind_T,ind_d,ind_TRF,,ind_L])
        
        print(" average_false_alarms_KS")
        print( average_false_alarms_KS[ind_T,ind_d,ind_TRF,,ind_L])
        print(" average_false_alarms_KS_star")
        print( average_false_alarms_KS_star[ind_T,ind_d,ind_TRF,,ind_L])
        print("average_false_alarms_KS_old")
        print(average_false_alarms_KS_old[ind_T,ind_d,ind_TRF,,ind_L] )
        print("average_false_alarms_Rnf")
        print(average_false_alarms_Rnf[ind_T,ind_d,ind_TRF,,ind_L] )
        print(" average_false_alarms_mle")
        print( average_false_alarms_mle[ind_T,ind_d,ind_TRF,,ind_L])
        
      }## close for Total rate
    }## close for d
  }## close for T
  
}### close for L



ind_TRF = 3
ind_d = 3
print("L")
print(L)
print("average v")
print(mean(v_array))
print("average_delay_time_KS")
print(average_delay_time_KS[ind_T,ind_d,ind_TRF,,ind_L])
print("average_delay_time_KS_star")
print(average_delay_time_KS_star[ind_T,ind_d,ind_TRF,,ind_L])
print("average_delay_time_KS_old")
print(average_delay_time_KS_old[ind_T,ind_d,ind_TRF,,ind_L] )
print("average_delay_time_Rnf")
print(average_delay_time_Rnf[ind_T,ind_d,ind_TRF,,ind_L] )
print("average_delay_time_mle")
print(average_delay_time_mle[ind_T,ind_d,ind_TRF,,ind_L])

print(" average_false_alarms_KS")
print( average_false_alarms_KS[ind_T,ind_d,ind_TRF,,ind_L])
print("average_false_alarms_KS_old")
print(average_false_alarms_KS_old[ind_T,ind_d,ind_TRF,,ind_L] )
print("average_false_alarms_Rnf")
print(average_false_alarms_Rnf[ind_T,ind_d,ind_TRF,,ind_L] )
print(" average_false_alarms_mle")
print( average_false_alarms_mle[ind_T,ind_d,ind_TRF,,ind_L])
