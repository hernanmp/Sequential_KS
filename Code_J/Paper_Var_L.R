

rm(list = ls())

## Code to create the plot of delay time varying the window size, Figure 7 in the main paper
## We start creating the pre and post change densities

working_directory = "/Users/user/Desktop/JASA_code_Anomaly/Experiments"
setwd(paste(working_directory,"/Code",sep=""))
source('utils.R')
setwd(paste(working_directory,"/Data",sep=""))

# How strong is the background rate, and how large is the source in milliCuries?
background_cps = 39  
source_size = 100  #mCi

# Distance to the source
dist_grid = c(50, 100,  150)
i = 2## index of distance used, distance is meassured in meters

# Read in background
background_train = as.numeric(read.csv('background_train_mean.csv', header=FALSE))
plot(background_train, type='l')
n_bins = length(background_train)

cs137_spectrum = as.numeric(read.csv('2013-cs137-5cm.csv', header=FALSE))

# Winsorize: assign all counts after bin 2048 to bin 2048, and then truncate
cs137_spectrum[2048] = sum(cs137_spectrum[2048:4096])
cs137_spectrum = head(cs137_spectrum, 2048)


# Normalize and plot
cs137_spectrum = cs137_spectrum/sum(cs137_spectrum)
plot(cs137_spectrum, type='l')


# Create composite spectrum from background + anomaly
this_dist =  dist_grid[i]
sim_spectrum = inject_source(background_train, background_cps, cs137_spectrum, source_size, this_dist)
sim_spectrum = sim_spectrum/sum(sim_spectrum)

# Plot the spectrum of background + anomaly: this is the "post-changepoint" distribution
this_title = paste0(this_dist, 'm')
plot(background_train, type='l',
     xlim=c(0,1250),
     main=this_title)
lines(sim_spectrum, col='red')


##################################################################################################################
########################################################################################################################
setwd(paste(working_directory,"/Code",sep=""))

##  pre and post change densities
theta0 = background_train
thetac = sim_spectrum
matplot(cbind(theta0,thetac),type ="l")

## number of bins
d = 11

## reading utilities
source("Generate_data.R")


####  setting parameters for simulations
##### Time horizon
T_grid   =   c(700)
##### 2^d is the number of bins or channels
d_grid = c(11)
#####   2^d *  Total_rate_factor_grid  is the average total number of photon counts per second
Total_rate_factor_grid = c(100/2^d) #c(2,3,4)
NMC = 5



L_grid = c(10,20,30,40,50,60,70,80)   ## window length 
k_alpha = 3 ## KS


size_prob  = .999 ### 1 - size_prob  is the probability of a false alarm

#####  Variables to store thresholds for different methods
k_alpha_grid = c(1.36) ## this is changed later
log_A_grid = c(0) ## this is changed later
k_alpha_grid_old = c(1.36)
tresd_hold_mle_grid = c(0)

tres_hold_KS = rep(0,length(k_alpha_grid))
tres_hold_KS_old = rep(0,length(k_alpha_grid_old))
tres_hold_log_Rnf = rep(0,length(log_A_grid))
tresd_hold_mle_grid = c(0)

#v = min(floor(runif(1)*T),T)
######################################################### ######################################   
ind_T  = 1
ind_d  = 1
ind_TRF = 1
T = T_grid[ind_T]
d = d_grid[ind_d]

average_false_alarms_KS_old = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid_old),length(L_grid)))
average_false_alarms_KS = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid),length(L_grid)))
average_false_alarms_Rnf = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(log_A_grid),length(L_grid)))
average_false_alarms_mle = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(tresd_hold_mle_grid),length(L_grid)))


average_delay_time_mle = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(tresd_hold_mle_grid),length(L_grid)))
average_delay_time_KS = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid),length(L_grid)))
average_delay_time_KS_old = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(k_alpha_grid_old),length(L_grid)))
average_delay_time_Rnf = array(0, c(length(T_grid),length(d_grid),length(Total_rate_factor_grid),length(log_A_grid),length(L_grid)))



for(ind_L in 1:length(L_grid))
{
  L = L_grid[ind_L]
  
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
        alpha = temp$alpha
        total_rate = temp$total_rate
        
        theta0_hat = temp$theta0
        F0 =  cumsum(theta0)
        ###########################
        
        ################################################################  
        
        #       N = 500
        kappa_grid = 10^seq(-5,-1,length= 20)
        
        ################################################################      
        ## KS tres_hold
        NT =  100  # number of MC simulations to estimate threshold
        
        k_alpha_grid[1]  = 0
        KSnew = matrix(0,NT,1000)
        for(i in 1:NT)
        {
          KSnew [i,] =   new_choose_tres_hold_KS(N = 1000,theta0,total_rate,d,F0, size_prob,L)
        }
        lower = min(KSnew )
        upper = max(KSnew )
        
        threshold_grid  =  seq(lower,upper,length = 1000)
        Expected_false_alarms = rep(0,length(threshold_grid ))
        
        for(j in 1:length(threshold_grid))
        {
          aux = rep(0,NT)
          for(i in 1:NT)
          {
            aux[i] =    length(which(KSnew[i,]> threshold_grid[j]))
          }
          Expected_false_alarms[j] =  mean(aux)     
        }
        
        j_best =  which.min(abs(Expected_false_alarms - 1))
        k_alpha_grid[1]  =  threshold_grid[j_best]
        
        ################################################################      
        ## Ks old
        
        KS = matrix(0,NT,1000)
        k_alpha_grid_old[1]  = 0
        for(i in 1:NT)
        {
          KS[i,] = new_choose_tres_hold_KS_old(1,N = 1000,theta0,total_rate,d,F0, size_prob,L)
        }
        
        lower = min(KS)
        upper = max(KS)
        
        threshold_grid  =  seq(lower,upper,length = 1000)
        Expected_false_alarms = rep(0,length(threshold_grid ))
        
        for(j in 1:length(threshold_grid))
        {
          aux = rep(0,NT)
          for(i in 1:NT)
          {
            aux[i] =    length(which(KS[i,]> threshold_grid[j]))
          }
          Expected_false_alarms[j] =  mean(aux)     
        }
        
        j_best =  which.min(abs(Expected_false_alarms - 1))
        
        k_alpha_grid_old[1]  =  threshold_grid[j_best]  
        
        
        ################################################################      
        ### log_Rnf tresd_hold
        
        log_A_grid[1]  =  0
        log_Rnf  = matrix(0,NT,1000)
        for(i in 1:NT)
        {
          log_Rnf[i,] = new_choose_tres_hold_log_Rnf(N = 1000,theta0,sum(lambda0),d, size_prob,L)
        }    
        
        lower = min( log_Rnf)
        upper = max( log_Rnf)
        
        threshold_grid  =  seq(lower,upper,length = 1000)
        Expected_false_alarms = rep(0,length(threshold_grid ))
        
        for(j in 1:length(threshold_grid))
        {
          aux = rep(0,NT)
          for(i in 1:NT)
          {
            aux[i] =    length(which(log_Rnf[i,]> threshold_grid[j]))
          }
          Expected_false_alarms[j] =  mean(aux)     
        }
        
        j_best =  which.min(abs(Expected_false_alarms - 1))
        
        log_A_grid[1] =  threshold_grid[j_best]  
        
        
        ################################################################        
        ### MLE tres_hold
        
        tresd_hold_mle_grid[1]   = 0
        mle  = matrix(0,NT,1000)
        for(i in 1:NT)
        {
          mle[i,] =    new_choose_tres_hold_mle(N = 1000,theta0,lambda0,d, size_prob,L)
        }   
        
        lower = min( mle)
        upper = max( mle)
        
        threshold_grid  =  seq(lower,upper,length = 1000)
        Expected_false_alarms = rep(0,length(threshold_grid ))
        
        for(j in 1:length(threshold_grid))
        {
          aux = rep(0,NT)
          for(i in 1:NT)
          {
            aux[i] =    length(which(mle[i,]> threshold_grid[j]))
          }
          Expected_false_alarms[j] =  mean(aux)     
        }
        
        j_best =  which.min(abs(Expected_false_alarms - 1))
        
        tresd_hold_mle_grid[1] =  threshold_grid[j_best]  
        
        
        ######################################    
        false_alarms_KS = matrix(0,NMC,length(k_alpha_grid))
        false_alarms_KS_old = matrix(0,NMC,length(k_alpha_grid_old))
        false_alarms_Rnf = matrix(0,NMC,length(log_A_grid))
        delay_time_KS =   matrix(0,NMC,length(k_alpha_grid)) 
        delay_time_KS_old =   matrix(0,NMC,length(k_alpha_grid_old)) 
        delay_time_Rnf = matrix(0,NMC,length(log_A_grid))
        false_alarms_mle = matrix(0,NMC,length(tresd_hold_mle_grid))
        delay_time_mle =   matrix(0,NMC,length(tresd_hold_mle_grid)) 
        
        v_array = rep(0,NMC )
        
        for(iter in 1:NMC)
        { 
          print("iter")
          print(iter) 
          
          KS = rep(0,T) 
          KS_old = rep(0,T) 
          log_Rnf = rep(0,T) 
          mle = rep(0,T) 
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
            } #rpois(length(theta0),lambdac_t)}
            #  ycounts =  gridCounts(x[t,],d)
            
            cum_x = cum_x + x[t,]
            
            KS_old[t] =  sqrt(sum(cum_x))*kolmogorov_statistic(cum_x,F0)
            
            #
            log_Rnf[t] = logRnf(t,1,1,x,L,lambda0)
            
            mle[t] = mle_detection(x,t,L,lambda0)
            
            KS[t] = sequential_KS_statistics(x,t,F0,L)
            
            
            
          }## close for 1:T      
          proc.time() - ptm
          
          for(j in 1:length(k_alpha_grid))
          {
            false_alarms_KS[iter,j] = length(which( KS[1:v]> k_alpha_grid[j]))
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
        for(j in 1:length( tresd_hold_mle_grid ))
        {
          average_delay_time_mle[ind_T,ind_d,ind_TRF,j,ind_L] = mean( delay_time_mle[,j])
          average_false_alarms_mle[ind_T,ind_d,ind_TRF,j,ind_L] = mean(false_alarms_mle[,j])
        }
        for(j in 1:length(log_A_grid))
        {
          average_delay_time_Rnf[ind_T,ind_d,ind_TRF,j,ind_L] = mean( delay_time_Rnf[,j])
          average_false_alarms_Rnf[ind_T,ind_d,ind_TRF,j,ind_L] = mean(false_alarms_Rnf[,j])
          #         prob_detect_T_Rnf[ind_T,ind_d,ind_TRF,j] = length(indices)/NMC
        }
        
        
        print("L")
        print(L)
        print("average v")
        print(mean(v_array))
        print("average_delay_time_KS")
        print(average_delay_time_KS[ind_T,ind_d,ind_TRF,,ind_L])
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
        
      }## close for Total rate
    }## close for d
  }## close for T
  
}### close for L





matplot(cbind(average_delay_time_KS[ind_T,ind_d,ind_TRF,1,] + 1,average_delay_time_Rnf[ind_T,ind_d,ind_TRF,1,]+ 1,average_delay_time_mle[ind_T,ind_d,ind_TRF,1,]+ 1),type = "l" ,ylab ="average delay time")


aa = average_delay_time_KS[ind_T,ind_d,ind_TRF,1,]
bb = rep(mean(average_delay_time_KS_old[ind_T,ind_d,ind_TRF,1,]),length(L_grid))
cc = average_delay_time_Rnf[ind_T,ind_d,ind_TRF,1,]
dd = average_delay_time_mle[ind_T,ind_d,ind_TRF,1,]

matplot(cbind(aa,bb,cc,dd),type= "l")



par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
matplot(cbind(aa+1,bb+1,cc+1,dd+1),type="l",xlab="L",ylab = "Delay time",main ="",xaxt = "n",cex.lab = 1.4,cex.main= 1.4,col=c(1,2,3,4,6))
axis(1, at = c(4,8,12,17 ), labels = L_grid[ c(4,8,12,17 )])# cex.axis = 0.7
legend("topright", inset=c(-0.3,0),cex = 1.4,legend=c("KS ", "PKS","EF","GLR"), pch=c(1,3),  col=c(1,2,3,4,6))
