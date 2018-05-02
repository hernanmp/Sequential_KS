
rm(list = ls()) 

###   Code to replicate the results in Table 3

##  working directory
working_directory = "/Users/user/Desktop/JASA_code_Anomaly/Experiments"
setwd(paste(working_directory,"/Code_J",sep=""))



##  break points of the different testing sets, all the testing data is one file
init_grid =  c(1,186,333,410,484,561,761,870)
end_grid =  c(185,332,409,483,560,760,869,996)

##  alarm thresholds for different methods, to be changed later
KS_detection_time = 0
PKS_detection_time = 0
mle_detection_time = 0
glr_detection_time = 0
D2_detection_time = 0
cusum_detection_time = 0

library(lubridate)


setwd(paste(working_directory,"/Data_J",sep=""))

checksource_test_hist =  as.matrix(read.table("testing_data.csv"))

x = as.matrix(read.table("trainingS2.txt"))
colnames(x) = NULL
rownames(x) = NULL
training_break_points = as.matrix(read.table("training_break_points.txt"))
colnames(training_break_points ) = NULL
rownames(training_break_points ) = NULL
drop(training_break_points )

num_tests = length(training_break_points) +1
labels =  rep(0,dim(x)[1])
labels[1:training_break_points[1]] = 1 

for(i in 1: (length(training_break_points) -1))
{
  labels[(training_break_points[i]+1):(training_break_points[i+1])] = i+1   
}
labels[(training_break_points[i+1]+1):dim(x)[1]] = i+2
ind =  1:num_tests
size_labels = rep(0,num_tests)
for(i in 1:num_tests)
{
  size_labels[i] = length(which( labels==ind[i]))
}


checksource_train_hist =   as.matrix(read.table("trainingS1.txt"))
colnames(checksource_train_hist) = NULL
rownames(checksource_train_hist) = NULL
lambda0 = colMeans(checksource_train_hist)
  #as.matrix(read.table("lambda0.txt"))

colnames(lambda0) = NULL
rownames(lambda0) = NULL
lambda0 = drop(lambda0)

setwd(paste(working_directory,"/Code_J",sep=""))
source("Generate_data.R")
theta0 =  lambda0/sum(lambda0)
F0 = empirical_cdf(lambda0)

num_training = 1000  
alpha = .2
L =  20 

###########################################################################################3
###########################################################################################3
##
###### D2 
lambda0_new =  rep(0,285)  # downsampling
tlambda_new =  rep(0,285)
x2 = donwsample(x)
lambda0_new = as.vector(donwsample(t(as.matrix(lambda0) )))+10^-6
D2 = matrix(0,num_tests,max(size_labels))
# scr_tresdhold  = 0
lamb =  0.008
for(i in 1:num_tests)
{
  indices = which(labels==ind[i])
  
  new_x = x2[indices,]
  
  lambda_t = lambda0_new
  omega_t =  diag(lambda0_new)  
  
  for(t in 1:dim(new_x)[1])
  {
    x_temp = new_x[t,]
    
    temp = scr_update(lambda_t,omega_t,lambda0_new,x_temp,lamb)
    
    D2[i,t] =  temp$D2
    lambda_t =  temp$lambda_t
    omega_t =   temp$omega_t
  }
}

lower_D2 = min( setdiff(as.vector(D2),0))
upper_D2 = max( setdiff(as.vector(D2),0))
lower_D2 
# 
threshold_grid_D2  =  seq(lower_D2,upper_D2,length = 1000)

Expected_false_alarms_D2 = rep(0,length(threshold_grid_D2 ))

for(j in 1:length(threshold_grid_D2))
{
  aux_D2 = rep(0,num_tests)
  
  for(i in 1:num_tests)
  {
    aux_D2[i] =    length(which(D2[i,1:size_labels[ind[i]]]> threshold_grid_D2[j]))
  }
  Expected_false_alarms_D2[j] =  mean(aux_D2)   
}

j =  which.min(abs(Expected_false_alarms_D2 - alpha))
tresd_hold_D2  =  threshold_grid_D2[j] 
tresd_hold_D2

# 
# ###  choosing threshold for other methods
######### KS, PKS, mle, glr
KS = matrix(0,num_tests, max(size_labels))
mle = matrix(0,num_tests, max(size_labels))
glr = matrix(0,num_tests, max(size_labels))
PKS = matrix(0,num_tests, max(size_labels))
cusum = matrix(0,num_tests, max(size_labels))

###########
d_prob  = c(colSums(checksource_train_hist),0)+1
d_prob = d_prob/sum(d_prob)

for(j in 1:num_tests)
{ 
  indices = which(labels==ind[j])
  
  new_x = x[indices,]
  
  cum_x = rep(0,length(lambda0))
  
  for(t in 1:dim(new_x)[1])
  { 
    cum_x =  cum_x + new_x[t,]
    
    KS[j,t] = sequential_KS_statistics(new_x,t,F0,L)
    mle[j,t] = mle_detection(new_x,t,L,lambda0)
    glr[j,t]  =  log(logRnf(t,1,1,new_x,L,lambda0  ))
    PKS[j,t] = sqrt(sum(cum_x))*kolmogorov_statistic(cum_x,F0)
    
    if(t==1)
    {  
      S1_prev =  rep(0,length(lambda0)+1)
      S2_prev =  rep(0,length(lambda0)+1)
    }## #
    temp  = cusum_statistic(c(new_x[t,]-lambda0,0),d_prob,S1_prev,S2_prev,k=0)
    cusum[j,t] =  temp$y            
    S1_prev = temp$S1_now
    S2_prev = temp$S2_now
  }######
}#####


lower_KS = min( setdiff(as.vector(KS),0))
upper_KS = max( setdiff(as.vector(KS),0))
lower_mle = min( setdiff(as.vector(mle),0))
upper_mle = max( setdiff(as.vector(mle),0))
lower_glr = min( setdiff(as.vector(glr),0))
upper_glr = max( setdiff(as.vector(glr),0))
lower_PKS = min( setdiff(as.vector(PKS),0))
upper_PKS = max( setdiff(as.vector(PKS),0))
lower_cusum = min( setdiff(as.vector(cusum),0))
upper_cusum = max( setdiff(as.vector(cusum),0))
# 
threshold_grid_KS  =  seq(lower_KS,upper_KS,length = 1000)
threshold_grid_mle  =  seq(lower_mle,upper_mle,length = 1000)
threshold_grid_glr  =  seq(lower_glr,upper_glr,length = 1000)
threshold_grid_PKS  =  seq(lower_PKS,upper_PKS,length = 1000)
threshold_grid_cusum  =  seq(lower_cusum,upper_cusum,length = 1000)

Expected_false_alarms_KS = rep(0,length(threshold_grid_KS ))
Expected_false_alarms_mle = rep(0,length(threshold_grid_mle ))
Expected_false_alarms_glr = rep(0,length(threshold_grid_glr ))
Expected_false_alarms_PKS = rep(0,length(threshold_grid_PKS ))
Expected_false_alarms_cusum = rep(0,length(threshold_grid_cusum ))

for(j in 1:length(threshold_grid_KS))
{
  aux_KS = rep(0,num_tests)
  aux_mle = rep(0,num_tests)
  aux_glr = rep(0,num_tests)
  aux_PKS = rep(0,num_tests)
  aux_cusum = rep(0,num_tests)
  
  for(i in 1:num_tests)
  {
    aux_KS[i]  =    length(which(KS[i,1:size_labels[ind[i]]]> threshold_grid_KS[j]))
    aux_mle[i] =    length(which(mle[i,1:size_labels[ind[i]]]> threshold_grid_mle[j]))
    aux_glr[i] =    length(which(glr[i,1:size_labels[ind[i]]]> threshold_grid_glr[j]))
    aux_PKS[i] =    length(which(PKS[i,1:size_labels[ind[i]]]> threshold_grid_PKS[j]))
    aux_cusum[i] =    length(which(cusum[i,1:size_labels[ind[i]]]> threshold_grid_cusum[j]))
  }
  Expected_false_alarms_KS[j] =  mean(aux_KS)   
  Expected_false_alarms_mle[j] =  mean(aux_mle)  
  Expected_false_alarms_glr[j] =  mean(aux_glr)  
  Expected_false_alarms_PKS[j] =  mean(aux_PKS)
  Expected_false_alarms_cusum[j] =  mean(aux_cusum)  
}

j =  which.min(abs(Expected_false_alarms_KS - alpha))
threshold_grid_KS[j] 
##  can change this to use the therotetical rule, KS^* in the paper
k_alpha_grid = threshold_grid_KS[j] #sqrt(log(2*L*min(size_labels)*1/alpha)/2)
#threshold_grid_KS[j] 
#sqrt(log(2*L*min(size_labels)*1/alpha)/2)
j =  which.min(abs(Expected_false_alarms_mle - alpha))
tresd_hold_mle  =  threshold_grid_mle[j] 
tresd_hold_mle

j =  which.min(abs(Expected_false_alarms_glr - alpha))
tresd_hold_glr  =  threshold_grid_glr[j] 
tresd_hold_glr

j =  which.min(abs(Expected_false_alarms_PKS - alpha))
tresd_hold_PKS  =  threshold_grid_PKS[j] 
tresd_hold_PKS

j =  which.min(abs(Expected_false_alarms_cusum - alpha))
tresd_hold_cusum  =  threshold_grid_cusum[j] 
tresd_hold_cusum


#############################################################################3
#############################################################################3
#############################################################################3
#############################################################################3
#############################################################################3
#############################################################################3
### Testing

training_set_dat = checksource_train_hist
for(j in 1:num_tests)
{
  indices = which(labels==ind[j])
  
  training_set_dat = rbind(training_set_dat, x[indices,])
}

lambda0 = colMeans(training_set_dat)
lambda0_new = as.vector(donwsample(t(as.matrix(lambda0) ))) + 10^-6

d_prob  =   c(colSums(training_set_dat),0)+1  ### posterior mean after dirichlet uniform prior
d_prob =  d_prob/sum(d_prob)

##################

for(test_ind in 1:length(init_grid))
{
  
  init = init_grid[test_ind]
  end  = end_grid[test_ind]
  
  test_set = checksource_test_hist[init:end,]
  test_set2  = donwsample(test_set)
  KS2 =  rep(0,dim(test_set)[1])
  mle2 = rep(0,dim(test_set)[1]) 
  glr2 = rep(0,dim(test_set)[1]) 
  PKS2 =  rep(0,dim(test_set)[1])
  D22 = rep(0,dim(test_set)[1])
  cusum2 = rep(0,dim(test_set)[1])
  
  cum_x = rep(0,length(lambda0))
  
  lambda_t = lambda0_new
  omega_t =  diag(lambda0_new )  
  
  for(t in 1:dim(test_set)[1])
  {
    
    #print(t) 
    
    cum_x =   cum_x + test_set[t,]
    
    KS2[t] = sequential_KS_statistics(test_set,t,F0,L)
    mle2[t] = mle_detection(test_set,t,L,lambda0)
    glr2[t]  = log( logRnf(t,1,1,test_set,L,lambda0  ))
    PKS2[t] = sqrt(sum(cum_x))*kolmogorov_statistic(cum_x,F0)
    
    x_temp = test_set2[t,]
    
    temp = scr_update(lambda_t,omega_t,lambda0_new,x_temp,lamb)
    
    D22[t] =  temp$D2
    lambda_t =  temp$lambda_t
    omega_t =   temp$omega_t
    
    if(t==1)
    {
      S1_prev =  rep(0,length(lambda0)+1)
      S2_prev =  rep(0,length(lambda0)+1)
    }###
    temp  = cusum_statistic(c(test_set[t,]-lambda0,0),d_prob,S1_prev,S2_prev,k=0)
    cusum2[t] =  temp$y            
    S1_prev = temp$S1_now
    
    
    S2_prev = temp$S2_now
    
  }## close for 1:T 
  
  init = 1
  end = dim(test_set)[1]
  matplot(cbind(KS2[init:end],rep(k_alpha_grid,length(KS2[init:end]))),type="l")
  matplot(cbind(mle2[init:end],rep(tresd_hold_mle,length(mle2[init:end]))),type="l")
  matplot(cbind(glr2[init:end],rep(tresd_hold_glr,length(glr2[init:end]))),type="l")
  matplot(cbind(D22[init:end],rep(tresd_hold_D2,length(D22[init:end]))),type="l")
  matplot(cbind(cusum2[init:end],rep(tresd_hold_cusum,length(cusum2[init:end]))),type="l")
  
  #   
  KS_detection_time    =  min(which(KS2>k_alpha_grid))
  mle_detection_time =   min(which(mle2>tresd_hold_mle))
  glr_detection_time =  min(which(glr2>tresd_hold_glr))
  PKS_detection_time =  min(which(PKS2>tresd_hold_PKS))
  D2_detection_time =  min(which(D22>tresd_hold_D2))
  cusum_detection_time =  min(which(cusum2>tresd_hold_cusum))
  
  
  print("init")
  print(init_grid[test_ind])
  print("KS")
  print(KS_detection_time)
  print("mle")
  print(mle_detection_time)
  print("glr")
  print(glr_detection_time)
  print("pks")
  print(PKS_detection_time )
  print("D2")
  print(D2_detection_time)
  print("cusum")
  print(cusum_detection_time)
  #  
}



