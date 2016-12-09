
rm(list = ls()) 

###   Code to replicate the results in Table 4

##  working directory
working_directory = "/Users/user/Desktop/JASA_code_Anomaly/Experiments"
setwd(paste(working_directory,"/Code",sep=""))

## Window size
L =  20 


##  break points of the different testing sets, all the testing data is one file
init_grid =  c(1,186,333,410,484,561,761,870)
end_grid =  c(185,332,409,483,560,760,869,996)

##  alarm thresholds for different methods, to be changed later
KS_detection_time = 0
PKS_detection_time = 0
mle_detection_time = 0
glr_detection_time = 0
D2_detection_time = 0

library(lubridate)


setwd(paste(working_directory,"/Data",sep=""))


checksource_train_hist = as.matrix(read.table("training_data.csv"))
checksource_test_hist =  as.matrix(read.table("testing_data.csv"))

dim(checksource_train_hist)




setwd(paste(working_directory,"/Code",sep=""))
source("Generate_data.R")


### We split training data in two parts, the first half to estimate the mean and the second one to 
###  calibrate the methods
num_training = 1000  
alpha = .2

###  mean of the first half of the training data
lambda0 =  colMeans(checksource_train_hist[1:num_training,])

## exploring if there is overdispersion
training_set = checksource_train_hist[1:num_training,]
mean_train = colMeans(training_set)
sd_train = apply(training_set,2,sd)
temp = which(mean_train > 0 )
mean(sd_train[temp]^2/mean_train[temp])


###  estimate background on training set
theta0 =  lambda0/sum(lambda0)
F0 = empirical_cdf(lambda0)

### held out set (testing set is different)
x = checksource_train_hist[(1+num_training):dim(checksource_train_hist)[1],]

##  split held out on a small number set to see if this  follow the same distribution from the trainig set
num_tests = 10

labels= cut(1:dim(x)[1],num_tests,labels = 1:num_tests)
size_labels = rep(0,num_tests)

for(jj in 1:num_tests)
{
  indices = which(labels==jj)
  size_labels[jj] = length(indices)
}

##  compute test statistics on held out sets
two_sample_test = rep(0,num_tests)
for(jj in 1:num_tests)
{
  indices = which(labels==jj)
  new_x = x[indices,]
  temp = colSums(checksource_train_hist[1:num_training,])
  n1 = sum(temp)
  n2 = sum(new_x)
  two_sample_test[jj] = max(abs(empirical_cdf(temp) - empirical_cdf(colSums(new_x)) ))*sqrt(   (n1*n2)/(n1+n2)) 
}
## compute p_values and adjusted p_values
p_values = 1-kolmogorov_distribution_cdf(two_sample_test)
adjusted_p_values = p.adjust(p_values, method = "BH", n = length(p_values)) 
ind = which(adjusted_p_values>.1)     
length(ind)

##disregard from checksource_train_hist   those subsets found to fail the KS two samples test
checksource_train_hist  = checksource_train_hist[1:num_training,]
for(jj in 1:length(ind))
{
  aux = x[which(labels ==ind[jj]),]
  checksource_train_hist = rbind(checksource_train_hist,aux)
}



#################################################################################3333
###  Constructing new held out sets
###  it is important that we avoid outliers as with small data set like this that would 
###   dramatically affect the selection of thresholds
x= checksource_train_hist[(1+num_training):dim(checksource_train_hist)[1],]

##  spliting new held out set in more subsets to estimate alarm thresdhold for different methods 
num_tests = 20#

labels= cut(1:dim(x)[1],num_tests,labels = 1:num_tests)
size_labels = rep(0,num_tests)

for(jj in 1:num_tests)
{
  indices = which(labels==jj)
  size_labels[jj] = length(indices)
}

two_sample_test = rep(0,num_tests)
for(jj in 1:num_tests)
{
  indices = which(labels==jj)
  new_x = x[indices,]
  temp = colSums(checksource_train_hist[1:num_training,])
  n1 = sum(temp)
  n2 = sum(new_x)
  two_sample_test[jj] = max(abs(empirical_cdf(temp) - empirical_cdf(colSums(new_x)) ))*sqrt(   (n1*n2)/(n1+n2)) 
}
# compute p_values
p_values = 1-kolmogorov_distribution_cdf(two_sample_test)
adjusted_p_values = p.adjust(p_values, method = "BH", n = length(p_values)) 
ind = which(adjusted_p_values>.1)
length(ind)
num_tests = length(ind)

###########################################################################################3
###########################################################################################3
##  Compute thresd hold but only using those held out sets that passed the test.

###### D2 
lambda0_new =  rep(0,285)  # downsampling
tlambda_new =  rep(0,285)
x2 = donwsample(x)
lambda0_new = as.vector(donwsample(t(as.matrix(lambda0) )))+10^-6



D2 = matrix(0,num_tests,max(size_labels))
# scr_tresdhold  = 0
lamb =  0.008
for(ii in 1:num_tests)
{
  indices = which(labels==ind[ii])
  
  new_x = x2[indices,]
  
  lambda_t = lambda0_new
  omega_t =  diag(lambda0_new)  
  
  for(t in 1:dim(new_x)[1])
  {
    xx = new_x[t,]
    
    temp = scr_update(lambda_t,omega_t,lambda0_new,xx,lamb)
    
    D2[ii,t] =  temp$D2
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

for(jj in 1:length(threshold_grid_D2))
{
  aux_D2 = rep(0,num_tests)
  
  for(ii in 1:num_tests)
  {
    aux_D2[ii] =    length(which(D2[ii,1:size_labels[ind[ii]]]> threshold_grid_D2[jj]))
  }
  Expected_false_alarms_D2[jj] =  mean(aux_D2)   
}

jjj =  which.min(abs(Expected_false_alarms_D2 - alpha))
tresd_hold_D2  =  threshold_grid_D2[jjj] 
tresd_hold_D2

# 
# ###  choosing threshold for other methods
######### KS, PKS, mle, glr
KS = matrix(0,num_tests, max(size_labels))
mle = matrix(0,num_tests, max(size_labels))
glr = matrix(0,num_tests, max(size_labels))
PKS = matrix(0,num_tests, max(size_labels))


for(jj in 1:num_tests)
{
  indices = which(labels==ind[jj])
  
  new_x = x[indices,]
  
  cum_x = rep(0,length(lambda0))
  
  for(t in 1:dim(new_x)[1])
  {
    cum_x =  cum_x + new_x[t,]
    
    KS[jj,t] = sequential_KS_statistics(new_x,t,F0,L)
    mle[jj,t] = mle_detection(new_x,t,L,lambda0)
    glr[jj,t]  =  log(logRnf(t,1,1,new_x,L,lambda0  ))
    PKS[jj,t] = sqrt(sum(cum_x))*kolmogorov_statistic(cum_x,F0)
  }
}


lower_KS = min( setdiff(as.vector(KS),0))
upper_KS = max( setdiff(as.vector(KS),0))
lower_mle = min( setdiff(as.vector(mle),0))
upper_mle = max( setdiff(as.vector(mle),0))
lower_glr = min( setdiff(as.vector(glr),0))
upper_glr = max( setdiff(as.vector(glr),0))
lower_PKS = min( setdiff(as.vector(PKS),0))
upper_PKS = max( setdiff(as.vector(PKS),0))
# 
threshold_grid_KS  =  seq(lower_KS,upper_KS,length = 1000)
threshold_grid_mle  =  seq(lower_mle,upper_mle,length = 1000)
threshold_grid_glr  =  seq(lower_glr,upper_glr,length = 1000)
threshold_grid_PKS  =  seq(lower_PKS,upper_PKS,length = 1000)


Expected_false_alarms_KS = rep(0,length(threshold_grid_KS ))
Expected_false_alarms_mle = rep(0,length(threshold_grid_mle ))
Expected_false_alarms_glr = rep(0,length(threshold_grid_glr ))
Expected_false_alarms_PKS = rep(0,length(threshold_grid_PKS ))



for(jj in 1:length(threshold_grid_KS))
{
  aux_KS = rep(0,num_tests)
  aux_mle = rep(0,num_tests)
  aux_glr = rep(0,num_tests)
  aux_PKS = rep(0,num_tests)
  
  for(ii in 1:num_tests)
  {
    aux_KS[ii]  =    length(which(KS[ii,1:size_labels[ind[ii]]]> threshold_grid_KS[jj]))
    aux_mle[ii] =    length(which(mle[ii,1:size_labels[ind[ii]]]> threshold_grid_mle[jj]))
    aux_glr[ii] =    length(which(glr[ii,1:size_labels[ind[ii]]]> threshold_grid_glr[jj]))
    aux_PKS[ii] =    length(which(PKS[ii,1:size_labels[ind[ii]]]> threshold_grid_PKS[jj]))
  }
  Expected_false_alarms_KS[jj] =  mean(aux_KS)   
  Expected_false_alarms_mle[jj] =  mean(aux_mle)  
  Expected_false_alarms_glr[jj] =  mean(aux_glr)  
  Expected_false_alarms_PKS[jj] =  mean(aux_PKS)  
}

jjj =  which.min(abs(Expected_false_alarms_KS - alpha))
threshold_grid_KS[jjj] 
##  can change this to use the therotetical rule, KS^* in the paper
k_alpha_grid = threshold_grid_KS[jjj] #sqrt(log(2*L*min(size_labels)*1/alpha)/2)
#threshold_grid_KS[jjj] 
#sqrt(log(2*L*min(size_labels)*1/alpha)/2)




jjj =  which.min(abs(Expected_false_alarms_mle - alpha))
tresd_hold_mle  =  threshold_grid_mle[jjj] 
tresd_hold_mle


jjj =  which.min(abs(Expected_false_alarms_glr - alpha))
tresd_hold_glr  =  threshold_grid_glr[jjj] 
tresd_hold_glr

jjj =  which.min(abs(Expected_false_alarms_PKS - alpha))
tresd_hold_PKS  =  threshold_grid_PKS[jjj] 
tresd_hold_PKS




#############################################################################3
#############################################################################3
#############################################################################3
#############################################################################3
#############################################################################3
#############################################################################3
### Testing

training_set_dat = checksource_train_hist[1:num_training,]
for(jj in 1:num_tests)
{
  indices = which(labels==ind[jj])
  
  training_set_dat = rbind(training_set_dat, x[indices,])
}

lambda0 = colMeans(training_set_dat)
lambda0_new = as.vector(donwsample(t(as.matrix(lambda0) ))) + 10^-6


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
    
    xx = test_set2[t,]
    
    temp = scr_update(lambda_t,omega_t,lambda0_new,xx,lamb)
    
    D22[t] =  temp$D2
    lambda_t =  temp$lambda_t
    omega_t =   temp$omega_t
    
  }## close for 1:T 
  
  init = 1
  end = dim(test_set)[1]
  matplot(cbind(KS2[init:end],rep(k_alpha_grid,length(KS2[init:end]))),type="l")
  matplot(cbind(mle2[init:end],rep(tresd_hold_mle,length(mle2[init:end]))),type="l")
  matplot(cbind(glr2[init:end],rep(tresd_hold_glr,length(glr2[init:end]))),type="l")
  matplot(cbind(D22[init:end],rep(tresd_hold_D2,length(D22[init:end]))),type="l")
  
  #   
  KS_detection_time    =  min(which(KS2>k_alpha_grid))
  mle_detection_time =   min(which(mle2>tresd_hold_mle))
  glr_detection_time =  min(which(glr2>tresd_hold_glr))
  PKS_detection_time =  min(which(PKS2>tresd_hold_PKS))
  D2_detection_time =  min(which(D22>tresd_hold_D2))
  
  
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
  #  
}


