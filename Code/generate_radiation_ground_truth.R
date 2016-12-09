






rm(list = ls())

### Code to generate an example  of normal background and anomlous densities
###  this can be used to replicate the results in Table 3
###  It requires to select the size of the type of radition source (celsium or cobalt),
### and its distance to the spetometer
## once that is selected this code will produce the respective pre an post change densities
working_directory = "/Users/user/Desktop/JASA_code_Anomaly/Experiments"

setwd(paste(working_directory,"/Code",sep=""))


source('utils.R')

setwd(paste(working_directory,"/Data",sep=""))



# How strong is the background rate, and how large is the source in milliCuries?
background_cps = 39  
source_size = 650  #mCi


# Distance to the source
dist_grid = c(50, 100,  150)
i = 3## index of distance used, if for example i = 1, then distance to source is 50m


# Read in background
background_train = as.numeric(read.csv('background_train_mean.csv', header=FALSE))
plot(background_train, type='l')
n_bins = length(background_train)

# choose anomaly
celsium  = 0 # if celsium = 1 then the anomaly is celsium is celsium-137. If set to zero we use cobalt

# Read in anomaly


if(celsium==1)
{
  cs137_spectrum = as.numeric(read.csv('2013-cs137-5cm.csv', header=FALSE))
}
if(celsium ==0)
{
  cs137_spectrum = as.numeric(read.csv('co60-4min-320cps-c7-20cm.csv', header=FALSE))
}





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
