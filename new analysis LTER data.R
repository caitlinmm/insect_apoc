# Here, I will try to take the data from Crossley and compare what the violin plot
# would look like without standardizing away the time trend/other odd pre-processing
# steps that seemed to distort the data unnecessarily.
#####
# Here is the file they were on before they started all sorts of weird distorting processing
data1 = read.csv('PerSpecies_Abundance_LTER_annotated.csv',as.is=T,check.names=F,header=T)

# I am going to try to sort out species and then for each species plot a lm and find the 
# slope of the regression line. Seems like at this point Species.code is where the most 
# accurate information about species is rather than the Species column
u.species = unique(data1$Species.code)
# This is an empty vector that we will fill with the slopes of the regressions of each species
slopes = rep(NA, length(u.species))

for(u in 1:length(u.species))
{
  sp_name = u.species[u]
  sp_ind = which(data1$Species.code == sp_name)
  sp_data = data1[sp_ind, ]
  
  # The Konza prairie subset of the data all meets the minimum filtering requirements (at least
  # 3 measurements per species/site and at least 1 non-zero measurement) so we can ignore the filtering
  # requirements.
  
  # Now that we have the indices of just our species of interest, we can call a linear
  # regression from that and evaluate the slope of it.
  slopes[u] = as.numeric(lm(Abundance ~ Year, data = data1[sp_ind, ])[[1]][2])
}
#####

# Trying to isolate the effect of the AR_reml() function here by doing all the other stretching they do but not that.
#####
LLS = paste(data1$LTER.site,data1$Locale,data1$Species.code,sep='_') #unique LTER by sub-site by species combinations
u.LLS = unique(LLS)
LLS.count = apply(array(u.LLS),1,function(x){d2=data1[which(LLS==x),];length(which(d2$Abundance!=0))})
u.LLS2 = u.LLS[which(LLS.count>0)]
time.trends = data.frame('LTER'=NA,'Site'=NA,'Species'=NA, 'coef_slope'=-999)
tx = 1
# this goes through each species in each site
for (l in 1:length(u.LLS2)){
  print(l)
  lls.data = data1[which(LLS==u.LLS2[l]),]
  ## This part I added in to try and keep the data to just the 2015 data used in Crossley.
  # lls.data = lls.data[which(lls.data$Year <= 2015), ]
  ## It seems like this is just recreating the site and study name? Not sure why
  LL = paste0(strsplit(u.LLS2[l],'_')[[1]][1:2],collapse='_') #trim "no sample" years
  ## Removed a harvard forest ants specific processing step
  ### Quality threshold
  ## Checking if there are more than 3 years where there are actual measurements for abundance and that of those non-NA
  ## years there are at least a few with non-zero values for abundance measured.
  if ( (length(which(!is.na(lls.data$Abundance))) > 3) & (length(which(lls.data$Abundance[!is.na(lls.data$Abundance)] > 0)) > 0) ){
    ### 
    ## This is measuring the non-zero minimum
    cmin = min(lls.data$Abundance[which(lls.data$Abundance>0)],na.rm=T)
    ## I don't understand how cmin would be NA after their preprocessing steps? Should be at least 1 non-zero values?
    ## I suppose this is their way of transforming so that they can take the log of the data
    if (is.na(cmin)){
      lls.data$Abundance[which(lls.data$Abundance==0)] = 0.5 
    } else {
      lls.data$Abundance[which(lls.data$Abundance==0)] = 0.5 * cmin #replace zeroes with 0.5*minimum abundance value in this time series
    }
    lls.data$Abundance = log(lls.data$Abundance) #log-transform abundances
    if (length(which(lls.data$Abundance==0))==length(lls.data$Abundance)){ #abundances are all = 1
      slope_plc = 0
    } else {
      lls.years = sort(lls.data$Year)
      lls.years.stretch = seq(lls.years[1],lls.years[length(lls.years)],1) #expand time series to include missing years
      llsy.match = match(lls.years.stretch,lls.years)
      X = llsy.match
      ## Here they are replacing the years with abundances and then having NAs for the years missing data in the expanded set
      X[which(!is.na(X))] = lls.data$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
      ## For some reason R is not allowing me to keep the format the same as they had it, so instead of X - mean() I am 
      ## trying -mean() + X
      Z = (-mean(X, na.rm=T) + X)/sd(X, na.rm=T) # z-transform
      t.scale = 1:length(lls.years.stretch)
      ## Here they are standardizing timescales, I assume to make so that it is the same across different lengths of monitoring.
      ## It says it is transforming it between 0 and 1 but I think that it would need to be divided by max(t.scale) - min(t.scale)
      # for that to be the case. Instead its transforming to ~almost~ 1.
      t.scale = (t.scale-min(t.scale))/max(t.scale) #original transform: scale between 0 and 1
      slope_plc = as.numeric(lm(Z ~ t.scale)[[1]][2])
    }
    time.trends[tx,1] = strsplit(u.LLS2[l],'_')[[1]][1] #LTER
    time.trends[tx,2] = strsplit(u.LLS2[l],'_')[[1]][2] #Site
    time.trends[tx,3] = paste0(strsplit(u.LLS2[l],'_')[[1]][3:length(strsplit(u.LLS2[l],'_')[[1]])],collapse='_') #Species
    time.trends[tx,4] =slope_plc
    tx = tx + 1
    # Plot high quality time series
    png(paste0(gsub('\\?','',u.LLS[l]),'.png'))
    plot(t.scale,Z,type='l',main=u.LLS2[l],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2)
    abline(a=add.coef.Z[3],b=add.coef.Z[4],lty=2,col='red',lwd=1.5)
    legend('topright',legend=c(paste0('slope = ',round(add.coef.Z[4],2)),paste0('length = ',length(which(!is.na(Z)))),paste0('autocor = ',round(add.coef.Z[2],2)),paste0('logLik = ',round(add.coef.Z[7],2))),ncol=1)
    dev.off()
  } else {
    #ignore low quality time series
  }
}



#####



# Evaluating ours with their pre-processing steps vs. without
######
# In my opinion, key things to bookmark from how they did it (that are sketchy) are that they log transformed everything even though
# we expect the trends to be basically linear without transformation, and they normalized the data over the standard deviation of time
# steps. They also then additionally used the AR_reml() function which is made to deal with autocorrelation, but they had already taken care
# of any autocorrelation in the data set in early processing steps.
m_dat = read.csv(file = "/Users/caitlinmiller/Desktop/Github/InsectTrends/data/Model2_sampledAbundance_RAW.csv", header = FALSE)

# They log transform all files so first we will do the same.
m_dat_l = log(m_dat)
slopes = rep(NA, ncol(m_dat_l))

# Next we will go through each column and transform values based on the mean and standard deviation of the species, and then use AR_reml()
# to derive a value for slope.
for(sp in 1:ncol(m_dat_l))
{
  X = m_dat_l[,sp]
  Z = (X - mean(X))/sd(X)
  t_scale = seq(from = 0, to = 1, length.out = nrow(m_dat_l))
  slopes[sp] = AR_reml(Z ~ t_scale)[[3]][2,1]
}

all_dat = data.frame(LTER = rep('mod', length(slopes)), slopes)


# Now trying to do an analysis where we just take a linear model of each column and then take the slope of the regression as the slope of change
slopes = rep(NA, ncol(m_dat))
times = 1:nrow(m_dat)

for(sp in 1:ncol(m_dat))
{
  abuns = m_dat[,sp]
  slopes[sp] = as.numeric(lm(abuns ~ times)[[1]][2])
}

all_dat = rbind(all_dat, cbind(LTER = rep('norm', length(slopes)), slopes))


# Now something that only log transforms data because our data does drop so suddenly it is exponential.
slopes = rep(NA, ncol(m_dat_l))
times = 1:nrow(m_dat_l)

for(sp in 1:ncol(m_dat_l))
{
  abuns = m_dat_l[,sp]
  slopes[sp] = as.numeric(lm(abuns ~ times)[[1]][2])
}

all_dat = rbind(all_dat, cbind(LTER = rep("log", length(slopes)), slopes))

all_dat$slopes = as.numeric(all_dat$slopes)

ggplot(data = all_dat, aes(LTER, slopes)) + geom_violin()
#####

# Calculating abundance over time in our model
#####
# For each row (which is equivalent to a time step in our model), we want to take the sum to evaluate abundance over
# time as a "biomass" metric
m_dat = read.csv(file = "/Users/caitlinmiller/Desktop/Github/InsectTrends/data/Model2_sampledAbundance_RAW.csv", header = FALSE)
abuns = rep(NA, nrow(m_dat))

for(r in 1:nrow(m_dat))
{
  abuns[r] = sum(m_dat[r, ])
}



#####


# Looking at abundance vs. slope to see if there is any sort of relationship
#####
# here, I am trying to create a data frame that brings together slopes, starting abundances, and 
# max abundances to look at the correlation between these measures in the real data vs. our simulated 
# data.
data1 = read.csv('PerSpecies_Abundance_LTER_annotated.csv',as.is=T,check.names=F,header=T)

u.species = unique(data1$Species.code)
no_species = length(u.species)

slopes = rep(NA, no_species)
max_abun = rep(NA, no_species)
start_abun = rep(NA, no_species)

for(spp in 1:no_species)
{
  spp_i = which(data1$Species.code == u.species[spp])
  slopes[spp] = as.numeric(lm(Abundance ~ Year, data = data1[spp_i, ])[[1]][2])
  max_abun[spp] = max(data1$Abundance[spp_i])
  # For start abundance, I will take an average of the first 5 measurements to account for the possibility
  # that the first measurement is an outlier.
  start_abun[spp] = mean(data1$Abundance[spp_i][1:5])
}










#####

