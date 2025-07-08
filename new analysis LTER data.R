# Here, I will try to take the data from Crossley and compare what the violin plot
# would look like without standardizing away the time trend/other odd pre-processing
# steps that seemed to distort the data unnecessarily.

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


# Evaluating ours with their pre-processing steps vs. without
######
# In my opinion, key things to bookmark from how they did it (that are sketchy) are that they log transformed everything even though
# we expect the trends to be basically linear without transformation, and they normalized the data over the standard deviation of time
# steps. They also then additionally used the AR_reml() function which is made to deal with autocorrelation, but they had already taken care
# of any autocorrelation in the data set in early processing steps.
m1_dat = read.csv(file = "/Users/caitlinmiller/Desktop/Github/InsectTrends/data/Model1_sampledAbundance_RAW.csv", header = FALSE)

# They log transform all files so first we will do the same.
m1_dat_l = log(m1_dat)
slopes = rep(NA, ncol(m1_dat_l))

# Next we will go through each column and transform values based on the mean and standard deviation of the species, and then use AR_reml()
# to derive a value for slope.
for(sp in 1:ncol(m1_dat_l))
{
  X = m1_dat_l[,sp]
  Z = (X - mean(X))/sd(X)
  t_scale = seq(from = 0, to = 1, length.out = nrow(m1_dat_l))
  slopes[sp] = AR_reml(Z ~ t_scale)[[3]][2,1]
}

all_dat = data.frame(LTER = rep('mod', length(slopes)), slopes)


# Now trying to do an analysis where we just take a linear model of each column and then take the slope of the regression as the slope of change
slopes = rep(NA, ncol(m1_dat))
times = 1:nrow(m1_dat)

for(sp in 1:ncol(m1_dat))
{
  abuns = m1_dat[,sp]
  slopes[sp] = as.numeric(lm(abuns ~ times)[[1]][2])
}

all_dat = rbind(all_dat, cbind(LTER = rep('norm', length(slopes)), slopes))


# Now something that only log transforms data because our data does drop so suddenly it is exponential.
slopes = rep(NA, ncol(m1_dat_l))
times = 1:nrow(m1_dat_l)

for(sp in 1:ncol(m1_dat_l))
{
  abuns = m1_dat_l[,sp]
  slopes[sp] = as.numeric(lm(abuns ~ times)[[1]][2])
}

all_dat = rbind(all_dat, cbind(LTER = rep("log", length(slopes)), slopes))

all_dat$slopes = as.numeric(all_dat$slopes)

ggplot(data = all.dat, aes(LTER, mod_slopes)) + geom_violin()
