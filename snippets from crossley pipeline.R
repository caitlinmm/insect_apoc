## These snippets from the Crossley et al. pipeline are just focused on processing Konza prairie information.
## Data from the gall insects isn't very interesting because there is just 1 single group so there is no distribution
## of slopes. An important caveat is that the current datasets include data up to 2020 but the paper only used data
## up until 2015 so more recent data needs to be filtered out for reproducibility. CSV names and locations may need to
## be edited depending on local information locations.

# Initial analysis of data sets: This section groups abundances of all observations from a single species in a single
# site in a single year.
#####
##### Konza Prairie gall insects
site.name = 'KonzaPrairie'
locale.name = 'gall'
data1 = read.csv("/Users/caitlinmiller/Desktop/Github/insect_apoc/CGP011.csv",as.is=T,check.names=F,header=T)
data1$Year = data1$RecYear
data1$Number = data1$GalledStems
u.years = sort(as.numeric(unique(data1$Year)))
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (y in 1:length(u.years)){
  y.data = data1[which(data1$Year==u.years[y]),]
  n.y1 = sum(y.data$SampledStems,na.rm=T)
  out.slopes[cx,1] = site.name #LTER.site
  out.slopes[cx,2] = locale.name #Locale
  out.slopes[cx,3] = 'gall insects' #Species.code
  out.slopes[cx,4] = u.years[y] #Year
  out.slopes[cx,5] = n.y1 #Number of observations
  out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
  cx = cx + 1
}
write.table(out.slopes, file = 'GallInsect_SpeciesSlopes.txt', sep='\t',quote=F,row.names=F)


##### Konza Prairie grasshoppers

# This analysis also had evaluations of diversity in the original Crossley pipeline (missing from gall insects 
# because only 1 subset was sampled) but we are not interested in that so I have pruned them out.

site.name = 'KonzaPrairie'
locale.name = 'grasshopper'
data1 = read.csv("/Users/caitlinmiller/Desktop/Github/insect_apoc/CGR022.csv",as.is=T,check.names=F,header=T)
data1$Species = data1$SPECIES
data1$Year = data1$RECYEAR
data1$Number = as.numeric(data1$TOTAL)
# Here, they screened out data where the species wasn't identified, including "unknown", that with a space after it, and a blank section
# for that entry. (only cut out ~100 entries)
data1 = data1[which(data1$Species!='unknown' & data1$Species!='unknown ' & data1$Species!=''),]
u.years = sort(as.numeric(unique(data1$Year)))
# I am not really sure why they have this line in here when they re-assign this object later with species names with proper spelling??
u.species = unique(data1$Species)
# Fix species names
# Merge duplicates that have different species name spelling
# This is a key that 
konza.duplicates = read.csv('/Users/caitlinmiller/Downloads/External_Database_S3/taxa_keys/KonzaPrairie_grasshoppers_duplicatenames.csv', 
                            as.is=T,check.names=F,header=F)
# Not sure why they are filtering for unknown again? Just because they left uppercase unknown out the first time?
data1 = data1[which(data1$Species!='Unknown'),]
# OK basically the taxa key matches up various alternative spellings of name with a standardized version of the same name.
# The first column is always the undesired potential alternative name, while the second column is the desired name. This function
# looks for which alternative (un-uniform) spelling is used from the first column, and replaces it with the standardized version from
# the second column. If there isn't a match in the first column, they instead look for a match in the second column which would indicate
# that it already matches the desired format. There are some columns that appear to have the proper name twice with no discernible different
# which seems redundant but probably has a reason? Also not positive what the else statement is doing here because it is empty...
# All species names are then converted to the standard.
for (i in 1:nrow(data1)){
  name.swap = konza.duplicates[which(konza.duplicates[,1]==data1$Species[i]),2]
  if (length(name.swap)==0){
    name.swap = konza.duplicates[match(data1$Species[i],konza.duplicates[,2]),2]
  } else {
  }
  data1$Species[i] = name.swap
  
}
# With the removal of duplicates/standardization, they reduce the number of species names from 128 to 55!
u.species = unique(data1$Species)


#Get values for species-centric table
# I believe this part here is the actually relevant bit for us!! Looking at slopes which I believe is what was used in fig. 2
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
  for (y in 1:length(u.years)){
    y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
    n.y1 = length(unique(paste(y.data$RECMONTH,y.data$RECDAY,y.data$REPSITE,sep='_')))
    out.slopes[cx,1] = site.name #LTER.site
    out.slopes[cx,2] = locale.name #Locale
    out.slopes[cx,3] = u.species[s] #Species.code
    out.slopes[cx,4] = u.years[y] #Year
    out.slopes[cx,5] = n.y1 #Number of observations
    out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
    cx = cx + 1
  }
}

# This is a line I added in to remove any data from years past 2015.
out.slopes = out.slopes[which(out.slopes$Year <= 2015), ]

# Altered these tables based on what I currently have
write.table(out.slopes, 'Grasshoppers_SpeciesSlopes.txt',sep='\t',quote=F,row.names=F)
#####

# Merging slope files, here we can make so it is only the 2 Konza Prairie datasets of interest that are being merged
#####
# Merge slopes files

# edited from their code to try and bring in just the 2 Konza prairie files I have in a different spot than what 
# their code was expecting
files = list.files(pattern = 'SpeciesSlopes.txt')[2]


output = c()
for (f in 1:length(files)){
  print(f)
  dat1 = read.table(files[f],sep='\t',as.is=T,check.names=F,header=T)
  output = data.frame(rbind(output,dat1))
}
LLSs = paste(output$LTER.site,output$Locale,output$Species.code,sep='_')
LLS = unique(LLSs)
LLS.years = apply(array(LLS),1,function(x){length(which(LLSs==x))})
LLS.keep = LLS[which(LLS.years>3)]

output2 = file('PerSpecies_Abundance_LTER.csv','w')
writeLines('LTER.site,Locale,Species.code,Year,N.Obs,Abundance',output2)
for (i in 1:length(LLS.keep)){
  add1 = output[which(LLSs==LLS.keep[i]),]
  add2 = apply(add1,1,function(x){paste0(x,collapse=',')})
  writeLines(add2,output2)
}
close(output2)

#####


# I believe this is adding more information to the table and merging all of the different LTER dataframes together
#####
#################################################################Section
# Classify species according to habitat and feeding guild

# Import species annotations
eco.key = read.csv('/Users/caitlinmiller/Downloads/External_Database_S3/taxa_keys/Ecologial Variables.Arthropods.Final_v2.csv',as.is=T,check.names=F,header=T)
colnames(eco.key)[1] = 'Code'
eco.key$Species = trimws(eco.key$Species,which='right')
eco.key$Family = trimws(eco.key$Family,which='right')
eco.key$Order = trimws(eco.key$Order,which='right')

# Import species count data
output = read.csv('PerSpecies_Abundance_LTER.csv',as.is=T,check.names=F,header=T) #file created in "analyze.LTER.arthropods.R"
# Check ambiguous species codes - I assume this is making sure that these don't have any hits, in our case they don't so I think 
# we're good
output[which(output$Species.code=='' | output$Species.code=='Unknown' | output$Species.code=='unk unk' | output$Species.code=='undet under' | output$Species.code=='na na' | output$Species.code=='na? na?' | output$Species.code=='none none'),]

# Check LTER names
lters = unique(output$LTER)
lters2 = unique(eco.key$LTER)
# Not sure how relevant this is for us because we are only using the 1 LTER site
match(lters,lters2) #check LTER names

# Start new dataframe that will contain species categorizations
output2 = data.frame(output,'Order'=rep(NA,nrow(output)),'Family'=rep(NA,nrow(output)),'Species'=rep(NA,nrow(output)),'Feeding'=rep(NA,nrow(output)),'Habitat'=rep(NA,nrow(output)),'Pollinator'=rep(NA,nrow(output)))

# Define function that adds species annotations to new output table
add.annot = function(sps,annot.pos,LTER,ecokey,output2){
  if (length(sps)!=length(annot.pos)){
    print('Species and positions lists do not match')
  } else {
    for (z in 1:length(annot.pos)){
      if (is.na(annot.pos[z])){
      } else {
        ecokey2 = ecokey[annot.pos[z],]
        out.pos = which(output2$LTER==LTER & output2$Species.code==sps[z])
        output2[out.pos,7] = ecokey2$Order ###WARNING: sensitive to variable order!!!
        output2[out.pos,8] = ecokey2$Family
        output2[out.pos,9] = ecokey2$Species
        output2[out.pos,10] = ecokey2$Feeding
        output2[out.pos,11] = ecokey2$Habitat
        output2[out.pos,12] = ecokey2$Pollinator
      }
    }
  }
  return(output2)
}

# I am not actually sure what the add.annot() sections do because they don't seem 
# to save the information anywhere but including them so I don't have to refind them 
# later if they are important.
LTER = 'KonzaPrairie'
species2 = output$Species[which(output$LTER==LTER)]
u.species = unique(species2)
eco.key2 = eco.key[which(eco.key$LTER==LTER),]
sp.match = match(u.species,eco.key2$Species)
output2 = add.annot(sps=u.species,annot.pos=sp.match,LTER=LTER,ecokey=eco.key2,output2=output2)
unique(output2$Feeding[which(output2$LTER==LTER)])
tmp = output2[which(output2$LTER==LTER),] #Check missing
miss = tmp[which(is.na(tmp$Habitat)),]
sp.miss = unique(miss$Species.code)
miss.annot = cbind(c(NA,rep('Orthoptera',length(u.species)-1)),rep('herbivore',length(u.species)),rep('terrestrial',length(u.species)),rep('n',length(u.species)))
for (i in 1:length(sp.miss)){ #Annotate missing
  pos1 = which(output2$LTER==LTER & output2$Species.code==sp.miss[i])
  output2[pos1,7] = miss.annot[i,1]
  output2[pos1,10] = miss.annot[i,2]
  output2[pos1,11] = miss.annot[i,3]
  output2[pos1,12] = miss.annot[i,4]
}
tmp = output2[which(output2$LTER==LTER),] #Check missing
miss = tmp[which(is.na(tmp$Habitat)),]
sp.miss = unique(miss$Species.code) #PERFECT



###
# QC
length(which(apply(output2,1,function(x){length(which(is.na(x)))})>6))
LLSs = paste(output2$LTER.site,output2$Locale,output2$Species.code,sep='_')
na.count = apply(output2[,7:12],1,function(x){length(which(is.na(x)))})
LLSs.na = LLSs[which(na.count==6)]
u.LLS.na = unique(LLSs.na)

# Prune & polish
output2 = output2[which(na.count<6),which(colnames(output2)!='Family')]
output2$Order[which(output2$Order=='Ephemerotpera')] = 'Ephemeroptera'
output2$Order[which(output2$Order=='Hymenotptera')] = 'Hymenoptera'
output2$Order[which(output2$Order=='Tricophera')] = 'Trichoptera'
write.csv(output2,'PerSpecies_Abundance_LTER_annotated.csv',quote=F,row.names=F)
#####

# Their section labeled "curate and visualize RANK abundance data"
#####
library('vegan')
library('sads')

rank.list = list()
rx=1
# Get rank abundance curves
calc.ranks = function(sp1,ab1,rank.list,rx,name1){
  #sp1 = a vector of species names (species records do not have to be unique)
  #ab1 = a vector of corresponding abundances of species in sp1
  usp = unique(sp1) #unique species names
  spc = apply(array(usp),1,function(x){sum(ab1[which(sp1==x)],na.rm=T)}) #species counts
  rich = length(which(spc>0))
  if (rich==0){
    return(list(rank.list,rx)) #if species are all zero-abundant, do not calculate evenness
  } else {		
    #Rank abundance: rate of decay, Fisher's alpha, dominance (Berger-Parker Index)
    spc2 = sort(spc[which(spc>0)])
    if (length(spc2)<5){
      return(list(rank.list,rx))
    } else {
      mod = rad.preempt(spc2)
      rank.list[[rx]] = rev(spc2)
      names(rank.list)[rx] = name1
      rx=rx+1
      return(list(rank.list,rx))
    }
  }
}


site.name = 'KonzaPrairie'
locale.name = 'grasshopper'
data1 = read.csv("CGR022.csv",as.is=T,check.names=F,header=T)
data1$Species = data1$SPECIES
data1$Year = data1$RECYEAR
data1$Number = as.numeric(data1$TOTAL)
## pruning out entries without information for species
data1 = data1[which(data1$Species!='unknown' & data1$Species!='unknown ' & data1$Species!=''),]
u.years = sort(as.numeric(unique(data1$Year)))
u.species = unique(data1$Species)
# Fix species names
# Merge duplicates that have different species name spelling
konza.duplicates = read.csv('KonzaPrairie_grasshoppers_duplicatenames.csv',as.is=T,check.names=F,header=F)
data1 = data1[which(data1$Species!='Unknown'),]
for (i in 1:nrow(data1)){
  name.swap = konza.duplicates[which(konza.duplicates[,1]==data1$Species[i]),2]
  if (length(name.swap)==0){
    name.swap = konza.duplicates[match(data1$Species[i],konza.duplicates[,2]),2]
  } else {
  }
  data1$Species[i] = name.swap
  
}
u.species = unique(data1$Species)
u.years = sort(as.numeric(unique(data1$Year)))
u.years = list(u.years[1:2],u.years[(length(u.years)-1):(length(u.years))])
for (y in 1:length(u.years)){
  y.data = data1[which(data1$Year==u.years[[y]][1] | data1$Year==u.years[[y]][2]),]
  rank.out = calc.ranks(y.data$Species,y.data$Number,rank.list,rx,paste(paste(u.years[[y]][1],u.years[[y]][2],sep='-'),site.name,locale.name,sep='_'))
  rank.list = rank.out[[1]]
  rx=rank.out[[2]]
}

# Print curves
lter.sites = apply(array(names(rank.list)),1,function(x){y=strsplit(x,'_')[[1]];paste(y[2],y[3],sep='_')})
uls = unique(lter.sites)
uls.count = apply(array(uls),1,function(x){length(which(lter.sites==x))})
uls2 = uls[which(uls.count>1)]


# Come back to lines 3359-3375 to visualize rank abundance stuff if that becomes important


#####


# Estimate arthropod taxa abundance time trends
#####
# AR_reml function is used and defined for this

# First trying to have a full pipeline of the relaxed inclusion criteria processing and visualization.
# Relaxed pipeline is what they used to recapitulate their results so that is what we are sticking with.

# Import species abundance data
data1 = read.csv('PerSpecies_Abundance_LTER_annotated.csv',as.is=T,check.names=F,header=T)
data1$LL = paste(data1$LTER.site,data1$Locale,sep='_')
sp.count = apply(array(unique(data1$LL)),1,function(x){length(unique(data1$Species.code[which(data1$LL==x)]))})


LLS = paste(data1$LTER.site,data1$Locale,data1$Species.code,sep='_') #unique LTER by sub-site by species combinations
u.LLS = unique(LLS)
LLS.count = apply(array(u.LLS),1,function(x){d2=data1[which(LLS==x),];length(which(d2$Abundance!=0))})
u.LLS2 = u.LLS[which(LLS.count>0)]
time.trends.Z = data.frame('LTER'=NA,'Site'=NA,'Species'=NA,'MSE'=-999,'b'=-999,'coef_int'=-999,'coef_slope'=-999,'Pr_int'=-999,'Pr_slope'=-999,'logLik'=-999,'length'=-999,'Y1'=NA,'Y2'=NA,'Feeding'=NA,'Habitat'=NA,'Pollinator'=NA,'Order'=NA)
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
      add.coef = add.coef.Z = c(NA,NA,0,0,NA,NA,NA)
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
      arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
      add.coef.Z = c(arr.Z[[1]],arr.Z[[2]],arr.Z[[3]][1,1],arr.Z[[3]][2,1],arr.Z[[5]][1],arr.Z[[5]][2],arr.Z[[6]])
    }
    time.trends.Z[tx,1] = strsplit(u.LLS2[l],'_')[[1]][1] #LTER
    time.trends.Z[tx,2] = strsplit(u.LLS2[l],'_')[[1]][2] #Site
    time.trends.Z[tx,3] = paste0(strsplit(u.LLS2[l],'_')[[1]][3:length(strsplit(u.LLS2[l],'_')[[1]])],collapse='_') #Species
    time.trends.Z[tx,4] = add.coef.Z[1] #MSE
    time.trends.Z[tx,5] = add.coef.Z[2] #b
    time.trends.Z[tx,6] = add.coef.Z[3] #coef_int
    time.trends.Z[tx,7] = add.coef.Z[4] #coef_slope
    time.trends.Z[tx,8] = add.coef.Z[5] #Pr_int
    time.trends.Z[tx,9] = add.coef.Z[6] #Pr_slope
    time.trends.Z[tx,10] = add.coef.Z[7] #logLik
    time.trends.Z[tx,11] = length(which(!is.na(Z))) #length of time series
    time.trends.Z[tx,12] = lls.years.stretch[1] #first year in time series
    time.trends.Z[tx,13] = lls.years.stretch[length(lls.years.stretch)] #last year in time series
    time.trends.Z[tx,14] = lls.data$Feeding[1]
    time.trends.Z[tx,15] = lls.data$Habitat[1]
    time.trends.Z[tx,16] = lls.data$Pollinator[1]
    time.trends.Z[tx,17] = lls.data$Order[1]
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
write.csv(time.trends.Z,'time_trends_arthropods_relaxed.csv',quote=F,row.names=F)


library(scales)

trends = read.csv('time_trends_arthropods_relaxed.csv',as.is=T,check.names=F,header=T)

# This will be significantly pared down because we aren't looking at most of the sites


# To check data in violin plot, remove the gall insects row and then:
ggplot(trends, aes(LTER, coef_slope)) + geom_violin()

######


# Trying to change from average of slopes to abundance in each year
#####
# Need to use a file from a little ways up the pipeline. This is fully added on afterwards.

data1 = read.csv('PerSpecies_Abundance_LTER_annotated.csv',as.is=T,check.names=F,header=T)
data1$LL = paste(data1$LTER.site,data1$Locale,sep='_')

# Here we are agnostic to the species and just want to aggregate abundance counts over time.
u.year = unique(data1$Year)
all_abun = rep(NA, length(u.year))

for(y in 1:length(u.year))
{
  y_ind = which(data1$Year == u.year[y])
  all_abun[y] = sum(data1$Abundance[y_ind])
}


#####




# Functions they used
#####
# Function for calculating species richness, evenness (Pielou's index), and rank abundance curve metrics
calc.diversity = function(sp1,ab1){
  #sp1 = a vector of species names (species records do not have to be unique)
  #ab1 = a vector of corresponding abundances of species in sp1
  usp = unique(sp1) #unique species names
  # converting usp to an array here allows use of the apply fxn. The 1 indicates that there is only 1 dimension
  # we could be considering values from. Then, a function is applied that sums all counts from indices that have
  # the same names as entries in the unique names grouping to group duplicate counts.
  spc = apply(array(usp),1,function(x){sum(ab1[which(sp1==x)],na.rm=T)}) #species counts
  # I am not sure how you could be 0 spc over zero, maybe if nothing at all was measured??
  rich = length(which(spc>0))
  if (rich==0){
    return(c(0,rep(NA,4))) #if species are all zero-abundant, do not calculate evenness
  } else {
    spcp = spc / sum(spc,na.rm=T) #proportion of each species
    sdi = 0
    # sdi is based on how evenly distributed species are. i.e. 2 pops of .5 will lead to a much larger value than a pop
    # of .9 and 1 of .1. Here, there is a psuedolog situation where zero values are made non-zero but in a way we don't 
    # expect to have a large impact.
    for (z in 1:length(usp)){ #for each species
      if (spcp[z] == 0){
        spcp[z] = spcp[z] + 0.0001 #add miniscule number to zeroes
      } else {	
      }
      # Specifically here, smaller proportions will have a larger log value but be weighted down by being multiplied by
      # a small number (why the miniscule psuedolog values don't have too much importance) while larger props have a smaller
      # log number but it is increased by being multipled by a large number.
      sdi = sdi + (spcp[z] * log(spcp[z])) #SDI
    }
    # All numbers are negative, this gives us the positive abundance.
    sdi = sdi * (-1)
    
    pei = sdi / log(length(spcp)) #PEI
    
    #Rank abundance: rate of decay, Fisher's alpha, dominance (Berger-Parker Index)
    spc2 = sort(spc[which(spc>0)])
    if (length(spc2)<5){
      print('Too few species to calculate rank abundance')
      return(c(rich,pei,NA,NA,NA))
    } else {
      mod = rad.preempt(spc2)
      decay.rate = coef(mod)[1]
      dominance = spc2[length(spc2)] / sum(spc2) #counts of most abundant species divided by total counts across all species
      sf = fitsad(spc2, sad="ls")
      fisher.alpha = coef(sf)[2]
      return(c(rich,pei,decay.rate,fisher.alpha,dominance))
    }
  }
}

# AR_reml that we are trying to figure out

AR_reml <- function(formula, data = list()) {
  
  AR_reml_funct <- function(par, x, u) {
    # b is a set parameter
    b <- par
    n.obs <- length(x)
    # q is the number of columns that are being considered in that data relevant to the model
    q <- dim(u)[2]
    # A n by n matrix with dimensions based on the number of observations. Populated with zeros and then
    # 1s on the diagonal.
    B <- diag(n.obs)
    # populates the diagonal under the main diagonal with negative the set parameter.
    diag(B[-1, ]) <- -b
    
    # New matrix of zeroes n by n dimensions with zeros on the principal diagonal.
    iS <- diag(n.obs)
    # Changes the first value based on the set parameter value fed in.
    iS[1, 1] <- (1 - b^2)
    # Ultimately the weirdness with converting the first value in iS just leads to it allowing that value to
    # be 1 when converted in this way.
    iV <- t(B) %*% iS %*% B
    logdetV <- -determinant(iV)$modulus[1]
    
    coef <- solve(t(u) %*% iV %*% u, t(u) %*% iV %*% x)
    H <- x - u %*% coef
    
    s2 <- (t(H) %*% iV %*% H)/(n.obs - q)
    LL <- 0.5 * ((n.obs - q) * log(s2) + logdetV + determinant(t(u) %*% iV %*% u)$modulus[1] + (n.obs - q))
    #show(c(LL,b))
    return(LL)
  }
  
  
  mf <- model.frame(formula = formula, data = data)
  u <- model.matrix(attr(mf, "terms"), data = mf)
  x <- model.response(mf)
  
  q <- dim(u)[2]
  
  opt <- optim(fn = AR_reml_funct, par = 0.2, method = "Brent", upper = 1, lower = -1, control = list(maxit = 10^4), x = x, u = u)
  b <- opt$par
  
  n.obs <- length(x)
  q <- dim(u)[2]
  B <- diag(n.obs)
  diag(B[-1, ]) <- -b
  
  iS <- diag(n.obs)
  iS[1, 1] <- (1 - b^2)
  iV <- t(B) %*% iS %*% B
  logdetV <- -determinant(iV)$modulus[1]
  
  coef <- solve(t(u) %*% iV %*% u, t(u) %*% iV %*% x)
  H <- x - u %*% coef
  
  MSE <- as.numeric((t(H) %*% iV %*% H)/(n.obs - q))
  s2coef <- MSE * solve(t(u) %*% iV %*% u)
  Pr <- 1:q
  for (i in 1:q) Pr[i] <- 2 * pt(abs(coef[i])/s2coef[i, i]^0.5, df = n.obs - q, lower.tail = F)
  
  logLik <- 0.5 * (n.obs - q) * log(2 * pi) + determinant(t(u) %*% u)$modulus[1] - opt$value
  
  return(list(MSE = MSE, b = b, coef = coef, s2coef = s2coef, Pr = Pr, logLik = logLik))
}
#####