
# R code used to curate, analyze, and visualize time trends in arthropod abundance and biodiversity, and to curate/analyze/visualize bird/fish abundance data.
# Raw data used for these analyses came from the Long Term Ecological Research sites described in Table S1.

# Code assumes the working directory contains the following directories:
# './raw_data/'
# './summary_tables/arthropods/Species-level/'
# './summary_tables/arthropods/Site-level/'
# './summary_tables/vertebrates/'
# './plots/rank_abundance/'
# './plots/time_trends/abundance/lineplots/'
# './plots/time_trends/diversity/'
# './plots/time_trends/diversity/lineplots/'
# './plots/time_trends/diversity/fitted/'
# './plots/time_trends/vertebrates/'
# './plots/time_trends/vertebrates/lineplots/'
# './missing_species/'
# './taxa_keys/'
# './diversity_trends/'
# './shapefiles/'

# Code also assumes that WorldClim and Human Footprint Index maps have been downloaded from the following:
# http://biogeo.ucdavis.edu/data/worldclim/v2.0/tif/base/wc2.0_10m_bio.zip
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.052q5

# Navigate sections by using Ctrl-F to find the text: ###Section

# Questions can be directed to Dr. Michael Crossley: michael.crossley@uga.edu 


### BEGIN pipeline

library('vegan')
library('sads')

# Function for calculating species richness, evenness (Pielou's index), and rank abundance curve metrics
calc.diversity = function(sp1,ab1){
	#sp1 = a vector of species names (species records do not have to be unique)
	#ab1 = a vector of corresponding abundances of species in sp1
	usp = unique(sp1) #unique species names
	spc = apply(array(usp),1,function(x){sum(ab1[which(sp1==x)],na.rm=T)}) #species counts
	rich = length(which(spc>0))
	if (rich==0){
		return(c(0,rep(NA,4))) #if species are all zero-abundant, do not calculate evenness
	} else {
		spcp = spc / sum(spc,na.rm=T) #proportion of each species
		sdi = 0
		for (z in 1:length(usp)){ #for each species
			if (spcp[z] == 0){
				spcp[z] = spcp[z] + 0.0001 #add miniscule number to zeroes
			} else {	
			}
			sdi = sdi + (spcp[z] * log(spcp[z])) #SDI
		}
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



##### Harvard Forest ants 1
# NOTE: diversity estimates limited to 47 ant species
data1 = read.csv("./raw_data/hf118-01-ants.csv",as.is=T,check.names=F,header=T)
data1$Species = data1$code
data1$Year = data1$year
data1$Number = data1$abundance
data1 = data1[which(data1$Year!=2007),]
u.species = unique(data1$Species)
u.methods = unique(data1$trap.type)
#fix method labels
methods1 = data1$trap.type
methods1[which(methods1=='Pitfall' | methods1=='pitfall' | methods1=='pit')]='pitfall'
methods1[which(methods1=='Litter' | methods1=='LItter' | methods1=='litter')]='litter'
methods1[which(methods1=='Hand sample' | methods1=='hand')]='hand'
methods1[which(methods1=='Bait' | methods1=='Biat' | methods1=='bait')]='bait'
u.methods1 = unique(methods1)
u.locales = u.methods1
data1$Site = methods1
#Get values for LTER site-centric table
out.diversity = data.frame('LTER.site'=NA,'Locale'=NA,'Year'=-999,'N.obs'=-999,'Total.abundance'=-999,'N.species'=-999,'Species.evenness'=-999,'Species.decay.rate'=-999,'Fishers.alpha'=-999,'Dominance'=-999)
lx = 1
for (l in 1:length(u.locales)){
	u.years = sort(as.numeric(unique(data1$Year[which(data1$Site==u.locales[l])])))
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Site==u.locales[l]),]
		n.obs = 1 #standard sampling protocol was followed across years
		divs = calc.diversity(y.data$Species,y.data$Number)
		out.diversity[lx,1] = site.name #LTER
		out.diversity[lx,2] = paste0('ants.',u.locales[l]) #Locale
		out.diversity[lx,3] = u.years[y] #Year
		out.diversity[lx,4] = n.obs #Number of observations
		out.diversity[lx,5] = sum(y.data$Number,na.rm=T) #Total abundance
		out.diversity[lx,6] = divs[1] #Total number of species
		out.diversity[lx,7] = divs[2] #Evenness
		out.diversity[lx,8] = divs[3] #rank species abundance decay rate
		out.diversity[lx,9] = divs[4] #Fisher's alpha (measure of rarity)
		out.diversity[lx,10] = divs[5] #Dominance
		lx = lx + 1
	}
}
# Rarefied alpha and beta diversity
comm.list = list()
clx = 1
for (l in 1:length(u.locales)){ #create matrices of species abundance over time, for each "site" and/or "observation method"
	u.years = sort(as.numeric(unique(data1$Year[which(data1$Site==u.locales[l])])))
	comm.mat = matrix(NA,nrow=length(u.years),ncol=length(u.species))
	rownames(comm.mat) = u.years
	colnames(comm.mat) = u.species
	for (s in 1:length(u.species)){
		for (y in 1:length(u.years)){
			y.dat = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s] & data1$Site==u.locales[l]),]
			comm.mat[y,s] = sum(y.dat$Number,na.rm=T)
		}
	}
	comm.list[[clx]] = comm.mat
	clx = clx + 1
}
names(comm.list) = paste('ants',u.locales,sep='.')
# Rarefy S and PIE
add.Srare = data.frame('N.species.rarefied'=rep(NA,nrow(out.diversity)))
for (l in 1:length(comm.list)){
	comm.mat = comm.list[[l]]
	S <- specnumber(comm.mat) # observed number of species
	comm.mat2 = comm.mat[which(S>0),] #remove years with 0 species
	raremax <- min(rowSums(comm.mat2))
	Srare <- rarefy(comm.mat2, raremax)
	for (y in 1:length(Srare)){ #create vector of rarefied S that matches order in out.diversity
		row.index = which(out.diversity$Locale==names(comm.list)[l] & out.diversity$Year==names(Srare)[y]) ###Varies by LTER
		add.Srare[row.index,1] = Srare[y]
	}
}
out.diversity = data.frame(out.diversity,add.Srare)
# Beta diversity
out.beta = data.frame('LTER.site'=NA,'Locale'=NA,'Year1'=-999,'Year2'=-999,'Beta.2'=-999,'Beta.j'=-999,'Beta.bray'=-999)
bx = 1
for (l in 1:length(comm.list)){
	comm.mat = comm.list[[l]]
	for (y in 1:(nrow(comm.mat)-1)){
		dat1 = comm.mat[c(y,y+1),]
		row.sums = rowSums(dat1)
		if (length(which(row.sums==0))>0){
			beta.2 = beta.j = vd.b = NA
		} else {
			beta.2 = betadiver(dat1,method="-2")	#presence/absence-based dissimilarity, according to Harrison et al. (1992)
			beta.j = 1-betadiver(dat1,method="j") #presence/absence-based, Jaccard (dis)similarity 
			vd.b = vegdist(dat1,method='bray') #abundance-based estimate of beta
		}
		out.beta[bx,1] = site.name #LTER
		out.beta[bx,2] = names(comm.list)[l] #locale
		out.beta[bx,3] = rownames(comm.mat)[y] #year1
		out.beta[bx,4] = rownames(comm.mat)[y+1] #year2
		out.beta[bx,5] = beta.2 #beta.2
		out.beta[bx,6] = beta.j #beta.j
		out.beta[bx,7] = vd.b #beta.bray
		bx = bx + 1
	}
}
#Get values for species-centric table
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (l in 1:length(u.locales)){
	u.years = sort(as.numeric(unique(data1$Year[which(data1$Site==u.locales[l])])))
	for (s in 1:length(u.species)){
		for (y in 1:length(u.years)){
			y.data = data1[which(data1$Year==u.years[y] & data1$Site==u.locales[l] & data1$Species==u.species[s]),]
			n.y1 = 1
			out.slopes[cx,1] = site.name #LTER.site
			out.slopes[cx,2] = paste0('ants.',u.locales[l]) #Locale
			out.slopes[cx,3] = u.species[s] #Species.code
			out.slopes[cx,4] = u.years[y] #Year
			out.slopes[cx,5] = n.y1 #Number of observations
			out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
			cx = cx + 1
		}
	}
}
write.table(out.diversity,paste0('./summary_tables/arthropods/Site-level/',site.name,'_ants1_TotalChange.txt'),sep='\t',quote=F,row.names=F)
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_ants1_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)
write.table(out.beta,paste0('./summary_tables/arthropods/Site-level/',site.name,'_ants1_BetaDiv.txt'),sep='\t',quote=F,row.names=F)



##### Harvard Forest ants 3 (Nantucket)
# NOTE: diversity estimates limited to 61 ant species
site.name = 'HarvardForest'
locale.name = 'ants.Nantucket'
data1 = read.csv("./raw_data/hf147-08-nantucket-ants-2004-09.csv",as.is=T,check.names=F,header=T)
data1$Species = data1$code
data1$Site = data1$site
data1$Year = data1$year
data1$Number = data1$qty
data1 = data1[which(!is.na(data1$Species)),]
u.species = unique(data1$Species)
u.years = sort(as.numeric(unique(data1$year)))
#Get values for LTER site-centric table
out.diversity = data.frame('LTER.site'=NA,'Locale'=NA,'Year'=-999,'N.obs'=-999,'Total.abundance'=-999,'N.species'=-999,'Species.evenness'=-999,'Species.decay.rate'=-999,'Fishers.alpha'=-999,'Dominance'=-999)
lx = 1
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[y]),]
	n.obs = length(unique(paste(y.data$month,y.data$day,y.data$community.type,sep='_')))
	divs = calc.diversity(y.data$Species,y.data$Number)
	out.diversity[lx,1] = site.name #LTER
	out.diversity[lx,2] = locale.name #Locale
	out.diversity[lx,3] = u.years[y] #Year
	out.diversity[lx,4] = n.obs #Number of observations
	out.diversity[lx,5] = sum(y.data$Number,na.rm=T) #Total abundance
	out.diversity[lx,6] = divs[1] #Total number of species
	out.diversity[lx,7] = divs[2] #Evenness
	out.diversity[lx,8] = divs[3] #rank species abundance decay rate
	out.diversity[lx,9] = divs[4] #Fisher's alpha (measure of rarity)
	out.diversity[lx,10] = divs[5] #Dominance
	lx = lx + 1
}
# Rarefied alpha and beta diversity
comm.list = list()
clx = 1
comm.mat = matrix(NA,nrow=length(u.years),ncol=length(u.species))
rownames(comm.mat) = u.years
colnames(comm.mat) = u.species
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.dat = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
		comm.mat[y,s] = sum(y.dat$Number,na.rm=T)
	}
}
comm.list[[clx]] = comm.mat
names(comm.list) = locale.name
# Alpha diversity
add.Srare = data.frame('N.species.rarefied'=rep(NA,nrow(out.diversity)))
for (l in 1:length(comm.list)){
	comm.mat = comm.list[[l]]
	S <- specnumber(comm.mat) # observed number of species
	if (sum(S)==0){
		Srare = rep(NA,length(S)) #do not rarefy communities with zero abundances
	} else {
		comm.mat2 = comm.mat[which(S>0),] #remove years with 0 species
		if (length(which(S>0))==1){ #do not rarefy when locale has less than 2 years with non-zero data
			Srare = rep(NA,length(S))
		} else {
			raremax <- min(rowSums(comm.mat2))
			Srare <- rarefy(comm.mat2, raremax)
		}
	}
	for (y in 1:length(Srare)){ #create vector of rarefied S that matches order in out.diversity
		row.index = which(out.diversity$Locale==names(comm.list)[l] & out.diversity$Year==names(Srare)[y]) ###Varies by LTER
		add.Srare[row.index,1] = Srare[y]
	}
}
out.diversity = data.frame(out.diversity,add.Srare)
# Beta diversity
out.beta = data.frame('LTER.site'=NA,'Locale'=NA,'Year1'=-999,'Year2'=-999,'Beta.2'=-999,'Beta.j'=-999,'Beta.bray'=-999)
bx = 1
for (l in 1:length(comm.list)){
	comm.mat = comm.list[[l]]
	for (y in 1:(nrow(comm.mat)-1)){
		dat1 = comm.mat[c(y,y+1),]
		row.sums = rowSums(dat1)
		if (length(which(row.sums==0))>0){
			beta.2 = beta.j = vd.b = NA
		} else {
			beta.2 = betadiver(dat1,method="-2")	#presence/absence-based dissimilarity, according to Harrison et al. (1992)
			beta.j = 1-betadiver(dat1,method="j") #presence/absence-based, Jaccard (dis)similarity 
			vd.b = vegdist(dat1,method='bray') #abundance-based estimate of beta
		}
		out.beta[bx,1] = site.name #LTER
		out.beta[bx,2] = names(comm.list)[l] ###varies with LTER
		out.beta[bx,3] = rownames(comm.mat)[y] #year1
		out.beta[bx,4] = rownames(comm.mat)[y+1] #year2
		out.beta[bx,5] = beta.2 #beta.2
		out.beta[bx,6] = beta.j #beta.j
		out.beta[bx,7] = vd.b #beta.bray
		bx = bx + 1
	}
}
#Get values for species-centric table
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
u.species = u.species[which(!is.na(u.species))]
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
		n.y1 = length(unique(paste(y.data$month,y.data$day,y.data$community.type,sep='_')))
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.diversity,paste0('./summary_tables/arthropods/Site-level/',site.name,'_ants3_TotalChange.txt'),sep='\t',quote=F,row.names=F)
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_ants3_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)
write.table(out.beta,paste0('./summary_tables/arthropods/Site-level/',site.name,'_ants3_BetaDiv.txt'),sep='\t',quote=F,row.names=F)



##### Midwest Aphid Suction Trap Network (formerly "Kellog Biological Station")
site.name = 'MidwestSTN'
data1 = read.table('C:/Users/mcros/Desktop/Postdoc UGA/Aphid_STN/STN_counts_curated_time-consistent.txt',sep='\t',as.is=T,check.names=F,header=T)
data1$Number = apply(data1[,6:7],1,function(x){if(is.na(x[1]) & is.na(x[2])){NA}else{sum(x,na.rm=T)}})
u.years = sort(as.numeric(unique(data1$Year)))
u.species = unique(data1$Species)
u.locales = unique(data1$Site)
data1$Site[grep('Urbana-Champaign',data1$Site)] = 'Urbana-Champaign' #merge trap sites that were moved trivial distances midway through the time series
data1$Site[grep('Saginaw',data1$Site)] = 'Saginaw'
u.locales = sort(unique(data1$Site))

#Get values for LTER site-centric table
out.diversity = data.frame('LTER.site'=NA,'Locale'=NA,'Year'=-999,'N.obs'=-999,'Total.abundance'=-999,'N.species'=-999,'Species.evenness'=-999,'Species.decay.rate'=-999,'Fishers.alpha'=-999,'Dominance'=-999)
lx = 1
for (l in 1:length(u.locales)){
	print(u.locales[l])
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Site==u.locales[l]),]
#		y.data$Number = ceiling(y.data$Number)
#		n.obs = length(unique(y.data$Sample.date))
		n.obs = length(unique(y.data$Date))
		divs = calc.diversity(y.data$Species,y.data$Number)
		out.diversity[lx,1] = site.name #LTER
		out.diversity[lx,2] = u.locales[l] #Locale
		out.diversity[lx,3] = u.years[y] #Year
		out.diversity[lx,4] = n.obs #Number of observations
		out.diversity[lx,5] = sum(y.data$Number,na.rm=T) #Total abundance
		out.diversity[lx,6] = divs[1] #Total number of species
		out.diversity[lx,7] = divs[2] #Evenness
		out.diversity[lx,8] = divs[3] #rank species abundance decay rate
		out.diversity[lx,9] = divs[4] #Fisher's alpha (measure of rarity)
		out.diversity[lx,10] = divs[5] #Dominance
		lx = lx + 1
	}
}
# Rarefied alpha and beta diversity
comm.list = list()
clx = 1
for (l in 1:length(u.locales)){ #create matrices of species abundance over time, for each "site" and/or "observation method"
	print(l)
	comm.mat = matrix(NA,nrow=length(u.years),ncol=length(u.species))
	rownames(comm.mat) = u.years
	colnames(comm.mat) = u.species
	for (s in 1:length(u.species)){
		for (y in 1:length(u.years)){
			y.dat = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s] & data1$Site==u.locales[l]),]
			comm.mat[y,s] = sum(y.dat$Number,na.rm=T)
		}
	}
	comm.list[[clx]] = comm.mat
	clx = clx + 1
}
names(comm.list) = u.locales
# Alpha diversity
add.Srare = data.frame('N.species.rarefied'=rep(NA,nrow(out.diversity)))
for (l in 1:length(comm.list)){
	comm.mat = comm.list[[l]]
	S <- specnumber(comm.mat) # observed number of species
	if (sum(S)==0){
		print('zero')
		Srare = rep(NA,length(S)) #do not rarefy communities with zero abundances
	} else {
		comm.mat2 = comm.mat[which(S>0),] #remove years with 0 species
		comm.mat2 = ceiling(comm.mat2)
		if (length(which(S>0))==1){ #do not rarefy when locale has less than 2 years with non-zero data
			print('Zero')
			Srare = rep(NA,length(S))
		} else {
			raremax <- min(rowSums(comm.mat2))
			Srare <- rarefy(comm.mat2, raremax)
		}
	}
	for (y in 1:length(Srare)){ #create vector of rarefied S that matches order in out.diversity
		row.index = which(out.diversity$Locale==names(comm.list)[l] & out.diversity$Year==names(Srare)[y]) ###Varies by LTER
		add.Srare[row.index,1] = Srare[y]
	}
}
out.diversity = data.frame(out.diversity,add.Srare)
#out.diversity2 = out.diversity[which(out.diversity$N.obs>0),]

# Beta diversity
out.beta = data.frame('LTER.site'=NA,'Locale'=NA,'Year1'=-999,'Year2'=-999,'Beta.2'=-999,'Beta.j'=-999,'Beta.bray'=-999)
bx = 1
for (l in 1:length(comm.list)){
	comm.mat = comm.list[[l]]
	for (y in 1:(nrow(comm.mat)-1)){
		dat1 = comm.mat[c(y,y+1),]
		dat1 = ceiling(dat1)
		row.sums = rowSums(dat1)
		if (length(which(row.sums==0))>0){
			beta.2 = beta.j = vd.b = NA
		} else {
			beta.2 = betadiver(dat1,method="-2")	#presence/absence-based dissimilarity, according to Harrison et al. (1992)
			beta.j = 1-betadiver(dat1,method="j") #presence/absence-based, Jaccard (dis)similarity 
			vd.b = vegdist(dat1,method='bray') #abundance-based estimate of beta
		}
		out.beta[bx,1] = site.name #LTER
		out.beta[bx,2] = names(comm.list)[l] ###varies with LTER
		out.beta[bx,3] = rownames(comm.mat)[y] #year1
		out.beta[bx,4] = rownames(comm.mat)[y+1] #year2
		out.beta[bx,5] = beta.2 #beta.2
		out.beta[bx,6] = beta.j #beta.j
		out.beta[bx,7] = vd.b #beta.bray
		bx = bx + 1
	}
}
#Get values for species-centric table
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (l in 1:length(u.locales)){
	print(l)
	l.data = data1[which(data1$Site==u.locales[l]),]
	u.species = sort(unique(l.data$Species))
	for (s in 1:length(u.species)){
		s.data = l.data[which(l.data$Species==u.species[s]),]
		u.years = sort(as.numeric(unique(s.data$Year)))
		if (length(u.years)>3){
			for (y in 1:length(u.years)){
				y.data = s.data[which(s.data$Year==u.years[y]),]
				n.y1 = length(unique(y.data$Date))
				out.slopes[cx,1] = site.name #LTER.site
				out.slopes[cx,2] = u.locales[l] #Locale
				out.slopes[cx,3] = u.species[s] #Species.code
				out.slopes[cx,4] = u.years[y] #Year
				out.slopes[cx,5] = n.y1 #Number of observations
				out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
				cx = cx + 1
			}
		} else { #skip species-sites that do not have at least 4 years of data
		}
	}
}
write.table(out.diversity,paste0('./summary_tables/arthropods/Site-level/',site.name,'_STN_TotalChange.txt'),sep='\t',quote=F,row.names=F)
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_STN_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)
write.table(out.beta,paste0('./summary_tables/arthropods/Site-level/',site.name,'_STN_BetaDiv.txt'),sep='\t',quote=F,row.names=F)



##### Cedar Creek Ecosystem - Aphid vac
# EXCLUDED: only 3 years of data for "aphids"
f=7
site.name = gsub(' ','',files.key$LTER[which(files.key$File==sub('./raw_data/','',files.arthropods[f]))])
files.arthropods[f]
data1 = read.table(files.arthropods[f],sep='\t',as.is=T,check.names=F,header=T)
data1$Number = data1$Aphids
u.years = sort(as.numeric(unique(data1$Year)))


##### Coweeta - aquatic invertebrates
# NOTE: Values reported in count/m2, assuming equal number of m2 observations across years and sites
site.name = 'Coweeta'
locale.name = 'aquatic'
data1 = read.csv("./raw_data/3023_invertebrate_1_6046.CSV",as.is=T,check.names=F,header=T)
data1 = data1[-c(1:2),c(2:4,grep('Abundance',colnames(data1)))]
data1 = data1[,-grep('Total',colnames(data1))]
data1$Year = data1$Starting_Year
u.species = colnames(data1)[4:13]
u.years = sort(as.numeric(unique(data1$Year)))
#u.locales = unique(data1$Site)
#Get values for species-centric table
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y]),c(3,14,which(colnames(data1)==u.species[s]))]
		n.y1 = 1
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(as.numeric(y.data[,3]),na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_AquaticInverts_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)



##### Central Arizona-Phoenix - Arthropod sweeps
# NOTE: insects are identified to various levels of taxonomy; estimates of diversity will be downward-biased
site.name = 'CentralArizona-Phoenix'
locale.name = 'sweep'
data1 = read.csv("./raw_data/652_arthropods_e9c22403e0f7243f241ed64862e22e05.csv",as.is=T,check.names=F,header=T)
data1$Year = apply(array(data1$sample_date),1,function(x){strsplit(x,'-')[[1]][1]})
data1$Site = data1$site_code
data1$Species = data1$arthropod_scientific_name
data1$Species = gsub('\\(immature\\)','',data1$Species)
data1$Species = gsub(' \\(immature\\)','',data1$Species)
data1$Species = gsub('Scale Insects \\(Coccoidea\\)','',data1$Species)
data1$Species = trimws(data1$Species,which='right')
data1 = data1[which(data1$Species!='' & data1$Species!='Unknown'),]
data1$Number = data1$number_of_arthropods
u.years = sort(as.numeric(unique(data1$Year)))
u.locales = unique(data1$Site)
u.species = unique(data1$Species)
#Get values for LTER site-centric table
out.diversity = data.frame('LTER.site'=NA,'Locale'=NA,'Year'=-999,'N.obs'=-999,'Total.abundance'=-999,'N.species'=-999,'Species.evenness'=-999,'Species.decay.rate'=-999,'Fishers.alpha'=-999,'Dominance'=-999)
lx = 1
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[y]),]
	n.obs = length(unique(y.data$sample_date))
	divs = calc.diversity(y.data$Species,y.data$Number)
	out.diversity[lx,1] = site.name #LTER
	out.diversity[lx,2] = locale.name #Locale
	out.diversity[lx,3] = u.years[y] #Year
	out.diversity[lx,4] = n.obs #Number of observations
	out.diversity[lx,5] = sum(y.data$Number,na.rm=T) #Total abundance
	out.diversity[lx,6] = divs[1] #Total number of species
	out.diversity[lx,7] = divs[2] #Evenness
	out.diversity[lx,8] = divs[3] #rank species abundance decay rate
	out.diversity[lx,9] = divs[4] #Fisher's alpha (measure of rarity)
	out.diversity[lx,10] = divs[5] #Dominance
	lx = lx + 1
}
# Rarefied alpha and beta diversity
comm.list = list()
comm.mat = matrix(NA,nrow=length(u.years),ncol=length(u.species))
rownames(comm.mat) = u.years
colnames(comm.mat) = u.species
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.dat = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
		comm.mat[y,s] = sum(y.dat$Number,na.rm=T)
	}
}
comm.list[[1]] = comm.mat
names(comm.list) = locale.name
# Alpha diversity
add.Srare = data.frame('N.species.rarefied'=rep(NA,nrow(out.diversity)))
for (l in 1:length(comm.list)){
	comm.mat = comm.list[[l]]
	S <- specnumber(comm.mat) # observed number of species
	if (sum(S)==0){
		print('zero')
		Srare = rep(NA,length(S)) #do not rarefy communities with zero abundances
	} else {
		comm.mat2 = comm.mat[which(S>0),] #remove years with 0 species
		if (length(which(S>0))==1){ #do not rarefy when locale has less than 2 years with non-zero data
			print('Zero')
			Srare = rep(NA,length(S))
		} else {
			raremax <- min(rowSums(comm.mat2))
			Srare <- rarefy(comm.mat2, raremax)
		}
	}
	for (y in 1:length(Srare)){ #create vector of rarefied S that matches order in out.diversity
		row.index = which(out.diversity$Locale==names(comm.list)[l] & out.diversity$Year==names(Srare)[y]) ###Varies by LTER
		add.Srare[row.index,1] = Srare[y]
	}
}
out.diversity = data.frame(out.diversity,add.Srare)
# Beta diversity
out.beta = data.frame('LTER.site'=NA,'Locale'=NA,'Year1'=-999,'Year2'=-999,'Beta.2'=-999,'Beta.j'=-999,'Beta.bray'=-999)
bx = 1
for (l in 1:length(comm.list)){
	comm.mat = comm.list[[l]]
	for (y in 1:(nrow(comm.mat)-1)){
		dat1 = comm.mat[c(y,y+1),]
		row.sums = rowSums(dat1)
		if (length(which(row.sums==0))>0){
			beta.2 = beta.j = vd.b = NA
		} else {
			beta.2 = betadiver(dat1,method="-2")	#presence/absence-based dissimilarity, according to Harrison et al. (1992)
			beta.j = 1-betadiver(dat1,method="j") #presence/absence-based, Jaccard (dis)similarity 
			vd.b = vegdist(dat1,method='bray') #abundance-based estimate of beta
		}
		out.beta[bx,1] = site.name #LTER
		out.beta[bx,2] = names(comm.list)[l] ###varies with LTER
		out.beta[bx,3] = rownames(comm.mat)[y] #year1
		out.beta[bx,4] = rownames(comm.mat)[y+1] #year2
		out.beta[bx,5] = beta.2 #beta.2
		out.beta[bx,6] = beta.j #beta.j
		out.beta[bx,7] = vd.b #beta.bray
		bx = bx + 1
	}
}
#Get values for species-centric table
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
		n.y1 = length(unique(y.data$sample_date))
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.diversity,paste0('./summary_tables/arthropods/Site-level/',site.name,'_Sweep_TotalChange.txt'),sep='\t',quote=F,row.names=F)
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_Sweep_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)
write.table(out.beta,paste0('./summary_tables/arthropods/Site-level/',site.name,'_Sweep_BetaDiv.txt'),sep='\t',quote=F,row.names=F)



##### Cedar Creek Ecosystem Arthropods Sweep 1
site.name = 'CedarCreek'
locale.name = 'sweep1'
data1 = read.table("./raw_data/e153_Arthropod sweepnet sampling.txt",sep='\t',as.is=T,check.names=F,header=T)
data1$Year = apply(array(data1$Date),1,function(x){strsplit(x,'/')[[1]][3]})
data1$Species = paste(data1$Genus,data1$Specific.epithet,sep=' ')
data1 = data1[which(data1$Species!='undet undet' & data1$Species!='undet under'),]
data1$Number = data1$Specimens
u.years = sort(as.numeric(unique(data1$Year)))
u.species = unique(data1$Species)
#Get values for LTER site-centric table
out.diversity = data.frame('LTER.site'=NA,'Locale'=NA,'Year'=-999,'N.obs'=-999,'Total.abundance'=-999,'N.species'=-999,'Species.evenness'=-999,'Species.decay.rate'=-999,'Fishers.alpha'=-999,'Dominance'=-999)
lx = 1
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[y]),]
	n.obs = length(unique(paste(y.data$Location,y.data$Plot,y.data$Date,sep='_')))
	divs = calc.diversity(y.data$Species,y.data$Number)
	out.diversity[lx,1] = site.name #LTER
	out.diversity[lx,2] = locale.name #Locale
	out.diversity[lx,3] = u.years[y] #Year
	out.diversity[lx,4] = n.obs #Number of observations
	out.diversity[lx,5] = sum(y.data$Number,na.rm=T) #Total abundance
	out.diversity[lx,6] = divs[1] #Total number of species
	out.diversity[lx,7] = divs[2] #Evenness
	out.diversity[lx,8] = divs[3] #rank species abundance decay rate
	out.diversity[lx,9] = divs[4] #Fisher's alpha (measure of rarity)
	out.diversity[lx,10] = divs[5] #Dominance
	lx = lx + 1
}
# Rarefied alpha and beta diversity
comm.list = list()
comm.mat = matrix(NA,nrow=length(u.years),ncol=length(u.species))
rownames(comm.mat) = u.years
colnames(comm.mat) = u.species
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		row.pos = which(data1$Year==u.years[y])
		if (length(row.pos)==0){
			comm.mat[y,s] = NA
		} else {
			y.dat = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
			comm.mat[y,s] = sum(y.dat$Number,na.rm=T)
		}
	}
}
comm.list[[1]] = comm.mat
names(comm.list) = locale.name
# Alpha diversity
add.Srare = data.frame('N.species.rarefied'=rep(NA,nrow(out.diversity)))
for (l in 1:length(comm.list)){
	print(names(comm.list)[l])
	comm.mat = comm.list[[l]]
	S <- specnumber(comm.mat) # observed number of species
	rkeep = which(!is.na(S))
	S = S[rkeep]
	comm.mat = comm.mat[rkeep,]
	if (sum(S)==0){
		print('Community has zero abundances')
		Srare = rep(NA,length(S)) #do not rarefy communities with zero abundances
	} else {
		comm.mat2 = comm.mat[which(S>0),] #remove years with 0 species
		if (length(which(S>0))==1){ #do not rarefy when locale has less than 2 years with non-zero data
			print('Community has zero abundances')
			Srare = rep(NA,length(S))
		} else {
			raremax <- min(rowSums(comm.mat2))
			Srare <- rarefy(comm.mat2, raremax)
		}
	}
	for (y in 1:length(Srare)){ #create vector of rarefied S that matches order in out.diversity
		row.index = which(out.diversity$Locale==names(comm.list)[l] & out.diversity$Year==names(Srare)[y]) ###Varies by LTER
		add.Srare[row.index,1] = Srare[y]
	}
}
out.diversity = data.frame(out.diversity,add.Srare)
# Beta diversity
out.beta = data.frame('LTER.site'=NA,'Locale'=NA,'Year1'=-999,'Year2'=-999,'Beta.2'=-999,'Beta.j'=-999,'Beta.bray'=-999)
bx = 1
for (l in 1:length(comm.list)){
	comm.mat = comm.list[[l]]
	for (y in 1:(nrow(comm.mat)-1)){
		dat1 = comm.mat[c(y,y+1),]
		row.sums = rowSums(dat1)
		if (length(which(row.sums==0))>0 | length(which(is.na(row.sums)))>0){
			beta.2 = beta.j = vd.b = NA
		} else {
			beta.2 = betadiver(dat1,method="-2")	#presence/absence-based dissimilarity, according to Harrison et al. (1992)
			beta.j = 1-betadiver(dat1,method="j") #presence/absence-based, Jaccard (dis)similarity 
			vd.b = vegdist(dat1,method='bray') #abundance-based estimate of beta
		}
		out.beta[bx,1] = site.name #LTER
		out.beta[bx,2] = names(comm.list)[l] ###varies with LTER
		out.beta[bx,3] = rownames(comm.mat)[y] #year1
		out.beta[bx,4] = rownames(comm.mat)[y+1] #year2
		out.beta[bx,5] = beta.2 #beta.2
		out.beta[bx,6] = beta.j #beta.j
		out.beta[bx,7] = vd.b #beta.bray
		bx = bx + 1
	}
}
#Get values for species-centric table
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
		n.y1 = length(unique(paste(y.data$SweepSampleEventID,y.data$DateSampled,sep='_')))
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.diversity,paste0('./summary_tables/arthropods/Site-level/',site.name,'_ArthropodSweep1_TotalChange.txt'),sep='\t',quote=F,row.names=F)
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_ArthropodSweep1_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)
write.table(out.beta,paste0('./summary_tables/arthropods/Site-level/',site.name,'_ArthropodSweep1_BetaDiv.txt'),sep='\t',quote=F,row.names=F)



##### Cedar Creek Ecosystem Arthropods Sweep 2
site.name = 'CedarCreek'
locale.name = 'sweep2'
data1 = read.table("./raw_data/e120_Main Plots All Arthropod Insect Sweepnet Sampling 1996-2006.txt",sep='\t',as.is=T,check.names=F,header=T)
data1$Species = paste(data1$Genus,data1$Specific.epithet,sep=' ')
data1 = data1[which(data1$Species!='undet_undet'),]
data1$Number = as.numeric(data1$Count)
data1 = data1[which(data1$Species!='undet undet' & data1$Species!='unk unk' & data1$Species!='none none' & data1$Species!='na? na?' & data1$Species!='na na'),]
u.years = sort(as.numeric(unique(data1$Year)))
u.species = unique(data1$Species)
#Get values for LTER site-centric table
out.diversity = data.frame('LTER.site'=NA,'Locale'=NA,'Year'=-999,'N.obs'=-999,'Total.abundance'=-999,'N.species'=-999,'Species.evenness'=-999,'Species.decay.rate'=-999,'Fishers.alpha'=-999,'Dominance'=-999)
lx = 1
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[y]),]
	n.obs = length(unique(paste(y.data$Month,y.data$Plot,sep='_')))
	divs = calc.diversity(y.data$Species,y.data$Number)
	out.diversity[lx,1] = site.name #LTER
	out.diversity[lx,2] = locale.name #Locale
	out.diversity[lx,3] = u.years[y] #Year
	out.diversity[lx,4] = n.obs #Number of observations
	out.diversity[lx,5] = sum(y.data$Number,na.rm=T) #Total abundance
	out.diversity[lx,6] = divs[1] #Total number of species
	out.diversity[lx,7] = divs[2] #Evenness
	out.diversity[lx,8] = divs[3] #rank species abundance decay rate
	out.diversity[lx,9] = divs[4] #Fisher's alpha (measure of rarity)
	out.diversity[lx,10] = divs[5] #Dominance
	lx = lx + 1
}
# Rarefied alpha and beta diversity
comm.list = list()
comm.mat = matrix(NA,nrow=length(u.years),ncol=length(u.species))
rownames(comm.mat) = u.years
colnames(comm.mat) = u.species
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		row.pos = which(data1$Year==u.years[y])
		if (length(row.pos)==0){
			comm.mat[y,s] = NA
		} else {
			y.dat = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
			comm.mat[y,s] = sum(y.dat$Number,na.rm=T)
		}
	}
}
comm.list[[1]] = comm.mat
names(comm.list) = locale.name
# Alpha diversity
add.Srare = data.frame('N.species.rarefied'=rep(NA,nrow(out.diversity)))
for (l in 1:length(comm.list)){
	print(names(comm.list)[l])
	comm.mat = comm.list[[l]]
	S <- specnumber(comm.mat) # observed number of species
	rkeep = which(!is.na(S))
	S = S[rkeep]
	comm.mat = comm.mat[rkeep,]
	if (sum(S)==0){ #do not rarefy communities with zero abundances for all years
		print('Community has all zero abundances')
		Srare = rep(NA,length(S)) 
	} else {
		comm.mat2 = comm.mat[which(S>0),] #remove years with 0 species
		if (length(which(S>0))==1){ #do not rarefy when locale has less than 2 years with non-zero data
			print('Community has zero abundances')
			Srare = rep(NA,length(S))
		} else {
			raremax <- min(rowSums(comm.mat2))
			Srare <- rarefy(comm.mat2, raremax)
		}
	}
	for (y in 1:length(Srare)){ #create vector of rarefied S that matches order in out.diversity
		row.index = which(out.diversity$Locale==names(comm.list)[l] & out.diversity$Year==names(Srare)[y]) ###Varies by LTER
		add.Srare[row.index,1] = Srare[y]
	}
}
#
out.diversity = data.frame(out.diversity,add.Srare)
# Beta diversity
out.beta = data.frame('LTER.site'=NA,'Locale'=NA,'Year1'=-999,'Year2'=-999,'Beta.2'=-999,'Beta.j'=-999,'Beta.bray'=-999)
bx = 1
for (l in 1:length(comm.list)){
	comm.mat = comm.list[[l]]
	for (y in 1:(nrow(comm.mat)-1)){
		dat1 = comm.mat[c(y,y+1),]
		row.sums = rowSums(dat1)
		if (length(which(row.sums==0))>0 | length(which(is.na(row.sums)))>0){
			beta.2 = beta.j = vd.b = NA
		} else {
			beta.2 = betadiver(dat1,method="-2")	#presence/absence-based dissimilarity, according to Harrison et al. (1992)
			beta.j = 1-betadiver(dat1,method="j") #presence/absence-based, Jaccard (dis)similarity 
			vd.b = vegdist(dat1,method='bray') #abundance-based estimate of beta
		}
		out.beta[bx,1] = site.name #LTER
		out.beta[bx,2] = names(comm.list)[l] ###varies with LTER
		out.beta[bx,3] = rownames(comm.mat)[y] #year1
		out.beta[bx,4] = rownames(comm.mat)[y+1] #year2
		out.beta[bx,5] = beta.2 #beta.2
		out.beta[bx,6] = beta.j #beta.j
		out.beta[bx,7] = vd.b #beta.bray
		bx = bx + 1
	}
}
#Get values for species-centric table
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
		n.y1 = length(unique(paste(y.data$Month,y.data$Plot,sep='_')))
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.diversity,paste0('./summary_tables/arthropods/Site-level/',site.name,'_ArthropodSweep2_TotalChange.txt'),sep='\t',quote=F,row.names=F)
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_ArthropodSweep2_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)
write.table(out.beta,paste0('./summary_tables/arthropods/Site-level/',site.name,'_ArthropodSweep2_BetaDiv.txt'),sep='\t',quote=F,row.names=F)


##### Harvard Forest - Arthropods Carnivorous Plants
# NOTE: used carnivorous plant species as a "site"/"method of observation"
# NOTE: most insects only classified to order. Estimates of diversity will be downward-biased
site.name = 'HarvardForest'
data1 = read.csv( "./raw_data/hf111-01-prey.csv",as.is=T,check.names=F,header=T)
data1 = data1[which(data1$units=='number_of_individuals'),c(1:26,29:31)]
data1$Year = apply(array(data1$study),1,function(x){y=strsplit(x,' ')[[1]];return(y[length(y)])})
data1$Site = data1$species
u.locales = unique(data1$Site) #list of carnivorous plant species sampled (use this as a "site")
u.species = colnames(data1)[5:29]
u.years = sort(as.numeric(unique(data1$Year)))
#Get values for LTER site-centric table
out.diversity = data.frame('LTER.site'=NA,'Locale'=NA,'Year'=-999,'N.obs'=-999,'Total.abundance'=-999,'N.species'=-999,'Species.evenness'=-999,'Species.decay.rate'=-999,'Fishers.alpha'=-999,'Dominance'=-999)
lx = 1
for (l in 1:length(u.locales)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Site==u.locales[l]),]
		n.obs = nrow(y.data)
		if (n.obs==0){
		} else {
			divs = calc.diversity(colnames(y.data)[5:29],apply(y.data[,5:29],2,function(x){sum(x,na.rm=T)}))
			out.diversity[lx,1] = site.name #LTER
			out.diversity[lx,2] = u.locales[l] #Locale
			out.diversity[lx,3] = u.years[y] #Year
			out.diversity[lx,4] = n.obs #Number of observations
			out.diversity[lx,5] = sum(as.numeric(unlist(y.data[,5:29])),na.rm=T) #Total abundance
			out.diversity[lx,6] = divs[1] #Total number of species
			out.diversity[lx,7] = divs[2] #Evenness
			out.diversity[lx,8] = divs[3] #rank species abundance decay rate
			out.diversity[lx,9] = divs[4] #Fisher's alpha (measure of rarity)
			out.diversity[lx,10] = divs[5] #Dominance
			lx = lx + 1
		}
	}
}
# Rarefied alpha and beta diversity
comm.list = list()
clx = 1
for (l in 1:length(u.locales)){ #create matrices of species abundance over time, for each "site" and/or "observation method"
	comm.mat = matrix(NA,nrow=length(u.years),ncol=length(u.species))
	rownames(comm.mat) = u.years
	colnames(comm.mat) = u.species
	for (s in 1:length(u.species)){
		for (y in 1:length(u.years)){
			row.pos = which(data1$Year==u.years[y] & data1$Site==u.locales[l])
			if (length(row.pos)==0){
				comm.mat[y,s] = NA
			} else {
				y.dat = sum(data1[row.pos,which(colnames(data1)==u.species[s])],na.rm=T)
				comm.mat[y,s] = as.numeric(y.dat)
			}
		}
	}
	comm.list[[clx]] = comm.mat
	clx = clx + 1
}
names(comm.list) = u.locales
# Alpha diversity
add.Srare = data.frame('N.species.rarefied'=rep(NA,nrow(out.diversity)))
for (l in 1:length(comm.list)){
	print(names(comm.list)[l])
	comm.mat = comm.list[[l]]
	S = specnumber(comm.mat) # observed number of species
	rkeep = which(!is.na(S))
	S = S[rkeep]
	comm.mat = comm.mat[rkeep,]
	if (sum(S)==0){
		print('Community has zero abundances')
		Srare = rep(NA,length(S)) #do not rarefy communities with zero abundances
	} else {
		if (length(which(S>0))<2){ #do not rarefy when locale has less than 2 years with non-zero data
			print('Community has < 2 years with non-zero abundances')
			Srare = rep(NA,length(S))
		} else {
			comm.mat2 = comm.mat[which(S>0),] #remove years with 0 species
			raremax <- min(rowSums(comm.mat2))
			Srare <- rarefy(comm.mat2, raremax)
		}
	}
	for (y in 1:length(Srare)){ #create vector of rarefied S that matches order in out.diversity
		row.index = which(out.diversity$Locale==names(comm.list)[l] & out.diversity$Year==names(Srare)[y]) ###Varies by LTER
		add.Srare[row.index,1] = Srare[y]
	}
}
#
out.diversity = data.frame(out.diversity,add.Srare)
# Beta diversity
out.beta = data.frame('LTER.site'=NA,'Locale'=NA,'Year1'=-999,'Year2'=-999,'Beta.2'=-999,'Beta.j'=-999,'Beta.bray'=-999)
bx = 1
for (l in 1:length(comm.list)){
	comm.mat = comm.list[[l]]
	for (y in 1:(nrow(comm.mat)-1)){
		dat1 = comm.mat[c(y,y+1),]
		row.sums = rowSums(dat1)
		if (length(which(row.sums==0))>0 | length(which(is.na(row.sums)))>0){
			beta.2 = beta.j = vd.b = NA
		} else {
			beta.2 = betadiver(dat1,method="-2")	#presence/absence-based dissimilarity, according to Harrison et al. (1992)
			beta.j = 1-betadiver(dat1,method="j") #presence/absence-based, Jaccard (dis)similarity 
			vd.b = vegdist(dat1,method='bray') #abundance-based estimate of beta
		}
		out.beta[bx,1] = site.name #LTER
		out.beta[bx,2] = names(comm.list)[l] ###varies with LTER
		out.beta[bx,3] = rownames(comm.mat)[y] #year1
		out.beta[bx,4] = rownames(comm.mat)[y+1] #year2
		out.beta[bx,5] = beta.2 #beta.2
		out.beta[bx,6] = beta.j #beta.j
		out.beta[bx,7] = vd.b #beta.bray
		bx = bx + 1
	}
}
#Get values for species-centric table
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (l in 1:length(u.locales)){
	for (s in 1:length(u.species)){
		for (y in 1:length(u.years)){
			y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s] & data1$Site==u.locales[l]),c(3,30:31,which(colnames(data1)==u.species[s]))]
			n.y1 = nrow(y.data)
			out.slopes[cx,1] = site.name #LTER.site
			out.slopes[cx,2] = u.locales[l] #Locale
			out.slopes[cx,3] = u.species[s] #Species.code
			out.slopes[cx,4] = u.years[y] #Year
			out.slopes[cx,5] = n.y1 #Number of observations
			out.slopes[cx,6] = sum(as.numeric(y.data[,4]),na.rm=T) #abundance
			cx = cx + 1
		}
	}
}
write.table(out.diversity,paste0('./summary_tables/arthropods/Site-level/',site.name,'_ArthroCarnPlant_TotalChange.txt'),sep='\t',quote=F,row.names=F)
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_ArthroCarnPlant_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)
write.table(out.beta,paste0('./summary_tables/arthropods/Site-level/',site.name,'_ArthroCarnPlant_BetaDiv.txt'),sep='\t',quote=F,row.names=F)



##### Bonanza Creek aspen leaf miners
site.name = 'BonanzaCreek'
data1 = read.csv("./raw_data/608_ALM_eggSurvey_2004-2015.txt",as.is=T,check.names=F,header=T)
data1$Number = data1$Total_esm
u.years = sort(as.numeric(unique(data1$Year)))
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[y]),]
	n.y1 = nrow(y.data)
	out.slopes[cx,1] = site.name #LTER.site
	out.slopes[cx,2] = 'leafminer' #Locale
	out.slopes[cx,3] = 'Phyllocnistis populiella' #Species.code
	out.slopes[cx,4] = u.years[y] #Year
	out.slopes[cx,5] = n.y1 #Number of observations
	out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
	cx = cx + 1
}
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_LeafMiner_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)



##### Bonanza Creek bark beetles
site.name = 'BonanzaCreek'
data1 = read.csv("./raw_data/35_BNZ_Beetles_Werner_1975-2012.txt",as.is=T,check.names=F,header=T)
d2010 = apply(data1[which(data1$Year==2010),2:4],2,function(x){mean(x,na.rm=T)})
data1 = data1[-which(data1$Year==2010),]
data1 = data.frame(rbind(data1,c(2010,d2010)))
u.years = sort(as.numeric(unique(data1$Year)))
u.species = colnames(data1)[2:4]
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y]),which(colnames(data1)==u.species[s])]
		n.y1 = 1
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = 'bark beetle' #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_BarkBeetles_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)



##### Konza Prairie gall insects
site.name = 'KonzaPrairie'
locale.name = 'gall'
data1 = read.csv("./raw_data/CGP011.csv",as.is=T,check.names=F,header=T)
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
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_GallInsect_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)


##### Konza Prairie grasshoppers
site.name = 'KonzaPrairie'
locale.name = 'grasshopper'
data1 = read.csv("./raw_data/CGR022.csv",as.is=T,check.names=F,header=T)
data1$Species = data1$SPECIES
data1$Year = data1$RECYEAR
data1$Number = as.numeric(data1$TOTAL)
data1 = data1[which(data1$Species!='unknown' & data1$Species!='unknown ' & data1$Species!=''),]
u.years = sort(as.numeric(unique(data1$Year)))
u.species = unique(data1$Species)
# Fix species names
# Merge duplicates that have different species name spelling
konza.duplicates = read.csv('./taxa_keys/KonzaPrairie_grasshoppers_duplicatenames.csv',as.is=T,check.names=F,header=F)
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
#Get values for LTER site-centric table
out.diversity = data.frame('LTER.site'=NA,'Locale'=NA,'Year'=-999,'N.obs'=-999,'Total.abundance'=-999,'N.species'=-999,'Species.evenness'=-999,'Species.decay.rate'=-999,'Fishers.alpha'=-999,'Dominance'=-999)
lx = 1
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[y]),]
	n.obs = length(unique(paste(y.data$RECMONTH,y.data$RECDAY,y.data$REPSITE,sep='_'))) #sweep sample events
	divs = calc.diversity(y.data$Species,y.data$Number)
	out.diversity[lx,1] = site.name #LTER
	out.diversity[lx,2] = locale.name #Locale
	out.diversity[lx,3] = u.years[y] #Year
	out.diversity[lx,4] = n.obs #Number of observations
	out.diversity[lx,5] = sum(y.data$Number,na.rm=T) #Total abundance
	out.diversity[lx,6] = divs[1] #Total number of species
	out.diversity[lx,7] = divs[2] #Evenness
	out.diversity[lx,8] = divs[3] #rank species abundance decay rate
	out.diversity[lx,9] = divs[4] #Fisher's alpha (measure of rarity)
	out.diversity[lx,10] = divs[5] #Dominance
	lx = lx + 1
}
# Rarefied alpha and beta diversity
comm.list = list()
comm.mat = matrix(NA,nrow=length(u.years),ncol=length(u.species))
rownames(comm.mat) = u.years
colnames(comm.mat) = u.species
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		row.pos = which(data1$Year==u.years[y])
		if (length(row.pos)==0){
			comm.mat[y,s] = NA
		} else {
			y.dat = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
			comm.mat[y,s] = sum(y.dat$Number,na.rm=T)
		}
	}
}
comm.list[[1]] = comm.mat
names(comm.list) = locale.name
# Alpha diversity
add.Srare = data.frame('N.species.rarefied'=rep(NA,nrow(out.diversity)))
for (l in 1:length(comm.list)){
	print(names(comm.list)[l])
	comm.mat = comm.list[[l]]
	S <- specnumber(comm.mat) # observed number of species
	rkeep = which(!is.na(S))
	S = S[rkeep]
	comm.mat = comm.mat[rkeep,]
	if (sum(S)==0){ #do not rarefy communities with zero abundances for all years
		print('Community has all zero abundances')
		Srare = rep(NA,length(S)) 
	} else {
		comm.mat2 = comm.mat[which(S>0),] #remove years with 0 species
		if (length(which(S>0))==1){ #do not rarefy when locale has less than 2 years with non-zero data
			print('Community has zero abundances')
			Srare = rep(NA,length(S))
		} else {
			raremax <- min(rowSums(comm.mat2))
			Srare <- rarefy(comm.mat2, raremax)
		}
	}
	for (y in 1:length(Srare)){ #create vector of rarefied S that matches order in out.diversity
		row.index = which(out.diversity$Locale==names(comm.list)[l] & out.diversity$Year==names(Srare)[y]) ###Varies by LTER
		add.Srare[row.index,1] = Srare[y]
	}
}
#
out.diversity = data.frame(out.diversity,add.Srare)
# Beta diversity
out.beta = data.frame('LTER.site'=NA,'Locale'=NA,'Year1'=-999,'Year2'=-999,'Beta.2'=-999,'Beta.j'=-999,'Beta.bray'=-999)
bx = 1
for (l in 1:length(comm.list)){
	comm.mat = comm.list[[l]]
	for (y in 1:(nrow(comm.mat)-1)){
		dat1 = comm.mat[c(y,y+1),]
		row.sums = rowSums(dat1)
		if (length(which(row.sums==0))>0 | length(which(is.na(row.sums)))>0){
			beta.2 = beta.j = vd.b = NA
		} else {
			beta.2 = betadiver(dat1,method="-2")	#presence/absence-based dissimilarity, according to Harrison et al. (1992)
			beta.j = 1-betadiver(dat1,method="j") #presence/absence-based, Jaccard (dis)similarity 
			vd.b = vegdist(dat1,method='bray') #abundance-based estimate of beta
		}
		out.beta[bx,1] = site.name #LTER
		out.beta[bx,2] = names(comm.list)[l] ###varies with LTER
		out.beta[bx,3] = rownames(comm.mat)[y] #year1
		out.beta[bx,4] = rownames(comm.mat)[y+1] #year2
		out.beta[bx,5] = beta.2 #beta.2
		out.beta[bx,6] = beta.j #beta.j
		out.beta[bx,7] = vd.b #beta.bray
		bx = bx + 1
	}
}
#Get values for species-centric table
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
write.table(out.diversity,paste0('./summary_tables/arthropods/Site-level/',site.name,'_Grasshoppers_TotalChange.txt'),sep='\t',quote=F,row.names=F)
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_Grasshoppers_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)
write.table(out.beta,paste0('./summary_tables/arthropods/Site-level/',site.name,'_Grasshoppers_BetaDiv.txt'),sep='\t',quote=F,row.names=F)


##### Cedar Creek Ecosystem - Old Field Grasshopper Sampling
site.name = 'CedarCreek'
locale.name = 'grasshopper'
data1 = read.table("./raw_data/e014_Core Old Field Grasshopper Sampling.txt",sep='\t',as.is=T,check.names=F,header=T)
data1$Species = paste(data1$Genus,data1$Specific.epithet,sep=' ')
data1$Number = data1$Specimens
u.years = sort(as.numeric(unique(data1$Year)))
u.species = unique(data1$Species)
#Get values for LTER site-centric table
out.diversity = data.frame('LTER.site'=NA,'Locale'=NA,'Year'=-999,'N.obs'=-999,'Total.abundance'=-999,'N.species'=-999,'Species.evenness'=-999,'Species.decay.rate'=-999,'Fishers.alpha'=-999,'Dominance'=-999)
lx = 1
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[y]),]
	n.obs = 1 #standardized sampling protocol
	divs = calc.diversity(y.data$Species,y.data$Number)
	out.diversity[lx,1] = site.name #LTER
	out.diversity[lx,2] = locale.name #Locale
	out.diversity[lx,3] = u.years[y] #Year
	out.diversity[lx,4] = n.obs #Number of observations
	out.diversity[lx,5] = sum(y.data$Number,na.rm=T) #Total abundance
	out.diversity[lx,6] = divs[1] #Total number of species
	out.diversity[lx,7] = divs[2] #Evenness
	out.diversity[lx,8] = divs[3] #rank species abundance decay rate
	out.diversity[lx,9] = divs[4] #Fisher's alpha (measure of rarity)
	out.diversity[lx,10] = divs[5] #Dominance
	lx = lx + 1
}
# Rarefied alpha and beta diversity
comm.list = list()
comm.mat = matrix(NA,nrow=length(u.years),ncol=length(u.species))
rownames(comm.mat) = u.years
colnames(comm.mat) = u.species
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		row.pos = which(data1$Year==u.years[y])
		if (length(row.pos)==0){
			comm.mat[y,s] = NA
		} else {
			y.dat = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
			comm.mat[y,s] = sum(y.dat$Number,na.rm=T)
		}
	}
}
comm.list[[1]] = comm.mat
names(comm.list) = locale.name
# Alpha diversity
add.Srare = data.frame('N.species.rarefied'=rep(NA,nrow(out.diversity)))
for (l in 1:length(comm.list)){
	print(names(comm.list)[l])
	comm.mat = comm.list[[l]]
	S <- specnumber(comm.mat) # observed number of species
	rkeep = which(!is.na(S))
	S = S[rkeep]
	comm.mat = comm.mat[rkeep,]
	if (sum(S)==0){ #do not rarefy communities with zero abundances for all years
		print('Community has all zero abundances')
		Srare = rep(NA,length(S)) 
	} else {
		comm.mat2 = comm.mat[which(S>0),] #remove years with 0 species
		if (length(which(S>0))==1){ #do not rarefy when locale has less than 2 years with non-zero data
			print('Community has zero abundances')
			Srare = rep(NA,length(S))
		} else {
			raremax <- min(rowSums(comm.mat2))
			Srare <- rarefy(comm.mat2, raremax)
		}
	}
	for (y in 1:length(Srare)){ #create vector of rarefied S that matches order in out.diversity
		row.index = which(out.diversity$Locale==names(comm.list)[l] & out.diversity$Year==names(Srare)[y]) ###Varies by LTER
		add.Srare[row.index,1] = Srare[y]
	}
}
#
out.diversity = data.frame(out.diversity,add.Srare)
# Beta diversity
out.beta = data.frame('LTER.site'=NA,'Locale'=NA,'Year1'=-999,'Year2'=-999,'Beta.2'=-999,'Beta.j'=-999,'Beta.bray'=-999)
bx = 1
for (l in 1:length(comm.list)){
	comm.mat = comm.list[[l]]
	for (y in 1:(nrow(comm.mat)-1)){
		dat1 = comm.mat[c(y,y+1),]
		row.sums = rowSums(dat1)
		if (length(which(row.sums==0))>0 | length(which(is.na(row.sums)))>0){
			beta.2 = beta.j = vd.b = NA
		} else {
			beta.2 = betadiver(dat1,method="-2")	#presence/absence-based dissimilarity, according to Harrison et al. (1992)
			beta.j = 1-betadiver(dat1,method="j") #presence/absence-based, Jaccard (dis)similarity 
			vd.b = vegdist(dat1,method='bray') #abundance-based estimate of beta
		}
		out.beta[bx,1] = site.name #LTER
		out.beta[bx,2] = names(comm.list)[l] ###varies with LTER
		out.beta[bx,3] = rownames(comm.mat)[y] #year1
		out.beta[bx,4] = rownames(comm.mat)[y+1] #year2
		out.beta[bx,5] = beta.2 #beta.2
		out.beta[bx,6] = beta.j #beta.j
		out.beta[bx,7] = vd.b #beta.bray
		bx = bx + 1
	}
}
#Get values for species-centric table
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
		n.y1 = 1
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.diversity,paste0('./summary_tables/arthropods/Site-level/',site.name,'_Grasshoppers_TotalChange.txt'),sep='\t',quote=F,row.names=F)
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_Grasshoppers_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)
write.table(out.beta,paste0('./summary_tables/arthropods/Site-level/',site.name,'_Grasshoppers_BetaDiv.txt'),sep='\t',quote=F,row.names=F)



##### Georgia Coastal Ecosystems - Grasshoppers
# NOTE: only 7 grasshopper species documented, skipping diversity statistics
site.name = 'GeorgiaCoastal'
locale.name = 'grasshopper'
files = list.files('./raw_data/GeorgiaCoast_grasshoppers',full.names=T)
data1 = c()
for (f in 1:length(files)){ #merge yearly data into one dataset
	add.data = read.table(files[f],sep='\t',as.is=T,check.names=F,header=T)[,1:7]
	data1 = data.frame(rbind(data1,add.data))
}
data1$Species = data1$Species_code
data1$Number = data1$Count
u.years = sort(as.numeric(unique(data1$Year)))
u.species = unique(data1$Species)
#Get values for species-centric table
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
		n.y1 = length(unique(paste(y.data$Month,y.data$Day,y.data$Transect,sep='_')))
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_Grasshoppers_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)



##### Sevilleta - Grasshoppers
site.name = 'Sevilleta'
locale.name = 'grasshopper'
data1 = read.csv("./raw_data/sev106_hopperdynamics_20150826.txt",as.is=T,check.names=F,header=T)
data1$Year = apply(array(data1$DATE),1,function(x){strsplit(x,'/')[[1]][3]})
data1$Site = data1$SITE
data1$Species = data1$SPECIES
data1$Number = data1$CNT
u.years = sort(as.numeric(unique(data1$Year)))
u.species = unique(data1$Species)
u.locales = unique(data1$Site)
#Get values for LTER site-centric table
out.diversity = data.frame('LTER.site'=NA,'Locale'=NA,'Year'=-999,'N.obs'=-999,'Total.abundance'=-999,'N.species'=-999,'Species.evenness'=-999,'Species.decay.rate'=-999,'Fishers.alpha'=-999,'Dominance'=-999)
lx = 1
#for (l in 1:length(u.locales)){
	l.data = data1
	u.years = sort(as.numeric(unique(l.data$Year)))
	for (y in 1:length(u.years)){
		y.data = l.data[which(l.data$Year==u.years[y] & l.data$Site==u.locales[l]),]
		n.obs = 1 #standardized protocol followed # n.obs = length(unique(paste(y.data$WEB,y.data$TRN,sep='_')))
		divs = calc.diversity(y.data$Species,y.data$Number)
		out.diversity[lx,1] = site.name #LTER
		out.diversity[lx,2] = locale.name #Locale
		out.diversity[lx,3] = u.years[y] #Year
		out.diversity[lx,4] = n.obs #Number of observations
		out.diversity[lx,5] = sum(y.data$Number,na.rm=T) #Total abundance
		out.diversity[lx,6] = divs[1] #Total number of species
		out.diversity[lx,7] = divs[2] #Evenness
		out.diversity[lx,8] = divs[3] #rank species abundance decay rate
		out.diversity[lx,9] = divs[4] #Fisher's alpha (measure of rarity)
		out.diversity[lx,10] = divs[5] #Dominance
		lx = lx + 1
	}
#}
# Rarefied alpha and beta diversity
comm.list = list()
clx = 1
comm.mat = matrix(NA,nrow=length(u.years),ncol=length(u.species))
rownames(comm.mat) = u.years
colnames(comm.mat) = u.species
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		row.pos = which(data1$Year==u.years[y] & data1$Species==u.species[s])
		if (length(row.pos)==0){
			comm.mat[y,s] = NA
		} else {
			y.dat = data1[row.pos,]
			comm.mat[y,s] = sum(y.dat$Number,na.rm=T)
		}
	}
}
comm.list[[clx]] = comm.mat
clx = clx + 1
names(comm.list) = locale.name
for (l in 1:length(comm.list)){ #replace NA with zero
	for (i in 1:nrow(comm.list[[l]])){
		comm.list[[l]][i,][which(is.na(comm.list[[l]][i,]))] = 0
	}
}
# Alpha diversity
add.Srare = data.frame('N.species.rarefied'=rep(NA,nrow(out.diversity)))
for (l in 1:length(comm.list)){
	print(names(comm.list)[l])
	comm.mat = comm.list[[l]]
	S <- specnumber(comm.mat) # observed number of species
	rkeep = which(!is.na(S))
	S = S[rkeep]
	comm.mat = comm.mat[rkeep,]
	if (sum(S)==0){ #do not rarefy communities with zero abundances for all years
		print('Community has all zero abundances')
		Srare = rep(NA,length(S)) 
	} else {
		comm.mat2 = comm.mat[which(S>0),] #remove years with 0 species
		if (length(which(S>0))==1){ #do not rarefy when locale has less than 2 years with non-zero data
			print('Community has zero abundances')
			Srare = rep(NA,length(S))
		} else {
			raremax <- min(rowSums(comm.mat2))
			Srare <- rarefy(comm.mat2, raremax)
		}
	}
	for (y in 1:length(Srare)){ #create vector of rarefied S that matches order in out.diversity
		row.index = which(out.diversity$Locale==names(comm.list)[l] & out.diversity$Year==names(Srare)[y]) ###Varies by LTER
		add.Srare[row.index,1] = Srare[y]
	}
}
#
out.diversity = data.frame(out.diversity,add.Srare)
# Beta diversity
out.beta = data.frame('LTER.site'=NA,'Locale'=NA,'Year1'=-999,'Year2'=-999,'Beta.2'=-999,'Beta.j'=-999,'Beta.bray'=-999)
bx = 1
for (l in 1:length(comm.list)){
	comm.mat = comm.list[[l]]
	for (y in 1:(nrow(comm.mat)-1)){
		dat1 = comm.mat[c(y,y+1),]
		row.sums = rowSums(dat1)
		if (length(which(row.sums==0))>0 | length(which(is.na(row.sums)))>0){
			beta.2 = beta.j = vd.b = NA
		} else {
			beta.2 = betadiver(dat1,method="-2")	#presence/absence-based dissimilarity, according to Harrison et al. (1992)
			beta.j = 1-betadiver(dat1,method="j") #presence/absence-based, Jaccard (dis)similarity 
			vd.b = vegdist(dat1,method='bray') #abundance-based estimate of beta
		}
		out.beta[bx,1] = site.name #LTER
		out.beta[bx,2] = names(comm.list)[l] ###varies with LTER
		out.beta[bx,3] = rownames(comm.mat)[y] #year1
		out.beta[bx,4] = rownames(comm.mat)[y+1] #year2
		out.beta[bx,5] = beta.2 #beta.2
		out.beta[bx,6] = beta.j #beta.j
		out.beta[bx,7] = vd.b #beta.bray
		bx = bx + 1
	}
}
#Get values for species-centric table
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
		n.y1 = 1
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.diversity,paste0('./summary_tables/arthropods/Site-level/',site.name,'_Grasshoppers_TotalChange.txt'),sep='\t',quote=F,row.names=F)
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_Grasshoppers_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)
write.table(out.beta,paste0('./summary_tables/arthropods/Site-level/',site.name,'_Grasshoppers_BetaDiv.txt'),sep='\t',quote=F,row.names=F)



##### Central Arizona-Phoenix - Arthropod Pitfall 1
# NOTE: taxon designations down to family level; diversity will be underestimated
site.name = 'CentralArizona-Phoenix'
locale.name = 'pitfall1'
data1 = read.csv("./raw_data/41_core_arthropods_1c843549f049757465771bf4c36e80ab.csv",as.is=T,check.names=F,header=T)
data1$Year = apply(array(data1$sample_date),1,function(x){strsplit(x,'-')[[1]][1]})
data1$Site = data1$site_code
data1$Species = data1$display_name
data1$Species = gsub('\\(adult\\)','',data1$Species)
data1$Species = gsub('\\(winged\\)','',data1$Species)
data1$Species = gsub('\\(larvae\\)','',data1$Species)
data1$Species = gsub('\\(nymph\\)','',data1$Species)
data1$Species = gsub(' \\(immature\\)','',data1$Species)
data1$Species = gsub(' \\(Winged Aphelinidae\\)','',data1$Species)
data1$Species = gsub(' \\(Aphelinidae\\)','',data1$Species)
data1$Species = gsub('_',' ',data1$Species)
data1$Species = gsub(' \\(Centipeds\\)','',data1$Species)
data1$Species = gsub('Scale Insects \\(Coccoidea\\)','',data1$Species)
data1$Species = gsub('Scale Insects \\(Coccoidea - Immature\\)','',data1$Species)
data1$Species = gsub(' \\(imm.\\)','',data1$Species)
data1$Species = gsub(' >10mm','',data1$Species)
data1$Species = gsub(' immature','',data1$Species)
data1$Species = trimws(data1$Species,which='right')
data1$Number = apply(data1[,13:17],1,function(x){sum(x,na.rm=T)})
data1 = data1[which(!is.na(data1$Species)),]
data1 = data1[which(data1$Species!='' & data1$Species!='Unknown'),]
u.locales = unique(data1$Site)
u.species = unique(data1$Species)
u.years = sort(as.numeric(unique(data1$Year)))
#Get values for LTER site-centric table
out.diversity = data.frame('LTER.site'=NA,'Locale'=NA,'Year'=-999,'N.obs'=-999,'Total.abundance'=-999,'N.species'=-999,'Species.evenness'=-999,'Species.decay.rate'=-999,'Fishers.alpha'=-999,'Dominance'=-999)
lx = 1
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[y]),]
	n.obs = 1 #length(unique(paste(y.data$sample_date,y.data$trap_name,sep='_')))
	divs = calc.diversity(y.data$Species,y.data$Number)
	out.diversity[lx,1] = site.name #LTER
	out.diversity[lx,2] = locale.name #Locale
	out.diversity[lx,3] = u.years[y] #Year
	out.diversity[lx,4] = n.obs #Number of observations
	out.diversity[lx,5] = sum(y.data$Number,na.rm=T) #Total abundance
	out.diversity[lx,6] = divs[1] #Total number of species
	out.diversity[lx,7] = divs[2] #Evenness
	out.diversity[lx,8] = divs[3] #rank species abundance decay rate
	out.diversity[lx,9] = divs[4] #Fisher's alpha (measure of rarity)
	out.diversity[lx,10] = divs[5] #Dominance
	lx = lx + 1
}
# Rarefied alpha and beta diversity
u.years = sort(as.numeric(unique(data1$Year)))
comm.list = list()
clx = 1
comm.mat = matrix(NA,nrow=length(u.years),ncol=length(u.species))
rownames(comm.mat) = u.years
colnames(comm.mat) = u.species
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		row.pos = which(data1$Year==u.years[y] & data1$Species==u.species[s])
		if (length(row.pos)==0){
			comm.mat[y,s] = NA
		} else {
			y.dat = data1[row.pos,]
			comm.mat[y,s] = sum(y.dat$Number,na.rm=T)
		}
	}
}
comm.list[[clx]] = comm.mat
clx = clx + 1
names(comm.list) = locale.name
for (l in 1:length(comm.list)){ #replace NA with zero
	for (i in 1:nrow(comm.list[[l]])){
		comm.list[[l]][i,][which(is.na(comm.list[[l]][i,]))] = 0
	}
}
# Alpha diversity
add.Srare = data.frame('N.species.rarefied'=rep(NA,nrow(out.diversity)))
for (l in 1:length(comm.list)){
	print(names(comm.list)[l])
	comm.mat = comm.list[[l]]
	S <- specnumber(comm.mat) # observed number of species
	rkeep = which(!is.na(S))
	S = S[rkeep]
	comm.mat = comm.mat[rkeep,]
	if (sum(S)==0){ #do not rarefy communities with zero abundances for all years
		print('Community has all zero abundances')
		Srare = rep(NA,length(S)) 
	} else {
		comm.mat2 = comm.mat[which(S>0),] #remove years with 0 species
		if (length(which(S>0))==1){ #do not rarefy when locale has less than 2 years with non-zero data
			print('Community has zero abundances')
			Srare = rep(NA,length(S))
		} else {
			raremax <- min(rowSums(comm.mat2))
			Srare <- rarefy(comm.mat2, raremax)
		}
	}
	for (y in 1:length(Srare)){ #create vector of rarefied S that matches order in out.diversity
		row.index = which(out.diversity$Locale==names(comm.list)[l] & out.diversity$Year==names(Srare)[y]) ###Varies by LTER
		add.Srare[row.index,1] = Srare[y]
	}
}
#
out.diversity = data.frame(out.diversity,add.Srare)
# Beta diversity
out.beta = data.frame('LTER.site'=NA,'Locale'=NA,'Year1'=-999,'Year2'=-999,'Beta.2'=-999,'Beta.j'=-999,'Beta.bray'=-999)
bx = 1
for (l in 1:length(comm.list)){
	comm.mat = comm.list[[l]]
	for (y in 1:(nrow(comm.mat)-1)){
		dat1 = comm.mat[c(y,y+1),]
		row.sums = rowSums(dat1)
		if (length(which(row.sums==0))>0 | length(which(is.na(row.sums)))>0){
			beta.2 = beta.j = vd.b = NA
		} else {
			beta.2 = betadiver(dat1,method="-2")	#presence/absence-based dissimilarity, according to Harrison et al. (1992)
			beta.j = 1-betadiver(dat1,method="j") #presence/absence-based, Jaccard (dis)similarity 
			vd.b = vegdist(dat1,method='bray') #abundance-based estimate of beta
		}
		out.beta[bx,1] = site.name #LTER
		out.beta[bx,2] = names(comm.list)[l] ###varies with LTER
		out.beta[bx,3] = rownames(comm.mat)[y] #year1
		out.beta[bx,4] = rownames(comm.mat)[y+1] #year2
		out.beta[bx,5] = beta.2 #beta.2
		out.beta[bx,6] = beta.j #beta.j
		out.beta[bx,7] = vd.b #beta.bray
		bx = bx + 1
	}
}
#Get values for species-centric table
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
		n.y1 = 1 #length(unique(paste(y.data$sample_date,y.data$trap_name,sep='_')))
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.diversity,paste0('./summary_tables/arthropods/Site-level/',site.name,'_ArthropodPitfall_TotalChange.txt'),sep='\t',quote=F,row.names=F)
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_ArthropodPitfall_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)
write.table(out.beta,paste0('./summary_tables/arthropods/Site-level/',site.name,'_ArthropodPitfall_BetaDiv.txt'),sep='\t',quote=F,row.names=F)


##### Central Arizona-Phoenix - Arthropod Pitfalls 2
site.name = 'CentralArizona-Phoenix'
locale.name = 'pitfall2'
data1 = read.csv("./raw_data/643_mcdowell_pitfall_arthropods_76003ef56b85e55df4e9f2605c3e8492.csv",as.is=T,check.names=F,header=T)
data1$Number = apply(data1[,13:17],1,function(x){sum(x,na.rm=T)})
data1$Year = apply(array(data1$sample_date),1,function(x){strsplit(x,'-')[[1]][1]})
data1$Species = data1$display_name
data1$Species = gsub('\\(adult\\)','',data1$Species)
data1$Species = gsub('\\(winged\\)','',data1$Species)
data1$Species = gsub('\\(larvae\\)','',data1$Species)
data1$Species = gsub('\\(nymph\\)','',data1$Species)
data1$Species = gsub(' \\(immature\\)','',data1$Species)
data1$Species = gsub(' \\(Winged Aphelinidae\\)','',data1$Species)
data1$Species = gsub(' \\(Aphelinidae\\)','',data1$Species)
data1$Species = gsub('_',' ',data1$Species)
data1$Species = gsub(' \\(Centipeds\\)','',data1$Species)
data1$Species = gsub('Scale Insects \\(Coccoidea\\)','',data1$Species)
data1$Species = gsub('Scale Insects \\(Coccoidea - Immature\\)','',data1$Species)
data1$Species = gsub(' \\(imm.\\)','',data1$Species)
data1$Species = gsub(' >10mm','',data1$Species)
data1$Species = gsub(' immature','',data1$Species)
data1$Species = trimws(data1$Species,which='right')
data1 = data1[which(data1$Species!=''),]
data1$Number = apply(data1[,13:17],1,function(x){sum(x,na.rm=T)})
data1 = data1[which(!is.na(data1$Species)),]
data1 = data1[which(data1$Species!='' & data1$Species!='Unknown'),]
u.species = unique(data1$Species)
u.years = sort(as.numeric(unique(data1$Year)))
#Get values for LTER site-centric table
out.diversity = data.frame('LTER.site'=NA,'Locale'=NA,'Year'=-999,'N.obs'=-999,'Total.abundance'=-999,'N.species'=-999,'Species.evenness'=-999,'Species.decay.rate'=-999,'Fishers.alpha'=-999,'Dominance'=-999)
lx = 1
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[y]),]
	n.obs = 1 #length(unique(paste(y.data$sample_date,y.data$trap_name,sep='_')))
	divs = calc.diversity(y.data$Species,y.data$Number)
	out.diversity[lx,1] = site.name #LTER
	out.diversity[lx,2] = locale.name #Locale
	out.diversity[lx,3] = u.years[y] #Year
	out.diversity[lx,4] = n.obs #Number of observations
	out.diversity[lx,5] = sum(y.data$Number,na.rm=T) #Total abundance
	out.diversity[lx,6] = divs[1] #Total number of species
	out.diversity[lx,7] = divs[2] #Evenness
	out.diversity[lx,8] = divs[3] #rank species abundance decay rate
	out.diversity[lx,9] = divs[4] #Fisher's alpha (measure of rarity)
	out.diversity[lx,10] = divs[5] #Dominance
	lx = lx + 1
}
# Rarefied alpha and beta diversity
u.years = sort(as.numeric(unique(data1$Year)))
comm.list = list()
clx = 1
comm.mat = matrix(NA,nrow=length(u.years),ncol=length(u.species))
rownames(comm.mat) = u.years
colnames(comm.mat) = u.species
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		row.pos = which(data1$Year==u.years[y] & data1$Species==u.species[s])
		if (length(row.pos)==0){
			comm.mat[y,s] = NA
		} else {
			y.dat = data1[row.pos,]
			comm.mat[y,s] = sum(y.dat$Number,na.rm=T)
		}
	}
}
comm.list[[clx]] = comm.mat
clx = clx + 1
names(comm.list) = locale.name
for (l in 1:length(comm.list)){ #replace NA with zero
	for (i in 1:nrow(comm.list[[l]])){
		comm.list[[l]][i,][which(is.na(comm.list[[l]][i,]))] = 0
	}
}
# Alpha diversity
add.Srare = data.frame('N.species.rarefied'=rep(NA,nrow(out.diversity)))
for (l in 1:length(comm.list)){
	print(names(comm.list)[l])
	comm.mat = comm.list[[l]]
	S <- specnumber(comm.mat) # observed number of species
	rkeep = which(!is.na(S))
	S = S[rkeep]
	comm.mat = comm.mat[rkeep,]
	if (sum(S)==0){ #do not rarefy communities with zero abundances for all years
		print('Community has all zero abundances')
		Srare = rep(NA,length(S)) 
	} else {
		comm.mat2 = comm.mat[which(S>0),] #remove years with 0 species
		if (length(which(S>0))==1){ #do not rarefy when locale has less than 2 years with non-zero data
			print('Community has zero abundances')
			Srare = rep(NA,length(S))
		} else {
			raremax <- min(rowSums(comm.mat2))
			Srare <- rarefy(comm.mat2, raremax)
		}
	}
	for (y in 1:length(Srare)){ #create vector of rarefied S that matches order in out.diversity
		row.index = which(out.diversity$Locale==names(comm.list)[l] & out.diversity$Year==names(Srare)[y]) ###Varies by LTER
		add.Srare[row.index,1] = Srare[y]
	}
}
#
out.diversity = data.frame(out.diversity,add.Srare)
# Beta diversity
out.beta = data.frame('LTER.site'=NA,'Locale'=NA,'Year1'=-999,'Year2'=-999,'Beta.2'=-999,'Beta.j'=-999,'Beta.bray'=-999)
bx = 1
for (l in 1:length(comm.list)){
	comm.mat = comm.list[[l]]
	for (y in 1:(nrow(comm.mat)-1)){
		dat1 = comm.mat[c(y,y+1),]
		row.sums = rowSums(dat1)
		if (length(which(row.sums==0))>0 | length(which(is.na(row.sums)))>0){
			beta.2 = beta.j = vd.b = NA
		} else {
			beta.2 = betadiver(dat1,method="-2")	#presence/absence-based dissimilarity, according to Harrison et al. (1992)
			beta.j = 1-betadiver(dat1,method="j") #presence/absence-based, Jaccard (dis)similarity 
			vd.b = vegdist(dat1,method='bray') #abundance-based estimate of beta
		}
		out.beta[bx,1] = site.name #LTER
		out.beta[bx,2] = names(comm.list)[l] ###varies with LTER
		out.beta[bx,3] = rownames(comm.mat)[y] #year1
		out.beta[bx,4] = rownames(comm.mat)[y+1] #year2
		out.beta[bx,5] = beta.2 #beta.2
		out.beta[bx,6] = beta.j #beta.j
		out.beta[bx,7] = vd.b #beta.bray
		bx = bx + 1
	}
}
#Get values for species-centric table
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
		n.y1 = 1 #length(unique(paste(y.data$sample_date,y.data$trap_name,sep='_')))
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.diversity,paste0('./summary_tables/arthropods/Site-level/',site.name,'_ArthropodPitfall2_TotalChange.txt'),sep='\t',quote=F,row.names=F)
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_ArthropodPitfall2_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)
write.table(out.beta,paste0('./summary_tables/arthropods/Site-level/',site.name,'_ArthropodPitfall2_BetaDiv.txt'),sep='\t',quote=F,row.names=F)


##### Sevilleta - Arthropod pitfalls
site.name = 'Sevilleta'
locale.name = 'pitfall'
data1 = read.csv("./raw_data/sev029_arthropop_02162009_0.txt",as.is=T,check.names=F,header=T)
for (i in 1:ncol(data1)){
	data1[which(data1[,i]==-888),i] = NA
}
data1 = data1[which(data1$Year>1994),]
data1$Species = paste(data1$Genus,data1$Species,sep='_')
data1$Number = data1$Count
data1 = data1[which(data1$Species!='_' & data1$Species!='NA_NA'),]
u.years = sort(as.numeric(unique(data1$Year)))
u.species = unique(data1$Species)
#Get values for LTER site-centric table
out.diversity = data.frame('LTER.site'=NA,'Locale'=NA,'Year'=-999,'N.obs'=-999,'Total.abundance'=-999,'N.species'=-999,'Species.evenness'=-999,'Species.decay.rate'=-999,'Fishers.alpha'=-999,'Dominance'=-999)
lx = 1
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[y]),]
	n.obs = 1 
	divs = calc.diversity(y.data$Species,y.data$Number)
	out.diversity[lx,1] = site.name #LTER
	out.diversity[lx,2] = locale.name #Locale
	out.diversity[lx,3] = u.years[y] #Year
	out.diversity[lx,4] = n.obs #Number of observations
	out.diversity[lx,5] = sum(y.data$Number,na.rm=T) #Total abundance
	out.diversity[lx,6] = divs[1] #Total number of species
	out.diversity[lx,7] = divs[2] #Evenness
	out.diversity[lx,8] = divs[3] #rank species abundance decay rate
	out.diversity[lx,9] = divs[4] #Fisher's alpha (measure of rarity)
	out.diversity[lx,10] = divs[5] #Dominance
	lx = lx + 1
}
# Rarefied alpha and beta diversity
u.years = sort(as.numeric(unique(data1$Year)))
comm.list = list()
clx = 1
comm.mat = matrix(NA,nrow=length(u.years),ncol=length(u.species))
rownames(comm.mat) = u.years
colnames(comm.mat) = u.species
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		row.pos = which(data1$Year==u.years[y] & data1$Species==u.species[s])
		if (length(row.pos)==0){
			comm.mat[y,s] = NA
		} else {
			y.dat = data1[row.pos,]
			comm.mat[y,s] = sum(y.dat$Number,na.rm=T)
		}
	}
}
comm.list[[clx]] = comm.mat
clx = clx + 1
names(comm.list) = locale.name
for (l in 1:length(comm.list)){ #replace NA with zero
	for (i in 1:nrow(comm.list[[l]])){
		comm.list[[l]][i,][which(is.na(comm.list[[l]][i,]))] = 0
	}
}
# Alpha diversity
add.Srare = data.frame('N.species.rarefied'=rep(NA,nrow(out.diversity)))
for (l in 1:length(comm.list)){
	print(names(comm.list)[l])
	comm.mat = comm.list[[l]]
	S <- specnumber(comm.mat) # observed number of species
	rkeep = which(!is.na(S))
	S = S[rkeep]
	comm.mat = comm.mat[rkeep,]
	if (sum(S)==0){ #do not rarefy communities with zero abundances for all years
		print('Community has all zero abundances')
		Srare = rep(NA,length(S)) 
	} else {
		comm.mat2 = comm.mat[which(S>0),] #remove years with 0 species
		if (length(which(S>0))==1){ #do not rarefy when locale has less than 2 years with non-zero data
			print('Community has zero abundances')
			Srare = rep(NA,length(S))
		} else {
			raremax <- min(rowSums(comm.mat2))
			Srare <- rarefy(comm.mat2, raremax)
		}
	}
	for (y in 1:length(Srare)){ #create vector of rarefied S that matches order in out.diversity
		row.index = which(out.diversity$Locale==names(comm.list)[l] & out.diversity$Year==names(Srare)[y]) ###Varies by LTER
		add.Srare[row.index,1] = Srare[y]
	}
}
#
out.diversity = data.frame(out.diversity,add.Srare)
# Beta diversity
out.beta = data.frame('LTER.site'=NA,'Locale'=NA,'Year1'=-999,'Year2'=-999,'Beta.2'=-999,'Beta.j'=-999,'Beta.bray'=-999)
bx = 1
for (l in 1:length(comm.list)){
	comm.mat = comm.list[[l]]
	for (y in 1:(nrow(comm.mat)-1)){
		dat1 = comm.mat[c(y,y+1),]
		row.sums = rowSums(dat1)
		if (length(which(row.sums==0))>0 | length(which(is.na(row.sums)))>0){
			beta.2 = beta.j = vd.b = NA
		} else {
			beta.2 = betadiver(dat1,method="-2")	#presence/absence-based dissimilarity, according to Harrison et al. (1992)
			beta.j = 1-betadiver(dat1,method="j") #presence/absence-based, Jaccard (dis)similarity 
			vd.b = vegdist(dat1,method='bray') #abundance-based estimate of beta
		}
		out.beta[bx,1] = site.name #LTER
		out.beta[bx,2] = names(comm.list)[l]
		out.beta[bx,3] = rownames(comm.mat)[y] #year1
		out.beta[bx,4] = rownames(comm.mat)[y+1] #year2
		out.beta[bx,5] = beta.2 #beta.2
		out.beta[bx,6] = beta.j #beta.j
		out.beta[bx,7] = vd.b #beta.bray
		bx = bx + 1
	}
}
#Get values for species-centric table
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
		n.y1 = 1
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.diversity,paste0('./summary_tables/arthropods/Site-level/',site.name,'_ArthropodPitfall_TotalChange.txt'),sep='\t',quote=F,row.names=F)
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_ArthropodPitfall_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)
write.table(out.beta,paste0('./summary_tables/arthropods/Site-level/',site.name,'_ArthropodPitfall_BetaDiv.txt'),sep='\t',quote=F,row.names=F)



##### Georgia Coastal Ecosystems - Invertebrates
# NOTE: data on 'grasshoppers' and Prokelisia marginata (a plant hopper)
site.name = 'GeorgiaCoastal'
locale.name = 'Salt Marsh'
data1 = read.table("./raw_data/PLT-GCEM-1610_Animals_3_0.TXT",sep='\t',as.is=T,check.names=F,header=T)
data1$Year = apply(array(data1$Date),1,function(x){strsplit(x,'-')[[1]][1]})
u.years = unique(data1$Year)
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[y]),]
	n.y1 = 1
	out.slopes[cx,1] = site.name #LTER.site
	out.slopes[cx,2] = locale.name #Locale
	out.slopes[cx,3] = 'Acrididae' #Species.code
	out.slopes[cx,4] = u.years[y] #Year
	out.slopes[cx,5] = n.y1 #Number of observations
	out.slopes[cx,6] = sum(y.data$Grasshopper_Abundance_Index,na.rm=T) #abundance
	cx = cx + 1
	out.slopes[cx,1] = site.name #LTER.site
	out.slopes[cx,2] = locale.name #Locale
	out.slopes[cx,3] = 'Prokelisia marginata' #Species.code
	out.slopes[cx,4] = u.years[y] #Year
	out.slopes[cx,5] = n.y1 #Number of observations
	out.slopes[cx,6] = sum(y.data$Prokelisia_Abundance_Index,na.rm=T) #abundance
	cx = cx + 1
}
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_Invertebrate_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)


##### Hubbard Brook - Lepidoptera 1
site.name = 'HubbardBrook'
locale.name = 'Lepidoptera1'
data1 = read.csv("./raw_data/lepidoptera_HB_1986-2018.csv",as.is=T,check.names=F,header=T)
data1$Number = data1$NumberIndividuals
data1$Species = data1$Taxon
data1 = data1[-which(data1$Species==0),]
data1 = data1[-which(is.na(data1$Species)),]
u.years = sort(as.numeric(unique(data1$Year)))
u.species = unique(data1$Species)
#Get values for species-centric table
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
		n.y1 = 1
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_Lepidoptera1_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)


##### Hubbard Brook - Lepidoptera 2
site.name = 'HubbardBrook'
locale.name = 'Lepidoptera2'
data1 = read.csv("./raw_data/lepidoptera_WMNFsites.csv",as.is=T,check.names=F,header=T)
data1$Number = data1$NumberIndividuals
data1$Species = data1$Taxon
u.years = sort(as.numeric(unique(data1$Year)))
u.species = unique(data1$Species)
#Get values for species-centric table
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
		n.y1 = 1
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_Lepidoptera2_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)



##### Baltimore
# NOTE: species counts sometimes occur as decimals; I suspect that instars were counted as fractions; I rounded-up to get integer counts
site.name = 'Baltimore'
locale.name = 'mosquito'
data1 = read.csv( "./raw_data/Biodiversity_-_Mosquito_-_ovitrap_mosquito_data.csv",as.is=T,check.names=F,header=T)
data1$Year = data1$year
u.years = unique(data1$Year)
u.species = colnames(data1)[6:17]
for (i in 1:nrow(data1)){
	for (j in 6:17){
		data1[i,j] = ceiling(data1[i,j])
	}
}
#Get values for LTER site-centric table
out.diversity = data.frame('LTER.site'=NA,'Locale'=NA,'Year'=-999,'N.obs'=-999,'Total.abundance'=-999,'N.species'=-999,'Species.evenness'=-999,'Species.decay.rate'=-999,'Fishers.alpha'=-999,'Dominance'=-999)
lx = 1
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[y]),]
	n.obs = length(unique(paste(y.data$site,y.data$week,sep='_'))) #site x week samples
	divs = calc.diversity(colnames(y.data)[6:17],apply(y.data[,6:17],2,function(x){sum(x,na.rm=T)}))
	out.diversity[lx,1] = site.name #LTER
	out.diversity[lx,2] = locale.name #Locale
	out.diversity[lx,3] = u.years[y] #Year
	out.diversity[lx,4] = n.obs #Number of observations
	out.diversity[lx,5] = sum(unlist(y.data[,6:17]),na.rm=T) #Total abundance
	out.diversity[lx,6] = divs[1] #Total number of species
	out.diversity[lx,7] = divs[2] #Evenness
	out.diversity[lx,8] = divs[3] #rank species abundance decay rate
	out.diversity[lx,9] = divs[4] #Fisher's alpha (measure of rarity)
	out.diversity[lx,10] = divs[5] #Dominance
	lx = lx + 1
}
# Rarefied alpha and beta diversity
u.years = sort(as.numeric(unique(data1$Year)))
comm.list = list()
clx = 1
comm.mat = matrix(NA,nrow=length(u.years),ncol=length(u.species))
rownames(comm.mat) = u.years
colnames(comm.mat) = u.species
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		row.pos = which(data1$Year==u.years[y])
		if (length(row.pos)==0){
			comm.mat[y,s] = NA
		} else {
			y.dat = data1[row.pos,which(colnames(data1)==u.species[s])]
			comm.mat[y,s] = sum(y.dat,na.rm=T)
		}
	}
}
comm.list[[clx]] = comm.mat
clx = clx + 1
names(comm.list) = locale.name
for (l in 1:length(comm.list)){ #replace NA with zero
	for (i in 1:nrow(comm.list[[l]])){
		comm.list[[l]][i,][which(is.na(comm.list[[l]][i,]))] = 0
	}
}
# Alpha diversity
add.Srare = data.frame('N.species.rarefied'=rep(NA,nrow(out.diversity)))
for (l in 1:length(comm.list)){
	print(names(comm.list)[l])
	comm.mat = comm.list[[l]]
	S <- specnumber(comm.mat) # observed number of species
	rkeep = which(!is.na(S))
	S = S[rkeep]
	comm.mat = comm.mat[rkeep,]
	if (sum(S)==0){ #do not rarefy communities with zero abundances for all years
		print('Community has all zero abundances')
		Srare = rep(NA,length(S)) 
	} else {
		comm.mat2 = comm.mat[which(S>0),] #remove years with 0 species
		if (length(which(S>0))==1){ #do not rarefy when locale has less than 2 years with non-zero data
			print('Community has zero abundances')
			Srare = rep(NA,length(S))
		} else {
			raremax <- min(rowSums(comm.mat2))
			Srare <- rarefy(comm.mat2, raremax)
		}
	}
	for (y in 1:length(Srare)){ #create vector of rarefied S that matches order in out.diversity
		row.index = which(out.diversity$Locale==names(comm.list)[l] & out.diversity$Year==names(Srare)[y]) ###Varies by LTER
		add.Srare[row.index,1] = Srare[y]
	}
}
#
out.diversity = data.frame(out.diversity,add.Srare)
# Beta diversity
out.beta = data.frame('LTER.site'=NA,'Locale'=NA,'Year1'=-999,'Year2'=-999,'Beta.2'=-999,'Beta.j'=-999,'Beta.bray'=-999)
bx = 1
for (l in 1:length(comm.list)){
	comm.mat = comm.list[[l]]
	for (y in 1:(nrow(comm.mat)-1)){
		dat1 = comm.mat[c(y,y+1),]
		row.sums = rowSums(dat1)
		if (length(which(row.sums==0))>0 | length(which(is.na(row.sums)))>0){
			beta.2 = beta.j = vd.b = NA
		} else {
			beta.2 = betadiver(dat1,method="-2")	#presence/absence-based dissimilarity, according to Harrison et al. (1992)
			beta.j = 1-betadiver(dat1,method="j") #presence/absence-based, Jaccard (dis)similarity 
			vd.b = vegdist(dat1,method='bray') #abundance-based estimate of beta
		}
		out.beta[bx,1] = site.name #LTER
		out.beta[bx,2] = names(comm.list)[l]
		out.beta[bx,3] = rownames(comm.mat)[y] #year1
		out.beta[bx,4] = rownames(comm.mat)[y+1] #year2
		out.beta[bx,5] = beta.2 #beta.2
		out.beta[bx,6] = beta.j #beta.j
		out.beta[bx,7] = vd.b #beta.bray
		bx = bx + 1
	}
}
#Get values for species-centric table
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y]),c(1,3,28,which(colnames(data1)==u.species[s]))]
		n.y1 = length(unique(paste(y.data$site,y.data$week,sep='_'))) #site x week samples
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(as.numeric(y.data[,4]),na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.diversity,paste0('./summary_tables/arthropods/Site-level/',site.name,'_Mosquito_TotalChange.txt'),sep='\t',quote=F,row.names=F)
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_Mosquito_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)
write.table(out.beta,paste0('./summary_tables/arthropods/Site-level/',site.name,'_Mosquito_BetaDiv.txt'),sep='\t',quote=F,row.names=F)


##### Arctic Tundra - Stream invertebrates
site.name = 'Arctic'
locale.name = 'aquatic'
data1 = read.csv("./raw_data/84-98hektot.csv",as.is=T,check.names=F,header=T)
data1$Year = apply(array(data1$Date),1,function(x){paste0(19,strsplit(x,'-')[[1]][3])})
data1 = data1[,which(colnames(data1)!='SNAILS')]
u.years = sort(as.numeric(unique(data1$Year)))
u.species = colnames(data1)[7:20]
#Get values for LTER site-centric table
out.diversity = data.frame('LTER.site'=NA,'Locale'=NA,'Year'=-999,'N.obs'=-999,'Total.abundance'=-999,'N.species'=-999,'Species.evenness'=-999,'Species.decay.rate'=-999,'Fishers.alpha'=-999,'Dominance'=-999)
lx = 1
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[y]),]
	n.obs = 1
	divs = calc.diversity(colnames(y.data)[7:20],apply(y.data[,7:20],2,function(x){sum(x,na.rm=T)}))
	out.diversity[lx,1] = site.name #LTER
	out.diversity[lx,2] = locale.name #Locale
	out.diversity[lx,3] = u.years[y] #Year
	out.diversity[lx,4] = n.obs #Number of observations
	out.diversity[lx,5] = sum(unlist(y.data[,7:20]),na.rm=T) #Total abundance
	out.diversity[lx,6] = divs[1] #Total number of species
	out.diversity[lx,7] = divs[2] #Evenness
	out.diversity[lx,8] = divs[3] #rank species abundance decay rate
	out.diversity[lx,9] = divs[4] #Fisher's alpha (measure of rarity)
	out.diversity[lx,10] = divs[5] #Dominance
	lx = lx + 1
}
# Rarefied alpha and beta diversity
u.years = sort(as.numeric(unique(data1$Year)))
comm.list = list()
clx = 1
comm.mat = matrix(NA,nrow=length(u.years),ncol=length(u.species))
rownames(comm.mat) = u.years
colnames(comm.mat) = u.species
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		row.pos = which(data1$Year==u.years[y])
		if (length(row.pos)==0){
			comm.mat[y,s] = NA
		} else {
			y.dat = data1[row.pos,which(colnames(data1)==u.species[s])]
			comm.mat[y,s] = sum(y.dat,na.rm=T)
		}
	}
}
comm.list[[clx]] = comm.mat
clx = clx + 1
names(comm.list) = locale.name
for (l in 1:length(comm.list)){ #replace NA with zero
	for (i in 1:nrow(comm.list[[l]])){
		comm.list[[l]][i,][which(is.na(comm.list[[l]][i,]))] = 0
	}
}
# Alpha diversity
add.Srare = data.frame('N.species.rarefied'=rep(NA,nrow(out.diversity)))
for (l in 1:length(comm.list)){
	print(names(comm.list)[l])
	comm.mat = comm.list[[l]]
	S <- specnumber(comm.mat) # observed number of species
	rkeep = which(!is.na(S))
	S = S[rkeep]
	comm.mat = comm.mat[rkeep,]
	if (sum(S)==0){ #do not rarefy communities with zero abundances for all years
		print('Community has all zero abundances')
		Srare = rep(NA,length(S)) 
	} else {
		comm.mat2 = comm.mat[which(S>0),] #remove years with 0 species
		if (length(which(S>0))==1){ #do not rarefy when locale has less than 2 years with non-zero data
			print('Community has zero abundances')
			Srare = rep(NA,length(S))
		} else {
			raremax <- min(rowSums(comm.mat2))
			Srare <- rarefy(comm.mat2, raremax)
		}
	}
	for (y in 1:length(Srare)){ #create vector of rarefied S that matches order in out.diversity
		row.index = which(out.diversity$Locale==names(comm.list)[l] & out.diversity$Year==names(Srare)[y]) ###Varies by LTER
		add.Srare[row.index,1] = Srare[y]
	}
}
#
out.diversity = data.frame(out.diversity,add.Srare)
# Beta diversity
out.beta = data.frame('LTER.site'=NA,'Locale'=NA,'Year1'=-999,'Year2'=-999,'Beta.2'=-999,'Beta.j'=-999,'Beta.bray'=-999)
bx = 1
for (l in 1:length(comm.list)){
	comm.mat = comm.list[[l]]
	for (y in 1:(nrow(comm.mat)-1)){
		dat1 = comm.mat[c(y,y+1),]
		row.sums = rowSums(dat1)
		if (length(which(row.sums==0))>0 | length(which(is.na(row.sums)))>0){
			beta.2 = beta.j = vd.b = NA
		} else {
			beta.2 = betadiver(dat1,method="-2")	#presence/absence-based dissimilarity, according to Harrison et al. (1992)
			beta.j = 1-betadiver(dat1,method="j") #presence/absence-based, Jaccard (dis)similarity 
			vd.b = vegdist(dat1,method='bray') #abundance-based estimate of beta
		}
		out.beta[bx,1] = site.name #LTER
		out.beta[bx,2] = names(comm.list)[l]
		out.beta[bx,3] = rownames(comm.mat)[y] #year1
		out.beta[bx,4] = rownames(comm.mat)[y+1] #year2
		out.beta[bx,5] = beta.2 #beta.2
		out.beta[bx,6] = beta.j #beta.j
		out.beta[bx,7] = vd.b #beta.bray
		bx = bx + 1
	}
}
#Get values for species-centric table
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y]),c(1:6,21,which(colnames(data1)==u.species[s]))]
		n.y1 = 1
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(as.numeric(y.data[,8]),na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.diversity,paste0('./summary_tables/arthropods/Site-level/',site.name,'_Aquatic_TotalChange.txt'),sep='\t',quote=F,row.names=F)
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_Aquatic_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)
write.table(out.beta,paste0('./summary_tables/arthropods/Site-level/',site.name,'_Aquatic_BetaDiv.txt'),sep='\t',quote=F,row.names=F)


##### Harvard Forest - Ticks
site.name = 'HarvardForest'
data1 = read.csv("./raw_data/hf299-01-survey.csv",as.is=T,check.names=F,header=T)
data1$Year = data1$year
data1$Number = apply(data1[,8:9],1,function(x){sum(x,na.rm=T)})
locales = gsub(' ','',gsub("'",'',data1$location.name))
u.locales = unique(locales)
u.locales = u.locales[-which(is.na(u.locales))]
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
sps = c(0,0,0,0,0,0,0,'wood.tick','deer.tick')
for (l in 1:length(u.locales)){
	l.data = data1[which(locales==u.locales[l]),]
	for (s in 8:9){
		s.data = data.frame('Year'=l.data$Year,'hours'=l.data$hours,'Number'=l.data[,s])
		u.years = sort(as.numeric(unique(s.data$Year)))
		for (y in 1:length(u.years)){
			y.data = s.data[which(s.data$Year==u.years[y]),]
			n.y1 = sum(y.data$hours,na.rm=T)
			out.slopes[cx,1] = site.name #LTER.site
			out.slopes[cx,2] = u.locales[l] #Locale
			out.slopes[cx,3] = sps[s] #Species.code
			out.slopes[cx,4] = u.years[y] #Year
			out.slopes[cx,5] = n.y1 #Number of observations
			out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
			cx = cx + 1
		}
	}
}
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_Tick_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)


##### Northern Temperate Lakes - Arthropods/Crustaceans
# Sampling effort accounted for by standardizing counts by volume sampled
# NOTE: only 4 species documented; skipping diversity stats
site.name = 'NorthTermperateLakes'
data1 = read.csv("./raw_data/ntl14_v8.csv",as.is=T,check.names=F,header=T)
data1$Year = data1$year4
data1$Site = data1$lakeid
data1$Species = data1$taxon
data1$Number = ceiling(data1$avg_ind_per_m3)
u.locales = unique(data1$Site)
u.years = sort(as.numeric(unique(data1$Year)))
u.species = unique(data1$Species)
#Get values for species-centric table
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (l in 1:length(u.locales)){
	for (s in 1:length(u.species)){
		for (y in 1:length(u.years)){
			y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s] & data1$Site==u.locales[l]),]
			n.y1 = length(unique(paste(y.data$sta,y.data$depth)))
			out.slopes[cx,1] = site.name #LTER.site
			out.slopes[cx,2] = u.locales[l] #Locale
			out.slopes[cx,3] = u.species[s] #Species.code
			out.slopes[cx,4] = u.years[y] #Year
			out.slopes[cx,5] = n.y1 #Number of observations
			out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
			cx = cx + 1
		}
	}
}
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_Arthropod_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)



##### Northern Temperate Lakes - Invertebrates 1
site.name = 'NorthTemperateLakes'
data1 = read.csv("./raw_data/ntl11_1_v8.csv",as.is=T,check.names=F,header=T)
data1$Year = data1$year4
data1$Site = data1$lakeid
data1$Species = data1$taxon_code
data1$Number = data1$number_indiv
u.locales = unique(data1$Site)
u.years = sort(as.numeric(unique(data1$Year)))
u.species = unique(data1$Species)
#Get values for LTER site-centric table
out.diversity = data.frame('LTER.site'=NA,'Locale'=NA,'Year'=-999,'N.obs'=-999,'Total.abundance'=-999,'N.species'=-999,'Species.evenness'=-999,'Species.decay.rate'=-999,'Fishers.alpha'=-999,'Dominance'=-999)
lx = 1
for (l in 1:length(u.locales)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Site==u.locales[l]),]
		n.obs = length(unique(paste(y.data$site,y.data$rep)))
		divs = calc.diversity(y.data$Species,y.data$Number)
		out.diversity[lx,1] = site.name #LTER
		out.diversity[lx,2] = u.locales[l] #Locale
		out.diversity[lx,3] = u.years[y] #Year
		out.diversity[lx,4] = n.obs #Number of observations
		out.diversity[lx,5] = sum(y.data$Number,na.rm=T) #Total abundance
		out.diversity[lx,6] = divs[1] #Total number of species
		out.diversity[lx,7] = divs[2] #Evenness
		out.diversity[lx,8] = divs[3] #rank species abundance decay rate
		out.diversity[lx,9] = divs[4] #Fisher's alpha (measure of rarity)
		out.diversity[lx,10] = divs[5] #Dominance
		lx = lx + 1
	}
}
# Rarefied alpha and beta diversity
comm.list = list()
clx = 1
for (l in 1:length(u.locales)){ #create matrices of species abundance over time, for each "site" and/or "observation method"
	comm.mat = matrix(NA,nrow=length(u.years),ncol=length(u.species))
	rownames(comm.mat) = u.years
	colnames(comm.mat) = u.species
	for (s in 1:length(u.species)){
		for (y in 1:length(u.years)){
			y.dat = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s] & data1$Site==u.locales[l]),]
			comm.mat[y,s] = sum(y.dat$Number,na.rm=T)
		}
	}
	comm.list[[clx]] = comm.mat
	clx = clx + 1
}
names(comm.list) = u.locales
# Alpha diversity
add.Srare = data.frame('N.species.rarefied'=rep(NA,nrow(out.diversity)))
for (l in 1:length(comm.list)){
	comm.mat = comm.list[[l]]
	S <- specnumber(comm.mat) # observed number of species
	comm.mat2 = comm.mat[which(S>0),] #remove years with 0 species
	raremax <- min(rowSums(comm.mat2))
	Srare <- rarefy(comm.mat2, raremax)
	for (y in 1:length(Srare)){ #create vector of rarefied S that matches order in out.diversity
		row.index = which(out.diversity$Locale==names(comm.list)[l] & out.diversity$Year==names(Srare)[y]) ###Varies by LTER
		add.Srare[row.index,1] = Srare[y]
	}
}
out.diversity = data.frame(out.diversity,add.Srare)
# Beta diversity
out.beta = data.frame('LTER.site'=NA,'Locale'=NA,'Year1'=-999,'Year2'=-999,'Beta.2'=-999,'Beta.j'=-999,'Beta.bray'=-999)
bx = 1
for (l in 1:length(comm.list)){
	comm.mat = comm.list[[l]]
	for (y in 1:(nrow(comm.mat)-1)){
		dat1 = comm.mat[c(y,y+1),]
		row.sums = rowSums(dat1)
		if (length(which(row.sums==0))>0){
			beta.2 = beta.j = vd.b = NA
		} else {
			beta.2 = betadiver(dat1,method="-2")	#presence/absence-based dissimilarity, according to Harrison et al. (1992)
			beta.j = 1-betadiver(dat1,method="j") #presence/absence-based, Jaccard (dis)similarity 
			vd.b = vegdist(dat1,method='bray') #abundance-based estimate of beta
		}
		out.beta[bx,1] = site.name #LTER
		out.beta[bx,2] = names(comm.list)[l] #locale
		out.beta[bx,3] = rownames(comm.mat)[y] #year1
		out.beta[bx,4] = rownames(comm.mat)[y+1] #year2
		out.beta[bx,5] = beta.2 #beta.2
		out.beta[bx,6] = beta.j #beta.j
		out.beta[bx,7] = vd.b #beta.bray
		bx = bx + 1
	}
}
out.beta2 = out.beta[which(!is.na(out.beta$Beta.2)),]
plot(out.beta2$Beta.2,out.beta2$Beta.j)
#Get values for species-centric table
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (l in 1:length(u.locales)){
	for (s in 1:length(u.species)){
		for (y in 1:length(u.years)){
			y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s] & data1$Site==u.locales[l]),]
			n.y1 = length(unique(paste(y.data$site,y.data$rep)))
			out.slopes[cx,1] = site.name #LTER.site
			out.slopes[cx,2] = u.locales[l] #Locale
			out.slopes[cx,3] = u.species[s] #Species.code
			out.slopes[cx,4] = u.years[y] #Year
			out.slopes[cx,5] = n.y1 #Number of observations
			out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
			cx = cx + 1
		}
	}
}
write.table(out.diversity,paste0('./summary_tables/arthropods/Site-level/',site.name,'_Invertebrate1_TotalChange.txt'),sep='\t',quote=F,row.names=F)
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_Invertebrate1_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)
write.table(out.beta,paste0('./summary_tables/arthropods/Site-level/',site.name,'_Invertebrate1_BetaDiv.txt'),sep='\t',quote=F,row.names=F)



##### Georgia Coastal Ecosystems
site.name = 'GeorgiaCoastal'
locale.name = 'Burrowing Crab'
data1 = read.table("./raw_data/INV-GCES-1609_1_0.TXT",sep='\t',as.is=T,check.names=F,header=T)
u.years = sort(as.numeric(unique(data1$Year)))
#Get values for species-centric table
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1	
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[y]),]
	n.y1 = length(unique(paste(y.data$Site,y.data$Zone,y.data$Plot)))
	out.slopes[cx,1] = site.name #LTER.site
	out.slopes[cx,2] = locale.name #Locale
	out.slopes[cx,3] = 'burrowing crab' #Species.code
	out.slopes[cx,4] = u.years[y] #Year
	out.slopes[cx,5] = n.y1 #Number of observations
	out.slopes[cx,6] = sum(as.numeric(y.data$Hole_Density),na.rm=T) #abundance
	cx = cx + 1
}
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_BurrowingCrab_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)



##### Temperate Lakes - Crayfish 1
site.name = 'NorthTemperateLakes'
data1 = read.csv("./raw_data/ntl3_v7.csv",as.is=T,check.names=F,header=T)
data1$Species = data1$spname
data1$Site = data1$lakeid
data1$Year = data1$year4
data1$Number = data1$total_caught
data1 = data1[which(data1$Species!='CRAYFISH'),]
u.locales = unique(data1$Site)
u.years = sort(as.numeric(unique(data1$Year)))
u.species = unique(data1$Species)
#Get values for species-centric table
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (l in 1:length(u.locales)){
	for (s in 1:length(u.species)){
		for (y in 1:length(u.years)){
			y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s] & data1$Site==u.locales[l]),]
			n.y1 = sum(y.data$effort,na.rm=T)
			out.slopes[cx,1] = site.name #LTER.site
			out.slopes[cx,2] = u.locales[l] #Locale
			out.slopes[cx,3] = u.species[s] #Species.code
			out.slopes[cx,4] = u.years[y] #Year
			out.slopes[cx,5] = n.y1 #Number of observations
			out.slopes[cx,6] = sum(as.numeric(y.data$Hole_Density),na.rm=T) #abundance
			cx = cx + 1
		}
	}
}	
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_Crayfish1_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)



##### Temperate Lakes - Crayfish 2
site.name = 'NorthTemperateLakes'
locale.name = 'Crayfish'
data1 = read.csv("./raw_data/sparkling_crayfish.csv",as.is=T,check.names=F,header=T)
data1$Year = data1$YEAR4
u.years = sort(as.numeric(unique(data1$Year)))
u.species = colnames(data1)[7:8]
#Get values for species-centric table
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y]),c(6,11,which(colnames(data1)==u.species[s]))]
		n.y1 = sum(y.data$TRAP_DAYS,na.rm=T)
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data[,3],na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_Crayfish2_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)


##### Georgia Coastal Ecosystems - Fiddler crabs
site.name = 'GeorgiaCoastal'
locale.name = 'Fiddler Crab'
data1 = read.table("./raw_data/INV-GCED-1406_1_0.TXT",sep='\t',as.is=T,check.names=F,header=T)
data1$Number = (data1$Adult_Fiddler_crab_burrow_count / data1$Quadrat_Area_Adult_Fiddler_Crab_Burrow_Count) + (data1$Juvenile_Fiddler_crab_burrow_count / data1$Quadrat_Area_Juvenile_Fiddler_Crab_Burrow_Count)
u.years = unique(data1$Year)
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[y]),]
	n.y1 = length(unique(paste(y.data$Month,y.data$Day,y.data$Experimental_Treatment,y.data$Plot_replicate)))
	out.slopes[cx,1] = site.name #LTER.site
	out.slopes[cx,2] = locale.name #Locale
	out.slopes[cx,3] = 'fiddler crab' #Species.code
	out.slopes[cx,4] = u.years[y] #Year
	out.slopes[cx,5] = n.y1 #Number of observations
	out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
	cx = cx + 1
}
write.table(out.slopes,paste0('./summary_tables/arthropods/Species-level/',site.name,'_FiddlerCrab_SpeciesSlopes.txt'),sep='\t',quote=F,row.names=F)


##########################################################Section
# Merge slopes files

files = list.files('./summary_tables/arthropods/Species-level',full.names=T)

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


#################################################################Section
# Classify species according to habitat and feeding guild

# Import species annotations
eco.key = read.csv('./taxa_keys/Ecologial Variables.Arthropods.Final_v2.csv',as.is=T,check.names=F,header=T)
colnames(eco.key)[1] = 'Code'
eco.key$Species = trimws(eco.key$Species,which='right')
eco.key$Family = trimws(eco.key$Family,which='right')
eco.key$Order = trimws(eco.key$Order,which='right')

# Import species count data
output = read.csv('PerSpecies_Abundance_LTER.csv',as.is=T,check.names=F,header=T) #file created in "analyze.LTER.arthropods.R"
# Check ambiguous species codes
output[which(output$Species.code=='' | output$Species.code=='Unknown' | output$Species.code=='unk unk' | output$Species.code=='undet under' | output$Species.code=='na na' | output$Species.code=='na? na?' | output$Species.code=='none none'),]

# Check LTER names
lters = unique(output$LTER)
lters2 = unique(eco.key$LTER)
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


# Arctic
LTER = 'Arctic'
species2 = output$Species[which(output$LTER==LTER)]
u.species = unique(species2)
eco.key2 = eco.key[which(eco.key$LTER==LTER),]
sp.match = match(u.species,eco.key2$Code)
output2 = add.annot(sps=u.species,annot.pos=sp.match,LTER=LTER,ecokey=eco.key2,output2=output2)
unique(output2$Feeding[which(output2$LTER==LTER)])
tmp = output2[which(output2$LTER=='Arctic'),] #Check missing
miss = tmp[which(is.na(tmp$Habitat)),] #PERFECT


# Baltimore
LTER = 'Baltimore'
species2 = output$Species[which(output$LTER==LTER)]
u.species = unique(species2)
eco.key2 = eco.key[which(eco.key$LTER==LTER),]
sp.match = match(u.species,eco.key2$Code)
output2 = add.annot(sps=u.species,annot.pos=sp.match,LTER=LTER,ecokey=eco.key2,output2=output2)
unique(output2$Feeding[which(output2$LTER==LTER)])
tmp = output2[which(output2$LTER=='Baltimore'),] #Check missing
miss = tmp[which(is.na(tmp$Habitat)),] #PERFECT


# Bonanza Creek
LTER = 'BonanzaCreek'
species2 = output$Species[which(output$LTER==LTER)]
u.species = unique(species2)
eco.key2 = eco.key[which(eco.key$LTER==LTER),]
sp.match = match(u.species,eco.key2$Code)
output2 = add.annot(sps=u.species,annot.pos=sp.match,LTER=LTER,ecokey=eco.key2,output2=output2)
unique(output2$Feeding[which(output2$LTER==LTER)])
tmp = output2[which(output2$LTER=='BonanzaCreek'),] #Check missing
miss = tmp[which(is.na(tmp$Habitat)),]
sp.miss = unique(miss$Species.code)
miss.annot = rbind(c('Coleoptera','herbivore','terrestrial','n'),c('Coleoptera','herbivore','terrestrial','n'),c('Coleoptera','herbivore','terrestrial','n'),c('Lepidoptera','herbivore','terrestrial','n'),c('Diptera','detritivore','terrestrial','n'))
for (i in 1:length(sp.miss)){ #Annotate missing
	pos1 = which(output2$LTER=='BonanzaCreek' & output2$Species.code==sp.miss[i])
	output2[pos1,7] = miss.annot[i,1]
	output2[pos1,10] = miss.annot[i,2]
	output2[pos1,11] = miss.annot[i,3]
	output2[pos1,12] = miss.annot[i,4]
}
tmp = output2[which(output2$LTER=='BonanzaCreek'),]
miss = tmp[which(is.na(tmp$Habitat)),] #PERFECT


# Cedar Creek
LTER = 'CedarCreek'
species2 = output$Species[which(output$LTER==LTER)]
u.species = unique(species2)
eco.key2 = eco.key[which(eco.key$LTER==LTER),]
sp.match = match(u.species,eco.key2$Species) #Code field is empty for this LTER. Use Species field instead
output2 = add.annot(sps=u.species,annot.pos=sp.match,LTER=LTER,ecokey=eco.key2,output2=output2)
unique(output2$Feeding[which(output2$LTER==LTER)])
tmp = output2[which(output2$LTER=='CedarCreek'),] #Check missing
miss = tmp[which(is.na(tmp$Habitat)),]
sp.miss = unique(miss$Species.code)
#write.table(sp.miss,'./missing_species/CedarCreek_missing_species.txt',sep='\t',row.names=F,col.names=F,quote=F)
sp.miss1 = read.csv('./missing_species/CedarCreek_missing_species2.csv',as.is=T,check.names=F,header=T)
sp.miss = sp.miss1[,1]
miss.annot = sp.miss1[,2:5]
for (i in 1:length(sp.miss)){ #Annotate missing
	pos1 = which(output2$LTER==LTER & output2$Species.code==sp.miss[i])
	output2[pos1,7] = miss.annot[i,1]
	output2[pos1,10] = miss.annot[i,2]
	output2[pos1,11] = miss.annot[i,3]
	output2[pos1,12] = miss.annot[i,4]
}
tmp = output2[which(output2$LTER==LTER),] #Check missing
miss = tmp[which(is.na(tmp$Habitat)),]
sp.miss = unique(miss$Species.code) #PERFECT (just missing some non-insect entries)


# Central Arizona-Phoenix
LTER = 'CentralArizona-Phoenix'
species2 = output$Species[which(output$LTER==LTER)]
u.species = unique(species2)
write.csv(u.species,'CentralArizona-Phoenix_duplicatenames.csv',quote=F,row.names=F)
eco.key2 = eco.key[which(eco.key$LTER==LTER),]
eco.key2$Species = trimws(eco.key2$Species,which='right')
sp.match = match(u.species,eco.key2$Species)
sp.match1 = match(u.species,eco.key2$Family) 
sp.match[which(is.na(sp.match))] = sp.match1[which(is.na(sp.match))] #fill missing species with family
sp.match2 = match(u.species,eco.key2$Order)
sp.match[which(is.na(sp.match))] = sp.match2[which(is.na(sp.match))] #Fill in missing Species with Family/Order
output2 = add.annot(sps=u.species,annot.pos=sp.match,LTER=LTER,ecokey=eco.key2,output2=output2)
unique(output2$Feeding[which(output2$LTER==LTER)])
tmp = output2[which(output2$LTER==LTER),] #Check missing
miss = tmp[which(is.na(tmp$Habitat)),]
sp.miss = unique(miss$Species.code)
#write.csv(sp.miss,'./missing_species/CentralArizona-Phoenix_missing_species.csv',row.names=F,quote=F)
sp.miss1 = read.csv('./missing_species/CentralArizona-Phoenix_missing_species2.csv',as.is=T,check.names=F,header=F)
sp.miss = sp.miss1[,1]
miss.annot = sp.miss1[,2:5]
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


# Coweeta
LTER = 'Coweeta'
output2$Species.code[which(output2$LTER=='Coweeta')] = gsub('_Abundance','',output2$Species.code[which(output2$LTER=='Coweeta')])
output$Species.code[which(output$LTER=='Coweeta')] = gsub('_Abundance','',output$Species.code[which(output$LTER=='Coweeta')])
species2 = output$Species[which(output$LTER==LTER)]
u.species = unique(species2)
eco.key2 = eco.key[which(eco.key$LTER==LTER),]
sp.match = match(u.species,eco.key2$Code)
sp.match[6:10] = sp.match[1:5]
output2 = add.annot(sps=u.species,annot.pos=sp.match,LTER=LTER,ecokey=eco.key2,output2=output2)
tmp = output2[which(output2$LTER==LTER),] #Check missing
miss = tmp[which(is.na(tmp$Habitat)),]
sp.miss = unique(miss$Species.code) #PERFECT


# Georgia Coastal Ecosystems
LTER = 'GeorgiaCoastal'
species2 = output$Species[which(output$LTER==LTER)]
u.species = unique(species2)
eco.key2 = eco.key[which(eco.key$LTER==LTER),]
sp.match = match(u.species,eco.key2$Code)
output2 = add.annot(sps=u.species,annot.pos=sp.match,LTER=LTER,ecokey=eco.key2,output2=output2)
unique(output2$Feeding[which(output2$LTER==LTER)])
tmp = output2[which(output2$LTER==LTER),] #Check missing
miss = tmp[which(is.na(tmp$Habitat)),]
sp.miss = unique(miss$Species.code)
miss.annot = rbind(c('Decapoda','omnivore','terrestrial','n'),c('Decapoda','omnivore','terrestrial','n'),c('Orthoptera','herbivore','terrestrial','n'),c('Orthoptera','herbivore','terrestrial','n'),c('Hemiptera','herbivore','terrestrial','n'))
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


# Harvard Forest
LTER = 'HarvardForest'
species2 = output$Species[which(output$LTER==LTER)]
u.species = unique(species2)
eco.key2 = eco.key[which(eco.key$LTER==LTER),]
sp.match = match(u.species,eco.key2$Code)
output2 = add.annot(sps=u.species,annot.pos=sp.match,LTER=LTER,ecokey=eco.key2,output2=output2)
unique(output2$Feeding[which(output2$LTER==LTER)])
tmp = output2[which(output2$LTER==LTER),] #Check missing
miss = tmp[which(is.na(tmp$Habitat)),]
sp.miss = unique(miss$Species.code)
#write.csv(sp.miss,'./missing_species/HarvardForest_missing_species.csv',row.names=F,quote=F)
sp.miss1 = read.csv('./missing_species/HarvardForest_missing_species2.csv',as.is=T,check.names=F,header=F)
sp.miss = sp.miss1[,1]
miss.annot = sp.miss1[,2:5]
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


# Hubbard Brook
LTER = 'HubbardBrook'
species2 = output$Species[which(output$LTER==LTER)]
u.species = unique(species2)
eco.key2 = eco.key[which(eco.key$LTER==LTER),]
sp.match = match(u.species,eco.key2$Code)
sp.miss = u.species
miss.annot = cbind(rep('Lepidoptera',5),rep('herbivore',5),rep('terrestrial',5),rep('n',5))
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


# Konza Prairie
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


# North Temperate Lakes
LTER = 'NorthTemperateLakes'
species2 = output$Species[which(output$LTER==LTER)]
u.species = unique(species2)
eco.key2 = eco.key[which(eco.key$LTER==LTER),]
sp.match = match(u.species,eco.key2$Code)
sp.match2 = match(u.species,eco.key2$Species)
sp.match[which(is.na(sp.match))] = sp.match2[which(is.na(sp.match))]
output2 = add.annot(sps=u.species,annot.pos=sp.match,LTER=LTER,ecokey=eco.key2,output2=output2)
unique(output2$Feeding[which(output2$LTER==LTER)])
tmp = output2[which(output2$LTER==LTER),] #Check missing
miss = tmp[which(is.na(tmp$Habitat)),]
sp.miss = unique(miss$Species.code)[1:2] 
miss.annot = rbind(c('Decapoda','omnivore','aquatic','n'),c('Decapoda','omnivore','aquatic','n'))
for (i in 1:length(sp.miss)){ #Annotate missing
	pos1 = which(output2$LTER==LTER & output2$Species.code==sp.miss[i])
	output2[pos1,7] = miss.annot[i,1]
	output2[pos1,10] = miss.annot[i,2]
	output2[pos1,11] = miss.annot[i,3]
	output2[pos1,12] = miss.annot[i,4]
}
tmp = output2[which(output2$LTER==LTER),] #Check missing
miss = tmp[which(is.na(tmp$Habitat)),]
sp.miss = unique(miss$Species.code) #GOOD; just missing zooplankton

# Remove crayfish data from 2018
output2 = output2[-which(output$LTER.site=="NorthTemperateLakes" & output$Year==2018 & (output$Species.code=="PROPINQUUS" | output$Species.code=="RUSTICUS" |  output$Species.code=="VIRILIS")),]


# Sevilleta
LTER = 'Sevilleta'
species2 = output$Species[which(output$LTER==LTER)]
u.species = unique(species2)
eco.key2 = eco.key[which(eco.key$LTER==LTER),]
sp.match = match(u.species,eco.key2$Code)
sp.match2 = match(u.species,gsub(' ','_',eco.key2$Species))
sp.match[which(is.na(sp.match))] = sp.match2[which(is.na(sp.match))]
output2 = add.annot(sps=u.species,annot.pos=sp.match,LTER=LTER,ecokey=eco.key2,output2=output2)
unique(output2$Feeding[which(output2$LTER==LTER)])
tmp = output2[which(output2$LTER==LTER),] #Check missing
miss = tmp[which(is.na(tmp$Habitat)),]
sp.miss = unique(miss$Species.code)
miss.annot = cbind(rep('Orthoptera',length(sp.miss)),rep('herbivore',length(sp.miss)),rep('terrestrial',length(sp.miss)),rep('n',length(sp.miss)))
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


# Midwest Suction Trap Network
LTER = 'MidwestSTN'
species2 = output$Species[which(output$LTER==LTER)]
u.species = unique(species2)
pos1 = which(output2$LTER==LTER)
output2[pos1,7] = 'Hemiptera'
output2[pos1,10] = 'herbivore'
output2[pos1,11] = 'terrestrial'
output2[pos1,12] = 'n'


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


########################################################################Section
# Curate and visualize rank abundance data

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


##### Harvard Forest ants 1
# NOTE: diversity estimates limited to 47 ant species
site.name = 'HarvardForest'
data1 = read.csv("./raw_data/hf118-01-ants.csv",as.is=T,check.names=F,header=T)
data1$Species = data1$code
data1$Year = data1$year
data1$Number = data1$abundance
data1 = data1[which(data1$Year!=2007),]
u.species = unique(data1$Species)
u.years = sort(as.numeric(unique(data1$Year)))
u.years = list(u.years[1:2],u.years[(length(u.years)-1):(length(u.years))])
u.methods = unique(data1$trap.type)
#fix method labels
methods1 = data1$trap.type
methods1[which(methods1=='Pitfall' | methods1=='pitfall' | methods1=='pit')]='pitfall'
methods1[which(methods1=='Litter' | methods1=='LItter' | methods1=='litter')]='litter'
methods1[which(methods1=='Hand sample' | methods1=='hand')]='hand'
methods1[which(methods1=='Bait' | methods1=='Biat' | methods1=='bait')]='bait'
u.methods1 = unique(methods1)
u.locales = u.methods1
data1$Site = methods1
for (l in 1:length(u.locales)){
	l.data = data1[which(data1$Site==u.locales[l]),]
	u.years = sort(as.numeric(unique(l.data$Year)))
	u.years = list(u.years[1:2],u.years[(length(u.years)-1):(length(u.years))])
	for (y in 1:length(u.years)){
		y.data = l.data[which(l.data$Year==u.years[[y]][1] | l.data$Year==u.years[[y]][2]),]
		if (nrow(y.data)==0){
		} else {
			rank.out = calc.ranks(y.data$Species,y.data$Number,rank.list,rx,paste(paste(u.years[[y]][1],u.years[[y]][2],sep='-'),site.name,paste0('ants.',u.locales[l]),sep='_'))
			rank.list = rank.out[[1]]
			rx=rank.out[[2]]
		}
	}
}


##### Harvard Forest ants 3 (Nantucket)
# NOTE: diversity estimates limited to 61 ant species
site.name = 'HarvardForest'
locale.name = 'ants.Nantucket'
data1 = read.csv("./raw_data/hf147-08-nantucket-ants-2004-09.csv",as.is=T,check.names=F,header=T)
data1$Species = data1$code
data1$Site = data1$site
data1$Year = data1$year
data1$Number = data1$qty
data1 = data1[which(!is.na(data1$Species)),]
u.species = unique(data1$Species)
u.years = sort(as.numeric(unique(data1$year)))
u.years = list(u.years[1:2],u.years[(length(u.years)-1):(length(u.years))])
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[[y]][1] | data1$Year==u.years[[y]][2]),]
	rank.out = calc.ranks(y.data$Species,y.data$Number,rank.list,rx,paste(paste(u.years[[y]][1],u.years[[y]][2],sep='-'),site.name,locale.name,sep='_'))
	rank.list = rank.out[[1]]
	rx=rank.out[[2]]
}


##### Midwest Suction Trap Network
site.name = 'MidwestSTN'
data1 = read.table('C:/Users/mcros/Desktop/Postdoc UGA/Aphid_STN/STN_counts_curated_time-consistent.txt',sep='\t',as.is=T,check.names=F,header=T)
data1$Number = apply(data1[,6:7],1,function(x){if(is.na(x[1]) & is.na(x[2])){NA}else{sum(x,na.rm=T)}})
u.years = sort(as.numeric(unique(data1$Year)))
u.species = unique(data1$Species)
data1$Site[grep('Urbana-Champaign',data1$Site)] = 'Urbana-Champaign' #merge trap sites that were moved trivial distances midway through the time series
data1$Site[grep('Saginaw',data1$Site)] = 'Saginaw'
u.locales = sort(unique(data1$Site))
#Merge species counts across sites for each year
data2 = data.frame('Site'=NA,'Year'=NA,'Species'=NA,'Number'=NA)
ix=1
for (i in 1:length(u.years)){
	print(u.years[i])
	y.data = data1[which(data1$Year==u.years[i]),]
	for (j in 1:length(u.species)){
		s.data = y.data[which(y.data$Species==u.species[j]),]
		data2[ix,1] = 'Midwest' #site
		data2[ix,2] = u.years[i] #year
		data2[ix,3] = u.species[j] #species
		data2[ix,4] = sum(s.data$Number,na.rm=T) #number
		ix=ix+1
	}
}
write.table(data2,'MidwestSTN_abundances_merged.txt',sep='\t',quote=F,row.names=F)
u.years = list(u.years[1:2],u.years[(length(u.years)-1):(length(u.years))])
for (y in 1:length(u.years)){
	y.data = data2[which(data2$Year==u.years[[y]][1] | data2$Year==u.years[[y]][2]),]
	rank.out = calc.ranks(y.data$Species,y.data$Number,rank.list,rx,paste(paste(u.years[[y]][1],u.years[[y]][2],sep='-'),site.name,'Midwest',sep='_'))
	rank.list = rank.out[[1]]
	rx=rank.out[[2]]
}


##### Central Arizona-Phoenix - Arthropod sweeps
# NOTE: insects are identified to various levels of taxonomy; estimates of diversity will be downward-biased
site.name = 'CentralArizona-Phoenix'
locale.name = 'sweep'
data1 = read.csv("./raw_data/652_arthropods_e9c22403e0f7243f241ed64862e22e05.csv",as.is=T,check.names=F,header=T)
data1$Year = apply(array(data1$sample_date),1,function(x){strsplit(x,'-')[[1]][1]})
data1$Site = data1$site_code
data1$Species = data1$arthropod_scientific_name
data1$Species = gsub('\\(immature\\)','',data1$Species)
data1$Species = gsub(' \\(immature\\)','',data1$Species)
data1$Species = gsub('Scale Insects \\(Coccoidea\\)','',data1$Species)
data1$Species = trimws(data1$Species,which='right')
data1 = data1[which(data1$Species!='' & data1$Species!='Unknown'),]
data1$Number = data1$number_of_arthropods
u.years = sort(as.numeric(unique(data1$Year)))
u.years = list(u.years[1:2],u.years[(length(u.years)-1):(length(u.years))])
u.species = unique(data1$Species)
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[[y]][1] | data1$Year==u.years[[y]][2]),]
	rank.out = calc.ranks(y.data$Species,y.data$Number,rank.list,rx,paste(paste(u.years[[y]][1],u.years[[y]][2],sep='-'),site.name,locale.name,sep='_'))
	rank.list = rank.out[[1]]
	rx=rank.out[[2]]
}


##### Cedar Creek Ecosystem Arthropods Sweep 1
site.name = 'CedarCreek'
locale.name = 'sweep1'
data1 = read.table("./raw_data/e153_Arthropod sweepnet sampling.txt",sep='\t',as.is=T,check.names=F,header=T)
data1$Year = apply(array(data1$Date),1,function(x){strsplit(x,'/')[[1]][3]})
data1$Species = paste(data1$Genus,data1$Specific.epithet,sep=' ')
data1$Number = data1$Specimens
u.years = sort(as.numeric(unique(data1$Year)))
u.years = list(u.years[1:2],u.years[(length(u.years)-1):(length(u.years))])
u.species = unique(data1$Species)
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[[y]][1] | data1$Year==u.years[[y]][2]),]
	rank.out = calc.ranks(y.data$Species,y.data$Number,rank.list,rx,paste(paste(u.years[[y]][1],u.years[[y]][2],sep='-'),site.name,locale.name,sep='_'))
	rank.list = rank.out[[1]]
	rx=rank.out[[2]]
}


##### Cedar Creek Ecosystem Arthropods Sweep 2
site.name = 'CedarCreek'
locale.name = 'sweep2'
data1 = read.table("./raw_data/e120_Main Plots All Arthropod Insect Sweepnet Sampling 1996-2006.txt",sep='\t',as.is=T,check.names=F,header=T)
data1$Species = paste(data1$Genus,data1$Specific.epithet,sep=' ')
data1 = data1[which(data1$Species!='undet undet' & data1$Species!='unk unk' & data1$Species!='none none' & data1$Species!='na? na?' & data1$Species!='na na'),]
data1$Number = as.numeric(data1$Count)
u.years = sort(as.numeric(unique(data1$Year)))
u.years = list(u.years[1:2],u.years[(length(u.years)-1):(length(u.years))])
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[[y]][1] | data1$Year==u.years[[y]][2]),]
	rank.out = calc.ranks(y.data$Species,y.data$Number,rank.list,rx,paste(paste(u.years[[y]][1],u.years[[y]][2],sep='-'),site.name,locale.name,sep='_'))
	rank.list = rank.out[[1]]
	rx=rank.out[[2]]
}


##### Harvard Forest - Arthropods Carnivorous Plants
# NOTE: used carnivorous plant species as a "site"/"method of observation"
# NOTE: most insects only classified to order. Estimates of diversity will be downward-biased
site.name = 'HarvardForest'
data1 = read.csv("./raw_data/hf111-01-prey.csv",as.is=T,check.names=F,header=T)
data1 = data1[which(data1$units=='number_of_individuals'),c(1:26,29:31)]
data1$Year = apply(array(data1$study),1,function(x){y=strsplit(x,' ')[[1]];return(y[length(y)])})
data1$Site = data1$species
u.locales = unique(data1$Site) #list of carnivorous plant species sampled (use this as a "site")
u.years = sort(as.numeric(unique(data1$Year)))
u.years = list(u.years[1:2],u.years[(length(u.years)-1):(length(u.years))])
for (l in 1:length(u.locales)){
	l.data = data1[which(data1$Site==u.locales[l]),]
	for (y in 1:length(u.years)){
		y.data = l.data[which(l.data$Year==u.years[[y]][1] | l.data$Year==u.years[[y]][2]),]
		n.obs = nrow(y.data)
		if (n.obs==0){
		} else {
			rank.out = calc.ranks(colnames(y.data)[5:29],apply(y.data[,5:29],2,function(x){sum(x,na.rm=T)}),rank.list,rx,paste(paste(u.years[[y]][1],u.years[[y]][2],sep='-'),site.name,u.locales[l],sep='_'))
			rank.list = rank.out[[1]]
			rx=rank.out[[2]]
		}
	}
}


##### Konza Prairie grasshoppers
site.name = 'KonzaPrairie'
locale.name = 'grasshopper'
data1 = read.csv("./raw_data/CGR022.csv",as.is=T,check.names=F,header=T)
data1$Species = data1$SPECIES
data1$Year = data1$RECYEAR
data1$Number = as.numeric(data1$TOTAL)
data1 = data1[which(data1$Species!='unknown' & data1$Species!='unknown ' & data1$Species!=''),]
u.years = sort(as.numeric(unique(data1$Year)))
u.species = unique(data1$Species)
# Fix species names
# Merge duplicates that have different species name spelling
konza.duplicates = read.csv('./taxa_keys/KonzaPrairie_grasshoppers_duplicatenames.csv',as.is=T,check.names=F,header=F)
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


##### Cedar Creek Ecosystem - Old Field Grasshopper Sampling
site.name = 'CedarCreek'
locale.name = 'grasshopper'
data1 = read.table("./raw_data/e014_Core Old Field Grasshopper Sampling.txt",sep='\t',as.is=T,check.names=F,header=T)
data1$Species = paste(data1$Genus,data1$Specific.epithet,sep=' ')
data1$Number = data1$Specimens
u.years = sort(as.numeric(unique(data1$Year)))
u.years = list(u.years[1:2],u.years[(length(u.years)-1):(length(u.years))])
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[[y]][1] | data1$Year==u.years[[y]][2]),]
	n.obs = 1 #standardized sampling protocol
	rank.out = calc.ranks(y.data$Species,y.data$Number,rank.list,rx,paste(paste(u.years[[y]][1],u.years[[y]][2],sep='-'),site.name,locale.name,sep='_'))
	rank.list = rank.out[[1]]
	rx=rank.out[[2]]
}



##### Sevilleta - Grasshoppers
site.name = 'Sevilleta'
locale.name = 'grasshopper'
data1 = read.csv("./raw_data/sev106_hopperdynamics_20150826.txt",as.is=T,check.names=F,header=T)
data1$Year = apply(array(data1$DATE),1,function(x){strsplit(x,'/')[[1]][3]})
data1$Site = data1$SITE
data1$Species = data1$SPECIES
data1$Number = data1$CNT
u.years = sort(as.numeric(unique(data1$Year)))
u.years = list(u.years[1:2],u.years[(length(u.years)-1):(length(u.years))])
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[[y]][1] | data1$Year==u.years[[y]][2]),]
	rank.out = calc.ranks(y.data$Species,y.data$Number,rank.list,rx,paste(paste(u.years[[y]][1],u.years[[y]][2],sep='-'),site.name,locale.name,sep='_'))
	rank.list = rank.out[[1]]
	rx=rank.out[[2]]
}



##### Central Arizona-Phoenix - Arthropod Pitfall
# NOTE: taxon designations down to family level; diversity will be underestimated
site.name = 'CentralArizona-Phoenix'
locale.name = 'pitfall1'
data1 = read.csv("./raw_data/41_core_arthropods_1c843549f049757465771bf4c36e80ab.csv",as.is=T,check.names=F,header=T)
data1$Year = apply(array(data1$sample_date),1,function(x){strsplit(x,'-')[[1]][1]})
data1$Site = data1$site_code
data1$Species = data1$display_name
data1$Species = gsub('\\(adult\\)','',data1$Species)
data1$Species = gsub('\\(winged\\)','',data1$Species)
data1$Species = gsub('\\(larvae\\)','',data1$Species)
data1$Species = gsub('\\(nymph\\)','',data1$Species)
data1$Species = gsub(' \\(immature\\)','',data1$Species)
data1$Species = gsub(' \\(Winged Aphelinidae\\)','',data1$Species)
data1$Species = gsub(' \\(Aphelinidae\\)','',data1$Species)
data1$Species = gsub('_',' ',data1$Species)
data1$Species = gsub(' \\(Centipeds\\)','',data1$Species)
data1$Species = gsub('Scale Insects \\(Coccoidea\\)','',data1$Species)
data1$Species = gsub('Scale Insects \\(Coccoidea - Immature\\)','',data1$Species)
data1$Species = gsub(' \\(imm.\\)','',data1$Species)
data1$Species = gsub(' >10mm','',data1$Species)
data1$Species = gsub(' immature','',data1$Species)
data1$Species = trimws(data1$Species,which='right')
data1 = data1[which(!is.na(data1$Species)),]
data1 = data1[which(data1$Species!='' & data1$Species!='Unknown'),]
data1$Number = apply(data1[,13:17],1,function(x){sum(x,na.rm=T)})
u.years = sort(as.numeric(unique(data1$Year)))
u.years = list(u.years[1:2],u.years[(length(u.years)-1):(length(u.years))])
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[[y]][1] | data1$Year==u.years[[y]][2]),]
	rank.out = calc.ranks(y.data$Species,y.data$Number,rank.list,rx,paste(paste(u.years[[y]][1],u.years[[y]][2],sep='-'),site.name,locale.name,sep='_'))
	rank.list = rank.out[[1]]
	rx=rank.out[[2]]
}


##### Central Arizona-Phoenix - Arthropod Pitfalls 2
site.name = 'CentralArizona-Phoenix'
locale.name = 'pitfall2'
data1 = read.csv("./raw_data/643_mcdowell_pitfall_arthropods_76003ef56b85e55df4e9f2605c3e8492.csv",as.is=T,check.names=F,header=T)
data1$Number = apply(data1[,13:17],1,function(x){sum(x,na.rm=T)})
data1$Year = apply(array(data1$sample_date),1,function(x){strsplit(x,'-')[[1]][1]})
data1$Species = data1$display_name
data1$Species = gsub('\\(adult\\)','',data1$Species)
data1$Species = gsub('\\(winged\\)','',data1$Species)
data1$Species = gsub('\\(larvae\\)','',data1$Species)
data1$Species = gsub('\\(nymph\\)','',data1$Species)
data1$Species = gsub(' \\(immature\\)','',data1$Species)
data1$Species = gsub(' \\(Winged Aphelinidae\\)','',data1$Species)
data1$Species = gsub(' \\(Aphelinidae\\)','',data1$Species)
data1$Species = gsub('_',' ',data1$Species)
data1$Species = gsub(' \\(Centipeds\\)','',data1$Species)
data1$Species = gsub('Scale Insects \\(Coccoidea\\)','',data1$Species)
data1$Species = gsub('Scale Insects \\(Coccoidea - Immature\\)','',data1$Species)
data1$Species = gsub(' \\(imm.\\)','',data1$Species)
data1$Species = gsub(' >10mm','',data1$Species)
data1$Species = gsub(' immature','',data1$Species)
data1$Species = trimws(data1$Species,which='right')
data1 = data1[which(data1$Species!=''),]
data1$Number = apply(data1[,13:17],1,function(x){sum(x,na.rm=T)})
data1 = data1[which(!is.na(data1$Species)),]
data1 = data1[which(data1$Species!='' & data1$Species!='Unknown'),]
u.years = sort(as.numeric(unique(data1$Year)))
u.years = list(u.years[1:2],u.years[(length(u.years)-1):(length(u.years))])
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[[y]][1] | data1$Year==u.years[[y]][2]),]
	rank.out = calc.ranks(y.data$Species,y.data$Number,rank.list,rx,paste(paste(u.years[[y]][1],u.years[[y]][2],sep='-'),site.name,locale.name,sep='_'))
	rank.list = rank.out[[1]]
	rx=rank.out[[2]]
}


##### Sevilleta - Arthropod pitfalls
site.name = 'Sevilleta'
locale.name = 'pitfall'
data1 = read.csv("./raw_data/sev029_arthropop_02162009_0.txt",as.is=T,check.names=F,header=T)
for (i in 1:ncol(data1)){
	data1[which(data1[,i]==-888),i] = NA
}
data1$Species = paste(data1$Genus,data1$Species,sep='_')
data1$Number = data1$Count
data1 = data1[which(data1$Species!='_' & data1$Species!='NA_NA'),]
u.years = sort(as.numeric(unique(data1$Year)))
u.years = list(u.years[1:2],u.years[(length(u.years)-1):(length(u.years))])
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[[y]][1] | data1$Year==u.years[[y]][2]),]
	rank.out = calc.ranks(y.data$Species,y.data$Number,rank.list,rx,paste(paste(u.years[[y]][1],u.years[[y]][2],sep='-'),site.name,locale.name,sep='_'))
	rank.list = rank.out[[1]]
	rx=rank.out[[2]]
}


##### Baltimore
# NOTE: species counts sometimes occur as decimals; I suspect that instars were counted as fractions; I rounded-up to get integer counts
site.name = 'Baltimore'
locale.name = 'mosquito'
data1 = read.csv("./raw_data/Biodiversity_-_Mosquito_-_ovitrap_mosquito_data.csv",as.is=T,check.names=F,header=T)
data1$Year = data1$year
u.years = unique(data1$Year)
u.years = list(u.years[1:2],u.years[(length(u.years)-1):(length(u.years))])
for (i in 1:nrow(data1)){
	for (j in 6:17){
		data1[i,j] = ceiling(data1[i,j])
	}
}
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[[y]][1] | data1$Year==u.years[[y]][2]),]
	rank.out = calc.ranks(colnames(y.data)[6:17],apply(y.data[,6:17],2,function(x){sum(x,na.rm=T)}),rank.list,rx,paste(paste(u.years[[y]][1],u.years[[y]][2],sep='-'),site.name,locale.name,sep='_'))
	rank.list = rank.out[[1]]
	rx=rank.out[[2]]
}


##### Arctic Tundra - Stream invertebrates
site.name = 'Arctic'
locale.name = 'aquatic'
data1 = read.csv("./raw_data/84-98hektot.csv",as.is=T,check.names=F,header=T)
data1$Year = apply(array(data1$Date),1,function(x){paste0(19,strsplit(x,'-')[[1]][3])})
data1 = data1[,which(colnames(data1)!='SNAILS')]
u.years = sort(as.numeric(unique(data1$Year)))
u.years = list(u.years[1:2],u.years[(length(u.years)-1):(length(u.years))])
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[[y]][1] | data1$Year==u.years[[y]][2]),]
	rank.out = calc.ranks(colnames(y.data)[7:20],apply(y.data[,7:20],2,function(x){sum(x,na.rm=T)}),rank.list,rx,paste(paste(u.years[[y]][1],u.years[[y]][2],sep='-'),site.name,locale.name,sep='_'))
	rank.list = rank.out[[1]]
	rx=rank.out[[2]]
}


##### Northern Temperate Lakes - Invertebrates 1
site.name = 'NorthTemperateLakes'
data1 = read.csv("./raw_data/ntl11_1_v8.csv",as.is=T,check.names=F,header=T)
data1$Year = data1$year4
data1$Site = data1$lakeid
data1$Species = data1$taxon_code
data1$Number = data1$number_indiv
u.locales = unique(data1$Site)
u.years = sort(as.numeric(unique(data1$Year)))
u.years = list(u.years[1:2],u.years[(length(u.years)-1):(length(u.years))])
for (l in 1:length(u.locales)){
	l.data = data1[which(data1$Site==u.locales[l]),]
	for (y in 1:length(u.years)){
		y.data = l.data[which(l.data$Year==u.years[[y]][1] | l.data$Year==u.years[[y]][2]),]
		rank.out = calc.ranks(y.data$Species,y.data$Number,rank.list,rx,paste(paste(u.years[[y]][1],u.years[[y]][2],sep='-'),site.name,u.locales[l],sep='_'))
		rank.list = rank.out[[1]]
		rx=rank.out[[2]]
	}
}

# Print curves
lter.sites = apply(array(names(rank.list)),1,function(x){y=strsplit(x,'_')[[1]];paste(y[2],y[3],sep='_')})
uls = unique(lter.sites)
uls.count = apply(array(uls),1,function(x){length(which(lter.sites==x))})
uls2 = uls[which(uls.count>1)]

for (i in 1:length(uls2)){
	spc1 = rank.list[[grep(uls2[i],names(rank.list))[1]]]
	spc2 = rank.list[[grep(uls2[i],names(rank.list))[2]]]
	year1 = strsplit(names(rank.list)[[grep(uls2[i],names(rank.list))[1]]],'_')[[1]][1]
	year2 = strsplit(names(rank.list)[[grep(uls2[i],names(rank.list))[2]]],'_')[[1]][1]
	ymax = max(c(spc1,spc2))
	xmax = max(c(length(spc1),length(spc2)))
	png(paste0('./plots/rank_abundance/',uls2[i],'.png'),width=400,height=400)
	par(oma=c(0,0,0,0),mar=c(5,5,2,0))
	plot(seq(1,xmax,1),c(rep(0,xmax-1),ymax),col='white',xlab='Species rank',ylab='Abundance',main=NULL,cex.axis=1.5,cex.lab=2,bty='n',xaxt='n',yaxt='n')
	axis(1,lwd=2,cex.lab=2,cex.axis=1.5)
	axis(2,lwd=2,cex.lab=2,cex.axis=1.5)
	lines(spc1,lwd=3)
	lines(spc2,lwd=3,col='red',lty=2)
	legend('topright',legend=c(year1,year2),col=c('black','red'),lty=1:2,cex=2,bty='n',lwd=2)
	dev.off()
}


########################################################################Section
# Estimate arthropod taxa abundance time trends

# Define time trend-finding function (created by A. R. Ives)
AR_reml <- function(formula, data = list()) {

	AR_reml_funct <- function(par, x, u) {
		b <- par
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


# Import species abundance data
data1 = read.csv('PerSpecies_Abundance_LTER_annotated.csv',as.is=T,check.names=F,header=T)
data1$LL = paste(data1$LTER.site,data1$Locale,sep='_')
sp.count = apply(array(unique(data1$LL)),1,function(x){length(unique(data1$Species.code[which(data1$LL==x)]))})

###
# relaxed inclusion criteria

LLS = paste(data1$LTER.site,data1$Locale,data1$Species.code,sep='_') #unique LTER by sub-site by species combinations
u.LLS = unique(LLS)
LLS.count = apply(array(u.LLS),1,function(x){d2=data1[which(LLS==x),];length(which(d2$Abundance!=0))})
u.LLS2 = u.LLS[which(LLS.count>0)]
time.trends.Z = data.frame('LTER'=NA,'Site'=NA,'Species'=NA,'MSE'=-999,'b'=-999,'coef_int'=-999,'coef_slope'=-999,'Pr_int'=-999,'Pr_slope'=-999,'logLik'=-999,'length'=-999,'Y1'=NA,'Y2'=NA,'Feeding'=NA,'Habitat'=NA,'Pollinator'=NA,'Order'=NA)
tx = 1
for (l in 1:length(u.LLS2)){
	print(l)
	lls.data = data1[which(LLS==u.LLS2[l]),]
	LL = paste0(strsplit(u.LLS2[l],'_')[[1]][1:2],collapse='_') #trim "no sample" years
	if (LL=='HarvardForest_ants.hand' | LL=='HarvardForest_ants.litter' | LL=='HarvardForest_ants.bait'){ #deal with HarvardForest cases where time series differed in length
		lls.data = lls.data[-which(lls.data$Year>2008),] #three ant datasets have no data past 2008
	}
	### Quality threshold
	if ( (length(which(!is.na(lls.data$Abundance))) > 3) & (length(which(lls.data$Abundance[!is.na(lls.data$Abundance)] > 0)) > 0) ){
	### 
		cmin = min(lls.data$Abundance[which(lls.data$Abundance>0)],na.rm=T)
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
			X[which(!is.na(X))] = lls.data$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
			Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
			t.scale = 1:length(lls.years.stretch)
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
		png(paste0('./plots/time_trends/abundance/lineplots/',gsub('\\?','',u.LLS[l]),'.png'))
		plot(t.scale,Z,type='l',main=u.LLS2[l],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2)
		abline(a=add.coef.Z[3],b=add.coef.Z[4],lty=2,col='red',lwd=1.5)
		legend('topright',legend=c(paste0('slope = ',round(add.coef.Z[4],2)),paste0('length = ',length(which(!is.na(Z)))),paste0('autocor = ',round(add.coef.Z[2],2)),paste0('logLik = ',round(add.coef.Z[7],2))),ncol=1)
		dev.off()
	} else {
		#ignore low quality time series
	}
}
write.csv(time.trends.Z,'time_trends_arthropods_relaxed.csv',quote=F,row.names=F)

###
# moderately strict inclusion criteria

LLS = paste(data1$LTER.site,data1$Locale,data1$Species.code,sep='_') #unique LTER by sub-site by species combinations
u.LLS = unique(LLS)
LLS.count = apply(array(u.LLS),1,function(x){d2=data1[which(LLS==x),];length(which(d2$Abundance!=0))})
u.LLS2 = u.LLS[which(LLS.count>0)]
time.trends.Z = data.frame('LTER'=NA,'Site'=NA,'Species'=NA,'MSE'=-999,'b'=-999,'coef_int'=-999,'coef_slope'=-999,'Pr_int'=-999,'Pr_slope'=-999,'logLik'=-999,'length'=-999,'Y1'=NA,'Y2'=NA,'Feeding'=NA,'Habitat'=NA,'Pollinator'=NA)
tx = 1
tx = 1
for (l in 1:length(u.LLS2)){
	print(l)
	lls.data = data1[which(LLS==u.LLS2[l]),]
	LL = paste0(strsplit(u.LLS2[l],'_')[[1]][1:2],collapse='_') #trim "no sample" years
	if (LL=='HarvardForest_ants.hand' | LL=='HarvardForest_ants.litter' | LL=='HarvardForest_ants.bait'){ #deal with HarvardForest cases where time series differed in length
		lls.data = lls.data[-which(lls.data$Year>2008),] #three ant datasets have no data past 2008
	}
	### Quality threshold
	if ( (length(which(!is.na(lls.data$Abundance))) > 7) & (length(which(lls.data$Abundance[!is.na(lls.data$Abundance)] > 0)) > 3) ){
	###
		cmin = min(lls.data$Abundance[which(lls.data$Abundance>0)],na.rm=T)
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
			X[which(!is.na(X))] = lls.data$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
			Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
			t.scale = 1:length(lls.years.stretch)
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
		tx = tx + 1
		# Plot high quality time series
		png(paste0('./plots/time_trends/abundance/lineplots/',gsub('\\?','',u.LLS[l]),'.png'))
		plot(t.scale,Z,type='l',main=u.LLS2[l],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2)
		abline(a=add.coef.Z[3],b=add.coef.Z[4],lty=2,col='red',lwd=1.5)
		legend('topright',legend=c(paste0('slope = ',round(add.coef.Z[4],2)),paste0('length = ',length(which(!is.na(Z)))),paste0('autocor = ',round(add.coef.Z[2],2)),paste0('logLik = ',round(add.coef.Z[7],2))),ncol=1)
		dev.off()
	} else {
		#ignore low quality time series
	}
}
write.csv(time.trends.Z,'time_trends_arthropods_moderate.csv',quote=F,row.names=F)



###
# strict inclusion criteria

LLS = paste(data1$LTER.site,data1$Locale,data1$Species.code,sep='_') #unique LTER by sub-site by species combinations
u.LLS = unique(LLS)
LLS.count = apply(array(u.LLS),1,function(x){d2=data1[which(LLS==x),];length(which(d2$Abundance!=0))})
u.LLS2 = u.LLS[which(LLS.count>0)]
time.trends.Z = data.frame('LTER'=NA,'Site'=NA,'Species'=NA,'MSE'=-999,'b'=-999,'coef_int'=-999,'coef_slope'=-999,'Pr_int'=-999,'Pr_slope'=-999,'logLik'=-999,'length'=-999,'Y1'=NA,'Y2'=NA,'Feeding'=NA,'Habitat'=NA,'Pollinator'=NA)
tx = 1
for (l in 1:length(u.LLS2)){
	print(l)
	lls.data = data1[which(LLS==u.LLS2[l]),]
	LL = paste0(strsplit(u.LLS2[l],'_')[[1]][1:2],collapse='_') #trim "no sample" years
	if (LL=='HarvardForest_ants.hand' | LL=='HarvardForest_ants.litter' | LL=='HarvardForest_ants.bait'){ #deal with HarvardForest cases where time series differed in length
		lls.data = lls.data[-which(lls.data$Year>2008),] #three ant datasets have no data past 2008
	}
	### Quality threshold
	if ( (length(which(!is.na(lls.data$Abundance))) > 14) & (length(which(lls.data$Abundance[!is.na(lls.data$Abundance)] > 0)) > 9) ){
	###
		cmin = min(lls.data$Abundance[which(lls.data$Abundance>0)],na.rm=T)
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
			X[which(!is.na(X))] = lls.data$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
			Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
			t.scale = 1:length(lls.years.stretch)
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
		tx = tx + 1
		# Plot high quality time series
		png(paste0('./plots/time_trends/abundance/lineplots/',gsub('\\?','',u.LLS[l]),'.png'))
		plot(t.scale,Z,type='l',main=u.LLS2[l],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2)
		abline(a=add.coef.Z[3],b=add.coef.Z[4],lty=2,col='red',lwd=1.5)
		legend('topright',legend=c(paste0('slope = ',round(add.coef.Z[4],2)),paste0('length = ',length(which(!is.na(Z)))),paste0('autocor = ',round(add.coef.Z[2],2)),paste0('logLik = ',round(add.coef.Z[7],2))),ncol=1)
		dev.off()
	} else {
		#ignore low quality time series
	}
}
time.trends.Z2 = time.trends.Z[-which(time.trends.Z$b>0.99),]
write.csv(time.trends.Z2,'time_trends_arthropods_strict.csv',quote=F,row.names=F)


###############################################################Section
# Visualize overall patterns in abundance trends

library('scales')

### Boxplots
# Relaxed criteria

trends.r = read.csv('/Users/caitlinmiller/Downloads/External_Database_S2_time_trends_arthropods_relaxed.csv',as.is=T,check.names=F,header=T)

# Group studies by site/taxa/method
trends.r$Site[which(trends.r$LTER=='MidwestSTN')] = 'aphids' #merge aphid records into one boxplot
trends.r$Site[which(trends.r$LTER=='HubbardBrook')] = 'butterflies/moths' #merge lepidoptera records into one boxplot
trends.r$Site[which(trends.r$LTER=='HarvardForest' & trends.r$Site!='ants' & trends.r$Site!='ants.pitfall' & trends.r$Site!='ants.bait' & trends.r$Site!='ants.hand' & trends.r$Site!='ants.litter')] = 'ticks' #merge tick records into one boxplot
trends.r$Site[grep('ants',trends.r$Site)] = 'ants' #merge tick records into one boxplot
trends.r$Site[which(trends.r$LTER=='NorthTemperateLakes' & trends.r$Site!='Crayfish')] = 'lake insects' #merge insect records into one boxplot
trends.r$Site[which(trends.r$LTER=='NorthTemperateLakes' & trends.r$Site=='Crayfish')] = 'crayfish' 
trends.r$Site[which(trends.r$LTER=='KonzaPrairie' & trends.r$Site=='gall')] = 'gall insects' 
trends.r$Site[which(trends.r$LTER=='GeorgiaCoastal' & (trends.r$Site=='Burrowing Crab' | trends.r$Site=='Fiddler Crab'))] = 'crabs'  
trends.r$Site[which(trends.r$LTER=='GeorgiaCoastal' & trends.r$Site=='Salt Marsh' & trends.r$Species=='Acrididae')] = 'grasshoppers'  
trends.r$Site[which(trends.r$LTER=='GeorgiaCoastal' & trends.r$Site=='Salt Marsh' & trends.r$Species=='Prokelisia marginata')] = 'planthoppers'  
trends.r$Site[which(trends.r$LTER=='CentralArizona-Phoenix' & (trends.r$Site=='pitfall1' | trends.r$Site=='pitfall2'))] = 'pitfalls' #merge pitfall records
trends.r$Site[which(trends.r$Site=='grasshopper')] = 'grasshoppers'
trends.r$Site[which(trends.r$Site=='aquatic')] = 'stream insects'
trends.r$Site[which(trends.r$Site=='leafminer')] = 'leaf miners'
trends.r$Site[which(trends.r$Site=='bark beetle')] = 'bark beetles'
trends.r$Site[which(trends.r$Site=='mosquito')] = 'mosquitoes'
trends.r$Site[grep('sweep',trends.r$Site)] = 'sweep net'
trends.r$Site[which(trends.r$LTER=='CedarCreek' & (trends.r$Site=='sweep1' | trends.r$Site=='sweep2'))] = 'sweep nets' 

trends.r$LL = as.factor(paste(trends.r$LTER,trends.r$Site,sep='_'))
lters = unique(apply(array(levels(trends.r$LL)),1,function(x){strsplit(x,'_')[[1]][1]}))
cols = c('lightblue',rep('white',7),'lightblue',rep('white',9),'lightblue','lightblue','white','white')
cbind(levels(trends.r$LL),cols)
# All taxa
png('ALL_time_trends_boxplot_relaxed_noannotations.png',width=500*2,height=500)
par(mfrow=c(1,1),cex.lab=2,cex.axis=1.5,lwd=2,oma=c(0,0,0,0),mar=c(15,8,1,1))
plt = plot(trends.r$LL,trends.r$coef_slope,col=cols,ylim=c(-8,8),xaxt='n',yaxt='n',ylab='Change in abundance',axes=F)
axis(2,at=seq(-8,8,2),labels=seq(-8,8,2),lwd=2)
abline(a=0,b=0,lty=1,col='black',lwd=2)
text(x=c(1,2,3.5,5.5,7.5,9,11,13.5,15,16,17.5,19.5,21.5),y=par("usr")[3]-0.2,labels=lters,srt=45,pos=2,xpd=T,cex=2)
tcols = c('blue',rep('black',7),'blue',rep('black',9),'blue','blue','black','black')
text(x=seq(1,length(levels(trends.r$LL)),1)-0.35,y=par("usr")[3]+0.5,col=tcols,labels=apply(array(levels(trends.r$LL)),1,function(x){strsplit(x,'_')[[1]][2]}),srt=90,adj=0,xpd=T,cex=1.2)
abline(a=0,b=0,lty=1,col='black',lwd=2)
abline(v=1.5,col='gray50');abline(v=2.5,col='gray50');abline(v=4.5,col='gray50');abline(v=6.5,col='gray50')
abline(v=8.5,col='gray50');abline(v=9.5,col='gray50');#abline(v=12.5,col='gray50');#abline(v=14.5,col='gray50')
abline(v=15.5,col='gray50');abline(v=16.5,col='gray50');abline(v=18.5,col='gray50');abline(v=20.5,col='gray50')
dev.off()

# Standard deviations of medians from zero for ms
sort(plt$stats[3,])

# Carnivores
png('./plots/time_trends/abundance/CARNIVORE_time_trends_boxplot_relaxed.png',width=500*2,height=500)
par(mfrow=c(1,1),cex.lab=2,cex.axis=1.5,lwd=2,oma=c(0,0,0,0),mar=c(15,8,1,1))
trends2 = trends.r[which(trends.r$Feeding=='carnivore'),]
plot(trends2$LL,trends2$coef_slope,col=cols,ylim=c(-8,8),xaxt='n',yaxt='n',ylab='Change in abundance',axes=F)
axis(2,at=seq(-8,8,2),labels=seq(-8,8,2),lwd=2)
abline(a=0,b=0,lty=1,col='black',lwd=2)
dev.off()

# Herbivores
png('./plots/time_trends/abundance/HERBIVORE_time_trends_boxplot_relaxed.png',width=500*2,height=500)
par(mfrow=c(1,1),cex.lab=2,cex.axis=1.5,lwd=2,oma=c(0,0,0,0),mar=c(15,8,1,1))
trends2 = trends.r[which(trends.r$Feeding=='herbivore'),]
plot(trends2$LL,trends2$coef_slope,col=cols,ylim=c(-8,8),xaxt='n',yaxt='n',ylab='Change in abundance',axes=F)
axis(2,at=seq(-8,8,2),labels=seq(-8,8,2),lwd=2)
abline(a=0,b=0,lty=1,col='black',lwd=2)
dev.off()

# Omnivores
png('./plots/time_trends/abundance/OMNIVORE_time_trends_boxplot_relaxed.png',width=500*2,height=500)
par(mfrow=c(1,1),cex.lab=2,cex.axis=1.5,lwd=2,oma=c(0,0,0,0),mar=c(15,8,1,1))
trends2 = trends.r[which(trends.r$Feeding=='omnivore'),]
plot(trends2$LL,trends2$coef_slope,col=cols,ylim=c(-8,8),xaxt='n',yaxt='n',ylab='Change in abundance',axes=F)
axis(2,at=seq(-8,8,2),labels=seq(-8,8,2),lwd=2)
abline(a=0,b=0,lty=1,col='black',lwd=2)
dev.off()

# Detritivores
png('./plots/time_trends/abundance/DETRITIVORE_time_trends_boxplot_relaxed.png',width=500*2,height=500)
par(mfrow=c(1,1),cex.lab=2,cex.axis=1.5,lwd=2,oma=c(0,0,0,0),mar=c(15,8,1,1))
trends2 = trends.r[which(trends.r$Feeding=='detritivore'),]
plot(trends2$LL,trends2$coef_slope,col=cols,ylim=c(-8,8),xaxt='n',yaxt='n',ylab='Change in abundance',axes=F)
axis(2,at=seq(-8,8,2),labels=seq(-8,8,2),lwd=2)
abline(a=0,b=0,lty=1,col='black',lwd=2)
dev.off()

# Parasite
png('./plots/time_trends/abundance/PARASITE_time_trends_boxplot_relaxed.png',width=500*2,height=500)
par(mfrow=c(1,1),cex.lab=2,cex.axis=1.5,lwd=2,oma=c(0,0,0,0),mar=c(15,8,1,1))
trends2 = trends.r[which(trends.r$Feeding=='parasite'),]
plot(trends2$LL,trends2$coef_slope,col=cols,ylim=c(-8,8),xaxt='n',yaxt='n',ylab='Change in abundance',axes=F)
axis(2,at=seq(-8,8,2),labels=seq(-8,8,2),lwd=2)
abline(a=0,b=0,lty=1,col='black',lwd=2)
dev.off()

# Aquatic
png('./plots/time_trends/abundance/AQUATIC_time_trends_boxplot_relaxed.png',width=500*2,height=500)
par(mfrow=c(1,1),cex.lab=2,cex.axis=1.5,lwd=2,oma=c(0,0,0,0),mar=c(15,8,1,1))
trends2 = trends.r[which(trends.r$Habitat=='aquatic'),]
plot(trends2$LL,trends2$coef_slope,col=cols,ylim=c(-8,8),xaxt='n',yaxt='n',ylab='Change in abundance',axes=F)
axis(2,at=seq(-8,8,2),labels=seq(-8,8,2),lwd=2)
abline(a=0,b=0,lty=1,col='black',lwd=2)
dev.off()

# Terrestrial
png('./plots/time_trends/abundance/TERRESTRIAL_time_trends_boxplot_relaxed.png',width=500*2,height=500)
par(mfrow=c(1,1),cex.lab=2,cex.axis=1.5,lwd=2,oma=c(0,0,0,0),mar=c(15,8,1,1))
trends2 = trends.r[which(trends.r$Habitat=='terrestrial'),]
plot(trends2$LL,trends2$coef_slope,col=cols,ylim=c(-8,8),xaxt='n',yaxt='n',ylab='Change in abundance',axes=F)
axis(2,at=seq(-8,8,2),labels=seq(-8,8,2),lwd=2)
abline(a=0,b=0,lty=1,col='black',lwd=2)
dev.off()

# Pollinator
png('./plots/time_trends/abundance/POLLINATOR_time_trends_boxplot_relaxed.png',width=500*2,height=500)
par(mfrow=c(1,1),cex.lab=2,cex.axis=1.5,lwd=2,oma=c(0,0,0,0),mar=c(15,8,1,1))
trends2 = trends.r[which(trends.r$Pollinator=='y'),]
plot(trends2$LL,trends2$coef_slope,col=cols,ylim=c(-8,8),xaxt='n',yaxt='n',ylab='Change in abundance',axes=F)
axis(2,at=seq(-8,8,2),labels=seq(-8,8,2),lwd=2)
abline(a=0,b=0,lty=1,col='black',lwd=2)
dev.off()

# Parasitoids
png('./plots/time_trends/abundance/PARASITOID_time_trends_boxplot_relaxed.png',width=500*2,height=500)
par(mfrow=c(1,1),cex.lab=2,cex.axis=1.5,lwd=2,oma=c(0,0,0,0),mar=c(15,8,1,1))
trends2 = trends.r[which(trends.r$Feeding=='parasitoid'),]
plot(trends2$LL,trends2$coef_slope,col=cols,ylim=c(-8,8),xaxt='n',yaxt='n',ylab='Change in abundance',axes=F)
axis(2,at=seq(-8,8,2),labels=seq(-8,8,2),lwd=2)
abline(a=0,b=0,lty=1,col='black',lwd=2)
dev.off()


##################################################################Section
# Estimate diversity time trends

# Alpha diversity and rank abundance curves

files = list.files('./summary_tables/arthropods/Site-level',full.names=T)
files2 = files[grep('TotalChange',files)]

variables = c('N.species','N.species.rarefied','Species.evenness','Species.decay.rate','Fishers.alpha','Dominance','Species.evenness.rarefied') #diversity metrics
thresh.na = 3
thresh.0 = 3
lost.sites = c()
all.diversity.metrics = c()

for (f in 1:length(files2)){
	print(f)
	dat1 = read.table(files2[f],sep='\t',as.is=T,check.names=F,header=T)
	all.diversity.metrics = data.frame(rbind(all.diversity.metrics,dat1))
	uL = unique(dat1$Locale)
	for (l in 1:length(uL)){
		if (files2[f]=="./summary_tables/arthropods/Site-level/HarvardForest_ants1_TotalChange.txt"){ #deal with HarvardForest cases where time series differed in length
			if (length(grep('pitfall',uL[l]))==0){ #trim artificially extended time series
				dat1 = dat1[-which(dat1$Locale==uL[l] & dat1$Year>2008),] #three ant datasets have no data past 2008
			} else {
			}
		}
		for (v in 1:length(variables)){
			vdat = dat1[which(dat1$Locale==uL[l]),which(colnames(dat1)==variables[v])]
			# Get time trends for site/method l and metric v
			trends = data.frame('LTER'=NA,'Site'=NA,'MSE'=-999,'b'=-999,'coef_int'=-999,'coef_slope'=-999,'Pr_int'=-999,'Pr_slope'=-999,'logLik'=-999,'length'=-999,'Y1'=NA,'Y2'=NA)
			tx = 1
			### quality filtering step ###
			if ( (length(which(!is.na(vdat))) > thresh.na) & (length(which(vdat[!is.na(vdat)] > 0)) > thresh.0) ){
			### quality filtering step ###
				years.stretch = seq(dat1$Year[1],dat1$Year[nrow(dat1)],1) #expand time series to include missing years
				y.match = match(years.stretch,dat1$Year)
				X = y.match
				X[which(!is.na(X))] = vdat[X[which(!is.na(X))]] #fill-in values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(years.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale) #original transform: scale between 0 and 1
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				add.coef.Z = c(arr.Z[[1]],arr.Z[[2]],arr.Z[[3]][1,1],arr.Z[[3]][2,1],arr.Z[[5]][1],arr.Z[[5]][2],arr.Z[[6]])
				trends[tx,1] = unique(dat1$LTER.site) #LTER
				trends[tx,2] = uL[l] #Site
				trends[tx,3] = add.coef.Z[1] #MSE
				trends[tx,4] = add.coef.Z[2] #b
				trends[tx,5] = add.coef.Z[3] #intercept
				trends[tx,6] = add.coef.Z[4] #slope
				trends[tx,7] = add.coef.Z[5] #pval intercept
				trends[tx,8] = add.coef.Z[6] #pval slope
				trends[tx,9] = add.coef.Z[7] #logLik
				trends[tx,10] = length(which(!is.na(Z))) #length
				trends[tx,11] = years.stretch[1] #first year
				trends[tx,12] = years.stretch[length(years.stretch)] #last year
				tx = tx + 1				
				# Plot Z transform on t.scale fit
				png(paste0('./plots/time_trends/diversity/fitted/',variables[v],'_',unique(dat1$LTER.site),'_',uL[l],'.png'))
				plot(Z ~ t.scale, main=paste0(variables[v], ": b = ",round(arr.Z$b, digits=3)))
				curve(arr.Z$coef[1] + arr.Z$coef[2] * x, from=0, to=max(t.scale), add=T)
				dev.off()
				#Plot "high quality" time series
				png(paste0('./plots/time_trends/diversity/lineplots/',variables[v],'_',unique(dat1$LTER.site),'_',uL[l],'.png'))
				plot(t.scale,Z,type='l',main=paste0(unique(dat1$LTER.site),'_',uL[l]),xlab='Scaled time',ylab=variables[v],lwd=2)
				abline(a=add.coef.Z[3],b=add.coef.Z[4],lty=2,col='red',lwd=1.5)
				legend('topright',legend=c(paste0('slope = ',round(add.coef.Z[4],2)),paste0('length = ',length(which(!is.na(Z)))),paste0('autocor = ',round(add.coef.Z[2],2)),paste0('logLik = ',round(add.coef.Z[7],2))),ncol=1)
				dev.off()
				#Write output file
				write.csv(trends,paste0('./diversity_trends/',variables[v],'_',unique(dat1$LTER.site),'_',uL[l],'.csv'),quote=F,row.names=F)
			} else {
				lost.sites = c(lost.sites,files2[f])
			}
		}
	}
}
str(all.diversity.metrics) #dataframe with all diversity trends (contrast with "trends" below which has only 'high quality' trends)
write.csv(all.diversity.metrics,'arthropods_ALPHAdiversity.csv',quote=F,row.names=F)

# Compile trends for each metric
files = list.files('./diversity_trends',full.names=T)
files = files[-grep('Beta',files)] #necessary if beta diversity metrics have already been calculated
files.split = lapply(files,function(x){strsplit(gsub('.csv','',strsplit(x,'/')[[1]][3]),'_')[[1]]})
LL = unlist(lapply(files.split,function(x){paste(x[2],x[3],sep='_')}))
uLL = unique(LL)
variables = c('N.species','N.species.rarefied','Species.evenness','Species.decay.rate','Fishers.alpha','Dominance','Species.evenness.rarefied') #diversity metrics
trends = data.frame('LTER.Locale'=uLL,'N.species'=rep(NA,length(uLL)),'N.species.rarefied'=rep(NA,length(uLL)),'Species.evenness'=rep(NA,length(uLL)),'Species.decay.rate'=rep(NA,length(uLL)),'Fishers.alpha'=rep(NA,length(uLL)),'Dominance'=rep(NA,length(uLL)),'Species.evenness.rarefied'=rep(NA,length(uLL)),stringsAsFactors=F)
for (f in 1:length(files)){
	split1 = strsplit(gsub('.csv','',strsplit(files[f],'/')[[1]][3]),'_')[[1]]
	LL1 = paste(split1[2],split1[3],sep='_')
	dat1 = read.csv(files[f],as.is=T,check.names=F,header=T)
	trends[which(uLL==LL1),which(variables==split1[1])+1] = dat1$coef_slope
}
write.csv(trends,'time_trends_arthropods_ALPHAdiversity.csv',quote=F,row.names=F)

trends = read.csv('time_trends_arthropods_ALPHAdiversity.csv',as.is=T,check.names=F,header=T)

# Compare relationships among diversity measures
trends$L = as.factor(apply(array(trends$LTER.Locale),1,function(x){strsplit(x,'_')[[1]][1]}))
library('RColorBrewer')
cols = brewer.pal(9,'Paired')
png('./plots/alphadiversity_metrics_change_comparison.png',width=480*3)
par(mfrow=c(1,3),cex=2,cex.lab=1.5,cex.axis=1.2,lwd=2,oma=c(0,0,0,0),mar=c(5,5,1,0))
plot(trends$N.species.rarefied,trends$Species.evenness,col=cols[trends$L],pch=16,cex=1.5,xlab='Richness (rarefied)',ylab='Evenness')
lm1 = lm(Species.evenness~N.species.rarefied,trends)
abline(lm1,lty=2)
summary(lm1)
plot(trends$Dominance,trends$Species.evenness,col=cols[trends$L],pch=16,cex=1.5,xlab='Dominance',ylab='Evenness')
lm2 = lm(Species.evenness~Dominance,trends)
abline(lm2,lty=2)
summary(lm2)
plot(trends$Dominance,trends$Species.evenness,type='n',ann=F,axes=F)
legend('topleft',legend=levels(trends$L),col=cols,pch=16,cex=1,bty='n')
dev.off()

# Group studies by site/taxa/method
trends$LTER = apply(array(trends$LTER.Locale),1,function(x){strsplit(x,'_')[[1]][1]})
lters = sort(unique(trends$LTER))
trends$Site = apply(array(trends$LTER.Locale),1,function(x){strsplit(x,'_')[[1]][2]})
trends$Site[which(trends$LTER=='MidwestSTN')] = 'aphids' #merge aphid records into one boxplot
trends$Site[which(trends$LTER=='HarvardForest' & (trends$Site=='ants' | trends$Site=='ants.pitfall' | trends$Site=='ants.hand' | trends$Site=='ants.Nantucket' | trends$Site=='ants.bait' | trends$Site=='ants.litter'))] = 'ants' #merge ant records into one boxplot
trends$Site[which(trends$LTER=='NorthTemperateLakes' & (trends$Site=='CR' | trends$Site=='SP' | trends$Site=='TR'))] = 'lake insects' #merge insect records into one boxplot
trends$Site[which(trends$LTER=='KonzaPrairie' & trends$Site=='gall')] = 'gall insects' 
trends$Site[which(trends$LTER=='NorthTemperateLakes' & trends$Site=='Crayfish')] = 'crayfish' 
trends$Site[which(trends$LTER=='CentralArizona-Phoenix' & (trends$Site=='pitfall1' | trends$Site=='pitfall2'))] = 'pitfall' #merge pitfall records
trends$Site[which(trends$Site=='grasshopper')] = 'grasshoppers'
trends$Site[which(trends$Site=='aquatic')] = 'stream insects'
trends$Site[grep('sweep',trends$Site)] = 'sweep net'
trends$Site[which(trends$Site=='mosquito')] = 'mosquitoes'
trends$Site[which(trends$LTER=='CedarCreek' & (trends$Site=='sweep1' | trends$Site=='sweep2'))] = 'sweep nets' 

# Boxplots
trends$LL = as.factor(paste(trends$LTER,trends$Site,sep='_'))
levels(trends$LL)
variables = c('N.species','N.species.rarefied','Species.evenness','Species.decay.rate','Fishers.alpha','Dominance','Species.evenness.rarefied')
var2 = c('richness (raw)','richness (rarefied)','evenness','rank abundance decay',"Fisher's alpha",'dominance','evenness (rarefied)')
for (v in 1:length(variables)){
	png(paste0('./plots/time_trends/diversity/',variables[v],'_trends_boxplot.png'),width=600,height=500)
	par(mfrow=c(1,1),cex.lab=2,cex.axis=1.5,lwd=1.5,oma=c(0,0,0,0),mar=c(15,8,1,1),lwd=2)
	plot(trends$LL,trends[,which(colnames(trends)==variables[v])],ylim=c(-10,10),yaxt='n',xaxt='n',axes=F,ylab=paste0('Change in ',var2[v]))
	axis(2,at=seq(-10,10,2),labels=seq(-10,10,2),lwd=2)
	text(x=c(1,2,3.5,5.5,7,8,9,10,11.5),y=par("usr")[3]-0.2,labels=lters,srt=45,pos=2,xpd=T,cex=2)
	text(x=seq(1,length(levels(trends$LL)),1)-0.35,y=par("usr")[3]+0.5,labels=apply(array(levels(trends$LL)),1,function(x){strsplit(x,'_')[[1]][2]}),srt=90,adj=0,xpd=T,cex=1.2)
	abline(a=0,b=0,lty=1,col='black',lwd=2)
	dev.off()	
}

# One Sample T-test
uL = unique(trends$LTER)
means.rS = apply(array(uL),1,function(x){mean(trends$N.species.rarefied[which(trends$LTER==x)],na.rm=T)})
means.PIE = apply(array(uL),1,function(x){mean(trends$Species.evenness[which(trends$LTER==x)],na.rm=T)})
t.rS = t.test(as.numeric(means.rS))
t.PIE = t.test(as.numeric(means.PIE))

# Barplots
# Richness (rarefied)
mean1 = t.rS[['estimate']]
uppers = t.rS[['conf.int']][2]
lowers = t.rS[['conf.int']][1]
png('./plots/t-test_richness-rarefied_barplot.png',width=200,height=500)
par(cex.lab=2,cex.axis=1.5,oma=c(0,0,0,0),mar=c(15,8,2,2),lwd=2)
bars = barplot(mean1,ylim=c(-2,2),names.arg='',ylab='Average change in richness (rarefied)',main=NULL,lwd=2)
abline(h=0)
axis(2,lwd=2)
text(x=bars,y=par("usr")[3]-0.05,srt=90,adj=1,labels='LTER',xpd=T,cex=2)
segments(bars,lowers,bars,uppers,lwd=2)
arrows(bars,lowers,bars,uppers,lwd=2,angle=90,code=3,length=(par("usr")[2]-par("usr")[1])/4)
dev.off()

# Evenness
mean1 = t.PIE[['estimate']]
uppers = t.PIE[['conf.int']][2]
lowers = t.PIE[['conf.int']][1]
png('./plots/t-test_evenness_barplot.png',width=200,height=500)
par(cex.lab=2,cex.axis=1.5,oma=c(0,0,0,0),mar=c(15,8,2,2),lwd=2)
bars = barplot(mean1,ylim=c(-2,2),names.arg='',ylab='Average change in evenness',main=NULL,lwd=2)
abline(h=0)
axis(2,lwd=2)
text(x=bars,y=par("usr")[3]-0.05,srt=90,adj=1,labels='LTER',xpd=T,cex=2)
segments(bars,lowers,bars,uppers,lwd=2)
arrows(bars,lowers,bars,uppers,lwd=2,angle=90,code=3,length=(par("usr")[2]-par("usr")[1])/4)
dev.off()

###
# BETA DIVERSITY

files = list.files('./summary_tables/arthropods/Site-level',full.names=T)
files2 = files[grep('Beta',files)]

variables = c('Beta.2','Beta.j','Beta.bray')
thresh.na = 2
thresh.0 = 0

all.metrics = c()
lost.sites=c()

for (f in 1:length(files2)){
	print(f)
	dat1 = read.table(files2[f],sep='\t',as.is=T,check.names=F,header=T)
	all.metrics = data.frame(rbind(all.metrics,dat1))
	uL = unique(dat1$Locale)
	for (l in 1:length(uL)){
		if (files2[f]=="./summary_tables/arthropods/Site-level/HarvardForest_ants1_TotalChange.txt"){ #deal with HarvardForest cases where time series differed in length
			if (length(grep('pitfall',uL[l]))==0){ #trim artificially extended time series
				dat1 = dat1[-which(dat1$Locale==uL[l] & dat1$Year>2008),] #three ant datasets have no data past 2008
			} else {
			}
		}
		for (v in 1:length(variables)){
			vdat = dat1[which(dat1$Locale==uL[l]),which(colnames(dat1)==variables[v])]
			# Get time trends for site/method l and metric v
			trends = data.frame('LTER'=NA,'Site'=NA,'MSE'=-999,'b'=-999,'coef_int'=-999,'coef_slope'=-999,'Pr_int'=-999,'Pr_slope'=-999,'logLik'=-999,'length'=-999,'Y1'=NA,'Y2'=NA)
			tx = 1
			### quality filtering step ###
			if ( (length(which(!is.na(vdat))) > thresh.na) & (length(which(vdat[!is.na(vdat)] > 0)) > thresh.0) ){
			### quality filtering step ###
				years.stretch = seq(dat1$Year1[1],dat1$Year1[nrow(dat1)],1) #expand time series to include missing years
				y.match = match(years.stretch,dat1$Year1)
				X = y.match
				X[which(!is.na(X))] = vdat[X[which(!is.na(X))]] #fill-in values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(years.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale) #original transform: scale between 0 and 1
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				add.coef.Z = c(arr.Z[[1]],arr.Z[[2]],arr.Z[[3]][1,1],arr.Z[[3]][2,1],arr.Z[[5]][1],arr.Z[[5]][2],arr.Z[[6]])
				trends[tx,1] = unique(dat1$LTER.site) #LTER
				trends[tx,2] = uL[l] #Site
				trends[tx,3] = add.coef.Z[1] #MSE
				trends[tx,4] = add.coef.Z[2] #b
				trends[tx,5] = add.coef.Z[3] #intercept
				trends[tx,6] = add.coef.Z[4] #slope
				trends[tx,7] = add.coef.Z[5] #pval intercept
				trends[tx,8] = add.coef.Z[6] #pval slope
				trends[tx,9] = add.coef.Z[7] #logLik
				trends[tx,10] = length(which(!is.na(Z))) #length
				trends[tx,11] = years.stretch[1] #first year
				trends[tx,12] = years.stretch[length(years.stretch)] #last year
				tx = tx + 1				
				# Plot Z transform on t.scale fit
				png(paste0('./plots/time_trends/diversity/fitted/',variables[v],'_',unique(dat1$LTER.site),'_',uL[l],'.png'))
				plot(Z ~ t.scale, main=paste0(variables[v], ": b = ",round(arr.Z$b, digits=3)))
				curve(arr.Z$coef[1] + arr.Z$coef[2] * x, from=0, to=max(t.scale), add=T)
				dev.off()
				#Plot "high quality" time series
				png(paste0('./plots/time_trends/diversity/lineplots/',variables[v],'_',unique(dat1$LTER.site),'_',uL[l],'.png'))
				plot(t.scale,Z,type='l',main=paste0(unique(dat1$LTER.site),'_',uL[l]),xlab='Scaled time',ylab=variables[v],lwd=2)
				abline(a=add.coef.Z[3],b=add.coef.Z[4],lty=2,col='red',lwd=1.5)
				legend('topright',legend=c(paste0('slope = ',round(add.coef.Z[4],2)),paste0('length = ',length(which(!is.na(Z)))),paste0('autocor = ',round(add.coef.Z[2],2)),paste0('logLik = ',round(add.coef.Z[7],2))),ncol=1)
				dev.off()
				#Write output file
				write.csv(trends,paste0('./diversity_trends/',variables[v],'_',unique(dat1$LTER.site),'_',uL[l],'.csv'),quote=F,row.names=F)
			} else {
				lost.sites = c(lost.sites,files2[f])
			}
		}
	}
}
write.csv(all.metrics,'arthropods_BETAdiversity.csv',quote=F,row.names=F)

# Compare relationships among diversity measures
files = list.files('./diversity_trends',full.names=T)
files = files[grep('Beta',files)]
files.split = lapply(files,function(x){strsplit(gsub('.csv','',strsplit(x,'/')[[1]][3]),'_')[[1]]})
LL = unlist(lapply(files.split,function(x){paste(x[2],x[3],sep='_')}))
uLL = unique(LL)
trends = data.frame('LTER.Locale'=uLL,'Beta.2'=rep(NA,length(uLL)),'Beta.j'=rep(NA,length(uLL)),'Beta.bray'=rep(NA,length(uLL)))
for (f in 1:length(files)){
	split1 = strsplit(gsub('.csv','',strsplit(files[f],'/')[[1]][3]),'_')[[1]]
	LL1 = paste(split1[2],split1[3],sep='_')
	dat1 = read.csv(files[f],as.is=T,check.names=F,header=T)
	trends[which(uLL==LL1),which(variables==split1[1])+1] = dat1$coef_slope
}
write.csv(trends,'time_trends_arthropods_BETAdiversity.csv',quote=F,row.names=F)

trends = read.csv('time_trends_arthropods_BETAdiversity.csv',as.is=T,check.names=F,header=T)

# Relationships among beta metrics
trends$L = as.factor(apply(array(trends$LTER.Locale),1,function(x){strsplit(x,'_')[[1]][1]}))
library('RColorBrewer')
cols = brewer.pal(9,'Paired')
png('./plots/betadiversity_metrics_change_comparison.png',width=480*2,height=480*2)
par(mfrow=c(2,2),cex=2,cex.lab=1.5,cex.axis=1.2,lwd=2,oma=c(0,0,0,0),mar=c(5,5,0,0))
plot(trends$Beta.2,trends$Beta.j,col=cols[trends$L],pch=16,xlab='Beta.2',ylab='Beta.j',cex=1.5)
abline(a=0,b=1,lty=2,col='gray50')
plot(trends$Beta.2,trends$Beta.bray,col=cols[trends$L],pch=16,xlab='Beta.2',ylab='Beta.bray',cex=1.5)
abline(a=0,b=1,lty=2,col='gray50')
plot(trends$Beta.j,trends$Beta.bray,col=cols[trends$L],pch=16,xlab='Beta.j',ylab='Beta.bray',cex=1.5)
abline(a=0,b=1,lty=2,col='gray50')
plot(trends$Beta.j,trends$Beta.bray,type='n',axes=F,ann=F)
legend('topleft',legend=levels(trends$L),col=cols,pch=16,cex=1.2,bty='n')
dev.off()

# Group studies by site/taxa/method
trends$LTER = apply(array(trends$LTER.Locale),1,function(x){strsplit(x,'_')[[1]][1]})
lters = sort(unique(trends$LTER))
trends$Site = apply(array(trends$LTER.Locale),1,function(x){strsplit(x,'_')[[1]][2]})
trends$Site[which(trends$LTER=='MidwestSTN')] = 'aphids' #merge aphid records into one boxplot
trends$Site[which(trends$LTER=='HarvardForest' & (trends$Site=='ants' | trends$Site=='ants.pitfall' | trends$Site=='ants.hand' | trends$Site=='ants.Nantucket' | trends$Site=='ants.bait' | trends$Site=='ants.litter'))] = 'ants' #merge ant records into one boxplot
trends$Site[which(trends$LTER=='NorthTemperateLakes' & (trends$Site=='CR' | trends$Site=='SP' | trends$Site=='TR'))] = 'lake insects' #merge insect records into one boxplot
trends$Site[which(trends$LTER=='KonzaPrairie' & trends$Site=='gall')] = 'gall insects' 
trends$Site[which(trends$LTER=='NorthTemperateLakes' & trends$Site=='Crayfish')] = 'crayfish' 
trends$Site[which(trends$LTER=='CentralArizona-Phoenix' & (trends$Site=='pitfall1' | trends$Site=='pitfall2'))] = 'pitfall' #merge pitfall records
trends$Site[which(trends$Site=='grasshopper')] = 'grasshoppers'
trends$Site[which(trends$Site=='aquatic')] = 'stream insects'
trends$Site[grep('sweep',trends$Site)] = 'sweep net'
trends$Site[which(trends$Site=='mosquito')] = 'mosquitoes'
trends$Site[which(trends$LTER=='CedarCreek' & (trends$Site=='sweep1' | trends$Site=='sweep2'))] = 'sweep nets' 

# Boxplots
trends$LL = as.factor(paste(trends$LTER,trends$Site,sep='_'))
levels(trends$LL)
variables = c('Beta.2','Beta.j','Beta.bray')
var2 = c('beta 2','beta (Jaccard)','beta (Bray-Curtis)')
for (v in 1:length(variables)){
	png(paste0('./plots/time_trends/diversity/',variables[v],'_trends_boxplot.png'),width=600,height=500)
	par(mfrow=c(1,1),cex.lab=2,cex.axis=1.5,lwd=2,oma=c(0,0,0,0),mar=c(15,8,1,1))
	plot(trends$LL,trends[,which(colnames(trends)==variables[v])],ylim=c(-15,15),xaxt='n',yaxt='n',axes=F,ylab=paste0('Change in ',var2[v]))
	axis(2,lwd=2,at=seq(-15,15,5),labels=seq(-15,15,5))
	text(x=c(1,2,3.5,5.5,7,8,9,10,11.5),y=par("usr")[3]-0.2,labels=lters,srt=45,pos=2,xpd=T,cex=2)
	text(x=seq(1,length(levels(trends$LL)),1)-0.35,y=par("usr")[3]+0.5,labels=apply(array(levels(trends$LL)),1,function(x){strsplit(x,'_')[[1]][2]}),srt=90,adj=0,xpd=T,cex=1.2)
	abline(a=0,b=0,lty=1,col='black',lwd=2)
	dev.off()	
}

# One Sample T-test
uL = unique(trends$LTER)
means.j = apply(array(uL),1,function(x){mean(trends$Beta.j[which(trends$LTER==x)],na.rm=T)})
t.j = t.test(as.numeric(means.j))
means.2 = apply(array(uL),1,function(x){mean(trends$Beta.2[which(trends$LTER==x)],na.rm=T)})
t.2 = t.test(as.numeric(means.2))
means.bray = apply(array(uL),1,function(x){mean(trends$Beta.bray[which(trends$LTER==x)],na.rm=T)})
t.bray = t.test(as.numeric(means.bray))

# Barplots
# Beta j
mean1 = t.j[['estimate']]
uppers = t.j[['conf.int']][2]
lowers = t.j[['conf.int']][1]
png('./plots/t-test_Beta.j_barplot.png',width=200,height=500)
par(cex.lab=2,cex.axis=1.5,oma=c(0,0,0,0),mar=c(15,8,2,2),lwd=2)
bars = barplot(mean1,ylim=c(-4,4),names.arg='',ylab='Average change in Beta.j',main=NULL,lwd=2)
abline(h=0)
axis(2,lwd=2)
text(x=bars,y=par("usr")[3]-0.05,srt=90,adj=1,labels='LTER',xpd=T,cex=2)
segments(bars,lowers,bars,uppers,lwd=2)
arrows(bars,lowers,bars,uppers,lwd=2,angle=90,code=3,length=(par("usr")[2]-par("usr")[1])/4)
dev.off()


#################################################Section
# One Sample T-tests of abundance trends

trends = read.csv('time_trends_arthropods_relaxed.csv',as.is=T,check.names=F,header=T)
trends$LTER.site = paste(trends$LTER,trends$Site,sep='_')

# t-test of means per LTER
uL = unique(trends$LTER)
means2 = apply(array(uL),1,function(x){mean(trends$coef_slope[which(trends$LTER==x)],na.rm=T)})
t2 = t.test(as.numeric(means2))

# Barplot
mean2 = t2[['estimate']]
uppers = t2[['conf.int']][2]
lowers = t2[['conf.int']][1]
png('./plots/t-test_abundance_lters_barplot.png',width=200,height=500)
par(cex.lab=2,cex.axis=1.5,oma=c(0,0,0,0),mar=c(15,8,2,2),lwd=2)
bars = barplot(mean2,ylim=c(-1.5,1.5),names.arg='',ylab='Average change in abundance',main=NULL,lwd=2)
abline(h=0)
axis(2,lwd=2)
text(x=bars,y=par("usr")[3]-0.05,srt=90,adj=1,labels='LTER',xpd=T,cex=2)
segments(bars,lowers,bars,uppers,lwd=2)
arrows(bars,lowers,bars,uppers,lwd=2,angle=90,code=3,length=(par("usr")[2]-par("usr")[1])/15)
dev.off()


# t-test of means per LTER in early vs. late time bins
trends$Site[which(trends$LTER=='MidwestSTN')] = 'aphids' #merge Kellogg records into one boxplot
trends$Site[which(trends$LTER=='HubbardBrook')] = 'butterflies/moths' #merge lepidoptera records into one boxplot
trends$Site[which(trends$LTER=='HarvardForest' & trends$Site!='ants' & trends$Site!='ants.pitfall' & trends$Site!='ants.bait' & trends$Site!='ants.hand' & trends$Site!='ants.litter')] = 'ticks' #merge tick records into one boxplot
trends$Site[grep('ants',trends$Site)] = 'ants' #merge tick records into one boxplot
trends$Site[which(trends$LTER=='NorthTemperateLakes' & trends$Site!='Crayfish')] = 'lake insects' #merge insect records into one boxplot
trends$Site[which(trends$LTER=='NorthTemperateLakes' & trends$Site=='Crayfish')] = 'crayfish' 
trends$Site[which(trends$LTER=='KonzaPrairie' & trends$Site=='gall')] = 'gall insects' 
trends$Site[which(trends$LTER=='GeorgiaCoastal' & (trends$Site=='Burrowing Crab' | trends$Site=='Fiddler Crab'))] = 'crabs'  
trends$Site[which(trends$LTER=='GeorgiaCoastal' & trends$Site=='Salt Marsh' & trends$Species=='Acrididae')] = 'grasshoppers'  
trends$Site[which(trends$LTER=='GeorgiaCoastal' & trends$Site=='Salt Marsh' & trends$Species=='Prokelisia marginata')] = 'planthoppers'  
trends$Site[which(trends$LTER=='CentralArizona-Phoenix' & (trends$Site=='pitfall1' | trends$Site=='pitfall2'))] = 'pitfalls' #merge pitfall records
trends$Site[which(trends$Site=='grasshopper')] = 'grasshoppers'
trends$Site[which(trends$Site=='aquatic')] = 'stream insects'
trends$Site[which(trends$Site=='leafminer')] = 'leaf miners'
trends$Site[which(trends$Site=='bark beetle')] = 'bark beetles'
trends$Site[which(trends$Site=='mosquito')] = 'mosquitoes'
trends$Site[grep('sweep',trends$Site)] = 'sweep net'
trends$Site[which(trends$LTER=='CedarCreek' & (trends$Site=='sweep1' | trends$Site=='sweep2'))] = 'sweep nets' 
trends$LTER.site = paste(trends$LTER,trends$Site,sep='_')
uLS = unique(trends$LTER.site)
year.ranges = apply(array(uLS),1,function(x){y=trends[which(trends$LTER.site==x),12:13];unique(paste(y[,1],y[,2],sep='_'))})
#take maximum range for combined sites
for (i in 1:length(year.ranges)){
	ranges1 = year.ranges[[i]]
	rsplit = t(apply(array(ranges1),1,function(x){strsplit(x,'_')[[1]]}))
	ymax = max(as.numeric(rsplit[,2]))
	ymin = min(as.numeric(rsplit[,1]))
	year.ranges[[i]] = paste(ymin,ymax,sep='_')
}
year.ranges = unlist(year.ranges)
yr = t(apply(array(year.ranges),1,function(x){as.numeric(strsplit(x,'_')[[1]])}))
early = which(yr[,1]<1990)
middle1 = which(yr[,1]>=1990 & yr[,1]<2000)
late = which(yr[,1]>=2000 & yr[,1]<2010)
late2 = which(yr[,1]>=2010)
means = apply(array(uLS),1,function(x){mean(trends$coef_slope[which(trends$LTER.site==x)],na.rm=T)})
t.early = t.test(as.numeric(means[early]))
t.middle1 = t.test(as.numeric(means[middle1]))
t.late = t.test(as.numeric(means[late]))
t.late2 = t.test(as.numeric(means[late2]))
mean1 = t.early[['estimate']]
uppers1 = t.early[['conf.int']][2]
lowers1 = t.early[['conf.int']][1]
mean2 = t.middle1[['estimate']]
uppers2 = t.middle1[['conf.int']][2]
lowers2 = t.middle1[['conf.int']][1]
mean3 = t.late[['estimate']]
uppers3 = t.late[['conf.int']][2]
lowers3 = t.late[['conf.int']][1]
mean4 = t.late2[['estimate']]
uppers4 = t.late2[['conf.int']][2]
lowers4 = t.late2[['conf.int']][1]
lowers = c(lowers1,lowers2,lowers3,lowers4)
uppers = c(uppers1,uppers2,uppers3,uppers4)
png('./plots/t-test_abundance_time_barplot.png',width=400,height=500)
par(cex.lab=2,cex.axis=1.5,oma=c(0,0,0,0),mar=c(15,8,3,2),lwd=2)
bars = barplot(c(mean1,mean2,mean3,mean4),ylim=c(-2,2),names.arg='',ylab='Average change in abundance',main=NULL,lwd=2)
abline(h=0)
axis(2,lwd=2)
text(x=bars,y=par("usr")[3]-0.05,srt=45,adj=1,labels=c('<1990','1990-2000','2000-2010','>2010'),xpd=T,cex=1.8)
title(xlab='Start year',line=7,cex.lab=2)
segments(bars,lowers,bars,uppers,lwd=2)
arrows(bars,lowers,bars,uppers,lwd=2,angle=90,code=3,length=(par("usr")[2]-par("usr")[1])/15)
text(x=bars,y=uppers+0.2,labels=paste('n',c(length(early),length(middle1),length(late),length(late2)),sep='='),cex=1.5)
dev.off()

# Summary statistics from time-binned t-tests
tlist = list(t.early,t.middle1,t.late,t.late2)
lapply(tlist,function(x){c(x$estimate,x$parameter,x$statistic,x$p.value)})


# t-test of means per LTER in early vs. late time bins (USING END YEAR)
early = which(yr[,2]<1990)
middle1 = which(yr[,2]>=1990 & yr[,2]<2000)
late = which(yr[,2]>=2000 & yr[,2]<2010)
late2 = which(yr[,2]>=2010)
means = apply(array(uLS),1,function(x){mean(trends$coef_slope[which(trends$LTER.site==x)],na.rm=T)})
#t.early = t.test(as.numeric(means[early]))
t.middle1 = t.test(as.numeric(means[middle1]))
t.late = t.test(as.numeric(means[late]))
t.late2 = t.test(as.numeric(means[late2]))
#mean1 = t.early[['estimate']]
#uppers1 = t.early[['conf.int']][2]
#lowers1 = t.early[['conf.int']][1]
mean2 = t.middle1[['estimate']]
uppers2 = t.middle1[['conf.int']][2]
lowers2 = t.middle1[['conf.int']][1]
mean3 = t.late[['estimate']]
uppers3 = t.late[['conf.int']][2]
lowers3 = t.late[['conf.int']][1]
mean4 = t.late2[['estimate']]
uppers4 = t.late2[['conf.int']][2]
lowers4 = t.late2[['conf.int']][1]
lowers = c(lowers2,lowers3,lowers4)
uppers = c(uppers2,uppers3,uppers4)
png('./plots/t-test_abundance_time_barplot_ENDYEAR.png',width=400,height=500)
par(cex.lab=2,cex.axis=1.5,oma=c(0,0,0,0),mar=c(15,8,3,2),lwd=2)
bars = barplot(c(mean2,mean3,mean4),ylim=c(-2,2),names.arg='',ylab='Average change in abundance',main=NULL,lwd=2)
abline(h=0)
axis(2,lwd=2)
text(x=bars,y=par("usr")[3]-0.05,srt=45,adj=1,labels=c('1990-2000','2000-2010','>2010'),xpd=T,cex=1.8)
title(xlab='End year',line=7,cex.lab=2)
segments(bars,lowers,bars,uppers,lwd=2)
arrows(bars,lowers,bars,uppers,lwd=2,angle=90,code=3,length=(par("usr")[2]-par("usr")[1])/15)
text(x=bars,y=uppers+0.2,labels=paste('n',c(length(middle1),length(late),length(late2)),sep='='),cex=1.5)
dev.off()

# Summary statistics from time-binned t-tests
tlist = list(t.middle1,t.late,t.late2)
lapply(tlist,function(x){c(x$estimate,x$parameter,x$statistic,x$p.value)})


# Trends among feeding guilds
# Herbivore
trends1 = trends[which(trends$Feeding=='herbivore'),]
uL = unique(trends1$LTER)
means = apply(array(uL),1,function(x){mean(trends1$coef_slope[which(trends1$LTER==x)],na.rm=T)})
t1 = t.test(as.numeric(means))
# Carnivore
trends2 = trends[which(trends$Feeding=='carnivore'),]
uL = unique(trends2$LTER)
means = apply(array(uL),1,function(x){mean(trends2$coef_slope[which(trends2$LTER==x)],na.rm=T)})
t2 = t.test(as.numeric(means))
# Omnivore
trends3 = trends[which(trends$Feeding=='omnivore'),]
uL = unique(trends3$LTER)
means = apply(array(uL),1,function(x){mean(trends3$coef_slope[which(trends3$LTER==x)],na.rm=T)})
t3 = t.test(as.numeric(means))
# Detritivore
trends4 = trends[which(trends$Feeding=='detritivore'),]
uL = unique(trends4$LTER)
means = apply(array(uL),1,function(x){mean(trends4$coef_slope[which(trends4$LTER==x)],na.rm=T)})
t4 = t.test(as.numeric(means))
# Parasite
trends5 = trends[which(trends$Feeding=='parasite'),]
uL = unique(trends5$LTER)
means = apply(array(uL),1,function(x){mean(trends5$coef_slope[which(trends5$LTER==x)],na.rm=T)})
t5 = t.test(as.numeric(means))
# Parasitoid
trends6 = trends[which(trends$Feeding=='parasitoid'),]
uL = unique(trends6$LTER)
means = apply(array(uL),1,function(x){mean(trends6$coef_slope[which(trends6$LTER==x)],na.rm=T)})
t6 = t.test(as.numeric(means))

# Barplot
tlist = list(t1,t2,t3,t4,t5,t6)
names1 = c('Herbivore','Carnivore','Omnivore','Detritivore','Parasite','Parasitoid')

for (i in 1:length(tlist)){
	t0 = tlist[[i]]
	mean1 = t0[['estimate']]
	uppers = t0[['conf.int']][2]
	lowers = t0[['conf.int']][1]
	png(paste0('./plots/t-test_abundance_lters_',names1[i],'_barplot.png'),width=200,height=500)
	par(cex.lab=2,cex.axis=1.5,oma=c(0,0,0,0),mar=c(15,8,2,2),lwd=2)
	if (i==6){
		bars = barplot(mean1,ylim=c(-6,6),names.arg='',ylab='Average change in abundance',main=NULL,lwd=2)
	} else {
		bars = barplot(mean1,ylim=c(-2,2),names.arg='',ylab='Average change in abundance',main=NULL,lwd=2)
	}
	abline(h=0)
	axis(2,lwd=2)
	text(x=bars,y=par("usr")[3]-0.05,srt=90,adj=1,labels='LTER/site',xpd=T,cex=2)
	segments(bars,lowers,bars,uppers,lwd=2)
	arrows(bars,lowers,bars,uppers,lwd=2,angle=90,code=3,length=(par("usr")[2]-par("usr")[1])/5)
	dev.off()
}

# Summary statistics from feeding guild t-tests
lapply(tlist,function(x){c(x$estimate,x$parameter,x$statistic,x$p.value)})


################################################Section
# Curate bird and fish data


# Arctic - fish
### WARNING: data records mass & length of fish, nothing about species or number. I calculated abundance from the number of fish measurements per year
### Data excluded


# Baltimore - birds
site.name = 'Baltimore'
locale.name = 'birds'
data1 = read.csv('./raw_data/bird-survey-2001-2015-birds.csv',as.is=T,check.names=F,header=T)
dates.key = read.csv('./raw_data/bird-survey-2001-2015-surveys.csv',as.is=T,check.names=F,header=T)
dates.key = data.frame('survey_id'=unique(dates.key$survey_id),'Year'=apply(array(unique(dates.key$survey_id)),1,function(x){y=unique(dates.key$survey_date[which(dates.key$survey_id==x)]);strsplit(y,'/')[[1]][3]}),stringsAsFactors=F)
data1$Year = apply(array(data1$survey_id),1,function(x){dates.key$Year[which(dates.key$survey_id==x)]})
data1$Site = data1$site_id
data1$Species = data1$species_id
data1$Number = data1$bird_count
data1$Number[which(is.na(data1$Number))] = 0
u.years = sort(as.numeric(unique(data1$Year)))
u.species = unique(data1$Species)
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
		n.y1 = 1
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.slopes,paste0('./summary_tables/vertebrates/',site.name,'_Bird_SpeciesSlopes.txt'),sep=';',quote=F,row.names=F)


# Central Arizona-Phoenix - birds 1
site.name = 'CentralArizona-Phoenix'
locale.name = 'birds1'
data1 = read.csv('./raw_data/46_core_birds_c217ccdaab8d94e2d6065aa301d8ea7e.csv',as.is=T,check.names=F,header=T)
data1$Year = apply(array(data1$survey_date),1,function(x){strsplit(x,'-')[[1]][1]})
data1$Species = data1$common_name
data1$Number = data1$bird_count
u.years = sort(as.numeric(unique(data1$Year)))
u.species = unique(data1$Species)
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
		n.y1 = 1
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.slopes,paste0('./summary_tables/vertebrates/',site.name,'_Bird1_SpeciesSlopes.txt'),sep=';',quote=F,row.names=F)


# Central Arizona-Phoenix - birds 2
site.name = 'CentralArizona-Phoenix'
locale.name = 'birds2'
data1 = read.csv('./raw_data/641_srbp_birds_7ffe7d716e64a4def128217851d6f2d0.csv',as.is=T,check.names=F,header=T)
str(data1)
data1$Year = apply(array(data1$survey_date),1,function(x){strsplit(x,'-')[[1]][1]})
data1$Site = data1$site_code
data1$Species = data1$common_name
data1$Number = data1$bird_count
u.species = unique(data1$Species)
u.years = sort(as.numeric(unique(data1$Year)))
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
		n.y1 = 1 #length(unique(paste(y.data$reach,y.data$survey_date,sep='_')))
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.slopes,paste0('./summary_tables/vertebrates/',site.name,'_Bird2_SpeciesSlopes.txt'),sep=';',quote=F,row.names=F)


# Coweeta - birds 1
site.name = 'Coweeta'
locale.name = 'birds1'
data1 = read.table('./raw_data/4046_bird_detections.csv',sep='\t',as.is=T,check.names=F,header=T)
str(data1)
data1$Year = data1$year
data1$Site = data1$point_num
data1$Species = data1$bird_name
data1$Number = apply(data1[,7:35],1,function(x){sum(x,na.rm=T)})
u.years = sort(as.numeric(unique(data1$Year)))
u.species = unique(data1$Species)
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
		n.y1 = length(unique(y.data$visit))
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.slopes,paste0('./summary_tables/vertebrates/',site.name,'_Bird1_SpeciesSlopes.txt'),sep=';',quote=F,row.names=F)


# Coweeta - fish and herptiles
 #has arthropods too
site.name = 'Coweeta'
locale.name = 'fish-herps'
data1 = read.csv('./raw_data/LTR_REG_4065_1_0.CSV',as.is=T,check.names=F,header=T)
data1 = data1[-c(1:2),]
str(data1)
data1$Site = data1$Site_Name
u.years = sort(as.numeric(unique(data1$Year)))
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 14:41){
	s.data = data1[,c(1:7,s)]
	s.data$Number = as.numeric(s.data[,8])
	for (y in 1:length(u.years)){
		y.data = s.data[which(s.data$Year==u.years[y]),]
		n.y1 = length(unique(paste(y.data$Month,y.data$Day,y.data$Plot_No,sep='_')))
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = colnames(data1)[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.slopes,paste0('./summary_tables/vertebrates/',site.name,'_FishHerp_SpeciesSlopes.txt'),sep=';',quote=F,row.names=F)


# Coweeta - fish 1
site.name = 'Coweeta'
locale.name = 'fish1'
data1 = read.csv('./raw_data/3047_fish_data_1_0.CSV',as.is=T,check.names=F,header=T)
str(data1)
data1$Site = data1$Creek_location
data1$Species = data1$Fish_species
u.years = sort(as.numeric(unique(data1$Year)))
u.species = unique(data1$Species)
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
		n.y1 = length(unique(paste0(y.data$Month,y.data$Day,sep='_')))
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = nrow(y.data) #abundance
		cx = cx + 1
	}
}
write.table(out.slopes,paste0('./summary_tables/vertebrates/',site.name,'_Fish1_SpeciesSlopes.txt'),sep=';',quote=F,row.names=F)


# Coweeta - fish 2
site.name = 'Coweeta'
locale.name = 'fish2'
data1 = read.csv('./raw_data/CWT_POP_3006_1_0.CSV',as.is=T,check.names=F,header=T)
str(data1)
data1$Year = paste0('19',data1$Year)
data1$Site = data1$Site
data1$Species = data1$Species
u.years = sort(as.numeric(unique(data1$Year)))
u.species = unique(data1$Species)
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
		n.y1 = length(unique(paste0(y.data$Month,y.data$Day,sep='_')))
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = nrow(y.data) #abundance
		cx = cx + 1
	}
}
write.table(out.slopes,paste0('./summary_tables/vertebrates/',site.name,'_Fish2_SpeciesSlopes.txt'),sep=';',quote=F,row.names=F)


# Hubbard Brook - birds 1
site.name = 'HubbardBrook'
locale.name = 'birds1'
data1 = read.csv('./raw_data/knb-lter-hbr.81.7/hb_bird.txt',as.is=T,check.names=F,header=T)
data2 = data.frame('Year'=NA,'Species'=NA,'Number'=-999)
u.species = unique(data1[,1])
years = colnames(data1)[-1]
ix=1
for (y in 1:length(years)){
	for (s in 1:length(u.species)){
		n = data1[which(data1[,1]==u.species[s]),which(colnames(data1)==years[y])]
		if (n=='t'){n=0}else{}
		data2[ix,1] = years[y]
		data2[ix,2] = u.species[s]
		data2[ix,3] = ceiling(as.numeric(n))
		ix=ix+1
	}
}
str(data2)
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data2[which(data2$Year==u.years[y] & data2$Species==u.species[s]),]
		n.y1 = 1
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.slopes,paste0('./summary_tables/vertebrates/',site.name,'_Bird1_SpeciesSlopes.txt'),sep=';',quote=F,row.names=F)


# Hubbard Brook - birds 2
site.name = 'HubbardBrook'
locale.name = 'birds2'
data1 = read.csv('./raw_data/knb-lter-hbr.81.7/mk_bird.txt',as.is=T,check.names=F,header=T)
data2 = data.frame('Year'=NA,'Species'=NA,'Number'=-999)
u.species = unique(data1[,1])
u.species = u.species[-c(35:36)]
years = colnames(data1)[-1]
ix=1
for (y in 1:length(years)){
	for (s in 1:length(u.species)){
		n = data1[which(data1[,1]==u.species[s]),which(colnames(data1)==years[y])]
		if (n=='t'){n=0}else{}
		data2[ix,1] = years[y]
		data2[ix,2] = u.species[s]
		data2[ix,3] = ceiling(as.numeric(n))
		ix=ix+1
	}
}
str(data2)
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data2[which(data2$Year==u.years[y] & data2$Species==u.species[s]),]
		n.y1 = 1
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.slopes,paste0('./summary_tables/vertebrates/',site.name,'_Bird2_SpeciesSlopes.txt'),sep=';',quote=F,row.names=F)


# Hubbard Brook - birds 3
site.name = 'HubbardBrook'
locale.name = 'birds3'
data1 = read.csv('./raw_data/knb-lter-hbr.81.7/rp_bird.txt',as.is=T,check.names=F,header=T)
data2 = data.frame('Year'=NA,'Species'=NA,'Number'=-999)
u.species = unique(data1[,1])
u.species = u.species[-c(35:36)]
years = colnames(data1)[-1]
ix=1
for (y in 1:length(years)){
	for (s in 1:length(u.species)){
		n = data1[which(data1[,1]==u.species[s]),which(colnames(data1)==years[y])]
		if (n=='t'){n=0}else{}
		data2[ix,1] = years[y]
		data2[ix,2] = u.species[s]
		data2[ix,3] = ceiling(as.numeric(n))
		ix=ix+1
	}
}
str(data2)
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data2[which(data2$Year==u.years[y] & data2$Species==u.species[s]),]
		n.y1 = 1
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.slopes,paste0('./summary_tables/vertebrates/',site.name,'_Bird3_SpeciesSlopes.txt'),sep=';',quote=F,row.names=F)


# Hubbard Brook - birds 4
site.name = 'HubbardBrook'
locale.name = 'birds4'
data1 = read.csv('./raw_data/knb-lter-hbr.81.7/sm_bird.txt',as.is=T,check.names=F,header=T)
data2 = data.frame('Year'=NA,'Species'=NA,'Number'=-999)
u.species = unique(data1[,1])
u.species = u.species[-c(36:37)]
years = colnames(data1)[-1]
ix=1
for (y in 1:length(years)){
	for (s in 1:length(u.species)){
		n = data1[which(data1[,1]==u.species[s]),which(colnames(data1)==years[y])]
		if (n=='t'){n=0}else{}
		data2[ix,1] = years[y]
		data2[ix,2] = u.species[s]
		data2[ix,3] = as.numeric(n)
		ix=ix+1
	}
}
str(data2)
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data2[which(data2$Year==u.years[y] & data2$Species==u.species[s]),]
		n.y1 = 1
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.slopes,paste0('./summary_tables/vertebrates/',site.name,'_Bird4_SpeciesSlopes.txt'),sep=';',quote=F,row.names=F)


# Konza Prairie - birds 1
site.name = 'KonzaPrairie'
locale.name = 'birds1'
data1 = read.csv('./raw_data/CBC011.csv',as.is=T,check.names=F,header=T)
data1$Year = data1$RecYear
u.species = unique(data1$Species)
u.years = sort(as.numeric(unique(data1$Year)))
data2 = data.frame('Year'=NA,'Species'=NA,'Number'=-999)
ix=1
for (y in 1:length(u.years)){
	for (s in 1:length(u.species)){
		n = data1[which(data1$Species==u.species[s] & data1$Year==years[y]),-c(1:4,57)]
		if (nrow(n)==0){
		} else {
			n[which(n=='x')] = 1 #x means bird was recorded as present
			n[which(n=='n')] = 1 #n means bird was recorded as nesting
			data2[ix,1] = years[y]
			data2[ix,2] = u.species[s]
			data2[ix,3] = mean(as.numeric(n),na.rm=T)
			ix=ix+1
		}
	}
}
str(data2)
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data2[which(data2$Year==u.years[y] & data2$Species==u.species[s]),]
		n.y1 = 1
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.slopes,paste0('./summary_tables/vertebrates/',site.name,'_Bird1_SpeciesSlopes.txt'),sep=';',quote=F,row.names=F)


# Konza Prairie - birds 2
site.name = 'KonzaPrairie'
locale.name = 'birds2'
data1 = read.csv('./raw_data/CBP011_0.csv',as.is=T,check.names=F,header=T)
str(data1)
data1$Species = data1$COMMONNAME
data1$Site = data1$WATERSHED
data1$Number = data1$COUNT
data1$Year = data1$RECYEAR
u.species = unique(data1$Species)
u.years = sort(as.numeric(unique(data1$Year)))
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
		n.y1 = length(unique(paste(y.data$RECMONTH,y.data$RECDAY,y.data$TRANSNUM,sep='_')))
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.slopes,paste0('./summary_tables/vertebrates/',site.name,'_Bird2_SpeciesSlopes.txt'),sep=';',quote=F,row.names=F)


# Konza Prairie - birds 3
# "Sparrow" counts
site.name = 'KonzaPrairie'
locale.name = 'sparrows'
data1 = read.csv('./raw_data/CBS031.csv',as.is=T,check.names=F,header=T)
str(data1)
data1$Site = data1$Watershed
data1$Year = apply(array(data1$SurveyDate),1,function(x){strsplit(x,'/')[[1]][3]})
u.years = sort(as.numeric(unique(data1$Year)))
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[y]),]
	n.y1 = 1
	out.slopes[cx,1] = site.name #LTER.site
	out.slopes[cx,2] = locale.name #Locale
	out.slopes[cx,3] = 'Ammodramus sp.' #Species.code
	out.slopes[cx,4] = u.years[y] #Year
	out.slopes[cx,5] = n.y1 #Number of observations
	out.slopes[cx,6] = nrow(y.data) #abundance
	cx = cx + 1
}
write.table(out.slopes,paste0('./summary_tables/vertebrates/',site.name,'_Bird3_SpeciesSlopes.txt'),sep=';',quote=F,row.names=F)


# Konza Prairie - birds 4
# Greater Prairie Chicken count
site.name = 'KonzaPrairie'
locale.name = 'prairie chickens'
data1 = read.csv('./raw_data/CPC011.csv',as.is=T,check.names=F,header=T)
str(data1)
data1$Year = data1$RECYEAR
data1$Number = data1$NUMBIRDS
u.years = sort(as.numeric(unique(data1$Year)))
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (y in 1:length(u.years)){
	y.data = data1[which(data1$Year==u.years[y]),]
	n.y1 = 1
	out.slopes[cx,1] = site.name #LTER.site
	out.slopes[cx,2] = locale.name #Locale
	out.slopes[cx,3] = 'Tympanuchus cupido' #Species.code
	out.slopes[cx,4] = u.years[y] #Year
	out.slopes[cx,5] = n.y1 #Number of observations
	out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
	cx = cx + 1
}
write.table(out.slopes,paste0('./summary_tables/vertebrates/',site.name,'_Bird4_SpeciesSlopes.txt'),sep=';',quote=F,row.names=F)


# Konza Prairie - fish
#has arthropods too
site.name = 'KonzaPrairie'
locale.name = 'fish'
data1 = read.csv('./raw_data/CFC011.csv',as.is=T,check.names=F,header=T)
str(data1)
data1$Site = data1$Watershed
data1$Number = data1$Count
data1$Year = data1$RecYear
data1$Site.Method = paste(data1$Site,data1$Method,sep='_')
u.years = sort(as.numeric(unique(data1$Year)))
u.species = unique(data1$Species)
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
		n.y1 = length(unique(paste(y.data$RecMonth,y.data$RecDay,y.data$Replicat,sep='_')))
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
		cx = cx + 1
	}
}
write.table(out.slopes,paste0('./summary_tables/vertebrates/',site.name,'_Fish_SpeciesSlopes.txt'),sep=';',quote=F,row.names=F)


# North Temperate Lakes - fish 1
site.name = 'NorthTemperateLakes'
locale.name = 'fish1'
data1 = read.csv('./raw_data/WifishAbundance.csv',as.is=T,check.names=F,header=T)
str(data1)
data1$Site = data1$GEARNAME
data1$Number = data1$CPUE
data1$Year = data1$YEAR
data1$Species = data1$taxon_id
u.years = sort(as.numeric(unique(data1$Year)))
u.years = u.years[which(u.years>1800)]
u.species = unique(data1$Species)
method = unique(data1$Site)
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (l in 1:length(method)){
	l.data = data1[which(data1$Site==method[l]),]
	for (s in 1:length(u.species)){
		for (y in 1:length(u.years)){
			y.data = l.data[which(l.data$Year==u.years[y] & l.data$Species==u.species[s]),]
			n.y1 = sum(y.data$N,na.rm=T)
			out.slopes[cx,1] = site.name #LTER.site
			out.slopes[cx,2] = paste(locale.name,method[l],sep='-') #Locale
			out.slopes[cx,3] = u.species[s] #Species.code
			out.slopes[cx,4] = u.years[y] #Year
			out.slopes[cx,5] = n.y1 #Number of observations
			out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
			cx = cx + 1
		}
	}
}
write.table(out.slopes,paste0('./summary_tables/vertebrates/',site.name,'_Fish1_SpeciesSlopes.txt'),sep=';',quote=F,row.names=F)


# North Temperate Lakes - fish 2
site.name = 'NorthTemperateLakes'
locale.name = 'fish2'
data1 = read.csv('./raw_data/ntl7_v7.csv',as.is=T,check.names=F,header=T)
str(data1)
data1$Site = data1$gearid
data1$Number = data1$total_caught
data1$Year = data1$year4
data1$Species = data1$spname
u.years = sort(as.numeric(unique(data1$Year)))
u.species = unique(data1$Species)
method = unique(data1$Site)
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (l in 1:length(method)){
	l.data = data1[which(data1$Site==method[l]),]
	for (s in 1:length(u.species)){
		for (y in 1:length(u.years)){
			y.data = l.data[which(l.data$Year==u.years[y] & l.data$Species==u.species[s]),]
			n.y1 = sum(y.data$effort,na.rm=T)
			out.slopes[cx,1] = site.name #LTER.site
			out.slopes[cx,2] = paste(locale.name,method[l],sep='-') #Locale
			out.slopes[cx,3] = u.species[s] #Species.code
			out.slopes[cx,4] = u.years[y] #Year
			out.slopes[cx,5] = n.y1 #Number of observations
			out.slopes[cx,6] = sum(y.data$Number,na.rm=T) #abundance
			cx = cx + 1
		}
	}
}
write.table(out.slopes,paste0('./summary_tables/vertebrates/',site.name,'_Fish2_SpeciesSlopes.txt'),sep=';',quote=F,row.names=F)


# Sevilleta - birds
site.name = 'Sevilleta'
locale.name = 'birds'
data1 = read.csv('./raw_data/sev017_birdpop_09251997_0.txt',as.is=T,check.names=F,header=T)
str(data1)
data1$Year = data1$year
data1$Species = data1$species_number
u.years = sort(as.numeric(unique(data1$Year)))
u.species = unique(data1$Species)
out.slopes = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999)
cx = 1
for (s in 1:length(u.species)){
	for (y in 1:length(u.years)){
		y.data = data1[which(data1$Year==u.years[y] & data1$Species==u.species[s]),]
		n.y1 = length(unique(y.data$day))
		out.slopes[cx,1] = site.name #LTER.site
		out.slopes[cx,2] = locale.name #Locale
		out.slopes[cx,3] = u.species[s] #Species.code
		out.slopes[cx,4] = u.years[y] #Year
		out.slopes[cx,5] = n.y1 #Number of observations
		out.slopes[cx,6] = nrow(y.data) #abundance
		cx = cx + 1
	}
}
write.table(out.slopes,paste0('./summary_tables/vertebrates/',site.name,'_Bird_SpeciesSlopes.txt'),sep=';',quote=F,row.names=F)



##############################################################Section
# Estimate time trends in bird and fish abundance

data1 = read.csv('./Bird_Abundance_LTER_annotated.csv',as.is=T,check.names=F,header=T)

# Get time trends
LLS = paste(data1$LTER.site,data1$Locale,data1$Species.code,sep='_') #unique LTER by sub-site by species combinations
u.LLS = unique(LLS)
LLS.count = apply(array(u.LLS),1,function(x){d2=data1[which(LLS==x),];length(which(d2$Abundance!=0))})
u.LLS2 = u.LLS[which(LLS.count>0)]
time.trends.Z = data.frame('LTER'=NA,'Site'=NA,'Species'=NA,'MSE'=-999,'b'=-999,'coef_int'=-999,'coef_slope'=-999,'Pr_int'=-999,'Pr_slope'=-999,'logLik'=-999,'length'=-999,'Y1'=NA,'Y2'=NA,'Feeding'=NA)
tx = 1
for (l in 1:length(u.LLS2)){
	print(l)
	lls.data = data1[which(LLS==u.LLS2[l]),]
	if ( (length(which(!is.na(lls.data$Abundance))) > 5) & (length(which(lls.data$Abundance[!is.na(lls.data$Abundance)] > 0)) > 4) ){
		#quality threshold: at least 5 non-zero years
		cmin = min(lls.data$Abundance[which(lls.data$Abundance>0)],na.rm=T)
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
			X[which(!is.na(X))] = lls.data$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
			Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
			t.scale = 1:length(lls.years.stretch)
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
		tx = tx + 1
		# Plot high quality time series
		png(paste0('./plots/time_trends/vertebrates/lineplots/',gsub('\\?','',u.LLS[l]),'.png'))
		plot(t.scale,Z,type='l',main=u.LLS2[l],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2)
		abline(a=add.coef.Z[3],b=add.coef.Z[4],lty=2,col='red',lwd=1.5)
		legend('topright',legend=c(paste0('slope = ',round(add.coef.Z[4],2)),paste0('length = ',length(which(!is.na(Z)))),paste0('autocor = ',round(add.coef.Z[2],2)),paste0('logLik = ',round(add.coef.Z[7],2))),ncol=1)
		dev.off()
	} else {
		#ignore low quality time series
	}
}
write.csv(time.trends.Z,'time_trends_birds.csv',quote=F,row.names=F)



########################################################################################Section
# Fish

# Import species abundance data
files = list.files('./summary_tables/vertebrates',full.names=T)
files = files[grep('Fish',files)]

abundances = read.table(files[1],sep=';',as.is=T,check.names=F,header=T)
for (f in 2:length(files)){
	print(f)
	dat1 = read.table(files[f],sep=';',as.is=T,check.names=F,header=T)
	abundances = data.frame(rbind(abundances,dat1))
}
data1 = abundances

# Get time trends
LLS = paste(data1$LTER.site,data1$Locale,data1$Species.code,sep='_') #unique LTER by sub-site by species combinations
u.LLS = unique(LLS)
LLS.count = apply(array(u.LLS),1,function(x){d2=data1[which(LLS==x),];length(which(d2$Abundance!=0))})
u.LLS2 = u.LLS[which(LLS.count>0)]
time.trends.Z = data.frame('LTER'=NA,'Site'=NA,'Species'=NA,'MSE'=-999,'b'=-999,'coef_int'=-999,'coef_slope'=-999,'Pr_int'=-999,'Pr_slope'=-999,'logLik'=-999,'length'=-999,'Y1'=NA,'Y2'=NA)
tx = 1
for (l in 1:length(u.LLS2)){
	print(l)
	lls.data = data1[which(LLS==u.LLS2[l]),]
	if ( (length(which(!is.na(lls.data$Abundance))) > 5) & (length(which(lls.data$Abundance[!is.na(lls.data$Abundance)] > 0)) > 4) ){
		#quality threshold: at least 5 non-zero years
		cmin = min(lls.data$Abundance[which(lls.data$Abundance>0)],na.rm=T)
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
			X[which(!is.na(X))] = lls.data$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
			Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
			t.scale = 1:length(lls.years.stretch)
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
		tx = tx + 1
		# Plot high quality time series
		png(paste0('./plots/time_trends/vertebrates/lineplots/',gsub('\\?','',u.LLS[l]),'.png'))
		plot(t.scale,Z,type='l',main=u.LLS2[l],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2)
		abline(a=add.coef.Z[3],b=add.coef.Z[4],lty=2,col='red',lwd=1.5)
		legend('topright',legend=c(paste0('slope = ',round(add.coef.Z[4],2)),paste0('length = ',length(which(!is.na(Z)))),paste0('autocor = ',round(add.coef.Z[2],2)),paste0('logLik = ',round(add.coef.Z[7],2))),ncol=1)
		dev.off()
	} else {
		#ignore low quality time series
	}
}
write.csv(time.trends.Z,'time_trends_fish.csv',quote=F,row.names=F)



###############################################################Section
# Visualize patterns of bird and fish abundance time trends

# Birds
trends = read.csv('time_trends_birds.csv',as.is=T,check.names=F,header=T)
trends$LL = as.factor(trends$LTER)

# All bird taxa
png('./plots/time_trends/vertebrates/ALL_Bird_time_trends_boxplot.png',width=480)
par(mfrow=c(1,1),cex.lab=1.5,cex.axis=1.2,lwd=1.5,oma=c(0,0,0,0),mar=c(12,8,1,1))
plot(trends$LL,trends$coef_slope,xaxt='n',ylab='Change in abundance')
text(x=seq(1,length(levels(trends$LL)),1),y=par("usr")[3]-0.2,labels=levels(trends$LL),srt=45,pos=2,xpd=T,cex=1.5)
abline(a=0,b=0,lty=1,col='black',lwd=1)
dev.off()

# Omnivore
trends2 = trends[which(trends$Feeding=='Omnivore'),]
png('./plots/time_trends/vertebrates/OMNIVORE_Bird_time_trends_boxplot.png',width=480)
par(mfrow=c(1,1),cex.lab=1.5,cex.axis=1.2,lwd=1.5,oma=c(0,0,0,0),mar=c(12,8,1,1))
plot(trends2$LL,trends2$coef_slope,xaxt='n',ylab='Change in abundance')
text(x=seq(1,length(levels(trends$LL)),1),y=par("usr")[3]-0.2,labels=levels(trends$LL),srt=45,pos=2,xpd=T,cex=1.5)
abline(a=0,b=0,lty=1,col='black',lwd=1)
dev.off()

# Non-arthropodivores
trends2 = trends[which(trends$Feeding=='Non-arthropodivore'),]
png('./plots/time_trends/vertebrates/NON-ARTHROPODIVORE_Bird_time_trends_boxplot.png',width=480)
par(mfrow=c(1,1),cex.lab=1.5,cex.axis=1.2,lwd=1.5,oma=c(0,0,0,0),mar=c(12,8,1,1))
plot(trends2$LL,trends2$coef_slope,xaxt='n',ylab='Change in abundance')
text(x=seq(1,length(levels(trends$LL)),1),y=par("usr")[3]-0.2,labels=levels(trends$LL),srt=45,pos=2,xpd=T,cex=1.5)
abline(a=0,b=0,lty=1,col='black',lwd=1)
dev.off()

### ms figure

# Birds - arthropodivores
trends = read.csv('time_trends_birds.csv',as.is=T,check.names=F,header=T)
trends$LL = as.factor(trends$LTER)
trends2 = trends[which(trends$Feeding=='Arthropodivore'),]
png('./plots/time_trends/vertebrates/ARTHROPODIVORE_Bird_time_trends_boxplot.png',width=480)
par(mfrow=c(1,1),cex.lab=2,cex.axis=1.5,lwd=2,oma=c(0,0,0,0),mar=c(12,8,1,1))
plot(trends2$LL,trends2$coef_slope,ylim=c(-4,4),xaxt='n',yaxt='n',axes=F,ylab='Change in abundance')
axis(2,at=seq(-4,4,2),labels=seq(-4,4,2),lwd=2)
text(x=seq(1,length(levels(trends$LL)),1),y=par("usr")[3]-0.2,labels=levels(trends$LL),srt=45,pos=2,xpd=T,cex=1.5)
abline(a=0,b=0,lty=1,col='black',lwd=2)
dev.off()

# Fish
trends = read.csv('time_trends_fish.csv',as.is=T,check.names=F,header=T)
#trends$Site[which(trends$LTER=='NorthTemperateLakes')] = 'fish' #merge NTL records into one boxplot
trends$L = as.factor(trends$LTER)
png('./plots/time_trends/vertebrates/ALL_Fish_time_trends_boxplot.png')
par(mfrow=c(1,1),cex.lab=2,cex.axis=1.5,lwd=2,oma=c(0,0,0,0),mar=c(12,8,1,1))
plot(trends$L,trends$coef_slope,ylim=c(-4,4),col='lightblue',xaxt='n',yaxt='n',axes=F,ylab='Change in abundance')
axis(2,at=seq(-4,4,2),labels=seq(-4,4,2),lwd=2)
text(x=seq(1,length(levels(trends$L)),1),y=par("usr")[3]-0.2,labels=levels(trends$L),srt=45,pos=2,xpd=T,cex=1.5)
abline(a=0,b=0,lty=1,col='black',lwd=2)
dev.off()


########################################################################
# Random Forest analysis

library(randomForest)
library(scales)
library(MASS)
library(tree)

# Tutorials:
# https://towardsdatascience.com/understanding-random-forest-58381e0602d2
# https://towardsdatascience.com/random-forest-in-r-f66adf80ec9
# https://www.rdocumentation.org/packages/randomForest/versions/4.6-14/topics/randomForest
# https://www.r-bloggers.com/how-to-implement-random-forests-in-r/
# https://www.blopig.com/blog/2017/04/a-very-basic-introduction-to-random-forests-using-r/
# http://rstudio-pubs-static.s3.amazonaws.com/156481_80ee6ee3a0414fd38f5d3ad33d14c771.html #Start with this one
# https://www.edureka.co/blog/random-forest-classifier/#Practical%20Implementation%20of%20Random%20Forest%20In%20R

# Import abundance trends
trends = read.csv('time_trends_arthropods_relaxed.csv',as.is=T,check.names=F,header=T)

# Append environmental variables
envars = read.csv('LTER_env.csv',as.is=T,check.names=F,header=T)
trends$Temp = unlist(apply(array(trends$LTER),1,function(x){return(envars$Temp[which(gsub(' ','',envars$LTER)==x)])}))
trends$Precip = unlist(apply(array(trends$LTER),1,function(x){return(envars$Precip[which(gsub(' ','',envars$LTER)==x)])}))
trends$HFI = unlist(apply(array(trends$LTER),1,function(x){return(envars$HFI[which(gsub(' ','',envars$LTER)==x)])}))
trends$Duration = trends$Y2 - trends$Y1
trends$LTER.site = paste(trends$LTER,trends$Site,sep='_')

# Quick visual of correlations among variables
png('./plots/trends_v_env_scatterplots.png',width=480*2,height=480*2)
par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(5,5,1,1),cex.lab=1.5,cex.axis=1.2,cex=2)
plot(trends$Temp,trends$coef_slope,xlab='Temperature',ylab='Abundance trend',pch=16,col=alpha('black',0.2))
lm1 = lm(coef_slope ~ Temp,data=trends)
lm1.s = summary(lm1)
abline(lm1,lwd=2)
legend('topright',legend=c(paste0('m = ',round(lm1.s$coefficients[2,1],4)),paste0('R2 = ',round(lm1.s$r.squared,4)),paste0('P = ',round(lm1.s$coefficients[2,4],4))),bty='n')
plot(trends$Precip,trends$coef_slope,xlab='Precipitation',ylab='Abundance trend',pch=16,col=alpha('black',0.2))
lm1 = lm(coef_slope ~ Precip,data=trends)
lm1.s = summary(lm1)
abline(lm1,lwd=2)
legend('topright',legend=c(paste0('m = ',round(lm1.s$coefficients[2,1],4)),paste0('R2 = ',round(lm1.s$r.squared,4)),paste0('P = ',round(lm1.s$coefficients[2,4],4))),bty='n')
plot(trends$HFI,trends$coef_slope,xlab='Human Footprint Index',ylab='Abundance trend',pch=16,col=alpha('black',0.2))
lm1 = lm(coef_slope ~ HFI,data=trends)
lm1.s = summary(lm1)
abline(lm1,lwd=2)
legend('topright',legend=c(paste0('m = ',round(lm1.s$coefficients[2,1],4)),paste0('R2 = ',round(lm1.s$r.squared,4)),paste0('P = ',round(lm1.s$coefficients[2,4],4))),bty='n')
plot(trends$Y1,trends$coef_slope,xlab='Start year',ylab='Abundance trend',pch=16,col=alpha('black',0.2))
lm1 = lm(coef_slope ~ Y1,data=trends)
lm1.s = summary(lm1)
abline(lm1,lwd=2)
legend('topright',legend=c(paste0('m = ',round(lm1.s$coefficients[2,1],4)),paste0('R2 = ',round(lm1.s$r.squared,4)),paste0('P = ',round(lm1.s$coefficients[2,4],4))),bty='n')
dev.off()

plot(trends$Duration,trends$coef_slope,xlab='Duration',ylab='Abundance trend',pch=16,col=alpha('black',0.2))
lm1 = lm(coef_slope ~ Duration,data=trends)
lm1.s = summary(lm1)
abline(lm1,lwd=2)
legend('topright',legend=c(paste0('m = ',round(lm1.s$coefficients[2,1],4)),paste0('R2 = ',round(lm1.s$r.squared,4)),paste0('P = ',round(lm1.s$coefficients[2,4],4))),bty='n')

plot(trends$Duration,trends$Y1)


##########################################################Section
# Random Forests analysis

# Exctract environmental variables
library('raster')
library('rgdal')
library('rgeos')

# FIX LTERs.shp
lters = readOGR(dsn='./shapefiles/LTERs.shp',layer='LTERs',stringsAsFactors=F,verbose=F)
lters@data$LTER[which(lters@data$LTER=='Goergia Coastal')] = 'Georgia Coastal'
writeOGR(lters,dsn='./shapefiles/LTERs.shp',layer='LTERs',driver='ESRI Shapefile',overwrite=T)

usa = readOGR(dsn='./shapefiles/us_alb/us_alb.shp',layer='us_alb',stringsAsFactors=F,verbose=F)
prj1 = proj4string(usa)
lters84 = spTransform(readOGR(dsn='./shapefiles/LTERs.shp',layer='LTERs',stringsAsFactors=F,verbose=F),CRS('+init=epsg:4326'))
lters = spTransform(readOGR(dsn='./shapefiles/LTERs.shp',layer='LTERs',stringsAsFactors=F,verbose=F),CRS(prj1))
lters = lters[which(lters@data$Included=='Y'),]
#Add Suction Trap Sites to LTER shapefile
stn = spTransform(readOGR(dsn='./shapefiles/STN.shp',layer='STN',stringsAsFactors=F,verbose=F),CRS(prj1)) #suction trap network sites
stn@data$Site[grep('Saginaw',stn@data$Site)] = 'Saginaw' #combine two proximal traps
stn@data$Site[grep('Urbana',stn@data$Site)] = 'Urbana-Champaign' #combine two proximal traps
stn@data$Site = paste('Midwest',stn@data$Site,sep='_')
stn@data = data.frame('LTER'=stn@data$Site,'Included'=rep('Y',length(stn)),stringsAsFactors=F)
lters2 = union(lters,stn)
lters84.2 = union(lters84,spTransform(stn,CRS('+init=epsg:4326')))

# Climate maps
temp = projectRaster(raster('./shapefiles/bioclim/wc2.0_bio_10m_01_meanTemp.tif'),crs=CRS(prj1))
precip = projectRaster(raster('./shapefiles/bioclim/wc2.0_bio_10m_12_Precip.tif'),crs=CRS(prj1))
# Human Footprint Index maps
HFI93 = raster('./shapefiles/HumanFootprintIndex/Maps/HFP1993_int.tif')
HFI09 = raster('./shapefiles/HumanFootprintIndex/Maps/HFP2009_int.tif')
prj2 = proj4string(HFI93)
lters_prj2 = spTransform(lters2,CRS(prj2))

# Add environmental variables to LTER shapefile
lters2@data$Temp = extract(temp,lters2,buffer=10000,fun=mean)
lters2@data$Precip = extract(precip,lters2,buffer=10000,fun=mean)
lters.HFI93 = extract(HFI93,lters_prj2,buffer=10000,fun=mean)
lters.HFI09 = extract(HFI09,lters_prj2,buffer=10000,fun=mean)
lters2@data$HFI = apply(cbind(lters.HFI93,lters.HFI09),1,mean)
lters2@data$Lat = lters84.2@coords[which(lters2@data$Included=='Y'),2]
# Write environmenal variables to file
write.csv(lters2@data,'LTER_env.csv',quote=F,row.names=F)

# Subset and classify factors
trends2 = trends[,c(1,7,12,14,15,18:20)]
trends2$Feeding = as.factor(trends2$Feeding)
trends2$Habitat = as.factor(trends2$Habitat)
#trends2$Order = as.factor(trends2$Order) #trends[,17] - too many levels (50 > max allowed = 32)
trends2$LTER = as.factor(trends2$LTER)
#trends2$LTER.site = as.factor(trends2$LTER.site) #trends[,22] - too many levels (104 > max allowed = 32)

###
# Random Forests
###
train = sample(1:nrow(trends2), nrow(trends2)/2)
trends.test=trends2[-train,"coef_slope"]
rf.trends=randomForest(coef_slope~.,data=trends2,subset=train,mtry=5,importance=TRUE)
yhat.rf = predict(rf.trends,newdata=trends2[-train,])
plot(yhat.rf, trends.test) # How well did our predictors predict trends?
abline(0,1)
mean((yhat.rf-trends.test)^2) #including 3 predictors in each tree minimized mean square error
rf.trends
importance(rf.trends)
varImpPlot(rf.trends)

# Bin high and low trends and treat them as categorical
thresh = 2
trends3 = trends2[which(abs(trends2$coef_slope)>thresh),]
slopes3 = trends3$coef_slope
trends3$coef_slope[which(slopes3<(-1*thresh))] = 'decrease'
trends3$coef_slope[which(slopes3>thresh)] = 'increase'
trends3$coef_slope = as.factor(trends3$coef_slope)
train3 = sample(1:nrow(trends3), nrow(trends3)/2)
trends.test=trends3[-train3,"coef_slope"]
# Find best number of predictors
rf.trends=randomForest(coef_slope~.,data=trends3,subset=train3,mtry=5,importance=TRUE)
yhat.rf = predict(rf.trends,newdata=trends3[-train3,])
importance(rf.trends)
varImpPlot(rf.trends)
rf.trends


###################################################################Section
# Order abundance trends boxplots by predictor variables
envars = read.csv('LTER_env.csv',as.is=T,check.names=F,header=T)
trends.r$Temp = unlist(apply(array(trends.r$LTER),1,function(x){return(envars$Temp[which(gsub(' ','',envars$LTER)==x)])}))
trends.r$Precip = unlist(apply(array(trends.r$LTER),1,function(x){return(envars$Precip[which(gsub(' ','',envars$LTER)==x)])}))
trends.r$HFI = unlist(apply(array(trends.r$LTER),1,function(x){return(envars$HFI[which(gsub(' ','',envars$LTER)==x)])}))
trends.r$Duration = trends.r$Y2 - trends.r$Y1

# Human Footprint Index
LL_L = apply(array(levels(trends.r$LL)),1,function(x){strsplit(x,'_')[[1]][1]})
L_HFI = apply(array(LL_L),1,function(x){envars$HFI[which(gsub(' ','',envars$LTER)==x)]})
trends.r$LL_sort = factor(as.character(trends.r$LL),levels=levels(trends.r$LL)[order(L_HFI)])
cols_sort = cols[order(L_HFI)]
png('./plots/time_trends/abundance/ALL_time_trends_boxplot_relaxed_noannotations_HFI.png',width=500*2,height=500)
par(mfrow=c(1,1),cex.lab=2,cex.axis=1.5,lwd=2,oma=c(0,0,0,0),mar=c(15,8,1,8))
plt = plot(trends.r$LL_sort,trends.r$coef_slope,col=cols_sort,ylim=c(-8,8),xaxt='n',yaxt='n',xlab=NULL,ylab='Change in abundance',axes=F)
axis(2,at=seq(-8,8,2),labels=seq(-8,8,2),lwd=2)
text(x=seq(1,length(levels(trends.r$LL_sort)),1),y=par("usr")[3]-0.2,labels=levels(trends.r$LL_sort),srt=45,pos=2,xpd=T,cex=1.2)
abline(a=0,b=0,lty=1,col='black',lwd=2)
# Add HFI line to boxplot
par(new=T)
plot(1:length(L_HFI),L_HFI,ylim=c(0,40),type='n',axes=F,xlab='',ylab='')
axis(4,at=seq(0,40,10),labels=seq(0,40,10),lwd=2,cex.axis=1.5)
mtext('Human Footprint Index',cex=2,side=4,line=3)
lines(1:length(L_HFI),sort(L_HFI),lwd=3,col='gray50')
dev.off()

# Start year
L_Y1 = apply(array(levels(trends.r$LL)),1,function(x){min(trends.r$Y1[which(as.character(trends.r$LL)==x)])})
trends.r$LL_sort = factor(as.character(trends.r$LL),levels=levels(trends.r$LL)[order(L_Y1)])
cols_sort = cols[order(L_Y1)]
png('./plots/time_trends/abundance/ALL_time_trends_boxplot_relaxed_noannotations_Y1.png',width=500*2,height=500)
par(mfrow=c(1,1),cex.lab=2,cex.axis=1.5,lwd=2,oma=c(0,0,0,0),mar=c(15,8,1,8))
plot(trends.r$LL_sort,trends.r$coef_slope,col=cols_sort,ylim=c(-8,8),xaxt='n',yaxt='n',xlab=NULL,ylab='Change in abundance',axes=F)
axis(2,at=seq(-8,8,2),labels=seq(-8,8,2),lwd=2)
text(x=seq(1,length(levels(trends.r$LL_sort)),1),y=par("usr")[3]-0.2,labels=levels(trends.r$LL_sort),srt=45,pos=2,xpd=T,cex=1.2)
abline(a=0,b=0,lty=1,col='black',lwd=2)
# Add env line to boxplot
par(new=T)
plot(1:length(L_Y1),L_Y1,ylim=c(1970,2020),type='n',axes=F,xlab='',ylab='')
axis(4,at=seq(1970,2020,10),labels=seq(1970,2020,10),lwd=2,cex.axis=1.5)
mtext('Start year',cex=2,side=4,line=3)
lines(1:length(L_Y1),sort(L_Y1),lwd=3,col='gray50')
dev.off()

L_temp = apply(array(LL_L),1,function(x){envars$Temp[which(gsub(' ','',envars$LTER)==x)]})
trends.r$LL_sort = factor(as.character(trends.r$LL),levels=levels(trends.r$LL)[order(L_temp)])
cols_sort = cols[order(L_temp)]
png('./plots/time_trends/abundance/ALL_time_trends_boxplot_relaxed_noannotations_Temp.png',width=500*2,height=500)
par(mfrow=c(1,1),cex.lab=2,cex.axis=1.5,lwd=2,oma=c(0,0,0,0),mar=c(15,8,1,8))
plot(trends.r$LL_sort,trends.r$coef_slope,col=cols_sort,ylim=c(-8,8),xaxt='n',yaxt='n',xlab=NULL,ylab='Change in abundance',axes=F)
axis(2,at=seq(-8,8,2),labels=seq(-8,8,2),lwd=2)
text(x=seq(1,length(levels(trends.r$LL_sort)),1),y=par("usr")[3]-0.2,labels=levels(trends.r$LL_sort),srt=45,pos=2,xpd=T,cex=1.2)
abline(a=0,b=0,lty=1,col='black',lwd=2)
# Add env line to boxplot
par(new=T)
plot(1:length(L_temp),L_temp,ylim=c(-8,24),type='n',axes=F,xlab='',ylab='')
axis(4,at=seq(-8,24,4),labels=seq(-8,24,4),lwd=2,cex.axis=1.5)
mtext('Mean annual temperature',cex=2,side=4,line=3)
lines(1:length(L_temp),sort(L_temp),lwd=3,col='gray50')
dev.off()

L_precip = apply(array(LL_L),1,function(x){envars$Precip[which(gsub(' ','',envars$LTER)==x)]})
trends.r$LL_sort = factor(as.character(trends.r$LL),levels=levels(trends.r$LL)[order(L_precip)])
cols_sort = cols[order(L_precip)]
png('./plots/time_trends/abundance/ALL_time_trends_boxplot_relaxed_noannotations_Precip.png',width=500*2,height=500)
par(mfrow=c(1,1),cex.lab=2,cex.axis=1.5,lwd=2,oma=c(0,0,0,0),mar=c(15,8,1,8))
plot(trends.r$LL_sort,trends.r$coef_slope,col=cols_sort,ylim=c(-8,8),xaxt='n',yaxt='n',xlab=NULL,ylab='Change in abundance',axes=F)
axis(2,at=seq(-8,8,2),labels=seq(-8,8,2),lwd=2)
text(x=seq(1,length(levels(trends.r$LL_sort)),1),y=par("usr")[3]-0.2,labels=levels(trends.r$LL_sort),srt=45,pos=2,xpd=T,cex=1.2)
abline(a=0,b=0,lty=1,col='black',lwd=2)
# Add env line to boxplot
par(new=T)
plot(1:length(L_precip),L_precip,ylim=c(200,1800),type='n',axes=F,xlab='',ylab='')
axis(4,at=seq(200,1800,200),labels=seq(200,1800,200),lwd=2,cex.axis=1.5)
mtext('Annual precipitation',cex=2,side=4,line=3)
lines(1:length(L_precip),sort(L_precip),lwd=3,col='gray50')
dev.off()


### END pipeline