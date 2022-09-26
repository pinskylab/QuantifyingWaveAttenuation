setwd("/Users/mpinsky/Documents/Stanford/NatCap/Coastal Protection Model")
source('Analysis/log10axis.R')
source('Analysis/p.values.lmer.R')

quant = read.csv("data/CoastProt Quant 110703.csv", header=TRUE)
sppsum = read.csv("output/SppSummary_2011-07-04.csv")


##################################
#### Data Prep before MatLab run #
#### (skip if already done)      #
##################################

# Create Hab2 habitat categories: combine bare, mud, and sand, rename artifical kelp and marsh
quant$Hab2 = as.character(quant$Habitat)
quant$Hab2[quant$Hab2 == "Sand"] = "Bare"
quant$Hab2[quant$Hab2 == "Mud"] = "Bare"
quant$Hab2[quant$Hab2 == "Artificial Kelp"] = "Kelp (Art.)"
quant$Hab2[quant$Hab2 == "Artificial Marsh"] = "Marsh (Art.)"
quant$Hab2[quant$Hab2 == "Artificial Seagrass"] = "Seagrass (Art.)"
quant$Hab2[quant$Hab2 == "Reef"] = "Coral Reef"
quant$Hab2 = as.factor(quant$Hab2)
levels(quant$Hab2)

# Create simpler habitat categories (don't differentiate artificial mimics): Hab3
quant$Hab3 = as.character(quant$Hab2)
quant$Hab3[quant$Hab3=="Kelp (Art.)"] = "Kelp"
quant$Hab3[quant$Hab3=="Marsh (Art.)"] = "Marsh"
quant$Hab3[quant$Hab3=="Seagrass (Art.)"] = "Seagrass"

# Add a simple lab/field column
quant$LF = NA
quant$LF[grep("Field", quant$Lab.field)] = "field"
quant$LF[grep("Lab", quant$Lab.field)] = "lab"
quant$Field = 0
quant$Field[quant$LF =='field'] = 1

# calc H2 (attenuated wave height) where we have atten% but not height data
rows = which(is.na(quant$H2) & !is.na(quant$Attenuation.....mean) & !is.na(quant$H1))
quant$H2[rows] = quant$H1[rows] - quant$H1[rows]*quant$Attenuation.....mean[rows]/100

# Rename T (period), h (mean depth), D (density), A (x-sectional area), hv (veg height) 
quant$x1 = quant$Distance.from.H1.to.seaward.bed.edge..m.
	quant$x1[is.na(quant$x1)] = 0 ## Would be better to verify each study
quant$T = quant$T.mean..s.
quant$h = quant$h.mean..m.
quant$D = as.numeric(as.character(quant$Density..shoots.m2..mean))
quant$A = quant$Stem.diameter..mm./1000 # convert mm to m
quant$hv = quant$Vegetation.height..m.

# calc %Attenuation for all studies: fraction lower that Hx is than H1, in terms of H1. >0 shows attenuation.
quant$atten = NA
for(i in 1:nrow(quant)){
	Hx = quant[i, c('H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8')]
	indx = which(!is.na(Hx)) # index of last wave height
	if(length(indx)>0 & !is.na(quant$H1[i])){
		indx = max(indx)+1
		Hx = quant[i, paste('H', indx, sep='')]
		quant$atten[i] = (quant$H1[i]-Hx)/quant$H1[i]
	}
}

# Fill in h1 and h2 where we have only mean depth
# assume flat, non-sloping bed (may not be true! need to verify in studies)
rows = !is.na(quant$h) & is.na(quant$h1) & is.na(quant$h2)
quant$h1[rows] = quant$h[rows]
quant$h2[rows] = quant$h[rows]

# Fill in h (average depth) where we have individual depths
rows = is.na(quant$h)
quant$h[rows] = rowMeans(quant[rows, c('h1', 'h2', 'h3', 'h4', 'h5','h6', 'h7','h8')], na.rm=T)

# Fill in h_hv (relative depth)
quant$h_hv = quant$h/quant$hv


# Calculate slope
quant$slope = NA
for(i in 1:nrow(quant)){
	x = as.numeric(c(0, quant[i,c('x2', 'x3', 'x4', 'x5', 'x6', 'x7', 'x8', 'x9')]))
	hv = as.numeric(quant[i,c('h1', 'h2', 'h3', 'h4', 'h5', 'h6', 'h7', 'h8', 'h9')])
	j = !is.na(x) & !is.na(hv)
	if(sum(j) >= 2){
		mod = lm(hv ~ x)
		quant$slope[i] = round(atan(-coef(mod)[2])*360/2/pi,2) # convert slope to degrees
	}
}

# Calculate maximum distance
quant$xmax = apply(quant[, c('x1', 'x2', 'x3', 'x4', 'x5','x6', 'x7','x8')], MARGIN=1, FUN=max, na.rm=T)

# Convert wave heights to Hrms
unique(quant$Wave.Height.Type)

	# Wavetype = 1 for Hrms and 0 for monochromatic
	quant$WaveType = NA
	quant$WaveTypeVerbose = quant$Wave.Height.Type

	i = which(quant$Wave.Height.Type=="Hrms" | quant$Wave.Height.Type=="Hrmsfreq" | quant$Wave.Height.Type=="Hrmstime")
	quant$WaveType[i] = 1
	quant$WaveTypeVerbose[i] = "Hrms"

	i = which(quant$Wave.Height.Type=="H1/3" | quant$Wave.Height.Type=="Hs")
	quant$H1[i] = quant$H1[i]/sqrt(log(1/(1/3)))
	quant$H2[i] = quant$H2[i]/sqrt(log(1/(1/3)))
	quant$H3[i] = quant$H3[i]/sqrt(log(1/(1/3)))
	quant$H4[i] = quant$H4[i]/sqrt(log(1/(1/3)))
	quant$H5[i] = quant$H5[i]/sqrt(log(1/(1/3)))
	quant$H6[i] = quant$H6[i]/sqrt(log(1/(1/3)))
	quant$H7[i] = quant$H7[i]/sqrt(log(1/(1/3)))
	quant$H8[i] = quant$H8[i]/sqrt(log(1/(1/3)))
	quant$WaveType[i] = 1
	quant$WaveTypeVerbose[i] = "Hrms"
	
	i = which(quant$Wave.Height.Type=="Hm0")
	quant$H1[i] = quant$H1[i]/sqrt(2)
	quant$H2[i] = quant$H2[i]/sqrt(2)
	quant$H3[i] = quant$H3[i]/sqrt(2)
	quant$H4[i] = quant$H4[i]/sqrt(2)
	quant$H5[i] = quant$H5[i]/sqrt(2)
	quant$H6[i] = quant$H6[i]/sqrt(2)
	quant$H7[i] = quant$H7[i]/sqrt(2)
	quant$H8[i] = quant$H8[i]/sqrt(2)
	quant$WaveType[i] = 1
	quant$WaveTypeVerbose[i] = "Hrms"
	
	i = which(quant$Wave.Height.Type=="Hmax" & quant$Study=="Allen et al. 2008") # assume 5 waves in a boat wake
	quant$H1[i] = quant$H1[i]/sqrt(log(1/(1/5)))
	quant$H2[i] = quant$H2[i]/sqrt(log(1/(1/5)))
	quant$H3[i] = quant$H3[i]/sqrt(log(1/(1/5)))
	quant$H4[i] = quant$H4[i]/sqrt(log(1/(1/5)))
	quant$H5[i] = quant$H5[i]/sqrt(log(1/(1/5)))
	quant$H6[i] = quant$H6[i]/sqrt(log(1/(1/5)))
	quant$H7[i] = quant$H7[i]/sqrt(log(1/(1/5)))
	quant$H8[i] = quant$H8[i]/sqrt(log(1/(1/5)))
	quant$WaveType[i] = 1
	quant$WaveTypeVerbose[i] = "Hrms"
	
	i = which(quant$Wave.Height.Type=="Hmax" & quant$Study=="Morgan et al. 2009 Estuaries and Coasts") # assume 10 waves in a boat wake from Morgan study
	quant$H1[i] = quant$H1[i]/sqrt(log(1/(1/10)))
	quant$H2[i] = quant$H2[i]/sqrt(log(1/(1/10)))
	quant$H3[i] = quant$H3[i]/sqrt(log(1/(1/10)))
	quant$WaveType[i] = 1
	quant$WaveTypeVerbose[i] = "Hrms"
	
	i = which(quant$Wave.Height.Type=="mean of top 5 waves" & quant$Study=="Bonham 1983") # Fig. 3 suggests ~ 11 waves observed
	quant$H1[i] = quant$H1[i]/sqrt(log(1/(5/11)))
	quant$H2[i] = quant$H2[i]/sqrt(log(1/(5/11)))
	quant$H3[i] = quant$H3[i]/sqrt(log(1/(5/11)))
	quant$H4[i] = quant$H4[i]/sqrt(log(1/(5/11)))
	quant$H5[i] = quant$H5[i]/sqrt(log(1/(5/11)))
	quant$H6[i] = quant$H6[i]/sqrt(log(1/(5/11)))
	quant$H7[i] = quant$H7[i]/sqrt(log(1/(5/11)))
	quant$H8[i] = quant$H8[i]/sqrt(log(1/(5/11)))
	quant$WaveType[i] = 1
	quant$WaveTypeVerbose[i] = "Hrms"

	i = which(quant$Wave.Height.Type=="Mean of three highest waves" & quant$Study == "Knutson et al 1982") # Assume a boat wake has about 10 waves in it
	quant$H1[i] = quant$H1[i]/sqrt(log(1/(3/10)))
	quant$H2[i] = quant$H2[i]/sqrt(log(1/(3/10)))
	quant$H3[i] = quant$H3[i]/sqrt(log(1/(3/10)))
	quant$H4[i] = quant$H4[i]/sqrt(log(1/(3/10)))
	quant$H5[i] = quant$H5[i]/sqrt(log(1/(3/10)))
	quant$H6[i] = quant$H6[i]/sqrt(log(1/(3/10)))
	quant$H7[i] = quant$H7[i]/sqrt(log(1/(3/10)))
	quant$H8[i] = quant$H8[i]/sqrt(log(1/(3/10)))
	quant$WaveType[i] = 1
	quant$WaveTypeVerbose[i] = "Hrms"
		
	i = which(quant$Wave.Height.Type=="monochromatic")
	quant$WaveType[i] = 0
		

# Fill in veg values where we can from SppData
# Ok that we're missing Rhizophora spp and Aegiceras spp: these are from Massel 1999 where we don't have H data
inds = which(is.na(quant$A))
for(i in inds){
	k = numeric(0)
	for(j in 1:nrow(sppsum)){ # find all spp listed in Species[i]
		if(length(grep(sppsum$Spp[j], quant$Species[i]))>0) k = c(k,j)
	}
	if(length(k)>0) quant$A[i] = mean(sppsum$A[k], na.rm=TRUE)
	else print(paste("No A match for", quant$Species[i]))
}

inds = which(is.na(quant$D))
for(i in inds){
	k = numeric(0)
	for(j in 1:nrow(sppsum)){ # find all spp listed in Species[i]
		if(length(grep(sppsum$Spp[j], quant$Species[i]))>0) k = c(k,j)
	}
	if(length(k)>0) quant$D[i] = mean(sppsum$D[k], na.rm=TRUE)
	else print(paste("No D match for", quant$Species[i]))
}

inds = which(is.na(quant$hv))
for(i in inds){
	k = numeric(0)
	for(j in 1:nrow(sppsum)){ # find all spp listed in Species[i]
		if(length(grep(sppsum$Spp[j], quant$Species[i]))>0) k = c(k,j)
	}
	if(length(k)>0) quant$hv[i] = mean(sppsum$hv[k], na.rm=TRUE)
	else print(paste("No hv match for", quant$Species[i]))
}


# calc exponential decay coeff from wave measurements (H1 and Hx)
# Don't eliminate any wave measurements outside the veg
# Use linear fit in log space because errors are symmetrically distributed (it's a log ratio)
	quant$k = NA
	for(i in 1:nrow(quant)){
		Hx = as.numeric(quant[i, c('H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8')])
		x = as.numeric(c(0, quant[i, c('x2', 'x3', 'x4', 'x5', 'x6', 'x7', 'x8')]))
		ind = which(!is.na(Hx))
		Hx = Hx[ind]/Hx[ind][1]
		x = x[ind]
		if(length(ind)>1){
			fit = lm(log(Hx) ~ x - 1) # linear fit through the origin
			quant$k[i] = -coef(fit)
		} else {
			print(paste('Did not calc for i=', i)) # if not enough Hx data
		}
	}

# Label unique species/study and habitat/study combinations
quant$spstudy = paste(quant$Species, quant$Study, sep="") 
quant$LFstudy = paste(quant$LF, quant$Hab3, quant$Study, sep='')
quant$LFHab3 = paste(quant$Hab3, quant$LF, sep='')
quant$hab2study = paste(quant$Hab2, quant$Study, sep='')
quant$habstudy = paste(quant$Hab3, quant$Study, sep='')

# Sort by habitat and lab/field
i = order(as.character(quant$Hab2), as.character(quant$LF))
quant = quant[i,]


## Write out
write.csv(quant, paste('output/quant_', Sys.Date(), '.csv', sep=''), row.names=FALSE)


###################
# Sample size
###################
setwd("/Users/mpinsky/Documents/Stanford/NatCap/Coastal Protection Model")
quant = read.csv('output/quant_2011-07-04.csv')

	table(quant$Habitat)
	table(quant$Habitat[quant$Study !="Coops et al. 1996"])
	
	# Count each species within each study as a unique data point
	length(unique(quant$spstudy)) # 60 total studies
		table(quant$Hab2[!duplicated(quant$spstudy)]) 
		
	length(unique(quant$spstudy[quant$Hab2 != "Bare"])) # 54 studies (excluding bare)
	length(unique(quant$spstudy[quant$Hab2 != "Bare" & quant$Hab2 != "Coral Reef"])) # 48 studies (excluding bare and coral reef)
	length((quant$spstudy[quant$Hab2 != "Bare" & quant$Hab2 != "Coral Reef"])) # 614 data points (excluding bare and coral reef)
	
	length(i<-unique(quant$spstudy[quant$Hab2 != "Bare" & quant$atten>0 & !is.na(quant$atten)])) # 49 studies show +atten over veg, plus XX studies that had attenuation (I checked) but don't have wave height measurements = XX studies
	setdiff(unique(quant$spstudy[quant$Hab2 !="Bare"]), i) # only Elwany doesn't show attenuation
	
	j = c('Mangrove', 'Kelp', 'Kelp (Art.)', 'Marsh', 'Marsh (Art.)', 'Seagrass', 'Seagrass (Art.)')
	j = c('Kelp', 'Kelp (Art.)', 'Marsh', 'Marsh (Art.)', 'Seagrass', 'Seagrass (Art.)')
	j = c('Mangrove', 'Marsh', 'Marsh (Art.)', 'Seagrass', 'Seagrass (Art.)')
	length(i<-unique(quant$spstudy[quant$Hab2 %in% j & quant$atten>0 & !is.na(quant$atten)])) 
	setdiff(unique(quant$spstudy[quant$Hab2 %in% j]), i) # only Elwany doesn't show attenuation
	
	table(quant$Hab2[!duplicated(quant$spstudy)])

	# Only lab studies
	table(quant$Hab2[!duplicated(quant$spstudy) & quant$Lab.field=="Lab"])
	
	spstudy = paste(quant$Species[i], quant$Study[i], sep="") # Studies without missing data
	table(quant$Hab2[i][!duplicated(spstudy)])

#######################################################
## Find studies missing full data for Cd calculation ##
#######################################################

i = (!is.na(quant$H1) & !is.na(quant$H2) & !is.na(quant$x2) & !is.na(quant$T) & !is.na(quant$h1) & !is.na(quant$h2) & !is.na(quant$A) & !is.na(quant$D) & !is.na(quant$hv))
miss = quant[!i, c("Study", "Lab.field", "Species", "Hab2", "H1", "H2", "x2", "T", "h1", "h2", "A", "D", "hv", 'spstudy')] 

write.csv(miss, paste("output/Missing_data", Sys.Date(), ".csv", sep="")) # write out


###########################################
## Output All Habitats for Greg's Fit_Cd ##
###########################################

outname = paste("output/QuantforGreg_", Sys.Date(), ".txt", sep="")
# Exclude Bonham because in a river channel with a steep bank: Greg says model won't work
# Also exclude Coops et al. 1996 because it suffers from resonant waves: model won't work
i =which(!is.na(quant$WaveType) & !is.na(quant$H1) & !is.na(quant$H2) & !is.na(quant$x2) & !is.na(quant$T) & !is.na(quant$h1) & !is.na(quant$h2) & !is.na(quant$A) & !is.na(quant$D) & !is.na(quant$hv) & quant$Study != "Bonham 1983" & quant$Study != "Coops et al. 1996")

#i2 =which(!is.na(quant$WaveType) & !is.na(quant$H1) & !is.na(quant$H2) & !is.na(quant$x2) & !is.na(quant$T) & !is.na(quant$h1) & !is.na(quant$h2) & quant$Hab2=="Bare" & quant$Study != "Bonham 1983" & quant$Study != "Coops et al. 1996") # don't include the bare studies: Greg's code can't do anything with them (but perhaps could estimate B?)
#i = c(i1, i2)

out = quant[i,]
out$kk = 1:nrow(out)

cols = c("kk", "PairID", "WaveType", "T", "D", "A", "hv", "H1", "x1", "h1", "H2", "x2", "h2", "H3", "x3", "h3", "H4", "x4", "h4", "H5", "x5", "h5", "H6", "x6", "h6", "H7", "x7", "h7", "H8", "x8", "h8", 'H9', 'x9', 'h9', 'Field', 'id', 'StudyID')
write.table(paste("%", paste(cols, collapse="\t"), sep=""), file=outname, row.names=F, quote=F, col.names=F)
write.table(out[,cols], file=outname, sep="\t", append=T, row.names=F, col.names=F, quote=F, na='NaN')

cols = c("kk", "id", "StudyID", "PairID", "Study", "LF", "Location", "Hab2", "Species", "Description", "WaveType", "WaveTypeVerbose", "Friction.Coef.mean", "Type.of.Friction.Coefficient", "T", "D", "A", "hv", "H1", "x1", "h1", "H2", "x2", "h2", "H3", "x3", "h3", "H4", "x4", "h4", "H5", "x5", "h5", "H6", "x6", "h6", "H7", "x7", "h7", "H8", "x8", "h8", 'H9', 'x9', 'h9', "Notes")
outname = paste("output/QuantforGregVerbose_", Sys.Date(), ".csv", sep="")
# have to use lowercase id for it to open in excel
write.table(paste(cols, collapse=","), file=outname, row.names=F, col.names=F, quote=F)
write.table(out[,cols], file=outname, sep=",", row.names=F, append=T, col.names=F)


########################################################################
## Skip if already done:                                              ##
##                       											  ##
## After Matlab run:                                                  ##
## Merge quant and fitcd, aggregate by hab3-study                     ##
## and trim to kelp, marsh, mangrove, seagrass                        ##
########################################################################

setwd("/Users/mpinsky/Documents/Stanford/NatCap/Coastal Protection Model")
quant = read.csv('output/quant_2011-07-04.csv')
	dim(quant)
	
fitcdfile = "analysis/Fit_Cd 110622/VegData2Out_25-Mar-2012.csv" # all observations, inside or outside the veg bed
fitcd = read.csv(fitcdfile, header=T) # output from MATLAB
	dim(fitcd)
KC2 = read.csv('analysis/Fit_Cd 110622/VegData2_Re_04-Mar-2013.csv', header=T)
	dim(KC2)
	
	# turn -999 or 999 or infinite to NA
	fitcd$Cdval2[fitcd$Cdval2 == -999] = NA
	fitcd$Difval2[fitcd$Difval2 == 999] = NA
	fitcd$slo[!is.finite(fitcd$slo)] = NA
	fitcd$ReMean[fitcd$ReMean==0] = NA
	fitcd$KCMean[!is.finite(fitcd$KCMean)] = NA
	KC2$KC2Mean[!is.finite(KC2$KC2Mean)] = NA
	KC2$KC2Sd[!is.finite(KC2$KC2Sd)] = NA

# Merge with fitcd
quantcd = merge(quant, fitcd[,c('id', 'X.kk', 'nHfit', 'Cdval2', 'Bval2', 'Difval2', 'R2val2', 'slo', 'flatbed', 'ReMean', 'ReSd', 'KCMean', 'KCSd')], by.x = 'id', by.y = 'id', all.x=T)
	dim(quantcd)
quantcd = merge(quantcd, KC2[,c('id', 'KC2Mean', 'KC2Sd')], by.x = 'id', by.y = 'id', all.x=T)
	dim(quantcd)

# Trim to vegetation and exclude problem studies
# Bonham 1983: slope is too steep in the river channel
# Coops et al. 1996: suffers from resonant standing waves in a wave tank
# Dubi & Torum 1995: same data as Mendez & Losada 2004, but harder to digitize
quantcdtrim = quantcd[quantcd$Hab3 %in% c('Marsh', 'Seagrass', 'Kelp', 'Mangrove'),]
quantcdtrim = quantcdtrim[!(quantcdtrim$Study %in% c('Bonham 1983', 'Coops et al. 1996', 'Dubi & Torum 1995')),]
quantcdtrim = droplevels(quantcdtrim)
quantcdtrim$Cdnozero = pmax(quantcdtrim$Cd, 0) # field with Cd no lower than zero
	dim(quantcdtrim)

# Calculate R2
	lreal = function(x) return(sum(!is.na(x))) # return the number of non-NA values in a vector
	sstot = function(x){ x = x[!is.na(x)]; mu = mean(x); return(sum((x-mu)^2))} # return sum square variation, after removing NAs
	L = apply(quantcdtrim[,c('H1','H2','H3','H4','H5','H6','H7','H8','H9')], MARGIN=1, FUN=lreal) # number of wave measurements
	SSerror = (quantcdtrim$Difval2^2)*L # back-calculate sum square error from RMSE
	SStotal = apply(quantcdtrim[,c('H1','H2','H3','H4','H5','H6','H7','H8','H9')], MARGIN=1, FUN=sstot) # total sum of squares
		i = which(L == quantcdtrim$nHfit+1) # find studies for which we removed the first measurement (outside bed or reflection problem)
		L[i] = L[i] - 1
		SSerror[i] = (quantcdtrim$Difval2[i]^2)*L[i]
		SStotal[i] = apply(quantcdtrim[i,c('H2','H3','H4','H5','H6','H7','H8','H9')], MARGIN=1, FUN=sstot)
	if(length(which(L>quantcdtrim$nHfit))>0) warning('Some studies had > 1 measurement removed')
	quantcdtrim$R2 = 1-SSerror/SStotal
		i = which(!is.na(quantcdtrim$R2val2)) # compare R2 to R2val2 from Matlab
		quantcdtrim[i, c('Study', 'Habitat', 'Cdval2', 'R2val2', 'R2')] # perfect match
		
	# Fix -Inf (no variation in observed waves) to NA
	quantcdtrim$R2[is.infinite(quantcdtrim$R2)] = NA

	range(quantcdtrim$R2, na.rm=T)
		i = which(quantcdtrim$R2<0.5)
		quantcdtrim[i, c('Study', 'Habitat', 'Lab.field', 'Cdval2', 'Difval2', 'R2val2', 'R2')]
			
# Add rest of nHfit (number of wave measurements per xsect)
	i = is.na(quantcdtrim$nHfit)
	quantcdtrim$nHfit[i] = L[i]

	cbind(quantcdtrim$nHfit[!i], L[!i]) # compare L with those cases we didn't replace: the few that Matlab calculated. Should now be equal, after fix above.
			
# Sort by habitat and lab/field
i = order(as.character(quantcdtrim$Hab2), as.character(quantcdtrim$LF))
quantcdtrim = quantcdtrim[i,]
	dim(quantcdtrim)

# Write out full data file
	write.csv(quantcdtrim, paste("output/quantcd_allH", Sys.Date(), ".csv", sep=''), row.names=F)


## Aggregate by LF-hab3-study

	quantcdlf = aggregate(list(nHfit = quantcdtrim$nHfit, H1 = quantcdtrim$H1, xmax=quantcdtrim$xmax, slo = quantcdtrim$slo, h = quantcdtrim$h, h_hv = quantcdtrim$h_hv, T= quantcdtrim$T, hv = quantcdtrim$hv, D = quantcdtrim$D, A = quantcdtrim$A, KC = quantcdtrim$KCMean, Re = quantcdtrim$ReMean, atten=quantcdtrim$atten, k=quantcdtrim$k, Cd = quantcdtrim$Cdval2, RMSE = quantcdtrim$Difval2, R2 = quantcdtrim$R2), by=list(LFstudy = quantcdtrim$LFstudy), FUN=mean, na.rm=T)
		dim(quantcdlf)
	
	quantlfsd = aggregate(list(nHfitsd = quantcdtrim$nHfit, H1sd = quantcdtrim$H1, xmaxsd=quantcdtrim$xmax, slosd = quantcdtrim$slo, hsd = quantcdtrim$h, h_hvsd = quantcdtrim$h_hv, Tsd= quantcdtrim$T, hvsd = quantcdtrim$hv, Dsd = quantcdtrim$D, Asd = quantcdtrim$A, KCsd = quantcdtrim$KCMean, Resd = quantcdtrim$ReMean, attensd = quantcdtrim$atten, ksd=quantcdtrim$k, Cdsd = quantcdtrim$Cdval2, RMSEsd = quantcdtrim$Difval2, R2sd = quantcdtrim$R2), by=list(LFstudy = quantcdtrim$LFstudy), FUN=sd, na.rm=T) # standard deviations
		dim(quantlfsd)
	
	quantlfn = aggregate(list(n = rep(1, nrow(quantcdtrim))), by=list(LFstudy = quantcdtrim$LFstudy), FUN=sum, na.rm=T)
		dim(quantlfn)
	
	# treat all Cd < 0 as 0
	quantlfnozero = aggregate(list(Cdnozero = pmax(quantcdtrim$Cdval2, 0)), by=list(LFstudy = quantcdtrim$LFstudy), FUN=mean, na.rm=T)
		dim(quantlfnozero)
	
	# logs ignore all values <= 0, but this is reasonable (waves shouldn't be able to gain energy)
	quantlfcdlog = aggregate(list(KClog = log10(quantcdtrim$KCMean), Relog = log10(quantcdtrim$ReMean), klog=log10(quantcdtrim$k), Cdlog = log10(quantcdtrim$Cdval2)), by=list(LFstudy = quantcdtrim$LFstudy), FUN=mean, na.rm=T) # means in log10 space after adjusting for negative values
		dim(quantlfcdlog)
	
	quantlfsdlog = aggregate(list(KCsdlog = log10(quantcdtrim$KCMean), Resdlog = log10(quantcdtrim$ReMean), ksdlog=log10(quantcdtrim$k), Cdsdlog = log10(quantcdtrim$Cdval2)), by=list(LFstudy = quantcdtrim$LFstudy), FUN=sd, na.rm=T) # standard deviations in log space after adjusting for negative values
		dim(quantlfsdlog)

	quantcdlf = merge(quantcdlf, quantlfn)
		dim(quantcdlf)
	quantcdlf = merge(quantcdlf, quantlfsd)
		dim(quantcdlf)
	quantcdlf = merge(quantcdlf, quantlfnozero)
		dim(quantcdlf)
	quantcdlf = merge(quantcdlf, quantlfcdlog)
		dim(quantcdlf)
	quantcdlf = merge(quantcdlf, quantlfsdlog)
		dim(quantcdlf)
	quantcdlf = merge(quantcdlf, quantcdtrim[!duplicated(quantcdtrim$LFstudy),c('LFstudy', 'LF', 'LFHab3', 'Hab3', 'Study')])
		dim(quantcdlf)

	# Write out aggregated data file
	write.csv(quantcdlf, paste("output/quantcdlf_allH", Sys.Date(), ".csv", sep=''), row.names=F)
	

## Aggregate by hab2-study

	quantcdhab2 = aggregate(list(H1 = quantcdtrim$H1, xmax=quantcdtrim$xmax, h = quantcdtrim$h, h_hv = quantcdtrim$h_hv, T= quantcdtrim$T, hv = quantcdtrim$hv, D = quantcdtrim$D, KC = quantcdtrim$KCMean, Re = quantcdtrim$ReMean, atten=quantcdtrim$atten, k=quantcdtrim$k, Cd = quantcdtrim$Cdval2), by=list(hab2study = quantcdtrim$hab2study), FUN=mean, na.rm=T)
		dim(quantcdhab2)
	
	quanthab2sd = aggregate(list(H1sd = quantcdtrim$H1, xmaxsd=quantcdtrim$xmax, hsd = quantcdtrim$h, h_hvsd = quantcdtrim$h_hv, Tsd= quantcdtrim$T, hvsd = quantcdtrim$hv, Dsd = quantcdtrim$D, KCsd = quantcdtrim$KCMean, Resd = quantcdtrim$ReMean, attensd = quantcdtrim$atten, ksd=quantcdtrim$k, Cdsd = quantcdtrim$Cdval2), by=list(hab2study = quantcdtrim$hab2study), FUN=sd, na.rm=T) # standard deviations
		dim(quanthab2sd)
	
	quanthab2n = aggregate(list(n = rep(1, nrow(quantcdtrim))), by=list(hab2study = quantcdtrim$hab2study), FUN=sum, na.rm=T)
		dim(quanthab2n)
	
	# treat all Cd < 0 as 0
	quanthab2nozero = aggregate(list(Cdnozero = pmax(quantcdtrim$Cdval2, 0)), by=list(hab2study = quantcdtrim$hab2study), FUN=mean, na.rm=T)
		dim(quanthab2nozero)
	
	# logs ignore all values <= 0, but this is reasonable (waves shouldn't be able to gain energy)
	quanthab2cdlog = aggregate(list(KClog = log10(quantcdtrim$KCMean), Relog = log10(quantcdtrim$ReMean), klog=log10(quantcdtrim$k), Cdlog = log10(quantcdtrim$Cdval2)), by=list(hab2study = quantcdtrim$hab2study), FUN=mean, na.rm=T) # means in log10 space after adjusting for negative values
		dim(quanthab2cdlog)
	
	quanthab2sdlog = aggregate(list(KCsdlog = log10(quantcdtrim$KCMean), Resdlog = log10(quantcdtrim$ReMean), ksdlog=log10(quantcdtrim$k), Cdsdlog = log10(quantcdtrim$Cdval2)), by=list(hab2study = quantcdtrim$hab2study), FUN=sd, na.rm=T) # standard deviations in log space after adjusting for negative values
		dim(quanthab2sdlog)

	quantcdhab2 = merge(quantcdhab2, quanthab2n)
		dim(quantcdhab2)
	quantcdhab2 = merge(quantcdhab2, quanthab2sd)
		dim(quantcdhab2)
	quantcdhab2 = merge(quantcdhab2, quanthab2nozero)
		dim(quantcdhab2)
	quantcdhab2 = merge(quantcdhab2, quanthab2cdlog)
		dim(quantcdhab2)
	quantcdhab2 = merge(quantcdhab2, quanthab2sdlog)
		dim(quantcdhab2)
	quantcdhab2 = merge(quantcdhab2, quantcdtrim[!duplicated(quantcdtrim$hab2study),c('hab2study', 'Hab2', 'Hab3', 'Study')])
		dim(quantcdhab2)
	quantcdhab2$Real = TRUE
		quantcdhab2$Real[grep('Art', quantcdhab2$Hab2)] = FALSE

	# Write out aggregated data file
	write.csv(quantcdhab2, paste("output/quantcdhab2_allH", Sys.Date(), ".csv", sep=''), row.names=F)
	

## Aggregate by hab3-study

	quantcdhab3 = aggregate(list(LF=as.numeric(quantcdtrim$LF=="field"), H1 = quantcdtrim$H1, xmax=quantcdtrim$xmax, h = quantcdtrim$h, h_hv = quantcdtrim$h_hv, slope = quantcdtrim$slope, T= quantcdtrim$T, hv = quantcdtrim$hv, D = quantcdtrim$D, KC = quantcdtrim$KCMean, Re = quantcdtrim$ReMean, atten=quantcdtrim$atten, k=quantcdtrim$k, Cd = quantcdtrim$Cdval2), by=list(habstudy = quantcdtrim$habstudy), FUN=mean, na.rm=T)
		dim(quantcdhab3)
	
	quantsd = aggregate(list(H1sd = quantcdtrim$H1, xmaxsd=quantcdtrim$xmax, hsd = quantcdtrim$h, h_hvsd = quantcdtrim$h_hv, Tsd= quantcdtrim$T, hvsd = quantcdtrim$hv, Dsd = quantcdtrim$D, slopesd = quantcdtrim$slop, KCsd = quantcdtrim$KCMean, Resd = quantcdtrim$ReMean, attensd = quantcdtrim$atten, ksd=quantcdtrim$k, Cdsd = quantcdtrim$Cdval2), by=list(habstudy = quantcdtrim$habstudy), FUN=sd, na.rm=T) # standard deviations
		dim(quantsd)
	
	quantn = aggregate(list(n = rep(1, nrow(quantcdtrim))), by=list(habstudy = quantcdtrim$habstudy), FUN=sum, na.rm=T)
		dim(quantn)
	
	# treat all Cd < 0 as 0
	quantnozero = aggregate(list(Cdnozero = pmax(quantcdtrim$Cdval2, 0)), by=list(habstudy = quantcdtrim$habstudy), FUN=mean, na.rm=T)
		dim(quantnozero)
	
	# logs ignore all values <= 0, but this is reasonable (waves shouldn't be able to gain energy)
	quantcdlog = aggregate(list(KClog = log10(quantcdtrim$KCMean), Relog = log10(quantcdtrim$ReMean), klog=log10(quantcdtrim$k), Cdlog = log10(quantcdtrim$Cdval2)), by=list(habstudy = quantcdtrim$habstudy), FUN=mean, na.rm=T) # means in log10 space after adjusting for negative values
		dim(quantcdlog)
	
	quantsdlog = aggregate(list(KCsdlog = log10(quantcdtrim$KCMean), Resdlog = log10(quantcdtrim$ReMean), ksdlog=log10(quantcdtrim$k), Cdsdlog = log10(quantcdtrim$Cdval2)), by=list(habstudy = quantcdtrim$habstudy), FUN=sd, na.rm=T) # standard deviations in log space after adjusting for negative values
		dim(quantsdlog)
	
		quantcdhab3 = merge(quantcdhab3, quantn)
			dim(quantcdhab3)
		quantcdhab3 = merge(quantcdhab3, quantsd)
			dim(quantcdhab3)
		quantcdhab3 = merge(quantcdhab3, quantnozero)
			dim(quantcdhab3)
		quantcdhab3 = merge(quantcdhab3, quantcdlog)
			dim(quantcdhab3)
		quantcdhab3 = merge(quantcdhab3, quantsdlog)
			dim(quantcdhab3)
		quantcdhab3 = merge(quantcdhab3, quantcdtrim[!duplicated(quantcdtrim$habstudy),c('habstudy', 'Hab3', 'Study')])
			dim(quantcdhab3)
	
	# Write out aggregated data file
		write.csv(quantcdhab3, paste("output/quantcdhab3_allH", Sys.Date(), ".csv", sep=''), row.names=F)


############################
## Output Table of Studies 
############################

setwd("/Users/mpinsky/Documents/Stanford/NatCap/Coastal Protection Model")
quantcd = read.csv('output/quantcd_allH2012-03-25.csv')

studies = vector(mode='list', length=4)
exclude = quantcd$Study == "Bonham 1983" | quantcd$Study == "Coops et al. 1996"
sum(exclude) # should be 0, since should have removed these studies before

habs = c('Kelpfield', 'Kelplab', 'Mangrovefield', 'Marshfield', 'Marshlab', 'Seagrassfield', 'Seagrasslab')
hablabels = c('Kelp', 'Kelp', 'Mangrove', 'Marsh', 'Marsh', 'Seagrass', 'Seagrass')
lablabels = c('Field', 'Lab', 'Field', 'Field', 'Lab', 'Field', 'Lab')
for(i in 1:length(habs)){
	studies[[i]] = unique(quantcd$Study[quantcd$LFHab3 == habs[i]])
}
len=0; for(i in 1:length(habs)) len = len + length(studies[[i]])

out = data.frame(Vegetation = character(len), Lab = character(len), Study = character(len), Location = character(len), Species = character(len), n = numeric(len), H1= numeric(len), Depth = numeric(len), Decay = character(len), Drag = character(len), R2 = numeric(len), stringsAsFactors=F)

allequalna = function(x){ # are all values equal, excluding NAs?
	x = x[!is.na(x)]
	return(all(x==x[1]))
}

i = 1
for(j in 1:length(habs)){
	for(k in 1:length(studies[[j]])){
		out$Vegetation[i] = hablabels[j]
		out$Lab[i] = lablabels[j]
		out$Study[i] = as.character(studies[[j]][k])
		ind = quantcd$Study==out$Study[i] & quantcd$LFHab3 == habs[j]
		out$Location[i] = paste(unique(as.character(quantcd$Location[ind])), collapse=", ")
		if(all(quantcd$LF[ind] == 'lab')){
			out$Location[i] = ""	
		}		
		out$Species[i] = paste(unique(quantcd$Species[ind]), collapse=", ")
		out$n[i] = sum(!is.na(quantcd$Cdnozero[ind]) & quantcd$Cdnozero[ind]>0)
		out$H1[i] = signif(mean(quantcd$H1[ind]),3)
		out$Depth[i] = signif(mean(quantcd$h1[ind]),3)
		f = apply(quantcd[ind,c('h1', 'h2', 'h3', 'h4', 'h5', 'h6', 'h7', 'h8', 'h9')], MARGIN=1, FUN=allequalna)
		out$Decay[i] = signif(mean(quantcd$k[ind][f], na.rm=T), 3) # only average the k's from flat transects

		out$Drag[i] = signif(mean(quantcd$Cdnozero[ind], na.rm=T), 3)
		
		out$R2[i] = signif(mean(quantcd$R2[ind], na.rm=T),2)
		
		i = i+1
	}	
}

i = order(out$Vegetation, out$Study, out$Lab)

write.csv(out[i,], paste("output/Studies", Sys.Date(), ".csv", sep=''), row.names=F)


##################################
##                              ##
## Analyze data and make graphs ##
##                              ##
##################################
setwd("/Users/mpinsky/Documents/Stanford/NatCap/Coastal Protection Model")
source('Analysis/log10axis.R')

quantcdold = read.csv('output/quantcd_allH2012-03-25.csv')
quantcd = read.csv('output/quantcd_allH2013-03-04.csv')
quantcdlf = read.csv('output/quantcdlf_allH2012-03-25.csv')

###################################
## Broad summaries
lunique = function(x) length(unique(x))
nreal = function(x) sum(!is.na(x))
se = function(x) return(sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x))))

nrow(quantcd) # number of observations
sum(!is.na(quantcd$Cdnozero)) # number of Cd observations
lunique(quantcd$LFstudy) # number of studies
length(unique(quantcd$LFstudy[!is.na(quantcd$Cdnozero)])) # number of studies with Cd values
aggregate(quantcdlf$Cdnozero, by=list(quantcdlf$Hab3), FUN=nreal) # num studies with Cd by habitat type
mean(quantcdlf$Cdnozero, na.rm=T) # mean Cd (including zero)
aggregate(quantcdlf$Cdnozero, by=list(quantcdlf$Hab3), FUN=mean, na.rm=T) # mean by habitat
aggregate(quantcdlf$Cdnozero, by=list(quantcdlf$Hab3), FUN=sd, na.rm=T) # SD by habitat
aggregate(quantcdlf$Cdnozero, by=list(quantcdlf$Hab3), FUN=se) # SE by habitat

sum(quantcd$Cdnozero >= 0 & quantcd$Cdnozero<=10, na.rm=T)/sum(!is.na(quantcd$Cdnozero)) # proportion of Cd observations 0-10
sum(quantcd$Cdnozero > 0 & quantcd$Cdnozero<1, na.rm=T)/sum(!is.na(quantcd$Cdnozero)) # proportion of Cd observations 0-10

mean(quantcdlf$k, na.rm=T) # mean k (decay)
aggregate(quantcdlf$k, by=list(quantcdlf$Hab3), FUN=mean, na.rm=T)
aggregate(quantcdlf$k, by=list(quantcdlf$Hab3), FUN=se)

aggregate(quantcdlf$k, by=list(quantcdlf$LF), FUN=mean, na.rm=T)
aggregate(quantcdlf$k, by=list(quantcdlf$LF), FUN=se)

############################################
# Between habitat variation in parameters

boxplot(A ~ LFHab3, data=quantcdlf) # mangrove has largest diameter, kelp is next largest

boxplot(H1 ~ LFHab3, data=quantcdlf) # kelp field has highest waves, kelplab/mangrovefield/marshfield have next highest. Seagrasslab is lowest. 

##################################################################################
# Between-habitat variation in k and Cd: average values for each lab/field study

# Violin plot of Cd by habitat type and lab vs. field (Figure 2)
	# Dot shows median, thick bar shows interquartile range. Set range to 0 so no 1.5x interquartile distance thin lines.
	require(vioplot)
	quartz(width=2.9, height=2.6)
	# pdf(width=2.9, height=2.6, file=paste('Figures/Cd vioplot LF ', Sys.Date(), '.pdf', sep=''))
	# tiff(width=2.9, height=2.6, filename=paste('Figures/Cd vioplot LF ', Sys.Date(), '.tif', sep=''), units='in', res=400)

	par(cex.lab=1.0, cex.axis = 0.8, mgp=c(1.5,0.6,0), mai=c(0.53,0.45,0.03,0.03), ps=8, tcl=-0.4, las=1)
	x1 = quantcdlf$Cdnozero[quantcdlf$LFHab3 == 'Kelpfield' & !is.na(quantcdlf$Cdnozero)]
	x2 = quantcdlf$Cdnozero[quantcdlf$LFHab3 == 'Kelplab' & !is.na(quantcdlf$Cdnozero)]
	x3 = quantcdlf$Cdnozero[quantcdlf$LFHab3 == 'Mangrovefield' & !is.na(quantcdlf$Cdnozero)]
	# x4 = quantcdlf$Cdnozero[quantcdlf$LFHab3 == 'Mangrovelab' & !is.na(quantcdlf$Cdnozero)] # no points
	x5 = quantcdlf$Cdnozero[quantcdlf$LFHab3 == 'Marshfield' & !is.na(quantcdlf$Cdnozero)]
	x6 = quantcdlf$Cdnozero[quantcdlf$LFHab3 == 'Marshlab' & !is.na(quantcdlf$Cdnozero)]
	x7 = quantcdlf$Cdnozero[quantcdlf$LFHab3 == 'Seagrassfield' & !is.na(quantcdlf$Cdnozero)]
	x8 = quantcdlf$Cdnozero[quantcdlf$LFHab3 == 'Seagrasslab' & !is.na(quantcdlf$Cdnozero)]
	cols = c('dark grey', 'light grey')
	plot(0,0, xlim=c(0.5,8.5), ylim=c(0,10), xlab='', ylab='', xaxt='n')
	vioplot(x1, x3, x5, x7, at = c(1.1, 3.1, 5.1, 7.1), col=cols[1], border=cols[1], add=TRUE, colMed='black', range=0)
	vioplot(x2, x6, x8, at = c(1.9, 5.9, 7.9), col=cols[2], border=cols[2], names=rep('L', 3), add=TRUE, colMed='black', range=0)
		mtext(expression(C[d]), side=2, line=par('mgp')[1], cex=1) # do this for better control over font size
		at = 1:8+c(0.1, -0.1)
		axis(1, at=at, labels=rep(c('F', 'L'), 4), cex.axis = par('cex.axis'))
		text(x=c(1.5, 3.5, 5.5, 7.5), y=-2.8, labels=c('Kelp', 'Mangrove', 'Marsh', 'Seagrass'), xpd=NA, cex= par('cex.lab'))

	dev.off()

	# Basic analysis of Cd by habitat type or lab/field
		aggregate(list(Cd = quantcdlf$Cdnozero), by=list(Hab = quantcdlf$Hab3), FUN=mean, na.rm=T)

		summary(lm(I(log(Cdnozero))~ Hab3, data=quantcdlf)) # p = 0.43, no variation by habitat type

		summary(mod1 <- lm(I(log(Cdnozero))~ LF*Hab3, data=quantcdlf, subset = Hab3 != 'Mangrove')) # p = 0.34, no variation by habitat type or LF. No lab studies for mangrove, so eliminate
		summary(mod2 <- lm(I(log(Cdnozero))~ LF+Hab3, data=quantcdlf, subset = Hab3 != 'Mangrove')) # p = 0.27, no variation by habitat type or LF. No lab studies for mangrove, so eliminate
		summary(mod3 <- lm(I(log(Cdnozero))~ LF, data=quantcdlf, subset = Hab3 != 'Mangrove')) # p = 0.22, no variation by habitat type or LF. No lab studies for mangrove, so eliminate
		summary(mod4 <- lm(I(log(Cdnozero))~ Hab3, data=quantcdlf, subset = Hab3 != 'Mangrove')) # p = 0.32, no variation by habitat type or LF. No lab studies for mangrove, so eliminate
		summary(mod5 <- lm(I(log(Cdnozero))~ 1, data=quantcdlf, subset = Hab3 != 'Mangrove')) # null
		anova(mod1, mod2, mod3, mod4, mod5) # mod 5 is fine, not significantly worse
	
		summary(mod1<-lm(Cdnozero ~ LF + Hab3 + LF:Hab3, data=quantcdlf, subset = Hab3 != 'Mangrove'))
		summary(mod2<-lm(Cdnozero ~ LF + Hab3, data=quantcdlf, subset = Hab3 != 'Mangrove'))
		summary(mod3<-lm(Cdnozero ~ Hab3, data=quantcdlf, subset = Hab3 != 'Mangrove')) # p=0.336, no variation among habitats
		summary(mod4<-lm(Cdnozero ~ LF, data=quantcdlf, subset=Hab3 != 'Mangrove')) # p = 0.6519
		mod5 = lm(Cdnozero ~ 1, data=quantcdlf, subset=Hab3 != 'Mangrove')
		anova(mod1, mod2, mod3, mod4, mod5) # nothing looks significant
		anova(mod1, mod2) # p = 0.7797, can drop interaction
		anova(mod2, mod3) # p = 0.5363, can drop LF if Hab3 present
		anova(mod2, mod4) # p = 0.3227, can drop Hab3 if LF present
		anova(mod3, mod5) # p = 0.336, can drop Hab3
		anova(mod4, mod5) # p = 0.6519, can drop LF
		mean(quantcdlf$Cdnozero, na.rm=T)
	
		summary(mod5 <- lm(Cdnozero ~ LFHab3, data=quantcdlf)) # p = 0.8219

		summary(mod1<-lm(log10(Cdnozero) ~ LF + Hab3 + LF:Hab3, data=quantcdlf, subset = Hab3 != 'Mangrove')) # p = 0.34
		summary(mod2<-lm(log10(Cdnozero) ~ LF + Hab3, data=quantcdlf, subset = Hab3 != 'Mangrove')) # p - 0.2749
		summary(mod3<-lm(log10(Cdnozero) ~ Hab3, data=quantcdlf, subset = Hab3 != 'Mangrove')) # p=0.3296, no variation among habitats
		summary(mod4<-lm(log10(Cdnozero) ~ LF, data=quantcdlf, subset=Hab3 != 'Mangrove')) # p = 0.1736
		mod5 = lm(Cdnozero ~ 1, data=quantcdlf, subset=Hab3 != 'Mangrove')

# Basic analysis of k by habitat type or lab/field
		aggregate(list(k = quantcdlf$k), by=list(Hab = quantcdlf$Hab3), FUN=mean, na.rm=T)

		summary(mod0 <- lm(k~ Hab3, data=quantcdlf)) # p = 0.33, no variation by habitat type

		summary(mod1 <- lm(k~ LF*Hab3, data=quantcdlf, subset = Hab3 != 'Mangrove')) # p = 0.11, no variation by habitat type or LF. No lab studies for mangrove, so eliminate
		summary(mod2 <- lm(k~ LF+Hab3, data=quantcdlf, subset = Hab3 != 'Mangrove')) # p = 0.039, LF is particularyl important. No lab studies for mangrove, so eliminate
		summary(mod3 <- lm(k~ LF, data=quantcdlf, subset = Hab3 != 'Mangrove')) # p = 0.035, so variation by LF. No lab studies for mangrove, so eliminate
		summary(mod4 <- lm(k~ Hab3, data=quantcdlf, subset = Hab3 != 'Mangrove')) # p = 0.44, no variation by habitat type. No lab studies for mangrove, so eliminate
		summary(mod5 <- lm(k~ 1, data=quantcdlf, subset = Hab3 != 'Mangrove')) # null
		anova(mod2, mod3) # p=0.13, so not bad to eliminate Hab3
		anova(mod1, mod2, mod3, mod4, mod5) # mod 5 is fine, not significantly worse
		AIC(mod1, mod2, mod3, mod4, mod5) # mod 2 is best




##################################################################################
# Between-habitat variation in Cd: all values for each lab/field study (not averages)

## One graph of Cd by Re within studies
	cols = as.numeric(quantcdtrim$Study)
	plot(quantcd$ReMean, quantcdtrim$Cdnozero, log='x', col=cols)



#######################################
# Mixed-Effected Models for Cd by Re
		require(nlme)
		require(languageR)
		require(AICcmodavg)

		temp = quantcd[!is.na(quantcd$Cdnozero) & quantcd$Cdnozero>0,c('Cdnozero', 'ReMean', 'KCMean', 'Hab3', 'Hab2', 'LFstudy', 'LF')]
		temp$logCdnozero = log10(temp$Cdnozero)
		temp$logRe = log10(temp$ReMean+1)
		adj = 3e-4; temp$logReadj = log10(temp$ReMean*adj) # adjustment factor chosen to minimize slope and intercept correlation
		temp = droplevels(temp)
		
			hist(temp$logCdnozero)

		# Finding error structure
		# include lab/field as a separate effect from Hab3? FULL MODEL
		mod1v0 <- gls(logCdnozero ~ logReadj + Hab3 + logReadj:Hab3 + LF + logReadj:LF, data=temp) # allow variance to vary with logReadj within habitats
		mod1v1 <- gls(logCdnozero ~ logReadj + Hab3 + logReadj:Hab3 + LF + logReadj:LF, data=temp, weights=varPower(form=~logReadj)) # allow variance to vary with logReadj within habitats
		mod1v2 <- gls(logCdnozero ~ logReadj + Hab3 + logReadj:Hab3 + LF + logReadj:LF, data=temp, weights=varPower(form=~logReadj|Hab3)) 
		mod1v3 <- gls(logCdnozero ~ logReadj + Hab3 + logReadj:Hab3 + LF + logReadj:LF, data=temp, weights=varExp(form=~logReadj)) 
		mod1v4 <- gls(logCdnozero ~ logReadj + Hab3 + logReadj:Hab3 + LF + logReadj:LF, data=temp, weights=varExp(form=~logReadj|Hab3)) 
		mod1v5 <- gls(logCdnozero ~ logReadj + Hab3 + logReadj:Hab3 + LF + logReadj:LF, data=temp, weights=varIdent(~1|Hab3)) 
		mod1v6 <- gls(logCdnozero ~ logReadj + Hab3 + logReadj:Hab3 + LF + logReadj:LF, data=temp, weights=varConstPower(form=~logReadj|Hab3)) 
		mod1v7 <- gls(logCdnozero ~ logReadj + Hab3 + logReadj:Hab3 + LF + logReadj:LF, data=temp, weights=varExp(form=~logReadj|LFstudy), control=list(msMaxIter=200)) 
		mod1v8 <- gls(logCdnozero ~ logReadj + Hab3 + logReadj:Hab3 + LF + logReadj:LF, data=temp, weights=varPower(form=~logReadj|LFstudy), control=list(msMaxIter=200)) 
		mod1v9 <- gls(logCdnozero ~ logReadj + Hab3 + logReadj:Hab3 + LF + logReadj:LF, data=temp, weights=varIdent(~1|LFstudy), control=list(msMaxIter=200)) 
		mod1v10 <- gls(logCdnozero ~ logReadj + Hab3 + logReadj:Hab3 + LF + logReadj:LF, data=temp, weights=varExp(form=~logRe|LFstudy), control=list(msMaxIter=200)) 
			summary(mod1v0)
			summary(mod1v1)
			summary(mod1v2)
			summary(mod1v3)
			summary(mod1v4)
			summary(mod1v5)
			anova(mod1v0, mod1v1, mod1v2, mod1v3, mod1v4, mod1v5, mod1v6, mod1v7, mod1v8, mod1v9, mod1v10)
			sapply(list(mod1v0, mod1v1, mod1v2, mod1v3, mod1v4, mod1v5, mod1v6, mod1v7, mod1v8, mod1v9, mod1v10), FUN=AICc)

			hist(resid(mod1)) # ok: reasonably normal
		
			plot(fitted(mod1), resid(mod1, type='response')); abline(h=0, col='grey') # bad: variance increases with mean
			plot(residuals(mod1, type='normalized') ~ fitted(mod1), data=temp)
			coplot(residuals(mod1, type='normalized') ~ fitted(mod1) | Hab3, data=temp)
			coplot(residuals(mod1v7, type='response') ~ fitted(mod1v7) | Hab3, data=temp, col=temp$LFstudy)
			coplot(residuals(mod1v7, type='normalized') ~ fitted(mod1v7) | Hab3, data=temp, col=temp$LFstudy)
			plot(temp$logReadj, resid(mod1)); abline(h=0, col='grey') # variance decreases with logReadj?

		# Now select the random effects
		#mod1 <- lmer(logCdnozero ~ logRe + Hab3 + logRe:Hab3 + (logRe|LFstudy), data=temp, REML=TRUE) # habitats vary by intercept and slope, plus lab and field studies vary by intercept and slope,
		mod1 <- lme(logCdnozero ~ logReadj + Hab3 + LF + logReadj:Hab3 + logReadj:LF, random = ~1+logReadj|LFstudy, data=temp, weights=varExp(form=~logRe|LFstudy), control=list(maxIter=200, msMaxIter=300, msMaxEval=500, msVerbose=F, opt='optim')) 
			summary(mod1)
			plot(fitted(mod1), resid(mod1, type='response')); abline(h=0, col='grey') # bad: variance decreases with mean
			coplot(residuals(mod1, type='response') ~ fitted(mod1) | Hab3, data=temp, col=temp$LF) # bad for Marsh and maybe Kelp, but these aren't the right resids to plot
			coplot(residuals(mod1, type='normalized') ~ fitted(mod1) | Hab3, data=temp, col=temp$LFstudy)
			plot(residuals(mod1, type='normalized') ~  Hab3, data=temp)
			plot(residuals(mod1, type='normalized') ~  LF, data=temp)
			coplot(residuals(mod1, type='normalized') ~ logReadj | Hab3, data=temp, col=temp$LFstudy)

		mod2 <- lme(logCdnozero ~ logReadj + Hab3 + logReadj:Hab3 + LF + logReadj:LF, random = ~1|LFstudy, data=temp, weights=varExp(form=~logRe|LFstudy), control=list(msMaxIter=200)) # only random intercepts

		mod3 <- gls(logCdnozero ~ logReadj + Hab3 +  LF + logReadj:Hab3 + logReadj:LF, data=temp, weights=varExp(form=~logRe|LFstudy)) # no random effects

			a=anova(mod1, mod2, mod3);a
		
			# corrected p-values for testing on the boundary (Zuur Ch. 4 p. 123)
			0.5 * (1-pchisq(a$L.Ratio[3],1)) # for mod2 vs. mod3 (include random intercept?)
			0.5 * ((1 - pchisq(a$L.Ratio[2], 1)) + (1 - pchisq(a$L.Ratio[2], 2))) # mod1 vs. mod2 (random slope?)

		# And finally, the fixed effects
		mod1 <- lme(logCdnozero ~ logReadj + Hab3 + LF + logReadj:Hab3 + logReadj:LF, random = ~1+logReadj|LFstudy, data=temp, weights=varExp(form=~logRe|LFstudy), control=list(msMaxIter=200, opt='optim'), method='ML') 
			summary(mod1)

				mod1a = lme(logCdnozero ~ logReadj + Hab3 + LF + logReadj:Hab3, random = ~1+logReadj|LFstudy, data=temp, weights=varExp(form=~logRe|LFstudy), control=list(msMaxIter=200, opt='optim'), method='ML') 
					summary(mod1a)
				mod1b = update(mod1, .~. - logReadj:Hab3, method='ML')
					summary(mod1b)
		
				anova(mod1, mod1a) # drop logReadj:LF
				anova(mod1, mod1b) # keep logReadj:Hab3 for the moment

		mod2 <- lme(logCdnozero ~ logReadj + Hab3 + LF + logReadj:Hab3, random = ~1+logReadj|LFstudy, data=temp, weights=varExp(form=~logRe|LFstudy), control=list(msMaxIter=200), method='ML') 
			summary(mod2)
	
			mod2a = update(mod2, .~. - logReadj:Hab3)
			mod2b = update(mod2, .~. - LF)

			anova(mod2, mod2a) # keep logReadj:Hab3
			anova(mod2, mod2b) # drop LF

		mod3 <- lme(logCdnozero ~ logReadj + Hab3 + logReadj:Hab3, random = ~1+logReadj|LFstudy, data=temp, weights=varExp(form=~logRe|LFstudy), control=list(msMaxIter=200), method='ML') 
			summary(mod3)
	
			mod3a = update(mod3, .~. - logReadj:Hab3) # keep logReadj:Hab3
			mod3b = update(mod3, .~. - Hab3)

			anova(mod3, mod3a)
			anova(mod3, mod3b)
		
			AIC(mod1, mod1a, mod1b, mod2, mod2a, mod2b, mod3, mod3a)
			BIC(mod1, mod1a, mod1b, mod2, mod2a, mod2b, mod3, mod3a)
			sapply(list(mod1, mod1a, mod1b, mod2, mod2a, mod2b, mod3, mod3a), FUN=AICc)
			
			coplot(residuals(mod2, type='normalized') ~ fitted(mod3) | Hab3, data=temp, col=temp$LFstudy)
			plot(residuals(mod1, type='normalized') ~  Hab3, data=temp)
			plot(residuals(mod1, type='normalized') ~  LF, data=temp)
			coplot(residuals(mod1, type='normalized') ~ logRe | Hab3, data=temp, col=temp$LFstudy)
			
		# Final model
		modfin <- lme(logCdnozero ~ logReadj + Hab3 + logReadj:Hab3 , random = ~1+logReadj|LFstudy, data=temp, weights=varExp(form=~logRe|LFstudy), method='REML', control = list(msMaxIter=200))
		modfin.ML <- lme(logCdnozero ~ logReadj + Hab3 + logReadj:Hab3 , random = ~1+logReadj|LFstudy, data=temp, weights=varExp(form=~logRe|LFstudy), method='ML', control = list(msMaxIter=200))
			summary(modfin)
	
			coplot(residuals(modfin, type='response') ~ fitted(modfin) | LFstudy, data=temp, col=temp$LFstudy)
			coplot(residuals(modfin, type='normalized') ~ fitted(modfin) | LFstudy, data=temp, col=temp$LFstudy)
			plot(residuals(modfin, type='response') ~ fitted(modfin), data=temp, col=temp$LFstudy)
			plot(residuals(modfin, type='normalized') ~ fitted(modfin), data=temp, col=temp$LFstudy)
			plot(residuals(modfin, type='normalized') ~ logReadj, data=temp, col=temp$LFstudy)
			plot(residuals(modfin, type='normalized') ~ Hab3, data=temp, col=temp$LFstudy)
			plot(residuals(modfin, type='normalized') ~ LF, data=temp, col=temp$LFstudy)
		
			# Print intercepts and slopes as a table
			params = data.frame(hab = c('Kelp', 'Mangrove', 'Marsh', 'Seagrass'), int = numeric(4), slope = numeric(4))
			params$int[1] = fixef(modfin)[1]			
			params$slope[1] = fixef(modfin)[2]			
			for(i in 2:4){
				params$int[i] = fixef(modfin)[1] + fixef(modfin)[1+i]
				params$slope[i] = fixef(modfin)[2] + fixef(modfin)[4+i]
			}	
			params
		
		# Plot the logRe*Hab3 model in five graphs: one for each habitat and one of all regressions
		quartz(width=6.1, height=4.5) # Ecosphere 2 columns is 6.1"
		# pdf(width=6.1, height=4.2, file=paste('Figures/Cd by Re ME models, 5 graphs_', Sys.Date(), '.pdf', sep=''))
		# tiff(filename=paste('Figures/Cd by Re ME models, 5 graphs_', Sys.Date(), '.tiff', sep=''), width=6.1, height=4.2, units='in', res = 400, compression='lzw')
		pchs = c(1,1,1,1)
		cols = c(rgb(166,206,227, maxColorValue=255), rgb(31,120,180, maxColorValue=255), rgb(178,223,138, maxColorValue=255), rgb(51,160,44, maxColorValue=255)) # only color-blind safe ColorBrewer2 set of 4 colors (qualitative)
		colsLF = c('black', 'grey') # for lab vs. field
		ltyLF = list(1,"12")
		cex.lab=1.0 # for axis labels
		cex.dot = 0.7 # for plotting
		cex.axis = 0.8 # for axis annotation
		habs = levels(temp$Hab3)
		ylabs = c(expression(C[d]), '', expression(C[d]), '') 
		xlabs = c('', '', 'Re', 'Re')
		sublabs = c('A', 'B', 'C', 'D', 'E')
		padj = 0.1
		adj = 0.05
		ps = 8; # point size for text
		tcl1 = -0.4 # tick length for big ticks (fraction of a line of text)
		tcl2 = -0.2 # small ticks

		layout(matrix(c(1,2,6,3,4,5), byrow=TRUE, nrow=2))
		par(mai=c(0.45, 0.45, 0.15, 0.028), mgp=c(1.25, 0.5 ,0), ps = ps, cex=1, las=1) # want 2 points (2/72") on all sides. Set par() after layout() because the latter adjusts par (especially cex)
		for(i in 1:4){ # Plot all points for each habitat in separate graphs
			k = temp$Hab3 == habs[i]
			xlims = range(temp$logRe[k]) # include storm Re
			if(floor(xlims[1]) == floor(xlims[2])){ # make sure xlims include at least one label so we can see the x-axis
				xlims[1] = floor(xlims[1])
				xlims[2] = ceiling(xlims[2])
			}
			plot(temp$logRe[k], temp$logCd[k], pch=pchs[i], ylab="", xlab="", main='', col=colsLF[temp$LF[k]], cex=cex.dot, xaxt='n', yaxt='n', xlim=xlims, cex.lab = cex.lab, ps=ps)
			mtext(ylabs[i], side=2, line=par('mgp')[1]) # do this for better control over font size
			mtext(xlabs[i], side=1, ps=ps, line=par('mgp')[1])
			log10axis(side=1, at=seq(1,8,by=1), log=FALSE, cex.axis=cex.axis, tcl1 = tcl1, tcl2=tcl2)
			log10axis(side=2, at=seq(-3,2,by=1), log=FALSE, cex.axis=cex.axis, tcl1 = tcl1, tcl2=tcl2)
			studs = unique(temp$LFstudy[k])
			for(j in 1:length(studs)){ # add the lines between points from the same study
				kk = temp$LFstudy == studs[j]
				x = temp$logRe[kk]
				y = temp$logCd[kk]
				ind = order(x)
				lines(x[ind], y[ind], col=colsLF[temp$LF[kk]], lwd=0.5)
			}
			#points(stormRe[i], max(temp$logCd[k], na.rm=T), pch=16, cex=2, col='red') # add storm conditions for each habitat
			mtext(sublabs[i], side=3, padj= padj, adj=adj, cex=cex.lab, ps=ps)
		}
		#legend('bottomleft', col=colsLF, lwd=1, legend=levels(temp$LF), cex=0.9, pch=1)

		# Plot fits in the bottom graph
		xlims = c(0.5,5.5)
		ylims = c(-3.5,3) # for log
		#ylims = c(-6,2.4) # w/out Hab3
			# Fixed effect fits
		f0 = fitted(modfin, level=0)
		plot(0,0, col='white', xlim=xlims, ylim=ylims, ylab='', xlab='', main='', xaxt='n', yaxt='n', cex=cex.dot)
			mtext('Re', side=1, ps=ps, line=par('mgp')[1])
			log10axis(side=1, at=seq(1,8,by=2), log=FALSE, cex.axis=cex.axis, tcl1 = tcl1, tcl2=tcl2)
			log10axis(side=2, at=seq(-3,6,by=2), log=FALSE, cex.axis=cex.axis, tcl1 = tcl1, tcl2=tcl2)
		mtext(sublabs[5], side=3, padj= padj, adj= adj)
		habs = sort(unique(temp$Hab3))
		lfs = sort(unique(temp$LF))
		for(i in 1:4){ # the habitat effects
			inds = temp$Hab3 == habs[i]
			if(sum(inds)>0){
				x = temp$logRe[inds]
				y = f0[inds]
				k = order(x)
				points(x[k],y[k],lwd=3, pch=pchs[i], col=cols[i], type='l', lty=1)
			}
		}
			# Random effect fits
		f1 = fitted(modfin, level=1)
		studs = unique(names(f1))
		fits = vector('list', length(studs))
		for(i in 1:length(studs)){
			inds = temp$LFstudy == studs[i]
			habind = which(habs %in% unique(temp$Hab3[inds]))
			if(sum(inds)>0){
				x = temp$logRe[inds]
				y = f1[inds]
				k = order(x)
				points(x[k],y[k],lwd=1, pch=pchs[habind], col=cols[habind], type='l')
			}
			fits[[i]] = list(Re = x, Cd = y)

		}
		names(fits) = studs
		
		# Plot legend for part E
		plot(1,1, bty='n', xlab='', ylab='', col='white', xaxt='n', yaxt='n')
		legend('bottom', col=cols, lwd=2, legend=levels(temp$Hab3), cex=cex.axis, seg.len=2, bty='n')
	
	
		dev.off()
		
		# compare ME fit to observed to get R2
		f = fixef(modfin)
			fint = f[c(1,3,4,5)]
			fslope = f[c(2,6,7,8)]
		ran = ranef(modfin)$LFstudy
		studs = levels(temp$LFstudy)
		par(mfrow=c(4,6), mai=c(0.5, 0.5, 0.15, 0.05), mgp=c(2,1,0))
		for(i in 1:length(studs)){
			j = temp$LFstudy == studs[i]
			thishab = unique(temp$Hab3[temp$LFstudy == studs[i]])
			f_ind = which(habs %in% thishab) # index into fixed effects
			r = temp$logReadj[j]
			ind = order(r)
			r = r[ind]
			cd = temp$logCdnozero[j][ind]
			y = fint[1] + fint[f_ind] + ran[i,1] + (fslope[1] + fslope[f_ind] + ran[i,2])*(r+log10(2e-4)) 
			if(f_ind == 1) y = fint[1] + ran[i,1] + (fslope[1] + ran[i,2])*(r- log10(3e4)) 
			sum = summary(lm(y ~ cd))
			print(paste(studs[i], ', r2=', signif(sum$r.squared,3), ', n=', sum(j), sep=''))
			plot(cd, y, xlab='Observed', ylab = 'Predicted', main=studs[i], cex.main=0.7)
		}
		

	detach(package:lme4)


#######################################
# Mixed-Effected Models for Cd by KC
		require(nlme)
		require(languageR)
		require(AICcmodavg)

		temp = quantcd[!is.na(quantcd$Cdnozero) & quantcd$Cdnozero>0,c('Cdnozero', 'ReMean', 'KCMean', 'KC2Mean', 'Hab3', 'Hab2', 'LFstudy', 'LF')]
		temp$logCdnozero = log10(temp$Cdnozero)
		temp$logKC = log10(temp$KC2Mean+1)
		adj = 3e-4; temp$logKCadj = log10(temp$KC2Mean*adj) # adjustment factor chosen to minimize slope and intercept correlation
		temp = droplevels(temp)
		
			hist(temp$logCdnozero)

		# Finding error structure
		# include lab/field as a separate effect from Hab3? FULL MODEL
		mod1v0 <- gls(logCdnozero ~ logKCadj + Hab3 + logKCadj:Hab3 + LF + logKCadj:LF, data=temp) # allow variance to vary with logKCadj within habitats
		mod1v1 <- gls(logCdnozero ~ logKCadj + Hab3 + logKCadj:Hab3 + LF + logKCadj:LF, data=temp, weights=varPower(form=~logKCadj)) # allow variance to vary with logKCadj within habitats
		mod1v2 <- gls(logCdnozero ~ logKCadj + Hab3 + logKCadj:Hab3 + LF + logKCadj:LF, data=temp, weights=varPower(form=~logKCadj|Hab3)) 
		mod1v3 <- gls(logCdnozero ~ logKCadj + Hab3 + logKCadj:Hab3 + LF + logKCadj:LF, data=temp, weights=varExp(form=~logKCadj)) 
		mod1v4 <- gls(logCdnozero ~ logKCadj + Hab3 + logKCadj:Hab3 + LF + logKCadj:LF, data=temp, weights=varExp(form=~logKCadj|Hab3)) 
		mod1v5 <- gls(logCdnozero ~ logKCadj + Hab3 + logKCadj:Hab3 + LF + logKCadj:LF, data=temp, weights=varIdent(~1|Hab3)) 
		mod1v6 <- gls(logCdnozero ~ logKCadj + Hab3 + logKCadj:Hab3 + LF + logKCadj:LF, data=temp, weights=varConstPower(form=~logKCadj|Hab3)) 
		mod1v7 <- gls(logCdnozero ~ logKCadj + Hab3 + logKCadj:Hab3 + LF + logKCadj:LF, data=temp, weights=varExp(form=~logKCadj|LFstudy), control=list(msMaxIter=200)) 
		mod1v8 <- gls(logCdnozero ~ logKCadj + Hab3 + logKCadj:Hab3 + LF + logKCadj:LF, data=temp, weights=varPower(form=~logKCadj|LFstudy), control=list(msMaxIter=200)) 
		mod1v9 <- gls(logCdnozero ~ logKCadj + Hab3 + logKCadj:Hab3 + LF + logKCadj:LF, data=temp, weights=varIdent(~1|LFstudy), control=list(msMaxIter=200)) 
		mod1v10 <- gls(logCdnozero ~ logKCadj + Hab3 + logKCadj:Hab3 + LF + logKCadj:LF, data=temp, weights=varExp(form=~logKC|LFstudy), control=list(msMaxIter=200)) 
			summary(mod1v0)
			summary(mod1v1)
			summary(mod1v2)
			summary(mod1v3)
			summary(mod1v4)
			summary(mod1v5)
			anova(mod1v0, mod1v1, mod1v2, mod1v3, mod1v4, mod1v5, mod1v6, mod1v7, mod1v8, mod1v9, mod1v10)
			sapply(list(mod1v0, mod1v1, mod1v2, mod1v3, mod1v4, mod1v5, mod1v6, mod1v7, mod1v8, mod1v9, mod1v10), FUN=AICc)

		mod1 = mod1v10

			hist(resid(mod1)) # ok: reasonably normal
		
			plot(fitted(mod1), resid(mod1, type='response')); abline(h=0, col='grey') # bad: trend in residuals with fitted
			plot(residuals(mod1, type='normalized') ~ fitted(mod1), data=temp)
			coplot(residuals(mod1, type='normalized') ~ fitted(mod1) | Hab3, data=temp)
			coplot(residuals(mod1v7, type='response') ~ fitted(mod1v7) | Hab3, data=temp, col=temp$LFstudy)
			coplot(residuals(mod1v7, type='normalized') ~ fitted(mod1v7) | Hab3, data=temp, col=temp$LFstudy)
			plot(temp$logKCadj, resid(mod1)); abline(h=0, col='grey') # variance increases with logKCadj?

		# Now select the random effects
		#mod1 <- lmer(logCdnozero ~ logKC + Hab3 + logKC:Hab3 + (logKC|LFstudy), data=temp, REML=TRUE) # habitats vary by intercept and slope, plus lab and field studies vary by intercept and slope,
		mod1 <- lme(logCdnozero ~ logKCadj + Hab3 + LF + logKCadj:Hab3 + logKCadj:LF, random = ~1+logKCadj|LFstudy, data=temp, weights=varExp(form=~logKC|LFstudy), control=list(maxIter=200, msMaxIter=300, msMaxEval=500, msVerbose=F)) 
			summary(mod1)
			plot(fitted(mod1), resid(mod1, type='response')); abline(h=0, col='grey') # bad: variance increases with mean
			coplot(residuals(mod1, type='response') ~ fitted(mod1) | Hab3, data=temp, col=temp$LF) # some are bad, but these aren't the right resids to plot
			coplot(residuals(mod1, type='normalized') ~ fitted(mod1) | Hab3, data=temp, col=temp$LFstudy)
			plot(residuals(mod1, type='normalized') ~  Hab3, data=temp)
			plot(residuals(mod1, type='normalized') ~  LF, data=temp)
			coplot(residuals(mod1, type='normalized') ~ logKCadj | Hab3, data=temp, col=temp$LFstudy)

		mod2 <- lme(logCdnozero ~ logKCadj + Hab3 + logKCadj:Hab3 + LF + logKCadj:LF, random = ~1|LFstudy, data=temp, weights=varExp(form=~logKC|LFstudy), control=list(msMaxIter=200)) # only random intercepts

		mod3 <- gls(logCdnozero ~ logKCadj + Hab3 +  LF + logKCadj:Hab3 + logKCadj:LF, data=temp, weights=varExp(form=~logKC|LFstudy)) # no random effects

			a=anova(mod1, mod2, mod3);a
		
			# corrected p-values for testing on the boundary (Zuur Ch. 4 p. 123)
			0.5 * (1-pchisq(a$L.Ratio[3],1)) # for mod2 vs. mod3 (include random intercept?)
			0.5 * ((1 - pchisq(a$L.Ratio[2], 1)) + (1 - pchisq(a$L.Ratio[2], 2))) # mod1 vs. mod2 (random slope?)

		# And finally, the fixed effects
		mod1 <- lme(logCdnozero ~ logKCadj + Hab3 + LF + logKCadj:Hab3 + logKCadj:LF, random = ~1+logKCadj|LFstudy, data=temp, weights=varExp(form=~logKC|LFstudy), control=list(msMaxIter=200), method='ML') 
			summary(mod1)

				mod1a = update(mod1, .~. - logKCadj:LF, method='ML')
					summary(mod1a)
				mod1b = update(mod1, .~. - logKCadj:Hab3, method='ML')
					summary(mod1b)
		
				anova(mod1, mod1a) # keep logKCadj:LF p = 3e-04
				anova(mod1, mod1b) # drop logKCadj:Hab3 p = 0.0425

		mod2 <- lme(logCdnozero ~ logKCadj + Hab3 + LF + logKCadj:LF, random = ~1+logKCadj|LFstudy, data=temp, weights=varExp(form=~logKC|LFstudy), control=list(msMaxIter=200), method='ML') 
			summary(mod2)
	
			mod2a = update(mod2, .~. - logKCadj:LF)
			mod2b = update(mod2, .~. - Hab3)

			anova(mod2, mod2a) # keep logKCadj:LF p = 0.0158
			anova(mod2, mod2b) # drop Hab3 p = 0.531

		mod3 <- lme(logCdnozero ~ logKCadj + LF + logKCadj:LF, random = ~1+logKCadj|LFstudy, data=temp, weights=varExp(form=~logKC|LFstudy), control=list(msMaxIter=200), method='ML') 
			summary(mod3)
	
			mod3a = update(mod3, .~. - logKCadj:LF) # keep logKCadj:LF
			mod3b = update(mod3, .~. - LF)

			anova(mod3, mod3a) # p = 0.02 keep logKCadj:LF
			anova(mod3, mod3b) # p = 0.0118 keep LF
		
			AIC(mod1, mod1a, mod1b, mod2, mod2a, mod2b, mod3, mod3a, mod3b)
			BIC(mod1, mod1a, mod1b, mod2, mod2a, mod2b, mod3, mod3a)
			sapply(list(mod1, mod1a, mod1b, mod2, mod2a, mod2b, mod3, mod3a, mod3b), FUN=AICc)
			
			coplot(residuals(mod2, type='normalized') ~ fitted(mod3) | Hab3, data=temp, col=temp$LFstudy)
			plot(residuals(mod1, type='normalized') ~  Hab3, data=temp)
			plot(residuals(mod1, type='normalized') ~  LF, data=temp)
			coplot(residuals(mod1, type='normalized') ~ logKC | Hab3, data=temp, col=temp$LFstudy)
			
		# Final model for KC2
		modfinKC <- lme(logCdnozero ~ logKCadj + LF + Hab3 + logKCadj:LF + logKCadj:Hab3, random = ~1+logKCadj|LFstudy, data=temp, weights=varExp(form=~logKC|LFstudy), method='REML', control = list(msMaxIter=200))
		modfinKC.ML <- lme(logCdnozero ~ logKCadj + LF + Hab3 + logKCadj:LF + logKCadj:Hab3, random = ~1+logKCadj|LFstudy, data=temp, weights=varExp(form=~logKC|LFstudy), method='ML', control = list(msMaxIter=200))
			summary(modfinKC)
	
			AIC(modfin, modfinKC)
			BIC(modfin, modfinKC)

			AIC(modfin, modfinKC, mod1) # include full model for KC
			BIC(modfin, modfinKC, mod1) # include full model for KC

			AIC(modfin.ML, modfinKC.ML)
			BIC(modfin, modfinKC)
	
			coplot(residuals(modfinKC, type='response') ~ fitted(modfinKC) | LFstudy, data=temp, col=temp$LFstudy)
			coplot(residuals(modfinKC, type='normalized') ~ fitted(modfinKC) | LFstudy, data=temp, col=temp$LFstudy)
			plot(residuals(modfinKC, type='response') ~ fitted(modfinKC), data=temp, col=temp$LFstudy)
			plot(residuals(modfinKC, type='normalized') ~ fitted(modfinKC), data=temp, col=temp$LFstudy)
			plot(residuals(modfinKC, type='normalized') ~ logKCadj, data=temp, col=temp$LFstudy)
			plot(residuals(modfinKC, type='normalized') ~ Hab3, data=temp, col=temp$LFstudy)
			plot(residuals(modfinKC, type='normalized') ~ LF, data=temp, col=temp$LFstudy)
		
			# Print intercepts and slopes as a table (if logKCadj:LF included)
			params = data.frame(hab = c('Field', 'Lab'), int = numeric(2), slope = numeric(2))
			params$int[1] = fixef(modfinKC)[1]
			params$slope[1] = fixef(modfinKC)[2]			
			params$int[2] = fixef(modfinKC)[1] + fixef(modfinKC)[3]
			params$slope[2] = fixef(modfinKC)[2] + fixef(modfinKC)[4]

			params
			
			
	# Alternative fixed-effects choice using AIC
	mod1 <- lme(logCdnozero ~ logKCadj + Hab3 + LF + logKCadj:Hab3 + logKCadj:LF, random = ~1+logKCadj|LFstudy, data=temp, weights=varExp(form=~logKC|LFstudy), control=list(msMaxIter=200), method='ML') 
	mod2 <- lme(logCdnozero ~ logKCadj + Hab3 + LF + logKCadj:Hab3, random = ~1+logKCadj|LFstudy, data=temp, weights=varExp(form=~logKC|LFstudy), control=list(msMaxIter=200), method='ML') 
	mod3 <- lme(logCdnozero ~ logKCadj + Hab3 + LF, random = ~1+logKCadj|LFstudy, data=temp, weights=varExp(form=~logKC|LFstudy), control=list(msMaxIter=200), method='ML') 
	mod4 <- lme(logCdnozero ~ logKCadj + Hab3, random = ~1+logKCadj|LFstudy, data=temp, weights=varExp(form=~logKC|LFstudy), control=list(msMaxIter=200), method='ML') 
	mod5 <- lme(logCdnozero ~ logKCadj + LF, random = ~1+logKCadj|LFstudy, data=temp, weights=varExp(form=~logKC|LFstudy), control=list(msMaxIter=200), method='ML') 
	mod6 <- lme(logCdnozero ~ Hab3 + LF, random = ~1+logKCadj|LFstudy, data=temp, weights=varExp(form=~logKC|LFstudy), control=list(msMaxIter=200), method='ML') 
	mod7 <- lme(logCdnozero ~ logKCadj, random = ~1+logKCadj|LFstudy, data=temp, weights=varExp(form=~logKC|LFstudy), control=list(msMaxIter=200), method='ML') 
	mod8 <- lme(logCdnozero ~ Hab3, random = ~1+logKCadj|LFstudy, data=temp, weights=varExp(form=~logKC|LFstudy), control=list(msMaxIter=200), method='ML') 
	mod9 <- lme(logCdnozero ~ LF, random = ~1+logKCadj|LFstudy, data=temp, weights=varExp(form=~logKC|LFstudy), control=list(msMaxIter=200), method='ML') 
	mod10 <- lme(logCdnozero ~ logKCadj + Hab3 + logKCadj:Hab3, random = ~1+logKCadj|LFstudy, data=temp, weights=varExp(form=~logKC|LFstudy), control=list(msMaxIter=200), method='ML') 
	mod11 <- lme(logCdnozero ~ logKCadj + Hab3 + LF + logKCadj:LF, random = ~1+logKCadj|LFstudy, data=temp, weights=varExp(form=~logKC|LFstudy), control=list(msMaxIter=200), method='ML') 
	mod12 <- lme(logCdnozero ~ logKCadj + LF + logKCadj:LF, random = ~1+logKCadj|LFstudy, data=temp, weights=varExp(form=~logKC|LFstudy), control=list(msMaxIter=200), method='ML') 
		
	AIC(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12)
	BIC(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12)


##########################
##### Helpful functions
##########################

printMeanSE = function(x, group){
	norig = table(group)
	n = table(group, is.finite(x))[,2]
	m = aggregate(x, by=list(group), mean, na.rm=T)
	names(m) = c("group", "mean")
	mn = aggregate(x, by=list(group), min, na.rm=T)
	mx = aggregate(x, by=list(group), max, na.rm=T)
	rg = round(log10(mx$x)-log10(mn$x),2)
	s = aggregate(x, by=list(group), sd, na.rm=T)
	names(s)[2] = "sd"
	cv = s$sd/abs(m$mean)
	out = cbind(norig, n, m$mean, mn$x, mx$x, rg, s$sd, cv)
	colnames(out) = c("n_orig", "n", "mean", "min", "max", "log10range", "sd", "cv")
	print(out)
}

# Errorbar function
errbar.xy <- function (x, y, xl=NA, xu=xl, yl=NA, yu = yl, length = 0.08, plot=FALSE, xlim=c(NA,NA), ylim=c(NA,NA),...){
	# xl, yl, xu, and yu are distances
	xlim2 = c(min(x-xl, na.rm=TRUE), 1.1*max(x+xl, na.rm=TRUE))
	ylim2 = c(min(y-yl, na.rm=TRUE), 1.1*max(y+yl, na.rm=TRUE))
	if(xlim2[1]<xlim[1] | is.na(xlim[1])) xlim[1] <- xlim2[1]
	if(xlim2[2]>xlim[2] | is.na(xlim[2])) xlim[2] <- xlim2[2]
	if(ylim2[1]<ylim[1] | is.na(ylim[1])) ylim[1] <- ylim2[1]
	if(ylim2[2]>ylim[2] | is.na(ylim[2])) ylim[2] <- ylim2[2]

	if(plot){
		plot(x,y, xlim=xlim, ylim=ylim, ...)	
	}
    if(!all(is.na(yl))){
    	arrows(x, y + yu, x, y - yl, angle = 90, code = 3, length = length, ...)
	}
	if(!all(is.na(xl))){
    	arrows(x+xu, y, x-xl, y, angle = 90, code = 3, length = length, ...)
	}
}

# Function to find f, given all the other params (use optimise(molleroptim, interval=c(_,), H1=_, ...))
# Shallow water approximation for wave celerity
molleroptim = function(f, H1, H2, x, vk, T, h, h1, h2, K, d, g, verbose = FALSE){
	# one of the inputs can be >1 for sensitivity analysis
	L = T*sqrt(g*h) # wave length (linear wave theory, shallow water approximation) ??? May not be accurate, since doesn't match molleroptim2 or 3
	c1 = sqrt(g*h1) # wave celerity (Airy wave theory, shallow water approximation), CERC 1984 Eq. 2-9
	c2 = sqrt(g*h2)

	# Viscous friction
	a1 = (2*pi/L)^2/((4*pi*h/L+sinh(4*pi*h/L))*sqrt(pi/(T*vk)))
	Kv = (1+a1*x)^-1

	# Shoaling
	n1 = (1/2)*(1+(4*pi*h1/L)/(sinh(4*pi*h1/L)))
	n2 = (1/2)*(1+(4*pi*h2/L)/(sinh(4*pi*h2/L)))
	Ks = sqrt(n1*c1/n2/c2)

	# Percolation
	a3 = K*(2*pi/T)/vk*(4*pi/L)/(4*pi*h/L+sinh(4*pi*h/L))*tanh(2*pi*d/L)
	Kp = (1+a3*x)^-1

	# Friction
	a2 = 64*pi^3*f*H1*h^2*Ks^2/(3*g^2*h^2*T^4*(sinh(2*pi*h/L))^3)
	Kf = (1+a2*x)^-1
	
	H2pred = H1*Kv*Ks*Kp*Kf	

	if(verbose){
		print(cbind(c("L", "c1", "c2", "a1", "Kv", "n1", "n2", "Ks", "a3", "Kp", "a2", "Kf", "H2pred"), c(L, c1, c2, a1, Kv, n1, n2, Ks, a3, Kp, a2, Kf, H2pred)))	
	}	
	
	return(abs(H2-H2pred))	
}

# Function to find f, given all the other params (use optimise(molleroptim, interval=c(_,), H1=_, ...))
# Uses full version of Airy equation to find celerity
# one of the inputs can be of length >1 for sensitivity analysis
molleroptim2 = function(f, H1, H2, x, vk, T, h, h1, h2, K, d, g, verbose=FALSE){
	L1 = g*T^2/(2*pi)*sqrt(tanh(4*pi^2*h1/(T^2*g))) # CERC 1984, Eq. 2-4b, accurate to w/in 5%
	L2 = g*T^2/(2*pi)*sqrt(tanh(4*pi^2*h2/(T^2*g)))
	L = (L1+L2)/2

	#c1 = sqrt(g*L1/(2*pi)*tanh((2*pi*h1)/L1)) # wave celerity (Airy wave theory), CERC 1984 Eq. 2-2
	#c2 = sqrt(g*L2/(2*pi)*tanh((2*pi*h2)/L2))

	c1 = L1/T # celerity at depth 1. From CERC 1984 Eq. 2-1
	c2 = L2/T

	# Viscous friction
	a1 = (2*pi/L)^2/((4*pi*h/L+sinh(4*pi*h/L))*sqrt(pi/(T*vk)))
	Kv = (1+a1*x)^-1

	# Shoaling
	n1 = (1/2)*(1+(4*pi*h1/L)/(sinh(4*pi*h1/L)))
	n2 = (1/2)*(1+(4*pi*h2/L)/(sinh(4*pi*h2/L)))
	Ks = sqrt(n1*c1/n2/c2)

	# Percolation
	a3 = K*(2*pi/T)/vk*(4*pi/L)/(4*pi*h/L+sinh(4*pi*h/L))*tanh(2*pi*d/L)
	Kp = (1+a3*x)^-1

	# Friction
	a2 = 64*pi^3*f*H1*h^2*Ks^2/(3*g^2*h^2*T^4*(sinh(2*pi*h/L))^3)
	Kf = (1+a2*x)^-1
	
	H2pred = H1*Kv*Ks*Kp*Kf		
	
	if(verbose){
		print(cbind(c("L", "c1", "c2", "a1", "Kv", "n1", "n2", "Ks", "a3", "Kp", "a2", "Kf", "H2pred"), c(L, c1, c2, a1, Kv, n1, n2, Ks, a3, Kp, a2, Kf, H2pred)))	
	}	
	
	return(abs(H2-H2pred))	
}

# Function to find f, given all the other params (use optimise(molleroptim, interval=c(_,), H1=_, ...))
# Uses full version of Airy equation to find celerity
# Uses linear dispersion relation and iterative Newton-Raphson method to find wave length
# one of the inputs can be of length >1 for sensitivity analysis
molleroptim3 = function(f, H1, H2, x, vk, T, h, h1, h2, K, d, g, verbose=FALSE){
	L1 = 2*pi/func_disp(T, h1, g) # wave length at depth 1 (linear dispersion relation for wavenumber k)
	L2 = 2*pi/func_disp(T, h2, g) 
	L = (L1+L2)/2

	#c1 = sqrt(g*L1/(2*pi)*tanh((2*pi*h1)/L1)) # wave celerity (Airy wave theory), CERC 1984 Eq. 2-2
	#c2 = sqrt(g*L2/(2*pi)*tanh((2*pi*h2)/L2))

	c1 = L1/T # celerity at depth 1. From CERC 1984 Eq. 2-1
	c2 = L2/T

	# Viscous friction
	a1 = (2*pi/L)^2/((4*pi*h/L+sinh(4*pi*h/L))*sqrt(pi/(T*vk)))
	Kv = exp(-a1*x) # Eq. 4 in Moller et al. 1999

	# Shoaling
	n1 = (1/2)*(1+(4*pi*h1/L)/(sinh(4*pi*h1/L))) # Eq. 3 in Moller et al. 1999
	n2 = (1/2)*(1+(4*pi*h2/L)/(sinh(4*pi*h2/L)))
	Ks = sqrt(n1*c1/(n2*c2))

	# Percolation
	a3 = K*(2*pi/T)/vk*(4*pi/L)/(4*pi*h/L+sinh(4*pi*h/L))*tanh(2*pi*d/L) # Eq. 7 in Moller et al. 1999
	Kp = (1+a3*x)^-1 # Eq. 7 in Moller et al. 1999

	# Friction
	a2 = 64*pi^3*f*H1*h^2*Ks^2/(3*g^2*h^2*T^4*(sinh(2*pi*h/L))^3) # Eq. 5 in Moller et al. 1999
	Kf = (1+a2*x)^-1
	
	H2pred = H1*Kv*Ks*Kp*Kf	
	
	Hb = L * 0.142*tanh(2*pi*h/L) # Miche's formula for limiting wave height due to breaking. Eq. 8 in Moller et al. 1999

	if(H2pred > Hb) H2pred = Hb
		
	if(verbose){
		print(cbind(c("L", "c1", "c2", "a1", "Kv", "n1", "n2", "Ks", "a3", "Kp", "a2", "Kf", "H2pred"), c(L, c1, c2, a1, Kv, n1, n2, Ks, a3, Kp, a2, Kf, H2pred)))	
	}	

	return(abs(H2-H2pred))	
}


#This function computes wave number based on water depth and wave period
#using the linear dispersion relation and an iterative Newton-Raphson
#technique (from Greg, via Katie Arkema)
func_disp = function(T,h,g){

	#Define Constants
	X = numeric(0)
	sigma=2*pi/T
	k=sigma^2/g*(tanh(sigma^2*h/g))^(-0.5) #rough approx of k
	X = c(X, k*h)

	cst=sigma^2*h/g

	#Initialize loop
	i=1
	f=cst-X[1]*tanh(X[1]);
	r=-tanh(X[1])-X[1]*(1-tanh(X[1])^2)

	#Loop to determine f

	while(abs(f)>=1e-3){
    	i=i+1
    	X[i]=X[1]-f/r;
    	f=cst-X[i]*tanh(X[i]);
    	r=-tanh(X[i])-X[i]*(1-tanh(X[i])^2)
	}

	k2=X[i]/h
	return(k2)
}

