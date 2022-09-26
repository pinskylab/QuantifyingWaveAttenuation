# Analyze SppData

setwd("/Users/mpinsky/Documents/Stanford/NatCap/Coastal Protection Model")

sppdata = read.csv("data/SppData 110704.csv")


# make a histogram
spp = sppdata[!duplicated(sppdata[c("Spp", "Habitat")]),c("Spp", "Habitat")]
	spp = spp[order(spp$Habitat, spp$Spp),]
habcols = list(Seagrass="green", Marsh="yellow", Kelp="brown", Mangrove="orange")
quartz(width=6, height=7.5)
par(mfrow=c(length(spp$Spp),3), mai=c(0.2, 0.2, 0.1, 0), mgp=c(1.7,.7,0), cex=0.3, cex.main=1.5)
maxA = max(sppdata$A, na.rm=T)
maxD = max(sppdata$D, na.rm=T)
maxhv = max(sppdata$hv, na.rm=T)

for(i in 1:length(spp$Spp)){
	if(sum(!is.na(sppdata$A[sppdata$Spp==spp$Spp[i]]))>0) hist(sppdata$A[sppdata$Spp==spp$Spp[i]], main=spp$Spp[i], xlab="stem width (m)", xlim=c(0,maxA), breaks=seq(0, maxA, length.out=40), col=habcols[[as.character(spp$Habitat[i])]])
	else plot(0,0, bty="n", xlab="", ylab="", xaxt="n", yaxt="n", col="white")
	if(sum(!is.na(sppdata$D[sppdata$Spp==spp$Spp[i]]))>0) hist(sppdata$D[sppdata$Spp==spp$Spp[i]], main=spp$Spp[i], xlab="density (shoots/m2)", xlim=c(0, maxD), breaks=seq(0, maxD, length.out=40), col=habcols[[as.character(spp$Habitat[i])]])
	else plot(0,0, bty="n", xlab="", ylab="", xaxt="n", yaxt="n", col="white")
	if(sum(!is.na(sppdata$hv[sppdata$Spp==spp$Spp[i]]))>0) hist(sppdata$hv[sppdata$Spp==spp$Spp[i]], main=spp$Spp[i], xlab="veg height (m)", xlim=c(0, maxhv), breaks=seq(0, maxhv, length.out=40), col=habcols[[as.character(spp$Habitat[i])]])
	else plot(0,0, bty="n", xlab="", ylab="", xaxt="n", yaxt="n", col="white")
}

## Aggregate by study and then by species

sppdata$SpRef = paste(sppdata$Spp, sppdata$Ref, sep="")

study = aggregate(list(A=sppdata$A, D=sppdata$D, hv=sppdata$hv, xmax = sppdata$Bed), by=list(SpRef = sppdata$SpRef, Spp = sppdata$Spp, Habitat=sppdata$Habitat), FUN=mean, na.rm=T)
dim(study)

sppsd = aggregate(list(Asd=study$A, Dsd=study$D, hvsd=study$hv, xmaxsd=study$xmax), by=list(Spp = study$Spp, Habitat=study$Habitat), FUN=sd, na.rm=T)
spp = aggregate(list(A=study$A, D=study$D, hv=study$hv, xmax=study$xmax), by=list(Spp = study$Spp, Habitat=study$Habitat), FUN=mean, na.rm=T)
dim(sppsd); dim(spp)

spp = merge(spp, sppsd, by=c("Habitat", 'Spp'))
names(spp)
dim(spp)
spp

## Add references
spp$Refs = NA
for(i in 1:nrow(spp)){
	j = which(as.character(sppdata$Spp) == as.character(spp$Spp[i]))
	spp$Refs[i] = paste(unique(sppdata$Ref[j]), collapse=';')
}

# Make aggregate entries (Spartina spp.,  Salicornia spp., and Rhizophora spp.)
i = grep("Spartina", spp$Spp)
temp1 = data.frame(Spp="Spartina spp.", Habitat='Marsh', A = mean(spp$A[i], na.rm=T), Asd = sd(spp$A[i], na.rm=T), D = mean(spp$D[i], na.rm=T), Dsd = sd(spp$D[i], na.rm=T),  hv = mean(spp$hv[i], na.rm=T), hvsd = sd(spp$hv[i], na.rm=T), xmax=mean(spp$xmax[i], na.rm=T), xmaxsd = sd(spp$xmax[i], na.rm=T), Refs=NA)

i = grep("Salicornia", spp$Spp)
temp2 = data.frame(Spp="Salicornia spp.", Habitat='Marsh', A = mean(spp$A[i], na.rm=T), Asd = sd(spp$A[i], na.rm=T), D = mean(spp$D[i], na.rm=T), Dsd = sd(spp$D[i], na.rm=T),  hv = mean(spp$hv[i], na.rm=T), hvsd = sd(spp$hv[i], na.rm=T), xmax=mean(spp$xmax[i], na.rm=T), xmaxsd = sd(spp$xmax[i], na.rm=T), Refs=NA)

i = grep("Rhizophora", spp$Spp)
j = grep("Avicennia", spp$Spp)
i = c(i,j)
temp3 = data.frame(Spp="Avicennia spp., Rhizophora spp.", Habitat='Mangrove', A = mean(spp$A[i], na.rm=T), Asd = sd(spp$A[i], na.rm=T), D = mean(spp$D[i], na.rm=T), Dsd = sd(spp$D[i], na.rm=T),  hv = mean(spp$hv[i], na.rm=T), hvsd = sd(spp$hv[i], na.rm=T), xmax=mean(spp$xmax[i], na.rm=T), xmaxsd = sd(spp$xmax[i], na.rm=T), Refs=NA)

spp = rbind(spp, temp1)
spp = rbind(spp, temp2)
spp = rbind(spp, temp3)
dim(spp)


## Write out

write.csv(spp, paste("output/SppSummary_", Sys.Date(), ".csv", sep=""))



#######################################################
############## Analyze by vegetation category #########
#######################################################

spp = read.csv('output/SppSummary_2011-01-02.csv')
source('Analysis/log10axis.R')

# pairs plot
pairs(~ log10(A) + log10(D) + log10(hv) + log10(xmax), data=spp, panel=panel.smooth)
	# A negatively correlated to D

# boxplot
par(mfrow=c(2,2))
	boxplot(A ~ Habitat, data=spp, log='y', ylab="Stem width")
	boxplot(D ~ Habitat, data=spp, log='y', ylab="Density")
	boxplot(hv ~ Habitat, data=spp, log='y', ylab="Height")
	boxplot(xmax ~ Habitat, data=spp, log='y', ylab="Bed width (m)")
	
# regressions
	pchs = c(0,16,4,17)[as.numeric(spp$Habitat)]
	cols = c('black', 'black', 'black', 'grey')[as.numeric(spp$Habitat)]

	plot(log10(A*100) ~ log10(D), data=spp, pch = pchs, col=cols, xaxt='n', yaxt='n', xlab=expression(paste("Density (shoots/", m^2, ")", sep="")), ylab='Diameter (cm)')
		log10axis(side=1, at=seq(-2,6,by=1), log=FALSE)
		log10axis(side=2, at=seq(-4,2,by=1), log=FALSE)
		legend('bottomleft', pch=c(0,16,4,17), col=c('black', 'black', 'black', 'grey'), legend=c('Kelp', 'Mangrove', 'Marsh', 'Seagrass'), cex=0.7)

		mod = lm(log10(A) ~ log10(D), data=spp)
		abline(mod)
		summary(mod)
			# log10(A) = -1.18782 - 0.42102*log10(D), p = 1.6e-5, r2 = 0.5787