## Analyze output from Fit_Cd (all spp) (as output by Matlab script)

########################################################
## Initial dataprep (can skip once done, see below)
 
# Read in and merge w/ verbose if needed
	setwd("/Users/mpinsky/Documents/Stanford/NatCap/Coastal Protection Lit Review")
#	fitcdfile = "analysis/Fit_Cd 101008/VegData2Out_04-Jan-2011_allH.csv" # all observations, inside or outside the veg bed
	fitcdfile = "analysis/Fit_Cd 101008/VegData2Out_04-Jan-2011_trimoutside.csv" # no obs outside the veg bed
	fitcd = read.csv(fitcdfile, header=T) # output from MATLAB
	quant = read.csv("output/QuantforGregVerbose_2011-01-04.csv") # verbose study data
	dim(fitcd)

	fitcd = merge(fitcd, quant[,c("id", "Study", "Location", "LF", "Hab2", "Species", "Description", "WaveTypeVerbose")], by.x="id", by.y="id")
	dim(fitcd)
	
	# turn -999 to NA
	fitcd$Cd[fitcd$Cd == -999] = NA
	fitcd$Cdf[fitcd$Cdf == -999] = NA
	fitcd$Cds[fitcd$Cds == -999] = NA
	fitcd$R2Cd[fitcd$R2Cd == -999] = NA
	fitcd$K1[fitcd$K1 == -999] = NA
	fitcd$R2K1[fitcd$R2K1 == -999] = NA
	fitcd$slo[fitcd$slo == -999] = NA
	
	# find average depth and max distance
	fitcd$have = rowMeans(fitcd[, c('h1', 'h2', 'h3', 'h4', 'h5','h6', 'h7','h8')], na.rm=T)
	fitcd$xmax = apply(fitcd[, c('x1', 'x2', 'x3', 'x4', 'x5','x6', 'x7','x8')], MARGIN=1, FUN=max, na.rm=T)

	# Create Hab3
	fitcd$Hab3 = fitcd$Hab2
	fitcd$Hab3[fitcd$Hab3=="Seagrass (Art.)"] = "Seagrass"
	fitcd$Hab3[fitcd$Hab3=="Marsh (Art.)"] = "Marsh"
	fitcd$Hab3[fitcd$Hab3=="Kelp (Art.)"] = "Kelp"

	# calc exponential decay coeff from endpoints (H1 and Hx) SHOULD TURN THIS INTO A BEST FIT LINE
	fitcd$k = NA
	for(i in 1:nrow(fitcd)){
		Hx = fitcd[i, c('H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8')]
		ind = which(!is.na(Hx))
		indx = max(ind)+1 # index of last wave height
		Hx = fitcd[i, paste('H', indx, sep='')] 
		x = fitcd[i, paste('x', indx, sep='')]
		if(fitcd$x1[i]>0 & length(ind)>1){ # don't use first wave height if it's outside of veg and we have >1 inside
			print(i)
			H1 = fitcd$H2[i]
			x = x-fitcd$x1[i]
		} else {
			H1 = fitcd$H1[i]	
		}
		fitcd$k[i] = -log(Hx/H1)/x # the exponent in H2/H1 = exp(-ax)
	}

	fitcd$SpStudy = paste(fitcd$Species, fitcd$Study, sep="")

	fitcd = fitcd[order(fitcd$X.kk),] # re-order by kk
	n = names(fitcd); n1 = which(n=="id"); n2 = which(n == "X.kk")
	newn = n; newn[n1] = "X.kk"; newn[n2] = "id" # switch position of id and kk (put kk first)
	fitcd = fitcd[,newn]
	
	if(length(grep("trimoutside", fitcdfile))>0){
		write.csv(fitcd, paste("output/Fit_Cd_verbose_trimoutside", Sys.Date(), ".csv", sep=''), row.names=F)
	}
	if(length(grep("allH", fitcdfile))>0){ 
		write.csv(fitcd, paste("output/Fit_Cd_verbose_allH", Sys.Date(), ".csv", sep=''), row.names=F)
	}	
	
	
############################################################
# Read in merged file and aggregate by hab3-study
############################################################

#	fitcd = read.csv('output/Fit_Cd_verbose_trimoutside2011-01-26.csv', header=T); fitcdfile = "trimoutside"
	fitcd = read.csv('output/Fit_Cd_verbose_allH2011-01-26.csv', header=T); fitcdfile = "allH"

# Average by hab3-study
fitcdag = aggregate(list(k=fitcd$k, H1 = quant$H1, xmax=quant$xmax, LF=as.numeric(quant$LF=="field")), by=list(habstudy = quant$habstudy), FUN=mean, na.rm=T)


# Sample size
	length(unique(fitcd$SpStudy)) # 36 total studies (by species)
	length(unique(paste(fitcd$Hab3, fitcd$Study))) # 26 total studies (by habitat-type)		
		sort(unique(paste(fitcd$Hab3, fitcd$Study)))
	table(fitcd$Hab2[!duplicated(fitcd$SpStudy)])
		# Bare            Kelp     Kelp (Art.)        Mangrove           Marsh    Marsh (Art.) 
        #    4               1               2               4              13               1 
        # Seagrass Seagrass (Art.) 
         #   7               4 
	length(unique(fitcd$SpStudy[fitcd$Hab2 != "Bare"])) # 32 studies (excluding bare)
	length((fitcd$SpStudy[fitcd$Hab2 != "Bare"])) # 592 data points (excluding bare)

	length(unique(fitcd$SpStudy[fitcd$Hab2 != "Bare" & !is.na(fitcd$Cd)])) # 28 studies (excluding bare)
	length((fitcd$SpStudy[fitcd$Hab2 != "Bare" & !is.na(fitcd$Cd)])) # 478 data points (excluding bare)

	j = c('Kelp', 'Kelp (Art.)', 'Marsh', 'Marsh (Art.)', 'Seagrass', 'Seagrass (Art.)')
	j = c('Mangrove', 'Marsh', 'Marsh (Art.)', 'Seagrass', 'Seagrass (Art.)')
	length(unique(fitcd$SpStudy[fitcd$Hab2 %in% j])) # 28 marsh/seagrass/kelp studies, 29 marsh/seagrass/mangrove

	# Studies missing from fitcd
	j = c('Mangrove', 'Marsh', 'Marsh (Art.)', 'Seagrass', 'Seagrass (Art.)')
	setdiff(quant$spstudy[quant$Hab2 %in% j], fitcd$SpStudy[fitcd$Hab2 %in% j])
		# NONE

###############
# Plots and Analysis

### Histograms
	# Cd: Plot histograms by habitat of all data points (no Inf)
i = is.finite(fitcd$Cd)
hab = as.character(unique(fitcd$Hab2[i]))
minCd = min(fitcd$Cd[i], na.rm=T)
maxCd = max(fitcd$Cd[i], na.rm=T)
bks = seq(minCd, maxCd, length.out=20)
quartz(width=8, height=3)
par(mfrow=c(2, 4), mai=c(0.4, 0.4, 0.3, 0.05), mgp=c(1.8, 0.6, 0))
hist(fitcd$Cd[i], breaks=bks, main="All Habitats", col="grey", xlab="Cd")
for(j in hab){
	hist(fitcd$Cd[i][fitcd$Hab2[i]==j], breaks=bks, main=j, col="grey", xlab="Cd")
}

	# log10(Cd): Plot histograms by habitat of all data points (no Inf and no < 0)
i = is.finite(fitcd$Cd) & fitcd$Cd>0
hab = sort(as.character(unique(fitcd$Hab2[i])))
minCd = min(fitcd$Cd[i], na.rm=T)
maxCd = max(fitcd$Cd[i], na.rm=T)
quartz(width=8, height=3)
par(mfrow=c(2, 4), mai=c(0.4, 0.4, 0.3, 0.05), mgp=c(1.8, 0.6, 0))
hist(log10(fitcd$Cd[i]), breaks=seq(log10(minCd), log10(maxCd), length.out=20), main="All Habitats", col="grey", xlab="log10 Cd")
for(j in hab){
	hist(log10(fitcd$Cd[i][fitcd$Hab2[i]==j]), breaks=seq(log10(minCd), log10(maxCd), length.out=20), main=j, col="grey", xlab="log10 Cd")
}

	# log10(Cd) Pretty plot for Arkema PPT: Plot histograms by habitat of all data points (no Inf)
i = is.finite(fitcd$Cd)
hab = c('Marsh', 'Seagrass', 'Mangrove')
minCd = min(fitcd$Cd[i & fitcd$Hab2 %in% hab], na.rm=T)
maxCd = max(fitcd$Cd[i & fitcd$Hab2 %in% hab], na.rm=T)
quartz(width=4, height=6)
# jpeg(width=4, height=6, units='in', res=300, file=paste('Figures/Cd hist_', Sys.Date(), '.jpg', sep=''))
par(mfrow=c(3, 1), mai=c(0.4, 0.4, 0.3, 0.05), mgp=c(1.8, 0.6, 0))
for(j in 1:length(hab)){
	hist(log10(fitcd$Cd[i & fitcd$Hab2 %in% hab[j]]), breaks=seq(log10(minCd), log10(maxCd), length.out=20), main=hab[j], 	col="dark grey", xlab="Cd", xaxt = 'n')
	axis(1, at=-4:1, labels = parse(text = paste("10^", -4:1, sep = "")))
}

dev.off()

	table(fitcd$Hab2[i & fitcd$Hab2 %in% hab])


	# K1: Plot histograms by habitat of all data points (no Inf)
i = is.finite(fitcd$K1)
hab = sort(as.character(unique(fitcd$Hab2[i])))
minK = min(fitcd$K1[i], na.rm=T)
maxK = max(fitcd$K1[i], na.rm=T)
par(mfrow=c(2, 4), mai=c(0.4, 0.4, 0.3, 0.05), mgp=c(1.8, 0.6, 0))
hist(fitcd$K1[i], breaks=seq(minK, maxK, length.out=20), main="All Habitats", col="grey", xlab="K")
for(j in hab){
	hist(fitcd$K1[i][fitcd$Hab2[i]==j], breaks=seq(minK, maxK, length.out=20), main=j, col="grey", xlab="K")
}

	# log10(K1+o): Plot histograms by habitat of all data points (no Inf)
i = is.finite(fitcd$K1)
hab = sort(as.character(unique(fitcd$Hab2[i])))
o = 0.001
o =1
minK = min(fitcd$K1[i]+o, na.rm=T)
maxK = max(fitcd$K1[i]+o, na.rm=T)
par(mfrow=c(2, 3), mai=c(0.4, 0.4, 0.3, 0.05), mgp=c(1.8, 0.6, 0))
hist(log10(fitcd$K1[i]+o), breaks=seq(log10(minK), log10(maxK), length.out=20), main="All Habitats", col="grey", xlab="K")
for(j in hab){
	hist(log10(fitcd$K1[i][fitcd$Hab2[i]==j]+o), breaks=seq(log10(minK), log10(maxK), length.out=20), main=j, col="grey", xlab="K")
}

	# log10(K1+o) Pretty plot for Arkema PPT: Plot histograms by habitat of all data points (no Inf)
i = is.finite(fitcd$K1)
o = 0.005 # offset
hab = c('Marsh', 'Seagrass', 'Mangrove')
minK = min(fitcd$K1[i & fitcd$Hab2 %in% hab], na.rm=T)
maxK = max(fitcd$K1[i & fitcd$Hab2 %in% hab], na.rm=T)
quartz(width=4, height=6)
# jpeg(width=4, height=6, units='in', res=300, file=paste('Figures/K1 hist_', Sys.Date(), '.jpg', sep=''))
par(mfrow=c(3, 1), mai=c(0.4, 0.4, 0.3, 0.05), mgp=c(1.8, 0.6, 0))
for(j in 1:length(hab)){
	hist(log10(fitcd$K1[i & fitcd$Hab2 %in% hab[j]]+o), breaks=seq(log10(minK+o), log10(maxK+o), length.out=10), main=hab[j], 	col="dark grey", xlab="K", xaxt = 'n')
	lab = c(-3, -2, -1)
	axis(1, at=log10(c(0,10^lab)+o), labels = parse(text = c("0", paste("10^", lab, sep = ""))))
}

dev.off()

	# log10(k+o) Pretty plot for Arkema PPT: Plot histograms by habitat of all data points (no Inf)
i = is.finite(fitcd$k)
o = 0.01 # offset
hab = c('Marsh', 'Seagrass', 'Mangrove', 'Kelp')
minK = min(fitcd$k[i & fitcd$Hab2 %in% hab], na.rm=T)
maxK = max(fitcd$k[i & fitcd$Hab2 %in% hab], na.rm=T)
quartz(width=4, height=6)
# jpeg(width=4, height=6, units='in', res=300, file=paste('Figures/k hist_', Sys.Date(), '.jpg', sep=''))
par(mfrow=c(4, 1), mai=c(0.4, 0.4, 0.3, 0.05), mgp=c(1.8, 0.6, 0))
for(j in 1:length(hab)){
	hist(log10(fitcd$k[i & fitcd$Hab2 %in% hab[j]]+o), breaks=seq(log10(minK+o), log10(maxK+o), length.out=10), main=hab[j], 	col="dark grey", xlab="K", xaxt = 'n')
	lab = c(-2, -1, 0)
	axis(1, at=log10(c(0,10^lab)+o), labels = parse(text = c("0", paste("10^", lab, sep = ""))))
}

mod = lm(k ~ Hab2, data=fitcd, subset=i & fitcd$Hab2 %in% hab)
summary(mod)
table(fitcd$Hab2[i & fitcd$Hab2 %in% hab])

dev.off()


##### Regressions

	# Cd vs. T, D, A, hv, H1, h, and x
par(mfrow=c(2,4))
i = is.finite(fitcd$Cd) # finite and not missing
habs = sort(as.character(unique(fitcd$Hab2[i])))
cols = c('black', 'red', 'blue', 'green', 'orange')
ylims = log10(range(fitcd$Cd[i], na.rm=T))
plot(fitcd$T[i][fitcd$Hab2[i]==habs[j]], log10(fitcd$Cd[i][fitcd$Hab2[i]==habs[j]]), col=cols[j], pch=16, xlab="T", ylab="log10 Cd", ylim=ylims, xlim=range(fitcd$T[i], na.rm=T))
	for(j in 2:length(habs)) points(fitcd$T[i][fitcd$Hab2[i]==habs[j]], log10(fitcd$Cd[i][fitcd$Hab2[i]==habs[j]]), col=cols[j], pch=16)
plot(fitcd$D[i][fitcd$Hab2[i]==habs[j]], log10(fitcd$Cd[i][fitcd$Hab2[i]==habs[j]]), col=cols[j], pch=16, xlab="D", ylab="log10 Cd", ylim=ylims, xlim=range(fitcd$D[i], na.rm=T))
	for(j in 2:length(habs)) points(fitcd$D[i][fitcd$Hab2[i]==habs[j]], log10(fitcd$Cd[i][fitcd$Hab2[i]==habs[j]]), col=cols[j], pch=16)
plot(fitcd$A[i][fitcd$Hab2[i]==habs[j]], log10(fitcd$Cd[i][fitcd$Hab2[i]==habs[j]]), col=cols[j], pch=16, xlab="A", ylab="log10 Cd", ylim=ylims, xlim=range(fitcd$A[i], na.rm=T))
	for(j in 2:length(habs)) points(fitcd$A[i][fitcd$Hab2[i]==habs[j]], log10(fitcd$Cd[i][fitcd$Hab2[i]==habs[j]]), col=cols[j], pch=16)
plot(fitcd$hv[i][fitcd$Hab2[i]==habs[j]], log10(fitcd$Cd[i][fitcd$Hab2[i]==habs[j]]), col=cols[j], pch=16, xlab="hv", ylab="log10 Cd", ylim=ylims, xlim = c(0,1.5)) # xlim=range(fitcd$hv[i], na.rm=T))
	for(j in 2:length(habs)) points(fitcd$hv[i][fitcd$Hab2[i]==habs[j]], log10(fitcd$Cd[i][fitcd$Hab2[i]==habs[j]]), col=cols[j], pch=16)
plot(fitcd$H1[i][fitcd$Hab2[i]==habs[j]], log10(fitcd$Cd[i][fitcd$Hab2[i]==habs[j]]), col=cols[j], pch=16, xlab="H1", ylab="log10 Cd", ylim=ylims, xlim = c(0,0.4)) # xlim=range(fitcd$H1[i], na.rm=T))
	for(j in 2:length(habs)) points(fitcd$H1[i][fitcd$Hab2[i]==habs[j]], log10(fitcd$Cd[i][fitcd$Hab2[i]==habs[j]]), col=cols[j], pch=16)
plot(fitcd$have[i][fitcd$Hab2[i]==habs[j]], log10(fitcd$Cd[i][fitcd$Hab2[i]==habs[j]]), col=cols[j], pch=16, xlab="h", ylab="log10 Cd", ylim=ylims, xlim=c(0,1.6)) # xlim=range(fitcd$have[i], na.rm=T))
	for(j in 2:length(habs)) points(fitcd$have[i][fitcd$Hab2[i]==habs[j]], log10(fitcd$Cd[i][fitcd$Hab2[i]==habs[j]]), col=cols[j], pch=16)
plot(fitcd$xmax[i][fitcd$Hab2[i]==habs[j]], log10(fitcd$Cd[i][fitcd$Hab2[i]==habs[j]]), col=cols[j], pch=16, xlab="x", ylab="log10 Cd", ylim=ylims, xlim=range(fitcd$xmax[i], na.rm=T))
	for(j in 2:length(habs)) points(fitcd$xmax[i][fitcd$Hab2[i]==habs[j]], log10(fitcd$Cd[i][fitcd$Hab2[i]==habs[j]]), col=cols[j], pch=16)
legend('topright', legend=habs, col=cols, pch=16)


# Cd vs. xmax, h, h/hv, D, Re, and KC
	require(RColorBrewer)
	# display.brewer.all()
	# piles data
	piles = read.csv("analysis/Fit_Cd 101008/Greg Piles data/DDp227Waves_CdRe.txt", skip=1, header=F, col.names=c('Re', 'Cd'))
	pilescem = read.csv("analysis/Fit_Cd 101008/Greg Piles data/CEMp571Waves_Smooth.txt", skip=1, header=F, col.names=c('Re', 'Cd'))

	# Data Prep
	fitcdtrim = fitcd[fitcd$Hab2 != "Bare" & !is.na(fitcd$Cd), ]
	fitcdtrim$Hab2 = fitcdtrim$Hab2[drop=TRUE]
	fitcdtrim$Hab3 = fitcdtrim$Hab3[drop=TRUE]
	fitcdtrim$habstudy = paste(fitcdtrim$Hab3, fitcdtrim$Study, sep='')
	
	y = fitcdtrim$Cd
	dropy = y<=0 & !is.na(y)
	which(dropy)
	fitcdtrim$X.kk[dropy]
	zero = min(y[y>0], na.rm=T)/2
	y[dropy] = zero
	fitcdtrim$Cdtrim = y
	
	habstdcd = aggregate(fitcdtrim[,c('k', 'Cd', 'ReMean', 'KCMean', 'xmax', 'have', 'hv', 'D')], by=list(habstudy = fitcdtrim$habstudy), FUN=mean, na.rm=T)
	dim(habstdcd)
	habstdcd = merge(habstdcd, fitcdtrim[!duplicated(fitcdtrim$habstudy), c("habstudy", "Species", "Hab2", "Hab3", "Study", "Field")])
	habstdcd$Cdtrim = habstdcd$Cd
	dropy = habstdcd$Cdtrim<=0 & !is.na(habstdcd$Cdtrim)
	habstdcd$Cdtrim[dropy] = zero

	# Function run regressions and plot lines
	regrcd = function(hab3, x, y, cols, newx, plot=F, xlog=F, zero){
		# fit linear models in log-log space (same as y = a*x^b)
		for(i in 1:length(unique(hab3))){
			hab = levels(hab3)[i]
			inds = hab3 == hab & !is.na(y) & y > zero
			data = data.frame(x = x[inds], y = y[inds])
			if(nrow(data)>0){
				mod = lm(I(log(y)) ~ x, data=data)
				if(xlog) mod = lm(I(log(y)) ~ I(log(x)), data=data)
				print(paste(hab, ": n = ", mod$df+3, ", r2 = ", signif(summary(mod)$r.squared,3), ", p = ", signif(summary(mod)$coefficients[2,4],3), sep=''))
				if(plot & !xlog){
					lines(newx, exp(predict(mod, newdata=list(x=newx))), col=cols[i], lty=3, lwd=1)
				}
				if(plot & xlog){
					lines(newx, exp(predict(mod, newdata=list(x=newx))), col=cols[i], lty=3, lwd=1)
				}
			}
		}
	}

	#cols = c('black', 'red4', 'seagreen4', 'sienna1', 'steelblue1', 'slategray', 'black', 'black') # use this 10/15/2010 and later
	cols = brewer.pal(4, name="Dark2")
	cols = brewer.pal(4, name="Paired")
	pchs = c(16, 1, 16, 16, 1, 16, 1)
	cex=0.8
	newx = exp(seq(1,log(max(fitcdtrim$ReMean, na.rm=T)), length.out=50))
	xlims = c(1, max(c(fitcdtrim$ReMean, piles$Re, pilescem$Re), na.rm=T))
	ylims = range(y, na.rm=T)


	# quartz(width=4, height=7.5)
	# pdf(file=paste("analysis/Fit_Cd 101008/Figures/CdvFactors_RealvArt", Sys.Date(), ".pdf", sep=""), width=6, height=6)
	# pdf(file=paste("analysis/Fit_Cd 101008/Figures/CdvFactors_RealvArt_wRegr", Sys.Date(), ".pdf", sep=""), width=5, height=7.5)
	layout(matrix(data=c(1,2,3,4,5,6,7,7,7,7), nrow=5, byrow=T), widths=rep(1.5,2), heights=rep(1.8,5))
	par(mai=c(0.5, 0.5, 0.05, 0.05), mgp = c(2.3,0.8,0))

	plot(habstdcd$xmax, habstdcd$Cdtrim, log='xy', pch=pchs[habstdcd$Hab2], col=cols[habstdcd$Hab3], xlab="Distance (m)", ylab="Cd", yaxt='n', cex=cex)
		at=c(log10(min(y, na.rm=T)), seq(-7, 3, by=4)); axis(2, at=10^at, labels=c("< 0", parse(text=paste("10^", at[2:length(at)], sep=""))), las=1)
		print("Distance")
		regrcd(hab3 = habstdcd$Hab3, x= habstdcd$xmax, y= habstdcd$Cdtrim, cols=cols, newx=newx, plot=F, xlog=T, zero=zero)
	plot(habstdcd$have, habstdcd$Cdtrim, log='xy', pch=pchs[habstdcd$Hab2], col=cols[habstdcd$Hab3], xlab="Depth (m)", ylab="", yaxt='n', cex=cex)
		at=c(log10(min(y, na.rm=T)), seq(-7, 3, by=4)); axis(2, at=10^at, labels=c("< 0", parse(text=paste("10^", at[2:length(at)], sep=""))), las=1)
		print("Depth")
		regrcd(hab3 = habstdcd$Hab3, x= habstdcd$have, y= habstdcd$Cdtrim, cols=cols, newx=newx, plot=F, xlog=T, zero=zero)
	plot(habstdcd$have/habstdcd$hv, habstdcd$Cdtrim, log='xy', pch=pchs[habstdcd$Hab2], col=cols[habstdcd$Hab3], xlab="Relative depth", ylab="Cd", yaxt='n', cex=cex)
		at=c(log10(min(y, na.rm=T)), seq(-7, 3, by=4)); axis(2, at=10^at, labels=c("< 0", parse(text=paste("10^", at[2:length(at)], sep=""))), las=1)
		print("Relative Depth")
		regrcd(hab3 = habstdcd$Hab3, x= habstdcd$have/habstdcd$hv, y= habstdcd$Cdtrim, cols=cols, newx=newx, plot=F, xlog=F, zero=zero)
	plot(habstdcd$D, habstdcd$Cdtrim, log='xy', pch=pchs[habstdcd$Hab2], col=cols[habstdcd$Hab3], xlab="Density (stems/m2)", ylab="", yaxt='n', xaxt='n', cex=cex)
		at=c(log10(min(y, na.rm=T)), seq(-7, 3, by=4)); axis(2, at=10^at, labels=c("< 0", parse(text=paste("10^", at[2:length(at)], sep=""))), las=1)
		at=seq(0,5,by=2); axis(1, at=10^at, labels=parse(text=paste("10^", at, sep="")))
		print("Density")
		regrcd(hab3 = habstdcd$Hab3, x= habstdcd$D, y= habstdcd$Cdtrim, cols=cols, newx=newx, plot=F, xlog=T, zero=zero)
	plot(habstdcd$KCMean, habstdcd$Cdtrim, log='xy', pch=pchs[habstdcd$Hab2], col=cols[habstdcd$Hab3], xlab="KC number", ylab="Cd", yaxt='n', xaxt='n', cex=cex)
		at=c(log10(min(y, na.rm=T)), seq(-7, 3, by=4)); axis(2, at=10^at, labels=c("< 0", parse(text=paste("10^", at[2:length(at)], sep=""))), las=1)
		at=seq(0,5,by=2); axis(1, at=10^at, labels=parse(text=paste("10^", at, sep="")))
		print("KC number")
		regrcd(hab3 = habstdcd$Hab3, x= habstdcd$KCMean, y= habstdcd$Cdtrim, cols=cols, newx=newx, plot=F, xlog=T, zero=zero)
	plot(habstdcd$ReMean, habstdcd$Cdtrim, log='xy', pch=pchs[habstdcd$Hab2], col=cols[habstdcd$Hab3], xlab="Reynolds Number", ylab="", xlim=xlims, yaxt='n', xaxt='n', cex=cex, ylim=ylims)
		at=c(log10(min(y, na.rm=T)), seq(-7, 3, by=4)); axis(2, at=10^at, labels=c("< 0", parse(text=paste("10^", at[2:length(at)], sep=""))), las=1)
		at=c(0,2,4,6); axis(1, at=10^at, labels=parse(text=paste("10^", at, sep="")))
		print("Re by study")
		regrcd(hab3 = habstdcd$Hab3, x=habstdcd$ReMean, y=habstdcd$Cdtrim, cols=cols, newx=newx, plot=F, xlog=T, zero=zero)


	# Main plot
	plot(fitcdtrim$ReMean, y, log='xy', pch=pchs[fitcdtrim$Hab2], col=cols[fitcdtrim$Hab3], xlab="Reynolds Number", ylab="Cd", xlim=xlims, yaxt='n', xaxt='n')

	points(pilescem$Re, pilescem$Cd, col="dark grey", pch=16) # piles data from CEM
	points(piles$Re, piles$Cd, col="black", pch=16) # piles data from D&D
	lines(c(3.3e4, 7.6e5), c(1000, 1000), col="black", lwd=2, lty=3) # range of Re for piles in storm waves of 0.1-2.1 m depth (H = 0.08-1.64)

	legend('bottomleft', legend=c(levels(fitcdtrim$Hab2), "Piles (D&D)", "Piles (CEM)"), pch=c(pchs, 16, 16), col = c(rep(cols, c(2,1,2,2)), "black", "dark grey"), bty='n', cex=1)
	at=c(log10(min(y, na.rm=T)), seq(-8, 3, by=2)); axis(2, at=10^at, labels=c("< 0", parse(text=paste("10^", at[2:length(at)], sep=""))), las=1)
	at=c(0,2,4,6); axis(1, at=10^at, labels=parse(text=paste("10^", at, sep="")))
	regrcd(hab3 = fitcdtrim$Hab3, x= fitcdtrim$ReMean, y= y, cols=cols, newx=newx, plot=T, xlog=T, zero=zero)

	dev.off()



	# K1 vs. T, D, A, hv, H1, h, x, H1/h, and hv/h
par(mfrow=c(2,5))
i = is.finite(fitcd$K1)
habs = sort(as.character(unique(fitcd$Hab2[i])))
cols = c('black', 'red', 'blue', 'green', 'orange')
ylims = range(fitcd$K1[i], na.rm=T)
plot(fitcd$T[i][fitcd$Hab2[i]==habs[1]], fitcd$K1[i][fitcd$Hab2[i]==habs[1]], col=cols[j], pch=16, xlab="T", ylab="K1", ylim=ylims, xlim=range(fitcd$T[i], na.rm=T))
	for(j in 2:length(habs)) points(fitcd$T[i][fitcd$Hab2[i]==habs[j]],fitcd$K1[i][fitcd$Hab2[i]==habs[j]], col=cols[j], pch=16)
plot(fitcd$D[i][fitcd$Hab2[i]==habs[1]],fitcd$K1[i][fitcd$Hab2[i]==habs[1]], col=cols[j], pch=16, xlab="D", ylab="K1", ylim=ylims, xlim=range(fitcd$D[i], na.rm=T))
	for(j in 2:length(habs)) points(fitcd$D[i][fitcd$Hab2[i]==habs[j]],fitcd$K1[i][fitcd$Hab2[i]==habs[j]], col=cols[j], pch=16)
plot(fitcd$A[i][fitcd$Hab2[i]==habs[1]],fitcd$K1[i][fitcd$Hab2[i]==habs[1]], col=cols[j], pch=16, xlab="A", ylab="K1", ylim=ylims, xlim=range(fitcd$A[i], na.rm=T))
	for(j in 2:length(habs)) points(fitcd$A[i][fitcd$Hab2[i]==habs[j]],fitcd$K1[i][fitcd$Hab2[i]==habs[j]], col=cols[j], pch=16)
plot(fitcd$hv[i][fitcd$Hab2[i]==habs[1]],fitcd$K1[i][fitcd$Hab2[i]==habs[1]], col=cols[j], pch=16, xlab="hv", ylab="K1", ylim=ylims, xlim=range(fitcd$hv[i], na.rm=T))
	for(j in 2:length(habs)) points(fitcd$hv[i][fitcd$Hab2[i]==habs[j]],fitcd$K1[i][fitcd$Hab2[i]==habs[j]], col=cols[j], pch=16)
plot(fitcd$H1[i][fitcd$Hab2[i]==habs[1]],fitcd$K1[i][fitcd$Hab2[i]==habs[1]], col=cols[j], pch=16, xlab="H1", ylab="K1", ylim=ylims, xlim=range(fitcd$H1[i], na.rm=T))
	for(j in 2:length(habs)) points(fitcd$H1[i][fitcd$Hab2[i]==habs[j]],fitcd$K1[i][fitcd$Hab2[i]==habs[j]], col=cols[j], pch=16)
plot(fitcd$have[i][fitcd$Hab2[i]==habs[1]],fitcd$K1[i][fitcd$Hab2[i]==habs[1]], col=cols[j], pch=16, xlab="h", ylab="K1", ylim=ylims, xlim=range(fitcd$have[i], na.rm=T))
	for(j in 2:length(habs)) points(fitcd$have[i][fitcd$Hab2[i]==habs[j]],fitcd$K1[i][fitcd$Hab2[i]==habs[j]], col=cols[j], pch=16)
plot(fitcd$xmax[i][fitcd$Hab2[i]==habs[1]],fitcd$K1[i][fitcd$Hab2[i]==habs[1]], col=cols[j], pch=16, xlab="x", ylab="K1", ylim=ylims, xlim=range(fitcd$xmax[i], na.rm=T))
	for(j in 2:length(habs)) points(fitcd$xmax[i][fitcd$Hab2[i]==habs[j]],fitcd$K1[i][fitcd$Hab2[i]==habs[j]], col=cols[j], pch=16)
plot(fitcd$H1[i][fitcd$Hab2[i]==habs[1]]/fitcd$have[i][fitcd$Hab2[i]==habs[1]],fitcd$K1[i][fitcd$Hab2[i]==habs[1]], col=cols[j], pch=16, xlab="H1/h", ylab="K1", ylim=ylims, xlim=range(fitcd$H1[i]/fitcd$have[i], na.rm=T))
	for(j in 2:length(habs)) points(fitcd$H1[i][fitcd$Hab2[i]==habs[j]]/fitcd$have[i][fitcd$Hab2[i]==habs[j]],fitcd$K1[i][fitcd$Hab2[i]==habs[j]], col=cols[j], pch=16)
plot(fitcd$hv[i][fitcd$Hab2[i]==habs[1]]/fitcd$have[i][fitcd$Hab2[i]==habs[1]],fitcd$K1[i][fitcd$Hab2[i]==habs[1]], col=cols[j], pch=16, xlab="hv/h", ylab="K1", ylim=ylims, xlim=range(fitcd$hv[i]/fitcd$have[i], na.rm=T))
	for(j in 2:length(habs)) points(fitcd$hv[i][fitcd$Hab2[i]==habs[j]]/fitcd$have[i][fitcd$Hab2[i]==habs[j]],fitcd$K1[i][fitcd$Hab2[i]==habs[j]], col=cols[j], pch=16)
legend('topright', legend=habs, col=cols, pch=16)

	
	# K1 Lab vs. Field boxplot
	quartz(width=4, height=4)
	# jpeg(width=4, height=4, res=300, units='in', file=paste('Figures/k lab vs. field_', Sys.Date(), '.jpg', sep=''))
	boxplot(I(log10(k)) ~ Field, data=fitcd, names=c('Lab', 'Field'), subset = fitcd$k>0, yaxt='n', ylab='k')
	axis(2, at= -4:1, labels = parse(text = paste("10^", -4:1, sep = "")))

	t.test(k ~ Field, data=fitcd, subset=fitcd$k>0) # p = 3e-6
	table(fitcd$Field[fitcd$k>0])

	dev.off()



### Cd vs. Re colored by habitat type
	# piles data
	piles = read.csv("analysis/Fit_Cd 101008/Greg Piles data/DDp227Waves_CdRe.txt", skip=1, header=F, col.names=c('Re', 'Cd'))
	pilescem = read.csv("analysis/Fit_Cd 101008/Greg Piles data/CEMp571Waves_Smooth.txt", skip=1, header=F, col.names=c('Re', 'Cd'))
	
	cols = c('royalblue4', 'red4', 'turquoise3', 'seagreen4', 'sienna1', 'slategray', 'steelblue1', 'wheat3')
	cols = c('black', 'red4', 'seagreen4', 'sienna1', 'steelblue1', 'slategray', 'black', 'black') # use this 10/15/2010 and later
	x = exp(seq(1,log(max(fitcd$ReMean, na.rm=T)), length.out=500))
	xlims = c(1, max(c(fitcd$ReMean, piles$Re, pilescem$Re), na.rm=T))

	y = fitcd$Cd
	y[y<=0] = min(y[y>0], na.rm=T)/2

	# quartz(width=5, height=5)
	# pdf(file=paste("analysis/Fit_Cd 101008/Figures/CdvRe_", Sys.Date(), ".pdf", sep=""), width=6, height=6)
	# pdf(file=paste("analysis/Fit_Cd 101008/Figures/CdvRe_wPiles", Sys.Date(), ".pdf", sep=""), width=6, height=6)
	# pdf(file=paste("analysis/Fit_Cd 101008/Figures/CdvRe_wRegr", Sys.Date(), ".pdf", sep=""), width=6, height=6)
	plot(fitcd$ReMean, y, log='xy', pch=16, col=cols[fitcd$Hab3], xlab="Reynolds Number", ylab="Cd", xlim=xlims, 	yaxt='n', xaxt='n')

	legend('bottomleft', legend=c(as.character(unique(fitcd$Hab3[!is.na(fitcd$Cd)])), "Piles (D&D)", "Piles (CEM)"), pch=16, col = c(cols[unique(as.numeric(fitcd$Hab3[!is.na(fitcd$Cd)]))], "black", "dark grey"), bty='n', cex=0.8)
	at=c(log10(min(y, na.rm=T)), seq(-10, 3, by=2)); axis(2, at=10^at, labels=c("< 0", parse(text=paste("10^", at[2:length(at)], sep=""))), las=1)
	at=c(0,2,4,6); axis(1, at=10^at, labels=parse(text=paste("10^", at, sep="")))

	points(pilescem$Re, pilescem$Cd, col="dark grey", pch=16) # piles data from CEM
	points(piles$Re, piles$Cd, col="black", pch=16) # piles data from D&D
	lines(c(3.3e4, 7.6e5), c(1000, 1000), col="black", lwd=2, lty=3) # range of Re for piles in storm waves of 0.1-2.1 m depth (H = 0.08-1.64)

	dev.off()

	# fit linear models in log-log space (same as y = a*x^b)
	for(i in 1:length(unique(fitcd$Hab3))){
		hab = unique(fitcd$Hab3)[i]
		ind = unique(as.numeric(fitcd$Hab3))[i]
		data = fitcd[fitcd$Hab3==hab & !is.na(fitcd$Cd),c('Cd', 'ReMean')]
		if(nrow(data)>0){
			mod = lm(I(log(Cd)) ~ I(log(ReMean)), data=data)
			lines(x, exp(predict(mod, newdata=list(ReMean=x))), col=cols[ind], lty=2, lwd=2)
		}
	}
	
	# fit nonlinear models: b1+(b2/Re)^b3
	for(i in 1:length(unique(fitcd$Hab3))){
		hab = unique(fitcd$Hab3)[i]
		ind = unique(as.numeric(fitcd$Hab3))[i]
		data = fitcd[fitcd$Hab3==hab & !is.na(fitcd$Cd) & fitcd$Cd>0,c('Cd', 'ReMean')]
		if(nrow(data)>0){
#			mod = nls(Cd ~ beta1+(beta2/ReMean)^beta3, data=data, start=list(beta1=1, beta2=1, beta3=1.1), algorithm="port", lower=list(beta1=min(y, na.rm=T), beta2=-1000, beta3=-100)) # doesn't converge for seagrass
			mod = nls(Cd ~ beta1+(beta2/ReMean)^beta3, data=data, start=list(beta1=-10, beta2=1000, beta3=0.63), control=list(maxiter=10000), trace=T) 	
				# Seagrass converges from -10, 1000, 0.63 (but not 0.61)
				# Mangrove converges from -10, 1000, 0.63 (but not 0.65)
			lines(x, predict(mod, newdata=list(ReMean=x)), col=cols[ind], lty=2, lwd=5)
		}
	}


	# fit nonlinear models: b1*Re^b2
	for(i in 1:length(unique(fitcd$Hab3))){
		hab = unique(fitcd$Hab3)[i]
		ind = unique(as.numeric(fitcd$Hab3))[i]
		data = fitcd[fitcd$Hab3==hab & !is.na(fitcd$Cd) & fitcd$Cd>0,c('Cd', 'ReMean')]
		if(nrow(data)>0){
			mod = nls(Cd ~ beta1*ReMean^beta2, data=data, start=list(beta1=10, beta2=1), control=list(maxiter=10000)) 
			lines(x, predict(mod, newdata=list(ReMean=x)), col=cols[ind], lty=2, lwd=5)
		}
	}	
	
	## Playing with nonlinear models to learn their behavior
	wavexp = function(x, b1, b2, b3){ y = b1+(b2/x)^b3; return(y) }
	plot(x, wavexp(x, 0, 10000, -0.6), col='red', log='xy')
	lines(x, wavexp(x, 1e-10, -500, 2), col='red')

	wavexp2 = function(x, b1, b2){ y = b1*x^b2; return(y) }
	plot(x, wavexp2(x, 1, -10), col='red', log='xy')
	lines(x, wavexp2(x, 200, -10), col='red')
	

### Cd vs. Re colored by habitat type (including artificial)
	
	# cols = c('purple', 'red3', 'red1', 'seagreen4', 'royalblue3', 'steelblue1', 'wheat4', 'wheat3')
	# cols = c('royalblue4', 'red4', 'turquoise3', 'seagreen4', 'sienna1', 'slategray', 'steelblue1', 'wheat3')
	cols = c('black', 'red4', 'seagreen4', 'sienna1', 'steelblue1', 'slategray', 'black', 'black') # use this 10/15/2010 and later
	pchs = c(16, 16, 1, 16, 16, 1, 16, 1)
	x = exp(seq(1,log(max(fitcd$ReMean, na.rm=T)), length.out=50))
	xlims = c(1, max(c(fitcd$ReMean, piles$Re, pilescem$Re), na.rm=T))

	y = fitcd$Cd
	dropy = y<=0 & !is.na(y)
	which(dropy)
	fitcd$X.kk[dropy]
	y[dropy] = min(y[y>0], na.rm=T)/2

	# quartz(width=5, height=5)
	# pdf(file=paste("analysis/Fit_Cd 101008/Figures/CdvRe_RealvArt", Sys.Date(), ".pdf", sep=""), width=6, height=6)
	# pdf(file=paste("analysis/Fit_Cd 101008/Figures/CdvRe_RealvArt_wRegr", Sys.Date(), ".pdf", sep=""), width=6, height=6)
	plot(fitcd$ReMean, y, log='xy', pch=pchs[fitcd$Hab2], col=cols[fitcd$Hab3], xlab="Reynolds Number", ylab="Cd", xlim=xlims, yaxt='n', xaxt='n')

	points(pilescem$Re, pilescem$Cd, col="dark grey", pch=16) # piles data from CEM
	points(piles$Re, piles$Cd, col="black", pch=16) # piles data from D&D
	lines(c(3.3e4, 7.6e5), c(1000, 1000), col="black", lwd=2, lty=3) # range of Re for piles in storm waves of 0.1-2.1 m depth (H = 0.08-1.64)

	legend('bottomleft', legend=c(as.character(unique(fitcd$Hab2[!is.na(fitcd$Cd)])), "Piles (D&D)", "Piles (CEM)"), pch=c(pchs[unique(as.numeric(fitcd$Hab2[!is.na(fitcd$Cd)]))], 16, 16), col = c(cols[c(5,5,4,4,3,2)], "black", "dark grey"), bty='n', cex=0.8)
	at=c(log10(min(y, na.rm=T)), seq(-10, 3, by=2)); axis(2, at=10^at, labels=c("< 0", parse(text=paste("10^", at[2:length(at)], sep=""))), las=1)
	at=c(0,2,4,6); axis(1, at=10^at, labels=parse(text=paste("10^", at, sep="")))

	# fit linear models in log-log space (same as y = a*x^b)
	for(i in 1:length(unique(fitcd$Hab3))){
		hab = unique(fitcd$Hab3)[i]
		ind = unique(as.numeric(fitcd$Hab3))[i]
		data = fitcd[fitcd$Hab3==hab & !is.na(fitcd$Cd),c('Cd', 'ReMean')]
		if(nrow(data)>0){
			mod = lm(I(log(Cd)) ~ I(log(ReMean)), data=data)
			lines(x, exp(predict(mod, newdata=list(ReMean=x))), col=cols[ind], lty=2, lwd=2)
		}
	}

	dev.off()



### Cd vs. Re colored by study (a complicated plot, could use tweaking colors/point styles)
	require(RColorBrewer)

	cols = c(brewer.pal(8, "Accent"), brewer.pal(8, "Dark2"), brewer.pal(9, "Set1"), brewer.pal(9, "Set3"), 'red', 'green')
	pchs = c(1,1,16,2,16,1,16,2,16,1,16,1,16,1,16,1,16,1,16,1,16,1,3,16,16,2,16,3,1,1,1,2,3,4,5,6)
	x = exp(seq(1,log(max(fitcd$ReMean, na.rm=T)), length.out=50))
	xlims = c(1, max(c(fitcd$ReMean, piles$Re, pilescem$Re), na.rm=T))
	y = fitcd$Cd
	y[y<=0] = min(y[y>0], na.rm=T)/2


	# quartz(width=10, height=5)
	# pdf(file=paste("analysis/Fit_Cd 101008/Figures/CdvRe_byStudy", Sys.Date(), ".pdf", sep=""), width=10, height=5)
	layout(mat = matrix(c(1,2), nrow=1), widths= c(1,1))
	plot(fitcd$ReMean, y, log='xy', pch=pchs[fitcd$SpStudy], col=cols[fitcd$SpStudy], xlab="Reynolds Number", ylab="Cd", xlim=xlims, 	yaxt='n', xaxt='n', cex=1.3)
	at=c(log10(min(y, na.rm=T)), seq(-10, 3, by=2)); axis(2, at=10^at, labels=c("< 0", parse(text=paste("10^", at[2:length(at)], sep=""))), las=1)
	at=c(0,2,4,6); axis(1, at=10^at, labels=parse(text=paste("10^", at, sep="")))

	opar = par(mar=c(0,0,0,0))
	plot(0,0, bty="n", xaxt='n', yaxt='n', xlab='', ylab='', col='white')
	i = match(unique(fitcd$SpStudy[!is.na(fitcd$Cd)]), fitcd$SpStudy)
	legend('left', legend=paste(fitcd$Study[i], " (", fitcd$Hab2[i], ")", sep=""), col = cols[unique(as.numeric(fitcd$SpStudy[!is.na(fitcd$Cd)]))], pch=pchs[unique(as.numeric(fitcd$SpStudy[!is.na(fitcd$Cd)]))], bty='n', cex=0.8)
	par(opar)

	dev.off()


### Cd vs. Re2 colored by habitat type
	# piles data
	piles = read.csv("analysis/Fit_Cd 101008/Greg Piles data/DDp227Waves_CdRe.txt", skip=1, header=F, col.names=c('Re', 'Cd'))
	pilescem = read.csv("analysis/Fit_Cd 101008/Greg Piles data/CEMp571Waves_Smooth.txt", skip=1, header=F, col.names=c('Re', 'Cd'))
	
	cols = c('royalblue4', 'red4', 'turquoise3', 'seagreen4', 'sienna1', 'slategray', 'steelblue1', 'wheat3')
	cols = c('black', 'red4', 'seagreen4', 'sienna1', 'steelblue1', 'slategray', 'black', 'black') # use this 10/15/2010 and later
	x = exp(seq(1,log(max(fitcd$Re2Mean, na.rm=T)), length.out=50))
	xlims = c(1, max(c(fitcd$Re2Mean, piles$Re, pilescem$Re), na.rm=T))
	y = fitcd$Cd
	y[y<=0] = min(y[y>0], na.rm=T)/2

	# quartz(width=5, height=5)
	# pdf(file=paste("analysis/Fit_Cd 101008/Figures/CdvRe2_", Sys.Date(), ".pdf", sep=""), width=5, height=5)
	# pdf(file=paste("analysis/Fit_Cd 101008/Figures/CdvRe2_wPiles", Sys.Date(), ".pdf", sep=""), width=5, height=5)
	plot(fitcd$Re2Mean, y, log='xy', pch=16, col=cols[fitcd$Hab3], xlab="Reynolds Number", ylab="Cd", xlim=xlims, 	yaxt='n', xaxt='n')

	legend('bottomleft', legend=c(as.character(unique(fitcd$Hab3[!is.na(fitcd$Cd)])), "Piles (D&D)", "Piles (CEM)"), pch=16, col = c(cols[unique(as.numeric(fitcd$Hab3[!is.na(fitcd$Cd)]))], "black", "dark grey"), bty='n', cex=0.8)
	at=c(log10(min(y, na.rm=T)), seq(-10, 3, by=2)); axis(2, at=10^at, labels=c("< 0", parse(text=paste("10^", at[2:length(at)], sep=""))), las=1)
	at=c(0,2,4,6); axis(1, at=10^at, labels=parse(text=paste("10^", at, sep="")))

	points(pilescem$Re, pilescem$Cd, col="dark grey", pch=16) # piles data from CEM
	points(piles$Re, piles$Cd, col="black", pch=16) # piles data from D&D
	lines(c(3.3e4, 7.6e5), c(1000, 1000), col="black", lwd=2, lty=3) # range of Re for piles in storm waves of 0.1-2.1 m depth (H = 0.08-1.64)

	dev.off()


### Cd vs. KC
	# piles data
	pileskc = read.csv("analysis/Fit_Cd 101008/Greg Piles data/DDp235WavesCdKC.txt", skip=1, header=F, col.names=c('KC', 'Cd'))
	pileskccem = read.csv("analysis/Fit_Cd 101008/Greg Piles data/CEMp570Waves_CdKcforRe10e4.txt", skip=1, header=F, col.names=c('KC', 'Cd'))

	cols = c('royalblue4', 'red4', 'turquoise3', 'seagreen4', 'sienna1', 'slategray', 'steelblue1', 'wheat3')
	cols = c('black', 'red4', 'seagreen4', 'sienna1', 'steelblue1', 'slategray', 'black', 'black') # use this 10/15/2010 and later
	x = exp(seq(1,log(max(fitcd$KCMean, na.rm=T)), length.out=50))
	xlims = c(1, max(c(fitcd$KCMean[is.finite(fitcd$KCMean)], pileskc$KC, pileskccem$KC), na.rm=T))
	y = fitcd$Cd
	y[y<=0] = min(y[y>0], na.rm=T)/2

	# quartz(width=5, height=5)
	# pdf(file=paste("analysis/Fit_Cd 101008/Figures/CdvKC_", Sys.Date(), ".pdf", sep=""), width=6, height=6)
	# pdf(file=paste("analysis/Fit_Cd 101008/Figures/CdvKC_wPiles", Sys.Date(), ".pdf", sep=""), width=6, height=6)
	plot(fitcd$KCMean, y, log='xy', pch=16, col=cols[fitcd$Hab3], xlab="KC Number", ylab="Cd", xlim=xlims, 	yaxt='n', xaxt='n')

	legend('bottomleft', legend=c(as.character(unique(fitcd$Hab3[!is.na(fitcd$Cd)])), "Piles (D&D)", "Piles (CEM)"), pch=16, col = c(cols	[unique(as.numeric(fitcd$Hab3[!is.na(fitcd$Cd)]))], "black", "dark grey"), bty='n', cex=0.8)
	at=c(log10(min(y, na.rm=T)), seq(-10, 3, by=2)); axis(2, at=10^at, labels=c("< 0", parse(text=paste("10^", at[2:length(at)], sep=""))), las=1)
	at=c(0,2,4); axis(1, at=10^at, labels=parse(text=paste("10^", at, sep="")))

	points(pileskc$KC, pileskc$Cd, col="black", pch=16) # piles data
	points(pileskccem$KC, pileskccem$Cd, col="dark grey", pch=16) # piles data from D&D
	lines(c(8, 177), c(1000, 1000), col="black", lwd=2, lty=3) # range of KC for piles in storm waves of 0.1-2.1 m depth (H = 0.08-1.64) and T = 10
	lines(c(12, 266), c(800, 800), col="black", lwd=2, lty=3) # range of KC for piles in storm waves of 0.1-2.1 m depth (H = 0.08-1.64) and T = 15

	 dev.off()



### K1 vs. Re
	cols = c('royalblue4', 'red4', 'turqoise3', 'seagreen4', 'sienna1', 'slategray', 'steelblue1', 'wheat3')
	cols = c('black', 'red4', 'seagreen4', 'sienna1', 'steelblue1', 'slategray', 'black', 'black') # use this 10/15/2010 and later
	x = exp(seq(1,log(max(fitcd$ReMean, na.rm=T)), length.out=50))
	xlims = c(1, max(fitcd$ReMean, na.rm=T))

	# quartz(width=5, height=5)
	# pdf(file=paste("analysis/Fit_Cd 101008/Figures/K1vRe_", Sys.Date(), ".pdf", sep=""), width=5, height=5)
	plot(fitcd$ReMean, fitcd$K1, log='xy', pch=16, col=cols[fitcd$Hab3], xlab="Reynolds Number", ylab="k (exponential decay)", xlim=xlims, yaxt='n', xaxt='n')

	legend('bottomleft', legend=as.character(unique(fitcd$Hab3[!is.na(fitcd$Cd)])), pch=16, col = cols	[unique(as.numeric(fitcd$Hab3[!is.na(fitcd$Cd)]))], bty='n', cex=0.8)
	axis(2, at=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c("0.0001", "0.001", "0.01", "0.1", "1"), las=1)
	at=c(0,2,4,6); axis(1, at=10^at, labels=parse(text=paste("10^", at, sep="")))

	# dev.off()



###############
## Exponential decay illustration

# All one graph
	x = seq(0, 200, by=0.1)
	ylims = c(0,0.8)
	xlims=c(0,100)
	cols = c('royalblue4', 'red4', 'turqoise3', 'seagreen4', 'sienna1', 'slategray', 'steelblue1', 'wheat3')
	cols = c('tan', 'red4', 'seagreen4', 'sienna1', 'steelblue1', 'slategray', 'black', 'black') # use this 10/15/2010 and later
	lwds=c(2,2,1,2,1,1,1,1)

	# quartz(width=5, height=7)
	# pdf(width=5, height=7, file=paste("analysis/Fit_Cd 101008/Figures/H_vs_x_", Sys.Date(), ".pdf", sep=""))
	x = seq(0, fitcd$xmax[1], by=0.1)
	plot(x, fitcd$H1[1]*exp(-fitcd$k[1]*x), type='l', xlab='Distance (m)', ylab='Wave height (m)', col=cols[fitcd$Hab3[1]], ylim=ylims, xlim=xlims, lwd=lwds[fitcd$Hab3[1]], log='')
	for(i in 2:nrow(fitcd)){
		x = seq(0, fitcd$xmax[i], by=0.1)
		lines(x, fitcd$H1[i]*exp(-fitcd$k[i]*x), col=cols[fitcd$Hab3[i]], lwd=lwds[fitcd$Hab3[i]])
	}
	
	legend('topright', legend=as.character(unique(fitcd$Hab3[!is.na(fitcd$K1)])), lwd = lwds[unique(as.numeric(fitcd$Hab3[!is.na(fitcd$K1)]))], col = cols	[unique(as.numeric(fitcd$Hab3[!is.na(fitcd$K1)]))], bty='o', cex=0.8, bg = "white")

	dev.off()


# Four graphs
	x = seq(0, 200, by=0.1)
	ylims = c(0.01,0.8)
	xlims=c(0.1,100)
	cols = c('royalblue4', 'red4', 'turqoise3', 'seagreen4', 'sienna1', 'slategray', 'steelblue1', 'wheat3')
	cols = c('tan', 'red4', 'seagreen4', 'sienna1', 'steelblue1', 'slategray', 'black', 'black') # use this 10/15/2010 and later
	lwds=c(2,2,1,2,1,1,1,1)
	habs = c('Kelp', 'Mangrove', 'Marsh', 'Seagrass')

	# quartz(width=5, height=7)
	# pdf(width=5, height=7, file=paste("analysis/Fit_Cd 101008/Figures/H_vs_x_", Sys.Date(), ".pdf", sep=""))
	par(mfrow=c(2,2))
	x = seq(0, fitcd$xmax[1], by=0.1)
	plot(x, fitcd$H1[1]*exp(-fitcd$k[1]*x), type='l', xlab='Distance (m)', ylab='Wave height (m)', col=cols[fitcd$Hab3[1]], ylim=ylims, xlim=xlims, lwd=lwds[fitcd$Hab3[1]], log='')
	for(i in 2:nrow(fitcd)){
		x = seq(0, fitcd$xmax[i], by=0.1)
		lines(x, fitcd$H1[i]*exp(-fitcd$k[i]*x), col=cols[fitcd$Hab3[i]], lwd=lwds[fitcd$Hab3[i]])
	}
	
	legend('topright', legend=as.character(unique(fitcd$Hab3[!is.na(fitcd$K1)])), lwd = lwds[unique(as.numeric(fitcd$Hab3[!is.na(fitcd$K1)]))], col = cols	[unique(as.numeric(fitcd$Hab3[!is.na(fitcd$K1)]))], bty='o', cex=0.8, bg = "white")

	dev.off()

#############################################
######### Aggregated by study ###############

# Find mean within studies
spstdcd = aggregate(fitcd[,c('k', 'Cd', 'ReMean')], by=list(SpStudy = fitcd$SpStudy), FUN=mean, na.rm=T)
dim(spstdcd)
spstdcd = merge(spstdcd, fitcd[!duplicated(fitcd$SpStudy), c("SpStudy", "Species", "Hab2", "Hab3", "Study", "Field")])


	## Plots
	# Plot histograms by habitat of aggregated data points
hab = as.character(unique(fitcd$Hab2))
par(mfrow=c(2, 4), mai=c(0.4, 0.4, 0.3, 0.05), mgp=c(1.8, 0.6, 0))
maxCd = max(spstdcd$Cd, na.rm=T)
hist(spstdcd$Cd, breaks=seq(0, maxCd, length.out=20), main="All Habitats", col="grey", xlab="Cd")
for(i in hab){
	hist(spstdcd$Cd[spstdcd$Hab2==i], breaks=seq(0, maxCd, length.out=20), main=i, col="grey", xlab="Cd")
}

	# Plot histograms by habitat of aggregated data points (lump real and artificial together)
hab = c("Seagrass", "Marsh", "Mangrove", "Kelp")
par(mfrow=c(2, 3), mai=c(0.4, 0.4, 0.3, 0.05), mgp=c(1.8, 0.6, 0))
maxCd = max(spstdcd$Cd, na.rm=T)
hist(spstdcd$Cd, breaks=seq(0, maxCd, length.out=20), main="All Habitats", col="grey", xlab="Cd")
for(i in hab){
	hist(spstdcd$Cd[grep(i, spstdcd$Hab2)], breaks=seq(0, maxCd, length.out=20), main=i, col="grey", xlab="Cd")
}

	# k Lab vs. Field boxplot (all habitats except bare)
	quartz(width=4, height=4)
	# jpeg(width=4, height=4, res=300, units='in', file=paste('Figures/k lab vs. field_spstudies_', Sys.Date(), '.jpg', sep=''))
	boxplot(I(log10(k)) ~ Field, data=spstdcd, names=c('Lab', 'Field'), subset = spstdcd$k>0, yaxt='n', ylab='k')
	axis(2, at= -4:1, labels = parse(text = paste("10^", -4:1, sep = "")))

	t.test(k ~ Field, data=spstdcd, subset=spstdcd$k>0) # p = 0.0022
	table(spstdcd$Field[spstdcd$k>0])
	table(spstdcd$Hab2[spstdcd$k>0])
	table(paste(spstdcd$Field[spstdcd$k>0], spstdcd$Hab2[spstdcd$k>0]))

	dev.off()

	# k Lab vs. Field boxplot (marsh/seagrass/mangrove)
	j = c('Marsh', 'Marsh (Art.)', 'Seagrass', 'Seagrass (Art.)', 'Mangrove')
	quartz(width=4, height=4)
	# jpeg(width=4, height=4, res=300, units='in', file=paste('Figures/k lab vs. field_spstudies_', Sys.Date(), '.jpg', sep=''))
	boxplot(I(log10(k)) ~ Field, data=spstdcd, names=c('Lab', 'Field'), subset = spstdcd$k>0 & spstdcd$Hab2 %in% j, yaxt='n', ylab='k')
	axis(2, at= -4:1, labels = parse(text = paste("10^", -4:1, sep = "")))

	t.test(k ~ Field, data=spstdcd, subset=spstdcd$k>0 & spstdcd$Hab2 %in% j) # p = 0.0022
	table(spstdcd$Field[i<-spstdcd$k>0 & spstdcd$Hab2 %in% j]) # 10 lab, 8 field
	table(paste(spstdcd$Field[i], spstdcd$Hab2[i])) # 3 marsh, 12 seagrass, 1 mangrove

	dev.off()

	# log10(k+o) Pretty plot for Arkema PPT: Plot histograms by habitat of all data points (no Inf)
	i = is.finite(spstdcd$k) & spstdcd$Field == 1
	# i = is.finite(spstdcd$k)
	o = 0.005 # offset
	hab = c('Marsh', 'Seagrass', 'Mangrove')
	minK = min(spstdcd$k[i & spstdcd$Hab2 %in% hab], na.rm=T)
	maxK = max(spstdcd$k[i & spstdcd$Hab2 %in% hab], na.rm=T)
	quartz(width=4, height=6)
	# jpeg(width=4, height=6, units='in', res=300, file=paste('Figures/k hist_spstudy_', Sys.Date(), '.jpg', sep=''))
	par(mfrow=c(3, 1), mai=c(0.4, 0.4, 0.3, 0.05), mgp=c(1.8, 0.6, 0))
	for(j in 1:length(hab)){
		hist(log10(spstdcd$k[i & spstdcd$Hab2 %in% hab[j]]+o), breaks=seq(log10(minK+o), log10(maxK+o), length.out=10), main=hab[j], col="dark grey", xlab="K", xaxt = 'n')
		lab = c(-3, -2, -1)
		axis(1, at=log10(c(0,10^lab)+o), labels = parse(text = c("0", paste("10^", lab, sep = ""))))
	}
	dev.off()

	mod = lm(k ~ Hab2, data=spstdcd, subset=i & spstdcd$Hab2 %in% hab)
	table(spstdcd$Hab2[i & spstdcd$Hab2 %in% hab])
	summary(mod)

	# log10(Cd) Pretty plot for Arkema PPT: Plot histograms by habitat of all data points (no Inf)
	i = is.finite(spstdcd$Cd)
	hab = c('Marsh', 'Seagrass', 'Mangrove')
	minCd = min(spstdcd$Cd[i & spstdcd$Hab2 %in% hab], na.rm=T)
	maxCd = max(spstdcd$Cd[i & spstdcd$Hab2 %in% hab], na.rm=T)
	quartz(width=4, height=6)
	# jpeg(width=4, height=6, units='in', res=300, file=paste('Figures/Cd hist_spstudy_', Sys.Date(), '.jpg', sep=''))
	par(mfrow=c(3, 1), mai=c(0.4, 0.4, 0.3, 0.05), mgp=c(1.8, 0.6, 0))
	for(j in 1:length(hab)){
		hist(log10(spstdcd$Cd[i & spstdcd$Hab2 %in% hab[j]]), breaks=seq(log10(minCd), log10(maxCd), length.out=20), main=hab[j], col="dark grey", xlab="Cd", xaxt = 'n')
		axis(1, at=-4:1, labels = parse(text = paste("10^", -4:1, sep = "")))
	}

	dev.off()

	mod = lm(Cd ~ Hab2, data=spstdcd, subset=i & spstdcd$Hab2 %in% hab)
	table(spstdcd$Hab2[i & spstdcd$Hab2 %in% hab])
	summary(mod)


	### Cd vs. Re colored by habitat type (including artificial)
	cols = c('black', 'red4', 'seagreen4', 'sienna1', 'steelblue1', 'slategray', 'black', 'black') # use this 10/15/2010 and later
	pchs = c(16, 16, 1, 16, 16, 1, 16, 1)
	x = exp(seq(1,log(max(spstdcd$ReMean, na.rm=T)), length.out=50))
	xlims = c(1, max(c(spstdcd$ReMean, piles$Re, pilescem$Re), na.rm=T))
	y = spstdcd$Cd
	y[y<=0] = min(y[y>0], na.rm=T)/2

	# quartz(width=5, height=5)
	# pdf(file=paste("analysis/Fit_Cd 101008/Figures/CdvRe_aggbyStudy_RealvArt", Sys.Date(), ".pdf", sep=""), width=6, height=6)
	# pdf(file=paste("analysis/Fit_Cd 101008/Figures/CdvRe_aggbyStudy_RealvArt_wRegr", Sys.Date(), ".pdf", sep=""), width=6, height=6)
	plot(spstdcd$ReMean, y, log='xy', pch=pchs[spstdcd$Hab2], col=cols[spstdcd$Hab3], xlab="Reynolds Number", ylab="Cd", xlim=xlims, yaxt='n', xaxt='n', cex=1.2)

	legend('bottomleft', legend=c(as.character(unique(spstdcd$Hab2[!is.na(spstdcd$Cd)])), "Piles (D&D)", "Piles (CEM)"), pch=c(pchs[unique(as.numeric(spstdcd$Hab2[!is.na(spstdcd$Cd)]))], 16, 16), col = c(cols[c(5,5,4,4,3,2)], "black", "dark grey"), bty='n', cex=0.8)
	at=c(log10(min(y, na.rm=T)), seq(-6, 3, by=2)); axis(2, at=10^at, labels=c("< 0", parse(text=paste("10^", at[2:length(at)], sep=""))), las=1)
	at=c(0,2,4,6); axis(1, at=10^at, labels=parse(text=paste("10^", at, sep="")))

	points(pilescem$Re, pilescem$Cd, col="dark grey", pch=16) # piles data from CEM
	points(piles$Re, piles$Cd, col="black", pch=16) # piles data from D&D
	lines(c(3.3e4, 7.6e5), c(100, 100), col="black", lwd=2, lty=3) # range of Re for piles in storm waves of 0.1-2.1 m depth (H = 0.08-1.64)

	# fit linear models in log-log space (same as y = a*x^b)
	for(i in 1:length(unique(spstdcd$Hab3))){
		hab = unique(spstdcd$Hab3)[i]
		ind = unique(as.numeric(spstdcd$Hab3))[i]
		data = spstdcd[spstdcd$Hab3==hab & !is.na(spstdcd$Cd),c('Cd', 'ReMean')]
		if(nrow(data)>0){
			mod = lm(I(log(Cd)) ~ I(log(ReMean)), data=data)
			lines(x, exp(predict(mod, newdata=list(ReMean=x))), col=cols[ind], lty=2, lwd=2)
		}
	}

	dev.off()



######## Aggregated by species ##############

spcd = aggregate(spstdcd$x, by=list(Species = spstdcd$Species), FUN=mean, na.rm=T)
spcdsd = aggregate(spstdcd$x, by=list(Species = spstdcd$Species), FUN=sd, na.rm=T)
spcd = merge(spcd, spcdsd, by="Species")
names(spcd) = c("Species", "Cd", "Cdsd")


