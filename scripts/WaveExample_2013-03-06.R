setwd("/Users/mpinsky/Documents/Stanford/NatCap/Coastal Protection Model")
source('Analysis/Wave_Baldock_GamKnown_H 2012-03-21.R')



#################################
## Storm conditions

# 3/21/2012 weak hurricane (8 sec period, max wave 5 m)
h = c(10,3,3,3) # depths for kelp, mangrove, marsh, and seagrass
cond = data.frame(Hab = c('Kelp', 'Mangrove', 'Marsh', 'Seagrass'), H1 = pmin(rep(5,4), 0.78*h), T = c(8,8,8,8), h = h, A = c(0.01, 0.2, 0.007, 0.005), hv = c(10, 5, 1.5, 0.33), D = c(6, 1, 500, 1500), a = c(-0.58, 3.16, -1.72, -0.36), ase = c(0.76, 1.1, 0.93, 0.96), b=c(-1.45, -3.57, -1.67, -0.46), bse = c(1.02, 1.37, 1.19, 1.24), Re = rep(NA, 4), Cd = rep(NA, 4), Cdlow = rep(NA, 4), Cdhigh = rep(NA, 4), Cddumb = c(0.4, 1.5, 2.6, 2.5), H2 = rep(NA,4), H2low = rep(NA,4), H2high = rep(NA,4), H2dumb = rep(NA,4), k = c(0.015, 0.0064, 0.11, 0.092), H2k = rep(NA,4)) # the storm conditions and Cd-Re relationships (a and b). Cddumb is from average across all studies within a habitat, without adjusting for Re. k is the exponential decay coefficient.


# 3/6/2013 calm conditions (5 sec period, max wave 1 m)
cond = data.frame(Hab = c('Kelp', 'Mangrove', 'Marsh', 'Seagrass'), H1 = c(1,1,1,1), T = c(5,5,5,5), h = c(10,3,3,3), A = c(0.01, 0.2, 0.007, 0.005), hv = c(10, 5, 1.5, 0.33), D = c(6, 1, 500, 1500), a = c(-0.58, 3.16, -1.72, -0.36), ase = c(0.76, 1.1, 0.93, 0.96), b=c(-1.45, -3.57, -1.67, -0.46), bse = c(1.02, 1.37, 1.19, 1.24), Re = rep(NA, 4), Cd = rep(NA, 4), Cdlow = rep(NA, 4), Cdhigh = rep(NA, 4), Cddumb = c(0.4, 1.5, 2.6, 2.5), H2 = rep(NA,4), H2low = rep(NA,4), H2high = rep(NA,4), H2dumb = rep(NA,4), k = c(0.015, 0.0064, 0.11, 0.092), H2k = rep(NA,4)) # the storm conditions and Cd-Re relationships (a and b). Cddumb is from average across all studies within a habitat, without adjusting for Re. k is the exponential decay coefficient.

# 3/21/2012 weak hurricane (8 sec period, max wave 5 m) with shallower marsh (2m)
h = c(10,3,2,3) # depths for kelp, mangrove, marsh, and seagrass
cond = data.frame(Hab = c('Kelp', 'Mangrove', 'Marsh', 'Seagrass'), H1 = pmin(rep(5,4), 0.78*h), T = c(8,8,8,8), h = h, A = c(0.01, 0.2, 0.007, 0.005), hv = c(10, 5, 1.5, 0.33), D = c(6, 1, 500, 1500), a = c(-0.58, 3.16, -1.72, -0.36), ase = c(0.76, 1.1, 0.93, 0.96), b=c(-1.45, -3.57, -1.67, -0.46), bse = c(1.02, 1.37, 1.19, 1.24), Re = rep(NA, 4), Cd = rep(NA, 4), Cdlow = rep(NA, 4), Cdhigh = rep(NA, 4), Cddumb = c(0.38, 1.5, 2.6, 2.5), H2 = rep(NA,4), H2low = rep(NA,4), H2high = rep(NA,4), H2dumb = rep(NA,4), k = c(0.015, 0.0064, 0.11, 0.092), H2k = rep(NA,4)) # the storm conditions and Cd-Re relationships (a and b). Cddumb is from average across all studies within a habitat, without adjusting for Re. k is the exponential decay coefficient.

xmax = 200 # width of bed
dx = 0.1 # step distance for calculation
B = 1 # breaking coefficient (irrelevant since a flat bed)
adj = 3e-4 # x-axis adjust from Cd-Re equation

for(i in 1:nrow(cond)){
	print(i)
	sig = 2*pi/cond$T[i]
	k=iterativek(sig,cond$h[i])$k
	ubm=0.5*cond$H1[i]*sig*cosh(k*(cond$h[i]-cond$hv[i]))/sinh(k*cond$h[i]); #Wave orb. velo. at top of veg
	cond$Re[i] = cond$A[i]*ubm/1.17e-6; # Reynolds number
	#cond$Cd[i] = 10^(cond$a[i] + cond$b[i]*log10(cond$Re[i])) 
	cond$Cd[i] = 10^(cond$a[i] + cond$b[i]*(log10(min(10^5, cond$Re[i])*adj))) #evaluate at Re no higher than fully turbulent (Cd stays constant after Re 10^5))
	cond$Cdlow[i] = 10^(cond$a[i] - cond$ase[i] + (cond$b[i]-cond$bse[i])*(log10(min(10^5, cond$Re[i])*adj))) # use lower estimates of slope and intercept, expanded to 95%CI
	cond$Cdhigh[i] = 10^(cond$a[i] + cond$ase[i] + (cond$b[i] + cond$bse[i])*(log10(min(10^5, cond$Re[i])*adj))) # use lower estimates of slope and intercept, expanded to 95%CI
	
	## Transform wave
	H = Wave_Baldock_GamKnown_H_flat(cond$H1[i],cond$T[i],cond$Cd[i],cond$A[i],cond$D[i],cond$hv[i],B,cond$h[i],xmax,dx);
	cond$H2[i] = H[length(H)]

	H = Wave_Baldock_GamKnown_H_flat(cond$H1[i],cond$T[i],cond$Cdlow[i],cond$A[i],cond$D[i],cond$hv[i],B,cond$h[i],xmax,dx);
	cond$H2high[i] = H[length(H)]

	H = Wave_Baldock_GamKnown_H_flat(cond$H1[i],cond$T[i],cond$Cdhigh[i],cond$A[i],cond$D[i],cond$hv[i],B,cond$h[i],xmax,dx);
	cond$H2low[i] = H[length(H)]

	H = Wave_Baldock_GamKnown_H_flat(cond$H1[i],cond$T[i],cond$Cddumb[i],cond$A[i],cond$D[i],cond$hv[i],B,cond$h[i],xmax,dx);
	cond$H2dumb[i] = H[length(H)]

	cond$H2k[i] = cond$H1[i]*exp(-cond$k[i]*xmax)
}
cond$atten = signif((cond$H1-cond$H2)/cond$H1*100,3) # in percent
cond$attenlow = signif((cond$H1-cond$H2high)/cond$H1*100,3) # in percent
cond$attenhigh = signif((cond$H1-cond$H2low)/cond$H1*100,3) # in percent
cond$attendumb = signif((cond$H1-cond$H2dumb)/cond$H1*100,3) # in percent
cond$attenk = signif((cond$H1-cond$H2k)/cond$H1*100,3) # in percent

cond[,c('Hab', 'H1', 'h', 'A', 'hv', 'D', 'Re', 'Cd', 'Cdlow', 'Cdhigh')] # compare the basic set-up of each scenario

cond[,c('Hab', 'attenlow', 'atten', 'attenhigh', 'attendumb', 'attenk')] # compare the 'dumb' attenuations

cond[,c('Hab', 'Re', 'Cdlow', 'Cd', 'Cdhigh', 'H1', 'H2', 'attenlow', 'atten', 'attenhigh')]

cond[,c('Hab', 'Re', 'Cdlow', 'Cd', 'Cdhigh', 'H1', 'H2')]

(cond$attendumb-cond$atten)/cond$atten # overestimation by 'attendumb' (average Cd)
(cond$attenk-cond$atten)/cond$atten # overestimation by 'attenk' (exponential decay)
