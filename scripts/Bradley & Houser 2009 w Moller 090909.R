# Fit Moller 1999 equations to Bradley & Houser 2009 results
# Extrapolate to erosion onshore using Stockdon et al. 2006 (run-up) and Komar et al. 1999 (erosion)

# Function to find f, given all the other params (use optimise(molleroptim, interval=c(_,), H1=_, ...))
molleroptim = function(f, H1, H2, x, vk, T, h, h1, h2, K, d, g){
	# one of the inputs can be >1 for sensitivity analysis
	L = T*sqrt(g*h) # wave length (linear wave theory, shallow water approximation)
	c1 = sqrt(g*h1) # wave celerity
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
	return(abs(H2-H2pred))	
}

# Function to compute H2 by Moller's 1999 formulas
moller = function(f, H1, x, vk, T, h, h1, h2, K, d, g){
	# one of the inputs can be >1 for sensitivity analysis
	L = T*sqrt(g*h) # wave length
	c1 = sqrt(g*h1) # wave celerity
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

	#print(paste("Ks: ", Ks, "  Kv: ", Kv, "  Kf: ", Kf, "  Kp: ", Kp))

	H2 = H1*Kv*Ks*Kp*Kf		
	return(H2)	
}

#This function computes wave number based on water depth and wave period#using the linear dispersion relation and an iterative Newton-Raphson#technique (from Greg, via Katie Arkema)
func_disp = function(T,h,g){	#Define Constants	X = numeric(0)
	sigma=2*pi/T	k=sigma^2/g*(tanh(sigma^2*h/g))^(-0.5) #rough approx of k	X = c(X, k*h)	cst=sigma^2*h/g	#Initialize loop	i=1	f=cst-X[1]*tanh(X[1]);	r=-tanh(X[1])-X[1]*(1-tanh(X[1])^2)	#Loop to determine f	while(abs(f)>=1e-3){    	i=i+1    	X[i]=X[1]-f/r;    	f=cst-X[i]*tanh(X[i]);    	r=-tanh(X[i])-X[i]*(1-tanh(X[i])^2)	}	k2=X[i]/h
	return(k2)
}

# Find H0 and L0 (deep water wave height and length from T) (from K. Arkema & Greg)
# Uses Greg's equation to find wave number (k)
deepH0 = function(T, h, H1, g){
	w = 2*pi/T # w = angular frequency (radians per second)
	k = func_disp(T,h, g) # wave number (2pi/L) (1/meters)
	c=w/k
	n= 1/2*(1+(2*k*h)/sinh(2*k*h))
	Cg=c*n
	Co=g*T/(2*pi)
	Cgo=1/2*Co
	H0=H1*sqrt(Cg/Cgo) # shoal back to deep water
	return(H0)
}

deepL0 = function(T, g){
	L0 = g*T^2/(2*pi)
	return(L0)
}

# Find H1 (shallow water wave height) from H0 (deep-water height) (from K. Arkema & Greg)
shallowH1 = function(T, h, H0, g){
	w = 2*pi/T
	k = func_disp(T,h,g)
	c = w/k
	n= 1/2*(1+(2*k*h)/sinh(2*k*h))
	Cg=c*n
	Co=g*T/(2*pi)
	Cgo=1/2*Co
	H1=H0/sqrt(Cg/Cgo) # shoal back to shallow water
	return(H1)
}

# Miche's formula for breaking waves, from Moller et al 1999 (Eq. 8)
breaking = function(H1, T, h, g){
	L = T*sqrt(g*h) # wave length (linear wave theory, shallow water approximation) SHOULD THIS BE GREG'S LINEAR DISPERSION CALC?
	Hb = L*0.142*tanh(2*pi*h/L) # breaking wave height
	if(length(Hb)==1){
		Hb = rep(Hb, length(H1))	
	}
	i = H1 > Hb
	H1[i] = Hb[i]	
	return(H1)
}

# Stockdon et al 2006 Eq. 19 for runup 2% exceedance value
stockdon = function(H0, L0, Bf){
	r2 = 1.1*(0.35*Bf*sqrt(H0*L0) + sqrt(H0*L0*(0.563*Bf^2+0.004))/2)
	return(r2)
}

# Modified Stockdon et al 2006 Eq. 19 for runup 2% exceedance value
# Uses deep-water wave height (H0) for wave set-up <n>
# Uses attenuated deep-water wave height (ratio*H0) for swash (S)
# Ratio is calculated as H2/H1, where H1 and H2 are wave height at 
#   offshore and inshore edge of biogenic habitat (respectively)
stockdon2 = function(H0, L0, Bf, ratio){
	r2 = 1.1*(0.35*Bf*sqrt(H0*L0) + sqrt(ratio*H0*L0*(0.563*Bf^2+0.004))/2)
	return(r2)
}

# Komar et al. 1999 Eq. 3 for maximum dune erosion
komar = function(r2, tide, beta, deltaBL, Ej){
	WL = tide + r2 # Total Water Level (tidal elevation + storm-wave runup)
	demax = ((WL-Ej)+deltaBL)/tan(beta)
	demax[demax<0] = 0
	return(demax)
}


#######################
## Set up parameters
#######################

# Moller 1999 params
#	H1 = 0.24 # wave height (m)
#	x = 180 # distance (m)
#	vk = 1.35e-6 # kinematic viscosity
#	T = 2.8 # wave period (sec)
#	h = 1.2 # water depth (m)
#	h1 = 1.3 # depth at outer edge
#	h2 = 0.8 # depth at inner edge
#	K = 1e-10 # specific permeability
#	d = 1 # sediment depth layer (m)
#	f = 0.2401 # friction factor marsh
#	fsand = 0.01 # friction factor sand
#	g = 9.81 # gravitational constant

# Bradley & Houser 2009 params
	g = 9.81 # gravity
	H1 = 0.095 # incident wave height (m)
	H2 = 0.068 # expect 40% attenuation
	x = 39 # distance (m)
	h = 1 # water depth (m)
	h1 = h # depth at offshore end of bed
	h2 = h # depth at inshore end of bed
	#D = 1*10^-2 # width of Thalassia testudinum (from Tussenbroek 1998 Aquatic Botany)
	#b = 1/sqrt(1100) # spacing of plants (m-1) (B&H report 1100 m-1, but this seems unreasonable, more likely m-2.
	T = 1/0.75 # wave period (sec) (peak period in Fig. 4 in Bradley & Houser 2009)
	#s =  0.225 # height of veg above bottom (m)
	#L = T*sqrt(g*h) # wave length (linear wave theory, shallow water)
	#k =  2*pi/L # wave number = 2pi/wavelength
	#w = 2*pi/T # radian wave frequency
	vk = 1.35e-6 # kinematic viscosity of water (m2/s)
	K = 1e-10 # specific permeability (from Moller 1999)
	d = 1 # sediment depth layer (m) (from Moller 1999)
	tide = 0
	Bf = 1.7*pi/180 # foreshore beach slope (radians) (value from Komar et al 1999)
	beta = Bf
	deltaBL = 0 # beach-level change (value from Komar example pg. 47)
	Ej = 3.2 # elevation of dune toe junction above mean water (value from Komar example pg. 47)

	# solve for f
	fsand = 0.01 # friction factor sand
	#f = optimize(molleroptim, interval = c(0,1), H1=H1, H2=H2, x=x, vk = vk, T=T, h=h, h1=h1, h2=h2, K=K, d=d, g=g, tol=1e-9)$minimum
	f = 0.47

	# a double check: should be H1*0.7
	H1*0.7
	H2 = moller(f=f, H1=H1, x=x, vk = vk, T=T, h=h, h1=h, h2=h, K=K, d=d, g=g)
	H2

# special Bradley & Houser params for large waves (e.g. Komar et al 1999)
	T = 17 # wave period (sec) (from Stockdon et al 2006, Duck82) (also Komar et al. 1999)
	#H0 = 7.7 # wave height (m) (from Stockdon et al 2006 (Duck82), similar to 10-yr extreme waves in Komar et al 1999)
	H0 = 7.8 # wave height (m)  10-yr extreme waves in Komar et al 1999)
	tide = 2.3 # 10-yr tides from Komar et al 1999)
	h = 3.3 # water depth with the tide (assumes 1 at mid-tide)
	h1 = 3.3
	h2 = 3.3
	H1 = H0

	H1 = shallowH1(T, h, H0, g)

########################################################
## Extrapolate to maximum dune erosion (point estimate)
########################################################

H1 = breaking(H1, T, h, g) # check to see if H1 breaks

H2 = breaking(moller(f=f, H1=H1, x=x, vk = vk, T=T, h=h, h1=h, h2=h, K=K, d=d, g=g), T, h, g) # attenuate H2 and again check for breaking
H2sand = breaking(moller(f=fsand, H1=H1, x=x, vk = vk, T=T, h=h, h1=h, h2=h, K=K, d=d, g=g), T, h, g) # same over sand

# Option #1: by reverse-shoaling broken waves
	H0prime = deepH0(T, h, H2, g)
	H0primesand = deepH0(T, h, H2sand, g)
	L0prime = deepL0(T, g)
	r2 = stockdon(H0prime, L0prime, Bf)
	r2sand = stockdon(H0primesand, L0prime, Bf)
	#r2
	#r2sand
	dmax = komar(r2, tide, beta, deltaBL, Ej)
	dmaxsand = komar(r2sand, tide, beta, deltaBL, Ej)
	dmax
	dmaxsand

# Option #2: by attenuating H0
	ratio = H2/H1
	ratiosand = H2sand/H1
	H0prime = ratio*H0
	H0primesand = ratiosand*H0
	L0prime = deepL0(T, g)
	r2 = stockdon(H0prime, L0prime, Bf)
	r2sand = stockdon(H0primesand, L0prime, Bf)
	r2
	r2sand
	dmax = komar(r2, tide, beta, deltaBL, Ej)
	dmaxsand = komar(r2sand, tide, beta, deltaBL, Ej)
	dmax
	dmaxsand

# Option #3: by attenuating H0 for swash only
	ratio = H2/H1
	ratiosand = H2sand/H1
	L0prime = deepL0(T, g)
	r2 = stockdon2(H0, L0prime, Bf, ratio)
	r2sand = stockdon2(H0, L0prime, Bf, ratiosand)
	r2
	r2sand
	dmax = komar(r2, tide, beta, deltaBL, Ej)
	dmaxsand = komar(r2sand, tide, beta, deltaBL, Ej)
	dmax
	dmaxsand

	


############################
## Sensitivity analysis
############################

H1sens = seq(0.1, 3.3, by=0.1)
xsens = seq(1, 300, by = 10)
Tsens = seq(1, 20, by = 1)
hsens = seq(0.31, 3, by=0.1)
fsens = seq(0.01, 1, by=0.01)
Bfsens = seq(.5, 30, by=.5)
betasens = seq(.5, 30, by=.5)


# special Bradley & Houser params for large waves (e.g. Komar et al 1999)
	H0sens = seq(0.1, 10, length.out = 50)
	xsens = seq(1, 300, by = 10)
	Tsens = seq(1, 25, by = 1)
	hsens = seq(0.31, 5, by=0.1)

	H1sens = H0sens

	# shoaling
	H1sens = shallowH1(T, h, H0sens, g)


## Sensitivity

# Calc wave attenuation
	H2height =moller(f, H1sens, x, vk,  T, h, h1, h2, K, d, g)
	H2dist = moller(f, H1, xsens, vk,  T, h, h1, h2, K, d, g)
	H2period = moller(f, H1, x, vk, Tsens, h, h1, h2, K, d, g)
	H2depth = moller(f, H1, x, vk, T, hsens, hsens, hsens, K, d, g)
	H2friction = moller(fsens, H1, x, vk, T, h, h1, h2, K, d, g)

	# Calcs w/ sand
	H2height_s = moller(fsand, H1sens, x, vk,  T, h, h1, h2, K, d, g)
	H2dist_s = moller(fsand, H1, xsens, vk,  T, h, h1, h2, K, d, g)
	H2period_s = moller(fsand, H1, x, vk, Tsens, h, h1, h2, K, d, g)
	H2depth_s = moller(fsand, H1, x, vk, T, hsens, hsens, hsens, K, d, g)

# with breaking
	H1 = breaking(shallowH1(T, g, H0, g), T, h, g)
	H1sens = breaking(shallowH1(T, g, H0sens, g), T, h, g)

	H2height = breaking(moller(f, H1sens, x, vk,  T, h, h1, h2, K, d, g), T, h, g)
	H2dist = breaking(moller(f, H1, xsens, vk,  T, h, h1, h2, K, d, g), T, h, g)
	H2period = breaking(moller(f, H1, x, vk, Tsens, h, h1, h2, K, d, g), T, h, g)
	H2depth = breaking(moller(f, H1, x, vk, T, hsens, hsens, hsens, K, d, g), T, hsens, g)
	H2friction = breaking(moller(fsens, H1, x, vk, T, h, h1, h2, K, d, g), T, h, g)

	# Calcs w/ sand
	H2height_s = breaking(moller(fsand, H1sens, x, vk,  T, h, h1, h2, K, d, g), T, h, g)
	H2dist_s = breaking(moller(fsand, H1, xsens, vk,  T, h, h1, h2, K, d, g), T, h, g)
	H2period_s = breaking(moller(fsand, H1, x, vk, Tsens, h, h1, h2, K, d, g), T, h, g)
	H2depth_s = breaking(moller(fsand, H1, x, vk, T, hsens, hsens, hsens, K, d, g), T, hsens, g)


# Plot sensitivity (% attenuation) w/ seagrass or sand
quartz(width=10, height=7)
ylim = c(-10,100)
par(mfrow=c(2,3), cex=1, mai=c(0.8, 0.6, 0.3, 0.3), mgp = c(2, 0.5, 0))
plot(H0sens, (H1sens-H2height)/H1sens*100, type="l", main="Incident wave height", ylim=ylim, ylab="Wave attenuation (%)", xlab="Incident wave height (m)")
	abline(v = H0, col="blue", lty=3)
	lines(H0sens, (H1sens-H2height_s)/H1sens*100, col="red")
plot(hsens, (H1-H2depth)/H1*100, type="l", main="Depth", ylim=ylim, ylab="Wave attenuation (%)", xlab="Depth (m)")
	abline(v = h, col="blue", lty=3)
	lines(hsens, (H1-H2depth_s)/H1*100, col="red")
plot(-1,-1, bty="n", col="white", xaxt="n", yaxt="n", xlab="", ylab="")
	legend("topright", legend=c(paste("Seagrass (f=", round(f, 3), ")", sep=""), "Sand (f=0.01)"), lty=1, col=c("black", "red"), cex=0.8, bty="n", title="B&H2009 w/ Moller1999 eq")
plot(xsens, (H1-H2dist)/H1*100, type="l", main="Distance", ylim=ylim, ylab="Wave attenuation (%)", xlab="Distance (m)")
	abline(v = x, col="blue", lty=3)
	lines(xsens, (H1-H2dist_s)/H1*100, col="red")
plot(Tsens, (H1-H2period)/H1*100, type="l", main="Period", ylim=ylim, ylab="Wave attenuation (%)", xlab="Period (sec)")
	abline(v = T, col="blue", lty=3)
	lines(Tsens, (H1-H2period_s)/H1*100, col="red")
plot(fsens, (H1-H2friction)/H1*100, type="l", main="Friction factor", ylim=ylim, ylab="Wave attenuation (%)", xlab="Friction factor")
	abline(v = f, col="blue", lty=3)

# Plot sensitivity: w/ attenuated wave height
ylim = c(0,0.15)
par(mfrow=c(2,3))
plot(H1sens, H2height, type="l", main="Incident wave height (m)", ylim=c(0,3), ylab="Attenuated wave height (m)")
	abline(a = 0, b = 1, lty=2)
	abline(v = H1, col="blue", lty=3)
	lines(H1sens, H2height_s, col="red")
plot(hsens, H2depth, type="l", main="Depth", ylim=ylim, ylab="Attenuated wave height (m)")
	abline(h = H1, lty=2)
	abline(v = h, col="blue", lty=3)
	lines(hsens, H2depth_s, col="red")
plot(-1,-1, bty="n", col="white", xaxt="n", yaxt="n", xlab="", ylab="")
	legend("topright", legend=c(paste("Seagrass (f=", round(f, 3), ")", sep=""), "Sand (f=0.01)"), lty=1, col=c("black", "red"), cex=1, bty="n", title="Bradley & Houser 2009 w/ Moller 1999 equations")
plot(xsens, H2dist, type="l", main="Distance", ylim=ylim, ylab="Attenuated wave height (m)")
	abline(h = H1, lty=2)
	abline(v = x, col="blue", lty=3)
	lines(xsens, H2dist_s, col="red")
plot(Tsens, H2period, type="l", main="Period (sec)", ylim=ylim, ylab="Attenuated wave height (m)")
	abline(h = H1, lty=2)
	abline(v = T, col="blue", lty=3)
	lines(Tsens, H2period_s, col="red")
plot(fsens, H2friction, type="l", main="Friction factor", ylim=ylim, ylab="Attenuated wave height (m)")
	abline(h = H1, lty=2)
	abline(v = f, col="blue", lty=3)


# Runup
r2height = stockdon(deepH0(T, h, H2height, g), deepL0(T, g), Bf)
r2dist = stockdon(deepH0(T, h, H2dist, g), deepL0(T, g), Bf)
r2period = stockdon(deepH0(T, h, H2period, g), deepL0(T, g), Bf)
r2depth = stockdon(deepH0(T, h, H2depth, g), deepL0(T, g), Bf)
r2friction = stockdon(deepH0(T, h, H2friction, g), deepL0(T, g), Bf)
r2slope = stockdon(deepH0(T, h, H2, g), deepL0(T, g), Bfsens*pi/180)

# Calcs w/ sand
r2height_s = stockdon(deepH0(T, h, H2height_s, g), deepL0(T, g), Bf)
r2dist_s = stockdon(deepH0(T, h, H2dist_s, g), deepL0(T, g), Bf)
r2period_s = stockdon(deepH0(T, h, H2period_s, g), deepL0(T, g), Bf)
r2depth_s = stockdon(deepH0(T, h, H2depth_s, g), deepL0(T, g), Bf)
r2slope_s = stockdon(deepH0(T, h, H2sand, g), deepL0(T, g), Bfsens*pi/180)


# Plot sensitivity of wave run-up
quartz(width=10, height=7)
ylim = c(0,2)
par(mfrow=c(2,3), cex=1, mai=c(0.8, 0.6, 0.3, 0.3), mgp = c(2, 0.5, 0))
plot(H0sens, r2height, type="l", main="Incident wave height", ylim=ylim, ylab="Run-up (m)", xlab="Incident wave height (m)")
	abline(v = H0, col="blue", lty=3)
	lines(H0sens, r2height_s, col="red")
plot(hsens, r2depth, type="l", main="Depth", ylim=ylim, ylab="Run-up (m)", xlab="Depth (m)")
	abline(v = h, col="blue", lty=3)
	lines(hsens, r2depth_s, col="red")
#plot(-1,-1, bty="n", col="white", xaxt="n", yaxt="n", xlab="", ylab="")
#	legend("topright", legend=c(paste("Seagrass (f=", round(f, 3), ")", sep=""), "Sand (f=0.01)"), lty=1, col=c("black", "red"), cex=0.8, bty="n", title="B&H2009 w/ Moller 1999 eq")
plot(xsens, r2dist, type="l", main="Distance", ylim=ylim, ylab="Run-up (m)", xlab="Distance (m)")
	abline(v = x, col="blue", lty=3)
	lines(xsens, r2dist_s, col="red")
plot(Tsens, r2period, type="l", main="Period", ylim=ylim, ylab="Run-up (m)", xlab="Period (sec)")
	abline(v = T, col="blue", lty=3)
	lines(Tsens, r2period_s, col="red")
plot(fsens, r2friction, type="l", main="Friction factor", ylim=ylim, ylab="Run-up (m)", xlab="Friction factor")
	abline(v = f, col="blue", lty=3)
plot(Bfsens, r2slope, type="l", main="Slope", ylim=c(0,10), ylab="Run-up (m)", xlab="Foreshore slope (°)")
	lines(Bfsens, r2slope_s, col="red")


# Erosion
dmaxheight = komar(r2height, tide, beta, deltaBL, Ej)
dmaxdist = komar(r2dist, tide, beta, deltaBL, Ej)
dmaxperiod = komar(r2period, tide, beta, deltaBL, Ej)
dmaxdepth = komar(r2depth, tide, beta, deltaBL, Ej)
dmaxfriction = komar(r2friction, tide, beta, deltaBL, Ej)
dmaxslope = komar(r2slope, tide, betasens*pi/180, deltaBL, Ej)

# Calcs w/ sand
dmaxheight_s = komar(r2height_s, tide, beta, deltaBL, Ej)
dmaxdist_s = komar(r2dist_s, tide, beta, deltaBL, Ej)
dmaxperiod_s = komar(r2period_s, tide, beta, deltaBL, Ej)
dmaxdepth_s = komar(r2depth_s, tide, beta, deltaBL, Ej)
dmaxslope_s = komar(r2slope_s, tide, betasens*pi/180, deltaBL, Ej)


# Plot sensitivity of dune erosion maximum
quartz(width=10, height=7)
ylim = c(0,25)
par(mfrow=c(2,3), cex=1, mai=c(0.8, 0.6, 0.3, 0.3), mgp = c(2, 0.5, 0))
plot(H0sens, dmaxheight, type="l", main="Incident wave height", ylim=ylim, ylab="Maximum erosion (m)", xlab="Incident wave height (m)")
	abline(v = H0, col="blue", lty=3)
	lines(H0sens, dmaxheight_s, col="red")
plot(hsens, dmaxdepth, type="l", main="Depth", ylim=ylim, ylab="Maximum erosion (m)", xlab="Depth (m)")
	abline(v = h, col="blue", lty=3)
	lines(hsens, dmaxdepth_s, col="red")
#plot(-1,-1, bty="n", col="white", xaxt="n", yaxt="n", xlab="", ylab="")
#	legend("topright", legend=c(paste("Seagrass (f=", round(f, 3), ")", sep=""), "Sand (f=0.01)"), lty=1, col=c("black", "red"), cex=0.8, bty="n", title="B&H2009 w/ Moller 1999 eq")
plot(xsens, dmaxdist, type="l", main="Distance", ylim=ylim, ylab="Maximum erosion (m)", xlab="Distance (m)")
	abline(v = x, col="blue", lty=3)
	lines(xsens, dmaxdist_s, col="red")
plot(Tsens, dmaxperiod, type="l", main="Period", ylim=ylim, ylab="Maximum erosion (m)", xlab="Period (sec)")
	abline(v = T, col="blue", lty=3)
	lines(Tsens, dmaxperiod_s, col="red")
plot(fsens, dmaxfriction, type="l", main="Friction factor", ylim=ylim, ylab="Maximum erosion (m)", xlab="Friction factor")
	abline(v = f, col="blue", lty=3)
plot(Bfsens, dmaxslope, type="l", main="Slope", ylim=ylim, ylab="Run-up (m)", xlab="Foreshore slope (°)")
	lines(Bfsens, dmaxslope_s, col="red")
	
	
######################################	
# Try longer-period waves

T = 17 # wave period (sec) (from Stockdon et al 2006, Duck82)
H1 = 7.7 # wave height (m) (from Stockdon et al 2006, Duck82)
tide = 2.3 # 10-yr tides from Komar et al 1999)

H2 = moller(f=f, H1=H1, x=x, vk = vk, T=T, h=h, h1=h, h2=h, K=K, d=d, g=g)
H2sand = moller(f=fsand, H1=H1, x=x, vk = vk, T=T, h=h, h1=h, h2=h, K=K, d=d, g=g)
r2 = stockdon(T,h,Bf, H2, g)
r2sand = stockdon(T,h,Bf, H2sand, g)
r2sand 
r2
demax = komar(r2, tide, beta, deltaBL, Ej)
demaxsand = komar(r2sand, tide, beta, deltaBL, Ej)
demax
demaxsand



# graph of a wave
quartz(width=12, height=5)
x = seq(0, 20*pi, length.out=300)
y = sin(x)*c(rep(1,100), exp(-0.05*x[1:200]))
plot(x,y, type="l", lwd=2, col="blue")

quartz(width=12, height=5)
x = seq(0, 20*pi, length.out=300)
y = sin(x)
plot(x,y, type="l", lwd=2, col="blue")