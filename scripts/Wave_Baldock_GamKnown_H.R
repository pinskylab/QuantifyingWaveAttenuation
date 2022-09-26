# May not work if bed is sloping (6/24/2011)

# A wrapper function to set up Wave_Baldock_GamKnown_H with a flat bed
Wave_Baldock_GamKnown_H_flat = function(H1,T,Cd,bv,N,hv,B,hd,xmax,dx){
	# Basic set-up
	x = seq(0, xmax, by=dx); # x-axis (m)
	h = rep(hd, length(x)); # depth
	hd = c(hd, hd)
	xd = c(0,xmax)
	N = rep(N, length(x)); bv = rep(bv, length(x)); hv = rep(hv, length(x)); # make veg params into 	vectors
	lx = length(x); L = length(xd);

	sig = 2*pi/T; fp = 1/T;

	H = Wave_Baldock_GamKnown_H(H1,fp,Cd,bv,N,hv,B,h,hd,x,xd)
	return(H)
}


# The function from Greg, modified to return H and to calculate flat/sloping inside the function
Wave_Baldock_GamKnown_H = function(H1,fp,Cd,bv,N,hv,B,h,hd,x,xd){
	require(VGAM)
	
	rho = 1024
	g = 9.81
	Cf = 0.003
	up=numeric(0)

	# ID Location in depth vector where bed is flat and sloping
	#segs gives loc where sloping bed starts and ends (row vector)
	#segf gives loc where flat bed starts and ends (row vector)
#	Dif=diff(hd); flat=which(Dif==0);flat=union(flat,flat+1);xflat=numeric(0); #ID where bed is flat is flagged by 1's
#	if(length(flat)==0){
#		flat=0;segf=0;segs=c(1, L);Lf=0; #Bed has no flat portions
#	} else {
#		temp=flat*0;Lf=length(flat);
#		for(ii in 1:(Lf-1)){
#			if(flat[ii+1] !=flat[ii+1]) temp[ii]=1
#		}
#		segf=which(temp>0);  #Locate where brk from 1 to 0 exists in the flat vector
#		if(length(segf)==0){
#			segf=c(1, Lf)
#		} else {
#			a=c(1, segf+1);b=c(segf, Lf);
#			segf=cbind(a, b); segf=flat[segf]; #segf gives loc where flat bed starts and ends (row vector)
#		}
#		if(mean(segf)>0) xflat=xd[segf] # get actual x-values for the beginning and end of the flat portion
#	
#		temp=setdiff(1:L,flat);
#		slopy=union(temp,temp+1); slopy=union(slopy,temp-1);
#		slopy[slopy<1]=NA;slopy[slopy>L]=NA; #slopy gives loc of sloping bed
#		if(length(slopy)>0){
#			temp=slopy*0;Ls=length(slopy);
#			for(ii in 1:(Ls-1)){
#				if(slopy[ii+1] !=slopy[ii]+1) temp[ii]=1
#			} 
#			segs=which(temp>0);  #Locate where brk from 1 to 0 exists
#			a=c(1, segs+1);b=c(segs, Ls);
#			segs=c(a, b); segs=slopy[segs];#segs gives loc where sloping bed starts and ends (row vector)
#		} else {
#			segs=0;
#		}
#	}
#	
#	temp=numeric(0);up=numeric(0);
#	
#	#Determine location of flat portion
#	if(any(segf !=0)){
#		for(ss in 1:dim(segf)[1]){
#			temp=c(temp, seq(xd[segf[ss,1]], xd[segf[ss,2]], by=dx))
#		}
#		xflat2=temp;
#		xflat2=round(xflat2*10)/10;x=round(x*10)/10;
#	} else {
#		xflat2=1;
#	}

	# Only works if whole bed is flat for now
	if(all(hd[1] == hd)){
		flat = c(1,2)
		xflat = xd
		xflat2 = x
	} else {
		warning('Bed was not flat!')
		stop()
	}

	
	# Wave Model
	sig=2*pi*fp;
	dx=x[2]-x[1];lx=length(x); L=length(hd);
	H=x*0;Ef=H;Db=H;Df=Db;Dv=Db;Cg=H;
	
	H[1]=H1;Ew=.125*rho*g*H[1]^2; 
	temp=iterativek(sig,hd[1]);
	Cg[1] = temp$Cg
	Ef[1]=Ew*Cg[1]; #Energy flux
	
	Lo=g/(2*pi*fp^2);Co=g/sig;Cgo=Co/2;#Deep water phase speed
	hi=h[1];ki=iterativek(sig,hi)$k;Li=2*pi/ki;#Wave # @ 1st grid pt
	ni=.5*(1+(2*ki*hi/sinh(2*ki*hi))); #to compute Cg @ 1st grid pt
	Ci=Li*fp;  Cgi=Ci*ni; #Group velocity @ 1st grid pt
	Ho=H[1]*sqrt(Cgi/Cgo);#We are in intermediate water. Assume no brkg occured
	So=Ho/Lo; #Transform significant wave height to rms wave height
	
	for(xx in 1:(lx-1)){ #Transform waves over sloping bottom
		Ef[xx]=0.125*rho*g*(H[xx]^2)*Cg[xx]; #Ef at (xx)      
		Ef[xx+1]=Ef[xx]-dx*(Db[xx]+Df[xx]+Dv[xx]); #Ef at (xx+1) 
	
		k=iterativek(sig,h[xx+1])$k;
		Cg[xx+1]=sig/k*0.5*(1+(2*k*h[xx+1]/sinh(2*k*h[xx+1]))); #Group velocity
	#     Gam=0.76*k*h(xx+1)+0.29;% Ruessink Gamma
		Gam=0.5+0.4*tanh(33*So);
		H[xx+1]=sqrt(8*Ef[xx+1]/(rho*g*Cg[xx+1]));#Wave height at (xx+1)
		Hb=0.88/k*tanh(Gam*k*h[xx+1]/0.88);
		A=0.25*rho*g*fp*B;
		temp1=((Hb/H[xx+1])^3+1.5*Hb/H[xx+1])*exp(-(Hb/H[xx+1])^2);
		if(is.finite(H[xx+1])){
			temp2=0.75*sqrt(pi)*(1-erf(Hb/H[xx+1]))
			nogo=FALSE
		} else {
			temp2=0
			nogo=TRUE;
		}
		Db[xx+1]=A*H[xx+1]^3/h[xx+1]*(temp1+temp2); #Dissipation due to brkg
		Dv[xx+1]=rho*Cd*bv[xx+1]*N[xx+1]*(k*g/(2*sig))^3/(2*sqrt(pi))*(sinh(k*hv[xx+1])^3+3*sinh(k*hv[xx+1]))/(3*k*cosh(k*h[xx+1])^3)*H[xx+1]^3;#Diss due to vegetation (MendezLosada2004)
		Df[xx+1]=rho*Cf/(16*sqrt(pi))*(2*pi*fp*H[xx+1]/sinh(k*h[xx+1]))^3;#Diss due to bot friction 
	
		#If there is a flat portion in profile, assume no breaking occurs
		if(x[xx+1] %in% xflat2) Db[xx+1]=0
		if(mean(flat)>0){ #if there is a flat portion
			temp=which(xflat<x[xx+1])
			if(length(temp)>0){
				up=temp[length(temp)]
			} else {
				up=1
			} #loc of last flat portion
		}
		if(h[xx+1]>hd[up]) Db[xx+1]=0  
	}

if(nogo) H=NA

return(H)
  
}



iterativek = function(sigma,hdummy,Uo=0){
	
	#function [k]=iterativek(sigma,hdummy)
	# The following loop performs an iterative solution for k 
	# The calculation is based on a given and constant angular wave frequency
		
	g=9.81;                                              # Gravity (m/s^2)
	b=((sigma^2)*hdummy)/g;                           # (dimensionless)
	kestimated=((sigma^2)/(g*(((tanh(b))^(1/2)))));  # Estimated wavenumber (m^-1) to begin the loop 
														 # (based onEckarts approximation.
	kprevious=0.0000001 
	count = 0;
	while ((abs(kestimated-kprevious) > 0.000005) & (count < 1000)){  # Designates 0.000005 as acceptable error
		count = count+1;
		kh=(kestimated*hdummy);                        # First calculates kh using the estimated k value
		cur = (sigma^2)*(1-Uo*kestimated/sigma)^2;
		kcalculated= cur/(tanh(kh)*g);      # Using the kh value the loop calculates k
		kprevious=kestimated;                        # The loop now sets the 'previous' k value to the 'estimated'
		kestimated=kcalculated;                      # And it sets the 'estimated' value to the 'calculated'
													# The loop then continues by subtracting the
													# two values and testing for the error between them.
	}
	
	if(is.nan(kestimated) | is.na(kestimated)) kcalculated=NA
	k=kcalculated;
		
	n=0.5*(1+(2*k*hdummy/sinh(2*k*hdummy)));
	C=sigma/k;  Cg=C*n; #Group velocity
	
	return(list(k = k, Cg = Cg, n = n))
	
}




