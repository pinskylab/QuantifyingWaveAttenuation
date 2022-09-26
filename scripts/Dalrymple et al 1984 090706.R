# Dalrymple et al 1984 equations: sensitivity analysis

dalrymple = function(H1, x, h, D, b, T, s, g, Cd){
	L = T*sqrt(g*h) # wave length (linear wave theory, shallow water)
	k =  2*pi/L # wave number = 2pi/wavelength
	w = 2*pi/T # radian frequency

	a = 2*Cd/(3*pi)*(D/b)*(H1/b)*((sinh(k*s))^3 + 3*(sinh(k*s))) * ((4*k)/(3*sinh(k*h)*(sinh(2*k*h) + 2*k*h)))
	H2 = H1/(1+a*x)
	atten = (H1-H2)/H1*100
	print(paste("%Atten: ", paste(round(atten, digits=3), collapse=" ")))
	return(H2)
}


g = 9.8 # gravity
H1 = 0.08 # incident wave height (m)
x = 43.3 # distance (m)
h = 1 # water depth (m)
D = 1*10^-2 # width of Thalassia testudinum (from Tussenbroek 1998 Aquatic Botany)
b = 1/sqrt(1100) # spacing of plants (m-1) (B&H report 1100 m-1, but this seems unreasonable, more likely m-2.
T = 1/0.75 # wave period (sec) (peak period in Fig. 4 in Bradley & Houser 2009)
s =  0.225 # height of veg above bottom (m)

L = T*sqrt(g*h) # wave length (linear wave theory, shallow water)
k =  2*pi/L # wave number = 2pi/wavelength
w = 2*pi/T # radian wave frequency

v = 1.35e-6 # kinematic viscosity of water (m2/s)

Cd = 1

#H1sens = seq(0.01, 3.01, by=0.1)
H1sens = seq(0.07, 3, length.out=40)
xsens = seq(1, 300, by = 10)
hsens = seq(0.1, 3, by=0.1)
Tsens = seq(0.5, 20, length.out=80)
Cdsens = 10^seq(-1, 3, length.out=80)
bsens = seq(0.0001, 0.1, length.out=40)


# Expect ~40% attenuation over 43m (B&H 2009)
height = dalrymple(H1=H1sens, x=x, h=h, D=D, b=b, T=T, s=s, g=g, Cd=Cd)
depth = dalrymple(H1=H1, x=x, h=hsens, D=D, b=b, T=T, s=s, g=g, Cd=Cd)
dist = dalrymple(H1=H1, x=xsens, h=h, D=D, b=b, T=T, s=s, g=g, Cd=Cd)
period = dalrymple(H1=H1, x=x, h=h, D=D, b=b, T=Tsens, s=s, g=g, Cd=Cd)
spacing = dalrymple(H1=H1, x=x, h=h, D=D, b=bsens, T=T, s=s, g=g, Cd=Cd)
friction = dalrymple(H1=H1, x=x, h=h, D=D, b=b, T=T, s=s, g=g, Cd=Cdsens)


# Plot sensitivity (% attenuation) w/ marsh or sand
ylim = c(-10,100)
par(mfrow=c(2,3))
plot(H1sens, (H1sens-height)/H1sens*100, type="l", main="Incident wave height", ylim=ylim, ylab="Wave attenuation (%)", xlab="Incident wave height (m)")
	abline(v = H1, col="blue", lty=3)
	abline(v = H1real, col="red", lty=3)
plot(hsens, (H1-depth)/H1*100, type="l", main="Depth", ylim=ylim, ylab="Wave attenuation (%)", xlab="Depth (m)")
	abline(v = h, col="blue", lty=3)
plot(xsens, (H1-dist)/H1*100, type="l", main="Distance", ylim=ylim, ylab="Wave attenuation (%)", xlab="Distance (m)")
	abline(v = x, col="blue", lty=3)
plot(Tsens, (H1-period)/H1*100, type="l", main="Period", ylim=ylim, ylab="Wave attenuation (%)", xlab="Period (sec)")
	abline(v = T, col="blue", lty=3)
plot(Cdsens, (H1-friction)/H1*100, type="l", main="Drag coefficient", ylim=ylim, ylab="Wave attenuation (%)", xlab="Drag coefficient (Cd)", log="x")
	abline(v = Cd, col="blue", lty=3)
plot(bsens, (H1-spacing)/H1*100, type="l", main="Spacing", ylim=ylim, ylab="Wave attenuation (%)", xlab="Spacing (m)")
	abline(v = b, col="blue", lty=3)
