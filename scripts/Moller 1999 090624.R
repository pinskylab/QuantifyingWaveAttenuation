# Moller 1999 equations: sensitivity analysis

moller = function(H1, x, vk, T, h, h1, h2, K, d, f, g){
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

	H2 = H1*Kv*Ks*Kp*Kf		
	return(H2)	
}


H1 = 0.24 # wave height (m)

x = 180 # distance (m)

vk = 1.35e-6 # kinematic viscosity

T = 2.8 # wave period (sec)

h = 1.2 # water depth (m)
h1 = 1.3 # depth at outer edge
h2 = 0.8 # depth at inner edge

K = 1e-10 # specific permeability

d = 1 # sediment depth layer (m)

f = 0.2401 # friction factor marsh
fsand = 0.01 # friction factor sand

g = 9.8 # gravitational constant


H1sens = seq(0.1, 3.3, by=0.1)
xsens = seq(1, 300, by = 10)
Tsens = seq(1, 20, by = 1)
hsens = seq(0.31, 3, by=0.1)
fsens = seq(0.01, 0.5, by=0.01)

height = moller(H1sens, x, vk,  T, h, h1, h2, K, d, f, g)
dist = moller(H1, xsens, vk,  T, h, h1, h2, K, d, f, g)
period = moller(H1, x, vk, Tsens, h, h1, h2, K, d, f, g)
depth = moller(H1, x, vk, T, hsens, hsens+0.1, hsens-0.3, K, d, f, g)
friction = moller(H1, x, vk, T, h, h1, h2, K, d, fsens, g)

# Calcs w/ sand
height_s = moller(H1sens, x, vk,  T, h, h1, h2, K, d, fsand, g)
dist_s = moller(H1, xsens, vk,  T, h, h1, h2, K, d, fsand, g)
period_s = moller(H1, x, vk, Tsens, h, h1, h2, K, d, fsand, g)
depth_s = moller(H1, x, vk, T, hsens, hsens+0.1, hsens-0.3, K, d, fsand, g)


# Plot sensitivity: w/ attenuated wave height
ylim = c(0,0.3)
par(mfrow=c(2,3))
plot(H1sens, height, type="l", main="Incident wave height (m)", ylim=ylim, ylab="Attenuated wave height (m)")
	abline(a = 0, b = 1, lty=2)
	abline(v = H1, col="blue", lty=3)
plot(hsens, depth, type="l", main="Depth", ylim=ylim, ylab="Attenuated wave height (m)")
	abline(h = H1, lty=2)
	abline(v = h, col="blue", lty=3)
plot(xsens, dist, type="l", main="Distance", ylim=ylim, ylab="Attenuated wave height (m)")
	abline(h = H1, lty=2)
	abline(v = x, col="blue", lty=3)
plot(Tsens, period, type="l", main="Period (sec)", ylim=ylim, ylab="Attenuated wave height (m)")
	abline(h = H1, lty=2)
	abline(v = T, col="blue", lty=3)
plot(fsens, friction, type="l", main="Friction factor", ylim=ylim, ylab="Attenuated wave height (m)")
	abline(h = H1, lty=2)
	abline(v = f, col="blue", lty=3)

# Plot sensitivity (% attenuation) w/ marsh or sand
ylim = c(-10,100)
par(mfrow=c(2,3))
plot(H1sens, (H1sens-height)/H1sens*100, type="l", main="Incident wave height", ylim=ylim, ylab="Wave attenuation (%)", xlab="Incident wave height (m)")
	abline(v = H1, col="blue", lty=3)
	lines(H1sens, (H1sens-height_s)/H1sens*100, col="red")
plot(hsens, (H1-depth)/H1*100, type="l", main="Depth", ylim=ylim, ylab="Wave attenuation (%)", xlab="Depth (m)")
	abline(v = h, col="blue", lty=3)
	lines(hsens, (H1-depth_s)/H1*100, col="red")
plot(-1,-1, bty="n", col="white", xaxt="n", yaxt="n", xlab="", ylab="")
	legend("topright", legend=c("Marsh (f=0.24)", "Sand (f=0.01)"), lty=1, col=c("black", "red"), cex=2, bty="n")
plot(xsens, (H1-dist)/H1*100, type="l", main="Distance", ylim=ylim, ylab="Wave attenuation (%)", xlab="Distance (m)")
	abline(v = x, col="blue", lty=3)
	lines(xsens, (H1-dist_s)/H1*100, col="red")
plot(Tsens, (H1-period)/H1*100, type="l", main="Period", ylim=ylim, ylab="Wave attenuation (%)", xlab="Period (sec)")
	abline(v = T, col="blue", lty=3)
	lines(Tsens, (H1-period_s)/H1*100, col="red")
plot(fsens, (H1-friction)/H1*100, type="l", main="Friction factor", ylim=ylim, ylab="Wave attenuation (%)", xlab="Friction factor")
	abline(v = f, col="blue", lty=3)
