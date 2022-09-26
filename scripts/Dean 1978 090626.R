# Dean 1978 equations: sensitivity analysis

dean = function(H1, x, h, CdDS2){
	A = CdDS2*x/(6*pi*h)
	H2 = H1/(1+A*H1)
	return(H2)
}

dean(.1, 20,.15,0.5) # Dean 1978 example = 7.4 cm

H1 = 0.24 # incident wave height (m)
x = 180 # distance (m)
h = 1.2 # water depth (m)
CdDS2 = 0.587 # Moller 1999 coef for marsh = Cd*D/s^2, where Cd is drag coef, D is stalk diam, and s is stalk spacing
CdDS2_sand = 0.0786 # Moller 1999 coef for sand
T = 2.8 # wave period (sec)



H1sens = seq(0.01, 3.01, by=0.1)
xsens = seq(1, 300, by = 10)
hsens = seq(0.1, 3, by=0.1)
Tsens = seq(1, 20, by = 1)
Cdsens = seq(0.01, 1, by=0.01)

height = dean(H1sens, x, h, CdDS2)
dist = dean(H1, xsens, h, CdDS2)
period = rep(dean(H1, x, h, CdDS2), length(Tsens))
depth = dean(H1, x, hsens, CdDS2)
friction = dean(H1, x, h, Cdsens)

# Calcs w/ sand
height_s = dean(H1sens, x, h, CdDS2_sand)
dist_s = dean(H1, xsens, h, CdDS2_sand)
period_s = rep(dean(H1, x, h, CdDS2_sand), length(Tsens))
depth_s = dean(H1, x, hsens, CdDS2_sand)


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
	legend("topright", legend=c("Marsh (CdDS2=0.587)", "Sand (CdDS2=0.0786)"), lty=1, col=c("black", "red"), cex=1, bty="n")
plot(xsens, (H1-dist)/H1*100, type="l", main="Distance", ylim=ylim, ylab="Wave attenuation (%)", xlab="Distance (m)")
	abline(v = x, col="blue", lty=3)
	lines(xsens, (H1-dist_s)/H1*100, col="red")
plot(Tsens, (H1-period)/H1*100, type="l", main="Period", ylim=ylim, ylab="Wave attenuation (%)", xlab="Period (sec)")
	abline(v = T, col="blue", lty=3)
	lines(Tsens, (H1-period_s)/H1*100, col="red")
plot(Cdsens, (H1-friction)/H1*100, type="l", main="Drag factor", ylim=ylim, ylab="Wave attenuation (%)", xlab="Drag factor (CdDS2)")
	abline(v = CdDS2, col="blue", lty=3)
