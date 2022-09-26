log10axis = function(side=1, at, cex.axis=1, log=FALSE, labs=TRUE, tcl1 = -0.5, tcl2 = -0.25){
	# at: exponents for 10^at: where to put labels
	# log: if TRUE, the axis is  automatically log-transformed by plot (e.g., plot(log='x'))?
	#		if FALSE, the values are log-transformed, but the axis is not
	# labs: whether to write numbers on the axis
	# tcl1: length of big ticks (fraction of a line of text)
	# tcl2: length of small ticks

	# check that all at values are whole numbers
	if(!all(is.wholenumber(at))) warning("This only works with whole numbers")
	
	# major tick marks
	tks = 10^(min(at):max(at))
	# build the intra-tick marks
	subtks = numeric(0)
	for(i in 1:length(tks)){
		subtks = c(subtks, 2:9*tks[i])
	}
	
	if(log){
		axis(side, at=subtks, labels=FALSE, tcl=tcl2) # short ticks
		axis(side, at=tks, labels=FALSE, tcl=tcl1) # big ticks
		if(labs) axis(side, at=10^at, labels=parse(text=paste("10^", at)), cex.axis=cex.axis, tcl=0)
	} else {
		axis(side, at=log10(subtks), labels=FALSE, tcl= tcl2)
		axis(side, at=log10(tks), labels=FALSE, tcl=tcl1)
		if(labs) axis(side, at=at, labels=parse(text=paste("10^", at)), cex.axis=cex.axis, tcl=0)
	}
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
