## Calculate relative depth for each of the Wave Quant studies


## Start Matlab and test it

system.file("externals", package="R.matlab") # get path to Matlab file, add this to Matlab search path

# then type "MatlabServer" in Matlab

# Create a Matlab client object used to communicate with Matlab
matlab <- Matlab()

# Confirm that the Matlab server is open, and running
isOpen <- open(matlab)
if (!isOpen)
  throw("Matlab server is not running: waited 30 seconds.")

#setVerbose(matlab, -2) # make Matlab verbose
#setVerbose(matlab, 0) # back to normal


evaluate(matlab, "A=1+2;", "B=ones(2,20)")

test <- getVariable(matlab, c("A", "B"))
cat("Received variables:\n")
str(test)
test$A

evaluate(matlab, "A = sin(1.2)")
getVariable(matlab, "A")


## Load data and calculate

setwd("C:/Documents and Settings/Lab/Desktop/Users/Malin/NatCap")
data = read.csv("Coastal Protection Quant Wave 010201.csv", header=T)

evaluate(matlab, "path('C:/Documents and Settings/Lab/Desktop/Users/Malin/NatCap',path);") # set Matlab working dir

kh = rep(NA, length(data[,1]))
for(i in 1:length(data[,1])){
	T = data$Wave.period.mean[i] # wave period
	sig = 2*pi/T # angular frequency
	h = data$Depth.at.incident[i] # depth at outer edge of bed
	if(is.na(h)) h = data$Depth.mean[i] # also try the depth.mean column
	if(!is.na(sig) & !is.na(h)){ # if we have freq and depth
		setVariable(matlab, sig=sig, h=h)
		evaluate(matlab, "a = iterativek(sig, h);")
		kh[i] = unlist(getVariable(matlab, "a"))*h
		#print(paste("sig:", sig, "h:", h, "k:", k[i]))
	}
}

# Define depth categories
data$RelativeDepth = NA
data$RelativeDepth[kh>pi] = "Deep"
data$RelativeDepth[kh<pi/10] = "Shallow"
data$RelativeDepth[kh<=pi & kh>=pi/10] = "Med"
data$RelativeDepth=factor(data$RelativeDepth)


## Summaries
summary(data$RelativeDepth)

aggregate(data$RelativeDepth, by=list(data$Study), FUN=function(x) length(unique(x))) # how many depth categories does each study span?

aggregate(data$RelativeDepth, by=list(data$Study), FUN=function(x) as.numeric(paste(as.numeric(unique(x)), collapse="0"))) # how many depth categories does each study span?

summary(data$RelativeDepth[!duplicated(data$Study)]) # only one line for each study




## When done, clean up and output files

close(matlab)
write.csv(data, paste("Coastal Protection Quant Wave RelDepth ", Sys.Date(), ".csv", sep=""))
