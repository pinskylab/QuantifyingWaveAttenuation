interp = function(a, b, x=NA){
	# Function to downsample and interpolate
	# a: dataframe with two columns (x and y values)
	# b: dataframe with two columns (x and y values)
	# x: optional vector of x-values to use
	# returns dataframe with three columns: x vals, y-vals from a, y-vals from b
	# higher frequency dataset downsampled to lower frequency

	# make sure high has the higher frequency (smaller timestep) time series
	if(any(is.na(x))){  # use lower frequency x-values if no x-vals specified
		stepa = diff(range(a[,1]))/length(a[,1])
		stepb = diff(range(b[,1]))/length(b[,1])
		if(stepb<=stepa){x=a[,1]; lowname="a"}
		if(stepb>stepa){x=b[,1]; lowname="b"}
	}
	
	# make sure range of a overlaps range of x, else drop some from x
	if(min(a[,1])>min(x) | max(a[,1])<max(x)){ 
		i = which(x<min(a[,1]) | x>max(a[,1]))
		x = x[-i]
		warning(paste("x does not span range of a. Trimmed", length(i), "value(s) off of x"))
	}

	# make sure range of b overlaps range of x, else drop some from x
	if(min(b[,1])>min(x) | max(b[,1])<max(x)){ 
		i = which(x<min(b[,1]) | x>max(b[,1]))
		x = x[-i]
		warning(paste("x does not span range of b. Trimmed", length(i), "value(s) off of x"))
	}
	
	len = length(x)
	new = data.frame(x = x, y1=numeric(len), y2=numeric(len))
	for(i in 1:len){ # interpolate a
		if(x[i] %in% a[,1]){ # if exact x-value is in a
			new$y1[i] = a[median(which(a[,1]==x[i]), na.rm=T),2]
		} else {
			ind1 = max(which(a[,1]<x[i])) # index into a just less than x-val
			ind2 = min(which(a[,1]>x[i])) # index into a just greater than x-val
			m = (a[ind2,2] - a[ind1,2])/(a[ind2,1]-a[ind1,1]) # slope of line from ind1 to ind2
			stepx = x[i] - a[ind1,1] # distance from x-value to a's x-value at ind1
			stepy = stepx*m # calculate the rise for this distance
			new$y1[i] = a[ind1,2]+stepy # add the rise to b's y-value to get predicted value
		}
	}
	for(i in 1:len){ # interpolate b
		if(x[i] %in% b[,1]){ # if exact x-value is in b
			new$y2[i] = b[median(which(b[,1]==x[i]), na.rm=T),2]
		} else {
			ind1 = max(which(b[,1]<x[i])) # index into b just less than x-val
			ind2 = min(which(b[,1]>x[i])) # index into b just greater than x-val
			m = (b[ind2,2] - b[ind1,2])/(b[ind2,1]-b[ind1,1]) # slope of line from ind1 to ind2
			stepx = x[i] - b[ind1,1] # distance from x-value to b's x-value at ind1
			stepy = stepx*m # calculate the rise for this distance
			new$y2[i] = b[ind1,2]+stepy # add the rise to b's y-value to get predicted value
		}
	}
	names(new) = paste(c(names(a), names(b)[2]), c("", "a", "b"), sep="")
	return(new)
}
