setwd('C:/Users/SXP/data/Penium/data_analysis/WGD')
rm(list=ls())
library(mclust)
library(ggplot2)
library(feature)
library("RColorBrewer")

myfunction <- function(output){
	densi <- densityMclust(data$ds)
	summary(densi, parameters = TRUE)
	sumData <-summary(densi, parameters = TRUE)
	#plot(densi, what = "BIC")
	#plot(densi, what = "diagnostic", type = "cdf")
	#plot(densi, what = "diagnostic", type = "qq")

	#br <- seq(min(data[,1]), max(data[,1]), length = 100)
	br <- seq(min(data[,1])-diff(range(data[,1]))/10,max(data[,1])+diff(range(data[,1]))/10, length=binNum)

	pdf(file=paste(output,".pdf"),useDingbats=FALSE)
	#par(mfrow=c(5,1))

	h <- hist(data$ds, breaks = br, plot = FALSE,right=F)
	plot(h, xlab="Ks", freq=TRUE, border=F, col="grey", xlim=range(br),main="Histogram")
	legend("topright", legend=paste("#bin=",binNum),box.lty=0)
	box()

	plot(densi, what="density",xlab="Ks",data=data$ds, breaks=br,col = "slategrey")
	legend("topright", legend=paste("#bin=",binNum),box.lty=0)
	title(main="Kernel density estimation")
	
	col <- adjustcolor(mclust.options("classPlotColors")[1:sumData$G], alpha = 0.3)
	cdens <- predict(densi, br, what = "cdens")
	cdens <- t(apply(cdens, 1, function(d) d*densi$parameters$pro))
	plot(densi, what="density", data=data$ds, breaks=br,col = "slategrey",xlab="Ks")
	matplot(br, cdens, type="l", lwd=1, add=TRUE, lty=1)
	legend("topright", legend=paste("#groups=",densi$G),box.lty=0)
	title(main="fitting with Gaussian distribution")
	
	SiZer(data$ds, signifLevel=0.05, logbw=FALSE, posLegend="bottomright",xlab="Ks")
	title(main="SiZer plot (Significant Zero crossings)")
	SiCon(data$ds, signifLevel=0.05, logbw=FALSE, posLegend="bottomright",xlab="Ks")
	title(main="SiCon plot (SiCon: Significant Convexity)")

	for(i in 1:length(bandwidth)){
		col2 <- brewer.pal(n = 8, name = "RdYlGn")
		data.fs <- featureSignif(data$ds, bw=bandwidth[i], signifLevel=0.05)
		## blue color is significant curvature regions
		plot(data.fs, xlab="Ks",densCol="black",gradCol=col2[1],curvCol=col2[8],addSignifGradRegion=TRUE,addSignifCurvRegion=TRUE,addData=TRUE,addSignifGradData=TRUE,addSignifCurvData=TRUE)
		legend("topright", legend=c("significant gradient region", "significant curvature region"),lty=1,lwd=3,col=c(col2[1],col2[8]),box.lty=0)
		box()
		title(main=paste("Feature significant density plot (bandwidth: ",bandwidth[i],")"))
	}
	
	dev.off()
}

data <- read.table('Penium_MaxOG4.txt',row.names=1,header=T)

#option 1
#data <- subset(data, ds > 0)
#output1 <- 'Penium_MaxOG4_Ks0_2'

#option 2
data <- subset(data, ds > 0.1)
output2 <- 'Penium_MaxOG4_Ks0.1_2'

binNum <- 100
bandwidth <- c(0.1,0.05) # bandwith for feature significant density plot  

myfunction(output2)








