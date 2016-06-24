# Functions Micha Rijkenberg (micharijkenberg@yahoo.com.au)
# To treat ICPMS data
# Version date: 07 June 2016


# Version log:
# 01/07/2015: first version for public use
# 06/07/2015: removed error message when no reference materials are used, implemented the option to show or not show outliers in the boxplot that brings together all checks in one figure
# 07/07/2015: added dotplots to the blanks to visualize the individual data points.
# 15/07/2015: function added allowing to read in files that have been produced in a previous use of the script for a certain ICP-MS run without rerunning the whole script.
# 04/08/2015: the ylim of the figures for the fStanAdd_HNO3 and fStanAdd_Seafast are adapted so that the scale adapts to Rh corrected data
# 07/08/2015: texts "Rh & Mo correction" etc has been added to the figures as text in a more consistent way
# 11/08/2015: units were changed from nmol/L to nmol/kg throughout the functions script.
# 21/10/2015: error messages have been added where it may help the user and reference samples will only be plotted as dotplots in the folder REFERENCE SAMPLES
# 26/10/2015: corrections on 21/10/2015 resulted in an error where the individual dotplots of the non-reference samples were not made for the seafast and ICPMS checks, this has been solved again
# 22/01/2016: Added an error message, you can find the location by using the code "# code 220116_1" with the find function
# 29/02/2016: when using the seafast online it is difficult to do the icpms checks. It is now possible to run the script without having icpms checks.
# 29/02/2015: rf is replaced by ref.mat as dataframe holding the reference consensus values as rf is used in R as a base function
# 09/03/2016: error messages were included to indicate missing or misspelled column names of input files.
# 09/03/2016: no.icpms.data <- 0 has been defined otherwise the missing object no.icpms.data would give an error
# 14/03/2016: bug removed concerning warning of column name index in file sample_list. This file does not have and should not have a column index
# 26/05/2016: the coding is adapted to allow elements indicated with only 1 letter, e.g. "V", "Y" and "U", to be included.
# 27/05/2016: included the possibility to indicate if your work on a mac. windows is than replaced by quartz in the fStElement.stand function
# 31/05/2016: included error messages to make it user friendlier, improved coding for fMo.interference.Cd, allowed variable use of standard.element to correct for evaporation and ICPMS variation in sensitivity
# 03/06/2016: more flexibility has been coded in the function fStElement.stand. You can now fill in any values in the line.high and line.low as long there is a number (at least a 0) and as long line.high > lin.low
# 06/06/2016: error vec.cd couldn't be found in function fMo.interference.Cd. This has been solved, see "# code 060616_1"
# 07/06/2016: when the Standard.element.correction was set on "no" the recoveries couldn't be calculated. This has been changed now, see "# code 070616_1"

######################################################################################################################
#######################################################################################################################
#### Function to make sure that created files in previous script runs can be accessed without rerunning the entire script

fAcces.prev.files <- function(mainDir, name.run,...){

	subDir.Rh <- paste("Figure",standard.element,"standard",name.run,sep="_")
	folder.rh <<- paste(mainDir,subDir.Rh,"/",sep="")

	subDir.fig.standadd <- paste("Figures standard addition HNO3",name.run,sep="_")
	folder.HNO3 <<- paste(mainDir,subDir.fig.standadd,"/",sep="")

	subDir.fig.standadd <- paste("Figures standard addition seafast",name.run,sep="_")
	folder.sf <<- paste(mainDir,subDir.fig.standadd,"/",sep="")

	subDir.fig.Mo <- paste("Figures & data Mo interference",name.run,sep="_")
	folder.Mo <<- paste(mainDir,subDir.fig.Mo,"/",sep="")

	subDir.fig.blank <- paste("Figures & data blank",name.run,sep="_")
	folder.blank <<- paste(mainDir,subDir.fig.blank,"/",sep="")

	subDir.fig.samples <- paste("Figures & data_samples",name.run,sep="_")
	folder.samples <<- paste(mainDir,subDir.fig.samples,"/",sep="")

	print("You can now load the data of a previous run of this dataset without having to rerun all previous functions.")

}

#############################################################################################################
# Function to make a figure for the Rh standard

fStElement.stand <- function(st.el, name.run, mainDir, resolution.figures, cex, cex.lab, cex.axis,...){

options(warn=-1) # to terminate all warning messages

rh <- st.el
rm(rh.result)

# Create new file folder for figures Rh
	subDir.Rh <- paste("Figure",standard.element,"standard",name.run,sep="_")
	dir.create(file.path(mainDir, subDir.Rh))			
	folder.rh <<- paste(mainDir,subDir.Rh,"/",sep="")
	name.fig.rh <- paste(standard.element," data_",name.run,sep="")

# Determine min and max values for Rh figure 1-3, <- rh[grepl("LR", names(rh))]
	rh1 <- rh[grepl("LR", names(rh))]; colnames(rh1) <- c("rh");rh2 <- rh[grepl("MR", names(rh))]; colnames(rh2) <- c("rh")
	rh_line.high.lr <-as.data.frame(line.high.lr); colnames(rh_line.high.lr) <- c("rh")
	rh_line.low.lr <-as.data.frame(line.low.lr); colnames(rh_line.low.lr) <- c("rh")
	rh_line.high.mr <-as.data.frame(line.high.mr); colnames(rh_line.high.mr) <- c("rh")
	rh_line.low.mr <-as.data.frame(line.low.mr); colnames(rh_line.low.mr) <- c("rh")

	# Figure 1
	rh.max <- max(rbind(rh1,rh2), na.rm=TRUE) + 0.1 * max(rbind(rh1,rh2), na.rm=TRUE)
	rh.min <- min(rbind(rh1,rh2), na.rm=TRUE)
	# Figure 2
	rh.max2 <- max(rbind(rh1,rh_line.high.lr,rh_line.low.lr), na.rm=TRUE) + 0.1 * max(rbind(rh1), na.rm=TRUE)
	rh.min2 <- min(rbind(rh1,rh_line.high.lr,rh_line.low.lr), na.rm=TRUE) 
	# Figure 3
	rh.max3 <- max(rbind(rh2,rh_line.high.mr,rh_line.low.mr), na.rm=TRUE) + 0.1 * max(rbind(rh2), na.rm=TRUE)
	rh.min3 <- min(rbind(rh2,rh_line.high.mr,rh_line.low.mr), na.rm=TRUE) 

# Check if index is written correctly
if(!("index" %in% colnames(rh))) {stop("ERROR: The column _index_ is missing or you wrote index with a capital I instead of i in file StandardElement.csv.")}

# Make figure
	tiff(file = paste(folder.rh, name.fig.rh, ".tiff", sep=""), width = 12, height = 12, units = "in", res=resolution.figures)

#windows(12,12)
	par(mfrow=c(3,1))
	par(mar=c(5, 8, 3, 2) + 0.1)		# c(bottom, left, top, right)
	
	plot(rh[,2] ~ rh[,1], 	las = 1, pch=16, ylab="", xlab="index", 
		cex=cex, cex.lab=cex.lab,cex.axis = cex.axis, col="black", main=name.run, ylim=c(rh.min,rh.max))
		title(ylab=paste(standard.element," icpms signal", sep=""), cex.lab=cex.lab, line = 7)
		points(rh[,3] ~ rh[,1],cex=cex,col="red",pch=16)
		mtext(paste(standard.element," at LR",sep=""), side=3, lin=1.6, adj=0, cex=1, col="black")
		mtext(paste(standard.element," at MR",sep=""), side=3, lin=0.2, adj=0, cex=1, col="red")
	
	plot(rh[,2] ~ rh[,1], 	las = 1, pch=16, ylab="", xlab="index", 
		cex=cex, cex.lab=cex.lab,cex.axis = cex.axis, col="black", main=name.run, ylim=c(rh.min2,rh.max2))
		title(ylab=paste(standard.element," icpms signal",sep=""), cex.lab=cex.lab, line = 7)
		mtext(paste(standard.element," at LR",sep=""), side=3, lin=0.5, adj=0, cex=1, col="black")
		lines(x=c(-100,1000),y=c(line.high.lr,line.high.lr), lty=3, lwd=2, col="red")
		lines(x=c(-100,1000),y=c(line.low.lr,line.low.lr), lty=3, lwd=2, col="red")

	plot(rh[,3] ~ rh[,1], 	las = 1, pch=16, ylab="", xlab="index", 
		cex=cex, cex.lab=cex.lab,cex.axis = cex.axis, col="red", main=name.run, ylim=c(rh.min3,rh.max3))
		title(ylab=paste(standard.element," icpms signal",sep=""), cex.lab=cex.lab, line = 7)
		mtext(paste(standard.element," at MR",sep=""), side=3, lin=0.5, adj=0, cex=1, col="red")
		lines(x=c(-100,1000),y=c(line.high.mr,line.high.mr), lty=3, lwd=2, col="black")
		lines(x=c(-100,1000),y=c(line.low.mr,line.low.mr), lty=3, lwd=2, col="black")

	dev.off()

# Same figure to be plotted in the window

if(Using.computer == "windows") {windows(12,12)}
if(Using.computer == "mac") {quartz(width=12,height=12)}
	par(mfrow=c(3,1))
	par(mar=c(5, 8, 3, 2) + 0.1)		# c(bottom, left, top, right)
	
	plot(rh[,2] ~ rh[,1], 	las = 1, pch=16, ylab="", xlab="index", 
		cex=cex, cex.lab=cex.lab,cex.axis = cex.axis, col="black", main=name.run, ylim=c(rh.min,rh.max))
		title(ylab=paste(standard.element," icpms signal", sep=""), cex.lab=cex.lab, line = 7)
		points(rh[,3] ~ rh[,1],cex=cex,col="red",pch=16)
		mtext(paste(standard.element," at LR",sep=""), side=3, lin=1.6, adj=0, cex=1, col="black")
		mtext(paste(standard.element," at MR",sep=""), side=3, lin=0.2, adj=0, cex=1, col="red")
	
	plot(rh[,2] ~ rh[,1], 	las = 1, pch=16, ylab="", xlab="index", 
		cex=cex, cex.lab=cex.lab,cex.axis = cex.axis, col="black", main=name.run, ylim=c(rh.min2,rh.max2))
		title(ylab=paste(standard.element," icpms signal",sep=""), cex.lab=cex.lab, line = 7)
		mtext(paste(standard.element," at LR",sep=""), side=3, lin=0.5, adj=0, cex=1, col="black")
		lines(x=c(-100,1000),y=c(line.high.lr,line.high.lr), lty=3, lwd=2, col="red")
		lines(x=c(-100,1000),y=c(line.low.lr,line.low.lr), lty=3, lwd=2, col="red")

	plot(rh[,3] ~ rh[,1], 	las = 1, pch=16, ylab="", xlab="index", 
		cex=cex, cex.lab=cex.lab,cex.axis = cex.axis, col="red", main=name.run, ylim=c(rh.min3,rh.max3))
		title(ylab=paste(standard.element," icpms signal",sep=""), cex.lab=cex.lab, line = 7)
		mtext(paste(standard.element," at MR",sep=""), side=3, lin=0.5, adj=0, cex=1, col="red")
		lines(x=c(-100,1000),y=c(line.high.mr,line.high.mr), lty=3, lwd=2, col="black")
		lines(x=c(-100,1000),y=c(line.low.mr,line.low.mr), lty=3, lwd=2, col="black")

## Select the different concentration ranges head(rh.result)
# LR

# Warning messages
if(mean(rh[,grepl("LR", names(rh))], na.rm=TRUE) == mean(rh[,2], na.rm=TRUE)) {} else {
		stop(paste("ERROR: The columns in the file StandardElement.csv should be in the order _index_ _",standard.element,"xxx(LR)_ _",standard.element,"xxx(MR)_.",sep=""))}  

if(line.low.lr > line.high.lr) 	{	stop("ERROR: the value of line.low.lr needs to be lower than line.high.lr")}

test1 <- subset(rh, rh[,grepl("LR", names(rh))] < line.high.lr & rh[,grepl("LR", names(rh))] > line.low.lr)
test2 <- subset(rh, rh[,grepl("MR", names(rh))] < line.high.mr & rh[,grepl("MR", names(rh))] > line.low.mr)


# First code in the case when someone selects two groups and put line.high.lr above all datapoints # source("Functions_SeafastICPMS_DataAnalysis_030616.r")


if(line.high.lr > max(rh1, na.rm=TRUE)){


	if(line.high.lr > 0 & line.low.lr > 0){
							
			if(dim(test1)[1] < 2) { stop("ERROR: there are not enough datapoints betweem line.high.lr and line.low.lr to calculate a mean.")}			

							Rh.lr.high <- subset(rh, rh[grepl("LR", names(rh))] < line.high.lr & rh[grepl("LR", names(rh))] > line.low.lr)
							Rh.lr.high[,paste("mean.",standard.element,".range.LR", sep="")] <- mean(Rh.lr.high[,2], na.rm=TRUE)
							Rh.lr.high <- Rh.lr.high[,c("index",paste("mean.",standard.element,".range.LR", sep=""))]

							Rh.lr.low <- subset(rh, rh[grepl("LR", names(rh))] < line.low.lr)
							Rh.lr.low[,paste("mean.",standard.element,".range.LR", sep="")] <- mean(Rh.lr.low[,2], na.rm=TRUE)
							Rh.lr.low <- Rh.lr.low[,c("index",paste("mean.",standard.element,".range.LR", sep=""))]

							Rh.lr <- rbind(Rh.lr.high, Rh.lr.low) 
							rh.result <- merge(rh, Rh.lr, by = "index", all.x=TRUE)

								}


	if(line.high.lr > 0 & line.low.lr == 0){
					
			if(dim(test1)[1] < 2) { stop("ERROR: there are not enough datapoints betweem line.high.lr and line.low.lr to calculate a mean.")}			
		
							Rh.lr.high <- subset(rh, rh[grepl("LR", names(rh))] < line.high.lr & rh[grepl("LR", names(rh))] > line.low.lr)
							Rh.lr.high[,paste("mean.",standard.element,".range.LR", sep="")] <- mean(Rh.lr.high[,2], na.rm=TRUE)
							Rh.lr.high <- Rh.lr.high[,c("index",paste("mean.",standard.element,".range.LR", sep=""))]

							rh.result <- merge(rh, Rh.lr.high, by = "index", all.x=TRUE)


								}

}

if(line.high.lr < max(rh1, na.rm=TRUE)){


	if(line.high.lr > 0 & line.low.lr > 0){
							
			if(dim(test1)[1] < 2) { stop("ERROR: there are not enough datapoints betweem line.high.lr and line.low.lr to calculate a mean.")}			

							Rh.lr.high <- subset(rh, rh[grepl("LR", names(rh))] > line.high.lr)
							Rh.lr.high[,paste("mean.",standard.element,".range.LR", sep="")] <- mean(Rh.lr.high[,2], na.rm=TRUE)
							Rh.lr.high <- Rh.lr.high[,c("index",paste("mean.",standard.element,".range.LR", sep=""))]

							Rh.lr.med <- subset(rh, rh[,2] < line.high.lr & rh[,2] > line.low.lr)
							Rh.lr.med[,paste("mean.",standard.element,".range.LR", sep="")] <- mean(Rh.lr.med[,2], na.rm=TRUE)
							Rh.lr.med <- Rh.lr.med[,c("index",paste("mean.",standard.element,".range.LR", sep=""))]

							Rh.lr.low <- subset(rh, rh[,2] < line.low.lr)
							Rh.lr.low[,paste("mean.",standard.element,".range.LR", sep="")] <- mean(Rh.lr.low[,2], na.rm=TRUE)
							Rh.lr.low <- Rh.lr.low[,c("index",paste("mean.",standard.element,".range.LR", sep=""))]

							Rh.lr <- rbind(Rh.lr.high, Rh.lr.med, Rh.lr.low) 
							rh.result <- merge(rh, Rh.lr, by = "index", all.x=TRUE)
							}

	if(line.high.lr > 0 & line.low.lr == 0){
							
			if(dim(test1)[1] < 2) { stop("ERROR: there are not enough datapoints betweem line.high.lr and line.low.lr to calculate a mean.")}			

							Rh.lr.high <- subset(rh, rh[,2] > line.high.lr)
							Rh.lr.high[,paste("mean.",standard.element,".range.LR", sep="")] <- mean(Rh.lr.high[,2], na.rm=TRUE)
							Rh.lr.high <- Rh.lr.high[,c("index",paste("mean.",standard.element,".range.LR", sep=""))]

							Rh.lr.low <- subset(rh, rh[,2] < line.high.lr)
							Rh.lr.low[,paste("mean.",standard.element,".range.LR", sep="")] <- mean(Rh.lr.low[,2], na.rm=TRUE)
							Rh.lr.low <- Rh.lr.low[,c("index",paste("mean.",standard.element,".range.LR", sep=""))]

							Rh.lr <- rbind(Rh.lr.high, Rh.lr.low) 
							rh.result <- merge(rh, Rh.lr, by = "index", all.x=TRUE)
							}

} #for if(line.high.lr < max(rh1, na.rm=TRUE))

	if(line.high.lr == 0 & line.low.lr == 0){
							Rh.lr.mean <- mean(rh[,2], na.rm=TRUE)
							rh.result <- rh
							rh.result[,paste("mean.",standard.element,".range.LR", sep="")] <- Rh.lr.mean
							}

# MR head(Rh.mr.high)

if(line.low.mr > line.high.mr) 	{	stop("ERROR: the value of line.low.mr needs to be lower than line.high.mr")}

if(line.high.mr > max(rh2, na.rm=TRUE)){


	if(line.high.mr > 0 & line.low.mr > 0){

			if(dim(test2)[1] < 2) { stop("ERROR: there are not enough datapoints betweem line.high.mr and line.low.mr to calculate a mean.")}			

							Rh.mr.high <- subset(rh, rh[grepl("MR", names(rh))] < line.high.mr & rh[grepl("MR", names(rh))] > line.low.mr)
							Rh.mr.high[,paste("mean.",standard.element,".range.MR", sep="")] <- mean(Rh.mr.high[,2], na.rm=TRUE)
							Rh.mr.high <- Rh.mr.high[,c("index",paste("mean.",standard.element,".range.MR", sep=""))]

							Rh.mr.low <- subset(rh, rh[grepl("MR", names(rh))] < line.low.mr)
							Rh.mr.low[,paste("mean.",standard.element,".range.MR", sep="")] <- mean(Rh.mr.low[,2], na.rm=TRUE)
							Rh.mr.low <- Rh.mr.low[,c("index",paste("mean.",standard.element,".range.MR", sep=""))]

							Rh.mr <- rbind(Rh.mr.high, Rh.mr.low) 
							rh.result <- merge(rh.result, Rh.mr, by = "index", all.x=TRUE)

								}


	if(line.high.mr > 0 & line.low.mr == 0){
	
			if(dim(test2)[1] < 2) { stop("ERROR: there are not enough datapoints betweem line.high.mr and line.low.mr to calculate a mean.")}			
							
							Rh.mr.high <- subset(rh, rh[grepl("MR", names(rh))] < line.high.mr & rh[grepl("MR", names(rh))] > line.low.mr)
							Rh.mr.high[,paste("mean.",standard.element,".range.MR", sep="")] <- mean(Rh.mr.high[,2], na.rm=TRUE)
							Rh.mr.high <- Rh.mr.high[,c("index",paste("mean.",standard.element,".range.MR", sep=""))]

							rh.result <- merge(rh.result, Rh.mr.high, by = "index", all.x=TRUE)



								}

}


if(line.high.mr < max(rh2, na.rm=TRUE)){

	if(line.high.mr > 0 & line.low.mr > 0){

			if(dim(test2)[1] < 2) { stop("ERROR: there are not enough datapoints betweem line.high.mr and line.low.mr to calculate a mean.")}			

							Rh.mr.high <- subset(rh, rh[,3] > line.high.mr)
							Rh.mr.high[,paste("mean.",standard.element,".range.MR", sep="")] <- mean(Rh.mr.high[,3], na.rm=TRUE)
							Rh.mr.high <- Rh.mr.high[,c("index",paste("mean.",standard.element,".range.MR", sep=""))]

							Rh.mr.med <- subset(rh, rh[,3] < line.high.mr & rh[,3] > line.low.mr)
							Rh.mr.med[,paste("mean.",standard.element,".range.MR", sep="")] <- mean(Rh.mr.med[,3], na.rm=TRUE)
							Rh.mr.med <- Rh.mr.med[,c("index",paste("mean.",standard.element,".range.MR", sep=""))]

							Rh.mr.low <- subset(rh, rh[,3] < line.low.mr)
							Rh.mr.low[,paste("mean.",standard.element,".range.MR", sep="")] <- mean(Rh.mr.low[,3], na.rm=TRUE)
							Rh.mr.low <- Rh.mr.low[,c("index",paste("mean.",standard.element,".range.MR", sep=""))]

							Rh.mr <- rbind(Rh.mr.high, Rh.mr.med, Rh.mr.low) 
							rh.result <- merge(rh.result, Rh.mr, by = "index", all.x=TRUE)
							}


	if(line.high.mr > 0 & line.low.mr == 0){

			if(dim(test2)[1] < 2) { stop("ERROR: there are not enough datapoints betweem line.high.mr and line.low.mr to calculate a mean.")}			

							Rh.mr.high <- subset(rh, rh[,3] > line.high.mr)
							Rh.mr.high[,paste("mean.",standard.element,".range.MR", sep="")] <- mean(Rh.mr.high[,3], na.rm=TRUE)
							Rh.mr.high <- Rh.mr.high[,c("index",paste("mean.",standard.element,".range.MR", sep=""))]

							Rh.mr.low <- subset(rh, rh[,3] < line.high.mr)
							Rh.mr.low[,paste("mean.",standard.element,".range.MR", sep="")] <- mean(Rh.mr.low[,3], na.rm=TRUE)
							Rh.mr.low <- Rh.mr.low[,c("index",paste("mean.",standard.element,".range.MR", sep=""))]

							Rh.mr <- rbind(Rh.mr.high, Rh.mr.low) 
							rh.result <- merge(rh.result, Rh.mr, by = "index", all.x=TRUE)
							}

} # for if(line.high.mr < max(rh2, na.rm=TRUE))

	if(line.high.mr == 0 & line.low.mr == 0){
							Rh.mr.mean <- mean(rh[,3], na.rm=TRUE)
							rh.result[,paste("mean.",standard.element,".range.MR", sep="")] <- Rh.mr.mean
							}

write.table(rh.result,file=paste(folder.rh,standard.element," standard data_",name.run,".csv", sep=""),sep=",",row.names=F)

print(head(rh.result))


}

#############################################################################################################
# Function to calculate the standard addition in HNO3 for the ICPMS


fStandAdd.HNO3 <- function(ic, ad, rh, name.run, mainDir, resolution.figures, cex, cex.lab, cex.axis,Standard.element.correction,...){

Rh.correction <- Standard.element.correction
rh <- st.el

# Check if column names are written correctly
if(!("index" %in% colnames(ic))) {stop("ERROR: The column _index_ is missing in one of the input files or you wrote index with a capital I instead of i in the file StanAdd_HNO3.csv.")}
if(!("sample" %in% colnames(ic))) {stop("ERROR: The column _sample_ is missing in the file StanAdd_HNO3.csv or is written incorrectly / don't use capital letters")}
if(!("repetition" %in% colnames(ic))) {stop("ERROR: The column _repetition_ is missing in the file StanAdd_HNO3.csv or is written incorrectly / don't use capital letters")}
if(!("conc.order" %in% colnames(ic))) {stop("ERROR: The column _conc.order_ is missing in the file StanAdd_HNO3.csv or is written incorrectly / don't use capital letters")}
if(!("conc.order" %in% colnames(ad))) {stop("ERROR: The column _conc.order_ is missing in the file Added_HNO3.csv or is written incorrectly / don't use capital letters")}
if(!("index" %in% colnames(rh))) {stop("ERROR: The column _index_ is missing in one of the input files or you wrote index with a capital I instead of i in file Rh_standard.csv.")}

library(stringr)
options(warn=-1) # to terminate all warning messages

vec.ad <- names(ad[2:dim(ad)[2]])
vec.repetition <- unique(ic$repetition)
	# to remove any NA's that excel may initiate
	vc <- is.na(vec.repetition)
	vec.repetition <- vec.repetition[vc==FALSE]
vec.ic <- names(ic)[5:dim(ic)[2]]

# Create table for data head(data.standadd)
data.standadd <- as.data.frame(matrix(NA,ncol=length(vec.ic)+1,nrow=10))
colnames(data.standadd) <- c("parameter",vec.ic)

#head(st)

for(i in 1:length(vec.ic)){
#i=5

			# select data
			st1 <- ic[,c("index","sample","repetition","conc.order",vec.ic[i])]
			nm <- str_sub(vec.ic[i],1,2)	#library(stringr)
			
			# to remove any row's with NA's that excel may initiate
				z <- is.na(st1[,c("index")])
				st1 <- st1[z==FALSE,]

if(Rh.correction == "yes"){
				# selection of LR or MR Rh standard data
				ne <- str_sub(vec.ic[i],-3,-2)

				if(exists("ne") == FALSE) {print("ERROR: you are working in an older version of R that does not support the stringr library!")}				
				if(ne == "") {print("ERROR: you are working in an older version of R that does not support the stringr library!")}				

				if(nm == paste(standard.element,sep="")) {next}				

				if(ne == "LR") {Rh.corr.a <- rh[grepl("LR", names(rh))]; index <- rh$index; Rh.corr <- cbind(index,Rh.corr.a); name.1 <- names(Rh.corr)[2]}
				if(ne == "MR") {Rh.corr.a <- rh[grepl("MR", names(rh))]; index <- rh$index; Rh.corr <- cbind(index,Rh.corr.a); name.1 <- names(Rh.corr)[2]} 

				st <- merge(st1, Rh.corr, by = "index", all.y=FALSE)
						
		# Correct data for Rh standard
		if(ne == "LR") {st$Me.rh.corrected <- (st[,c(vec.ic[i])]*st[,c(paste("mean.",standard.element,".range.LR",sep=""))])/st[,c(name.1)]}		
		if(ne == "MR") {st$Me.rh.corrected <- (st[,c(vec.ic[i])]*st[,c(paste("mean.",standard.element,".range.MR",sep=""))])/st[,c(name.1)]}
				} else {st <- st1}
	
			# add the added standard concentrations
			if(nm %in% colnames(ad)) {} else {
			nm <- str_sub(vec.ic[i],1,1)     # it is an element indicated with one letter only, e.g. V
			if(nm %in% colnames(ad)) {} else {next}
			}
			tad <- ad[,c("conc.order",nm)]      
			st <- merge(st, tad, by = "conc.order", all.y=FALSE)

			# units for the different metals
				if(nm %in% units.nmolKg) {units <- "nmol/kg"}
				if(nm %in% units.pmolKg) {units <- "pmol/kg"}
			
			# Create new file folder for figures standard addition
			subDir.fig.standadd <- paste("Figures standard addition HNO3",name.run,sep="_")
			dir.create(file.path(mainDir, subDir.fig.standadd))			
			folder.HNO3 <<- paste(mainDir,subDir.fig.standadd,"/",sep="")

			# make figures standard addition
			name.fig <- paste(vec.ic[i],name.run,sep="_")
			tiff(file = paste(folder.HNO3, name.fig, ".tiff", sep=""), width = 5, height = 5, units = "in", res=resolution.figures)
			
			#windows(5,5)
			par(mfrow=c(1,1))
			par(mar=c(5, 6, 3, 2) + 0.1)		# c(bottom, left, top, right)
		
			if(Rh.correction == "yes")	{
				y.max <- max(c(st[,c("Me.rh.corrected")],st[,c(vec.ic[i])]), na.rm=TRUE)
				y.min <- min(c(st[,c("Me.rh.corrected")],st[,c(vec.ic[i])]), na.rm=TRUE)
					}    else 	{
				y.max <- max(c(st[,c(vec.ic[i])]), na.rm=TRUE)
				y.min <- min(c(st[,c(vec.ic[i])]), na.rm=TRUE)
							}
	
			plot(st[,c(vec.ic[i])] ~ st[,nm], 	las = 1, pch=1, ylab="", xlab=paste("added ",nm," (",units,")",sep=""), 
									cex=cex, cex.lab=cex.lab,cex.axis = cex.axis, col="black", main=name.run, 
									ylim=c(y.min,y.max))
					title(ylab=vec.ic[i], cex.lab=cex.lab, line = 4.5)
					model <- lm(st[,c(vec.ic[i])] ~ st[,nm])
					abline(model)
							slope <- round(summary(model)[4][[1]][[2]],2)
							slope.StdError <- round(summary(model)[4][[1]][[4]],2)
							intercept <- round(summary(model)[4][[1]][[1]],2)
							intercept.StdError <- round(summary(model)[4][[1]][[3]],2)
							Rsq <- round(summary(model)[9][[1]],3)
			
				mtext("Standard addition in HNO3", side=3, lin=0.1, adj=0.5, cex=0.8, col="black")
				mtext(paste("y = ",slope," x + ",intercept,sep=""), side=3, lin=-1, adj=0.1, cex=0.8, col="black")
				mtext(paste("R^2 = ",Rsq, sep=""), side=3, lin=-2, adj=0.07, cex=0.8, col="black")

if(Rh.correction != "yes"){
					legend("bottomright", legend = c("uncorrected for Rh"), 
					col=c("black"),pch=1, cex=1, bty="n")
					}


if(Rh.correction == "yes"){
					points(st[,c("Me.rh.corrected")] ~ st[,nm], pch=1, col="red")
					model.rh <- lm(st[,c("Me.rh.corrected")] ~ st[,nm])
					abline(model.rh,col="red")
							slope.rh <- round(summary(model.rh)[4][[1]][[2]],2)
							slope.StdError.rh <- round(summary(model.rh)[4][[1]][[4]],2)
							intercept.rh <- round(summary(model.rh)[4][[1]][[1]],2)
							intercept.StdError.rh <- round(summary(model.rh)[4][[1]][[3]],2)
							Rsq.rh <- round(summary(model.rh)[9][[1]],3)
					mtext(paste("y = ",slope.rh," x + ",intercept.rh,sep=""), side=3, lin=-3, adj=0.1, cex=0.8, col="red")
					mtext(paste("R^2 = ",Rsq.rh, sep=""), side=3, lin=-4, adj=0.07, cex=0.8, col="red")
					legend("bottomright", legend = c(paste("uncorrected for ",standard.element,sep=""), paste("corrected for ",standard.element,sep="")), 
					col=c("black","red"),pch=1, cex=0.8, bty="n")

}

			dev.off() 	

			# Bring data together in table

			data.standadd[1,1] <- "slope"
			data.standadd[2,1] <- "slope.StdError"
			data.standadd[3,1] <- "intercept"
			data.standadd[4,1] <- "intercept.StdError"
			data.standadd[5,1] <- "R^2"
			data.standadd[6,1] <- paste("slope.",standard.element,".corr", sep="")
			data.standadd[7,1] <- paste("slope.StdError.",standard.element,".corr", sep="")
			data.standadd[8,1] <- paste("intercept.",standard.element,".corr", sep="")
			data.standadd[9,1] <- paste("intercept.StdError.",standard.element,".corr", sep="")
			data.standadd[10,1] <- paste("R^2.",standard.element,".corr", sep="")

			data.standadd[1,i+1] <- slope
			data.standadd[2,i+1] <- slope.StdError
			data.standadd[3,i+1] <- intercept
			data.standadd[4,i+1] <- intercept.StdError
			data.standadd[5,i+1] <- Rsq
			
if(Rh.correction == "yes"){

			data.standadd[6,i+1] <- slope.rh
			data.standadd[7,i+1] <- slope.StdError.rh
			data.standadd[8,i+1] <- intercept.rh
			data.standadd[9,i+1] <- intercept.StdError.rh
			data.standadd[10,i+1] <- Rsq.rh
					}

graphics.off()

print(paste(i," & ", vec.ic[i],sep=""))


}

write.table(data.standadd,file=paste(folder.HNO3,"StandAdd_HNO3_results_",name.run,".csv", sep=""),sep=",",row.names=F)

print(data.standadd)

} # end function

#############################################################################################################
# Function to calculate the standard addition in seawater extracted by seafast and measured on the icpms


fStandAdd.seafast <- function(sf, ad.sf, st.el, name.run, mainDir, resolution.figures, cex, cex.lab, cex.axis,Standard.element.correction,...){

library(stringr)
options(warn=-1) # to terminate all warning messages

Rh.correction <- Standard.element.correction
rh <- st.el

# Check if column names are written correctly
if(!("index" %in% colnames(sf))) {stop("ERROR: The column _index_ is missing in one of the input files or you wrote index with a capital I instead of i in the file StanAdd_seafast.csv.")}
if(!("sample" %in% colnames(sf))) {stop("ERROR: The column _sample_ is missing in the file StanAdd_seafast.csv or is written incorrectly / don't use capital letters")}
if(!("repetition" %in% colnames(sf))) {stop("ERROR: The column _repetition_ is missing in the file StanAdd_seafast.csv or is written incorrectly / don't use capital letters")}
if(!("conc.order" %in% colnames(sf))) {stop("ERROR: The column _conc.order_ is missing in the file StanAdd_seafast.csv or is written incorrectly / don't use capital letters")}
if(!("conc.order" %in% colnames(ad.sf))) {stop("ERROR: The column _conc.order_ is missing in the file Added_seafast.csv or is written incorrectly / don't use capital letters")}
if(!("index" %in% colnames(rh))) {stop("ERROR: The column _index_ is missing in one of the input files or you wrote index with a capital I instead of i in file Rh_standard.csv.")}

vec.ad.sf <- names(ad.sf[2:dim(ad.sf)[2]])
vec.repetition <- unique(sf$repetition)
	# to remove any NA's that excel may initiate
	vc <- is.na(vec.repetition)
	vec.repetition <- vec.repetition[vc==FALSE]
vec.sf <- names(sf)[5:dim(sf)[2]]

# Create table for data head(data.standadd)
data.standadd <- as.data.frame(matrix(NA,ncol=length(vec.sf)+1,nrow=10))
colnames(data.standadd) <- c("parameter",vec.sf)

# Create new file folder for figures standard addition
			subDir.fig.standadd <- paste("Figures standard addition seafast",name.run,sep="_")
			dir.create(file.path(mainDir, subDir.fig.standadd))			
			folder.sf <<- paste(mainDir,subDir.fig.standadd,"/",sep="")

for(i in 1:length(vec.sf)){
#i=5

			# select data
			st1 <- sf[,c("index","sample","repetition","conc.order",vec.sf[i])]
			nm <- str_sub(vec.sf[i],1,2)	#library(stringr)
			
			# to remove any row's with NA's that excel may initiate
				z <- is.na(st1[,c("index")])
				st1 <- st1[z==FALSE,]

if(Rh.correction == "yes"){
				# selection of LR or MR Rh standard data
				ne <- str_sub(vec.sf[i],-3,-2)

				if(exists("ne") == FALSE) {print("ERROR: you are working in an older version of R that does not support the stringr library!")}				
				if(ne == "") {print("ERROR: you are working in an older version of R that does not support the stringr library!")}				

				if(nm == paste(standard.element,sep="")) {next}				

				if(ne == "LR") {Rh.corr.a <- rh[grepl("LR", names(rh))]; index <- rh$index; Rh.corr <- cbind(index,Rh.corr.a); name.1 <- names(Rh.corr)[2]}
				if(ne == "MR") {Rh.corr.a <- rh[grepl("MR", names(rh))]; index <- rh$index; Rh.corr <- cbind(index,Rh.corr.a); name.1 <- names(Rh.corr)[2]} 

				st <- merge(st1, Rh.corr, by = "index", all.y=FALSE)
						
		# Correct data for Rh standard
		if(ne == "LR") {st$Me.rh.corrected <- (st[,c(vec.sf[i])]*st[,c(paste("mean.",standard.element,".range.LR",sep=""))])/st[,c(name.1)]}		
		if(ne == "MR") {st$Me.rh.corrected <- (st[,c(vec.sf[i])]*st[,c(paste("mean.",standard.element,".range.MR",sep=""))])/st[,c(name.1)]}
				} else {st <- st1}
	
			# add the added standard concentrations
			if(nm %in% colnames(ad)) {} else {
			nm <- str_sub(vec.sf[i],1,1)     # it is an element indicated with one letter only, e.g. V
			if(nm %in% colnames(ad)) {} else {next}
			}
			tad <- ad.sf[,c("conc.order",nm)]      
			st <- merge(st, tad, by = "conc.order", all.y=FALSE)

			# units for the different metals
				if(nm %in% units.nmolKg) {units <- "nmol/kg"}
				if(nm %in% units.pmolKg) {units <- "pmol/kg"}

			# make figures standard addition
			name.fig <- paste(vec.sf[i],name.run,sep="_")
			tiff(file = paste(folder.sf, name.fig, ".tiff", sep=""), width = 5, height = 5, units = "in", res=resolution.figures)
			
			#windows(5,5)
			par(mfrow=c(1,1))
			par(mar=c(5, 6, 3, 2) + 0.1)		# c(bottom, left, top, right)

			if(Rh.correction == "yes")	{
				y.max <- max(c(st[,c("Me.rh.corrected")],st[,c(vec.sf[i])]), na.rm=TRUE)
				y.min <- min(c(st[,c("Me.rh.corrected")],st[,c(vec.sf[i])]), na.rm=TRUE)
					}    else 	{
				y.max <- max(c(st[,c(vec.sf[i])]), na.rm=TRUE)
				y.min <- min(c(st[,c(vec.sf[i])]), na.rm=TRUE)
							}

			plot(st[,c(vec.sf[i])] ~ st[,nm], 	las = 1, pch=1, ylab="", xlab=paste("added ",nm," (",units,")",sep=""), 
									cex=cex, cex.lab=cex.lab,cex.axis = cex.axis, col="black", main=name.run, 
									ylim=c(y.min,y.max))
					title(ylab=vec.sf[i], cex.lab=cex.lab, line = 4.5)
					model <- lm(st[,c(vec.sf[i])] ~ st[,nm])
					abline(model)
							slope <- round(summary(model)[4][[1]][[2]],2)
							slope.StdError <- round(summary(model)[4][[1]][[4]],2)
							intercept <- round(summary(model)[4][[1]][[1]],2)
							intercept.StdError <- round(summary(model)[4][[1]][[3]],2)
							Rsq <- round(summary(model)[9][[1]],3)
			
				mtext("Standard addition seafast", side=3, lin=0.1, adj=0.5, cex=0.8, col="black")
				mtext(paste("y = ",slope," x + ",intercept,sep=""), side=3, lin=-1, adj=0.1, cex=0.8, col="black")
				mtext(paste("R^2 = ",Rsq, sep=""), side=3, lin=-2, adj=0.07, cex=0.8, col="black")

if(Rh.correction != "yes"){
					legend("bottomright", legend = c("uncorrected for Rh"), 
					col=c("black"),pch=1, cex=1, bty="n")
					}


if(Rh.correction == "yes"){
					points(st[,c("Me.rh.corrected")] ~ st[,nm], pch=1, col="red")
					model.rh <- lm(st[,c("Me.rh.corrected")] ~ st[,nm])
					abline(model.rh,col="red")
							slope.rh <- round(summary(model.rh)[4][[1]][[2]],2)
							slope.StdError.rh <- round(summary(model.rh)[4][[1]][[4]],2)
							intercept.rh <- round(summary(model.rh)[4][[1]][[1]],2)
							intercept.StdError.rh <- round(summary(model.rh)[4][[1]][[3]],2)
							Rsq.rh <- round(summary(model.rh)[9][[1]],3)
					mtext(paste("y = ",slope.rh," x + ",intercept.rh,sep=""), side=3, lin=-3, adj=0.1, cex=0.8, col="red")
					mtext(paste("R^2 = ",Rsq.rh, sep=""), side=3, lin=-4, adj=0.07, cex=0.8, col="red")
					legend("bottomright", legend = c(paste("uncorrected for ",standard.element,sep=""), paste("corrected for ",standard.element,sep="")), 
					col=c("black","red"),pch=1, cex=0.8, bty="n")

}

			dev.off() 	

			# Bring data together in table

			data.standadd[1,1] <- "slope"
			data.standadd[2,1] <- "slope.StdError"
			data.standadd[3,1] <- "intercept"
			data.standadd[4,1] <- "intercept.StdError"
			data.standadd[5,1] <- "R^2"
			data.standadd[6,1] <- paste("slope.",standard.element,".corr", sep="")
			data.standadd[7,1] <- paste("slope.StdError.",standard.element,".corr", sep="")
			data.standadd[8,1] <- paste("intercept.",standard.element,".corr", sep="")
			data.standadd[9,1] <- paste("intercept.StdError.",standard.element,".corr", sep="")
			data.standadd[10,1] <- paste("R^2.",standard.element,".corr", sep="")


			data.standadd[1,i+1] <- slope
			data.standadd[2,i+1] <- slope.StdError
			data.standadd[3,i+1] <- intercept
			data.standadd[4,i+1] <- intercept.StdError
			data.standadd[5,i+1] <- Rsq
			
if(Rh.correction == "yes"){

			data.standadd[6,i+1] <- slope.rh
			data.standadd[7,i+1] <- slope.StdError.rh
			data.standadd[8,i+1] <- intercept.rh
			data.standadd[9,i+1] <- intercept.StdError.rh
			data.standadd[10,i+1] <- Rsq.rh
					}

graphics.off()

print(paste(i," & ", vec.sf[i],sep=""))


}

write.table(data.standadd,file=paste(folder.sf,"StandAdd_seafast_results_",name.run,".csv", sep=""),sep=",",row.names=F)

print(data.standadd)

} # end function

#####################################################################################################################################
##### Mo interference Cd function ###################################


fMo.interference.Cd <- function(mo, st.el, resolution.figures, cex, cex.lab, cex.axis, Standard.element.correction,...){

# Create new file folder for figures standard addition
			subDir.fig.Mo <- paste("Figures & data Mo interference",name.run,sep="_")
			dir.create(file.path(mainDir, subDir.fig.Mo))			
			folder.Mo <<- paste(mainDir,subDir.fig.Mo,"/",sep="")

Rh.correction <- Standard.element.correction
rh <- st.el
mi <- mo

test1 <- mi[grepl("Cd", names(mi))]
if(dim(test1)[2] > 1) {stop("ERROR: use only 1 Cd isotope (e.g. only Cd110 or only Cd111) in the Mo_interference.csv file.")}

if(!("index" %in% colnames(mi))) {stop("ERROR: The column _index_ is missing in the file Mo_interference.csv or you wrote index with a capital I instead of i.")}
if(!("sample" %in% colnames(mi))) {stop("ERROR: The column _sample_ is missing in the file Mo_interference.csv or is written incorrectly / don't use capital letters")}
if(!("index" %in% colnames(rh))) {stop("ERROR: The column _index_ is missing in one of the input files or you wrote index with a capital I instead of i in file Rh_standard.csv.")}


## Correct Mo and Cd for Rh

vec.mi <- names(mi[3:dim(mi)[2]])

mi.rh <- mi[,c("index","sample")]

	for(i in 1:length(vec.mi)){
#i=1

			# select data
			st1 <- mi[,c("index","sample",vec.mi[i])]
			nm <- str_sub(vec.mi[i],1,2)	#library(stringr)

			# to remove any row's with NA's that excel may initiate
				z <- is.na(st1[,c("index")])
				st1 <- st1[z==FALSE,]

			if(Rh.correction == "yes"){
				# selection of LR or MR Rh standard data
				ne <- str_sub(vec.mi[i],-3,-2)

				if(exists("ne") == FALSE) {print("ERROR: you are working in an older version of R that does not support the stringr library!")}				
				if(ne == "") {print("ERROR: you are working in an older version of R that does not support the stringr library!")}				

				if(nm == paste(standard.element,sep="")) {next}				

				if(ne == "LR") {Rh.corr.a <- rh[grepl("LR", names(rh))]; index <- rh$index; Rh.corr <- cbind(index,Rh.corr.a); name.1 <- names(Rh.corr)[2]}
				if(ne == "MR") {Rh.corr.a <- rh[grepl("MR", names(rh))]; index <- rh$index; Rh.corr <- cbind(index,Rh.corr.a); name.1 <- names(Rh.corr)[2]} 

				st <- merge(st1, Rh.corr, by = "index", all.y=FALSE)
						
		# Correct data for Rh standard
		if(ne == "LR") {st$Me.rh.corrected <- (st[,c(vec.mi[i])]*st[,c(paste("mean.",standard.element,".range.LR",sep=""))])/st[,c(name.1)]}		
		if(ne == "MR") {st$Me.rh.corrected <- (st[,c(vec.mi[i])]*st[,c(paste("mean.",standard.element,".range.MR",sep=""))])/st[,c(name.1)]}
				} else {st <- st1}


			if(Rh.correction == "yes"){

					temp <- st[,c("index", "Me.rh.corrected")]
					names(temp)[names(temp) == "Me.rh.corrected"] <- paste(vec.mi[i],standard.element,".corr",sep="")
					mi.rh <- merge(mi.rh, temp, by = "index", all.y=FALSE)			
								}

			if(Rh.correction == "no"){mi.noRh <- mi}


} # end loop i

# make figure Mo interference
			name.fig.Mo <- paste("Figure Mo interference",name.run,sep="_")
			tiff(file = paste(folder.Mo, name.fig.Mo, ".tiff", sep=""), width = 5, height = 5, units = "in", res=resolution.figures)


	# Define vec.cd	
	if(Rh.correction == "yes") {
	if(paste("Cd111.LR.",standard.element,".corr",sep="") %in% names(mi.rh)) {vec.cd <- paste("Cd111.LR.",standard.element,".corr",sep=""); vec.lim.b <- mi.rh[,paste("Cd111.LR.",standard.element,".corr",sep="")]}
	if(paste("Cd110.LR.",standard.element,".corr",sep="") %in% names(mi.rh)) {vec.cd <- paste("Cd110.LR.",standard.element,".corr",sep=""); vec.lim.b <- mi.rh[,paste("Cd110.LR.",standard.element,".corr",sep="")]}
	vec.lim.c <- mi.rh[,paste("Mo95.LR.",standard.element,".corr",sep="")]
	} else 	{  
		if("Cd111.LR." %in% names(mi)) {vec.cd <- "Cd111.LR."}		# code 060616_1
		if("Cd110.LR." %in% names(mi)) {vec.cd <- "Cd110.LR."}		# code 060616_1
		vec.lim.b <- NA
	 	vec.lim.c <- NA
			}

	# define the axis scales for Cd
	if("Cd111.LR." %in% names(mi)) {vec.lim.a <- mi[,"Cd111.LR."]; nb <- "Cd111.LR."}
	if("Cd110.LR." %in% names(mi)) {vec.lim.a <- mi[,"Cd110.LR."]; nb <- "Cd110.LR."}

	vec.lim.cd <- c(vec.lim.a,vec.lim.b)

	x.a <- min(vec.lim.cd, na.rm=TRUE)
	x.b <- max(vec.lim.cd, na.rm=TRUE)

	# define the axis scales for Mo
	vec.lim.d <- mi[,"Mo95.LR."]

	vec.lim.mo <- c(vec.lim.c,vec.lim.d)

	y.a <- min(vec.lim.mo, na.rm=TRUE)
	y.b <- max(vec.lim.mo, na.rm=TRUE)


			# select for Cd only				
			nm <- str_sub(vec.cd,1,2)	#library(stringr)
			
			#windows(7,7)
			par(mfrow=c(1,1))
			par(mar=c(5, 6, 3, 2) + 0.1)		# c(bottom, left, top, right)
			
			plot(vec.lim.d ~ vec.lim.a, las = 1, pch=16, xlab="Cd (signal)", ylab="", 
											cex=cex, cex.lab=cex.lab,cex.axis = cex.axis, col="black", main=name.run, xlim=c(x.a,x.b))
			title(ylab="Mo (signal)", cex.lab=cex.lab, line = 4.5)
			
			# Uncorrected for standard.element		
			if("Cd111.LR." %in% names(mi))  {
			points(mi[,"Mo95.LR."] ~ mi[,"Cd111.LR."],pch=1,cex=cex, cex.lab=cex.lab,cex.axis = cex.axis, col="black")
			modeluncor <- lm(mi[,c("Mo95.LR.")] ~ mi[,c("Cd111.LR.")])
			an <- "Cd111.LR." 
			}

			if("Cd110.LR." %in% names(mi))  {modeluncor <- lm(mi[,c("Mo95.LR.")] ~ mi[,c("Cd110.LR.")])
			points(mi[,"Mo95.LR."] ~ mi[,"Cd110.LR."],pch=1,cex=cex, cex.lab=cex.lab,cex.axis = cex.axis, col="black")
			modeluncor <- lm(mi[,c("Mo95.LR.")] ~ mi[,c("Cd110.LR.")])
			an <- "Cd110.LR."
			}

			abline(modeluncor, col="black") 

					slope.unc <- round(summary(modeluncor)[4][[1]][[2]],2)
					slope.StdError.unc <- round(summary(modeluncor)[4][[1]][[4]],2)
					intercept.unc <- round(summary(modeluncor)[4][[1]][[1]],2)
					intercept.StdError.unc <- round(summary(modeluncor)[4][[1]][[3]],2)
					Rsq.unc <- round(summary(modeluncor)[9][[1]],3)

			mtext("Mo interference of Cd signal", side=3, lin=0.1, adj=0.5, cex=0.8, col="black")
			mtext(paste(an," = ",slope.unc," x + ",intercept.unc,sep=""), side=3, lin=-1, adj=0.1, cex=0.8, col="black")
			mtext(paste("R^2 = ",Rsq.unc, sep=""), side=3, lin=-2, adj=0.07, cex=0.8, col="black")

			# Corrected for standard.element
		if(Rh.correction == "yes") {
			
			points(vec.lim.c ~ vec.lim.b, pch=1, col="red", cex=cex)

			model <- lm(vec.lim.c ~ vec.lim.b)
			abline(model, col="red", lty=3)
							slope <- round(summary(model)[4][[1]][[2]],2)
							slope.StdError <- round(summary(model)[4][[1]][[4]],2)
							intercept <- round(summary(model)[4][[1]][[1]],2)
							intercept.StdError <- round(summary(model)[4][[1]][[3]],2)
							Rsq <- round(summary(model)[9][[1]],3)
			
				mtext(paste(vec.cd," = ",slope," x + ",intercept,sep=""), side=3, lin=-3, adj=0.1, cex=0.8, col="red")
				mtext(paste("R^2 = ",Rsq, sep=""), side=3, lin=-4, adj=0.07, cex=0.8, col="red")
				
			legend("bottomright", legend = c(paste("uncorrected for ",standard.element,sep=""),paste("corrected for ",standard.element,sep="")), 
			col=c("black","red"),pch=c(16,1), lty = c(1,3), cex=1, bty="n")
			} else {
			legend("bottomright", legend = c(paste("uncorrected for ",standard.element,sep="")), 
			col=c("black"),pch=c(16), lty = c(1), cex=1, bty="n")}

dev.off()

# Bring data together in table
# Create table for data head(data.standadd)

			data.standadd <- as.data.frame(matrix(NA,ncol=3,nrow=5))
			colnames(data.standadd) <- c("parameter",nb,paste(nb,standard.element,".corr",sep=""))
			
			data.standadd[1,1] <- "slope"
			data.standadd[2,1] <- "slope.StdError"
			data.standadd[3,1] <- "intercept"
			data.standadd[4,1] <- "intercept.StdError"
			data.standadd[5,1] <- "R^2"

			data.standadd[1,2] <- slope.unc
			data.standadd[2,2] <- slope.StdError.unc
			data.standadd[3,2] <- intercept.unc
			data.standadd[4,2] <- intercept.StdError.unc
			data.standadd[5,2] <- Rsq.unc
	
if(Rh.correction == "yes") 	{
		
			data.standadd[1,3] <- slope
			data.standadd[2,3] <- slope.StdError
			data.standadd[3,3] <- intercept
			data.standadd[4,3] <- intercept.StdError
			data.standadd[5,3] <- Rsq
					}

write.table(data.standadd,file=paste(folder.Mo,"Mo_interferenceCd_results_",name.run,".csv", sep=""),sep=",",row.names=F)

} # end of function

###########################################################################################
######### Blank determination #############################################################

fBlank.Conc.Stats <- function(blnk, sa.sf, st.el, mor, name.run, mainDir,resolution.figures, cex, cex.lab, cex.axis, Standard.element.correction,Mo.correction.Cd,...){

options(warn=-1) # to terminate all warning messages

Rh.correction <- Standard.element.correction 
rh <- st.el

if(Rh.correction == "yes"){
				tst2 <- sa.sf[6,]
				tst2a <- cbind(tst2[2:dim(tst2)[2]])
				if(is.na(rowMeans(tst2a, na.rm=TRUE)[1][[1]]) == TRUE) {stop(paste("You need to set Standard.element.correction at no as you did not ",standard.element," correct your SEAFAST standard additions.", sep=""))}
				}

if(!("index" %in% colnames(blnk))) {stop("ERROR: The column _index_ is missing in the file blanks.csv or you wrote index with a capital I instead of i.")}
if(!("sample" %in% colnames(blnk))) {stop("ERROR: The column _sample_ is missing in the file blanks.csv or was written incorrectly / don't use capital letters")}
if(!("blank" %in% colnames(blnk))) {stop("ERROR: The column _blank_ is missing in the file blanks.csv or was written incorrectly / don't use capital letters")}
if(!("index" %in% colnames(rh))) {stop("ERROR: The column _index_ is missing in one of the input files or you wrote index with a capital I instead of i in file Rh_standard.csv.")}

if(Rh.correction == "yes" & Mo.correction.Cd == "no") {print(paste("All data are ",standard.element," corrected. Cd has not been corrected for Mo interference.", sep=""))} 
if(Rh.correction == "yes" & Mo.correction.Cd == "yes") {print(paste("All data are ",standard.element," corrected. Cd has been corrected for Mo interference. Cd and Mo were ",standard.element," corrected.", sep=""))} 
if(Rh.correction == "no" & Mo.correction.Cd == "yes") {print(paste("All data are not ",standard.element," corrected. Cd has been corrected for Mo interference. Cd and Mo were not ",standard.element," corrected.", sep=""))} 
if(Rh.correction == "no" & Mo.correction.Cd == "no") {print(paste("All data are not ",standard.element," corrected. Cd has not been corrected for Mo interference.", sep=""))} 

vec.bl <- names(blnk)[4:dim(blnk)[2]]
vec.type <- unique(blnk$blank)
bl <- blnk

# Create new file folder for figures and data of the blanks
	subDir.fig.blank <- paste("Figures & data blank",name.run,sep="_")
	dir.create(file.path(mainDir, subDir.fig.blank))			
	folder.blank <<- paste(mainDir,subDir.fig.blank,"/",sep="")

	for(i in 1:length(vec.bl)){
#i=3

				# select data & safe original data frame
				st1 <- bl[,c("index","sample","blank",vec.bl[i])]
				nm <- str_sub(vec.bl[i],1,2)	#library(stringr)
				if(nm == "Yb"){nm <- "Yb"} else {				
				nm <- str_sub(vec.bl[i],1,1)
				if(nm == "V" | nm == "Y" | nm == "U"){} else { 
				nm <- str_sub(vec.bl[i],1,2)	#library(stringr)
				}
				}

				# to remove any row's with NA's that excel may initiate
				z <- is.na(st1[,c("index")])
				st1 <- st1[z==FALSE,]

			if(Rh.correction == "yes"){
				# selection of LR or MR Rh standard data
				ne <- str_sub(vec.bl[i],-3,-2)

				if(exists("ne") == FALSE) {print("ERROR: you are working in an older version of R that does not support the stringr library!")}				
				if(ne == "") {print("ERROR: you are working in an older version of R that does not support the stringr library!")}				

				if(nm == paste(standard.element,sep="")) {next}				

				if(ne == "LR") {Rh.corr.a <- rh[grepl("LR", names(rh))]; index <- rh$index; Rh.corr <- cbind(index,Rh.corr.a); name.1 <- names(Rh.corr)[2]}
				if(ne == "MR") {Rh.corr.a <- rh[grepl("MR", names(rh))]; index <- rh$index; Rh.corr <- cbind(index,Rh.corr.a); name.1 <- names(Rh.corr)[2]} 

				st <- merge(st1, Rh.corr, by = "index", all.y=FALSE)
						
		# Correct data for Rh standard
		if(ne == "LR") {st$Me.rh.corrected <- (st[,c(vec.bl[i])]*st[,c(paste("mean.",standard.element,".range.LR",sep=""))])/st[,c(name.1)]}		
		if(ne == "MR") {st$Me.rh.corrected <- (st[,c(vec.bl[i])]*st[,c(paste("mean.",standard.element,".range.MR",sep=""))])/st[,c(name.1)]}
				} else {st <- st1}
		
				# The metal is skipped if no slope was determined
				if(is.na(sa.sf[1,c(vec.bl[i])]) == TRUE) {next}

				# To correct Cd signal for Mo interference
				if(nm == "Cd" & Mo.correction.Cd == "yes" & Rh.correction == "yes"){
					if(Rh.correction == "yes" & is.na(mor[1,3]) == TRUE){stop("ERROR: you have run fMo.interference.Cd without standard.element correction. 
You have 3 options: 	
i) run fMo.interference.Cd with correction,
ii) run fBlank.Conc.Stats with Standard.element.correction <- no
iii) run fBlank.Conc.Stats with Mo.correction.Cd <- no")}

									# Select the correct slope
									if(vec.bl[i] == "Cd110.LR.") {slopeMo <- mor[1,paste("Cd110.LR.",standard.element,".corr", sep="")]}	
									if(vec.bl[i] == "Cd111.LR.") {slopeMo <- mor[1,paste("Cd111.LR.",standard.element,".corr", sep="")]}			
									
									if(is.null(slopeMo) == TRUE) {
										print("ERROR: You used _Standard.element.correction <- no_ with the function _fMo.interference.Cd_.")}

									# Bring Mo into st
									MoIn <- bl[,c("index","Mo95.LR.")]
									st <- merge(st, MoIn, by = "index", all.y=FALSE)
									st$Me.rh.corrected <- st$Me.rh.corrected - (st$Mo95.LR./slopeMo)
									st <- cbind(st[match("index",names(st)):match("Me.rh.corrected",names(st))])
									print(paste(vec.bl[i]," blank was corrected for ",standard.element," and Mo interference",sep=""))
											}

				if(nm == "Cd" & Mo.correction.Cd == "yes" & Rh.correction == "no"){

									# Select the correct slope
									if(vec.bl[i] == "Cd110.LR.") {slopeMo <- mor[1,"Cd110.LR."]}	
									if(vec.bl[i] == "Cd111.LR.") {slopeMo <- mor[1,"Cd111.LR."]}			
									
									if(is.null(slopeMo) == TRUE) {
										print("ERROR: You have not used the function _fMo.interference.Cd_.")}

									# Bring Mo into st
									MoIn <- bl[,c("index","Mo95.LR.")]
									st <- merge(st, MoIn, by = "index", all.y=FALSE)
									st[,vec.bl[i]] <- st[,vec.bl[i]] - (st$Mo95.LR./slopeMo)
									st <- cbind(st[match("index",names(st)):match(vec.bl[i],names(st))])
									#print(paste(vec.bl[i]," blank was not corrected for ",standard.element," but was for the (",standard.element," uncorrected) Mo interference",sep=""))
											}

				# Calculate concentrations

				if(Rh.correction == "yes"){
					st$slope.rh.corr <- sa.sf[6,c(vec.bl[i])]
					st$slope.rh.uncorr <- sa.sf[1,c(vec.bl[i])]
				
					st$Me.conc.rh.corr <- st$Me.rh.corrected/st$slope.rh.corr
					st$Me.conc.rh.uncorr <- st[,c(vec.bl[i])]/st$slope.rh.uncorr

				} else {
					
					st$slope.rh.uncorr <- sa.sf[1,c(vec.bl[i])]
					st$Me.conc.rh.uncorr <- st[,c(vec.bl[i])]/st$slope.rh.uncorr

				}
						
if(Rh.correction == "yes"){

					temp <- st[,c("index", "Me.conc.rh.corr")]
					names(temp)[names(temp) == "Me.conc.rh.corr"] <- paste(vec.bl[i],"conc",sep="")
					bl <- merge(bl, temp, by = "index", all.y=FALSE)
							
				} else {

					temp <- st[,c("index","Me.conc.rh.uncorr")]
					names(temp)[names(temp) == "Me.conc.rh.uncorr"] <- paste(vec.bl[i],"conc",sep="")
					bl <- merge(bl, temp, by = "index", all.y=FALSE)
				
				}}

write.table(bl,file=paste(folder.blank,"Blank concentrations_",name.run,".csv", sep=""),sep=",",row.names=F)

#if(Rh.correction == "yes"){print(paste("Blank metal concentrations are calculated: ",standard.element," corrected", sep=""))} else {
#				   print(paste("Blank metal concentrations are calculated: not ",standard.element," corrected", sep=""))} 


# Create separate files with the values for the means, medians, standard errors etc. of the different blanks

# Create vecor with metal blank concentrations head(bl)
vec.conc <- names(bl)[(dim(blnk)[2]+1):dim(bl)[2]]

#data.blank <- as.data.frame(matrix(NA,ncol=(length(vec.conc)+1),nrow=5*length(unique(bl$blank))))
#colnames(data.blank) <- c("parameter",vec.conc)


for(j in 1:length(vec.type)){
#j=1					

		data.blank <- as.data.frame(matrix(NA,ncol=(length(vec.conc)+1),nrow=5))
		colnames(data.blank) <- c("parameter",vec.conc)
		m <- seq(1, dim(data.blank)[1], by =5)

		bld3 <- subset(bl, blank == vec.type[j])
						
						z <- m[j]

						data.blank[1,1] 	<- paste("1st_Qu.",vec.type[j],sep="_")
						data.blank[2,1] <- paste("median",vec.type[j],sep="_")
						data.blank[3,1] <- paste("3rd_Qu.",vec.type[j],sep="_")
						data.blank[4,1] <- paste("mean",vec.type[j],sep="_")
						data.blank[5,1] <- paste("StDev",vec.type[j],sep="_")
						data.blank[6,1] <- paste("RSD%",vec.type[j],sep="_")

							
			for(k in 1:length(vec.conc)+1) 	{
#k=5
							
							res1 <- summary(bld3[,c(vec.conc[k-1])])
							
							data.blank[1,k] <- res1[2][[1]]
							data.blank[2,k] <- res1[3][[1]]
							data.blank[3,k] <- res1[5][[1]]
							data.blank[4,k] <- res1[4][[1]]
							data.blank[5,k] <- sd(bld3[,c(vec.conc[k-1])], na.rm=TRUE)
							data.blank[6,k] <- (100*sd(bld3[,c(vec.conc[k-1])], na.rm=TRUE))/res1[4][[1]]
							if(summary(bld3[,c(vec.conc[k-1])])[2][[1]] == "try-error") next

									} # k loop

write.table(data.blank,file=paste(folder.blank,"StatsData_blank_",vec.type[j],"_",name.run,".csv", sep=""),sep=",",row.names=F)

} # j loop

print("Separate data files with concentrations of the different blanks have been created.")

library(plyr)

#################
# The function needed to label the outliers in a boxplot was taken from the internet
#https://raw.githubusercontent.com/talgalili/R-code-snippets/master/boxplot.with.outlier.label.r

# last updated: 31.10.2011
# 		This is instead of the 20.6.11 version...

boxplot.with.outlier.label <- function(y, label_name, ..., spread_text = T, data, plot = T, range = 1.5, label.col = "blue", push_text_right = 1.3, # enlarge push_text_right in order to push the text labels further from their point
segement_width_as_percent_of_label_dist = .45, # Change this if you want to have the line closer to the label (range should be between 0 to 1
	jitter_if_duplicate = T, jitter_only_positive_duplicates = F)
{	
	# notes - this functions doesn't work if there are any missing values in the data.
	#		You must pre-process the data to make sure it is "complete".


	# change log:
	# 19.04.2011 - added support to "names" and "at" parameters.


	# jitter_if_duplicate - will jitter (Actually just add a bit of numbers) so to be able to decide on which location to plot the label when having identical variables...
	require(plyr) # for is.formula and ddply

	# a function to jitter data in case of ties in Y's
	jitter.duplicate <- function(x, only_positive = F)
	{
		if(only_positive) {
			ss <- x > 0
		} else {
			ss <- T
		}	
		ss_dup <- duplicated(x[ss])
		# ss <- ss & ss_dup
		temp_length <- length(x[ss][ss_dup])	
		x[ss][ss_dup] <- x[ss][ss_dup] + seq(from = 0.00001, to = 0.00002, length.out = temp_length)
		x
	}
		
	# handle cases where 
	if(jitter_if_duplicate) {
		# warning("duplicate jutter of values in y is ON")
		if(!missing(data)) {	#e.g: we DO have data
			# if(exists("y") && is.formula(y)) {		# F && NULL # F & NULL
			y_name <- as.character(substitute(y))	# I could have also used as.list(match.call())
												# credit to Uwe Ligges and Marc Schwartz for the help
												# https://mail.google.com/mail/?shva=1#inbox/12dd7ca2f9bfbc39
			if(length(y_name) > 1) {	# then it is a formula (for example: "~", "y", "x"
				model_frame_y <- model.frame(y, data = data)
				temp_y <- model_frame_y[,1]
				temp_y  <- jitter.duplicate(temp_y, jitter_only_positive_duplicates)	# notice that the default of the function is to work only with positive values...
				# the_txt <- paste(names(model_frame_y)[1], "temp_y", sep = "<<-") # wrong...
				the_txt <- paste("data['",names(model_frame_y)[1],"'] <- temp_y", sep = "")				
				eval(parse(text = the_txt))	# jutter out y var so to be able to handle identical values.
			} else {	# this isn't a formula
				data[,y_name] <- jitter.duplicate(data[,y_name], jitter_only_positive_duplicates)
				y <- data[,y_name]	# this will make it possible for boxplot(y, data) to work later (since it is not supposed to work with data when it's not a formula, but now it does :))
			}		
		} else {	# there is no "data"		 
			if(is.formula(y)) { # if(exists("y") && is.formula(y)) {		# F && NULL # F & NULL
				temp_y <- model.frame(y)[,1]
				temp_y  <- jitter.duplicate(temp_y, jitter_only_positive_duplicates)	# notice that the default of the function is to work only with positive values...
				temp_y_name <- names(model.frame(y))[1]	# we must extract the "names" before introducing a new enbironment (or there will be an error)
				environment(y) <- new.env()
				assign(temp_y_name, temp_y, environment(y))
					# Credit and thanks for doing this goes to Niels Richard Hansen (2 Jan 30, 2011)
					# http://r.789695.n4.nabble.com/environment-question-changing-variables-from-a-formula-through-model-frame-td3246608.html
				# warning("Your original variable (in the global environemnt) was just jittered.")	# maybe I should add a user input before doing this....
				# the_txt <- paste(names(model_frame_y)[1], "temp_y", sep = "<<-")
				# eval(parse(text = the_txt))	# jutter out y var so to be able to handle identical values.
			} else {
				y <- jitter.duplicate(y, jitter_only_positive_duplicates)
			}		
		}
	}
	
	
	# y should be a formula of the type: y~x, y~a*b
	# or it could be simply y
	if(missing(data)) {
			boxdata <- boxplot(y, plot = plot,range = range ,...)
		} else {
			boxdata <- boxplot(y, plot = plot,data = data, range = range ,...)
		}
	if(length(boxdata$names) == 1 && boxdata$names =="") boxdata$names <- 1	# this is for cases of type: boxplot(y) (when there is no dependent group)
	if(length(boxdata$out) == 0 ) {
		warning("No outliers detected for this boxplot")
		return(invisible())
		}
	
	if(!missing(data)) attach(data)	# this might lead to problams I should check out for alternatives for using attach here...
	

	# creating a data.frame with information from the boxplot output about the outliers (location and group)
	boxdata_group_name <- factor(boxdata$group)
	levels(boxdata_group_name) <- boxdata$names[as.numeric(levels(boxdata_group_name))]	# the subseting is for cases where we have some sub groups with no outliers
	if(!is.null(list(...)$at))	{	# if the user chose to use the "at" parameter, then we would like the function to still function (added on 19.04.2011)
		boxdata$group <- list(...)$at[boxdata$group]		
		}
	boxdata_outlier_df <- data.frame(group = boxdata_group_name, y = boxdata$out, x = boxdata$group)
	

	# Let's extract the x,y variables from the formula:
	if(is.formula(y))
	{
		model_frame_y <- model.frame(y)
			# old solution: (which caused problems if we used the names parameter when using a 2 way formula... (since the order of the names is different then the levels order we get from using factor)
			# y <- model_frame_y[,1]
			# x <- model_frame_y[,-1]

		y <- model_frame_y[,1]
		x <- model_frame_y[,-1]
		if(!is.null(dim(x))) {	# then x is a matrix/data.frame of the type x1*x2*..and so on - and we should merge all the variations...
			x <- apply(x,1, paste, collapse = ".")
		}
	} else {
		# if(missing(x)) x <- rep(1, length(y))
		x <- rep(1, length(y))	# we do this in case y comes as a vector and without x
	}	
	
	# and put all the variables (x, y, and outlier label name) into one data.frame
	DATA <- data.frame(label_name, x ,y)
	
	if(!is.null(list(...)$names))	{	# if the user chose to use the names parameter, then we would like the function to still function (added on 19.04.2011)
		DATA$x <- factor(DATA$x, levels = unique(DATA$x))
		levels(DATA$x) = list(...)$names	# enable us to handle when the user adds the "names" parameter # fixed on 19.04.11	# notice that DATA$x must be of the "correct" order (that's why I used split above
		# warning("Careful, the use of the 'names' parameter is experimental.  If you notice any errors please e-mail me at: tal.galili@gmail.com")
		}

	if(!missing(data)) detach(data)	# we don't need to have "data" attached anymore.

	# let's only keep the rows with our outliers 
	boxplot.outlier.data <- function(xx, y_name = "y")
	{
		y <- xx[,y_name]
		boxplot_range <- range(boxplot.stats(y, coef = range )$stats)
		ss <- (y < boxplot_range[1]) | (y > boxplot_range[2])
		return(xx[ss,])	
	}
	outlier_df <-ddply(DATA, .(x), boxplot.outlier.data)
	

	# create propor x/y locations to handle over-laping dots...
	if(spread_text) {
		# credit: Greg Snow
		require(TeachingDemos)		
		temp_x <- boxdata_outlier_df[,"x"]
		temp_y1 <- boxdata_outlier_df[,"y"]
		temp_y2 <- temp_y1
		for(i in unique(temp_x))
		{
			tmp <- temp_x == i
			temp_y2[ tmp ] <- spread.labs( temp_y2[ tmp ], 1.3*strheight('A'), maxiter=6000, stepsize = 0.05) #, min=0 )
		}
		
	}
	
	# plotting the outlier labels :)  (I wish there was a non-loop wise way for doing this)
	for(i in seq_len(dim(boxdata_outlier_df)[1]))
	{
		# ss <- (outlier_df[,"x"]  %in% boxdata_outlier_df[i,]$group) & (outlier_df[,"y"] %in% boxdata_outlier_df[i,]$y)

		# if(jitter_if_duplicate) {
			# ss <- (outlier_df[,"x"]  %in% boxdata_outlier_df[i,]$group) & closest.number(outlier_df[,"y"]  boxdata_outlier_df[i,]$y)
		# } else {
		ss <- (outlier_df[,"x"]  %in% boxdata_outlier_df[i,]$group) & (outlier_df[,"y"] %in% boxdata_outlier_df[i,]$y)
		# }

		current_label <- outlier_df[ss,"label_name"]
		temp_x <- boxdata_outlier_df[i,"x"]
		temp_y <- boxdata_outlier_df[i,"y"]		
		# cbind(boxdata_outlier_df,		temp_y2)
		# outlier_df

		
		
		if(spread_text) {
			temp_y_new <- temp_y2[i] # not ss			
			move_text_right <- strwidth(current_label) * push_text_right
			text( temp_x+move_text_right, temp_y_new, current_label, col = label.col)			
			# strwidth
			segments( temp_x+(move_text_right/6), temp_y, temp_x+(move_text_right*segement_width_as_percent_of_label_dist), temp_y_new )
		} else {
			text(temp_x, temp_y, current_label, pos = 4, col = label.col)
		}		
	}

	# outputing some of the information we collected
	invisible(list(boxdata = boxdata, boxdata_outlier_df = boxdata_outlier_df, outlier_df=outlier_df))
}

######################

# Create boxplots for the blanks of the metals

# create directory for boxplots
				subDir.boxplots.blank <- paste("Boxplots_blanks",name.run,sep="_")
				dir.create(file.path(mainDir, subDir.fig.blank,subDir.boxplots.blank))

for(k in 1:length(vec.conc)){
# k=1
				# select the data
				bld2 <- as.data.frame(bl[,c("index", "sample", "blank", vec.conc[k])])
				nm <- str_sub(vec.conc[k],1,2)	#library(stringr)
				if(nm == "Yb"){nm <- "Yb"} else {				
				nm <- str_sub(vec.conc[k],1,1)
				if(nm == "V" | nm == "Y" | nm == "U"){} else { 
				nm <- str_sub(vec.conc[k],1,2)	#library(stringr)
				}
				}

				# remove rows with NA's
				z <- is.na(bld2[,c(vec.conc[k])])
				bld2 <- bld2[z==FALSE,]
				
				# The metal is skipped if no concentrations at all were determined
				nod <- mean(bld2[,c(vec.conc[k])],na.rm=TRUE)
				if(is.na(nod) == TRUE) {next}

				# units for the different metals
				if(nm %in% units.nmolKg) {units <- "nmol/kg"}
				if(nm %in% units.pmolKg) {units <- "pmol/kg"}

				# Make figure
				name.fig.blank <- paste(vec.conc[k],name.run,sep="_")
				tiff(file = paste(folder.blank, subDir.boxplots.blank,"/", name.fig.blank, ".tiff", sep=""), width = 5, height = 6, units = "in", res=resolution.figures)
	
				#windows(5,6)
				par(mfrow=c(1,1))
				par(mar=c(10, 5, 3, 2) + 0.1)		# c(bottom, left, top, right)

				lab_y <- bld2$index
				boxplot.with.outlier.label(bld2[,4] ~ bld2[,"blank"], lab_y, las = 2, ylab="", xlab="", 
						cex=cex, cex.lab=cex.lab,cex.axis = cex.axis, main=name.run,
						label.col = "blue",
						push_text_right = 1.5,
						segement_width_as_percent_of_label_dist = .4)
						title(ylab=paste("Blank conc. of ",nm," (",units,")",sep=""), cex.lab=cex.lab, line = 4)
						means <- tapply(bld2[,4],bld2[,3],function(x) mean(x,na.rm=TRUE))
						points(means,col="red",pch=18,cex=cex)
						mtext("index numbers belonging to the outliers", side=3, lin=0.1, adj=0.0, cex=0.8, col="blue")
						legend("topright", legend = c("median","mean"), 
								col=c("black","red"),pch=c(NA,18), lty = c(1,NA), lwd=c(2,NA), cex=1, bty="n")
						
						if(Rh.correction == "no") {
						mtext(paste("not ",standard.element," corrected", sep=""), side=3, lin=0.1, adj=0.98, cex=0.8, col="black")
							} else {
						mtext(paste(standard.element," corrected", sep=""), side=3, lin=0.1, adj=0.98, cex=0.8, col="black")}
						
						if(Mo.correction.Cd == "no" & nm == "Cd") {
						mtext("not Mo corrected", side=3, lin=1, adj=0.98, cex=0.8, col="black")}
											 
						if(Mo.correction.Cd == "yes" & nm == "Cd") {
						mtext("Mo corrected", side=3, lin=1, adj=0.98, cex=0.8, col="black")}

						if(length(unique(bld2[,3])) == 1) {title(xlab=paste(unique(bld2[,3])), cex.lab=1, line = 1)}



				dev.off()

				graphics.off()

				rm(bld2)

} # end loop k

print("Boxplots of the blanks have been created.")

# Create dot plots for the blanks of the metals

# create directory for dotplots
				subDir.dotplots.blank <- paste("Dotplots_blanks",name.run,sep="_")
				dir.create(file.path(mainDir, subDir.fig.blank,subDir.dotplots.blank))

for(k in 1:length(vec.conc)){
# k=1
				# select the data
				bld2 <- as.data.frame(bl[,c("index", "sample", "blank", vec.conc[k])])
				nm <- str_sub(vec.conc[k],1,2)	#library(stringr)
				if(nm == "Yb"){nm <- "Yb"} else {				
				nm <- str_sub(vec.conc[k],1,1)
				if(nm == "V" | nm == "Y" | nm == "U"){} else { 
				nm <- str_sub(vec.conc[k],1,2)	#library(stringr)
				}
				}

				# remove rows with NA's
				z <- is.na(bld2[,c(vec.conc[k])])
				bld2 <- bld2[z==FALSE,]
				
				# The metal is skipped if no concentrations at all were determined
				nod <- mean(bld2[,c(vec.conc[k])],na.rm=TRUE)
				if(is.na(nod) == TRUE) {next}

				# add arbitrary index
				bld2$arbitrary.index <- seq(1,dim(bld2)[1],1)

				# units for the different metals
				if(nm %in% units.nmolKg) {units <- "nmol/kg"}
				if(nm %in% units.pmolKg) {units <- "pmol/kg"}

				# Make figure
				name.fig.blank <- paste(vec.conc[k],name.run,sep="_")
				tiff(file = paste(folder.blank, subDir.dotplots.blank,"/", name.fig.blank, ".tiff", sep=""), width = 6, height = 4, units = "in", res=resolution.figures)
	
#windows(6,4)
						par(mfrow=c(1,1))
						par(mar=c(5, 5, 3, 2) + 0.1)		# c(bottom, left, top, right)

						plot(bld2[,4] ~ bld2[,"arbitrary.index"],las = 1, ylab="", xlab="arbitrary index", pch=16,
						cex=cex, cex.lab=cex.lab,cex.axis = cex.axis, main=name.run, xlim=c(0,(dim(bld2)[1]+2)))
						title(ylab=paste(nm," (",units,")",sep=""), cex.lab=cex.lab, line = 4)
						text(bld2$arbitrary.index,bld2[,4], label=bld2$index, pos=2, col="blue", cex=0.7)
						text(bld2$arbitrary.index,bld2[,4], label=bld2$blank, pos=4, col="red",cex=0.7) 
						mtext("data index",side=3, lin=1, adj=0.02, cex=0.8, col="blue")
						mtext("Blank",side=3, lin=-0.1, adj=0.5, cex=1, col="black")

						if(Rh.correction == "no") {
						mtext(paste("not ",standard.element," corrected", sep=""), side=3, lin=0.1, adj=0.02, cex=0.8, col="black")
							} else {
						mtext(paste(standard.element," corrected", sep=""), side=3, lin=0.1, adj=0.02, cex=0.8, col="black")}
						
						if(Mo.correction.Cd == "no" & nm == "Cd") {
						mtext("not Mo corrected", side=3, lin=0.1, adj=0.98, cex=0.8, col="black")}
											 
						if(Mo.correction.Cd == "yes" & nm == "Cd") {
						mtext("Mo corrected", side=3, lin=0.1, adj=0.98, cex=0.8, col="black")}


				dev.off()

				graphics.off()

				rm(bld2)

				

} # end loop k

print("Dotplots of the blanks have been created.")

} # end function

####################################################################################################################################

fChecks.balances <- function(ck, sa.sf, sa.HNO3, st.el, blank, ref.mat, mor, coding, sample.naming, name.run, mainDir, use.Blank, Seafast.conc.factor, outliers.boxplot,resolution.figures, cex, cex.lab, cex.axis, Standard.element.correction,Mo.correction.Cd,...){

options(warn=-1) # to terminate all warning messages

Rh.correction <- Standard.element.correction
rh <- st.el

# code 070616_1
if(Rh.correction == "yes"){
				tst2 <- sa.sf[6,]
				tst2a <- cbind(tst2[2:dim(tst2)[2]])
				if(is.na(rowMeans(tst2a, na.rm=TRUE)[1][[1]]) == TRUE) {stop(paste("You need to set Standard.element.correction at no as you did not ",standard.element," correct your SEAFAST standard additions.", sep=""))}
				}

if(Rh.correction == "yes"){
				tst3 <- sa.HNO3[6,]
				tst3a <- cbind(tst3[2:dim(tst3)[2]])
				if(is.na(rowMeans(tst3a, na.rm=TRUE)[1][[1]]) == TRUE) {
				print(paste("WARNING: You you did not ",standard.element," correct your ICPMS standard additions therefore no recoveries can be calculated.", sep="")) 
				print("Change Standard.element.correction into no.")}
				}

if(!("index" %in% colnames(ck))) {stop("ERROR: The column _index_ is missing in the file checks.csv or you wrote index with a capital I instead of i.")}
if(!("sample" %in% colnames(ck))) {stop("ERROR: The column _sample_ is missing in the file checks.csv or was written incorrectly / don't use capital letters")}
if(!("check" %in% colnames(ck))) {stop("ERROR: The column _check_ is missing in the file checks.csv or was written incorrectly / don't use capital letters")}
if(!("reference" %in% colnames(ck))) {stop("ERROR: The column _reference_ is missing in the file checks.csv or you wrote index with a capital I instead of i.")}
if(!("where" %in% colnames(ck))) {stop("ERROR: The column _where_ is missing in the file checks.csv or you wrote index with a capital I instead of i.")}
if(!("index" %in% colnames(rh))) {stop("ERROR: The column _index_ is missing in one of the input files or you wrote index with a capital I instead of i in file Rh_standard.csv.")}
if(!("seafast.name" %in% colnames(coding))) {stop("ERROR: The column _seafast.name_ is missing in the file sample_list.csv or was written incorrectly / don't use capital letters")}
if(!("sample.name" %in% colnames(coding))) {stop("ERROR: The column _sample.name_ is missing in the file sample_list.csv or was written incorrectly / don't use capital letters")}

if(Rh.correction == "yes" & Mo.correction.Cd == "no") {print(paste("All data are ",standard.element," corrected. Cd has not been corrected for Mo interference.", sep=""))} 
if(Rh.correction == "yes" & Mo.correction.Cd == "yes") {print(paste("All data are ",standard.element," corrected. Cd has been corrected for Mo interference. Cd and Mo were ",standard.element," corrected.", sep=""))} 
if(Rh.correction == "no" & Mo.correction.Cd == "yes") {print(paste("All data are not ",standard.element," corrected. Cd has been corrected for Mo interference. Cd and Mo were not ",standard.element," corrected.", sep=""))} 
if(Rh.correction == "no" & Mo.correction.Cd == "no") {print(paste("All data are not ",standard.element," corrected. Cd has not been corrected for Mo interference.", sep=""))} 

no.icpms.data <- 0

vec.ck <- names(ck)[6:dim(ck)[2]]
vec.ck.conc <- names(blank)[2:dim(blank)[2]]
vec.check.icpms <- unique(ck$check[ck$where == "icpms"])
vec.check.seafast <- unique(ck$check[ck$where == "seafast"])

vec.check.icpms.dotplot <- unique(ck$check[ck$where == "icpms" & is.na(ck$reference) == TRUE])
vec.check.seafast.dotplot <- unique(ck$check[ck$where == "seafast" & is.na(ck$reference) == TRUE])


cks <- ck

# Create new file folder for figures and data of the blanks
	subDir.fig.check <- paste("Figures_checks & balances",name.run,sep="_")
	dir.create(file.path(mainDir, subDir.fig.check))			
	folder.check <<- paste(mainDir,subDir.fig.check,"/",sep="")
	
	subDir.fig.ref <- "REFERENCE SAMPLES"
	dir.create(file.path(mainDir, subDir.fig.check, subDir.fig.ref))			
	folder.ref.met <- paste(mainDir,subDir.fig.check,"/",subDir.fig.ref,"/",sep="")

#### Calculation of concentrations for seafast checks

	for(i in 1:length(vec.ck)){
#i=11

				# select data & safe original data frame
				st1 <- cks[,c("index","sample","check","where",vec.ck[i])]
				st1 <- subset(st1, where == "seafast")
				nm <- str_sub(vec.ck[i],1,2)	#library(stringr)
				if(nm == "Yb"){nm <- "Yb"} else {				
				nm <- str_sub(vec.ck[i],1,1)
				if(nm == "V" | nm == "Y" | nm == "U"){} else { 
				nm <- str_sub(vec.ck[i],1,2)	#library(stringr)
				}
				}

				# remove rows with NA's
				z <- is.na(st1[,c("index")])
				st1 <- st1[z==FALSE,]
			
				if(Rh.correction == "yes"){
				# selection of LR or MR Rh standard data
				ne <- str_sub(vec.ck[i],-3,-2)

				if(exists("ne") == FALSE) {print("ERROR: you are working in an older version of R that does not support the stringr library!")}				
				if(ne == "") {print("ERROR: you are working in an older version of R that does not support the stringr library!")}				

				if(nm == paste(standard.element,sep="")) {next}				

				if(ne == "LR") {Rh.corr.a <- rh[grepl("LR", names(rh))]; index <- rh$index; Rh.corr <- cbind(index,Rh.corr.a); name.1 <- names(Rh.corr)[2]}
				if(ne == "MR") {Rh.corr.a <- rh[grepl("MR", names(rh))]; index <- rh$index; Rh.corr <- cbind(index,Rh.corr.a); name.1 <- names(Rh.corr)[2]} 

				st <- merge(st1, Rh.corr, by = "index", all.y=FALSE)
						
		# Correct data for Rh standard
		if(ne == "LR") {st$Me.rh.corrected <- (st[,c(vec.ck[i])]*st[,c(paste("mean.",standard.element,".range.LR",sep=""))])/st[,c(name.1)]}		
		if(ne == "MR") {st$Me.rh.corrected <- (st[,c(vec.ck[i])]*st[,c(paste("mean.",standard.element,".range.MR",sep=""))])/st[,c(name.1)]}
				} else {st <- st1}

				# The metal is skipped if no slope was determined
				if(is.na(sa.sf[1,c(vec.ck[i])]) == TRUE) {next}

				# To correct Cd signal for Mo interference
				if(nm == "Cd" & Mo.correction.Cd == "yes" & Rh.correction == "yes"){
						
				if(Rh.correction == "yes" & is.na(mor[1,3]) == TRUE){stop("ERROR: you have run fMo.interference.Cd without standard.element correction. 
You have 3 options: 	
i) run fMo.interference.Cd with correction,
ii) run fChecks.balances with Standard.element.correction <- no
iii) run fChecks.balances with Mo.correction.Cd <- no")}

									# Select the correct slope
									if(vec.ck[i] == "Cd110.LR.") {slopeMo <- mor[1,paste("Cd110.LR.",standard.element,".corr", sep="")]}	
									if(vec.ck[i] == "Cd111.LR.") {slopeMo <- mor[1,paste("Cd111.LR.",standard.element,".corr", sep="")]}			
									
									#if(is.null(slopeMo) == TRUE) {
									#	print("ERROR: you confused the identity of the Cd istopes as used for Mo interference and as measured by ICPMS.")
									#	print("ERROR: alternatively, you used _Standard.element.correction <- no_ with the function _fMo.interference.Cd_.")}

									# Bring Mo into st
									MoIn <- cks[,c("index","Mo95.LR.")]
									st <- merge(st, MoIn, by = "index", all.y=FALSE)
									st$Me.rh.corrected <- st$Me.rh.corrected - (st$Mo95.LR./slopeMo)
									st <- cbind(st[match("index",names(st)):match("Me.rh.corrected",names(st))])
									print(paste(vec.ck[i]," blank was corrected for ",standard.element," and Mo interference in seafast checks",sep=""))
											}

				if(nm == "Cd" & Mo.correction.Cd == "yes" & Rh.correction == "no"){

									# Select the correct slope
									if(vec.ck[i] == "Cd110.LR.") {slopeMo <- mor[1,"Cd110.LR."]}	
									if(vec.ck[i] == "Cd111.LR.") {slopeMo <- mor[1,"Cd111.LR."]}			
									
									if(is.null(slopeMo) == TRUE) {
										#print("ERROR: you confused the identity of the Cd istopes as used for Mo interference and as measured by ICPMS.")
										print("ERROR: you have not run the function _fMo.interference.Cd_.")}

									# Bring Mo into st
									MoIn <- cks[,c("index","Mo95.LR.")]
									st <- merge(st, MoIn, by = "index", all.y=FALSE)
									st[,vec.ck[i]] <- st[,vec.ck[i]] - (st$Mo95.LR./slopeMo)
									st <- cbind(st[match("index",names(st)):match(vec.ck[i],names(st))])
									print(paste(vec.ck[i]," blank was not corrected for ",standard.element," but was for the (",standard.element," uncorrected) Mo interference in seafast checks",sep=""))
											}

				# Calculate concentrations
				
				# determine the blank
				name <- paste(vec.ck[i],"conc", sep="")
				if(use.Blank == "mean") {seafast.blank <- blank[4,c(name)]}
				if(use.Blank == "median") {seafast.blank <- blank[2,c(name)]}
															

				if(Rh.correction == "yes"){
					st$slope.rh.corr <- sa.sf[6,c(vec.ck[i])]
					st$slope.rh.uncorr <- sa.sf[1,c(vec.ck[i])]

					st$seafast.blank <- seafast.blank
								
					st$Me.conc.rh.corr <- (st$Me.rh.corrected/st$slope.rh.corr)-seafast.blank
					st$Me.conc.rh.uncorr <- (st[,c(vec.ck[i])]/st$slope.rh.uncorr)-seafast.blank

				} else {
					
					st$slope.rh.uncorr <- sa.sf[1,c(vec.ck[i])]
					st$seafast.blank <- seafast.blank
					st$Me.conc.rh.uncorr <- (st[,c(vec.ck[i])]/st$slope.rh.uncorr)-seafast.blank
				}
						
if(Rh.correction == "yes"){

					temp <- st[,c("index", "Me.conc.rh.corr")]
					names(temp)[names(temp) == "Me.conc.rh.corr"] <- paste(vec.ck[i],"conc",sep="")
					cks <- merge(cks, temp, by = "index", all.y=FALSE)
				
				
				} else {

					temp <- st[,c("index","Me.conc.rh.uncorr")]
					names(temp)[names(temp) == "Me.conc.rh.uncorr"] <- paste(vec.ck[i],"conc",sep="")
					cks <- merge(cks, temp, by = "index", all.y=FALSE)
				
				}
} # end i loop

ckkk.seafast <- cbind(cks[1:5],cks[(length(vec.ck)+6):dim(cks)[2]])

write.table(ckkk.seafast,file=paste(folder.check,"Checks concentrations_seafast_",name.run,".csv", sep=""),sep=",",row.names=F)

print("Calculation of concentrations Seafast checks: completed")

#### Calculation of concentrations for ICPMS checks

cki <- ck

#####
	for(i in 1:length(vec.ck)){
#i=11

				# select data & safe original data frame
				st1 <- cki[,c("index","sample","check","where",vec.ck[i])]
				st1 <- subset(st1, where == "icpms")
				if(dim(st1)[1] == 0) {no.icpms.data <- 1; break}
				nm <- str_sub(vec.ck[i],1,2)	#library(stringr)
				if(nm == "Yb"){nm <- "Yb"} else {				
				nm <- str_sub(vec.ck[i],1,1)
				if(nm == "V" | nm == "Y" | nm == "U"){} else { 
				nm <- str_sub(vec.ck[i],1,2)	#library(stringr)
				}
				}

				# remove rows with NA's
				z <- is.na(st1[,c("index")])
				st1 <- st1[z==FALSE,]
			
				if(Rh.correction == "yes"){
				# selection of LR or MR Rh standard data
				ne <- str_sub(vec.ck[i],-3,-2)

				if(exists("ne") == FALSE) {print("ERROR: you are working in an older version of R that does not support the stringr library!")}				
				if(ne == "") {print("ERROR: you are working in an older version of R that does not support the stringr library!")}				

				if(nm == paste(standard.element,sep="")) {next}				

				if(ne == "LR") {Rh.corr.a <- rh[grepl("LR", names(rh))]; index <- rh$index; Rh.corr <- cbind(index,Rh.corr.a); name.1 <- names(Rh.corr)[2]}
				if(ne == "MR") {Rh.corr.a <- rh[grepl("MR", names(rh))]; index <- rh$index; Rh.corr <- cbind(index,Rh.corr.a); name.1 <- names(Rh.corr)[2]} 

				st <- merge(st1, Rh.corr, by = "index", all.y=FALSE)
						
		# Correct data for Rh standard
		if(ne == "LR") {st$Me.rh.corrected <- (st[,c(vec.ck[i])]*st[,c(paste("mean.",standard.element,".range.LR",sep=""))])/st[,c(name.1)]}		
		if(ne == "MR") {st$Me.rh.corrected <- (st[,c(vec.ck[i])]*st[,c(paste("mean.",standard.element,".range.MR",sep=""))])/st[,c(name.1)]}
				} else {st <- st1}
		
				# The metal is skipped if no slope was determined
				if(is.na(sa.HNO3[1,c(vec.ck[i])]) == TRUE) {next}

				# To correct Cd signal for Mo interference
				if(nm == "Cd" & Mo.correction.Cd == "yes" & Rh.correction == "yes"){
						
									# Select the correct slope
									if(vec.ck[i] == "Cd110.LR.") {slopeMo <- mor[1,paste("Cd110.LR.",standard.element,".corr", sep="")]}	
									if(vec.ck[i] == "Cd111.LR.") {slopeMo <- mor[1,paste("Cd111.LR.",standard.element,".corr", sep="")]}			
									
									#if(is.null(slopeMo) == TRUE) {
									#	print("ERROR: you confused the identity of the Cd istopes as used for Mo interference and as measured by ICPMS.")
									#	print("ERROR: alternatively, you used _Standard.element.correction <- no_ with the function _fMo.interference.Cd_.")}

									# Bring Mo into st
									MoIn <- cki[,c("index","Mo95.LR.")]
									st <- merge(st, MoIn, by = "index", all.y=FALSE)
									st$Me.rh.corrected <- st$Me.rh.corrected - (st$Mo95.LR./slopeMo)
									st <- cbind(st[match("index",names(st)):match("Me.rh.corrected",names(st))])
									print(paste(vec.ck[i]," checks were corrected for Rh and Mo interference in ICPMS checks",sep=""))
											}

				if(nm == "Cd" & Mo.correction.Cd == "yes" & Rh.correction == "no"){

									# Select the correct slope
									if(vec.ck[i] == "Cd110.LR.") {slopeMo <- mor[1,"Cd110.LR."]}	
									if(vec.ck[i] == "Cd111.LR.") {slopeMo <- mor[1,"Cd111.LR."]}			
									
									if(is.null(slopeMo) == TRUE) {
										#print("ERROR: you confused the identity of the Cd istopes as used for Mo interference and as measured by ICPMS.")
										print("ERROR: you have not run the function _fMo.interference.Cd_.")}

									# Bring Mo into st
									MoIn <- cki[,c("index","Mo95.LR.")]
									st <- merge(st, MoIn, by = "index", all.y=FALSE)
									st[,vec.ck[i]] <- st[,vec.ck[i]] - (st$Mo95.LR./slopeMo)
									st <- cbind(st[match("index",names(st)):match(vec.ck[i],names(st))])
									print(paste(vec.ck[i]," checks were not corrected for Rh but was for Mo interference in ICPMS checks",sep=""))
											}


				# Calculate concentrations
				
				# no blank subtracted for the ICPMS checks 
				seafast.blank <- 0
																	
				if(Rh.correction == "yes"){
					st$slope.rh.corr <- sa.HNO3[6,c(vec.ck[i])]
					st$slope.rh.uncorr <- sa.HNO3[1,c(vec.ck[i])]

					st$seafast.blank <- seafast.blank
								
					st$Me.conc.rh.corr <- (st$Me.rh.corrected/st$slope.rh.corr)-seafast.blank
					st$Me.conc.rh.uncorr <- (st[,c(vec.ck[i])]/st$slope.rh.uncorr)-seafast.blank

				} else {
					
					st$slope.rh.uncorr <- sa.HNO3[1,c(vec.ck[i])]
					st$seafast.blank <- seafast.blank
					st$Me.conc.rh.uncorr <- (st[,c(vec.ck[i])]/st$slope.rh.uncorr)-seafast.blank
				}
						
if(Rh.correction == "yes"){

					temp <- st[,c("index", "Me.conc.rh.corr")]
					names(temp)[names(temp) == "Me.conc.rh.corr"] <- paste(vec.ck[i],"conc",sep="")
					cki <- merge(cki, temp, by = "index", all.y=FALSE)
				
				
				} else {

					temp <- st[,c("index","Me.conc.rh.uncorr")]
					names(temp)[names(temp) == "Me.conc.rh.uncorr"] <- paste(vec.ck[i],"conc",sep="")
					cki <- merge(cki, temp, by = "index", all.y=FALSE)
				
				}
} # end i loop


if(no.icpms.data == 1) {print("ICPMS checks were not calculated: no data were given")} else {

ckkk.icpms <- cbind(cki[1:5],cki[(length(vec.ck)+6):dim(cki)[2]])

write.table(ckkk.icpms,file=paste(folder.check,"Checks concentrations_icpms_",name.run,".csv", sep=""),sep=",",row.names=F)

print("Calculation of concentrations ICPMS checks: completed")
	}

#### creating figures for the different checks of the Seafast together in a boxplot

vec.ch <- names(ckkk.seafast)[6:dim(ckkk.seafast)[2]]

# Create figures for the different checks & balances
for(k in 1:length(vec.ch)){
# k=6
				# select the data
				ch <- as.data.frame(ckkk.seafast[,c("index", "sample", "check","where", vec.ch[k])])
				nm <- str_sub(vec.ch[k],1,2)	#library(stringr)
				if(nm == "Yb"){nm <- "Yb"} else {				
				nm <- str_sub(vec.ch[k],1,1)
				if(nm == "V" | nm == "Y" | nm == "U"){} else { 
				nm <- str_sub(vec.ch[k],1,2)	#library(stringr)
				}
				}

				# for boxplot drop levels
				ch$check <- droplevels(ch$check)
				
				# remove rows with NA's
				z <- is.na(ch[,c(vec.ch[k])])
				ch <- ch[z==FALSE,]
				
				# The metal is skipped if no concentrations at all were determined
				nod <- mean(ch[,c(vec.ch[k])],na.rm=TRUE)
				if(is.na(nod) == TRUE) {next}

				# units for the different metals
				if(nm %in% units.nmolKg) {units <- "nmol/kg"}
				if(nm %in% units.pmolKg) {units <- "pmol/kg"}

				# Make figure for all checks as 1 boxplot
				name.fig.check <- paste(vec.ch[k],"boxplot_Checks_Seafast_",name.run,sep="_")
				tiff(file = paste(folder.check, name.fig.check, ".tiff", sep=""), width = 8, height = 5, units = "in", res=resolution.figures)
	
				#windows(8,5)
				par(mfrow=c(1,1))
				par(mar=c(8, 5, 3, 2) + 0.1)		# c(bottom, left, top, right)

				if(outliers.boxplot == "yes") {boxplot(ch[,5] ~ ch[,3],las = 2, ylab="", xlab="", 
						cex=cex, cex.lab=cex.lab,cex.axis = 0.8, main=name.run, outline=TRUE)} else {	
									 boxplot(ch[,5] ~ ch[,3],las = 2, ylab="", xlab="", 
						cex=cex, cex.lab=cex.lab,cex.axis = 0.8, main=name.run, outline=FALSE)}

						title(ylab=paste(nm," (",units,")",sep=""), cex.lab=cex.lab, line = 4)
						means <- tapply(ch[,5],ch[,3],function(x) mean(x,na.rm=TRUE))
						points(means,col="red",pch=18,cex=cex)
						legend("topright", legend = c("median","mean"), 
								col=c("black","red"),pch=c(NA,18), lty = c(1,NA), lwd=c(2,NA), cex=1, bty="n")

				# To include information in the boxplot only
				# determine the blank
				name <- paste(vec.ch[k], sep="")
				if(use.Blank == "mean") {seafast.blank <- blank[4,c(name)]}
				if(use.Blank == "median") {seafast.blank <- blank[2,c(name)]}

				mtext(paste("Seafast checks: i) ",use.Blank," blank correction of: ",seafast.blank,
					" ii) use of seafast standard addition",sep=""), side=3, lin=0.0, adj=0.0, cex=0.8, col="blue")

				if(Rh.correction == "no") {
						mtext(paste("not ",standard.element," corrected", sep=""), side=3, lin=0.1, adj=0.98, cex=0.8, col="black")
							} else {
						mtext(paste(standard.element," corrected", sep=""), side=3, lin=0.1, adj=0.98, cex=0.8, col="black")}
						
						if(Mo.correction.Cd == "no" & nm == "Cd") {
						mtext("not Mo corrected", side=3, lin=1, adj=0.98, cex=0.8, col="black")}
											 
						if(Mo.correction.Cd == "yes" & nm == "Cd") {
						mtext("Mo corrected", side=3, lin=1, adj=0.98, cex=0.8, col="black")}
			

	
				dev.off()

				graphics.off()

				
} # end loop k

print("Figures of Seafast checks: completed")

#### creating figures for the different checks of the ICPMS together in a boxplot

if(no.icpms.data == 1) {print("No figures of ICPMS checks have been created: no data were given")} else {

vec.ch <- names(ckkk.icpms)[6:dim(ckkk.icpms)[2]]

# Create figures for the different checks & balances
for(k in 1:length(vec.ch)){
# k=6
				# select the data
				ch <- as.data.frame(ckkk.icpms[,c("index", "sample", "check","where", vec.ch[k])])
				nm <- str_sub(vec.ch[k],1,2)	#library(stringr)
				if(nm == "Yb"){nm <- "Yb"} else {				
				nm <- str_sub(vec.ch[k],1,1)
				if(nm == "V" | nm == "Y" | nm == "U"){} else { 
				nm <- str_sub(vec.ch[k],1,2)	#library(stringr)
				}
				}

				# for boxplot drop levels
				ch$check <- droplevels(ch$check)
				
				# remove rows with NA's
				z <- is.na(ch[,c(vec.ch[k])])
				ch <- ch[z==FALSE,]
				
				# The metal is skipped if no concentrations at all were determined
				nod <- mean(ch[,c(vec.ch[k])],na.rm=TRUE)
				if(is.na(nod) == TRUE) {next}

				# units for the different metals
				if(nm %in% units.nmolKg) {units <- "nmol/kg"}
				if(nm %in% units.pmolKg) {units <- "pmol/kg"}

				# Make figure for all checks as 1 boxplot
				name.fig.check <- paste(vec.ch[k],"boxplot_Checks_ICPMS_",name.run,sep="_")
				tiff(file = paste(folder.check, name.fig.check, ".tiff", sep=""), width = 8, height = 5, units = "in", res=resolution.figures)
	
				#windows(8,5)
				par(mfrow=c(1,1))
				par(mar=c(7, 5, 3, 2) + 0.1)		# c(bottom, left, top, right)

				if(outliers.boxplot == "yes") {boxplot(ch[,5] ~ ch[,3],las = 2, ylab="", xlab="", 
						cex=cex, cex.lab=cex.lab,cex.axis = 0.8, main=name.run, outline = TRUE)} else {
									 boxplot(ch[,5] ~ ch[,3],las = 2, ylab="", xlab="", 
						cex=cex, cex.lab=cex.lab,cex.axis = 0.8, main=name.run, outline = FALSE)}

						title(ylab=paste(nm," (",units,")",sep=""), cex.lab=cex.lab, line = 4)
						means <- tapply(ch[,5],ch[,3],function(x) mean(x,na.rm=TRUE))
						points(means,col="red",pch=18,cex=cex)
						legend("topright", legend = c("median","mean"), 
								col=c("black","red"),pch=c(NA,18), lty = c(1,NA), lwd=c(2,NA), cex=1, bty="n")
						
				mtext("ICP-MS checks: i) uncorrected for blank, ii) use of ICPMS standard addition", side=3, lin=0.0, adj=0.0, cex=0.8, col="blue")
				
				if(Rh.correction == "no") {
						mtext(paste("not ",standard.element," corrected", sep=""), side=3, lin=0.1, adj=0.98, cex=0.8, col="black")
							} else {
						mtext(paste(standard.element," corrected", sep=""), side=3, lin=0.1, adj=0.98, cex=0.8, col="black")}
						
						if(Mo.correction.Cd == "no" & nm == "Cd") {
						mtext("not Mo corrected", side=3, lin=1, adj=0.98, cex=0.8, col="black")}
											 
						if(Mo.correction.Cd == "yes" & nm == "Cd") {
						mtext("Mo corrected", side=3, lin=1, adj=0.98, cex=0.8, col="black")}



				dev.off()

				graphics.off()

				
} # end loop k

print("Figures of ICPMS checks: completed")
		
		}  #if(no.icpms.data == 1)

# Create individual dotplots for the different checks & balances

for(k in 1:length(vec.check.seafast.dotplot)){
# k =1

				if(length(vec.check.seafast.dotplot) == 0) {break}				
				if(k > length(vec.check.seafast.dotplot)) {break}			


								# select the data
				ch2 <- subset(ckkk.seafast, check == vec.check.seafast.dotplot[k])

				# Create a separate folder for a certain check
				subDir.fig.indv.check <- paste(vec.check.seafast.dotplot[k],"_",name.run,sep="")
				dir.create(file.path(mainDir, subDir.fig.check, subDir.fig.indv.check))			
				folder.individual.check <- paste(mainDir,subDir.fig.check,"/",subDir.fig.indv.check,"/",sep="")
		

		for(j in 1:length(vec.ch)){
#j=2
				ch <- as.data.frame(ch2[,c("index", "sample", "check", "where",vec.ch[j])])
				nm <- str_sub(vec.ch[j],1,2)	#library(stringr)
				if(nm == "Yb"){nm <- "Yb"} else {				
				nm <- str_sub(vec.ch[j],1,1)
				if(nm == "V" | nm == "Y" | nm == "U"){} else { 
				nm <- str_sub(vec.ch[j],1,2)	#library(stringr)
				}
				}

				# remove rows with NA's
				z <- is.na(ch[,c(vec.ch[j])])
				ch <- ch[z==FALSE,]
				
				# The metal is skipped if no concentrations at all were determined
				nod <- mean(ch[,c(vec.ch[j])],na.rm=TRUE)
				if(is.na(nod) == TRUE) {next}

				# units for the different metals
				if(nm %in% units.nmolKg) {units <- "nmol/kg"}
				if(nm %in% units.pmolKg) {units <- "pmol/kg"}

				# The used blank for in the figures
				name <- paste(vec.ch[j], sep="")
				if(use.Blank == "mean") {seafast.blank <- blank[4,c(name)]}
				if(use.Blank == "median") {seafast.blank <- blank[2,c(name)]}


						# Make figure for all checks separately
						name.fig.check2 <- paste(vec.ch[j],vec.check.seafast.dotplot[k],name.run,sep="_")
						tiff(file = paste(folder.individual.check, name.fig.check2,".tiff", sep=""), width = 6, height = 4, units = "in", res=resolution.figures)
	
						# ylim values
						y.a <- min(c(ch[,5]), na.rm=TRUE) - (0.3 * mean(ch[,5],na.rm=TRUE))
						y.b <- max(c(ch[,5]), na.rm=TRUE) + (0.3 * mean(ch[,5],na.rm=TRUE))


						#windows(6,4)
						par(mfrow=c(1,1))
						par(mar=c(5, 5, 3, 2) + 0.1)		# c(bottom, left, top, right)

						plot(ch[,5] ~ ch[,1],las = 1, ylab="", xlab="index", pch=16,
						cex=cex, cex.lab=cex.lab,cex.axis = cex.axis, main=name.run,
						ylim=c(y.a,y.b))
						title(ylab=paste(nm," (",units,")",sep=""), cex.lab=cex.lab, line = 4)
						mtext(paste("Seafast checks & balances: ",vec.check.seafast.dotplot[k],", blank corrected",sep=""), side=3, lin=0.1, adj=0.0, cex=0.8, col="blue")
						mtext(paste(use.Blank," blank: ",round(seafast.blank,3),sep=""), 
											side=1, lin=-1, adj=0.01, cex=0.8, col="blue")
		
						if(Rh.correction == "no") {
						mtext(paste("not ",standard.element," corrected", sep=""), side=3, lin=0.1, adj=0.98, cex=0.8, col="black")
							} else {
						mtext(paste(standard.element," corrected", sep=""), side=3, lin=0.1, adj=0.98, cex=0.8, col="black")}
						
						if(Mo.correction.Cd == "no" & nm == "Cd") {
						mtext("not Mo corrected", side=3, lin=1, adj=0.98, cex=0.8, col="black")}
											 
						if(Mo.correction.Cd == "yes" & nm == "Cd") {
						mtext("Mo corrected", side=3, lin=1, adj=0.98, cex=0.8, col="black")}



				dev.off()

				graphics.off()

#print(paste("k = ",k," : ",vec.check.seafast.dotplot[k]," and j = ",j," : ",vec.ch[j]), sep="")

						} # end loop j
} # end loop k

print("Figures of the individual Seafast checks: completed")


# Create individual dotplots for the different checks & balances of the ICPMS

if(no.icpms.data == 1) {print("No figures of the individual ICPMS checks have been created: no data were given")} else {

for(k in 1:length(vec.check.icpms.dotplot)){
# k =1
				if(length(vec.check.icpms.dotplot) == 0) {break}
				if(k > length(vec.check.icpms.dotplot)) {break}			

				# select the data
				ch2 <- subset(ckkk.icpms, check == vec.check.icpms.dotplot[k])

				# Create a separate folder for a certain check
				subDir.fig.indv.check <- paste(vec.check.icpms.dotplot[k],"_",name.run,sep="")
				dir.create(file.path(mainDir, subDir.fig.check, subDir.fig.indv.check))			
				folder.individual.check <- paste(mainDir,subDir.fig.check,"/",subDir.fig.indv.check,"/",sep="")
		

		for(j in 1:length(vec.ch)){
#j=1
				ch <- as.data.frame(ch2[,c("index", "sample", "check", "where",vec.ch[j])])
				nm <- str_sub(vec.ch[j],1,2)	#library(stringr)
				if(nm == "Yb"){nm <- "Yb"} else {				
				nm <- str_sub(vec.ch[j],1,1)
				if(nm == "V" | nm == "Y" | nm == "U"){} else { 
				nm <- str_sub(vec.ch[j],1,2)	#library(stringr)
				}
				}

				# remove rows with NA's
				z <- is.na(ch[,c(vec.ch[j])])
				ch <- ch[z==FALSE,]
				
				# The metal is skipped if no concentrations at all were determined
				nod <- mean(ch[,c(vec.ch[j])],na.rm=TRUE)
				if(is.na(nod) == TRUE) {next}

				# units for the different metals
				if(nm %in% units.nmolKg) {units <- "nmol/kg"}
				if(nm %in% units.pmolKg) {units <- "pmol/kg"}

						# Make figure for all checks separately
						name.fig.check2 <- paste(vec.ch[j],vec.check.icpms.dotplot[k],name.run,sep="_")
						tiff(file = paste(folder.individual.check, name.fig.check2,".tiff", sep=""), width = 6, height = 4, units = "in", res=resolution.figures)
	
						# ylim values
						y.a <- min(c(ch[,5]), na.rm=TRUE) - (0.3 * mean(ch[,5],na.rm=TRUE))
						y.b <- max(c(ch[,5]), na.rm=TRUE) + (0.3 * mean(ch[,5],na.rm=TRUE))


						#windows(6,4)
						par(mfrow=c(1,1))
						par(mar=c(5, 5, 3, 2) + 0.1)		# c(bottom, left, top, right)

						plot(ch[,5] ~ ch[,1],las = 1, ylab="", xlab="index", pch=16,
						cex=cex, cex.lab=cex.lab,cex.axis = cex.axis, main=name.run,
						ylim=c(y.a,y.b))
						title(ylab=paste(nm," (",units,")",sep=""), cex.lab=cex.lab, line = 4)
						mtext(paste("ICPMS checks & balances: ",vec.check.icpms.dotplot[k],", no blank correction",sep=""), side=3, lin=0.1, adj=0.0, cex=0.8, col="blue")
						mtext("Warning: no preconc. factor used, don't compare with seafast data", 
											side=1, lin=-1, adj=0.01, cex=0.8, col="red")
		
						if(Rh.correction == "no") {
						mtext(paste("not ",standard.element," corrected", sep=""), side=3, lin=0.1, adj=0.98, cex=0.8, col="black")
							} else {
						mtext(paste(standard.element," corrected", sep=""), side=3, lin=0.1, adj=0.98, cex=0.8, col="black")}
						
						if(Mo.correction.Cd == "no" & nm == "Cd") {
						mtext("not Mo corrected", side=3, lin=1, adj=0.98, cex=0.8, col="black")}
											 
						if(Mo.correction.Cd == "yes" & nm == "Cd") {
						mtext("Mo corrected", side=3, lin=1, adj=0.98, cex=0.8, col="black")}


				dev.off()

				graphics.off()

#print(paste("k = ",k," : ",vec.check.icpms.dotplot[k]," and j = ",j," : ",vec.ch[j]), sep="")

						} # end loop j
} # end loop k

print("Figures for the individual ICPMS checks: completed")

} #if(no.icpms.data == 1)

## Calculate the recovery for the different metals

# Function to properly designate the header in the output file.
# From internet: http://stackoverflow.com/questions/2478352/write-table-in-r-screws-up-header-when-has-rownames
		my.write <- function(x, file, header, f = write.csv, ...){
						# create and open the file connection
		datafile <- file(file, open = 'wt')
						# close on exit 
		on.exit(close(datafile))
						# if a header is defined, write it to the file (@CarlWitthoft's suggestion)
		if(!missing(header)) {
		writeLines(header,con=datafile, sep='\t')
		writeLines('', con=datafile, sep='\n')
		}
						# write the file using the defined function and required addition arguments  
		f(x, datafile,...)
		}


if(no.icpms.data == 1) {print("Recoveries have not been calculated: no standard addition in eluent directly measured by ICPMS (no seafast) has been performed")} else {


# code 070616_1
if(Rh.correction == "yes"){
				rec.HNO3 <- sa.HNO3[6,]
				rec.sf <- sa.sf[6,]/Seafast.conc.factor
				rec <- rbind(rec.HNO3,rec.sf)

				rec <- setNames(data.frame(t(rec[,-1])), rec[,1])
				colnames(rec) <- c("slope.Rh.corr.HNO3","slope.Rh.corr.sf")
				rec$perc.recovery <- round((rec$slope.Rh.corr.sf*100)/rec$slope.Rh.corr.HNO3,2)	
				recfinal <- cbind(rec[3])

				print(paste("Calculation of recoveries based on ",standard.element," corrected standard addition slopes : completed", sep=""))

				}

if(Rh.correction == "no"){
				rec.HNO3 <- sa.HNO3[1,]
				rec.sf <- sa.sf[1,]/Seafast.conc.factor
				rec <- rbind(rec.HNO3,rec.sf)

				rec <- setNames(data.frame(t(rec[,-1])), rec[,1])
				colnames(rec) <- c("slope.HNO3","slope.sf")
				rec$perc.recovery <- round((rec$slope.sf*100)/rec$slope.HNO3,2)	
				recfinal <- cbind(rec[3])

				print("Calculation of recoveries based on uncorrected standard addition slopes : completed")

				}

my.write(recfinal,file=paste(folder.check,"Percentage_recovery_",name.run,".csv", sep=""),sep=",",row.names=TRUE)


} #if(no.icpms.data == 1)

# Create a csv file with only the results of the reference materials

# check if reference values have been used for the SEAFAST

x <- subset(ckkk.seafast, !is.na(reference))

if(dim(x)[1] > 0) {

			ref.res <- ckkk.seafast 

			# remove rows with NA's reference
			z <- is.na(ref.res[,"reference"])
			ref.res <- ref.res[z==FALSE,]

			# order the data frame names(ref.res)
			ref.res <- ref.res[order(ref.res$reference, ref.res$check),]
			vec.metalnames <- names(ref.res[6:dim(ref.res)[2]])

			# Add run identifier useful for later collation of all reference material results
			ref.res$name.run <- name.run

			# reorder the columns
			ref.res <- ref.res[c("index","sample","check","reference","where","name.run",vec.metalnames)]

			# create the csv file
			my.write(ref.res,file=paste(folder.check,"Results_ReferenceMaterials_",name.run,".csv", sep=""),sep=",",row.names=FALSE)

			# Report
			print("A csv file reporting the results of the Reference materials has been created.")

	# Make figures with the reference consensus values

if(exists("ref.mat") == FALSE) {
	stop("your file TraceMetalReferenceConsensusValues.csv has not been placed in your working directory or the name is incorrectly spelled")}				

	for(k in 1:length(vec.ch)){
	# k=1
				# select the data
				cr <- as.data.frame(ckkk.seafast[,c("index", "sample", "check","reference", vec.ch[k])])
				nm <- str_sub(vec.ch[k],1,2)	#library(stringr)
				if(nm == "Yb"){nm <- "Yb"} else {				
				nm <- str_sub(vec.ch[k],1,1)
				if(nm == "V" | nm == "Y" | nm == "U"){} else { 
				nm <- str_sub(vec.ch[k],1,2)	#library(stringr)
				}
				}

				# do next if the metal has no reference value
				if(nm == "Cd" | nm == "Co" | nm == "Cu" | nm == "Fe" | nm == "Mn" | nm == "Ni" | nm == "Pb" | nm == "Zn") {} else {next}

				# remove rows with NA's in metals
				z <- is.na(cr[,c(vec.ch[k])])
				cr <- cr[z==FALSE,]

				# remove rows with NA's in reference
				z <- is.na(cr[,c("reference")])
				cr <- cr[z==FALSE,]

				# The metal is skipped if no concentrations at all were determined
				nod <- mean(cr[,c(vec.ch[k])],na.rm=TRUE)
				if(is.na(nod) == TRUE) {next}

				vec.rf <- unique(cr$reference)

				# remove NA's in vec.rf
				z <- is.na(vec.rf)
				vec.rf <- vec.rf[z==FALSE]

				# The used blank for in the figures
				name <- paste(vec.ch[k], sep="")
				if(use.Blank == "mean") {seafast.blank <- blank[4,c(name)]}
				if(use.Blank == "median") {seafast.blank <- blank[2,c(name)]}

				# units for the different metals
				if(nm %in% units.nmolKg) {units <- "nmol/kg"}
				if(nm %in% units.pmolKg) {units <- "pmol/kg"}

				
		for(j in 1:length(vec.rf)){
		# j=1; ref.mat
						ctr <- subset(cr, reference == vec.rf[j])
						ctr$newIndex <- seq(1,dim(ctr)[1],1)						

						# include sample code in ctr data.frame
						if(sample.naming == "yes")	{
								ctr <- merge(ctr, coding, by.x = "sample", by.y = "seafast.name", all = FALSE) 
											}
						if(dim(ctr)[1] == 0) {
		stop("Did you fill in the file sample_list.csv in a correct way? If you don't want to include sample names define sample.names <- no")}				

						# to select the reference value
						vec.ref.code <- ref.mat$code
						if(vec.rf[j] %in% vec.ref.code) 	{
						   		sel.ref <- ref.mat[ref.mat$code == as.character(vec.rf[j]),]
								Me.mean.ref <- sel.ref[,c(paste(nm,"mean",sep="_"))] 	
								Me.stdev.ref <- sel.ref[,c(paste(nm,"stdev",sep="_"))]			
								Ref.mat <- as.character(sel.ref[,c("Ref")])
																	

						# to calculate the % of between our reference value and the consensus value
						ctr$perc.difference <- (((ctr[,5] - Me.mean.ref)*100)/Me.mean.ref)
						mean.perc.difference <- mean(ctr[,"perc.difference"], na.rm=TRUE)
						stdev.perc.difference <- sd(ctr[,"perc.difference"], na.rm=TRUE)
												}
						
						# Make figure 
						name.fig.ref <- paste(vec.ch[k],vec.rf[j],"RefMat",name.run,sep="_")
						tiff(file = paste(folder.ref.met, name.fig.ref, ".tiff", sep=""), width = 6, height = 4, units = "in", res=resolution.figures)

						# xlim values
						x.a <- 0
						x.b <- dim(ctr)[1]+1

						# ylim values
						if(vec.rf[j] %in% vec.ref.code) 	{
							y.a <- min(c(ctr[,5],Me.mean.ref), na.rm=TRUE) - (0.3 * mean(ctr[,5],na.rm=TRUE))
							y.b <- max(c(ctr[,5],Me.mean.ref), na.rm=TRUE) + (0.3 * mean(ctr[,5],na.rm=TRUE))
												} else {
							y.a <- min(c(ctr[,5]), na.rm=TRUE) - (0.3 * mean(ctr[,5],na.rm=TRUE))
							y.b <- max(c(ctr[,5]), na.rm=TRUE) + (0.3 * mean(ctr[,5],na.rm=TRUE))
												}
			
						#windows(6,4)
						par(mfrow=c(1,1))
						par(mar=c(5, 5, 3, 2) + 0.1)		# c(bottom, left, top, right)

						plot(ctr[,5] ~ ctr$newIndex,las = 1, ylab="", xlab="arbitrary index", pch=16,
						cex=cex, cex.lab=cex.lab,cex.axis = cex.axis, main=name.run, xlim=c(x.a,x.b), ylim=c(y.a,y.b))
						title(ylab=paste(nm," (",units,")",sep=""), cex.lab=cex.lab, line = 4)
						text(ctr$newIndex,ctr[,5], label=round(ctr[,5],2), pos=2) 
						mtext(paste(use.Blank," blank: ",round(seafast.blank,3),sep=""), 
											side=1, lin=-1, adj=0.01, cex=0.8, col="blue")
						if(sample.naming == "yes")	{
								text(ctr$newIndex,ctr[,5], label=ctr$sample.name, pos=3, cex=0.5, offset=1.5) 				
						}
						
						if(vec.rf[j] %in% vec.ref.code) 	{
						
						mtext(paste("Reference material: ",Ref.mat,sep=""), side=3, lin=0.1, adj=0.0, cex=0.8, col="red")
						points(Me.mean.ref ~ x.b, col="red", pch=16, cex=cex)
						plotCI(x.b, Me.mean.ref, err="y", uiw=Me.stdev.ref, liw=Me.stdev.ref, add=TRUE,col="red", gap=0, sfrac=0.005)		
						text(x.b,Me.mean.ref, label=round(Me.mean.ref,2), pos=2) 				
						mtext(paste("% off from consensus, mean: ",round(mean.perc.difference,2)," +/- ", round(stdev.perc.difference,2),
								sep=""), side=1, lin=-1, adj=0.99, cex=0.8, col="blue")
						}

						if(Rh.correction == "no") {
						mtext(paste("not ",standard.element," corrected", sep=""), side=3, lin=0.1, adj=0.98, cex=0.8, col="black")
							} else {
						mtext(paste(standard.element," corrected", sep=""), side=3, lin=0.1, adj=0.98, cex=0.8, col="black")}
						
						if(Mo.correction.Cd == "no" & nm == "Cd") {
						mtext("not Mo corrected", side=3, lin=1, adj=0.98, cex=0.8, col="black")}
											 
						if(Mo.correction.Cd == "yes" & nm == "Cd") {
						mtext("Mo corrected", side=3, lin=1, adj=0.98, cex=0.8, col="black")}


	
				dev.off()

				graphics.off()

				} # end loop j
		
} # end loop k

print("Figures of the reference samples with consencus values: completed")

} else {print("No reference materials have been used with the Seafast.")} 

} # end function


#################################################################################################################################################

fsample.concentrations <- function(fr, sa.sf, st.el, blank, mor, name.run, mainDir, use.Blank, resolution.figures, cex, cex.lab, cex.axis, Standard.element.correction,Mo.correction.Cd,...){

options(warn=-1) # to terminate all warning messages

Rh.correction <- Standard.element.correction
rh <- st.el

if(Rh.correction == "yes"){
				tst2 <- sa.sf[6,]
				tst2a <- cbind(tst2[2:dim(tst2)[2]])
				if(is.na(rowMeans(tst2a, na.rm=TRUE)[1][[1]]) == TRUE) {stop(paste("You need to set Standard.element.correction at no as you did not ",standard.element," correct your SEAFAST standard additions.", sep=""))}
				}

if(!("index" %in% colnames(fr))) {stop("ERROR: The column _index_ is missing in the file samples.csv or you wrote index with a capital I instead of i.")}
if(!("sample" %in% colnames(fr))) {stop("ERROR: The column _sample_ is missing in the file samples.csv or was written incorrectly / don't use capital letters")}
if(!("index" %in% colnames(rh))) {stop("ERROR: The column _index_ is missing or you wrote index with a capital I instead of i in file Rh_standard.csv.")}

if(Rh.correction == "yes" & Mo.correction.Cd == "no") {print(paste("All data are ",standard.element," corrected. Cd has not been corrected for Mo interference.", sep=""))} 
if(Rh.correction == "yes" & Mo.correction.Cd == "yes") {print(paste("All data are ",standard.element," corrected. Cd has been corrected for Mo interference. Cd and Mo were ",standard.element," corrected.", sep=""))} 
if(Rh.correction == "no" & Mo.correction.Cd == "yes") {print(paste("All data are not ",standard.element," corrected. Cd has been corrected for Mo interference. Cd and Mo were not ",standard.element," corrected.", sep=""))} 
if(Rh.correction == "no" & Mo.correction.Cd == "no") {print(paste("All data are not ",standard.element," corrected. Cd has not been corrected for Mo interference.", sep=""))} 


vec.fr <- names(fr)[3:dim(fr)[2]]
frr <- fr

# Create new file folder for figures and data of the blanks
	subDir.fig.samples <- paste("Figures & data_samples",name.run,sep="_")
	dir.create(file.path(mainDir, subDir.fig.samples))			
	folder.samples <<- paste(mainDir,subDir.fig.samples,"/",sep="")
	
	
	for(i in 1:length(vec.fr)){
#i=3

				# select data & safe original data frame
				st1 <- frr[,c("index","sample",vec.fr[i])]
				nm <- str_sub(vec.fr[i],1,2)	#library(stringr)
				if(nm == "Yb"){nm <- "Yb"} else {				
				nm <- str_sub(vec.fr[i],1,1)
				if(nm == "V" | nm == "Y" | nm == "U"){} else { 
				nm <- str_sub(vec.fr[i],1,2)	#library(stringr)
				}
				}

				if(Rh.correction == "yes"){
				# selection of LR or MR Rh standard data
				ne <- str_sub(vec.fr[i],-3,-2)

				if(exists("ne") == FALSE) {print("ERROR: you are working in an older version of R that does not support the stringr library!")}				
				if(ne == "") {print("ERROR: you are working in an older version of R that does not support the stringr library!")}				

				if(nm == paste(standard.element,sep="")) {next}				

				if(ne == "LR") {Rh.corr.a <- rh[grepl("LR", names(rh))]; index <- rh$index; Rh.corr <- cbind(index,Rh.corr.a); name.1 <- names(Rh.corr)[2]}
				if(ne == "MR") {Rh.corr.a <- rh[grepl("MR", names(rh))]; index <- rh$index; Rh.corr <- cbind(index,Rh.corr.a); name.1 <- names(Rh.corr)[2]} 

				st <- merge(st1, Rh.corr, by = "index", all.y=FALSE)
						
		# Correct data for Rh standard
		if(ne == "LR") {st$Me.rh.corrected <- (st[,c(vec.fr[i])]*st[,c(paste("mean.",standard.element,".range.LR",sep=""))])/st[,c(name.1)]}		
		if(ne == "MR") {st$Me.rh.corrected <- (st[,c(vec.fr[i])]*st[,c(paste("mean.",standard.element,".range.MR",sep=""))])/st[,c(name.1)]}
				} else {st <- st1}
		
	
				# The metal is skipped if no slope was determined
				if(is.na(sa.sf[1,c(vec.fr[i])]) == TRUE) {next}

				# To correct Cd signal for Mo interference
				if(nm == "Cd" & Mo.correction.Cd == "yes" & Rh.correction == "yes"){
					if(Rh.correction == "yes" & is.na(mor[1,3]) == TRUE){stop("ERROR: you have run fMo.interference.Cd without standard.element correction. 
You have 3 options: 	
i) run fMo.interference.Cd with correction,
ii) run fsample.concentrations with Standard.element.correction <- no
iii) run fsample.concentrations with Mo.correction.Cd <- no")}
	
									# Select the correct slope
									if(vec.fr[i] == "Cd110.LR.") {slopeMo <- mor[1,paste("Cd110.LR.",standard.element,".corr", sep="")]}	
									if(vec.fr[i] == "Cd111.LR.") {slopeMo <- mor[1,paste("Cd111.LR.",standard.element,".corr", sep="")]}			
																
									#if(is.null(slopeMo) == TRUE) {
									#	print("ERROR: you confused the identity of the Cd istopes as used for Mo interference and as measured by ICPMS.")
									#	print("ERROR: alternatively, you used _Standard.element.correction <- no_ with the function _fMo.interference.Cd_.")}


									# Bring Mo into st
									MoIn <- frr[,c("index","Mo95.LR.")]
									st <- merge(st, MoIn, by = "index", all.y=FALSE)
									st$Me.rh.corrected <- st$Me.rh.corrected - (st$Mo95.LR./slopeMo)
									st <- cbind(st[match("index",names(st)):match("Me.rh.corrected",names(st))])
									print(paste(vec.fr[i]," samples were corrected for Rh and Mo interference in ICPMS checks",sep=""))
											}

				if(nm == "Cd" & Mo.correction.Cd == "yes" & Rh.correction == "no"){

									# Select the correct slope
									if(vec.fr[i] == "Cd110.LR.") {slopeMo <- mor[1,"Cd110.LR."]}	
									if(vec.fr[i] == "Cd111.LR.") {slopeMo <- mor[1,"Cd111.LR."]}			
									
									if(is.null(slopeMo) == TRUE) {
										#print("ERROR: you confused the identity of the Cd istopes as used for Mo interference and as measured by ICPMS.")
										print("ERROR: You have not run the function _fMo.interference.Cd_.")}

									# Bring Mo into st
									MoIn <- frr[,c("index","Mo95.LR.")]
									st <- merge(st, MoIn, by = "index", all.y=FALSE)
									st[,vec.fr[i]] <- st[,vec.fr[i]] - (st$Mo95.LR./slopeMo)
									st <- cbind(st[match("index",names(st)):match(vec.fr[i],names(st))])
									#print(paste(vec.fr[i]," samples were not corrected for ",standard.element," but was for Mo interference in ICPMS checks",sep=""))
											}


				# Calculate concentrations

				name <- paste(vec.fr[i],"conc", sep="")
				if(use.Blank == "mean") {seafast.blank <- blank[4,c(name)]}
				if(use.Blank == "median") {seafast.blank <- blank[2,c(name)]}

				if(Rh.correction == "yes"){
					st$slope.rh.corr <- sa.sf[6,c(vec.fr[i])]
					st$slope.rh.uncorr <- sa.sf[1,c(vec.fr[i])]

					st$seafast.blank <- seafast.blank
								
					st$Me.conc.rh.corr <- (st$Me.rh.corrected/st$slope.rh.corr)-seafast.blank
					st$Me.conc.rh.uncorr <- (st[,c(vec.fr[i])]/st$slope.rh.uncorr)-seafast.blank

				} else {
					
					st$slope.rh.uncorr <- sa.sf[1,c(vec.fr[i])]
					st$seafast.blank <- seafast.blank
					st$Me.conc.rh.uncorr <- (st[,c(vec.fr[i])]/st$slope.rh.uncorr)-seafast.blank
				}
						
if(Rh.correction == "yes"){

					temp <- st[,c("index", "Me.conc.rh.corr")]
					names(temp)[names(temp) == "Me.conc.rh.corr"] <- paste(vec.fr[i],"conc",sep="")
					frr <- merge(frr, temp, by = "index", all.y=FALSE)
				
				
				} else {

					temp <- st[,c("index","Me.conc.rh.uncorr")]
					names(temp)[names(temp) == "Me.conc.rh.uncorr"] <- paste(vec.fr[i],"conc",sep="")
					frr <- merge(frr, temp, by = "index", all.y=FALSE)
				
				}
} # end i loop

# head(frr); names(frr)

frrr <- cbind(frr[1:2],frr[(length(vec.fr)+3):dim(frr)[2]])

write.table(frrr,file=paste(folder.samples,"Sample concentrations_",name.run,".csv", sep=""),sep=",",row.names=F)

#### creating figures for the sample concentrations

vec.frr <- names(frrr)[3:dim(frrr)[2]]

# Create figures for the different metals of the samples

for(k in 1:length(vec.frr)){
# k=1
				# select the data
				ch <- as.data.frame(frrr[,c("index", "sample", vec.frr[k])])
				nm <- str_sub(vec.frr[k],1,2)	#library(stringr)
				if(nm == "Yb"){nm <- "Yb"} else {				
				nm <- str_sub(vec.frr[k],1,1)
				if(nm == "V" | nm == "Y" | nm == "U"){} else { 
				nm <- str_sub(vec.frr[k],1,2)	#library(stringr)
				}
				}

				# remove rows with NA's
				z <- is.na(ch[,c(vec.frr[k])])
				ch <- ch[z==FALSE,]
				
				# The metal is skipped if no concentrations at all were determined
				nod <- mean(ch[,c(vec.frr[k])],na.rm=TRUE)
				if(is.na(nod) == TRUE) {next}

				# units for the different metals
				if(nm %in% units.nmolKg) {units <- "nmol/kg"}
				if(nm %in% units.pmolKg) {units <- "pmol/kg"}

				# Make figure for all checks as 1 boxplot
				name.fig.samples <- paste(vec.frr[k],"plot_samples",name.run,sep="_")
				tiff(file = paste(folder.samples, name.fig.samples, ".tiff", sep=""), width = 12, height = 6, units = "in", res=resolution.figures)
	
						#windows(12,6)
						par(mfrow=c(1,1))
						par(mar=c(5, 5, 3, 2) + 0.1)		# c(bottom, left, top, right)

						plot(ch[,3] ~ ch[,1],las = 1, ylab="", xlab="index", pch=16,
						cex=cex, cex.lab=cex.lab,cex.axis = cex.axis, main=name.run)
						title(ylab=paste(nm," (",units,")",sep=""), cex.lab=cex.lab, line = 4)
						if(Rh.correction == "yes") {
							if(nm == "Cd" & Mo.correction.Cd == "yes"){mtext(paste("The final concentration of the samples according to their given index: ",standard.element," corrected, corrected for Mo interference",sep=""), side=3, lin=0.1, adj=0.0, cex=0.8, col="blue")} else {
							mtext(paste("The final concentration of the samples according to their given index: ",standard.element," corrected", sep=""), side=3, lin=0.1, adj=0.0, cex=0.8, col="blue")}
											}
						if(Rh.correction == "no") {
							if(nm == "Cd" & Mo.correction.Cd == "yes"){mtext(paste("The final concentration of the samples according to their given index: not ",standard.element," corrected, corrected for Mo interference", sep=""), side=3, lin=0.1, adj=0.0, cex=0.8, col="blue")} else {
							mtext(paste("The final concentration of the samples according to their given index: not ",standard.element," corrected", sep=""), side=3, lin=0.1, adj=0.0, cex=0.8, col="blue")}
											}
						
						#if(Rh.correction == "no") {
						#mtext("The final concentration of the samples according to their given index: Rh corrected", side=3, lin=0.1, adj=0.0, cex=0.8, col="blue")
						#
						#mtext("no Rh correction", side=3, lin=-1, adj=0.01, cex=0.8, col="red")
						#					}
						if(Mo.correction.Cd == "no" & nm == "Cd") {
						mtext("no Mo correction", side=3, lin=-1, adj=0.01, cex=0.8, col="red")
											}


				dev.off()

				graphics.off()

} # end loop k

print("The metal concentrations for the samples have been calculated")

} # end function


##################################################################################################################################
######## Function to include original sample names to the final concentrations results ###########################################

fsample.codes <- function(samples, coding,...){

options(warn=-1) # to terminate all warning messages

if(!("index" %in% colnames(samples))) {stop("ERROR: The column _index_ is missing in the file samples.csv or you wrote index with a capital I instead of i.")}
if(!("sample" %in% colnames(samples))) {stop("ERROR: The column _sample_ is missing in the file samples.csv or was written incorrectly / don't use capital letters")}
if(!("seafast.name" %in% colnames(coding))) {stop("ERROR: The column _seafast.name_ is missing in the file sample_list.csv or was written incorrectly / don't use capital letters")}
if(!("sample.name" %in% colnames(coding))) {stop("ERROR: The column _sample.name_ is missing in the file sample_list.csv or was written incorrectly / don't use capital letters")}

# Work within file folder for figures and data of the blanks
	subDir.fig.samples <- paste("Figures & data_samples",name.run,sep="_")
	dir.create(file.path(mainDir, subDir.fig.samples))			
	folder.samples <<- paste(mainDir,subDir.fig.samples,"/",sep="")

# remove index from coding
	cod <- cbind(coding[match("seafast.name",names(coding)):match("sample.name",names(coding))])			

# megre the files
	sampcod <- merge(samples, cod, by.x = "sample", by.y = "seafast.name", all.y=FALSE)
	
# replace original file with concentrations with file that includes sample names
	write.table(sampcod,file=paste(folder.samples,"Sample concentrations_",name.run,".csv", sep=""),sep=",",row.names=F)

print("Original sample coding added to the .csv file with final metal concentrations")

} # end function


################################################################################################################################






		