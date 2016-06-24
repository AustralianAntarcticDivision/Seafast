# Script for the analysis of SEAFAST ICP-MS data
# Micha Rijkenberg (micharijkenberg@yahoo.com.au)
# Version: 07 June 2016
# Version log is heading the functions script

### TO SET UP R FOR USE OF THE SCRIPT
# Activate the line below (remove "#" and run) to download the necessary packages
# you only have to do this the first time that you want to use a package
# install.packages("stringr", dependencies = TRUE)
# install.packages("TeachingDemos", dependencies = TRUE)
# install.packages("plyr", dependencies = TRUE)
# install.packages("gplots", dependencies = TRUE)

#### TO USE THIS R SCRIPT
# To start with a clean sheet (is not necessary, just neat)
rm(list = ls())

# To activate the necessary libraries
library(stringr)
library(TeachingDemos)
library(plyr)
library(gplots)

##### This script is easiest to use with the buffered output turned off, i.e.:
# 1) click with your mouse on the R Console
# 2) use one time "cntrl + W" to turn buffered output off, repeating turns it on again

##### NEEDS YOUR INPUT
# Are you using a windows computer or a mac: i.e. "windows" or "mac"
Using.computer <- "windows"

##### NEEDS YOUR INPUT, SETTING THE R WORKING DIRECTORY
# setwd("E:/NIOZ_this computer/Patrick_Seafast/public version/"): KEEP EXAMPLE
setwd("E:/NIOZ_this computer/Public R_code/Seafast ICPMS R_code/")

##### THE GENERAL PARAMETERS
##### TO ACTIVATE THE FUNCTIONS, DON'T CHANGE THE CODE HERE
source("Functions_SeafastICPMS_DataAnalysis_080616.r")

##### The units used for the different parameters
units.nmolKg <- c("Fe", "Mn", "Ni", "Zn","Cd", "Cu")
units.pmolKg <- c("Pb","Co","Ti")

##### NEEDS YOUR INPUT: The element used as standard to correct for variability in the ICPMS sensitivity, e.g. "Rh" of "Ga" or....
standard.element <- "Rh"

##### NEEDS YOUR INPUT: Defines the directory where the data can be found and the results will be placed
# mainDir: The main directory in which will be worked
# name.run: the name of the icpms run with the series of sample, it will appear in file names and figures
mainDir <- 	"E:/NIOZ_this computer/Public R_code/Seafast ICPMS R_code/test case/"	
name.run <- 20160304

# To allow to access files created in previous script use for a certain ICPM-MS run activate below function
fAcces.prev.files(mainDir, name.run)

################ Standard Element ##############################################################
# With Standard Element the elment like Rh or Ga is meant which is used 
# among others to correct for variation in the ICPMS sensitivity.

# To read in the results of the StandardElement.csv, DON'T CHANGE THE CODE HERE
st.el <- read.csv(paste(mainDir,"StandardElement.csv",sep=""),header=TRUE)
head(st.el)

# Graphing parameters, CAN BE CHANGED
resolution.figures <- 150
cex=1.5
cex.lab=1.5
cex.axis = 1.5

# NEEDS YOUR INPUT, These lines separate the different Standard Element concentrations used.
# With the line.high and line.low you can divide the Standard.Element in a maximum of three distinct groups.
# You can fill in any number for the lines as long line.high > line.low and as long there is a number, at least a 0. 
# LR
line.high.lr <- 2.78e6
line.low.lr <- 1.5e6

# MR
line.high.mr <- 4e4	
line.low.mr <- 1e4

# Function, DON'T CHANGE THE CODE HERE
fStElement.stand(st.el, name.run, mainDir, resolution.figures, cex, cex.lab, cex.axis)
		
######### Standard addition in HNO3 ##########################################################

# To read in the results of the ICPMS, DON'T CHANGE THE CODE HERE 
ic <- read.csv(paste(mainDir,"StanAdd_HNO3.csv",sep=""),header=TRUE)
head(ic)

# To read in the added metal concentrations for the standard addition, DON'T CHANGE THE CODE HERE
ad <- read.csv(paste(mainDir,"Added_HNO3.csv",sep=""),header=TRUE)
head(ad)

# To read in the results of the standard.element, DON'T CHANGE THE CODE HERE
st.el <- read.csv(paste(folder.rh,standard.element," standard data_",name.run,".csv",sep=""),header=TRUE)
head(st.el)

# Graphing parameters, CAN BE CHANGED
resolution.figures <- 150
cex=1.3
cex.lab=1.2
cex.axis = 1

# Other parameters, CAN BE CHANGED
# Standard.element correction
Standard.element.correction <- "yes"

# Function, DON'T CHANGE THE CODE HERE
fStandAdd.HNO3(ic, ad, st.el, name.run, mainDir, resolution.figures, cex, cex.lab, cex.axis, Standard.element.correction)

###########################################################################################
######### Standard addition in seawater over the seafast ##################################

# To read in the results of the ICPMS, DON'T CHANGE THE CODE HERE 
sf <- read.csv(paste(mainDir,"StanAdd_seafast.csv",sep=""),header=TRUE)
head(sf)

# To read in the added metal concentrations for the standard addition, DON'T CHANGE THE CODE HERE
ad.sf <- read.csv(paste(mainDir,"Added_seafast.csv",sep=""),header=TRUE)
head(ad.sf)

# To read in the results of the standard.element, DON'T CHANGE THE CODE HERE
st.el <- read.csv(paste(folder.rh,standard.element," standard data_",name.run,".csv",sep=""),header=TRUE)
head(st.el)

# Graphing parameters, CAN BE CHANGED
resolution.figures <- 150
cex=1.3
cex.lab=1.2
cex.axis = 1

# Other parameters, CAN BE CHANGED
# Standard element correction
Standard.element.correction <- "yes"

# Function, DON'T CHANGE THE CODE HERE
fStandAdd.seafast(sf, ad.sf, st.el, name.run, mainDir, resolution.figures, cex, cex.lab, cex.axis, Standard.element.correction)

##########################################################################################
### Correct Cd concentrations for Mo interference #######################################

# To read in the results of the ICPMS for the checks and balances, DON'T CHANGE THE CODE HERE 
mo <- read.csv(paste(mainDir,"Mo_interference.csv",sep=""),header=TRUE)
head(mo)

# To read in the results of the standard.element, DON'T CHANGE THE CODE HERE
st.el <- read.csv(paste(folder.rh,standard.element," standard data_",name.run,".csv",sep=""),header=TRUE)
head(st.el)

# Graphing parameters, CAN BE CHANGED
resolution.figures <- 150
cex=1.3
cex.lab=1.2
cex.axis = 1

# Other parameters, CAN BE CHANGED
# Standard element correction
Standard.element.correction <- "yes"

# Function, DON'T CHANGE THE CODE HERE
fMo.interference.Cd(mo, st.el, resolution.figures, cex, cex.lab, cex.axis, Standard.element.correction)

###########################################################################################
######### Blank determination #############################################################

# To read in the results of the ICPMS, DON'T CHANGE THE CODE HERE 
blnk <- read.csv(paste(mainDir,"blanks.csv",sep=""),header=TRUE)
head(blnk)
table(blnk$blank)

# To read in the results of the seafast standard additions, DON'T CHANGE THE CODE HERE
sa.sf <- read.csv(paste(folder.sf,"StandAdd_seafast_results_",name.run,".csv", sep=""),header=TRUE)
head(sa.sf)

# To read in the results of the standard.element, DON'T CHANGE THE CODE HERE
st.el <- read.csv(paste(folder.rh,standard.element," standard data_",name.run,".csv",sep=""),header=TRUE)
head(st.el)

# To read in the results of Mo interference, DON'T CHANGE THE CODE HERE
mor <- read.csv(paste(folder.Mo,"Mo_interferenceCd_results_",name.run,".csv",sep=""),header=TRUE)
head(mor)

# Graphing parameters, CAN BE CHANGED
resolution.figures <- 150
cex=1.3
cex.lab=1.2
cex.axis = 1

# Other parameters, CAN BE CHANGED
# Standard element correction
Standard.element.correction <- "yes"
Mo.correction.Cd <- "yes"

# The function, DON'T CHANGE THE CODE HERE
fBlank.Conc.Stats(blnk, sa.sf, st.el, mor, name.run, mainDir, resolution.figures, cex, cex.lab, cex.axis, Standard.element.correction, Mo.correction.Cd)

###########################################################################################
######### Checks & balances #############################################################

# To read in the results of the ICPMS for the checks and balances, DON'T CHANGE THE CODE HERE 
ck <- read.csv(paste(mainDir,"checks.csv",sep=""),header=TRUE)
head(ck)
table(ck$check)

# To read in the results of the seafast standard additions, DON'T CHANGE THE CODE HERE
sa.sf <- read.csv(paste(folder.sf,"StandAdd_seafast_results_",name.run,".csv", sep=""),header=TRUE)
head(sa.sf)

# To read in the results of the HNO3 standard additions, DON'T CHANGE THE CODE HERE
sa.HNO3 <- read.csv(paste(folder.HNO3,"StandAdd_HNO3_results_",name.run,".csv", sep=""),header=TRUE)
head(sa.HNO3)

# To read in the results of the standard.element, DON'T CHANGE THE CODE HERE
st.el <- read.csv(paste(folder.rh,standard.element," standard data_",name.run,".csv",sep=""),header=TRUE)
head(st.el)

# To read in the results of the seafast blank, ONLY CHANGE "blank.type.name" to the name that you gave to the blank that you will use
blank.type.name <- "MQ_blank"
blank <- read.csv(paste(folder.blank,"StatsData_blank_",blank.type.name,"_",name.run,".csv",sep=""),header=TRUE)
head(blank)

# To read in the reference material consensus values, DON'T CHANGE THE CODE HERE
ref.mat <- read.csv("TraceMetalReferenceConsensusValues.csv",header=TRUE)
head(ref.mat)

# To read in the results of Mo interference, DON'T CHANGE THE CODE HERE
mor <- read.csv(paste(folder.Mo,"Mo_interferenceCd_results_",name.run,".csv",sep=""),header=TRUE)
head(mor)

# To read the sample coding, DON'T CHANGE THE CODE HERE
coding <- read.csv(paste(mainDir,"sample_list.csv",sep=""),header=TRUE)
head(coding)

# Graphing parameters, CAN BE CHANGED
resolution.figures <- 150
cex=1.3
cex.lab=1.2
cex.axis = 1

# Other parameters, CAN BE CHANGED
Standard.element.correction <- "yes"
Mo.correction.Cd <- "yes"
use.Blank <- "median"		# either use "mean" or "median"
Seafast.conc.factor <- 10/0.65
outliers.boxplot <- "no"

# includes sample names in the graphs of the reference values, if you don't have a sample_list say "no"
sample.naming <- "yes"		

# the function, DON'T CHANGE THE CODE HERE
fChecks.balances(ck, sa.sf, sa.HNO3, st.el, blank, ref.mat, mor, coding, sample.naming, name.run, mainDir, use.Blank, Seafast.conc.factor, outliers.boxplot,resolution.figures, cex, cex.lab, cex.axis, Standard.element.correction,Mo.correction.Cd)

###########################################################################################
######### The final concentrations of the samples #########################################

# To read in the results of the ICPMS for the checks and balances, DON'T CHANGE THE CODE HERE 
fr <- read.csv(paste(mainDir,"samples.csv",sep=""),header=TRUE)
head(fr)

# To read in the results of the seafast standard additions, DON'T CHANGE THE CODE HERE
sa.sf <- read.csv(paste(folder.sf,"StandAdd_seafast_results_",name.run,".csv", sep=""),header=TRUE)
head(sa.sf)

# To read in the results of the standard.element, DON'T CHANGE THE CODE HERE
st.el <- read.csv(paste(folder.rh,standard.element," standard data_",name.run,".csv",sep=""),header=TRUE)
head(st.el)

# To read in the results of the seafast blank, DON'T CHANGE THE CODE HERE
blank.type.name <- "MQ_blank_SW"
blank <- read.csv(paste(folder.blank,"StatsData_blank_",blank.type.name,"_",name.run,".csv",sep=""),header=TRUE)
head(blank)

# To read in the results of Mo interference, DON'T CHANGE THE CODE HERE
mor <- read.csv(paste(folder.Mo,"Mo_interferenceCd_results_",name.run,".csv",sep=""),header=TRUE)
head(mor)

# Graphing parameters, CAN BE CHANGED
resolution.figures <- 150
cex=1.3
cex.lab=1.2
cex.axis = 1

# Other parameters, CAN BE CHANGED
Standard.element.correction <- "yes"
use.Blank <- "median"				# either use "mean" or "median"
Mo.correction.Cd <- "yes"

# the function, DON'T CHANGE THE CODE HERE
fsample.concentrations(fr, sa.sf, st.el, blank, mor, name.run, mainDir, use.Blank, resolution.figures, cex, cex.lab, cex.axis, Standard.element.correction,Mo.correction.Cd)

##################################################################################################################################
######## Function to include original sample names to the final concentrations results ###########################################

# To read the final concentration results, DON'T CHANGE THE CODE HERE
samples <- read.csv(paste(folder.samples,"Sample concentrations_",name.run,".csv",sep=""),header=TRUE)
head(samples)

# To read the sample coding, DON'T CHANGE THE CODE HERE
coding <- read.csv(paste(mainDir,"sample_list.csv",sep=""),header=TRUE)
head(coding)

# Function, DON'T CHANGE THE CODE HERE
fsample.codes(samples, coding)

##################################################################################################################################






