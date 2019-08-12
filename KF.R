
rm(list = ls())

# libraries for compiling C++ codes used in this work
library(Rcpp) 
library(RcppArmadillo)

# library for checking positive semi-definite
library(matrixcalc)
 
# libraries for calAUPR
library(MESS)
library(pracma)
library(ROCR) 
library(Bolstad2)

# sourceCpp
sourceCpp("fastKgipMat.cpp")
sourceCpp("fastKF.cpp")
sourceCpp("fastSolve.cpp")

# source
source("RLS_KF.R")
#source("calAUPR.R")

###################################################################################################
# You just modify the partfn to different data sets
# file name to be used: nr, gpcr, ic, e
#partfn = "nr"
# Take long time for Enzyme dataset
#if (partfn == "e") cat("Need several hours to finish the big [Enzymes] data set, please be patient!\n")
###################################################################################################


	  yFn <- paste0( "adj.txt")
	  #yFn <- paste0("case3_adj.txt")
	  
	  y <- read.table(yFn)
	  y=t(y)
	  # simmatCompd
	  #simCompdFn <- paste0(partfn, "_simmat_dc.txt")
	  simCompdFn <- paste0( "miRNA_sim.txt")
	  #simCompdFn <- paste0("case3_miRNAsimilarity.txt")
	  
	  #simCompdFn <- paste0( "dis_sim.txt")
	  simmatCompd <- read.table(simCompdFn)
	  # simmatTarget
	  #simTargetFn <- paste0(partfn, "_simmat_dg.txt")
	  simTargetFn <- paste0( "dis_sim.txt")
	  #simTargetFn <- paste0("case3_diseasesimilarity.txt")
	  
	  #simTargetFn <- paste0( "miRNA_sim.txt")
	  simmatTarget <- read.table(simTargetFn)
	  # convert into matrix
	  y <- as.matrix(y)
	  simmatCompd <- as.matrix(simmatCompd)
	  simmatTarget <- as.matrix(simmatTarget)
	  # check matrix symmetric
	  if (!isSymmetric(simmatCompd)) simmatCompd <- (simmatCompd + t(simmatCompd))/2
	  # check matrix positive semi-definite 
	  epsilon <- 0.1
	  while (!is.positive.semi.definite(simmatCompd)) simmatCompd <- simmatCompd + epsilon * diag(nrow(simmatCompd))
	  # check matrix symmetric
	  if (!isSymmetric(simmatTarget)) simmatTarget <- (simmatTarget + t(simmatTarget))/2
	  # check matrix positive semi-definite 
	  epsilon <- 0.1
	  while (!is.positive.semi.definite(simmatTarget)) simmatTarget <- simmatTarget + epsilon * diag(nrow(simmatTarget))
	
# parameters
gamma0   = 1    
lambda   = 1    
nfold    = 5
# k for KF
numNeig  = 4    
# t for KF
numIter  = 2
	flush.console()

	# (1) Prediction based on the target similarity
#y=t(y)
numRows2 <- nrow(y)
numCols2 <- ncol(y)
#row miRNA
cat(numRows2,"\n")
#while(1)
id=0

#for (ii in 1:numRows2){
#	for(jj in 1:numCols2){
#	if(y[ii,jj]==0) next
#	y[ii,jj]=0
	#cat("i ",ii," j ",jj ,"\n")
	#a=c(50,91,92,126,145,159,59,236,240,252,254,260,327,205)
	#for(ii in 1:14){
	#temp=a[ii]
	#temp_y=y
	#y[temp,1:495]=0
	id=id+1
	cat(id,"\n")
	yCompd <- t(y)
	ytrMat <- yCompd
	k4simmat <- simmatCompd
	Kgip <- fastKgipMat(ytrMat, gamma0)
  # Kernel Fusion (KF)
	K <- fastKF(Kgip, k4simmat, numNeig, numIter)
	
	ytrMat <- y
	k4simmat <- simmatTarget
	Kgip <- fastKgipMat(ytrMat, gamma0)
  # Kernel Fusion (KF)
	K2 <- fastKF(Kgip, k4simmat, numNeig, numIter)
	
	#currY <- y[,jj]
	#currYsum <- sum(ytrMat[, i])
	#write.table(K,paste("K_",id,".txt",sep=""),quote = FALSE,row.names = FALSE,
    #           col.names = FALSE)
	#write.table(K2,paste("K2_",id,".txt",sep=""),quote = FALSE,row.names = FALSE,
    #           col.names = FALSE)
    write.table(K,paste("K.txt",sep=""),quote = FALSE,row.names = FALSE,
               col.names = FALSE)
	write.table(K2,paste("K2.txt",sep=""),quote = FALSE,row.names = FALSE,
               col.names = FALSE)
	#write.table(K,paste("case2_K_",id,".txt",sep=""),quote = FALSE,row.names = FALSE,
    #           col.names = FALSE)
	#write.table(K2,paste("case2_K2_",id,".txt",sep=""),quote = FALSE,row.names = FALSE,
    #           col.names = FALSE)
    #write.table(K,paste("case3_K.txt",sep=""),quote = FALSE,row.names = FALSE,
    #           col.names = FALSE)
	#write.table(K2,paste("case3_K2.txt",sep=""),quote = FALSE,row.names = FALSE,
    #           col.names = FALSE)
    #y=temp_y
#}
#	y[ii,jj]=1
#}
#}