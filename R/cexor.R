
cexor <- function( bam, chrN, chrL, p=1e-9, dpeaks=c(0,150), dpairs=100, idr=0.01, N=5e6, bedfile=TRUE, mu=2.6, sigma=1.3, rho=0.8, prop=0.7 )
{
  options(digits=10)

  # Number of treatment samples
  NT <- length(bam)
  # Check chr length and chr names
  LL <- length(chrL)
  #LN <- length(chrN)

  # Read treatment M samples in BAM format (M biological replicates)
  listLen <- 0
  SampleCovPlus <- list()
  SampleCovMinus <- list()
  sumCovPlusMinus <- c()
  totalSum <- c()

  # Read BAM files
  for (j in seq(from=1, to=NT, by=1) ){

      for (i in 1:LL){
              listLen <-  listLen + 1
              which <- GRanges(seqnames = chrN[i], ranges = IRanges(c(1), c(chrL[i])))

              what <- c("rname", "strand", "pos", "qwidth","seq")
              param <- ScanBamParam(which = which, what = what)
              SampleBam <- scanBam(file=bam[j], param=param)
              sapply(SampleBam[[1]], class);

              # + strand
              SamplelstPlus <- lapply(names(SampleBam[[1]]), function(elt) { do.call(c, unname(lapply(SampleBam, "[[", elt)))})
              names(SamplelstPlus) <- names(SampleBam[[1]])
              SampledfPlus <- do.call("DataFrame", SamplelstPlus)
              # correction for reverse strands
              SampledfPlus$pos[which(SampledfPlus$strand=='+')] <- SampledfPlus$pos[which(SampledfPlus$strand=='+')] - 1
              # Modify read length (exonuclease stop sites, 1bp)
              SampledfPlus$qwidth <- 1
              # strand selection
              SampledfPlus <- SampledfPlus[which(SampledfPlus$strand=='+'),]
              SampleCovPlus[[listLen]]  <- coverage(IRanges(SampledfPlus[["pos"]], width = SampledfPlus[["qwidth"]]))
              if (length(SampleCovPlus[[listLen]])   < chrL[i] ) { SampleCovPlus[[listLen]]   <- c(SampleCovPlus[[listLen]],rep(0,chrL[i] - length(SampleCovPlus[[listLen]]) ) )  }
              if (length(SampleCovPlus[[listLen]])   > chrL[i] ) { SampleCovPlus[[listLen]]   <-  SampleCovPlus[[listLen]][1:chrL[i]]     }

              # - strand
              rm(SampledfPlus,SamplelstPlus)
              SamplelstMinus <- lapply(names(SampleBam[[1]]), function(elt) { do.call(c, unname(lapply(SampleBam, "[[", elt)))})
              names(SamplelstMinus) <- names(SampleBam[[1]])
              SampledfMinus  <- do.call("DataFrame", SamplelstMinus )
              # correction for reverse strands
              SampledfMinus$pos[which(SampledfMinus$strand=='-')] <- SampledfMinus$pos[which(SampledfMinus$strand=='-')] -1 + SampledfMinus$qwidth[which(SampledfMinus$strand=='-')]
              # Modify read length  (exonuclease stop sites, 1bp)
              SampledfMinus$qwidth <- 1
              # strand selection
              SampledfMinus <- SampledfMinus[which(SampledfMinus$strand=='-'),]
              SampleCovMinus[[listLen]]  <- coverage(IRanges(SampledfMinus[["pos"]], width = SampledfMinus[["qwidth"]]))
              if (length(SampleCovMinus[[listLen]])   < chrL[i] ) { SampleCovMinus[[listLen]]   <- c(SampleCovMinus[[listLen]],rep(0,chrL[i] - length(SampleCovMinus[[listLen]]) ) )  }
              if (length(SampleCovMinus[[listLen]])   > chrL[i] ) { SampleCovMinus[[listLen]]   <-  SampleCovMinus[[listLen]][1:chrL[i]]     }

             # All reads are considered (information from both strands is taken into account)

              sumCovPlusMinus[listLen]  <- sum(SampleCovPlus[[listLen]])  +  sum(SampleCovMinus[[listLen]])

              rm(SampleBam,SamplelstMinus,SampledfMinus,param,what,which)
      }

      # get total lambda exonuclease cuts at each sample
      totalSum[j] <- sum(sumCovPlusMinus[(((j-1)*LL) + 1):(j*LL)])
  }

 rm(i,j)
 rm(sumCovPlusMinus)
 
 
 # Normalization
 # To the smallest sample  - to make the samples 'comparable'
 # and estimation of Skellam Lambda1 and Lambda2 for each sample
 listLen <- 0
 meanCovPlus <- c()
 meanCovMinus <- c()
 lambda1 <- c()
 lambda2 <- c()

 factorNorm <-  min(totalSum) / totalSum

 for (j in 1:NT){

      for (i in 1:LL){
              listLen <-  listLen + 1
              SampleCovPlus[[listLen]]   <- round( factorNorm[j] * SampleCovPlus[[listLen]]   )      # it is an approximation
              SampleCovMinus[[listLen]]  <- round( factorNorm[j] * SampleCovMinus[[listLen]]  )      # it is an approximation
              meanCovPlus[listLen]  <- mean(SampleCovPlus[[listLen]])
              meanCovMinus[listLen] <- mean(SampleCovMinus[[listLen]])
      }
      # estimating Poisson means of the distributions
      lambda1_ini <- sum((chrL/sum(chrL))  * meanCovPlus[(((j-1)*LL) + 1):(j*LL)])
      lambda2_ini <- sum((chrL/sum(chrL))  * meanCovMinus[(((j-1)*LL) + 1):(j*LL)])

      # lambda1 > lambda2
      lambda1[j] <- max(lambda1_ini,lambda2_ini)
      lambda2[j] <- min(lambda1_ini,lambda2_ini)

      rm(lambda1_ini, lambda2_ini)

   }


  # Get paired peaks at each sample
  # using the Skellam cumulative distribution function
  pairedPeaks  <- list()
  listLen <- 0
  chrLwithPeaks <- list()

  for (w in 1:NT){
      	chrLwithPeaksCount <- 0        
      	filteredScoreFinalChr  <- data.frame()

  		for (i in 1:LL){  
      		filteredScoreFinal     <- data.frame() 
      		listLen <-  listLen + 1
      		mL <- max(length(SampleCovPlus[[listLen]]), length(SampleCovMinus[[listLen]]))
      		#print(LL)
      		#print(mL)
      		#frames <- floor(mL / N) +1 
      		frames <- floor(mL / N) #+1 ###NEW
      		#print(frames)
      		for (m in 1:frames){

         		if (m != frames){
            		lFrame <-  (((m-1)*N)+1):(m*N)
         		}
         		if (m == frames){
            		lFrame<- (((m-1)*N)+1):mL   
         		}
         	#print(m)
         	#print(lFrame)	
            	sample_forward <- as.integer(SampleCovPlus[[listLen]] [  lFrame ])
            	sample_reverse <- as.integer(SampleCovMinus[[listLen]][  lFrame ])

         		# Difference of exonuclease start sites at each strand
         		k <- as.integer(sample_forward - sample_reverse)

         		rm(sample_forward,sample_reverse)

         		# Skellam cumulative distribution function
         		# Calculate p-value (two-sided test)
         		score <- c()
         		iPlus  <- which(k>=0)   # '='can be at both
         		iMinus <- which(k<0)
         		score[iPlus]  <- pskellam(q=k[iPlus],    lambda1= lambda1[w], lambda2 = lambda2[w], lower.tail = FALSE, log.p = FALSE)
         		score[iMinus] <- pskellam(q=k[iMinus]-1, lambda1= lambda1[w], lambda2 = lambda2[w], lower.tail = TRUE, log.p = FALSE)

         		# p-value threshold
         		positions <- which(score <= p)   # score has peak p-values, F(x)

         		# which ones are close to each other
         		maxD <-  dpeaks[2]  #150 #30
         		minD <-  dpeaks[1]  #0   #8

         		filteredScore <- data.frame()
         		idx <- 0

         		if (length(positions)>1){
              		# keep just Max and Min in a region of length dpeaks[2]+1
              		W <- dpeaks[2]+1   #151
              		new_positions <- c()
              		for (ww in 1:length(positions)){
                		if ( positions[ww]-W < ceiling( W)  ){  vectorAux <- k[1:(positions[ww]+W) ]  }
                   		if ( positions[ww]-W >= ceiling( W)  &  positions[ww]+W <= N-floor(W)  ){  vectorAux <- k[(positions[ww]-W):(positions[ww]+W) ]  }
                   		if ( positions[ww]+W > N-floor(W)  ){  vectorAux <- k[(positions[ww]-W):N ]  }
                   		maximo <- max(vectorAux)
                   		minimo <- min(vectorAux)
  #               		print(ww)
  #                 		print(maximo)
  #                 		print(minimo)
                   		if( k[positions[ww]]==maximo && maximo >0 ){ new_positions <- c( new_positions,positions[ww]) }
                   		if( k[positions[ww]]==minimo && minimo <0 ){ new_positions <- c( new_positions,positions[ww]) }
              	}
        		positions <- new_positions
        		rm(new_positions)
   				if (length(positions) >= 2){ #######################NEW
        			for (j in 2:length(positions)){
                		A <- positions[j-1]
	                	B <- positions[j]
	                	C <- sign(k[positions[j-1]])
	                	D <- sign(k[positions[j]])
	                	# region for max and min
                		if (C != 0  && D != 0 ){
	             			if (B< A+maxD && B> A+minD  && (C == 1  && D== -1 || C== -1  && D == 1)  ){
		            			idx=idx+1
		                   		filteredScore[idx,1] <- chrN[i]
		                   		filteredScore[idx,2] <- round((lFrame[positions[j-1]] + lFrame[positions[j]])/2) -dpairs
		                   		filteredScore[idx,3] <- round((lFrame[positions[j-1]] + lFrame[positions[j]])/2) +dpairs
		                   		filteredScore[idx,4] <- lFrame[positions[j]]-lFrame[positions[j-1]]
		                   		filteredScore[idx,5] <- signif(score[positions[j-1]], digits=30)
		                   		filteredScore[idx,6] <- signif(score[positions[j]] , digits=30)
				           	pvalVector <- c(score[positions[j-1]] , score[positions[j]])
				   		filteredScore[idx,7] <- sum(pvalVector)  #pchisq(-2*sum(log(pvalVector)), df=2*length(pvalVector), lower=FALSE )
			           	   	filteredScore[idx,8] <- -1*log10(filteredScore[idx,7])
			        		}
		        	 	}
        	 		}
 				}########################NEW             
 		}

		if ( length(filteredScoreFinal) ==0  &&  length(filteredScore) > 0) { filteredScoreFinal <- filteredScore  }
                   
			if ( length(filteredScoreFinal) > 0 &&  length(filteredScore) > 0) { filteredScoreFinal <- rbind(filteredScoreFinal,filteredScore)}

    }     #end frame
    
    

   if ( length(filteredScoreFinal) > 0  ){                                                
        if (length(filteredScoreFinalChr) ==0 ) {  filteredScoreFinalChr <- filteredScoreFinal   }
    	if (length(filteredScoreFinalChr) > 0 ) {  filteredScoreFinalChr <- rbind(filteredScoreFinalChr, filteredScoreFinal )  }
        #save the chromosome
        chrLwithPeaksCount <- chrLwithPeaksCount + 1                                        
        if(chrLwithPeaksCount == 1) { chrLwithPeaks[[w]] <- chrL[i]  }                      
        if(chrLwithPeaksCount > 1)  { chrLwithPeaks[[w]][chrLwithPeaksCount] <- chrL[i]  }  
	}                                                                                      
}         #end chromosome
                                     

 colnames(filteredScoreFinalChr) <- c("chr","start","end","length","pval.peak.forward","pval.peak.reverse", "pval.binding.event","log10.pval.binding.event")
 pairedPeaks[[w]] <- filteredScoreFinalChr

 rm(filteredScoreFinalChr)
}

#print(head(filteredScoreFinalChr))



##############################################################################

#create a GRanges structure for each biol. replicate
repl <- list()

for (w in 1:NT){

      repl[[w]] <- GRanges(seqnames=pairedPeaks[[w]]$chr,
                  ranges=IRanges(pairedPeaks[[w]]$start, end=pairedPeaks[[w]]$end),
                  score=pairedPeaks[[w]]$log10.pval.binding.event)

       seqlengths(repl[[w]])   <- chrLwithPeaks[[w]]    ##chrL
       
       repl[[w]] <- trim(repl[[w]], use.names=TRUE)
}


# (intersect) final list of chip-exo peaks
finalset  <- repl[[1]]

for (w in 2:NT){
         finalsetini <-  finalset
         finalset <- intersect(repl[[w]],finalsetini )

}

       finalset <- trim(finalset, use.names=TRUE)
       

# prepare matrix MM for IDR assessment
MM <- matrix(data = NA, nrow = length(finalset), ncol = NT)

for (w in 1:NT){
    for (i in 1: length(finalset) ){
      MM[i,w] <-   max(score( subsetByOverlaps(repl[[w]],finalset[i]) ))   #MIN is very conservative
    }
}

# send matrix MM to IDR analysis
#  mu <- 2.6
#  sigma <- 1.3
#  rho <- 0.8
#  p <- 0.7

 # library(idr)
 # x = a m by n numeric matrix, where m= num of replicates, n=num of observations. Numerical values representing the significance
 # of the observations, where signals are expected to have large values, for example, -log(p-value).
 idr.out <- est.IDR(x=MM, mu, sigma, rho, prop, eps=0.00001)
 #idr.out$IDR

# FINAL  sites
# values(finalset) <- idr.out$IDR
 passingIDR <- which(idr.out$IDR< idr)

# prepare final GRanges
PVALMATRIX=10^-MM

MM2 <- as.data.frame(MM)
colnames(MM2) <- paste(paste("rep", as.character(1:NT), sep=""), "neg.log10pvalue",sep=".")
# get p-values



#http://en.wikipedia.org/wiki/Fisher%27s_method
Stouffer.test <- function(p, w) { # p is a vector of p-values
  if (missing(w)) {
    w <- rep(1, length(p))/length(p)
  } else {
    if (length(w) != length(p))
      stop("Length of p and w must equal!")
  }
  Zi <- qnorm(1-p) 
  Z  <- sum(w*Zi)/sqrt(sum(w^2))
  p.val <- 1-pnorm(Z)
  return(c(p.value = p.val))
}
#http://en.wikipedia.org/wiki/Fisher%27s_method
Fisher.test <- function(p) {
  Xsq <- -2*sum(log(p))
  p.val <- 1-pchisq(Xsq, df = 2*length(p))
  return(c(p.value = p.val))
}

Stouffer <-c()
Fisher <-c()
for (i in 1:dim(PVALMATRIX)[1]) {

	Stouffer[i] <- Stouffer.test(p=PVALMATRIX[i,])[[1]]
	Fisher[i] <- Fisher.test(p=PVALMATRIX[i,])[[1]]
 

}

# Fisher's
 mcols(finalset) <- data.frame(IDR=idr.out$IDR, MM2, Stouffer.pvalue=Stouffer, Fisher.pvalue=Fisher)
 

###############################################################
 #centres
 finalsetcentres <- finalset

 start(finalsetcentres) <- round( (start(finalset) + end(finalset) )/2 )
 end(finalsetcentres)  <-  round( (start(finalset) + end(finalset) )/2 ) + 1

 if (bedfile==TRUE){
    export.bedGraph(object=finalset[passingIDR], con=paste("cexor_peak_events_IDR_",".bed", sep=as.character(idr)))
    export.bedGraph(object=finalsetcentres[passingIDR], con=paste("cexor_peak_centres_IDR_",".bed", sep=as.character(idr)))
 }

  return(list( bindingEvents=finalset[passingIDR], bindingCentres =finalsetcentres[passingIDR],  pairedPeaksRepl = repl ))

}

