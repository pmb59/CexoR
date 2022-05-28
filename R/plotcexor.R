plotcexor <- function( bam, peaks, EXT=500 ){  
  
  N <- length(bam)
  colo <- colorRampPalette(brewer.pal(N,"Set2"))(N)
  
  # extend peak centre Upstream and Dowstream EXT bp
  peaks <- peaks$bindingCentres 
  start(peaks) <- start(peaks) - EXT
  end(peaks)   <- end(peaks)   + EXT
  # genomation
  nBins <- (2*EXT)+1
  scaleData <- FALSE
  # read bam
  sm <- list()
  smF <- list()
  smR <- list()
  for (j in seq(from=1, to=N, by=1) ){
    sm[[j]]  <- ScoreMatrixBin(target = bam[j], bin.num = nBins, windows = peaks, type="bam",strand.aware = TRUE, extend=0,rpm=F, param = ScanBamParam(which=reduce(peaks, ignore.strand=TRUE), flag=scanBamFlag(isMinusStrand=NA)))
    smF[[j]] <- ScoreMatrixBin(target = bam[j], bin.num = nBins, windows = peaks, type="bam",strand.aware = TRUE, extend=1,rpm=F, param = ScanBamParam(which=reduce(peaks, ignore.strand=TRUE), flag=scanBamFlag(isMinusStrand=FALSE)))
    smR[[j]] <- ScoreMatrixBin(target = bam[j], bin.num = nBins, windows = peaks, type="bam",strand.aware = TRUE, extend=1,rpm=F, param = ScanBamParam(which=reduce(peaks, ignore.strand=TRUE), flag=scanBamFlag(isMinusStrand=TRUE)))
  }

  # plot Ylim 
  YMAX <-0; YMAXF <-0; YMAXR <-0;
  for (j in 1:N){
    YMAX  <- max(YMAX, colMeans(sm[[j]],na.rm=TRUE) )
    YMAXF <- max(YMAXF, colMeans(smF[[j]],na.rm=TRUE) )
    YMAXR <- max(YMAXR, colMeans(smR[[j]],na.rm=TRUE) )
  }
  
  par(mfrow=c(2,1))
 
  L <- -EXT:EXT
  
  # upper panel
  for (j in seq(from=1, to=N, by=1) ){
    ym <- max(YMAXF,YMAXR)
    if (j==1){
      plot(L,colMeans(smF[[j]],na.rm=TRUE),type="l", col="red", lwd=1, frame.plot=F, xlab="Distance to \nbinding centre (bp)", ylab=expression(paste("Average ",lambda," exonuclease stop sites")), ylim=c(0,ym), lty=j)
      points(L,colMeans(smR[[j]],na.rm=TRUE),type="l", col="blue", lwd=1, lty=j) 
    }
    if (j!=1){
      points(L,colMeans(smF[[j]],na.rm=TRUE),type="l", col="red", lwd=1, lty=j) 
      points(L,colMeans(smR[[j]],na.rm=TRUE),type="l", col="blue", lwd=1, lty=j)   
    }
    
  }
  abline(v=0, lty=2)
  legend("topleft", text.col=c("red","blue"), legend=c("+ strand","- strand"), bty="n")
  legend("topright", col="black", legend=paste("rep", 1:N), lty=1:N, bty="n", lwd=1)
  
  # lower panel
  for (j in seq(from=1, to=N, by=1) ){
    if (j==1) plot(L,colMeans(sm[[j]],na.rm=TRUE),type="l", col=colo[j], lwd=2, frame.plot=F, xlab="Distance to \nbinding centre (bp)", ylab="Average ChIP-exo reads", ylim=c(0,YMAX))
    if (j!=1) points(L,colMeans(sm[[j]],na.rm=TRUE),type="l", col=colo[j], lwd=2)
  }
  abline(v=0, lty=2)
  legend("topright", col=colo, legend=paste("rep", 1:N), bty="n", lty=1, lwd=2)
  
}
