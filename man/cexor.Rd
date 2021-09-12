\name{cexor}
\alias{cexor}

\title{ ChIP-exo peak-pair calling with replicates}
\description{ ChIP-exo peak-pair calling with replicates. }

\usage{
cexor(bam, chrN, chrL, p=1e-9, dpeaks=c(0,150), dpairs=100, idr=0.01, 
N=5e6, bedfile=TRUE, mu=2.6, sigma=1.3, rho=0.8, prop=0.7)
}

\arguments{
  \item{bam}{ BAM alignment files of biological replicates. }
  \item{chrN}{ Vector of chromosome names. }
  \item{chrL}{ Vector of chromosome sizes (bp). }
  \item{p}{ P-value cutoff (should be relaxed, e.g. 1e-3, to allow the correct estimation of the irreproducible discovery rate (idr). However, this depends on the sequencing depth. For datasets with high number of tag counts, 1e-9 can be appropriate. See the vignette for more information.) }
  \item{dpeaks}{ Min. and max. allowed distance between peak pairs located at opposed strands in a replicate (bp). }
  \item{dpairs}{ Max. allowable distance between peak-pair centres across replicates (bp). }
  \item{idr}{ Irreproducible discovery rate cutoff [0-1]. }
  \item{N}{ Genome is divided in blocks of N bp. for processing. N must be not higher than the size of the smallest chromosome. }
  \item{bedfile}{ Generate BED files of ChIP-exo reproducible peak pairs.}
  \item{mu}{ A starting value for the mean of the reproducible component (see 'idr' package).}
  \item{sigma}{ A starting value for the standard deviation of the reproducible component (see 'idr' package).}
  \item{rho}{ A starting value for the correlation coefficient of the reproducible component (see 'idr' package).}
  \item{prop}{ A starting value for the proportion of reproducible component (see 'idr' package).}

  
}

\details{ Strand specific peak-pair calling in ChIP-exo replicates. The cumulative Skellam distribution function (package 'skellam') is used to detect significant normalized count differences of opposed sign at each DNA strand (peak-pairs). Irreproducible discovery rate for overlapping peak-pairs across biological replicates is estimated using the package 'idr'. 
The internal functions \code{pskellam} and \code{pskellam.sp} from the Jerry W. Lewis' 'skellam' R package (version 0.0-8-7) are used to calculate the cumulative Skellam distribution (see LICENSE file).}

\value{
  A list containing the following elements:
  \item{bindingEvents}{ A GRanges object with reproducible peak pair locations. The metadata 'value' indicates the Irreproducible discovery rate (IDR) estimated at this region, while 'repI.neg.log10pvalue' indicates -log10(p-value) for the replicate I. 'Stouffer.pvalue' and 'Fisher.pvalue' report the combined p-value considering they come from from independent significance tests. }
  \item{bindingCentres}{ A GRanges object with centre position of reproducible peak pair locations. The metadata 'value' indicates the Irreproducible discovery rate (IDR) estimated at this region, while 'repI.neg.log10pvalue' indicates -log10(p-value) for the replicate I. 'Stouffer.pvalue' and 'Fisher.pvalue' report the combined p-value considering they come from from independent significance tests.}
  \item{pairedPeaksRepl}{ A GRangesList object with the location of peak pairs retrieved at each replicate.  The metadata 'score' indicates -log10(p-value).}
  
}

\references{ Madrigal P (2015) CexoR: an R/Bioconductor package to uncover high-resolution protein-DNA interactions in ChIP-exo replicates. EMBnet.journal  21: e837.  }

\author{ Pedro Madrigal, \email{pmadrigal@ebi.ac.uk} }

\seealso{ \code{\link{CexoR-package}} }

\examples{

## hg19. chr2:1-1,000,000. CTCF data from Rhee and Pugh (2011)

owd <- setwd(tempdir())

rep1 <- "CTCF_rep1_chr2_1-1e6.bam"
rep2 <- "CTCF_rep2_chr2_1-1e6.bam"
rep3 <- "CTCF_rep3_chr2_1-1e6.bam"
r1 <- system.file("extdata", rep1, package="CexoR",mustWork = TRUE)
r2 <- system.file("extdata", rep2, package="CexoR",mustWork = TRUE)
r3 <- system.file("extdata", rep3, package="CexoR",mustWork = TRUE)

chipexo <- cexor(bam=c(r1,r2,r3), chrN="chr2", chrL=1e6, idr=0.01, p=1e-12, N=3e4)

plotcexor(bam=c(r1,r2,r3), peaks=chipexo, EXT=500)

setwd(owd)

}
