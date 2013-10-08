\name{CexoR-package}
\alias{CexoR-package}
\alias{CexoR}
\docType{package}

\title{ An R package to uncover high-resolution protein-DNA interactions in ChIP-exo replicates }
\description{ Strand specific peak-pair calling in ChIP-exo replicates. The cumulative Skellam distribution function (package 'skellam') is used to detect significant normalized count differences of opposed sign at each DNA strand (peak-pairs). Irreproducible discovery rate for overlapping peak-pairs across biological replicates is estimated using the package 'idr'.}
\details{
\tabular{ll}{
Package: \tab CexoR\cr
Type: \tab Package\cr
Version: \tab 0.99.5\cr
Date: \tab 2013-09-09\cr
License: \tab Artistic-2.0 | GPL-2 + file LICENSE\cr
LazyLoad: \tab yes\cr
}
}

\author{
Pedro Madrigal,

Maintainer: Pedro Madrigal \email{pm@engineering.com}

}
\references{ Madrigal P, et al. (in preparation).\cr 
Skellam JG (1946) The frequency distribution of the difference between two Poisson variates belonging to different populations. J R Stat Soc Ser A 109: 296. \cr 
Li Q, Brown J, Huang H, Bickel P (2011) Measuring reproducibility of high-throughput experiments. Ann Appl Stat 5: 1752-1779. \cr 
Rhee HS, Pugh BF (2011) Comprehensive genome-wide protein-DNA interactions detected at single-nucleotide resolution. Cell 147: 1408-1419. }

\examples{

## hg19. chr2:1-1,000,000. CTCF data from Rhee and Pugh (2011)

owd <- setwd(tempdir())

rep1 <- "CTCF_rep1_chr2_1-1e6.bam"
rep2 <- "CTCF_rep2_chr2_1-1e6.bam"
rep3 <- "CTCF_rep3_chr2_1-1e6.bam"
r1 <- system.file("extdata", rep1, package="CexoR",mustWork = TRUE)
r2 <- system.file("extdata", rep2, package="CexoR",mustWork = TRUE)
r3 <- system.file("extdata", rep3, package="CexoR",mustWork = TRUE)

cexor(bam=c(r1,r2,r3), chrN="chr2", chrL=1e6, idr=0.01, p=1e-12, N=3e4)

setwd(owd)
}
