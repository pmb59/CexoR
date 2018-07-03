test_cexor <- function() {

	owd <- setwd(tempdir())

	rep1 <- "CTCF_rep1_chr2_1-1e6.bam"
	rep2 <- "CTCF_rep2_chr2_1-1e6.bam"
	rep3 <- "CTCF_rep3_chr2_1-1e6.bam"

	r1 <- system.file("extdata", rep1, package="CexoR",mustWork = TRUE)
	r2 <- system.file("extdata", rep2, package="CexoR",mustWork = TRUE)
	r3 <- system.file("extdata", rep3, package="CexoR",mustWork = TRUE)
	
	out <- cexor(bam=c(r1,r2,r3), chrN="chr2", chrL=1e6, N=3e4, p=1e-12)

	#check that the example data outputs a GRanges with 13 ranges 
    checkEquals(length(out$bindingEvents), 13)

        #check that if minimal input is not introduced the package recognizes an error situation
	#bam
	checkException(cexor(bam=c(), chrN="chr2", chrL=1e6, N=3e4))
	#chrN
	checkException(cexor(bam=c(r1,r2,r3), chrN=c(), chrL=1e6, N=3e4))
	#chrL
	checkException(cexor(bam=c(r1,r2,r3), chrN="chr2", chrL=c(), N=3e4))

	setwd(owd)

}
