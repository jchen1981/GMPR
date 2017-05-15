# Title: Geometric Mean of Pairwise Ratios (GMPR) for Microbiome Sequencing data normalization
# Version: 0.1
# Authors: Jun Chen (chen.jun2@mayo.edu)
# Date: 2017/02/07
# Description: The function calculates the normalizing factors for microbiome sequencing data or generally zeroinflated sequencing data. 
# The size factors can be used as offsets in count-based regression models or as devisors to produce normalized data


require(matrixStats)

GMPR <- function (comm, intersect.no=4, ct.min=5) {
	# Computes the GMPR size factor
	#
	# Args:
	#   comm: a matrix of counts, row - features (OTUs, genes, etc) , column - sample
	#   intersect.no: the minimum number of shared features between sample pair, where the ratio is calculated
	#   ct.min: the minimum number of counts required to calculate ratios （Empirical study found ct.min=5 is suitable)

	#
	# Returns:
	#   a list that contains:
	#      gmpr： the GMPR size factors for all samples; Samples with distinct sets of features will be output as NA.
	#      nss:   number of samples with significant sharing (> intersect.no) including itself
	
	# mask counts < ct.min
	comm[comm < ct.min] <- 0
	
	if (is.null(colnames(comm))) {
		colnames(comm) <- paste0('S', 1:ncol(comm))
	}
	
	cat('Begin GMPR size factor calculation ...\n')
	
	comm.no <- numeric(ncol(comm))
	gmpr <- sapply(1:ncol(comm),  function(i) {		
				if (i %% 50 == 0) {
					cat(i, '\n')
				}
				x <- comm[, i]
				# Compute the pairwise ratio
				pr <- x / comm
				# Handling of the NA, NaN, Inf
				pr[is.nan(pr) | !is.finite(pr) | pr == 0] <- NA
				# Counting the number of non-NA, NaN, Inf
				incl.no <- colSums(!is.na(pr))		
				# Calculate the median of PR
				pr.median <- colMedians(pr, na.rm=TRUE)
				# Record the number of samples used for calculating the GMPR
				comm.no[i] <<- sum(incl.no >= intersect.no)
				# Geometric mean of PR median
				if (comm.no[i] > 1) {
					return(exp(mean(log(pr.median[incl.no >= intersect.no]))))
				} else {
					return(NA)
				}
			}
	)
	
	if (sum(is.na(gmpr))) {
		warning(paste0('The following samples\n ', paste(colnames(comm)[is.na(gmpr)], collapse='\n'), 
				'\ndo not share at least ', intersect.no, ' common taxa with the rest samples! ',
				'For these samples, their size factors are set to be NA! \n', 
				'You may consider removing these samples since they are potentially outliers or negative controls!\n',
				'You may also consider decreasing the minimum number of intersecting taxa and rerun the procedure!\n'))
	}
	
	cat('Completed!\n')
	cat('Please watch for the samples with limited sharing with other samples based on NSS! They may be outliers! \n')
	names(gmpr) <- names(comm.no) <- colnames(comm)
	return(list(gmpr=gmpr, nss=comm.no))
}
