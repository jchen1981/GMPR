# Demonstrate the use of GMPR factor
#setwd('path2GMPR')

require(GUniFrac)
require(vegan)
require(DESeq2)
source('GMPR.R')

data(throat.otu.tab)
data(throat.meta)

###########################################################################################################
# Calculate GMPR size factor
# Row - features, column - samples
otu.tab <- t(throat.otu.tab)
gmpr.size.factor <- GMPR(otu.tab)$gmpr
###########################################################################################################

# Two potential applications of GMPR size factors

###########################################################################################################
# Application 1: Count are normalized by size factors to reduce the variation due to different library sizes
# The normalized counts are subject to further downstream analysis such as ordination (PCA, PCoA), clustering,
# and other multivariate methods. Note that further data transformation such as VST transformation (DESeq2)
# may be needed in order to reveal patterns. Here shows an example of BC distance based ordination 
otu.tab.norm <- t(t(otu.tab) / gmpr.size.factor)
dist.mat <- vegdist(t(otu.tab.norm))
PCs <- cmdscale(dist.mat, k=2)
plot(PCs[, 1], PCs[, 2], col=factor(throat.meta$SmokingStatus))
###########################################################################################################

###########################################################################################################
# Application 2:  Differential abundance analysis using DESeq2 using GMPR size factors instead of the default
dds <- DESeqDataSetFromMatrix(countData = otu.tab,
		colData = throat.meta,
		design= ~ SmokingStatus)

# Replace size factor
# dds <- estimateSizeFactors(dds)
sizeFactors(dds) <-  gmpr.size.factor
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
results(dds, contrast=c("SmokingStatus", "NonSmoker", "Smoker"))
###########################################################################################################


