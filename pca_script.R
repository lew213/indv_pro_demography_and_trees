#===============================================================================
#                               Welcome!
# This adgenet script was provided by Levi Yant and altered by our group (R3G2)
# From an initial VCF file (from GATK) it does a PCA, K-means clustering
# and also calculates a matrix of genetic distances and does an AMOVA.
#
# It was written for plants so it can deal with different ploidies (useful!)
#
# If a comment doesn't have name, it was made by Levi, otherwise it has the
# name of the person that has done it between hashtags
#
# Script by: Filip Kolar 2017, further edits by Sian Bray 2018 & Levi Yant 2022
# altered by: Lewis Wood 2022
# Date: july 2022
#===============================================================================
rm(list = ls())

#installation of packages
install.packages("adegenet", dep=TRUE)
install.packages("StAMPP")
install.packages("phangorn")
install.packages("BiocManager")
BiocManager::install("msa")
#print warnings as they occur
options(warn=1)

#call libraries:
library(adegenet)
library(StAMPP)
library(vcfR)
library(ggplot2)
library(MASS)
library(adegraphics)
library(ape)
library(phangorn)
library(msa)
######################=========MODIFIED FUNCTIONS=========######################

# a function for conversion from vcfR object to genlight in tetraploids
vcfR2genlight.tetra <- function (x, n.cores = 1) 
{
  bi <- is.biallelic(x)
  if (sum(!bi) > 0) {
    msg <- paste("Found", sum(!bi), "loci with more than two alleles.")
    msg <- c(msg, "\n", paste("Objects of class genlight only support loci with two alleles."))
    msg <- c(msg, "\n", paste(sum(!bi), "loci will be omitted from the genlight object."))
    warning(msg)
    x <- x[bi, ]
  }
  x <- addID(x)
  CHROM <- x@fix[, "CHROM"]
  POS <- x@fix[, "POS"]
  ID <- x@fix[, "ID"]
  x <- extract.gt(x)
  x[x == "0|0"] <- 0    
  x[x == "0|1"] <- 1    
  x[x == "1|0"] <- 1     
  x[x == "1|1"] <- 2     
  x[x == "0/0"] <- 0     
  x[x == "0/1"] <- 1     
  x[x == "1/0"] <- 1    
  x[x == "1/1"] <- 2    
  x[x == "1/1/1/1"] <- 4
  x[x == "0/1/1/1"] <- 3
  x[x == "0/0/1/1"] <- 2
  x[x == "0/0/0/1"] <- 1
  x[x == "0/0/0/0"] <- 0
  x[x == "0/0/0/0/0/0"] <- 0
  x[x == "0/0/0/0/0/1"] <- 1
  x[x == "0/0/0/0/1/1"] <- 2
  x[x == "0/0/0/1/1/1"] <- 3
  x[x == "0/0/1/1/1/1"] <- 4
  x[x == "0/1/1/1/1/1"] <- 5
  x[x == "1/1/1/1/1/1"] <- 6
  if (requireNamespace("adegenet")) {
    x <- new("genlight", t(x), n.cores = n.cores)
  }
  else {
    warning("adegenet not installed")
  }
  adegenet::chromosome(x) <- CHROM
  adegenet::position(x) <- POS
  adegenet::locNames(x) <- ID
  return(x)
}
# a patch for MUCH MUCH faster PCA calculation on genlight objects
# see https://github.com/thibautjombart/adegenet/pull/150
glPcaFast <- function(x,
                      center=TRUE,
                      scale=FALSE,
                      nf=NULL,
                      loadings=TRUE,
                      alleleAsUnit=FALSE,
                      returnDotProd=FALSE){
  
  if(!inherits(x, "genlight")) stop("x is not a genlight object")
  # keep the original mean / var code, as it's used further down
  # and has some NA checks...
  if(center) {
    vecMeans <- glMean(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecMeans))) stop("NAs detected in the vector of means")
  }
  if(scale){
    vecVar <- glVar(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecVar))) stop("NAs detected in the vector of variances")
  }
  # convert to full data, try to keep the NA handling as similar to the original as possible
  # - dividing by ploidy keeps the NAs
  mx <- t(sapply(x$gen, as.integer)) / ploidy(x)
  # handle NAs
  NAidx <- which(is.na(mx), arr.ind = T)
  if (center) {
    mx[NAidx] <- vecMeans[NAidx[,2]]
  } else {
    mx[NAidx] <- 0
  }
  # center and scale
  mx <- scale(mx,
              center = if (center) vecMeans else F,
              scale = if (scale) vecVar else F)
  # all dot products at once using underlying BLAS to support thousands of samples,
  # this could be replaced by 'Truncated SVD', but it would require more changes
  # in the code around
  allProd <- tcrossprod(mx) / nInd(x) # assume uniform weights
  
  ## PERFORM THE ANALYSIS ## ---------------------------------------------------
  # eigen analysis
  eigRes <- eigen(allProd, symmetric=TRUE, only.values=FALSE)
  rank <- sum(eigRes$values > 1e-12)
  eigRes$values <- eigRes$values[1:rank]
  eigRes$vectors <- eigRes$vectors[, 1:rank, drop=FALSE]
  # scan nb of axes retained
  if(is.null(nf)){
    barplot(eigRes$values, main="Eigenvalues", col=heat.colors(rank))
    cat("Select the number of axes: ")
    nf <- as.integer(readLines(n = 1))
  }
  # rescale PCs
  res <- list()
  res$eig <- eigRes$values
  nf <- min(nf, sum(res$eig>1e-10))
  ##res$matprod <- allProd # for debugging
  ## use: li = XQU = V\Lambda^(1/2)
  eigRes$vectors <- eigRes$vectors * sqrt(nInd(x)) # D-normalize vectors
  res$scores <- sweep(eigRes$vectors[, 1:nf, drop=FALSE],2, sqrt(eigRes$values[1:nf]), FUN="*")
  
  ## GET LOADINGS ## -----------------------------------------------------------
  # need to decompose X^TDV into a sum of n matrices of dim p*r
  # but only two such matrices are represented at a time
  if(loadings){
    if(scale) {
      vecSd <- sqrt(vecVar)
    }
    res$loadings <- matrix(0, nrow=nLoc(x), ncol=nf) # create empty matrix
    ## use: c1 = X^TDV
    ## and X^TV = A_1 + ... + A_n
    ## with A_k = X_[k-]^T v[k-]
    myPloidy <- ploidy(x)
    for(k in 1:nInd(x)){
      temp <- as.integer(x@gen[[k]]) / myPloidy[k]
      if(center) {
        temp[is.na(temp)] <- vecMeans[is.na(temp)]
        temp <- temp - vecMeans
      } else {
        temp[is.na(temp)] <- 0
      }
      if(scale){
        temp <- temp/vecSd
      }
      res$loadings <- res$loadings + matrix(temp) %*% eigRes$vectors[k, 1:nf, drop=FALSE]
    }
    res$loadings <- res$loadings / nInd(x) # don't forget the /n of X_tDV
    res$loadings <- sweep(res$loadings, 2, sqrt(eigRes$values[1:nf]), FUN="/")
  }
  ## FORMAT OUTPUT ## ----------------------------------------------------------
  colnames(res$scores) <- paste("PC", 1:nf, sep="")
  if(!is.null(indNames(x))){
    rownames(res$scores) <- indNames(x)
  } else {
    rownames(res$scores) <- 1:nInd(x)
  }
  if(!is.null(res$loadings)){
    colnames(res$loadings) <- paste("Axis", 1:nf, sep="")
    if(!is.null(locNames(x)) & !is.null(alleles(x))){
      rownames(res$loadings) <- paste(locNames(x),alleles(x), sep=".")
    } else {
      rownames(res$loadings) <- 1:nLoc(x)
    }
  }
  if(returnDotProd){
    res$dotProd <- allProd
    rownames(res$dotProd) <- colnames(res$dotProd) <- indNames(x)
  }
  res$call <- match.call()
  class(res) <- "glPca"
  return(res)
}
######################====================================######################

# IMPORT SNP data from VCF
filename <- readline(prompt="Enter file full path and name file (example ~/Desktop/intermidiate.vcf: ")
vcf <- read.vcfR(filename, nrows=10000)

# convert to genlight 	
aa.genlight <- vcfR2genlight.tetra(vcf)
locNames(aa.genlight) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_")  
## ----
pop(aa.genlight) <-substr(indNames(aa.genlight),1,5)  

#check   
aa.genlight
indNames(aa.genlight) 
ploidy(aa.genlight)

######################====================================######################
#   PCA     --------------------------------------------------------------------

#remove missing values:
toRemove <- is.na(glMean(aa.genlight, alleleAsUnit = FALSE)) 
which(toRemove) # position of entirely non-typed loci
aa.genlight_correct <- aa.genlight[, !toRemove]

#run PCA
pca.1 <- glPcaFast(aa.genlight_correct, nf=300)

#Plotting PCA
scatter(pca.1, posi="topleft") # plot scatter with the individual labels


# proportion of explained variance by each axes
pca.1$eig[1]/sum(pca.1$eig) # proportion of variation explained by 1st axis
pca.1$eig[2]/sum(pca.1$eig) # proportion of variation explained by 2nd axis
pca.1$eig[3]/sum(pca.1$eig) # proportion of variation explained by 3rd axis
pca.1$eig[4]/sum(pca.1$eig) # proportion of variation explained by 4th axis
pca.1$eig[5]/sum(pca.1$eig) # proportion of variation explained by 5th axis
pca.1$eig[6]/sum(pca.1$eig) # proportion of variation explained by 6th axis
pca.1$eig[7]/sum(pca.1$eig) # proportion of variation explained by 7th axis

## ----
#there are only 3 variables (genotype: 00 / 01-10 / 11) when doing PCA genomics you only focus on PCA1 and PCA2
## ----

col <- funky(41) #total locations n 
#compare pcas 1 and 2 and check graph settings 
s.class(pca.1$scores, pop(aa.genlight),  xax=1, yax=2, col=transp(col,1),
        ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=T, plabels.cex=F, ylab="PCA1", xlab="PCA2",)
        
  # save nice figs configured for my needs
pdf ("fin_no_gap_all_pop_ac30", width=14, height=7)
g1 <- s.class(pca.1$scores, pop(aa.genlight),  xax=1, yax=2, col=transp(col,.6),
              ellipseSize=0, starSize=0, ppoints.cex=2, paxes.draw=T, plegend.size=0.4, 
              ylab="PCA1", xlab="PCA2", pgrid.draw =F, plot = T, plabels.cex = F)


g2 <- s.label(pca.1$scores, xax=1, yax=2, ppoints.col = F, 
              plabels = list(box = list(draw = F), optim = TRUE),
              ylab="PCA1", xlab="PCA2", paxes.draw=T, pgrid.draw=F, plabels.cex=0.2, plot = T) 




ADEgS(c(g1, g2), layout = c(1, 2))
dev.off()

#making nj trees
NJtree <- nj(dist(as.matrix(aa.genlight_correct)))
NJtree

#plot unrooted tree
plot(NJtree, typ="unrooted", cex=0.4, hjust=-0.3)
title(expression("Unrooted Neighbour-joining tree original cohort "))
#save the tree using:
write.tree((NJtree),file="NJ.113_ac30.txt")

#plot rooted tree
plot(NJtree, cex=0.5)
title(expression("Rooted Neighbour-joining tree original cohort "))
write.tree((NJtree),file="NJ.rooted_all_ind_no_alignment2_ac30.txt")


