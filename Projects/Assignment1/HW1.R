rm(list=ls())
setwd('/Users/krishna/MIRI/BSG/R_BSG')

# 1 Load Chinese Dataset
load('CHBChr2.rda')
head(X)

# 2 Replace NN with NA
X[X == 'NN'] = NA
mean(is.na(X)) * 100
# 36.428

# Transposing X
X = t(X)

# 3 SNPs the genotype information is completely missing?
length(colMeans(is.na(X)))
missing_snp = sum(colMeans(is.na(X)) == 1)
# 2975

X.clean = X[,colMeans(is.na(X)) != 1]

# 4 What is, on the average, the percentage of missing information per individual?
mean(unname(unlist(rowMeans(is.na(X.clean))))) * 100
#  9.50 %
head(X.clean)

# 5 How many markers are monomorphic?
X.table = apply(X.clean, 2, table)
sum(lapply(X.table,length) == 1)
# 1860

# 6 Write a function to compute the minor allele frequency?
maf <- function(x){
  if(mean(is.na(x))==1){
    return(NA)
  }
  alleles <- table(unlist(strsplit(as.character(unlist(unname(x))),split="")))
  alleles_freq <- alleles/sum(alleles)
  m <- min(alleles_freq)
  if(m == 1){
    return(0)
  }else{
    return(m)
  }
}

maf = apply(X.clean , 2, maf) 


# 7
# Compute the minor allele frequencies for all markers, and make a histogram of it.
hist(maf)

# 8
# What percentage of the markers have a maf below 0.05? And below 0.01?
mean(maf < 0.05) * 100
# 34.51957
mean(maf < 0.01) * 100
# 26.03

# 9
# Make a histogram of the expected heterozygosity.
heterozygosity <- function(x){
  alleles <- table(unlist(strsplit(as.character(unlist(unname(x))),split="")))   
  af <-  alleles/sum(alleles)
  het <- 1- sum(af^2)
  return(het)
}

heter <- apply(X.clean,2,heterozygosity)
mean(heter)
hist(heter)



