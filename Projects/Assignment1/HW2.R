# STR Dataset
rm(list=ls())
setwd('/Users/krishna/MIRI/BSG/R_BSG')

# 1
# Load Data
data = read.table("FrenchSTRs.dat")
dim(data)

# 2
# How many individuals and how many STRs contains the database?
X <- data[,2:ncol(data)]
n <- nrow(X)/2
p <- ncol(X)
paste("Individuals", n)
paste("STRs ", p)

# 3
# The value 􀀀9 indicates a missing value. Replace all missing values by NA. What percentage
# of the total amount of datavalues is missing?
sum(X==-9)
X[X==-9] <- NA
mean(is.na(X)) * 100
# 4.20

# 4
# (mean, standard deviation, median, minimum, maximum) for each STR in the database.
n.alleles <- function(x) {
  y <- length(unique(x[!is.na(x)]))
  return(y)
}
alleles.str = as.data.frame(apply(X, 2, n.alleles))
names(alleles.str) = c('alleles')
# mean
lapply(alleles.str, mean)
# 6.374631

# sd
lapply(alleles.str, sd)
# 1.823385

# medain
lapply(alleles.str, median)
# 6

# max
lapply(alleles.str, max)
# 16

# min
lapply(alleles.str, min)
# 3

# 5
# Make a boxplot and a histogram of the number of alleles per STR. What is the most common
# number of alleles for an STR?
par(mfrow=c(1,2))
hist(alleles.str$alleles, ylab="Frequency",xlab="Number of Alleles", main="Alleles Histogram")
boxplot(alleles.str$alleles,ylab="Number of Alleles", main="Alleles Boxplot")

# 6
# Compute the expected heterozygosity for each STR. Make a histogram of the expected 
# heterozygosity over all STRS. Compute the average expected heterozygosity over all STRs.
heterozygosity <- function(x){
  freq <- table(x)/sum(table(x))
  1 - sum(freq^2)
}

allels.het = apply(X, 2, heterozygosity)
m
ean(allels.het)
# 0.7172662

hist(allels.het)
#
 7
# Compare the results you obtained for the SNP database with those you obtained for the STR
# database. What differences do you observe between these two types of genetic markers?
# The expected heterozygosity of the STRs is more than twice 
# the expected heterozygosity of the SNPs

