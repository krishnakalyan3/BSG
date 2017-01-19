
library(HardyWeinberg)
library(genetics)
library(LDheatmap)


# 2) Load data

load(url("http://www-eio.upc.es/~jan/data/bsg/CHBChr2-2000.rda"))


# 3) Calculate D, D', R2 and Chi2 for SNP12 and SNP13

X[1:10,1:10]
X <- t(X)
dim(X)

X[,12]

Alleles[12]
snp12 <- genotype(X[,12],alleles=c("C","G"),sep="")

X[,13]
Alleles[13]
snp13 <- genotype(X[,13],alleles=c("C","T"),sep="")

summary(snp12)
summary(snp13)

Z <- data.frame(snp12,snp13) # It does not work for monomorphic
out <- LD(Z)


# 4) SNP12 and SNP1000

X[,1000]
Alleles[1000]
snp1000 <- genotype(X[,1000],alleles=c("A","G"),sep="")

Z <- data.frame(snp12,snp1000)
out <- LD(Z)
print(out)


# 5) Select the first 100 SNPs no missings and no monomorphics
nmis <- function(x) {
  nmis <- sum(is.na(x))
  return(nmis)
}

mono <- function(x) {
  l <- length(table(x))
  if(l==1) y <- TRUE else y <- FALSE
}    

n.missings <- apply(X,2,nmis)
n.missings

monomorphic <- apply(X,2,mono)

X2 <- X[,n.missings==0 & !monomorphic]
dim(X2)

Snp1 <- genotype(X2[,1],sep="")
Z <- data.frame(Snp1)
for(i in 2:100) {
  snp <- genotype(X2[,i],sep="")
  Z <- cbind(Z,snp)
}
Z[1:5,1:5]


# 6) Compute 4 matrices of D, D', R2 and Chi2 

out <- LD(Z)
attributes(out)
D <- out$D
Dp <- out$"D'"
R2 <- out$"R^2"
X2 <- out$"X^2"

D[1:10,1:10]


# 7) Extract the subdiagonal part

Dv <- D[upper.tri(D)]
Dpv <- Dp[upper.tri(Dp)]
R2v <- R2[upper.tri(R2)]
X2v <- X2[upper.tri(X2)]


# 8) Scatterplot of D, D', R2 and Chi2 

windows()
pairs(cbind(Dv,Dpv,R2v,X2v))


# 9) LDheatmap of D, D', R2 and Chi2 

LDheatmap(Z[,1:100],LDmeasure="D'")

windows()
rgb.palette <- colorRampPalette(rev(c("blue", "orange", "red")), space = "rgb")

LDheatmap(Z, LDmeasure="r",color=rgb.palette(18))

LDheatmap(Z, LDmeasure="D'", color=rgb.palette(18) )


