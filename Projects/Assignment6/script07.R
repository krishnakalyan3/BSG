#install.packages("HardyWeinberg")
library(HardyWeinberg)
load(url("http://www-eio.upc.es/~jan/data/bsg/Yoruba10000.rda"))
ls()
head(X.Fam)

X.Geno[1:5,1:5]

nrow(X.Geno)
ncol(X.Geno)

n.alleles <- function(x) {
  y <- length(unique(unlist(strsplit(names(table(x)),""))))
  return(y)
}

n.alleles(X.Geno[,1])

numalleles <- apply(X.Geno,2,n.alleles)
table(numalleles)

get.alleles <- function(x) {
  y <- unique(unlist(strsplit(names(table(x)),"")))
  if(length(y)==1) als <- paste(y,"X",sep="/") else
    als <- paste(y[1],y[2],sep="/")
  return(als)
}

get.alleles(X.Geno[,2])

alleles <- apply(X.Geno,2,get.alleles)
head(alleles)

X <- recode(X.Geno,alleles)

data.frame(X.Geno[,1:3],X[,1:3])

ibs.mean <- function(x,y) {
  y <- mean(2 - abs(x - y),na.rm=TRUE)
  return(y)
}

ibs.sd <- function(x,y) {
  y <- sd(abs(x-y),na.rm=TRUE)
  return(y)
}

ibs.mean(X[1,],X[2,])
ibs.sd(X[1,],X[2,])
n <- nrow(X)
n

Dmean <- matrix(NA,nrow=n,ncol=n)
Dsd <- matrix(NA,nrow=n,ncol=n)

for(i in 1:n) {
  for(j in 1:n) {
    Dmean[i,j] <- ibs.mean(X[i,],X[j,])
    Dsd[i,j] <- ibs.sd(X[i,],X[j,])
  }
}
Dmean[1:5,1:5]
Dsd[1:5,1:5]

ibs.m <- Dmean[lower.tri(Dmean)]
ibs.s <- Dsd[lower.tri(Dsd)]

plot(ibs.m,ibs.s,xlim=c(1,2),ylim=c(0,1),xlab="Mean",ylab="Standard deviation")

head(X.Fam)

Relationship <- matrix(NA,nrow=n,ncol=n)

is.po <- function(x,y) {
  xchildy  <- x[1]==y[1] & (x[3]==y[2] | x[4]==y[2]) # x child of y
  xparenty <- x[1]==y[1] & (x[2]==y[3] | x[2]==y[4]) # x parent of y
  y <- xchildy | xparenty
  return(y)
}

is.po(X.Fam[1,],X.Fam[2,])
is.po(X.Fam[1,],X.Fam[3,])
is.po(X.Fam[2,],X.Fam[3,])
is.po(X.Fam[2,],X.Fam[1,])
is.po(X.Fam[3,],X.Fam[1,])
is.po(X.Fam[3,],X.Fam[2,])

for(i in 1:n) {
  for(j in 1:n) {
    Relationship[i,j] <- is.po(X.Fam[i,],X.Fam[j,])
  }
}

Relationship[1:5,1:5]

rel.pair <- Relationship[lower.tri(Relationship)]

colvec <- rep("blue",length(rel.pair))
colvec[rel.pair] <- "red"

plot(ibs.m,ibs.s,xlim=c(1,2),ylim=c(0,1),xlab="Mean",ylab="Standard deviation",
     col=colvec)







