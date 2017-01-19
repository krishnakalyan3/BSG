######################################################
##
##  Script session 5: Population substructure
##


# 1) Load data

library(MASS)
library(HardyWeinberg)

load(url("http://www-eio.upc.es/~jan/data/bsg/CHBChr2-200.rda"))

ls()

# 2) Convert the genotype data into an n x n distance matrix.

dim(X)
X <- t(X)

X[1:5,1:5]

all.missing <- function(x) {
  n <- length(x)
  n.mis <- sum(is.na(x))
  y <- n==n.mis
  return(y)
}

all.missing(X[,1])

index.all.missing <- apply(X,2,all.missing)
sum(index.all.missing)

X <- X[,!index.all.missing]
Alleles <- Alleles[!index.all.missing]

dim(X)
length(Alleles)

sum(is.na(X))

X[1:5,1:5]

?recode
table(X[,1],useNA="always")
y <- recode(X[,1:2],Alleles[1:2])
Y <- recode(X,Alleles)


library(SNPassoc)
additive(X[,1])


Y = apply(X,2,additive)


De <- as.matrix(dist(Y))
sum(is.na(De))


# 3) Produce a map of the individuals by metric multidimensional scaling 

out.mds <- cmdscale(De,eig=TRUE)
attributes(out.mds)
dim(out.mds$points)
X.sol <- out.mds$points[,1:2]
plot(X.sol[,1],X.sol[,2],asp=1,xlab="First principal axis",
     ylab="Second principal axis")
abline(v=0,lty=3)
abline(h=0,lty=3)
text(X.sol[,1],X.sol[,2],1:nrow(X),pos=1)
out.mds$GOF
lamb <- out.mds$eig
lamb
sum(lamb[1:2])/sum(abs(lamb))
sum(lamb[1:2])/sum(lamb[lamb>0])



# 4) Make a graph of the fitted against the observed distances,

D.fitted <- as.matrix(dist(X.sol))

D.obs <- De[lower.tri(De)]
D.fit <- D.fitted[lower.tri(D.fitted)]

plot(D.obs,D.fit,asp=1,xlab="Observed distance",ylab="Fitted distance")
abline(0,1,col="red",lw=2)
out.lm <- lm(D.fit~D.obs)
summary(out.lm)$coef
abline(summary(out.lm)$coef[1,1],summary(out.lm)$coef[2,1],col="green",lw=2)



# 5) Produce a map of the individuals by non-metric multidimensional scaling.


?isoMDS
out.nmds.1 <- isoMDS(De,k=2)
attributes(out.nmds.1)
set.seed(1234)
x.init <- scale(matrix(runif(2*nrow(X)),ncol=2))
out.nmds.2 <- isoMDS(De,y=x.init)
X.nmds <- out.nmds.1$points
plot(X.nmds[,1],X.nmds[,2],asp=1)
text(X.nmds[,1],X.nmds[,2],1:nrow(X),pos=1)

dim(X.sol)
dim(X.nmds)
pairs(cbind(X.sol,X.nmds))
round(cor(cbind(X.sol,X.nmds)),digits=2)

cor(cbind(X.sol,X.nmds))




