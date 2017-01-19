######################################################
##
##  Script 2: Exercise on Hardy-Weinberg equilibrium
##


# 1) Install HardyWeinberg package


#install.packages("HardyWeinberg")
library(HardyWeinberg)



# 2) HWE for a C/G polymorphism

x <- c(CC=23,CG=48,GG=29)
hwe.x <- HWChisq(x,cc=0,verbose=TRUE)
hwe.x.cc <- HWChisq(x,cc=0.5,verbose=TRUE)
hwe.x.exa <- HWExact(x,verbose=TRUE)


# 3) HWE for a C/T polymorphism

y <- c(CC=0,CT=7,TT=93)
hwe.y <- HWChisq(y,cc=0,verbose=TRUE)
hwe.y.cc <- HWChisq(y,cc=0.5,verbose=TRUE)
hwe.y.exa <- HWExact(y,verbose=TRUE)


# 4) Ternary plot representation

X <- rbind(x,y)
X
colnames(X) <- c("AA","AB","BB")
out <- HWTernaryPlot(X)


# 5) Function for a permutation test based on

# chi-square with or without continuity correction

# given the genotype counts, return the p-value, after N permutations

permutation_test = function(x,N,c){
  
  chi2_0 = HWChisq(y,cc=c,verbose=F)$chi
  
  nA = 2*x[1] + x[2]
  nB = 2*x[3] + x[2]
  
  
  
  chi2_perm = rep(0,N)
  
  
  for(i in 1:N){
    
    alleles = c(rep("A",nA),rep("B",nB))
    
    new_alleles = sample(alleles)
    
    new_genotypes = matrix(new_alleles,ncol=2,byrow=F)
    
    new_genotypes = do.call(paste0,list(new_genotypes[,1],new_genotypes[,2]))
    
    new_genotypes[new_genotypes == "BA"] = "AB"
    
    new_counts = c(AA = sum(new_genotypes=="AA"),
                   AB = sum(new_genotypes=="AB"),
                   BB = sum(new_genotypes=="BB"))
    
    chi2_perm[i] = HWChisq(new_counts,cc=c,verbose=F)$chi
    
  }
  return(sum(chi2_perm>chi2_0)/N)
}

permutation_test(x,N=1000,c=0)
permutation_test(y,N=1000,c=0.5)

# with HardyWeinberg package

HWPerm
hwe.x.per <- HWPerm(x,verbose=TRUE)
hwe.y.per <- HWPerm(y,verbose=TRUE)



# 7) Simulate 100 SNPs

Z <- HWData(100)
Z
out <- HWTernaryPlot(Z)


# chi square test

out_chi = HWChisq(Z[1,])

out_chi # we will save the $chisq and $pval 

out.chi <- matrix(unlist(apply(Z,1,HWChisq)),ncol=5,byrow=TRUE)
head(out.chi)
dim(out.chi)

X2 <- out.chi[,1]
pval <- out.chi[,2]

sum(is.na(X2))
Z[is.na(X2),]

sum(X2[!is.na(X2)] >= qchisq(0.95,1)) # comparing with X2 distribution
sum(pval <= 0.05) # comparing with the significance level


# exact test

exact.pval <- numeric(nrow(Z))
for(i in 1:nrow(Z)) {
  exact.pval[i] <- HWExact(Z[i,])$pval
}

sum(exact.pval <= 0.05) # more conservative



# 8) Histogram

hist(X2)

Z <- HWData(10000)
out.chi <- matrix(unlist(apply(Z,1,HWChisq)),ncol=5,byrow=TRUE)
X2 <- out.chi[,1]
hist(X2)

boxplot(X2)
hist(X2[X2<10],breaks=30)


hist(rchisq(10000,1)) # random chisquare sample of N = 10000

# It follows a Chi square distribution of one degree of freedom
