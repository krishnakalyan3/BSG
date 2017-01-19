rm(list=ls())
library(HardyWeinberg)
library(genetics)
library(LDheatmap)
setwd('/Users/krishna/MIRI/BSG/R_BSG/Assignment3')

# 1 Load
par(mfrow=c(1,1))
# Load data
load('ABO-CHB.rda')

ls()
# 2 Missing data
head(Z)
ind = ncol(Z)
snps = nrow(Z)
Z = t(Z)
missing = apply(Z,2,function(x)(is.na(x)))
per_missing = sum(missing)/length(missing)
per_missing * 100


# 3. ternary plot
snpg <- genotype(Z,sep="")
geno_count = MakeCounts(Z,alleles)
ug <- UniqueGenotypeCounts(unlist(geno_count[,1:3]),verbose=F)
ternary <- HWTernaryPlot(ug[,c(1:3)])
# Interpretation 
# We see that most of the observation are in equilibrium
# there are 2 observations that are significant
# 5% of the observations might lost by chance
# 45 * 0.05
# Yes I belive Hardy-Weinberg equilibrium is tenable

# 4 LD 
ls()
SNP1 = genotype(Z[,1],sep="")
SNP2 = genotype(Z[,2],sep="")
snp12 = data.frame(SNP1,SNP2)
LD(snp12)
# Ho : No Assoiciation
# Low Cor
# Aceept Ho 
# D` is close to 1 hence its assoiciated

# 5
# By Hand comput vals


# 6
# which haplotye is commom

# 7
geno <- function(x){
  return(genotype(x, sep=""))
}
gt = apply(Z,2,geno)
test = LD(makeGenotypes(gt))
vecD <- test$D[upper.tri(test$D)]
vecD_p <- test$`D'`[upper.tri(test$`D'`)]
vec_chi <- test$`X^2`[upper.tri(test$`X^2`)]
vec_rsq <- test$`R^2`[upper.tri(test$`R^2`)]
vectors <- as.data.frame(cbind("D" = vecD, "D'" = vecD_p, "X^2" = vec_chi, "R^2" = vec_rsq))
plot(vectors)

# 8
plot(vectors$`R^2`, matrix(dist(pos)))

# 9
LDheatmap(test$`R^2`)
LDheatmap(test$`D'`)

# 10

