rm(list=ls())
library(HardyWeinberg)
library(genetics)

setwd('/Users/krishna/MIRI/BSG/R_BSG/Assignment2')
load('CHBChr2-2000.rda')
ls()
X = t(X)
X
missing_col = which(colMeans(is.na(X)) == 1)
remove_mono = which(apply(X,2,function(x) length(unique(na.omit(unlist(x))))) == 1)
remove_cols =  c(missing_col,remove_mono)
Alleles = t(Alleles)
Alleles = Alleles[-remove_cols]
X = X[,-remove_cols]
paste("snps remaining", length(Alleles))
X[,1]

k = c(rep("TT",6),"TM")
snpg <- genotype(k,sep="")
HWE.chisq(snpg)

hw <- function(x){
  snpg <- genotype(x,sep="")
  hw <- HWE.chisq(snpg)$p.value
  return(hw)
}
hw <- apply(X,2,hw)
paste("significant snps",sum(hw > 0.05))
ncol(X) * 0.050 

ncol(X) - sum(hw > 0.05)

hws <- function(x){
  snpg <- genotype(x,sep="")
  hw = HWExact(snpg,pvaluetype="selome",verbose=F)$pval
  return(hw)
}
hw <- apply(X,2,hws)

hw = HWExact(12,pvaluetype="selome",verbose=F)$pval
pvalues_exact <- apply(genotypes,2,function(x){
  HWExact("12",pvaluetype="selome",verbose=F)$pval
})


snpg <- genotype(X,sep="")
snpg
Y <- MakeCounts(X,X[,1])
Y
hw = HWExact(unlist(Y[,1:3]),verbose=F)$pval
hwe.x <- HWChisq((na.omit(Y)),cc=0,verbose=TRUE)

k = unlist(Y[,1:3])
k
unlist(k)


Y[1,]

?HWExact
X[,1]
snpg <- genotype(X,sep=",")
snpg
MakeCounts(X[,1],c("CC/CG"))
?HWExact
SNP1 <- c("GG","GG","GG","GG","GG","GG","GG","GG","GG")
SNP2 <- c("CG","GG","CC","GG","GG","CG","CG","CG","CG")
SNP3 <- c("AA","AA","AA","AG","AA","AG","AA","AA","AA")
SNP4 <- c("GG","GG","GG","GG","GG","GG","GG","GG","GG")
SNP5 <- c("CC","CC","CC","CC","CC","CC","CT","CT","CT")
Z <- cbind(SNP1,SNP2,SNP3,SNP4,SNP5)
Y <- MakeCounts(Z,c("A/G","C/G","A/G","A/G","C/T"))
print(Y)
X

Z <- cbind(X[,1])
Z2 <- cbind(X[,2])
Y <- MakeCounts(Z,c("A/G","C/G","A/G","A/G","C/T","G/G"))
Y <- MakeCounts(Z2,Alleles[2,])
Y <- MakeCounts(X,Alleles,na.rm=TRUE)
Y
geno_count
ml <- apply(snpg,2,function(x){
  HWLratio(x, verbose = F, x.linked = FALSE)$pval
})
ml <- sum(ml > 0.05)

unique_genotypes <- UniqueGenotypeCounts(t(Y),verbose=F)
unique_genotypes
unique_genotypes <- UniqueGenotypeCounts((Y),verbose=F)
unique_genotypes
Y

missing = which(colMeans(is.na(X))==1)
num_missing = length(missing)
if (num_missing>0){ X <- X[,-missing]}

# Remove monomporhic
monomorphic = which(apply(X,2,function(x) length(unique(na.omit(unlist(x))))) == 1)
num_monomorphic = length(monomorphic)
if (num_monomorphic>0){X <- X[,-monomorphic]}

#Remaining SNPs
paste("SNPs remaining =", ncol(X))
identical_letters <- function(x){
  identical(unlist(strsplit(x,split=""))[1],unlist(strsplit(x,split=""))[2])
}
#Function to get the genotypes
get_genotype <- function(x) { 
  genotype <- table(x)
  if(length(genotype)==2){
    if(identical_letters(names(genotype)[1]) & identical_letters(names(genotype)[2])){
      genotype <- c(genotype[1], 0, genotype[2])
    }else if(identical_letters(names(genotype)[1]) & !identical_letters(names(genotype)[2])){
      genotype <- c(genotype[1], genotype[2],0)
    }else{
      genotype <- c(0, genotype[1], genotype[2])
    }
    
  }else{
    genotype <- genotype
  }
  return(unname(genotype))
}

#Generate Genotypes
genotypes <- apply(X,2,get_genotype)


#Significant SNPs without continuity correction
pvalues <- apply(genotypes,2,function(x){
  HWChisq(x,cc = 0,verbose=F)$pval
})
paste("Significant SNPs =", sum(pvalues > 0.05))
hist(pvalues,20)

chisel = apply(unlist(geno_count[,1:3]), 1, function(x){
  HWExact(x,verbose=F)$pval
})
hist(chisel,20)

geno <- unlist(geno_count[,1:3])
#Sample Size
n = unname(rowSums(unlist(geno_count[,1:3])))
#Allele frequency
p = unname(apply(geno,1,maf))

#Simulated value for each genotype
simul <-c()
for(i in 1:length(n)){
  s <- HWData(1,n=n[i],p=p[i])
  simul <- rbind(simul,s)
}
length(simul)
length(unique(simul))

chisq_simul <- apply((simul),1,function(x){
  HWChisq(x,cc = 0,verbose=F)$chi
})

n = unname(unlist(rowSums(Alleles)))
p = unname(apply(Alleles,1,maf))

simul <-c()
for(i in 1:length(n)){
  s <- HWData(1,n=n[i],p=p[i])
  simul <- rbind(simul,s)
}

chisq_simul <- apply(t(simul),2,function(x){
  HWChisq(x,cc = 0,verbose=F)$chi
})
qqplot(x=chi_sq,y=chisq_simul,xlab = "Chi square Database", 
       ylab="Chi square Simulated", main="QQ plot",
       xlim=c(0,50),ylim=c(0,50))
