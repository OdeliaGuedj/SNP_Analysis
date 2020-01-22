### -------------------------------------------------------------------
### TP1
### -------------------------------------------------------------------
setwd("C:/Users/odeli/Desktop/statisticsForGeneticsAndGenomics_M2/tp1")

### 3. Loading                                  
### genotypes
vgeno<-read.table("data_matrix.txt",comment.char="!",sep="\t",header=TRUE,row.names=1,na.strings = "NoCall",nrows=1000)

### phenotypes
vpheno <- as.factor(t(read.table("data_matrix.txt", skip = 36, nrows = 1 ))[-1])


### 4. Sample sizes
table(vpheno)


### 5. Call Rates
### If "NoCall==NA"
CallRateCalc <- function(v){
  res <- 1-sum(is.na(v))/length(v)
  # equivalently:  res <- sum(!is.na(v))/length(v)
  return(res)
}

SNPCallRate <- apply(vgeno,1,CallRateCalc) # Call Rates calculation for SNPs
SampleCallRate <- apply(vgeno,2,CallRateCalc) # Call Rates calcultation for samples/individuals

par(mfrow=c(1,2))
hist(SNPCallRate); hist(SampleCallRate)

sum(SNPCallRate<0.95) # Number of filtered SNPs at a 0.95 level
vgeno2 <- vgeno[SNPCallRate>=0.95,] # Filtered genotypes matrix 


### 6. 
# a) MAF calculation
maf.calc <- function(vsnp){
  vsnp <- na.omit(unlist(vsnp))
  fA <- ( 2*sum(vsnp=="AA")+sum(vsnp=="AB") ) /(2*length(vsnp))
  res <- min(fA,1-fA)
  return(res)
}

vi <- sample(1:1000,1)
maf.calc(vgeno2[vi,])

resmaf <- apply(vgeno2,1,maf.calc)
hist(resmaf)
sum(resmaf<0.05,na.rm=TRUE)
sum(is.na(resmaf))
vgeno3 <- vgeno2[ !is.na(resmaf)&resmaf>=0.05,] # Filtered genotypes matrix 

# b) Hardy Weinberg equilibrium test
HWtest <- function(vsnp){
  vsnp <- na.omit(unlist(vsnp))
  restable <- table(vsnp)
  vn <- length(vsnp)
  hatp <- ( 2*sum(vsnp=="AA")+sum(vsnp=="AB") ) /(2*length(vsnp))
  Obs <- c(ifelse("AA"%in%names(restable),restable["AA"],0),
           ifelse("AB"%in%names(restable),restable["AB"],0),
           ifelse("BB"%in%names(restable),restable["BB"],0))
  Exp <- length(vsnp)*c(hatp^2,2*hatp*(1-hatp),(1-hatp)^2)
  vstat <- sum((Obs-Exp)^2/Exp)
  vpval <- 1-pchisq(vstat,1)
  return(vpval)
}

HWtest(vgeno3[2,])
levels(vpheno) <- c("cases","controls")
HWtest(vgeno3[2,vpheno=="controls"])
resHW <- apply(vgeno3[,vpheno=="controls"],1,HWtest)

sum(resHW<0.001,na.rm=TRUE) ## Number of snps that would be removed at 0.001 level
vgeno3 <- vgeno3[!is.na(resHW)&resHW>=0.001,] # Filtered genotypes matrix 

# c) $\chi^2$ test for association for SNP i
### Contingency tables "tageno" (genotypes) and "taballel" (alles) for SNP i=7
i <- 7 
restable <- tapply(unlist(vgeno3[i,]),vpheno,table)
for(k in 1:2){
  for(j in 1:3){
    restable[[k]][c("AA","AB","BB")[j]] <- ifelse(c("AA","AB","BB")[j]%in%names(restable[[k]]),restable[[k]][c("AA","AB","BB")[j]],0)
  }
  restable[[k]] <- restable[[k]][order(names(restable[[k]]))]
}

tabgeno <- matrix(unlist(restable),ncol=2)
taballel <- matrix(c(2*restable$cases["AA"]+restable$cases["AB"],2*restable$cases["BB"]+restable$cases["AB"],2*restable$controls["AA"]+restable$controls["AB"],2*restable$controls["BB"]+restable$controls["AB"]),2,2)
print(t(tabgeno))
print(taballel)

### allelic
resa <- chisq.test(taballel)
resa$expected
resa
### genotypic
## resg <- chisq.test(matrix(unlist(tapply(unlist(vgeno3[i,]),vpheno,table)),ncol=2))
resg <- chisq.test(tabgeno)
resg$expected
resg

# d) For SNP "SNP_A-1780342"
fisher.test(taballel)$p.value
#chisq.test(taballel)$p.value
#chisq.test(tabgeno)$p.value

# e) Logistic regression
vX <- as.factor(unlist(vgeno3[i,]))
resglm <- glm(vpheno~vX,family=binomial(link="logit"))
summary(resglm)

# OR estimate :
exp(resglm$coefficients[2])
exp(resglm$coefficients[3])

# f) Logistic regression
vX <- relevel(vX,ref="BB")
resglm2 <- glm(vpheno~vX,family=binomial(link="logit"))
summary(resglm2)
exp(resglm2$coefficients[2])
exp(resglm2$coefficients[3])

# alternatively:
vX <- as.numeric(vX)
resglm2 <- glm(vpheno~vX,family=binomial(link="logit"))
summary(resglm2)

# g) Cochran Armitage
library(DescTools)
CochranArmitageTest(tabgeno)

