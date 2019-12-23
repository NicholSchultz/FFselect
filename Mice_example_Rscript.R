library(lrgpr)
as.character.formula<-formula.tools:::as.character.formula
library(plyr)


###Mouse example
###load mouse data from BGLR packages
library(BGLR)
data(mice, package="BGLR")
colnames(mice.X)<-gsub("_A","", colnames(mice.X))
colnames(mice.X)<-gsub("_G","", colnames(mice.X))
colnames(mice.X)<-gsub("_C","", colnames(mice.X))
colnames(mice.X)<-gsub("_T","", colnames(mice.X))
rm(mice.A)

###load map file for mouse data
mouse.map<-read.table("http://mtweb.cs.ucl.ac.uk/HSMICE/GENOTYPES/mapfile.txt", header=TRUE)
#save(mouse.map, file="mouse.map.Rdat")
#load("mouse.map.Rdat")
table(colnames(mice.X)%in%mouse.map[,'marker'])
mouse.map<-mouse.map[match(colnames(mice.X), mouse.map[,'marker']),]
mouse.map[,'SNP.name']<-paste('Chr', mouse.map[,2], mouse.map[,3], sep="_")
x<-!duplicated(mouse.map[,'SNP.name'])
mouse.map<-mouse.map[x,]
mice.X<-mice.X[,x]

###rename genotype matrix columns in the format "Chr_1_3187380"
colnames(mice.X)<-mouse.map[,'SNP.name']

###match(phenotype matrix to genotype matrix order)
table(rownames(mice.X)%in%mice.pheno[,'SUBJECT.NAME'])
length(unique(mice.pheno[,'SUBJECT.NAME']))
mice.pheno<-mice.pheno[match(rownames(mice.X) ,mice.pheno[,'SUBJECT.NAME']),]

###remove X chromosome markers
mice.X<-mice.X[,!seq(1:ncol(mice.X))%in%grep("Chr_X_", colnames(mice.X))]

###remove duplicate genotypes
set.seed(1234)
mice.X<-mice.X[ ,sample(1:ncol(mice.X), ncol(mice.X))]
x<-which(duplicated(t(mice.X[,1:ncol(mice.X)])))
duplicated.snps<-colnames(mice.X)[x]
mice.X<-mice.X[,!colnames(mice.X)%in%duplicated.snps]

###order genotype matrix by genome position
chr<-as.numeric(ldply(strsplit(colnames(mice.X),split= "_"))[[2]])
bp<-as.numeric(ldply(strsplit(colnames(mice.X),split= "_"))[[3]])
M<-mice.X[ ,order(chr, bp)]
rm(chr, bp, mice.X)

###recode to minor allele
x<-which((colMeans(M[,1:ncol(M)])*0.5)>0.5)
M[,x]<-(M[,x]*-1)+2
#hist(colMeans(M[,1:ncol(M)])*0.5)
#range(colMeans(M[,1:ncol(M)])*0.5)
M<-M[,(colMeans(M[,1:ncol(M)])*0.5)>=0.05]

###Create binary file of genotype matrix to conserve memory in R for large genotype matrices
z <- filebacked.big.matrix(nrow(M), ncol(M),backingfile="M.bin", descriptorfile="M.desc", dimnames=list(rownames(M), colnames(M)))
z[1:nrow(M),1:ncol(M)]<-  M[1:nrow(M),1:ncol(M)]
rm(z, M)
invisible(gc())
M <- attach.big.matrix("M.desc")

mouse.map<-mouse.map[match(colnames(M), mouse.map[,'SNP.name']),]
rownames(mouse.map)<-mouse.map[,'SNP.name']

###create shared environment variable
shared.env<-as.factor(as.numeric(mice.pheno[,'cage']))

###Simulate trait with 50% heritability (30 QTL) and 10% environment and 40% residual
set.seed(1234)
n_loci<-30
causal_loci<-sample(colnames(M), n_loci)
loci_effects<-rnorm(n_loci,0, 1)
names(loci_effects)<-causal_loci
additive<-M[,match(causal_loci, colnames(M))]%*%(loci_effects)
rnd_int<-rnorm(length(unique(shared.env)), 0, sqrt(var(additive)*0.2))
names(rnd_int)<-unique(shared.env)
env_effect<-rnd_int[match(shared.env, names(rnd_int))]
residuals<-rnorm(nrow(M), 0, sqrt(var(additive)*0.8))
y<-additive + env_effect + residuals	
y<-as.vector(scale(y))

###Simulate covariate
simulated.covariate<-rnorm(length(y))

### y=numeric phenotype vector
### M=n x p genotype matrix (genotypes coded 0,1,2, for copies of minor allele and columns ordered according to genome position and colnames in the format "Chr_1_3187380", "Chr_1_3397373", etc)
### covariate.mat=numeric covariate matrix (first column must be column of one's for the intercept)
### env=factor variable for shared environments
### allSnp.weight=0.05; select features until genetic variance estimate is reduced by 95%
### max.features=50 maximum number of loci selected for the featured SNP GRM
### adjacent.snps= number of SNPs upstream and downstream of selected SNP at each feature locus to be included in the featured SNP GRM
### LD.window=base pair distance from test SNP for removing select SNPs from feature GRM
### ncores=number of cores for parallel processing

#load FFselect functions
source("FFselect_functions.R")

#Example with shared environment effect, no covariates
Result<-FFselect(y=y, M=M, covariate.mat=NULL, env=shared.env, allSnp.weight=0.05, max.features=50, adjacent.snps=3, LD.window=500000, ncores=1)

#Example with shared environment effect and a single covariate
#Result<-FFselect(y=y, M=M, covariate.mat=cbind(1, simulated.covariate), env=shared.env, allSnp.weight=0.05, max.features=50, adjacent.snps=3, LD.window=500000, ncores=1)

#Example with no shared environment effect and no ccovariates
#Result<-FFselect(y=y, M=M, covariate.mat=NULL, env=NULL, allSnp.weight=0.05, max.features=50, adjacent.snps=3, LD.window=500000, ncores=1)

names(Result)
print(Result$featureSNP)
print(c(Result$RM1weight,Result$RM2weight,Result$RM3weight))

##Plot Results		
result<-cbind(colnames(M), ldply(strsplit(colnames(M),split= "_"))[[2]], ldply(strsplit(colnames(M),split= "_"))[[3]], Result$pvals)
colnames(result)<-c('SNP', 'Chr', 'Pos', 'P.value')
result<-as.data.frame(result)
result[,3]<-as.numeric(as.character(result[,3]))
result[,4]<-as.numeric(as.character(result[,4]))
head(result)

bonferonni.correction<--log10(0.05/ncol(M))
plot(-log10(result$P.value) ~ seq(1:length(result$P.value)), col=ifelse(as.numeric(gsub("X", "20", result$Chr)) %% 2 == 0,"gray17","turquoise4"), xaxt="n", main='', xlab='Genome position', ylab=expression("-Log"[10]*'(p-value)'), cex=0.6, )
axis(1, at=sapply(unique(result$Chr), function (x) which(result$Chr==x & result$Pos==result$Pos[result$Chr==x][trunc(length(result$Pos[result$Chr==x])/2)])), labels=seq(1:length(unique(result$Chr))))
abline(h=bonferonni.correction, v=which(rownames(result)%in%causal_loci))


