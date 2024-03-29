if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager",repos = "http://cran.us.r-project.org")
BiocManager::install("gdsfmt")
BiocManager::install("GWASTools")
library(gdsfmt)
# to update Matrix package:
# https://cran.r-project.org/bin/macosx/tools/
# sudo cp -r /opt/gfortran /usr/local
# vim ~/.R/Makevars
# FLIBS=-L/usr/local/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin20.0/12.2.0
library(GWASTools)
library(SNPRelate)
library(GENESIS)
library(SeqVarTools)

directory <- "/Users/wletsou/Library/CloudStorage/OneDrive-NewYorkInstituteofTechnology/Courses/BIOL 350 Spring 2024/bioGWAS/data"

# convert continuous phenotype to binary
pheno <- read.table(sprintf("%s/pat_snps_sim_phenos.tsv",directory),header = FALSE) # computed phenotypes
colnames(pheno) <- c("FID","IID","Phenotype")
if (!file.exists(sprintf("%s/pat_filt_sim_continuous.fam",directory))) {
  fam <- read.table(sprintf("%s/pat_filt_sim.fam",directory),header = FALSE) # plink fam file
  write.table(fam,sprintf("%s/pat_filt_sim_continuous.fam",directory),col.names = FALSE,row.names = FALSE,quote = FALSE)
} else {
  fam <- read.table(sprintf("%s/pat_filt_sim_continuous.fam",directory),header = FALSE)
}
colnames(fam) <- c("FID","IID","MID","PID","Sex","Phenotype")
set.seed(566)
fam$Phenotype <- as.numeric(runif(length(pheno$Phenotype)) <= exp(pheno$Phenotype) / (1 + exp(pheno$Phenotype))) + 1 # convert betas to disease states
write.table(fam,sprintf("%s/pat_filt_sim.fam",directory),col.names = FALSE,row.names = FALSE,quote = FALSE)
write.table(fam,"/Users/wletsou/Library/CloudStorage/OneDrive-NewYorkInstituteofTechnology/Courses/BIOL 350 Spring 2024/EUR_Bca.fam",col.names = FALSE,row.names = FALSE,quote = FALSE) # phenotypes fam file

vcf.fn <- sprintf("%s/pat_filt_sim.vcf",directory)
gds.fn <- sprintf("%s/pat_filt_sim.gds",directory)
snpgdsVCF2GDS(vcf.fn,gds.fn) # create gds file from vcf
genofile <- snpgdsOpen(gds.fn) # import gds object

pca <- snpgdsPCA(genofile) # principal components analysis

# LD pruning https://uw-gac.github.io/SISG_2021/ancestry-and-relatedness-inference.html#ld-pruning
set.seed(566)
snpset <- snpgdsLDpruning(genofile,method = "corr", slide.max.bp = 10e6,ld.threshold = sqrt(0.1), verbose = FALSE) 
pruned <- unlist(snpset, use.names = FALSE)

# IBD and kinship calculation (1)
ibd <- snpgdsIBDKING(genofile,snp.id = pruned) # KING kinship estimation
colnames(ibd$kinship) <- ibd$sample.id
rownames(ibd$kinship) <- ibd$sample.id

par(mar = c(5.1,5.1,4.1,2.1) ) # left default plus one
plot(ibd$IBS0,ibd$kinship,ylab = "Kinship coeffecient",xlab = "IBS0",main = "KING relatedness estimation", ylim = c(0,1),pch = 19,cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5) # plot kinship coefficient vs. proportion of SNPs not shared IBS

closefn.gds(genofile)

geno <- GdsGenotypeReader(sprintf("%s/pat_filt_sim.gds",directory))
genoData <- GenotypeData(geno)

vcfWrite(genoData,vcf.file = "/Users/wletsou/Library/CloudStorage/OneDrive-NewYorkInstituteofTechnology/Courses/BIOL 350 Spring 2024/EUR_Bca.vcf",sample.col = "sample.id",id.col = "snp.id")

mypcair <- pcair(genoData,kinobj = ibd$kinship,divobj = ibd$kinship,snp.include = pruned) # genotype principal components based on a subset of unrelated individuals

genoData.iterator <- GenotypeBlockIterator(genoData,snpInclude = pruned) 
mypcrel <- pcrelate(genoData.iterator,pcs = mypcair$vectors[,1:2],training.set = mypcair$unrels) # kinship based on unrelated individuals
par(mar = c(5.1,5.1,4.1,2.1) ) # left default plus one
plot(mypcrel$kinBtwn$k0,mypcrel$kinBtwn$kin,xlab = "IBD0 proportion",ylab = "Kinship coefficient",main = "PC-relate relatedness estimation",pch = 19,cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5) # kinship vs. IBD0

myGRM <- pcrelateToMatrix(mypcrel,scaleKin = 2) # genomic relationship matrix, multiplying each phi by 2

mydat <- data.frame(scanID = mypcair$sample.id,pc1 = mypcair$vectors[,1],pc2 = mypcair$vectors[,2],pc3 = mypcair$vectors[,3],pc4 = mypcair$vectors[,4],pc5 = mypcair$vectors[,5],pc6 = mypcair$vectors[,6],pc7 = mypcair$vectors[,7],pc8 = mypcair$vectors[,8],pc9 = mypcair$vectors[,9],pc10 = mypcair$vectors[,10],pheno = fam$Phenotype - 1) # data frame of fixed effects

scanAnnot <- ScanAnnotationDataFrame(mydat)
nullmod <- fitNullModel(scanAnnot, outcome = "pheno", covars = c("pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10"), cov.mat = myGRM, family = "binomial") # null model
assoc <- assocTestSingle(genoData.iterator,null.model = nullmod) # model including SNPs

par(mar = c(5.1,5.1,4.1,2.1) ) # left default plus one
plot(assoc$pos,-log10(assoc$Score.pval),xlab = "Position",ylab = "-log10(p)",pch = 19,cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5) # Manhattan plot
abline(h = 8 - log10(5),col = 'red',lty = 2) # genome-wide significance threshold
123331690
close(genoData)

# compare to continuous phenotype

pat_snps_sim_gwas <- read.delim("~/Library/CloudStorage/OneDrive-NewYorkInstituteofTechnology/Courses/BIOL 350 Spring 2024/bioGWAS/data/pat_snps_sim_gwas.tsv")
par(mar = c(5.1,5.1,4.1,2.1) ) # left default plus one
plot(pat_snps_sim_gwas$pos,-log10(pat_snps_sim_gwas$pval),xlab = "Position",ylab = "-log10(p)",pch = 19,cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5) # Manhattan plot
abline(h = 8 - log10(5),col = 'red',lty = 2) # genome-wide significance threshold
