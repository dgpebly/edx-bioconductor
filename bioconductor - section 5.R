BiocManager::install(c("SpikeInSubset",
                       "genefilter",
                       "qvalue",
                       "limma",
                       "genomicsclass/maPooling",
                       "leukemiasEset",
                       "preprocessCore"))
## biological versus technical variability:
library(devtools)
install_github("genomicsclass/maPooling")
library(Biobase)
library(maPooling)
data("maPooling")
head(maPooling)
head(pData(maPooling))
# illustrating which mice included in which samples
library(rafalib)
mypar()
flipt <- function(m) t(m[nrow(m):1,])
myimage <- function(m,...) {
  image(flipt(m), xaxt="n", yaxt="n",...)
}
myimage(as.matrix(pData(maPooling)), col=c("white", "black"),
        xlab="experiments",
        ylab="individuals",
        main="phenoData")
# we want to detect genes differentially expressed between two mice strains:
# we can apply tests to pooled samples. 
# can identify pooled samples bc all mice represented in these samples and sum of rows of experimental design matrix adds to 12
data("maPooling")
pd = pData(maPooling)
# determine strains from column names:
factor(as.numeric(grepl("b", names(pooled))))
pooled = which(rowSums(pd)==12)
# compare mean expression between groups:
i = 11425; j = 11878
pooled_y = exprs(maPooling[,pooled])
pooled_g = factor(as.numeric(grepl("b", names(pooled))))
mypar(1,2)
stripchart(split(pooled_y[i,], pooled_g), vertical=TRUE, method="jitter",
           col=c(1,2), main="Gene 1", xlab ="Group", pch=15)
stripchart(split(pooled_y[j,], pooled_g), vertical=TRUE, method="jitter",
           col=c(1,2), main="Gene 2", xlab="Group", pch=15)
# compute a t-test from these values:
library(genefilter)
pooled_tt = rowttests(pooled_y, pooled_g)
pooled_tt$p.value[i]
pooled_tt$p.value[j]
# we are replicating experimental protocol here
# created four technical replicates for each pooled sample
# for each strain we have 12 biological replicates:
individuals=which(rowSums(pd)==1)
# some technical replicates included so we remove them:
individuals = individuals[-grep("tr", names(individuals))]
y = exprs(maPooling)[,individuals]
g = factor(as.numeric(grepl("b", names(individuals))))
# compute sample variance for each gene and compare to standard deviation from technical replicates:
technicalsd <- rowSds(pooled_y[,pooled_g==0])
biologicalsd <- rowSds(y[,g==0])
LIM = range(c(technicalsd, biologicalsd))
mypar(1,1)
boxplot(technicalsd, biologicalsd, names=c("technical", "biological"),
        ylab="standard deviation")
# biological variance is larger than technical variance
# variability of variances is larger for biological variance
mypar(1,2)
stripchart(split(y[i,], g), vertical=TRUE, method="jitter", col=c(1,2),
           xlab="Gene 1", pch=15)
points(c(1,2), tapply(y[i,], g, mean), pch=4, cex=1.5)
stripchart(split(y[j,], g), vertical=TRUE, method="jitter", col=c(1,2),
           xlab="Gene 2", pch=15)
points(c(1,2), tapply(y[j,], g, mean), pch=4, cex=1.5)
# now the p values tell a different story:
library(genefilter)
tt = rowttests(y,g)
tt$p.value[i]
tt$p.value[j]
