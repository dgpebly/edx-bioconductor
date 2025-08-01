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


## t tests in genomics
# load pooling experiment data:
library(Biobase)
library(maPooling)
data("maPooling")
pd = pData(maPooling)
individuals = which(rowSums(pd)==1)
# extract individual mice and their strain:
individuals = which(rowSums(pd)==1)
individuals = individuals[-grep("tr", names(individuals))]
y = exprs(maPooling)[,individuals]
g = factor(as.numeric(grepl("b", names(individuals))))
 # apply a t-test to each gene using rowttest function in genefilter
library(genefilter)
tt = rowttests(y,g)
NsigAt01 = sum(tt$p.value<0.01)
NsigAt01
NsigAt05 = sum(tt$p.value<0.05)
NsigAt05
# split first group into 2 and force null to be true
set.seed(0)
shuffledIndex <- factor(sample(c(0,1), sum(g==0), replace=TRUE))
nulltt <- rowttests(y[,g==0], shuffledIndex)
NfalselySigAt01 = sum(nulltt$p.value<0.01)
NfalselySigAt01
NfalselySigAt05 = sum(nulltt$p.value<0.05)
NfalselySigAt05
library(qvalue)
# adjust for false positives reported:
qvals = qvalue(tt$p.value)$qvalue
sum(qvals<0.05)
sum(qvals<0.01)
# now null case generates no false positives:
library(qvalue)
nullqvals = qvalue(nulltt$p.value)$qvalue
sum(nullqvals<0.05)
sum(nullqvals<0.01)


## moderated t tests with limma package
library(SpikeInSubset)
data("rma95")
fac <- factor(rep(1:2, each=3))
# 16 mRNA species, fixed concentration samples prepared and mixed using hgu95 array
# subset is such that for each of spiked in mRNA, first and second trio samples have fixed distinct values
pData(rma95)
# get a feel for response of array quantifications to this design
# rma (robust multi array average) quantifications on a log2 scale:
par(mfrow=c(2,2))
for (i in 1:4){
  spg = names(pData(rma95))
  plot(1:6, exprs(rma95)[spg[i+6],], main=spg[i+6], ylab="RMA",
       xlab="nominal", axes=FALSE)
  axis(2)
  axis(1, at=1:6, labels=pData(rma95)[[spg[i+6]]])
}
# perform simple t tests:
library(genefilter)
rtt <- rowttests(exprs(rma95), fac)
# define colors depending on if p value is small, absolute diff. in means is large, and whether feature is spike in value
mask <- with(rtt, abs(dm)<0.2 & p.value<0.01)
spike <- rownames(rma95) %in% colnames(pData(rma95))
cols <- ifelse(mask, "red", ifelse(spike, "dodgerblue", "black"))
# plot results:
with(rtt, plot(-dm, -log10(p.value), cex=0.8, pch=16,
               xlim=c(-1,1), ylim=c(0,5),
               xlab="difference in means",
               col=cols))
abline(h=2, v=c(-0.2,0.2), lty=2)
# we see red genes have mostly low estimates of std deviation
rtt$s <- apply(exprs(rma95), 1, function(row) sqrt(.5 * (var(row[1:3])+var(row[4:6]))))
with(rtt, plot(s, -log10(p.value), cex=0.8, pch=16,
               log="x", xlab="estimate of standard deviation",
               col=cols))
# perform basic limma analysis
library(limma)
options(digits=3)
fit <- lmFit(rma95, design=model.matrix(~ fac)) # least squares estimate
colnames(coef(fit))
fit <- eBayes(fit) # moderate t statistics
tt <- topTable(fit, coef=2) # report
tt
# topTable return top genes ranked by value defined
# by default method of Behamini-Hochberg is used
dim(topTable(fit, coef=2, number=Inf, sort.by="none"))
# compare previous volcano plot with limma results:
# note that red points are under line where -log10(p.value) = 2
# also blue points represent real differences have p-values higher than before
limmares <- data.frame(dm=coef(fit)[, "fac2"], p.value=fit$p.value[, "fac2"])
with(limmares, plot(dm, -log10(p.value), cex=0.8, pch=16,
                    col=cols, xlab="difference in means",
                    xlim=c(-1,1), ylim=c(0,5)))
abline(h=2, v=c(-0.2,0.2), lty=2)
# construct plot to show how limma shrinks variance est towards common value, eliminating false positives
# we pick for each of 40 bins of different variance estimates, one gene that falls in that bin and remove bins without such genes
n <- 40
qs <- seq(from=0, to=0.2, length=n)
idx <- sapply(seq_len(n), function(i) which(as.integer(cut(rtt$s^2,qs))==i)[1])
idx <- idx[!is.na(idx)]
# plot a line from initial estimate of variance for these genes to estimate after running limma
par(mar=c(5,5,2,2))
plot(1,1, xlim=c(0,0.21), ylim=c(0,1), type="n",
     xlab="variance estimates", ylab="", yaxt="n")
axis(2,at=c(0.1,0.9), c("before", "after"), las=2)
segments((rtt$s^2)[idx], rep(0.1,n),
         fit$s2.post[idx], rep(0.9,n))

## NOTE - WASN'T ABLE TO GET PORTION OF CODE TO WORK. possibly bc it is outdated.
## gene sets (summary statistics)
# preparing data
BiocManager::install(c("sva", "hgfocus.db"))
install.packages("GSA")
library(rafalib)
library(GSEABase)
library(GSE5859Subset)
library(sva)
library(limma)
library(matrixStats)
library(devtools)
data(GSE5859Subset)
X = sampleInfo$group
mod <- model.matrix(~X)
svafit <- sva(geneExpression, mod)
svaX <- model.matrix(~X + svafit$sv)
lmfit <- lmFit(geneExpression, svaX)
tt <- lmfit$coef[,2]*sqrt(lmfit$df.residual)/(2*lmfit$sigma)
pval <- 2*(1-pt(abs(tt), lmfit$df.residual[1]))
qval <- p.adjust(pval, "BH")
library(GSA)
gfiles <- GSA.read.gmt("C:/Users/dgpeb/Downloads/c1.all.v2025.1.Hs.entrez.gmt") # original data not available
length(gfiles$genesets)
head(names(gfiles))
gfiles[["geneset.names"]]
gfiles[["genesets"]][[301]] # gene ids
mapGMT2Affy <- function(object, gfiles){
  ann <- annotation(object)
  dbname <- paste(ann, "db", sep=".")
  require(dbname, character.only=TRUE)
  gns <- featureNames(object)

  map <- select(get(dbname), keys = gns,columns=c("ENTREZID", "PROBEID"))
  map <- split(map[,1],map[,2])
  indexes <- sapply(gfiles, function(ids){
    gns2 <- unlist(map[geneIds(ids)])
    match(gns2, gns)
  })
  names(indexes) <- names(gfiles)
  return(indexes)
}
rownames(sampleInfo) <- colnames(geneExpression)
e = ExpressionSet(assay=geneExpression,
                  phenoData = AnnotatedDataFrame(sampleInfo),
                  annotation = "hgfocus")
gsids <- mapGMT2Affy(e, gfiles[["genesets"]])
class(gfiles[["genesets"]])
# approaches based on association tests
tab <- table(ingeneset=1:nrow(e) %in% gsids[["chrYq11"]], signif=qval<0.05)
tab # object gsids not found - tab won't work.
# gene set summary statistics
## END CODE THAT DOESN'T WORK

## gene testing in r - part 1 and 2
library(GEOquery)
g <- getGEO("GSE34313")
e <- g[[1]]
e$condition <- e$characteristics_ch1.2
e$condition
table(e$condition)
# examine boxplots - subset control samples and samples treated with dexamethasone
boxplot(exprs(e), range=0)
names(fData(e))
levels(e$condition) <- c("treatment: dexamethasone for 24 hr", "treatment: dexamethasone for 4 hr", "treatment: none")
lvls <- c("treatment: dexamethasone for 4 hr", "treatment: none")
es <- e[,e$condition %in% lvls]
es$condition <- factor(es$condition, levels=lvls)
es$condition
# run linear model in limma
library(limma)
design <- model.matrix(~ es$condition)
fit <- lmFit(es, design=design)
fit <- eBayes(fit)
tt <- topTable(fit, coef=2, genelist=fData(es)$GENE_SYMBOL)
tt
# immune response:
idx <- grep("GO:0006955", fData(es)$GO_ID)
length(idx)
r1 <- roast(es, idx, design)
r1
# testing multiple gene sets:
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
org.Hs.egGO2EG
go2eg <- as.list(org.Hs.egGO2EG)
head(go2eg)
# unlist the list and match Entrez ID to index in ExpressionSet:
govector <- unlist(go2eg)
golengths <- sapply(go2eg, length)
head(fData(es)$GENE)
idxvector <- match(govector, fData(es)$GENE)
table(is.na(idxvector))
idx <- split(idxvector, rep(names(go2eg), golengths))
go2eg[[1]]
fData(es)$GENE[idx[[1]]]
# clean lists:
idxclean <- lapply(idx, function(x) x[!is.na(x)])
idxlengths <- sapply(idxclean, length)
idxsub <- idxclean[idxlengths > 10]
length(idxsub)
# run multiple ROAST test
r2 <- mroast(es, idxsub, design) # doesn't work - no residual degrees of freedom
head(r2)
r2 <- r2[order(r2$PValue.Mixed),]
# extract GO terms for top results by mixed test
BiocManager::install("GO.db")
library(GO.db)
columns(GO.db)
keytypes(GO.db)
GOTERM[[rownames(r2)[1]]]
r2tab <- select(GO.db, keys=rownames(r2)[1:10],
                columns=c("GOID", "TERM", "DEFINITION"),
                keytype="GOID")
r2tab[,1:2]
# look for top results using standard p value in up direction:
r2 <- r2[order(r2$PValue),]
r2tab <- select(GO.db, keys=rownames(r2)[r2$Direction=="Up"][1:10],
                columns = c("GOID", "TERM", "DEFINITION"),
                keytype="GOID")
r2tab[,1:2]
# for down direction:
r2tab <- select(GO.db, keys=rownames(r2)[r2$Direction=="Down"][1:5],
                columns = c("GOID", "TERM", "DEFINITION"),
                keytype="GOID")
r2tab[,1:2]