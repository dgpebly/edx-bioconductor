BiocManager::install(c("Biobase",
                       "GEOquery",
                       "genomicsclass/GSE5859Subset",
                       "affy",
                       "hgu95acdf",
                       "genefilter",
                       "parathyroidSE",
                       "airway",
                       "pasillaBamSubset",
                       "Rsamtools",
                       "GenomicAlignments",
                       "ArrayExpress",
                       "NGScopyData",
                       "AnnotationDbi"))
## expressionset class:
library(knitr)
library(Biobase)
library(GEOquery)
geoq <- getGEO("GSE9514")
names(geoq)
e <- geoq[[1]] # extract ExpressionSet
# exprs gives expression matrix:
dim(e) # number of features and samples in ExpressionSet
exprs(e)[1:3,1:3]
dim(exprs(e)) # rows = features, columns = samples
#pData gives phenotype data (sample info):
pData(e)[1:3,1:6]
dim(pData(e)) # rows of pData correspond to columns of exprs
names(pData(e))
pData(e)$characteristics_ch1
# fData gives feature data (probe info):
fData(e)[1:3,1:3]
dim(fData(e)) # rows correspond to rows of exprs
names(fData(e))
head(fData(e)$"Gene Symbol")
head(rownames(e))
# additional annotation tied to ExpressionSet:
experimentData(e)
annotation(e)


## reading microarray raw data single color arrays:
wd <- getwd()
datadir <- paste0(wd, "/rawdata-master")
basedir <- paste0(datadir, "/celfiles")
setwd(basedir)
library(affy)
tab <- read.delim("sampleinfo.txt", check.names=FALSE, as.is=TRUE)
rownames(tab) <- tab$filenames
tab
fns <- list.celfiles(basedir)
fns
fns %in% tab[,1] ## check
ab <- ReadAffy(phenoData=tab)
# creates affybatch object which contains needed info
dim(pm(ab))
dim(pData(ab))
rownames(pData(ab))
annotation(ab)
# preprocess rma:
e <- rma(ab)
# go back to previous working directory:
setwd(wd)
# if not interested in probe level data:
setwd(basedir)
ejust <- justRMA(filenames=tab[,1], phenoData=tab)
dim(ejust)
# agilent data:
library(limma)
install.packages("rafalib")
library(rafalib)
basedir <- paste0(datadir, "/agilent")
setwd(basedir)
targets <- readTargets("TargetBeta7.txt")
RG <- read.maimages(targets$FileName, source="genepix")
MA <- MA.RG(RG, bc.method="none")
mypar(1,1)
imageplot(MA$M[,2], RG$printer, zlim=c(-3,3))
dev.off()
# can use oligo to read affy arrays:
detach("package:affy")
library(oligo)
basedir <- paste0(datadir, "/celfiles")
setwd(basedir)
tab <- read.delim("sampleinfo.txt", check.names=FALSE, as.is=TRUE)
fns <- list.celfiles(listGzipped=TRUE)
fns %in% tab[,1] ## check
pd <- as(tab, "AnnotatedDataFrame")
efs <- read.celfiles(filenames=tab[,1], phenoData=pd,
                     sampleNames=sampleNames(pd))
e <- rma(efs)


## reading microarray raw data agilent two color arrays:
library(limma)
library(rafalib)
basedir <- paste0(datadir, "/agilent")
setwd(basedir)
targets <- readTargets("TargetBeta7.txt")
RG <- read.maimages(targets$FileName, source="genepix")
MA <- MA.RG(RG, bc.method="none")
dim(RG$R)
dim(RG$G)
dim(MA$M)
dim(MA$A)
plot(MA$A[,1], MA$M[,1]) # MA plot for first sample
# microarray image
mypar(1,1)
imageplot(MA$M[,2], RG$printer, zlim=c(-3,3))
dev.off()


## the summarizedexperiment class
library(parathyroidSE)
data("parathyroidGenesSE")
se <- parathyroidGenesSE
se
dim(se)
assay(se)[1:3,1:3]
dim(assay(se))
# column data - equivalent to pdata in expressionset
colData(se)[1:3,1:6]
dim(colData(se))
names(colData(se))
colData(se)[1]
colData(se)$treatment
# row data 
rowRanges(se)[1]
# double brackets extract out single GRanges object:
rowRanges(se)[[1]]
class(rowRanges(se))
length(rowRanges(se)) # number of genes
length(rowRanges(se)[[1]]) # number of exons for first gene
# metadata tells us how this GRanges list was constructed
head(rowRanges(se))
metadata(rowRanges(se))
# more info about experiment
metadata(se)$MIAME
abstract(metadata(se)$MIAME)


## importing ngs data in r
library(pasillaBamSubset)
library(Rsamtools)
filename <- untreated1_chr4()
# create BamFile object to allow other functions to process file:
bf <- BamFile(filename)
# ask for info on chromosomes declared in header of BAM file:
seqinfo(bf)
sl <- seqlengths(bf)
# summary of alignment types in file:
quickBamFlagSummary(bf)
# count number of reads on chr. 4
gr <- GRanges("chr4", IRanges(1, sl["chr4"]))
countBam(bf, param = ScanBamParam(which=gr))
# specify new bamfile and limit reads to 5 at a time:
reads <- scanBam(BamFile(filename, yieldSize=5))
# examining output of scanbam:
class(reads)
names(reads[[1]])
reads[[1]]$pos # aligned start position
reads[[1]]$rname # chromosome
reads[[1]]$strand # the strand
reads[[1]]$qwidth # width of the read
reads[[1]]$seq # seq of read
# example of specifiying what and which:
gr <- GRanges("chr4", IRanges(500000,700000))
reads <- scanBam(bf, param=ScanBamParam(what=c("pos", "strand"), 
                                        which=gr))
## examine output of readGAlignments():
library(GenomicAlignments)
ga <- readGAlignments(bf)
length(ga)
# can extract GRanges object within GAlignments object
granges(ga[1])
# can use familiar GenomicRanges functions on GAlignments:
gr <- GRanges("chr4", IRanges(700000,800000))
(fo <- findOverlaps(ga,gr)) # reads over this range
countOverlaps(gr,ga) # count overlaps in range with reads
table(ga %over% gr) # logical vector of read overlaps in range
# integer vector with overlaps for each read in range in gr:
countOverlaps(ga,gr)


## creating a count table from a BAM file
# load packages:
BiocManager::install("TxDb.Dmelanogaster.UCSC.dm3.ensGene")
library(pasillaBamSubset)
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
# pull out exons for genes:
grl <- exonsBy(txdb, by="gene")
grl[100] #GRangesList of exons for 100th gene
grl[[100]] #Granges with exons of 100th gene
grl[[100]][1] # first exon of 100th gene
# paths to BAM files:
fl1 <- untreated1_chr4()
fl2 <- untreated3_chr4()
# libraries for importing BAM files:
library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
# specify files with BamFileList
fls <- BamFileList(c(fl1,fl2))
names(fls) <- c("first", "second")
# find reads that overlap exons
so1 <- summarizeOverlaps(features=grl,
                         read=fls,
                         ignore.strand=TRUE)
so1
# examine count matrix:
head(assay(so1))
colSums(assay(so1))
# examine rest of SummarizedExperiment components:
rowRanges(so1)
colData(so1)
colData(so1)$sample <- c("one", "two") # add sample info
colData(so1)
metadata(rowRanges(so1))
# exploratory data analysis of counts:
x <- assay(so1)[,1]
hist(x[x>0], col="grey")
hist(x[x>0 & x<1000], col="grey")
plot(assay(so1)+1, log="xy")
# count second file as paired-end reads
fls <- BamFileList(fl2)
so2 <- summarizeOverlaps(features = grl,
                         read=fls,
                         ignore.strand=TRUE,
                         singleEnd = FALSE,
                         fragments=TRUE)
colSums(assay(so2))
colSums(assay(so1))
# show there are half as many reads in so2 as so1:
plot(assay(so1)[,2], assay(so2)[,1], xlim=c(0,5000), ylim=c(0,5000),
     xlab="single end counting", ylab="paired end counting")
abline(0,1)
abline(0,.5)


## accessing GEOquery microarray experiments
library(GEOquery)
glioMA = getGEO("GSE78703")[[1]]
glioMA
names(pData(glioMA)) # variable names
glioMA$molecule_ch1 # molecule being assayed (RNA)
table(glioMA$'treated with:ch1', glioMA$'cell type:ch1') # experimental variables


## array oriented annotation
BiocManager::install(c("hgu133plus2.db",
                     "hgu133plus2probe"))
library(hgu133plus2.db)
hgu133plus2.db
library(hgu133plus2probe)
head(hgu133plus2probe)
dim(hgu133plus2probe)
select(hgu133plus2.db, keytype="PROBEID",
       columns=c("SYMBOL", "GENENAME", "PATH", "MAP"), 
       keys="1007_s_at")
