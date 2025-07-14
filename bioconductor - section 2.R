BiocManager::install(c("Homo.sapiens",
                       "GenomicFeatures",
                       "genomicsclass/ERBS",
                       "genomicsclass/ph525x"))
library(base)

## motivation sample code
install.packages("devtools")
library(devtools)
install_github("genomicsclass/ERBS",
               force=TRUE)
library(ERBS)
data(HepG2)
HepG2 # one of the cell lines
# HepG2 is summarized in Bioconductor using a GRanges object

## introduction to using genomicranges
library(ERBS)
library(GenomicRanges)
data(HepG2)
data(GM12878)
# inspect HepG2
class(HepG2)
HepG2
values(HepG2)
# extract chromosome names
seqnames(HepG2)
chr = seqnames(HepG2)
as.character(chr)
# table of numbers of sequences on each chromosome
table(chr)
table(chr)[1:24] # restrict to autosomes, X and Y
#subset and order GRanges
HepG2[chr=="chr20",]
x = HepG2[order(HepG2),]
x
seqnames(x)
as.character(seqnames(x))