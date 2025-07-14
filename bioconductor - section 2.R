BiocManager::install(c("Homo.sapiens",
                       "GenomicFeatures",
                       "genomicsclass/ERBS",
                       "genomicsclass/ph525x"))
library(base)

## motivation sample code
install.packages("devtools")
install.packages("RTools45")
library(devtools)
install_github("genomicsclass/ERBS")
library(ERBS)
data(HepG2)
HepG2 # one of the cell lines
# HepG2 is summarized in Bioconductor using a GRanges object

## introduction to using genomicranges
library(ERBS)
data(HepG2)
data(GM12878)
# inspect HepG2
class(HepG2)
HepG2
values(HepG2)
# extract chromosome names
seqnames(HepG2)
