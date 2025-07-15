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



## introduction to using GenomicRanges
library(ERBS)
library(GenomicRanges)
data(HepG2)
data(GM12878)
# inspect HepG2
class(HepG2)
HepG2
values(HepG2)
# can treat this data like a matrix
HepG2[1:10,]
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



## interval ranges: IRanges
library(IRanges)
ir <- IRanges(5,10)
ir # gives us start, end, and width
IRanges(5, width=6) # identical result, tells us start is 5, it is 6 bp long
start(ir)
end(ir)
width(ir)
# can specify more than one range at a time
ir <- IRanges(start=c(3,5,17), end=c(10,8,20))
ir
length(ir)
start(ir)
end(ir)
width(ir)
# intra range methods = operation will occur for each range you have,
# doesn't depend on other ranges contained in IRanges object
ir <- IRanges(5,10)
shift(ir, -2)
ir
# narrow = relative to start, instead start this range at second base pair
narrow(ir, start=2)
# relative to 5, should end on 5th base pair, i.e. 9
narrow(ir, end=5)
# flank = gets flanking sequence
# flanking sequence 3 bp from start
flank(ir, width=3, start=TRUE, both=FALSE)
# flanking sequence 3 bp from end
flank(ir, width=3, start=FALSE, both=FALSE)
# bidirectional flanking sequence
flank(ir, width=3, start=TRUE, both=TRUE)

# setting up a plot to look at range operations
# here we are drawing base pairs
# number indicates the base pair
# we are drawing the original range
plotir <- function(ir, i){arrows(start(ir)-.5, i, end(ir)+.5, i, code=3,
                                 angle=90, lwd=3)}
plot(0, 0, xlim=c(0,15), ylim=c(0,8), type="n", xlab="", ylab="", xaxt="n")
axis(1,0:15)
abline(v=0:30+.5, col=rgb(0, 0, 0, .5))
# plot original IRange:
plotir(ir,1)
# draw red shadow box for original IRange
polygon(c(start(ir)-.5, start(ir)-.5, end(ir)+.5, end(ir)+.5),
        c(-1,9,9,-1), col=rgb(1,0,0,.2), border=NA)
plotir(shift(ir,-2), 2)
plotir(narrow(ir, start=2), 3) # start at 2nd bp not 1st
plotir(narrow(ir, end=5), 4) # specify end at 5
plotir(flank(ir, width=3, start=TRUE, both=FALSE), 5) # 3 bp from start to left
plotir(flank(ir, width=3, start=FALSE, both=FALSE), 6) # flank sequence downstream of range
plotir(flank(ir, width=3, start=TRUE, both=TRUE), 7) # 3 bp to left and right of start
# these are the most basic intra range methods

# inter range methods = functions depend on other ranges in object
ir <- IRanges(start=c(3,5,17), end=c(10,8,20)) # create IRange object w 3 ranges
ir
range(ir) # gives beginning of IRanges to end, including gaps
reduce(ir) # gives basepairs covered by original ranges
# with this we don't get the gap at the end of 10 and beginning of 17
# we can then ask for the gaps:
gaps(ir)
# gap therefore begins at 11, ends at 16
# get set of ranges that has same coverage as original IRange object:
disjoin(ir)
# same coverage but they don't overlap
# also contains union of all endpoints of original range



## genomic ranges: GRanges
library(GenomicRanges)
# IRange on Chromosome Z, can contain strand info and seq lengths:
gr <- GRanges("chrZ", IRanges(start=c(5,10), end=c(35,45)),
              strand="+", seqlengths=c(chrZ=100L))
gr
# we have two ranges, zero metadata columns
# we specified chromosome z is 100 bp long
shift(gr, 10)
# moves starts and ends by 10 bp to right
shift(gr, 80)
# this goes off end of chromosome, gives us an error
# wrap this in a trim function, makes sure it ends at chromosome end: 
trim(shift(gr, 80))
mcols(gr) #access metadata columns
mcols(gr)$value <- c(-1,4) # add columns with $
gr
# additional class in GRanges package: GRangesList
# GRanglesList groups GRanges together
gr2 <- GRanges("chrZ", IRanges(11:13,51:53))
mcols(gr)$value <- NULL
# create new GRangesList:
grl <- GRangesList(gr, gr2)
# contains 2 GRanges - first GRanges has 2 ranges, second has 3 ranges:
grl
length(grl)
grl[[1]]
mcols(grl)$value <- c(5,7)
grl
mcols(grl)
# finding overlaps:
gr1 <- GRanges("chrZ", IRanges(c(1,11,21,31,41), width = 5))
gr2 <- GRanges("chrZ", IRanges(c(19,33), c(38,35)))
# both GRanges are on same sequence, chromosome Z
gr1
gr2
fo <- findOverlaps(gr1, gr2)
fo
# hits tells us how many overlaps occurred
# [1] 3rd element of query (gr1) intersected with 1st element of subject (gr2)
# [2] 4th element of query intersected with 1st element of subject
# [3] 4th element of query intersected with 2nd element of subject
queryHits(fo)
subjectHits(fo)
# another way to find overlap:
gr1 %over% gr2 # gives logical vector
gr1[gr1 %over% gr2] # logical subsetting returns ranges in gr1 with overlap

# introductioon to Rle: 
r <- Rle(c(1,1,1,0,0,-2,-2,-2,rep(-1,20)))
r
str(r)
as.numeric(r)
# create 2 views, one that starts at 4 and ends at 7, other starts at 2 and ends at 6
Views(r, start=c(4,2), end=c(7,6))



## operating on GRanges
BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg19",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "ph525x",
  "Gviz"))
library(Biostrings)
library(knitr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Homo.sapiens)
library(BSgenome)
library(Gviz)

library(IRanges)
ir <- IRanges(c(3,8,14,15,19,34,40), 
              width=c(12,6,6,15,6,2,7))
# set up plotRanges function
plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),
                       col="black", sep = 0.5, ...)
{
  height <- 1
  if (is(xlim, "Ranges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col,...)
  title(main)
  axis(1)
}

plotRanges(ir)
par(mfrow=c(4,1), mar=c(4,2,2,2))
plotRanges(ir, xlim=c(0,60))
plotRanges(reduce(ir), xlim=c(0,60))
plotRanges(disjoin(ir), xlim=c(0,60))
plotRanges(gaps(ir), xlim=c(0,60))

library(GenomicRanges)
gir = GRanges(seqnames="chr1")
gir
# we need to supply metadata
strand(gir) = c(rep("+", 4), rep("-", 3))
gir
genome(gir) = "hg19"
seqinfo(gir)

gir = GRanges(seqnames="chr1", ir, strand=c(rep("+", 4),
                                            rep("-", 3)))
# set up plotGRanges function
plotGRanges <- function(x, xlim=x, col="black", sep=0.5,
                        xlimits=c(0,60),...)
{
  main=deparse(substitute(x))
  ch = as.character(seqnames(x)[1])
  x = ranges(x)
  height <- 1
  if(is(xlim, "Ranges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim = xlimits, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x) - 0.5, ybottom, end(x) + 0.5, ybottom + height,
       col = col, ...)
  title(main, xlab=ch)
  axis(1)
}
par(mfrow=c(4,1), mar=c(4,2,2,2))
plotGRanges(gir, xlim=c(0,60))
plotGRanges(resize(gir, 1), xlim=c(0,60), col="green")
plotGRanges(flank(gir, 3), xlim=c(0,60), col="purple")
plotGRanges(flank(gir, 2, start=FALSE), xlim=c(0,60), col="brown")



## finding overlaps
# load libraries in
library(GenomicFeatures)
library(GenomicRanges)
library(IRanges)
library(ERBS)
# cell lines:
data(HepG2)
data(GM12878)
browseVignettes("GenomicRanges")
# investigate sites in both cell lines:
res = findOverlaps(HepG2, GM12878)
res
class(res)
HepG2[1,]
GM12878[12,]
queryHits(res)
# index this
index = queryHits(res)
# subset HepG2, create new function called erbs
erbs = HepG2[index,]
erbs
# we just want location:
granges(erbs)



## genes as granges
# load up gene info:
library(Homo.sapiens)
ghs = genes(Homo.sapiens)
ghs
# precede function:
index = precede(erbs, ghs)
index
# find closest preceding range in ghs:
ghs[index[1:3]]
erbs[1:3]



## finding the nearest gene
# distance between binding sites and nearest preceding genes:
distance(erbs, ghs[index])
# find transcription start site nearest to each binding site:
tssgr = resize(ghs, 1)
tssgr
# distance between binding site and nearest tss:
d = distanceToNearest(erbs, tssgr)
d
queryHits(d)
dists = values(d)$distance
dists
# histograms of distances:
hist(dists, nc=100)
hist(dists, nc=1000)
hist(dists, nc=100, xlim=c(0,100000))
hist(dists, nc=1000, xlim=c(0,100000))
# take a closer look at the genes close to our tss
index = subjectHits(d)
index
index = subjectHits(d)[dists<1000]
index



## annotating genes and querying with select
# use select function to query hono sapiens database:
tssgr[index,]
keytypes(Homo.sapiens)
keys = as.character(values(tssgr[index])$GENEID)
columns(Homo.sapiens)
res = select(Homo.sapiens, keys = keys,
             columns = c("SYMBOL", "GENENAME"), keytype="GENEID")
res
res[1:2,]



## dnastring object
library(Biostrings)
# define a dnastring
dna <- DNAString("TCGAGCAAT")
dna
length(dna)
# invalid base:
DNAString("JQX")
# valid sequence with unknowns and gaps:
DNAString("NNNACGCGC-TTA-CGGGCTANN")
# index into a dnastring:
dna[4:6]
# convert dnastring back to character string:
as.character(dna)
# combine multiple dna sequences:
set1 <- DNAStringSet(c("TCA", "AAATCG", "ACGTGCCTA", 
                       "CGCGCA", "GTT", "TCA"))
set1
set1[2:3] # extract sequence subset
set1[[4]] # extract one sequence as a single dnastring
length(set1) # returns set length not size of each sequence
width(set1) # returns size of each individual sequence
duplicated(set1) # tells us which sequences are duplicated
unique(set1) # keeps only unique sequences
sort(set1) # sort sequences alphabetically 
# consider the dna sequence:
dna_seq <- DNAString("ATCGCGCGCGGCTCTTTTAAAAAAACGCTACTACCATGTGTGTCTATC")
letterFrequency(dna_seq, "A") # count times specific letter appears
letterFrequency(dna_seq, "GC") # counts times letters appear together
# can see frequency of all dinucleotides and trinucleotides
dinucleotideFrequency(dna_seq)
trinucleotideFrequency(dna_seq)
# reverse complement of dnastring:
reverseComplement(dna_seq)
# amino acid translation of dnastring:
translate(dna_seq)
# count number of occurences of pattern and find location of those patterns
dna_seq <- DNAString("ATCGCGCGCGGCTCTTTTAAAAAAACGCTACTACCATGTGT")
dna_seq
countPattern("CG", dna_seq)
matchPattern("CG", dna_seq)
start(matchPattern("CG", dna_seq))
matchPattern("CTCTTTTAAAAAAACGCTACTACCATGTG", dna_seq)
# dna is double stranded so we may want to check for reverse complement:
countPattern("TAG", dna_seq)
countPattern(reverseComplement(DNAString("TAG")), dna_seq)
# can also count and locate patterns in dnastring objects:
set2 <- DNAStringSet(c("AACCGGTTTCGA", "CATGCTGCTACA",
                       "CGATCGCGCCGG", "TACAACCGTACA"))
set2
# count number of occurences of pattern across many biostrings:
vcountPattern("CG", set2)
# return locations of pattern across many biostrings:
vmatchPattern("CG", set2)
# work with matches from a single string in a set:
vmatchPattern("CG", set2)[[1]]

eco <- DNAString("GGTTTCACCGCCGGTAATGAAAAAGGCGAACTGGTGGTGCTTGGACGCAACGGTTCCGACTACTCTGCTGCGGTGCTGGCTGCCTGTTTACGCGCCGATTGTTGCGAGATTTGGACGGACGTTGACGGGGTCTATACCTGCGACCCGCGTCAGGTGCCCGATGCGAGGTTGTTGAAGTCGA")
eco
dinucleotideFrequency(eco)
trinucleotideFrequency(eco)
translate(eco)
reverseComplement(eco)
countPattern("ATG", eco)
countPattern("TATA", eco)



## getting the sequence of regions
library(ERBS)
data(HepG2)
HepG2
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
Hsapiens
c17 = Hsapiens$chr17
c17
class(Hsapiens)
hepseq = getSeq(Hsapiens, HepG2)
hepseq
length(HepG2)
width(HepG2)[1:5]
# shift binding regions so there is no principled relationship to binding peaks:
rhepseq = getSeq(Hsapiens, shift(HepG2, 2500)) 
rhepseq
mot = "TCAAGGTCA" # the motif we are investigating
vcountPattern(mot, hepseq)
# frequency of motif appearance:
sum(vcountPattern(mot, hepseq))
sum(vcountPattern(mot, reverseComplement(hepseq)))
sum(vcountPattern(mot, reverseComplement(hepseq)))+
  sum(vcountPattern(mot, hepseq))
sum(vcountPattern(mot, reverseComplement(rhepseq)))+
  sum(vcountPattern(mot, rhepseq))
# there is roughly a 9 fold enrichment of sequences under the binding peaks of ER protein
