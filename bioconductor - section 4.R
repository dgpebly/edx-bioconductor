BiocManager::install(c("BSgenome",
                       "BSgenome.Hsapiens.UCSC.hg19.masked",
                       "org.Hs.eg.db",
                       "ensembldb",
                       "EnsDb.Hsapiens.v75",
                       "AnnotationHub",
                       "rtracklayer",
                       "TxDb.Hsapiens.UCSC.hg38.knownGene",
                       "KEGGREST",
                       "rols",
                       "GSEABase"))
install.packages(c("R.utils",
                   "png",
                   "DT"))


## discovering availible reference genomes:
library(Biostrings)
library(BSgenome)
ag = available.genomes()
length(ag)
head(ag) # get names of reference genomic sequences
# test commands
grep("Scerev", ag, value=TRUE)
grep("Hsap", ag, value=TRUE)
# inspect human genome:
library(BSgenome.Hsapiens.UCSC.hg19)
length(Hsapiens)
class(Hsapiens)
methods(class="BSgenome")
Hsapiens$chrX
substr(Hsapiens$chrX, 5e6, 5.1e6) # moving in 5 megabases and we see familiar nucleotide codes
nchar(Hsapiens$chrY)
nchar(Hsapiens[[24]])
sum(unlist(lapply(18:24, function(x) nchar(Hsapiens[[x]]))))
system.time(sum(unlist(lapply(18:24, function(x) nchar(Hsapiens[[x]])))))
# improve response time:
library(parallel)
detectCores()
options(mc.cores=16)
system.time(sum(unlist(mclapply(18:24, function(x) nchar(Hsapiens[[x]])))))


## packages for gene and transcript catalogs:
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
class(txdb)
methods(class="TxDb")
# extract and inspect genes from TxDb
genes(txdb)
table(strand(genes(txdb)))
summary(width(genes(txdb)))
# inspect largest gene in genome
id = which.max(width(genes(txdb)))
genes(txdb)[id]
library(org.Hs.eg.db)
select(org.Hs.eg.db, keys="286297", keytype="ENTREZID", 
       columns=c("SYMBOL","GENENAME"))
# compare total size of exons to total size of genes:
ex = exons(txdb)
rex = reduce(ex)
ex_width = sum(width(rex)) # bases in exons
gene_width = sum(width(genes(txdb))) # bases in genes
ex_width/gene_width


## ensembldb, EnsDb: annotation from EMBL:
# inspect data available from Ensembl
library(ensembldb)
library(EnsDb.Hsapiens.v75)
names(listTables(EnsDb.Hsapiens.v75))
# extract Ensembl transcripts
edb = EnsDb.Hsapiens.v75
txs <- transcripts(edb, filter=GeneNameFilter("ZBTB16"),
                   columns=c("protein_id", "uniprot_id", "tx_biotype"))
txs
# compare Ensembl and UCSC transcripts
alltx = transcripts(edb) # Ensembl is larger
utx = transcripts(txdb) # UCSC is smaller
# table of biological types of transcripts
table(alltx$tx_biotype)


## annotationhub: finding and caching important annotation:
library(AnnotationHub)
ah <- AnnotationHub()
ah
length(unique(ah$species))
ah_human <- subset(ah, species == "Homo sapiens")
ah_human
query(ah, "HepG2")
query(ah, c("HepG2", "H3K4me3"))
hepg2 <- query(ah, "HepG2")
hepg2_h3k4me3 <- query(ah, "H3K4me3")
hepg2_h3k4me3
hepg2_h3k4me3$tags
e118_broadpeak <- query(hepg2_h3k4me3, c("E118", "broadPeak"))
id <- e118_broadpeak$ah_id
id
hepg2_h3k4me3_broad <- ah[["AH29728"]]
hepg2_h3k4me3_broad
alt_format <- ah[[id]]
identical(hepg2_h3k4me3_broad, alt_format)


## liftover: translating between reference builds
library(rtracklayer)
# chromosome 1 gene locations in hg38:
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
tx38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
seqlevels(tx38, pruning.mode="coarse")="chr1"
g1_38 = genes(tx38)
# download hg38 to hg19 chain file:
library(AnnotationHub)
ah <- AnnotationHub()
ah.chain <- subset(ah, rdataclass == "ChainFile" 
                   & species == "Homo sapiens")
query(ah.chain, c("hg19", "hg38"))
ch <- ah[["AH14108"]]
# perform the liftover:
g1_19L <- liftOver(g1_38, ch)
g1_19L


## data import and export with rtracklayer
library(devtools)
install_github("genomicsclass/ERBS") # install ERBS package
f1 = dir(system.file("extdata", package="ERBS"), full=TRUE)[1] # access data
readLines(f1, 4) # look at a few lines
library(rtracklayer)
imp = import(f1, format="bedGraph") # import as bedGraph format
imp
genome(imp) # gene identifier tag not set, but should be set
genome(imp) = "hg19" # set genome
genome(imp)
export(imp, "demoex.bed") # export as BED format
cat(readLines("demoex.bed", n=5), sep="\n") # check output file
