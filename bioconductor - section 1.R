## part one
## installing bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("3.15")
BiocManager::version()
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install(c("genefilter", "geneplotter"))
library(BSgenome.Hsapiens.UCSC.hg19)
library(genefilter)
library(geneplotter)
## part two
## getting help in r and bioconductor
help.start()
?mean
help(mean)
help(package = "genefilter")
# inspect objects, classes, and methods
library(Biobase) #loading one of core Bioconductor package
?ExpressionSet
?"ExpressionSet-class"
methods(class = ExpressionSet)
# inspect source code for functions and methods
read.csv
plotMA
showMethods("plotMA")
getMethod("plotMA","data.frame")
# vignettes teach how to use functions in a package
vignette(package = "Biobase")
vignette("ExpressionSetIntroduction")
browseVignettes(package = "Biobase")
#report key details about R session
sessionInfo()

## section 1
BiocManager::install(version="3.21")
BiocManager::install(c("genefu",
                       "COPDSexualDimorphism",
                       "gwascat",
                       "hgu133a.db",
                       "genomicsclass/tissuesGeneExpression"))
library(genefu)
data(sig.gene70)
dim(sig.gene70)
head(sig.gene70)[,1:6]

library(gwascat)
data("ebicat_2020_04_30")
ebicat_2020_04_30

BiocManager::install("genomicsclass/tissuesGeneExpression", 
                     force=TRUE)
BiocManager::install("remotes")
library(tissuesGeneExpression)
data(tissuesGeneExpression)
head(e[,1:5])
table(tissue)