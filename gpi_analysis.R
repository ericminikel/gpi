options(stringsAsFactors=FALSE)
setwd('~/d/sci/src/gpi')

genes = read.table('data/gpi_anchored.tsv',sep='\t',header=F)
colnames(genes) = c('symbol')

# random spot checks - Google a few & see what is known
set.seed(1)
sample(genes$symbol, size=10)
# result:
# "FBLN1"    "GPC1"     "DNASE1L1" "VDAC3"    "SLC7A5"   "NCAM1"    "BASP1"    "TMPO"     "GPC6"     "SLC25A5" 
