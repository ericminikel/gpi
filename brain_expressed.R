options(stringsAsFactors=FALSE)
library(gplots)
setwd('~/d/sci/src/gpi')
gtex = read.table('data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct',sep='\t',header=T, skip=2)
dim(gtex)
orig_colnames = read.table('data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct',sep='\t',header=F, skip=2)[1,]
colnames(gtex) = gsub('[^a-z0-9_]','_',tolower(colnames(gtex)))
brain_cols = colnames(gtex)[grepl('^brain',colnames(gtex))]
head(gtex[,1:3])
range(gtex[1:2500,3:25]) # 0 to 7739, so this appears to be raw tpm not log scale
gtex$brain_mean = apply(gtex[,brain_cols],1,FUN=mean)
head(gtex[,c(brain_cols,'brain_mean')])
sum(gtex$brain_mean > 10) # 7660

brain_expressed_genes_10tpm = data.frame(symbol=sort(unique(gtex[gtex$brain_mean > 10, 'description'])))

source('~/d/sci/src/exac_papers/exac_constants.R')
kstraint = load_constraint_data(reload=FALSE) # constraint is a SQLite keyword, so I call the table kstraint
kstraint$lof_obs_exp = kstraint$n_lof / kstraint$exp_lof

gpiaps = read.table('../gene_lists/lists/gpi_anchored.tsv',sep='\t',header=F)
colnames(gpiaps) = c('symbol')

gpiaps_empiric = read.table('data/gpi_anchored_empirical.tsv',sep='\t',header=F)
colnames(gpiaps_empiric) = c('symbol')


kstraint$gpiap = kstraint$gene %in% gpiaps$symbol
kstraint$gpiap_empiric = kstraint$gene %in% gpiaps_empiric$symbol

ctable = table(kstraint[,c('gpiap','gpiap_empiric')])
fisher.test(as.matrix(ctable), alternative='two.sided')
fisher.test(as.matrix(ctable), alternative='two.sided')$p.value

colnames(gtex)
colnames(kstraint)

h_breaks = 10^(seq(1,6,by=.25))
all_bp_hist = hist(kstraint$bp, breaks=h_breaks, plot=FALSE)
gpi_bp_hist = hist(kstraint$bp[kstraint$gene %in% gpiaps$symbol], breaks=h_breaks, plot=FALSE)

plot(all_bp_hist)


png('~/d/j/cureffi/media/2018/07/gpi-vs-all-length.png',width=800,height=500,res=150)
par(mar=c(4,5,4,5))
plot(all_bp_hist$mids, all_bp_hist$counts, type='h', lwd=40, lend=1, log='x', xlim=c(100, 10000), col='#77777788', xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
axis(side=1, at=10^(2:4))
axis(side=2, col='#777777', col.axis='#777777', at=(0:6)*1000, las=2, font=2)
par(new=TRUE)
plot(gpi_bp_hist$mids, gpi_bp_hist$counts, type='h', lwd=40, lend=1, log='x', xlim=c(100, 10000), col='#FF991288', xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
axis(side=4, col='#FF9912', col.axis='#FF9912', font=2, at=(0:5)*10, las=2)
mtext(side=1, line=2.5, text='coding sequence length (bp)')
mtext(side=2, line=3.0, text='count (all genes)', col='#777777', font=2)
mtext(side=4, line=2.5, text='count (GPI-anchored proteins)', col='#FF9912', font=2)
mtext(side=3, line=1, text='length distribution\nGPI-anchored proteins vs. all genes', font=2, cex=1.2)
dev.off()

ks.test(kstraint$bp[kstraint$gene %in% gpiaps$symbol], kstraint$bp[!(kstraint$gene %in% gpiaps$symbol)], alternative='two.sided')
mean(kstraint$bp[kstraint$gene %in% gpiaps$symbol])
mean(kstraint$bp[!(kstraint$gene %in% gpiaps$symbol)])


mean(kstraint$n_exons[kstraint$gene %in% gpiaps$symbol])
mean(kstraint$n_exons[!(kstraint$gene %in% gpiaps$symbol)])

sum(kstraint$n_exons[kstraint$gene %in% gpiaps$symbol]==1)
kstraint$gene[(kstraint$n_exons==1 & kstraint$gene %in% gpiaps$symbol)]
mean(kstraint$n_exons[kstraint$gene %in% gpiaps$symbol]==1)
mean(kstraint$n_exons[!(kstraint$gene %in% gpiaps$symbol)]==1)


ks.test(kstraint$lof_obs_exp[kstraint$gene %in% gpiaps$symbol], kstraint$lof_obs_exp[!(kstraint$gene %in% gpiaps$symbol)], alternative='two.sided')
mean(kstraint$lof_obs_exp[kstraint$gene %in% gpiaps$symbol])
mean(kstraint$lof_obs_exp[!(kstraint$gene %in% gpiaps$symbol)])
m = lm(lof_obs_exp ~ bp + gpiap, data=kstraint)
summary(m)
kstraint[kstraint$gpiap & kstraint$lof_obs_exp == 0 & kstraint$exp_lof > 10,]
kstraint[kstraint$gpiap & kstraint$lof_obs_exp >.9 & kstraint$exp_lof > 10,]


ks.test(gtex$brain_mean[gtex$description %in% gpiaps$symbol], gtex$brain_mean[!(gtex$description %in% gpiaps$symbol)], alternative='two.sided')
mean(gtex$brain_mean[gtex$description %in% gpiaps$symbol])
mean(gtex$brain_mean[!(gtex$description %in% gpiaps$symbol)])
mean(gtex$brain_mean[gtex$description %in% kstraint$gene & !(gtex$description %in% gpiaps$symbol)])

gtex$in_universe = gtex$description %in% kstraint$gene
gtex$expressed_in_brain = gtex$brain_mean > 10
gtex$gpiap = gtex$description %in% gpiaps$symbol
ctable = table(gtex[gtex$in_universe,c('expressed_in_brain','gpiap')]) # subset to kstraint$gene to match the "universe" these GPI-APs were required to be part of
fisher.test(as.matrix(ctable), alternative='two.sided')

gpiap_gtex_matrix = as.matrix(log10(gtex[gtex$gpiap,3:55]))
rownames(gpiap_gtex_matrix) = gtex$description[gtex$gpiap]
colnames(gpiap_gtex_matrix) = colnames(gtex)[3:55]
gpiap_gtex_matrix[gpiap_gtex_matrix < 0] = 0
gpiap_gtex_matrix[gpiap_gtex_matrix > 4] = 4
allzero_rows = which(rowSums(gpiap_gtex_matrix)==0)
gpiap_gtex_matrix = gpiap_gtex_matrix[-allzero_rows,]

colors = c('#fcfbfd','#efedf5','#dadaeb','#bcbddc','#9e9ac8','#807dba','#6a51a3','#4a1486')

pdf('~/d/j/cureffi/media/2018/07/gpi-gtex-heatmap-new.pdf',width=10,height=20)
ph = pheatmap(gpiap_gtex_matrix, treeheight_row=0, treeheight_col=0, col=colors, border_color = NA, labels_col=orig_colnames[3:55])
dev.off()