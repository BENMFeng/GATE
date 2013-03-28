### R code from vignette source 'vignettes/cummeRbund/inst/doc/cummeRbund-example-workflow.Rnw'

###################################################
### code chunk number 1: init
###################################################
options(width=65)


###################################################
### code chunk number 2: loadLib
###################################################
library(cummeRbund)


###################################################
### code chunk number 3: read
###################################################
cuff <- readCufflinks(dir=system.file("extdata", package="cummeRbund"))
cuff


###################################################
### code chunk number 4: model_fit_1
###################################################
d<-dispersionPlot(genes(cuff))
d


###################################################
### code chunk number 5: model_fit_1_plot
###################################################
d<-dispersionPlot(genes(cuff))
d
print(d)


###################################################
### code chunk number 6: rep_boxplot_1
###################################################
pBoxRep<-csBoxplot(genes(cuff),replicates=T)
pBoxRep


###################################################
### code chunk number 7: rep_dendro_1
###################################################
pDendro<-csDendro(genes(cuff),replicates=T)
pDendro


###################################################
### code chunk number 8: rep_boxplot_1_plot
###################################################
pBoxRep<-csBoxplot(genes(cuff),replicates=T)
pBoxRep
print(pBoxRep)


###################################################
### code chunk number 9: rep_dendro_1_plot
###################################################
pDendro<-csDendro(genes(cuff),replicates=T)
pDendro
print(pDendro)


###################################################
### code chunk number 10: boxplot_1
###################################################
pBox<-csBoxplot(genes(cuff))
pBox


###################################################
### code chunk number 11: boxplot_1_plot
###################################################
pBox<-csBoxplot(genes(cuff))
pBox
print(pBox)


###################################################
### code chunk number 12: diff_exp_genes_1
###################################################
sigGeneIds<-getSig(cuff,alpha=0.05,level="genes")
head(sigGeneIds)
length(sigGeneIds)


###################################################
### code chunk number 13: diff_exp_genes_2
###################################################
hESCvsFibroblast.sigGeneIds<-getSig(cuff,"hESC","Fibroblasts",alpha=0.05,level="genes")
head(hESCvsFibroblast.sigGeneIds)
length(hESCvsFibroblast.sigGeneIds)


###################################################
### code chunk number 14: diff_exp_genes_3
###################################################
sigGenes<-getGenes(cuff,sigGeneIds)
sigGenes


###################################################
### code chunk number 15: diff_exp_feat_1
###################################################
sigGeneIds<-getSig(cuff,alpha=0.05,level="isoforms")
head(sigGeneIds)
length(sigGeneIds)


###################################################
### code chunk number 16: ind_gene_1
###################################################
Pink1<-getGene(cuff,'PINK1')
Pink1


