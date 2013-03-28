### R code from vignette source 'vignettes/cummeRbund/inst/doc/cummeRbund-manual.Rnw'

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
myDir<-system.file("extdata", package="cummeRbund") #You can leave blank if cwd or replace with your own directory path.
gtfFile<-system.file("extdata/chr1_snippet.gtf",package="cummeRbund") #path to .gtf file used in cuffdiff analysis.
cuff <- readCufflinks(dir=myDir,gtfFile=gtfFile,genome="hg19",rebuild=T)


###################################################
### code chunk number 4: read2 (eval = FALSE)
###################################################
## cuff<-readCufflinks()


###################################################
### code chunk number 5: read3
###################################################
cuff


###################################################
### code chunk number 6: add_features
###################################################
#annot<-read.table("gene_annotation.tab",sep="\t",header=T,na.string="-")
#addFeatures(cuff,annot,level="genes")


###################################################
### code chunk number 7: global_dispersion
###################################################
disp<-dispersionPlot(genes(cuff))
disp


###################################################
### code chunk number 8: global_dispersion_plot
###################################################
	print(disp)


###################################################
### code chunk number 9: SCV_visualization
###################################################
genes.scv<-fpkmSCVPlot(genes(cuff))
isoforms.scv<-fpkmSCVPlot(isoforms(cuff))


###################################################
### code chunk number 10: global_plots_1
###################################################
dens<-csDensity(genes(cuff))
dens
densRep<-csDensity(genes(cuff),replicates=T)
densRep


###################################################
### code chunk number 11: global_plots_dens
###################################################
dens<-csDensity(genes(cuff))
dens
densRep<-csDensity(genes(cuff),replicates=T)
densRep
	print(dens)


###################################################
### code chunk number 12: global_plots_dens_rep
###################################################
	print(densRep)


###################################################
### code chunk number 13: global_plots_2
###################################################
b<-csBoxplot(genes(cuff))
b
brep<-csBoxplot(genes(cuff),replicates=T)
brep


###################################################
### code chunk number 14: global_plots_box
###################################################
b<-csBoxplot(genes(cuff))
b
brep<-csBoxplot(genes(cuff),replicates=T)
brep
	print(b)


###################################################
### code chunk number 15: global_plots_box_rep
###################################################
	print(brep)


###################################################
### code chunk number 16: global_plots_3.1
###################################################
s<-csScatterMatrix(genes(cuff))



###################################################
### code chunk number 17: global_plots_scatter_1
###################################################
s<-csScatterMatrix(genes(cuff))

	print(s)


###################################################
### code chunk number 18: global_plots_3.2
###################################################
s<-csScatter(genes(cuff),"hESC","Fibroblasts",smooth=T)
s


###################################################
### code chunk number 19: global_plots_scatter_2
###################################################
s<-csScatter(genes(cuff),"hESC","Fibroblasts",smooth=T)
s
	print(s)


###################################################
### code chunk number 20: global_plots_6
###################################################
dend<-csDendro(genes(cuff))
dend.rep<-csDendro(genes(cuff),replicates=T)


###################################################
### code chunk number 21: global_plots_dendro
###################################################
dend<-csDendro(genes(cuff))
dend.rep<-csDendro(genes(cuff),replicates=T)
	plot(dend)


###################################################
### code chunk number 22: global_plots_dendro_rep
###################################################
	plot(dend.rep)


###################################################
### code chunk number 23: global_plots_4
###################################################
m<-MAplot(genes(cuff),"hESC","Fibroblasts")
m
mCount<-MAplot(genes(cuff),"hESC","Fibroblasts",useCount=T)
mCount


###################################################
### code chunk number 24: global_plots_MA
###################################################
m<-MAplot(genes(cuff),"hESC","Fibroblasts")
m
mCount<-MAplot(genes(cuff),"hESC","Fibroblasts",useCount=T)
mCount
	print(m)


###################################################
### code chunk number 25: global_plots_MA_count
###################################################
	print(mCount)


###################################################
### code chunk number 26: global_plots_5_1
###################################################
v<-csVolcanoMatrix(genes(cuff))
v


###################################################
### code chunk number 27: global_plots_volcano_1
###################################################
v<-csVolcanoMatrix(genes(cuff))
v
print(v)


###################################################
### code chunk number 28: global_plots_5_2
###################################################
v<-csVolcano(genes(cuff),"hESC","Fibroblasts")
v


###################################################
### code chunk number 29: global_plots_volcano_2
###################################################
v<-csVolcano(genes(cuff),"hESC","Fibroblasts")
v
print(v)


###################################################
### code chunk number 30: data_access_0
###################################################
runInfo(cuff)
replicates(cuff)


###################################################
### code chunk number 31: data_access_1
###################################################
gene.features<-annotation(genes(cuff))
head(gene.features)

gene.fpkm<-fpkm(genes(cuff))
head(gene.fpkm)

gene.repFpkm<-repFpkm(genes(cuff))
head(gene.repFpkm)

gene.counts<-count(genes(cuff))
head(gene.counts)

isoform.fpkm<-fpkm(isoforms(cuff))
head(isoform.fpkm)

gene.diff<-diffData(genes(cuff))
head(gene.diff)


###################################################
### code chunk number 32: data_access_2
###################################################
sample.names<-samples(genes(cuff))
head(sample.names)
gene.featurenames<-featureNames(genes(cuff))
head(gene.featurenames)


###################################################
### code chunk number 33: data_access_3
###################################################
gene.matrix<-fpkmMatrix(genes(cuff))
head(gene.matrix)


###################################################
### code chunk number 34: data_access_4
###################################################
gene.rep.matrix<-repFpkmMatrix(genes(cuff))
head(gene.rep.matrix)


###################################################
### code chunk number 35: data_access_5
###################################################
gene.count.matrix<-countMatrix(genes(cuff))
head(gene.count.matrix)


###################################################
### code chunk number 36: create_geneset_1
###################################################
data(sampleData)
myGeneIds<-sampleIDs
myGeneIds
myGenes<-getGenes(cuff,myGeneIds)
myGenes


###################################################
### code chunk number 37: create_geneset_2
###################################################
#FPKM values for genes in gene set
head(fpkm(myGenes))

#Isoform-level FPKMs for gene set
head(fpkm(isoforms(myGenes)))

#Replicate FPKMs for TSS groups within gene set
head(repFpkm(TSS(myGenes)))



###################################################
### code chunk number 38: geneset_plots_1
###################################################
h<-csHeatmap(myGenes,cluster='both')
h
h.rep<-csHeatmap(myGenes,cluster='both',replicates=T)
h.rep


###################################################
### code chunk number 39: geneset_plots_heatmap
###################################################
h<-csHeatmap(myGenes,cluster='both')
h
h.rep<-csHeatmap(myGenes,cluster='both',replicates=T)
h.rep
print(h)


###################################################
### code chunk number 40: geneset_plots_heatmap_rep
###################################################
print(h.rep)


###################################################
### code chunk number 41: geneset_plots_1.5
###################################################
b<-expressionBarplot(myGenes)
b


###################################################
### code chunk number 42: geneset_plots_barplot
###################################################
b<-expressionBarplot(myGenes)
b
print(b)


###################################################
### code chunk number 43: geneset_plots_2
###################################################
s<-csScatter(myGenes,"Fibroblasts","hESC",smooth=T)
s


###################################################
### code chunk number 44: geneset_plots_scatter
###################################################
s<-csScatter(myGenes,"Fibroblasts","hESC",smooth=T)
s
print(s)


###################################################
### code chunk number 45: geneset_plots_3
###################################################
v<-csVolcano(myGenes,"Fibroblasts","hESC")
v


###################################################
### code chunk number 46: geneset_plots_volcano
###################################################
v<-csVolcano(myGenes,"Fibroblasts","hESC")
v
print(v)


###################################################
### code chunk number 47: geneset_plots_4
###################################################
ih<-csHeatmap(isoforms(myGenes),cluster='both',labRow=F)
ih
th<-csHeatmap(TSS(myGenes),cluster='both',labRow=F)
th


###################################################
### code chunk number 48: geneset_plots_isoform_heatmap
###################################################
ih<-csHeatmap(isoforms(myGenes),cluster='both',labRow=F)
ih
th<-csHeatmap(TSS(myGenes),cluster='both',labRow=F)
th
print(ih)


###################################################
### code chunk number 49: geneset_plots_TSS_heatmap
###################################################
print(th)


###################################################
### code chunk number 50: geneset_plots_5
###################################################
den<-csDendro(myGenes)


###################################################
### code chunk number 51: geneset_plots_dendro
###################################################
den<-csDendro(myGenes)
plot(den)


###################################################
### code chunk number 52: gene_level_1
###################################################
myGeneId<-"PINK1"
myGene<-getGene(cuff,myGeneId)
myGene
head(fpkm(myGene))
head(fpkm(isoforms(myGene)))


###################################################
### code chunk number 53: gene_plots_1
###################################################
gl<-expressionPlot(myGene)
gl

gl.rep<-expressionPlot(myGene,replicates=TRUE)
gl.rep

gl.iso.rep<-expressionPlot(isoforms(myGene),replicates=T)
gl.iso.rep

gl.cds.rep<-expressionPlot(CDS(myGene),replicates=T)
gl.cds.rep


###################################################
### code chunk number 54: gene_plots_line
###################################################
gl<-expressionPlot(myGene)
gl

gl.rep<-expressionPlot(myGene,replicates=TRUE)
gl.rep

gl.iso.rep<-expressionPlot(isoforms(myGene),replicates=T)
gl.iso.rep

gl.cds.rep<-expressionPlot(CDS(myGene),replicates=T)
gl.cds.rep
	print(gl)


###################################################
### code chunk number 55: gene_plots_replicate_line
###################################################

	print(gl.rep)


###################################################
### code chunk number 56: gene_plots_iso_replicate_line
###################################################

	print(gl.iso.rep)


###################################################
### code chunk number 57: gene_plots_cds_replicate_line
###################################################

	print(gl.cds.rep)


###################################################
### code chunk number 58: gene_plots_2
###################################################
gb<-expressionBarplot(myGene)
gb
gb.rep<-expressionBarplot(myGene,replicates=T)
gb.rep


###################################################
### code chunk number 59: gene_plots_bar
###################################################
gb<-expressionBarplot(myGene)
gb
gb.rep<-expressionBarplot(myGene,replicates=T)
gb.rep
print(gb)


###################################################
### code chunk number 60: gene_plots_bar_rep
###################################################
print(gb.rep)


###################################################
### code chunk number 61: gene_plots_3
###################################################
igb<-expressionBarplot(isoforms(myGene),replicates=T)
igb


###################################################
### code chunk number 62: gene_plots_bar_isoforms
###################################################
igb<-expressionBarplot(isoforms(myGene),replicates=T)
igb
print(igb)


###################################################
### code chunk number 63: features_1
###################################################
head(features(myGene))


###################################################
### code chunk number 64: features_2
###################################################
genetrack<-makeGeneRegionTrack(myGene)
plotTracks(genetrack)


###################################################
### code chunk number 65: features_3
###################################################
trackList<-list()
myStart<-min(features(myGene)$start)
myEnd<-max(features(myGene)$end)
myChr<-unique(features(myGene)$seqnames)
genome<-'hg19'

ideoTrack <- IdeogramTrack(genome = genome, chromosome = myChr)
trackList<-c(trackList,ideoTrack)

axtrack<-GenomeAxisTrack()
trackList<-c(trackList,axtrack)
 
genetrack<-makeGeneRegionTrack(myGene)
genetrack
 
trackList<-c(trackList,genetrack)
 
biomTrack<-BiomartGeneRegionTrack(genome=genome,chromosome=as.character(myChr),
		start=myStart,end=myEnd,name="ENSEMBL",showId=T)
 
trackList<-c(trackList,biomTrack)
 
conservation <- UcscTrack(genome = genome, chromosome = myChr,
		track = "Conservation", table = "phyloP46wayPlacental",
		from = myStart-2000, to = myEnd+2000, trackType = "DataTrack",
		start = "start", end = "end", data = "score",
		type = "hist", window = "auto", col.histogram = "darkblue",
		fill.histogram = "darkblue", ylim = c(-3.7, 4),
		name = "Conservation")

trackList<-c(trackList,conservation)
 
plotTracks(trackList,from=myStart-2000,to=myEnd+2000)



###################################################
### code chunk number 66: sig_mat_1
###################################################
mySigMat<-sigMatrix(cuff,level='genes',alpha=0.05)



###################################################
### code chunk number 67: sig_mat_plot_1
###################################################
mySigMat<-sigMatrix(cuff,level='genes',alpha=0.05)

print(mySigMat)


###################################################
### code chunk number 68: get_sig_1
###################################################
mySigGeneIds<-getSig(cuff,alpha=0.05,level='genes')
head(mySigGeneIds)
length(mySigGeneIds)


###################################################
### code chunk number 69: get_sig_2
###################################################
hESC_vs_iPS.sigIsoformIds<-getSig(cuff,x='hESC',y='iPS',alpha=0.05,level='isoforms')
head(hESC_vs_iPS.sigIsoformIds)
length(hESC_vs_iPS.sigIsoformIds)


###################################################
### code chunk number 70: get_sig_3
###################################################
mySigGenes<-getGenes(cuff,mySigGeneIds)
mySigGenes



###################################################
### code chunk number 71: get_sig_4
###################################################
mySigTable<-getSigTable(cuff,alpha=0.01,level='genes')
head(mySigTable,20)


###################################################
### code chunk number 72: dist_heat_1
###################################################
myDistHeat<-csDistHeat(genes(cuff))



###################################################
### code chunk number 73: dist_heat_plot_1
###################################################
myDistHeat<-csDistHeat(genes(cuff))

print(myDistHeat)


###################################################
### code chunk number 74: dist_heat_2
###################################################
myRepDistHeat<-csDistHeat(genes(cuff),replicates=T)



###################################################
### code chunk number 75: dist_heat_plot_2
###################################################
myRepDistHeat<-csDistHeat(genes(cuff),replicates=T)

print(myRepDistHeat)


###################################################
### code chunk number 76: dim_reduction_1
###################################################
genes.PCA<-PCAplot(genes(cuff),"PC1","PC2")
genes.MDS<-MDSplot(genes(cuff))

genes.PCA.rep<-PCAplot(genes(cuff),"PC1","PC2",replicates=T)
genes.MDS.rep<-MDSplot(genes(cuff),replicates=T)


###################################################
### code chunk number 77: gene_PCA
###################################################
genes.PCA<-PCAplot(genes(cuff),"PC1","PC2")
genes.MDS<-MDSplot(genes(cuff))

genes.PCA.rep<-PCAplot(genes(cuff),"PC1","PC2",replicates=T)
genes.MDS.rep<-MDSplot(genes(cuff),replicates=T)
	print(genes.PCA)


###################################################
### code chunk number 78: gene_MDS
###################################################

	print(genes.MDS)


###################################################
### code chunk number 79: gene_PCA_rep
###################################################

	print(genes.PCA.rep)


###################################################
### code chunk number 80: gene_MDS_rep
###################################################

	print(genes.MDS.rep)


###################################################
### code chunk number 81: geneset_cluster_1
###################################################
ic<-csCluster(myGenes,k=4)
head(ic$cluster)
icp<-csClusterPlot(ic)
icp


###################################################
### code chunk number 82: geneset_plots_cluster
###################################################
print(icp)


###################################################
### code chunk number 83: specificity_1
###################################################
myGenes.spec<-csSpecificity(myGenes)
head(myGenes.spec)


###################################################
### code chunk number 84: similar_1
###################################################
mySimilar<-findSimilar(cuff,"PINK1",n=20)
mySimilar.expression<-expressionPlot(mySimilar,logMode=T,showErrorbars=F)


###################################################
### code chunk number 85: similar_plots_1
###################################################
mySimilar<-findSimilar(cuff,"PINK1",n=20)
mySimilar.expression<-expressionPlot(mySimilar,logMode=T,showErrorbars=F)
print(mySimilar.expression)


###################################################
### code chunk number 86: similar_2
###################################################
myProfile<-c(500,0,400)
mySimilar2<-findSimilar(cuff,myProfile,n=10)
mySimilar2.expression<-expressionPlot(mySimilar2,logMode=T,showErrorbars=F)


###################################################
### code chunk number 87: similar_plots_2
###################################################
myProfile<-c(500,0,400)
mySimilar2<-findSimilar(cuff,myProfile,n=10)
mySimilar2.expression<-expressionPlot(mySimilar2,logMode=T,showErrorbars=F)
print(mySimilar2.expression)


###################################################
### code chunk number 88: close_connection
###################################################
end<-sqliteCloseConnection(cuff@DB)


###################################################
### code chunk number 89: session
###################################################
sessionInfo()


