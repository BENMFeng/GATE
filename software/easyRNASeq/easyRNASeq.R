### R code from vignette source 'vignettes/easyRNASeq/inst/doc/easyRNASeq.Rnw'

###################################################
### code chunk number 1: Load the library and create the object
###################################################
library(easyRNASeq)
library(RnaSeqTutorial)
library(BSgenome.Dmelanogaster.UCSC.dm3)

count.table <- easyRNASeq(filesDirectory=system.file(
                            "extdata",
                            package="RnaSeqTutorial"),
                          pattern="[A,C,T,G]{6}\\.bam$",
                          readLength=30L,
                          organism="Dmelanogaster",
                          chr.sizes=seqlengths(Dmelanogaster),
                          annotationMethod="rda",
                          annotationFile=system.file(
                            "data",
                            "gAnnot.rda",
                            package="RnaSeqTutorial"),
                          count="exons"
                          )


###################################################
### code chunk number 2: Show the results
###################################################
head(count.table)
dim(count.table)


###################################################
### code chunk number 3: Create the object with only 2 files (eval = FALSE)
###################################################
## count.table <- easyRNASeq(system.file(
##                                       "extdata",
##                                       package="RnaSeqTutorial"),
##                           organism="Dmelanogaster",
##                           readLength=30L,
##                           annotationMethod="rda",
##                           annotationFile=system.file(
##                             "data",
##                             "gAnnot.rda",
##                             package="RnaSeqTutorial"),
##                           count="exons",
##                           filenames=c("ACACTG.bam", "ACTAGC.bam"))


###################################################
### code chunk number 4: Show (eval = FALSE)
###################################################
## head(count.table)
## dim(count.table)


###################################################
### code chunk number 5: read an export (eval = FALSE)
###################################################
## count.table <- easyRNASeq(system.file(
##                                       "extdata",
##                                       package="RnaSeqTutorial"),
##                           organism="Dmelanogaster",
##                           chr.sizes=seqlengths(Dmelanogaster),
##                           readLength=36L,
##                           annotationMethod="rda",
##                           annotationFile=system.file(
##                             "data",
##                             "gAnnot.rda",
##                             package="RnaSeqTutorial"),
##                           format="aln",
##                           type="SolexaExport",
##                           count="exons",
##                           pattern="subset_export",
##                           filter=compose(
##                             chromosomeFilter(regex="chr"),
##                             chastityFilter()))


###################################################
### code chunk number 6: read a gapped  bam file (eval = FALSE)
###################################################
## count.table <- easyRNASeq(system.file(
##                                       "extdata",
##                                       package="RnaSeqTutorial"),
##                           organism="Dmelanogaster",
##                           chr.sizes="auto",
##                           readLength=36L,
##                           annotationMethod="rda",
##                           annotationFile=system.file(
##                             "data",
##                             "gAnnot.rda",
##                             package="RnaSeqTutorial"),
##                           gapped=TRUE,
##                           count="exons",
##                           filenames="gapped.bam")


###################################################
### code chunk number 7: GFF3 (eval = FALSE)
###################################################
## system.file(
##             "extdata",
##             "annot.gff",
##             package="RnaSeqTutorial")


###################################################
### code chunk number 8: load annotation from biomaRt (eval = FALSE)
###################################################
## count.table <- easyRNASeq(system.file(
##                                       "extdata",
##                                       package="RnaSeqTutorial"),
##                           organism="Dmelanogaster",
##                           readLength=36L,
##                           annotationMethod="biomaRt",
##                           gapped=TRUE,
##                           count="exons",
##                           filenames="gapped.bam")


###################################################
### code chunk number 9: load the annotation from a gff (eval = FALSE)
###################################################
## count.table <- easyRNASeq(system.file(
##                                       "extdata",
##                                       package="RnaSeqTutorial"),
##                           organism="Dmelanogaster",
##                           readLength=36L,
##                           annotationMethod="gff",
##                           annotationFile=system.file(
##                             "extdata",
##                             "annot.gff",
##                             package="RnaSeqTutorial"),
##                           gapped=TRUE,
##                           count="exons",
##                           filenames="gapped.bam")


###################################################
### code chunk number 10: RangedData structure
###################################################
gAnnot <- RangedData(
                     IRanges(
                             start=c(10,30,100),
                             end=c(21,53,123)),
                          space=c("chr01","chr01","chr02"),
                          strand=c("+","+","-"),
                          transcript=c("trA1","trA2","trB"),
                          gene=c("gA","gA","gB"),
                          exon=c("e1","e2","e3"),
                          universe = "Hs19"
                          )
gAnnot


###################################################
### code chunk number 11: look at the gAnnot (eval = FALSE)
###################################################
## data(gAnnot)
## gAnnot


###################################################
### code chunk number 12: GRangesList coercion
###################################################
grngs <- as(gAnnot,"GRanges")
grngsList<-split(grngs,seqnames(grngs))
grngsList


###################################################
### code chunk number 13: easyRNASeq output
###################################################
rnaSeq <- easyRNASeq(system.file(
                                 "extdata",
                                 package="RnaSeqTutorial"),
                     organism="Dmelanogaster",
                     readLength=30L,
                     annotationMethod="rda",
                     annotationFile=system.file(
                       "data",
                       "gAnnot.rda",
                       package="RnaSeqTutorial"),
                     count="exons",
                     pattern="[A,C,T,G]{6}\\.bam$",
                     outputFormat="RNAseq")


###################################################
### code chunk number 14: Transcript counts
###################################################
rnaSeq <- transcriptCounts(rnaSeq)
head(readCounts(rnaSeq,'transcripts'))


###################################################
### code chunk number 15: geneCounts fails (eval = FALSE)
###################################################
## try(geneCounts(rnaSeq,summarization='geneModels'))


###################################################
### code chunk number 16: current access to geneCounts (eval = FALSE)
###################################################
## rnaSeq2 <- easyRNASeq(system.file(
##                                   "extdata",
##                                   package="RnaSeqTutorial"),
##                       organism="Dmelanogaster",
##                       readLength=30L,
##                       annotationMethod="rda",
##                       annotationFile=system.file(
##                         "data",
##                         "gAnnot.rda",
##                         package="RnaSeqTutorial"),
##                       count="genes",
##                       summarization="geneModels",
##                       pattern="[A,C,T,G]{6}\\.bam$",
##                       outputFormat="RNAseq")
## head(readCounts(rnaSeq2,'genes','geneModels'))


###################################################
### code chunk number 17: RPKM normalization (eval = FALSE)
###################################################
## count.table <- easyRNASeq(system.file(
##                                       "extdata",
##                                       package="RnaSeqTutorial"),
##                           organism="Dmelanogaster",
##                           readLength=30L,
##                           annotationMethod="rda",
##                           annotationFile=system.file(
##                             "data",
##                             "gAnnot.rda",
##                             package="RnaSeqTutorial"),
##                           count="exons",
##                           filenames=c("ACACTG.bam", "ACTAGC.bam", 
##                             "ATGGCT.bam", "TTGCGA.bam"),
##                           normalize=TRUE
##                           )


###################################################
### code chunk number 18: RPKM counts
###################################################
feature.size = width(genomicAnnotation(rnaSeq))
names(feature.size) = genomicAnnotation(rnaSeq)$exon
lib.size=c(
  "ACACTG.bam"=56643,
  "ACTAGC.bam"=42698,
  "ATGGCT.bam"=55414,
  "TTGCGA.bam"=60740)
head(RPKM(readCounts(rnaSeq,summarization="exons")$exons,
          NULL,
          lib.size=lib.size,
          feature.size=feature.size))


###################################################
### code chunk number 19: RPKM on RNAseq
###################################################
head(RPKM(rnaSeq,from="transcripts"))


###################################################
### code chunk number 20: RPKM on RNAseq (eval = FALSE)
###################################################
## head(RPKM(rnaSeq2,from="geneModels"))


###################################################
### code chunk number 21: define the conditions
###################################################
conditions=c("A","A","B","B")
names(conditions) <- c("ACACTG.bam", "ACTAGC.bam", 
                       "ATGGCT.bam", "TTGCGA.bam")


###################################################
### code chunk number 22: DESeq normalization (eval = FALSE)
###################################################
## countDataSet <- easyRNASeq(system.file(
##                                       "extdata",
##                                       package="RnaSeqTutorial"),
##                           organism="Dmelanogaster",
##                           readLength=30L,
##                           annotationMethod="rda",
##                           annotationFile=system.file(
##                             "data",
##                             "gAnnot.rda",
##                             package="RnaSeqTutorial"),
##                           count="exons",
##                           filenames=c("ACACTG.bam", "ACTAGC.bam", 
##                             "ATGGCT.bam", "TTGCGA.bam"),
##                           normalize=TRUE,
##                           outputFormat="DESeq",
##                           conditions=conditions,
##                           fitType="local"
##                           )


###################################################
### code chunk number 23: edgeR normalization (eval = FALSE)
###################################################
## dgeList <- easyRNASeq(system.file(
##                                   "extdata",
##                                   package="RnaSeqTutorial"),
##                       organism="Dmelanogaster",
##                       readLength=30L,
##                       annotationMethod="rda",
##                       annotationFile=system.file(
##                         "data",
##                         "gAnnot.rda",
##                         package="RnaSeqTutorial"),
##                       count="exons",
##                       filenames=c("ACACTG.bam", "ACTAGC.bam",
##                         "ATGGCT.bam", "TTGCGA.bam"),
##                       normalize=TRUE,
##                       outputFormat="edgeR",
##                       conditions=conditions
##                       )


###################################################
### code chunk number 24: cleanup
###################################################
file2clean=c(
  "ACACTG.fastq",
  "ACTAGC.fastq",
  "ATGGCT.fastq",
  "TTGCGA.fastq")
silent <- sapply(
                 file2clean,
                 function(f2c){
                   if(file.exists(f2c)){file.remove(f2c)}
                 })


###################################################
### code chunk number 25: Demultiplexing sample
###################################################
alns <- readAligned(
                    system.file(
                                "extdata",
                                package="RnaSeqTutorial"),
                    pattern="multiplex_export",
                    filter=compose(
                      chastityFilter(),
                      nFilter(2),
                      chromosomeFilter(regex="chr")),
                    type="SolexaExport",
                    withAll=TRUE)


###################################################
### code chunk number 26: Illumina example (eval = FALSE)
###################################################
## alignData(alns)$multiplexIndex


###################################################
### code chunk number 27: barcodes
###################################################
barcodes=c("ACACTG","ACTAGC","ATGGCT","TTGCGA")


###################################################
### code chunk number 28: barcodePlot
###################################################
barcodePlot(alns,
            barcodes=barcodes,
            type="within",
            barcode.length=6,
            show.barcode=20,
            main="All samples",
            xlim=c(0,0.5))


###################################################
### code chunk number 29: demultiplex
###################################################
dem.alns <- demultiplex(alns,
                        barcodes=barcodes,
                        edition.dist=2,
                        barcodes.qty=4,
                        type="within")


###################################################
### code chunk number 30: check the objects (eval = FALSE)
###################################################
## dem.alns$reads[[1]]
## dem.alns$barcodes[[1]]


###################################################
### code chunk number 31: demPlot
###################################################
par(mfrow=c(2,2))
barcode.frequencies <- lapply(
                              names(dem.alns$barcodes),
                              function(barcode,alns){                                
                                barcodePlot(
                                            alns$barcodes[[barcode]],
                                            barcodes=barcode,
                                            type="within",barcode.length=6,
                                            show.barcode=20,
                                            main=paste(
                                              "Expected barcode:",
                                              barcode))
                              },dem.alns)


###################################################
### code chunk number 32: write the file
###################################################
status <- lapply(
                 names(dem.alns$barcodes),
                 function(barcode,alns){         
                   writeFastq(
                              alns$reads[[barcode]],
                              file=paste(
                                barcode,
                                "fastq",sep="."))
                 },dem.alns)


###################################################
### code chunk number 33: cleanup
###################################################
file2clean=c(
  "ACACTG.fastq",
  "ACTAGC.fastq",
  "ATGGCT.fastq",
  "TTGCGA.fastq",
  "Rplots.pdf")
silent <- sapply(
                 file2clean,
                 function(f2c){
                   if(file.exists(f2c)){file.remove(f2c)}
                 })


###################################################
### code chunk number 34: count function example (eval = FALSE)
###################################################
## ## creating a SummarizedExperiment from 4 bam files
## sumExp <- count(filesDirectory=system.file(
##                   "extdata",
##                   package="RnaSeqTutorial"),
##                 pattern="[A,C,T,G]{6}\\.bam$",
##                 readLength=30L,
##                 organism="Dmelanogaster",
##                 chr.sizes=seqlengths(Dmelanogaster),
##                 annotationMethod="rda",
##                 annotationFile=system.file(
##                   "data",
##                   "gAnnot.rda",
##                   package="RnaSeqTutorial"),
##                 count="exons"
##                 )
## ## the counts
## assays(sumExp)
## ## the sample info
## colData(sumExp)
## ## the 'features' info
## rowData(sumExp)
## 


###################################################
### code chunk number 35: load the lib (eval = FALSE)
###################################################
## library(easyRNASeq)


###################################################
### code chunk number 36: bam prep (eval = FALSE)
###################################################
## library(BSgenome.Hsapiens.UCSC.hg19)
## library(GEOquery)
## library(SRAdb)
## library(Rsamtools)
## library(Rsubread)
## 
## ## create a temp dir
## dir.create(file.path(getwd(),"tmp1234"))
## 
## ## change the working directory
## setwd(file.path(getwd(),"tmp1234"))
## 
## ## init SRA
## sqlfile <- getSRAdbFile()
## 
## ## init a connection
## sra_con <- dbConnect(SQLite(),sqlfile)
## 
## ## list the files
## acc <- c("SRR490224","SRR490225")
## getFASTQinfo( in_acc = acc, srcType = 'ftp' )
## 
## ## get the read files
## getSRAfile( in_acc=acc, sra_con, destDir = getwd(),
##            fileType = 'fastq', srcType = 'ftp' )
## 
## ## close the connection
## dbDisconnect( sra_con )
## 
## ## write the human genome sequences
## writeXStringSet(Reduce(append,
##                        lapply(seqnames(Hsapiens),
##                               function(nam){
##                                 dss <- DNAStringSet(unmasked(Hsapiens[[nam]]))
##                                 names(dss) <- nam
##                                 dss
##                               })),
##                 file="hg19.fa")
## 
## ## create the indexes
## dir.create("indexes")
## buildindex(basename=file.path("indexes","hg19"),
##            reference="hg19.fa")
## 
## ## align the reads
## sapply(dir(pattern="*\\.gz$"),function(fil){
##   ## decompress the files
##   gunzip(fil)
##   
##   ## align
##   align(index=file.path("indexes","hg19"),
##         readfile1=sub("\\.gz$","",fil),
##         nsubreads=2,TH1=1,
##         output_file=sub("\\.fastq\\.gz$","\\.sam",fil)
##         )
##   
##   ## create bam files
##   asBam(file=sub("\\.fastq\\.gz$","\\.sam",fil),
##         destination=sub("\\.fastq\\.gz$","",fil),
##         indexDestination=TRUE)
## })
## 


###################################################
### code chunk number 37: chromosome size (eval = FALSE)
###################################################
## library(BSgenome.Hsapiens.UCSC.hg19)
## chr.sizes=seqlengths(Hsapiens)


###################################################
### code chunk number 38: bam files name (eval = FALSE)
###################################################
## bamfiles=dir(getwd(),pattern="*\\.bam$")


###################################################
### code chunk number 39: Fetching with biomaRt (eval = FALSE)
###################################################
## rnaSeq <- easyRNASeq(filesDirectory=getwd(),
##                     organism="Hsapiens",
##                     chr.sizes=chr.sizes,
##                     readLength=58L,
##                     annotationMethod="biomaRt",
##                     count="exons",
##                     filenames=bamfiles[1],
##                     outputFormat="RNAseq"
##                     )
## gAnnot <- genomicAnnotation(rnaSeq)


###################################################
### code chunk number 40: filtering the annot (eval = FALSE)
###################################################
## gAnnot <- gAnnot[space(gAnnot) %in% paste("chr",c(1:22,"X","Y","M"),sep=""),]
## save(gAnnot,file="gAnnot.rda")


###################################################
### code chunk number 41: get a count table (eval = FALSE)
###################################################
## countTable <- easyRNASeq(filesDirectory=getwd(),
##                         organism="Hsapiens",
##                         chr.sizes=chr.sizes,
##                         readLength=58L,
##                         annotationMethod="rda",
##                         annotationFile="gAnnot.rda",
##                         count="exons",
##                         filenames=bamfiles[1]
##                         )


###################################################
### code chunk number 42: using different annotation (eval = FALSE)
###################################################
## library(GenomicFeatures)
## hg19.tx <- makeTranscriptDbFromUCSC(
##                                    genome="hg19",
##                                    tablename="refGene")
## 
## gAnnot <- exons(hg19.tx)
## 
## colnames(elementMetadata(gAnnot)) <- "exon"
## 
## gAnnot <- split(gAnnot,seqnames(gAnnot))


###################################################
### code chunk number 43: sub setting (eval = FALSE)
###################################################
## countTable <- easyRNASeq(filesDirectory=getwd(),
##                         organism="Hsapiens",
##                         readLength=58L,
##                         annotationMethod="env",
##                         annotationObject=gAnnot,
##                         count="exons",
##                         filenames=bamfiles[1],
##                         chr.sel="chr1"
##                         )


###################################################
### code chunk number 44: using the count function (eval = FALSE)
###################################################
## sumExp <- count(filesDirectory=getwd(),
##                     organism="Hsapiens",
##                     readLength=58L,
##                     annotationMethod="env",
##                     annotationObject=gAnnot,
##                     count="exons",
##                     filenames=bamfiles[1],
##                     chr.sel="chr1"
##                     )


###################################################
### code chunk number 45: DESeq run (eval = FALSE)
###################################################
## conditions <- c("A","B")
## names(conditions) <- bamfiles
## 
## countDataSet <- easyRNASeq(filesDirectory=getwd(),
##                           organism="Hsapiens",
##                           chr.sizes=chr.sizes,
##                           annotationMethod="env",
##                           annotationObject=gAnnot,
##                           count="exons",
##                           filenames=bamfiles,
##                           chr.sel="chr1",
##                           outputFormat="DESeq",
##                           normalize=TRUE,
##                           conditions=conditions,
##                           fitType="local",
##                           method="blind"
##                         )


###################################################
### code chunk number 46: read the data (eval = FALSE)
###################################################
## aln <- readAligned("data",type="SolexaExport",pattern="*.txt.gz")
## gc()
## levels(chromosome(aln))


###################################################
### code chunk number 47: fake the data
###################################################
c("chr1.fa","chr10.fa","...","chrY.fa")


###################################################
### code chunk number 48: fetch the annotation (eval = FALSE)
###################################################
## obj <- fetchAnnotation(new('RNAseq',
##                           organismName="Mmusculus"
##                           ),
##                       method="biomaRt")
## 
## gAnnot <- genomicAnnotation(obj)
## 
## length(grep("NT_",space(gAnnot)))
## 


###################################################
### code chunk number 49: fake number 2
###################################################
1181


###################################################
### code chunk number 50: modify the annotation (eval = FALSE)
###################################################
## names(gAnnot) <- paste("chr",names(gAnnot),".fa",sep="")


###################################################
### code chunk number 51: clean (eval = FALSE)
###################################################
## rm(aln,obj)
## gc()


###################################################
### code chunk number 52: chrsize (eval = FALSE)
###################################################
## library(BSgenome.Mmusculus.UCSC.mm9)
## chr.sizes<- seqlengths(Mmusculus)


###################################################
### code chunk number 53: create the chr.map (eval = FALSE)
###################################################
## chr.map <- data.frame(
##                       from=paste("chr",c(1:19,"X","Y"),".fa",sep=""),
##                       to=paste("chr",c(1:19,"X","Y"),sep=""))


###################################################
### code chunk number 54: no filter (eval = FALSE)
###################################################
## rnaSeq <- easyRNASeq(filesDirectory="data",
##                     organism="custom",
##                     chr.map=chr.map,
##                     chr.sizes=chr.sizes,
##                     filter=compose(
##                       naPositionFilter(),
##                       chastityFilter()),
##                     readLength=50L,
##                     annotationMethod="env",
##                     annotationObject=gAnnot,
##                     format="aln",
##                     count="genes",
##                     summarization= "geneModels",
##                     filenames="1-Feb_ATCACG_L003_R1_001_export.txt.gz",
##                     outputFormat="RNAseq",
##                     nbCore=2
##                     )


###################################################
### code chunk number 55: selected chr (eval = FALSE)
###################################################
## rnaSeq <- easyRNASeq(filesDirectory="data",
##                     organism="custom",
##                     chr.map=chr.map,
##                     chr.sizes=chr.sizes,
##                     chr.sel=chr.map$from,
##                     filter=compose(
##                       naPositionFilter(),
##                       chastityFilter()),
##                     readLength=50L,
##                     annotationMethod="env",
##                     annotationObject=gAnnot,
##                     format="aln",
##                     count="genes",
##                     summarization= "geneModels",
##                     filenames="1-Feb_ATCACG_L003_R1_001_export.txt.gz",
##                     outputFormat="RNAseq",
##                     nbCore=2
##                     )


###################################################
### code chunk number 56: remove annotation (eval = FALSE)
###################################################
## sel <- grep("NT_",names(gAnnot))
## gAnnot <- RangedData(ranges=ranges(gAnnot)[-sel,],values=values(gAnnot)[-sel,])
## colnames(gAnnot) <- gsub("values\\.","",colnames(gAnnot))


###################################################
### code chunk number 57: final call (eval = FALSE)
###################################################
## rnaSeq <- easyRNASeq(filesDirectory="data",
##                     organism="custom",
##                     chr.map=chr.map,
##                     chr.sizes=chr.sizes,
##                     chr.sel=chr.map$from,
##                     filter=compose(
##                       naPositionFilter(),
##                       chastityFilter()),
##                     readLength=50L,
##                     annotationMethod="env",
##                     annotationObject=gAnnot,
##                     format="aln",
##                     count="genes",
##                     summarization= "geneModels",
##                     filenames="1-Feb_ATCACG_L003_R1_001_export.txt.gz",
##                     outputFormat="RNAseq",
##                     nbCore=2
##                     )


###################################################
### code chunk number 58: SessionInfo
###################################################
library(easyRNASeq)
library(BSgenome.Dmelanogaster.UCSC.dm3)
library(RnaSeqTutorial)
sessionInfo()


###################################################
### code chunk number 59: some test
###################################################
grngs = GRanges(seqnames=c("chr1", "chr2"),
  ranges=IRanges(start=1:2, end=2:3), 
  strand=c("+","-"))
silent <- strand(grngs)
silent <- reduce(grngs)


