# cat all.htt.vcf | grep -v ^## | sed 's/^#CHROM/CHROM/' > headerless.vcf

library(reshape) # for melt
options(stringsAsFactors=FALSE)
setwd('$allele-spec-dir')  # need to modify

qual_threshold = 3000 # for me this was the lowest I could set it without including any no-call genotypes
vcf = read.table('headerless.vcf',header=TRUE) 
dim(vcf[vcf$QUAL > qual_threshold, 1:9]) # check the number of remaining SNPs after filtering

mat = as.matrix(vcf[vcf$QUAL > qual_threshold, 10:33]) # cast the value portion of the VCF to a matrix
dim(mat)
rownames(mat) = vcf$ID[vcf$QUAL > qual_threshold] # re-apply the row names

rel = melt(mat) # melt matrix to a 3-column relational table
dim(rel)
colnames(rel) = c("ID","SAMPLE","VALUE")

rel$AD = sapply(strsplit(rel$VALUE,":"),"[[",2) # parse allelic depth (AD) from value field
rel$REFAD = as.integer(sapply(strsplit(rel$AD,","),"[[",1)) # parse REF allele depth from AD
rel$ALTAD = as.integer(sapply(strsplit(rel$AD,","),"[[",2)) # parse ALT allele depth from AD
rel$ARATIO = rel$ALTAD / (rel$REFAD + rel$ALTAD) # calculate allelic ratio (% reads supporting ALT allele)

rel$BP = as.integer(vcf$POS[match(rel$ID,vcf$ID)])
plot(rel$BP, rel$ARATIO, pch=19, xlab='chr4 BP', ylab='ALT allele depth as fraction of total', main='Allelic balance at HTT SNPs')
abline(h=.5)