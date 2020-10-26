require(DESeq2)


#define function for geometric mean
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


#counts table geometric mean transformation
phyloseq_gm_mean_trans <- function(physeq){
  
Dor_ps.dds <- phyloseq_to_deseq2(physeq, ~1)
geoMeans = apply(counts(Dor_ps.dds), 1, gm_mean)
Dor_ps.dds = estimateSizeFactors(Dor_ps.dds, geoMeans = geoMeans)
Dor_ps.dds <- estimateDispersions(Dor_ps.dds)
otu.vst <- getVarianceStabilizedData(Dor_ps.dds)

Dor_ps.prev.vst<-physeq
otu_table(Dor_ps.prev.vst)<- otu_table(otu.vst, taxa_are_rows = TRUE)

Dor_ps.prev.vst
}


