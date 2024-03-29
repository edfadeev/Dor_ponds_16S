require(DESeq2)
library(phyloseq); packageVersion("phyloseq")

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

#scale parameters
scale_par <- function(x) scale(x, center = FALSE, scale = TRUE)[,1]

#define presence-absence function #https://github.com/hms-dbmi/UpSetR/issues/85
pres_abs_matrix <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}

#calculate standard error
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
