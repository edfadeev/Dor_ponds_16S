#load libraries and set random seed
library(ggplot2); packageVersion("ggplot2")
library(dada2); packageVersion("dada2")
library(gridExtra); packageVersion("gridExtra")
library(DECIPHER); packageVersion("DECIPHER")
library(phangorn); packageVersion("phangorn")
set.seed(123)


# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files("Clipped", pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files("Clipped", pattern="_R2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_clip_R1.fastq"), `[`, 1)


# Make directories and filenames for the filtered fastqs
filt_path <- file.path("Report")
if(!file_test("-d", filt_path)) dir.create(filt_path)

filt_path <- file.path("Filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filt_path <- file.path("Seq.Tables")
if(!file_test("-d", filt_path)) dir.create(filt_path)

filtFs <- file.path("Filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("Filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#separate the different runs
#Hiseq
fnFs_hi<- sort(file.path("Clipped",paste(c(1:69), "_clip_R1.fastq", sep = "")))
fnRs_hi<- sort(file.path("Clipped",paste(c(1:69), "_clip_R2.fastq", sep = ""))) 
filtFs_hi<- sort(file.path("Filtered",paste(c(1:69), "_F_filt.fastq.gz", sep = "")))
filtRs_hi<- sort(file.path("Filtered",paste(c(1:69), "_R_filt.fastq.gz", sep = ""))) 

# quality check run1
QualityProfileFs <- list()
for(i in 1:length(fnFs_hi)) {
  QualityProfileFs[[i]] <- list()
  QualityProfileFs[[i]][[1]] <- plotQualityProfile(fnFs_hi[i])
}
pdf(file.path("Report","RawProfileForward_run1.pdf"))
for(i in 1:length(fnFs_hi)) {
  do.call("grid.arrange", QualityProfileFs[[i]])  
}
dev.off()
rm(QualityProfileFs)

QualityProfileRs <- list()
for(i in 1:length(fnRs_hi)) {
  QualityProfileRs[[i]] <- list()
  QualityProfileRs[[i]][[1]] <- plotQualityProfile(fnRs_hi[i])
}
pdf(file.path("Report","RawProfileReverse_run1.pdf"))
for(i in 1:length(fnRs_hi)) {
  do.call("grid.arrange", QualityProfileRs[[i]])  
}
dev.off()
rm(QualityProfileRs)

#generate aggregated quality overview
ggsave(file.path("Report","RawProfileForward_run1_agg.pdf"),
       plot = plotQualityProfile(fnFs_hi, aggregate = TRUE),
       width = 20, height = 20, units = "cm")

ggsave(file.path("Report","RawProfileReverse_run1_agg.pdf"),
       plot = plotQualityProfile(fnRs_hi, aggregate = TRUE),
       width = 20, height = 20, units = "cm")


#Filter and trim
out_hi <- filterAndTrim(fnFs_hi, filtFs_hi, fnRs_hi, filtRs_hi, truncLen=c(235,235),
                        maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
                        compress=TRUE, multithread=TRUE, verbose = TRUE)

#remove samples that were filtered out completele
filtFs_hi<-sort(file.path("Filtered",row.names(as.data.frame(out_hi))[as.data.frame(out_hi)$reads.out >0]))
filtFs_hi<-gsub("_clip_R1.fastq","_F_filt.fastq.gz", filtFs_hi)
filtRs_hi<- gsub("_F","_R",filtFs_hi)

#generate aggregated quality overview
ggsave(file.path("Report","FiltProfileForward_run1_agg.pdf"),
       plot = plotQualityProfile(filtFs_hi, aggregate = TRUE),
       width = 20, height = 20, units = "cm")

ggsave(file.path("Report","FiltProfileReverse_run1_agg.pdf"),
       plot = plotQualityProfile(filtRs_hi, aggregate = TRUE),
       width = 20, height = 20, units = "cm")

# Learn errors 
errF_hi <- learnErrors(filtFs_hi, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
errR_hi <- learnErrors(filtRs_hi, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
# Sample Inference 
dadaFs_hi <- dada(filtFs_hi, err=errF_hi, multithread=TRUE, verbose = TRUE)
dadaRs_hi <- dada(filtRs_hi, err=errR_hi, multithread=TRUE, verbose = TRUE)

#Merge paired reads
mergers_hi <- mergePairs(dadaFs_hi, filtFs_hi, dadaRs_hi, filtRs_hi, verbose=TRUE, minOverlap = 10)
#generate sequence table and save it
seqtab_hi <- makeSequenceTable(mergers_hi)
saveRDS(seqtab_hi, file.path("Seq.Tables","seqtab_hi.rds"))

#MiSeq 1 - 493:307
fnFs_mi1 <- sort(file.path("Clipped",paste(c(70:105), "_clip_R1.fastq", sep = "")))
fnRs_mi1 <- sort(file.path("Clipped",paste(c(70:105), "_clip_R2.fastq", sep = "")))
filtFs_mi1 <- sort(file.path("Filtered",paste(c(70:105), "_F_filt.fastq.gz", sep = "")))
filtRs_mi1 <- sort(file.path("Filtered",paste(c(70:105), "_R_filt.fastq.gz", sep = "")))

# quality check run2
QualityProfileFs <- list()
for(i in 1:length(fnFs_mi1)) {
  QualityProfileFs[[i]] <- list()
  QualityProfileFs[[i]][[1]] <- plotQualityProfile(fnFs_mi1[i])
}
pdf(file.path("Report","RawProfileForward_run2.pdf"))
for(i in 1:length(fnFs_mi1)) {
  do.call("grid.arrange", QualityProfileFs[[i]])  
}
dev.off()
rm(QualityProfileFs)

QualityProfileRs <- list()
for(i in 1:length(fnRs_mi1)) {
  QualityProfileRs[[i]] <- list()
  QualityProfileRs[[i]][[1]] <- plotQualityProfile(fnRs_mi1[i])
}
pdf(file.path("Report","RawProfileReverse_run2.pdf"))
for(i in 1:length(fnRs_mi1)) {
  do.call("grid.arrange", QualityProfileRs[[i]])  
}
dev.off()
rm(QualityProfileRs)

#generate aggregated quality overview
ggsave(file.path("Report","RawProfileForward_run2_agg.pdf"),
       plot = plotQualityProfile(fnFs_mi1, aggregate = TRUE),
       width = 20, height = 20, units = "cm")

ggsave(file.path("Report","RawProfileReverse_run2_agg.pdf"),
       plot = plotQualityProfile(fnRs_mi1, aggregate = TRUE),
       width = 20, height = 20, units = "cm")

#Filter and trim
out_mi1 <- filterAndTrim(fnFs_mi1, filtFs_mi1, fnRs_mi1, filtRs_mi1, truncLen=c(235,235),
                         maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE, verbose = TRUE)


#remove samples that were filtered out completely
filtFs_mi1<-sort(file.path("Filtered",row.names(as.data.frame(out_mi1))[as.data.frame(out_mi1)$reads.out >0]))
filtFs_mi1<-gsub("_clip_R1.fastq","_F_filt.fastq.gz", filtFs_mi1)
filtRs_mi1<- gsub("_F","_R",filtFs_mi1)

#generate aggregated quality overview
ggsave(file.path("Report","FiltProfileForward_run2_agg.pdf"),
       plot = plotQualityProfile(filtFs_mi1, aggregate = TRUE),
       width = 20, height = 20, units = "cm")

ggsave(file.path("Report","FiltProfileReverse_run2_agg.pdf"),
       plot = plotQualityProfile(filtRs_mi1, aggregate = TRUE),
       width = 20, height = 20, units = "cm")

# Learn errors
errF_mi1 <- learnErrors(filtFs_mi1, nbases = 2e+08, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
errR_mi1 <- learnErrors(filtRs_mi1, nbases = 2e+08, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
# Sample Inference 
dadaFs_mi1 <- dada(filtFs_mi1, err=errF_mi1, multithread=TRUE, verbose = TRUE)
dadaRs_mi1 <- dada(filtRs_mi1, err=errR_mi1, multithread=TRUE, verbose = TRUE)
#Merge paired reads
mergers_mi1 <- mergePairs(dadaFs_mi1, filtFs_mi1, dadaRs_mi1, filtRs_mi1, verbose=TRUE, minOverlap = 10)
#generate sequence table and save it
seqtab_mi1 <- makeSequenceTable(mergers_mi1)
saveRDS(seqtab_mi1, file.path("Seq.Tables","seqtab_mi1.rds"))

# Plot error profiles
pdf(file.path("Report","ErrorProfiles.pdf"))
plotErrors(errF_hi, nominalQ = TRUE)+ggtitle("Run1-Forward")
plotErrors(errR_hi, nominalQ = TRUE)+ggtitle("Run1-Reverse")
plotErrors(errF_mi1, nominalQ = TRUE)+ggtitle("Run2-Forward")
plotErrors(errR_mi1, nominalQ = TRUE)+ggtitle("Run2-Reverse")
dev.off()

#write out filtered read counts
write.csv(rbind(out_hi,out_mi1), 
          file= file.path("Report","dada2_filterAndTrim_output.csv"))



#merge the two runs into a single sequence table
seqtab<- mergeSequenceTables(table1= makeSequenceTable(mergers_hi), 
                             table2 = makeSequenceTable(mergers_mi1))

#Combine together sequences that are identical 
seqtab1 <- collapseNoMismatch(seqtab, verbose = TRUE)

dim(seqtab1)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab1)))

seqtab.nochim <- removeBimeraDenovo(seqtab1, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#proportion of chimeras
sum(seqtab.nochim)/sum(seqtab1)

# inspect output: remove singletons and 'junk' sequences
# read lengths modified for V3V4 amplicons / based upon output table where majority of reads occurs
seqtab.nochim2 <- seqtab.nochim[, nchar(colnames(seqtab.nochim)) %in% c(400:480) & colSums(seqtab.nochim) > 1]
dim(seqtab.nochim2)
summary(rowSums(seqtab.nochim2)/rowSums(seqtab.nochim))

#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim2, "./tax/silva_nr99_v138_train_set.fa.gz", multithread=TRUE, tryRC = TRUE, verbose = TRUE)
taxa <- addSpecies(taxa, "./tax/silva_species_assignment_v138.fa.gz", tryRC = TRUE, verbose = TRUE)

# get summary tables 
sample.order <- c(basename(filtFs_mi1),basename(filtFs_hi))

getN <- function(x) sum(getUniques(x))
track <- cbind(rbind(as.data.frame(out_mi1)[as.data.frame(out_mi1)$reads.out >0,],as.data.frame(out_hi)[as.data.frame(out_hi)$reads.out >0,]), 
               sapply(c(dadaFs_mi1,dadaFs_hi), getN),
               sapply(c(mergers_mi1,mergers_hi), getN),
               rowSums(seqtab.nochim[sample.order,]),
               rowSums(seqtab.nochim2[sample.order,]))
colnames(track) <- c("input", "filtered", "denoised", "merged", "nochim", "tabled")
rownames(track) <- gsub("_F_filt.fastq.gz","",sample.order)
track <- data.frame(track)

#add unclassified levels of taxonomy 
TAX <- taxa
k <- ncol(TAX) - 1
for (i in 2:k) {
  if (sum(is.na(TAX[, i])) > 1) {
    test <- TAX[is.na(TAX[, i]), ]
    for (j in 1:nrow(test)) {
      if (sum(is.na(test[j, i:(k + 1)])) == length(test[j, i:(k + 1)])) {
        test[j, i] <- paste(test[j, (i - 1)], "_uncl", sep = "")
        test[j, (i + 1):(k + 1)] <- test[j, i]
      }
    }
    TAX[is.na(TAX[, i]), ] <- test
  }
  if (sum(is.na(TAX[, i])) == 1) {
    test <- TAX[is.na(TAX[, i]), ]
    if (sum(is.na(test[i:(k + 1)])) == length(test[i:(k + 1)])) {
      test[i] <- paste(test[(i - 1)], "_uncl", sep = "")
      test[(i + 1):(k + 1)] <- test[i]
    }
    TAX[is.na(TAX[, i]),] <- test
  }
}
TAX[is.na(TAX[, (k + 1)]), (k + 1)] <- paste(TAX[is.na(TAX[, (k + 1)]), k], "_uncl", sep = "")

# write output
write.table(t(seqtab.nochim2), "dada2_seqtab_nochim2.txt", quote = F, sep = "\t")
write.table(TAX, "dada2_taxonomy_table.txt", sep = "\t", quote = F)
write.table(track, "libs_summary_table.txt", sep = "\t", quote = F)
save.image("V3V4_dada2_output_object.Rdata")
