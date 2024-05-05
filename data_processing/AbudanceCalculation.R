if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.18")
library(dada2); packageVersion("dada2")
# Liu
path<-"~/Downloads/CompGenomicsProject/Liu/Liu"
list.files(path)
fnFs <- sort(list.files(path, pattern="_1.fastq$", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq$", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
print(fnFs)
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
print(sample.names)
plotQualityProfile(fnFs[1])
plotQualityProfile(fnRs[1])
filtFs <- file.path(path, "filtered1", paste0(sample.names, "_1_filt.fastq.gz"))
filtRs <- file.path(path, "filtered1", paste0(sample.names, "_2_filt.fastq.gz"))
path <- "~/Downloads/CompGenomicsProject/Liu/Liu/filtered1"
filtFs<- sort(list.files(path, pattern="1_filt.fastq\\.gz$", full.names = TRUE))
filtRs<- sort(list.files(path, pattern="2_filt.fastq\\.gz$", full.names = TRUE))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,maxN=0, maxEE=c(2,2), truncLen=c(165,135), truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=TRUE, verbose=TRUE)
errF <- learnErrors(filtFs, multithread=TRUE, verbose=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE, verbose=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaFs[[1]]
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaRs[[1]]
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
taxa <- assignTaxonomy(seqtab.nochim, "~/Downloads/CompGenomicsProject/Liu/Liu/taxa/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "~/Downloads/CompGenomicsProject/Liu/Liu/taxa/silva_species_assignment_v138.1.fa.gz")
BiocManager::install("phyloseq")
BiocManager::install("ggplot2")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample.names, tax_table(taxa))
tax_table <- tax_table(ps)
tax_info <- apply(tax_table, 1, function(x) paste(na.omit(x), collapse = ";"))
colnames(abundance_matrix) <- tax_info

#Lloyd-Price
path<-"~/Downloads/CompGenomicsProject/Lloyd-Price(HMP2)/Lloyd"
list.files(path)
fnFs2.2<- sort(list.files(path, pattern="^SRR67.*\\_1.fastq\\.gz$", full.names = TRUE))
fnRs2.2<- sort(list.files(path, pattern="^SRR67.*\\_2.fastq\\.gz$", full.names = TRUE))
sample2.2.names <- sapply(strsplit(basename(fnFs2.2), "_"), `[`, 1)
filtFs2.2 <- file.path(path, "filtered1.2", paste0(sample2.2.names, "_1_filt.fastq.gz"))
filtRs2.2 <- file.path(path, "filtered1.2", paste0(sample2.2.names, "_2_filt.fastq.gz"))
path1.2<-"~/Downloads/CompGenomicsProject/Lloyd-Price(HMP2)/Lloyd/filtered1.2"
out2.2 <- filterAndTrim(fnFs2.2, filtFs2.2, fnRs2.2, filtRs2.2,maxN=0, maxEE=c(2,2), truncLen=c(240,160), truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=TRUE, verbose=TRUE)

errF2.2 <- learnErrors(filtFs2.2, multithread=TRUE, verbose=TRUE)
errR2.2 <- learnErrors(filtRs2.2, multithread=TRUE, verbose=TRUE)

derepFs2.2 <- derepFastq(filtFs2.2, verbose=TRUE)
derepRs2.2 <- derepFastq(filtRs2.2, verbose=TRUE)

names(derepFs2.2) <- sample2.2.names
names(derepRs2.2) <- sample2.2.names

dadaFs2.2 <- dada(derepFs2.2, err=errF2.2, multithread=TRUE)
dadaRs2.2 <- dada(derepRs2.2, err=errR2.2, multithread=TRUE)
mergers2.2 <- mergePairs(dadaFs2.2, derepFs2.2, dadaRs2.2, derepRs2.2, verbose=TRUE)
head(mergers2.2[[1]])
seqtab2.2 <- makeSequenceTable(mergers2.2)
dim(seqtab2.2)
table(nchar(getSequences(seqtab2.2)))
seqtab2.2.nochim <- removeBimeraDenovo(seqtab2.2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab2.2.nochim)
sum(seqtab2.2.nochim)/sum(seqtab2.2)
getN <- function(x) sum(getUniques(x))
track2.2 <- cbind(out2.2, sapply(dadaFs2.2, getN), sapply(dadaRs2.2, getN), sapply(mergers2.2, getN), rowSums(seqtab2.2.nochim))
colnames(track2.2) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track2.2) <- sample2.2.names
head(track2.2)
taxa2.2 <- assignTaxonomy(seqtab2.2.nochim, "~/Downloads/CompGenomicsProject/Liu/Liu/taxa/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxaWSpecies2.2 <- addSpecies(taxa2.2, "~/Downloads/CompGenomicsProject/Liu/Liu/taxa/silva_species_assignment_v138.1.fa.gz")
ps2.2 <- phyloseq(otu_table(seqtab2.2.nochim, taxa_are_rows=FALSE), sample2.2.names, tax_table(taxa2.2))

taxa_names_physeq1 <- taxa_names(ps)
taxa_names_physeq2 <- taxa_names(ps2.2)

# Identify overlapping and unique taxa
overlapping_taxa <- intersect(taxa_names_physeq1, taxa_names_physeq2)
unique_taxa_physeq1 <- setdiff(taxa_names_physeq1, taxa_names_physeq2)
unique_taxa_physeq2 <- setdiff(taxa_names_physeq2, taxa_names_physeq1)

merged_physeq <- merge_phyloseq(ps, ps2.2, taxa_are_rows = FALSE, common_by_taxa = TRUE)

# Add unique taxa from physeq1
if (length(unique_taxa_physeq1) > 0) {
  physeq1_unique_taxa <- prune_taxa(unique_taxa_physeq1, ps)
  merged_physeq <- merge_phyloseq(merged_physeq, physeq1_unique_taxa, taxa_are_rows = FALSE, sample = FALSE)
}

# Add unique taxa from physeq2
if (length(unique_taxa_physeq2) > 0) {
  physeq2_unique_taxa <- prune_taxa(unique_taxa_physeq2, ps2.2)
  merged_physeq <- merge_phyloseq(merged_physeq, physeq2_unique_taxa, taxa_are_rows = FALSE, sample = FALSE)
}

abundance_matrix_merged <- as.data.frame(otu_table(merged_physeq))
tax_table <- tax_table(merged_physeq)
tax_info <- apply(tax_table, 1, function(x) paste(na.omit(x), collapse = ";"))
colnames(abundance_matrix_merged) <- tax_info
write.csv(abundance_matrix_merged, file = "abundance_data_with_taxonomy.csv", row.names = TRUE)



psWSpec <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample.names, tax_table(taxaWSpecies))
psWSpec2.2 <- phyloseq(otu_table(seqtab2.2.nochim, taxa_are_rows=FALSE), sample2.2.names, tax_table(taxaWSpecies2.2))
taxa_names_physeq1_spec <- taxa_names(psWSpec)
taxa_names_physeq2_spec <- taxa_names(psWSpec2.2)

# Identify overlapping and unique taxa
overlapping_taxa_spec <- intersect(taxa_names_physeq1_spec, taxa_names_physeq2_spec)
unique_taxa_physeq1_spec <- setdiff(taxa_names_physeq1_spec, taxa_names_physeq2_spec)
unique_taxa_physeq2_spec <- setdiff(taxa_names_physeq2_spec, taxa_names_physeq1_spec)

merged_physeq_spec <- merge_phyloseq(psWSpec, psWSpec2.2, taxa_are_rows = FALSE, common_by_taxa = TRUE)

# Add unique taxa from physeq1
if (length(unique_taxa_physeq1_spec) > 0) {
  physeq1_unique_taxa_spec <- prune_taxa(unique_taxa_physeq1_spec, psWSpec)
  merged_physeq_spec <- merge_phyloseq(merged_physeq_spec, physeq1_unique_taxa_spec, tax = FALSE, sample = FALSE)
}

# Add unique taxa from physeq2
if (length(unique_taxa_physeq2_spec) > 0) {
  physeq2_unique_taxa_spec <- prune_taxa(unique_taxa_physeq2_spec, psWSpec2.2)
  merged_physeq_spec <- merge_phyloseq(merged_physeq_spec, physeq2_unique_taxa_spec, tax = FALSE, sample = FALSE)
}

abundance_matrix_merged_spec <- as.data.frame(otu_table(merged_physeq_spec))
tax_table_spec <- tax_table(merged_physeq_spec)
tax_info_spec <- apply(tax_table_spec, 1, function(x) paste(na.omit(x), collapse = ";"))
colnames(abundance_matrix_merged_spec) <- tax_info_spec
write.csv(abundance_matrix_merged_spec, file = "abundance_data_with_taxonomy_wspecies.csv", row.names = TRUE)


taxa_names_physeq1_WSpec <- taxa_names(psWSpec)
taxa_names_physeq2_WSpec <- taxa_names(psWSpec2.2)
shared_taxa <- intersect(taxa_names_physeq1_WSpec, taxa_names_physeq2_WSpec)
# Subset phyloseq objects to keep only shared taxa
merged_physeq_spec_onlyShared <- merge_phyloseq(physeq1_shared, physeq2_shared, taxa_are_rows = FALSE, common_by_taxa = TRUE)
abundance_matrix_merged_spec_onlyShared <- as.data.frame(otu_table(merged_physeq_spec_onlyShared))
tax_table_spec_onlyShared <- tax_table(merged_physeq_spec_onlyShared)
tax_info_spec_onlyShared <- apply(tax_table_spec_onlyShared, 1, function(x) paste(na.omit(x), collapse = ";"))
colnames(abundance_matrix_merged_spec_onlyShared) <- tax_info_spec_onlyShared
write.csv(abundance_matrix_merged_spec_onlyShared, file = "abundance_data_with_taxonomy_wspecies_onlyShared.csv", row.names = TRUE)

abundance_matrix_psWspec <- as.data.frame(otu_table(psWSpec))
abundance_matrix_psWspec2.2 <- as.data.frame(otu_table(psWSpec2.2))
tax_tablepsWSpec <- tax_table(psWSpec)
tax_infopsWSpec <- apply(tax_tablepsWSpec, 1, function(x) paste(na.omit(x), collapse = ";"))
colnames(abundance_matrix_psWspec) <- tax_infopsWSpec







file_contents <- readLines("/Users/harrisonfried/Downloads/CompGenomicsProject/SRAsUsed.txt")
path<-"/Users/harrisonfried/Downloads/CompGenomicsProject/Gevers1/Gevers1/"
list.files(path)
fnFs3 <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs3 <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
sample3.names <- sapply(strsplit(basename(fnFs3), "_"), `[`, 1)
plotQualityProfile(fnFs3[1])
plotQualityProfile(fnRs3[1])
filtFs3 <- file.path(path, "filtered", paste0(sample3.names, "_1_filt.fastq.gz"))
filtRs3 <- file.path(path, "filtered", paste0(sample3.names, "_2_filt.fastq.gz"))
names(filtFs3) <- sample3.names
names(filtRs3) <- sample3.names
out3 <- filterAndTrim(fnFs3, filtFs3, fnRs3, filtRs3,maxN=0, maxEE=c(2,2), truncLen=c(165,135), truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=TRUE, verbose=TRUE)
errF3 <- learnErrors(filtFs3, multithread=TRUE, verbose=TRUE)
errR3 <- learnErrors(filtRs3, multithread=TRUE, verbose=TRUE)
path<-"/Users/harrisonfried/Downloads/CompGenomicsProject/Gevers1/Gevers1/filtered/"
list.files(path)
filtFs3 <- sort(list.files(path, pattern="_1_filt.fastq.gz", full.names = TRUE))
filtRs3 <- sort(list.files(path, pattern="_2_filt.fastq.gz", full.names = TRUE))
derepFs3.1 <- derepFastq(filtFs3[0:150], verbose=TRUE)
rm(derepRs3.1)
derepRs3.1 <- derepFastq(filtRs3[0:150], verbose=TRUE)
names(derepFs3.1) <- sample3.names[0:150]
names(derepRs3.1) <- sample3.names[0:150]
dadaFs3.1 <- dada(derepFs3.1, err=errF3, multithread=TRUE)
dadaRs3.1 <- dada(derepRs3.1, err=errR3, multithread=TRUE)
mergers3.1 <- mergePairs(dadaFs3.1, derepFs3.1, dadaRs3.1, derepRs3.1, verbose=TRUE)
head(mergers3.1[[1]])
seqtab3.1 <- makeSequenceTable(mergers3.1)
dim(seqtab3.1)
table(nchar(getSequences(seqtab3.1)))
seqtab3.1.nochim <- removeBimeraDenovo(seqtab3.1, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab3.1.nochim)
sum(seqtab3.1.nochim)/sum(seqtab3.1)
getN <- function(x) sum(getUniques(x))
track3.1 <- cbind(out3[0:150], sapply(dadaFs3.1, getN), sapply(dadaRs3.1, getN), sapply(mergers3.1, getN), rowSums(seqtab3.1.nochim))
colnames(track3.1) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track3.1) <- sample.names3.1
head(track3.1)

derepFs3.2 <- derepFastq(filtFs3[151:300], verbose=TRUE)
derepRs3.2 <- derepFastq(filtRs3[151:300], verbose=TRUE)
names(derepFs3.2) <- sample3.names[151:300]
names(derepRs3.2) <- sample3.names[151:300]
dadaFs3.2 <- dada(derepFs3.2, err=errF3, multithread=TRUE)
dadaRs3.2 <- dada(derepRs3.2, err=errR3, multithread=TRUE)
mergers3.2 <- mergePairs(dadaFs3.2, derepFs3.2, dadaRs3.2, derepRs3.2, verbose=TRUE)
head(mergers3.2[[1]])
seqtab3.2 <- makeSequenceTable(mergers3.2)
dim(seqtab3.2)
table(nchar(getSequences(seqtab3.2)))
seqtab3.2.nochim <- removeBimeraDenovo(seqtab3.2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab3.2.nochim)
sum(seqtab3.2.nochim)/sum(seqtab3.2)
getN <- function(x) sum(getUniques(x))
track3.2 <- cbind(out3[151:300], sapply(dadaFs3.2, getN), sapply(dadaRs3.2, getN), sapply(mergers3.2, getN), rowSums(seqtab3.2.nochim))
rownames(track3.2) <- sample3.names[151:300]
head(track3.2)
rm(derepRs3.2)
rm(derepFs3.2)

derepFs3.3 <- derepFastq(filtFs3[301:450], verbose=TRUE)
derepRs3.3 <- derepFastq(filtRs3[301:450], verbose=TRUE)
names(derepFs3.3) <- sample3.names[301:450]
names(derepRs3.3) <- sample3.names[301:450]
dadaFs3.3 <- dada(derepFs3.3, err=errF3, multithread=TRUE)
dadaRs3.3 <- dada(derepRs3.3, err=errR3, multithread=TRUE)
mergers3.3 <- mergePairs(dadaFs3.3, derepFs3.3, dadaRs3.3, derepRs3.3, verbose=TRUE)
head(mergers3.3[[1]])
seqtab3.3 <- makeSequenceTable(mergers3.3)
dim(seqtab3.3)
table(nchar(getSequences(seqtab3.3)))
seqtab3.3.nochim <- removeBimeraDenovo(seqtab3.3, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab3.3.nochim)
sum(seqtab3.3.nochim)/sum(seqtab3.3)
getN <- function(x) sum(getUniques(x))
track3.3 <- cbind(out3[301:450], sapply(dadaFs3.3, getN), sapply(dadaRs3.3, getN), sapply(mergers3.3, getN), rowSums(seqtab3.3.nochim))
rownames(track3.3) <- sample3.names[301:450]
head(track3.3)
rm(derepRs3.3)
rm(derepFs3.3)

derepFs3.4 <- derepFastq(filtFs3[451:600], verbose=TRUE)
derepRs3.4 <- derepFastq(filtRs3[451:600], verbose=TRUE)
names(derepFs3.4) <- sample3.names[451:600]
names(derepRs3.4) <- sample3.names[451:600]
dadaFs3.4 <- dada(derepFs3.4, err=errF3, multithread=TRUE)
dadaRs3.4 <- dada(derepRs3.4, err=errR3, multithread=TRUE)
mergers3.4 <- mergePairs(dadaFs3.4, derepFs3.4, dadaRs3.4, derepRs3.4, verbose=TRUE)
head(mergers3.4[[1]])
seqtab3.4 <- makeSequenceTable(mergers3.4)
dim(seqtab3.4)
table(nchar(getSequences(seqtab3.4)))
seqtab3.4.nochim <- removeBimeraDenovo(seqtab3.4, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab3.4.nochim)
sum(seqtab3.4.nochim)/sum(seqtab3.4)
getN <- function(x) sum(getUniques(x))
track3.4 <- cbind(out3[451:600], sapply(dadaFs3.4, getN), sapply(dadaRs3.4, getN), sapply(mergers3.4, getN), rowSums(seqtab3.4.nochim))
rownames(track3.4) <- sample3.names[451:600]
head(track3.4)
rm(derepRs3.4)
rm(derepFs3.4)

derepFs3.5 <- derepFastq(filtFs3[601:750], verbose=TRUE)
derepRs3.5 <- derepFastq(filtRs3[601:750], verbose=TRUE)
names(derepFs3.5) <- sample3.names[601:750]
names(derepRs3.5) <- sample3.names[601:750]
dadaFs3.5 <- dada(derepFs3.5, err=errF3, multithread=TRUE)
dadaRs3.5 <- dada(derepRs3.5, err=errR3, multithread=TRUE)
mergers3.5 <- mergePairs(dadaFs3.5, derepFs3.5, dadaRs3.5, derepRs3.5, verbose=TRUE)
head(mergers3.5[[1]])
seqtab3.5 <- makeSequenceTable(mergers3.5)
dim(seqtab3.5)
table(nchar(getSequences(seqtab3.5)))
seqtab3.5.nochim <- removeBimeraDenovo(seqtab3.5, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab3.5.nochim)
sum(seqtab3.5.nochim)/sum(seqtab3.5)
getN <- function(x) sum(getUniques(x))
track3.5 <- cbind(out3[601:750], sapply(dadaFs3.5, getN), sapply(dadaRs3.5, getN), sapply(mergers3.5, getN), rowSums(seqtab3.5.nochim))
rownames(track3.5) <- sample3.names[601:750]
head(track3.5)
rm(derepRs3.5)
rm(derepFs3.5)

derepFs3.6 <- derepFastq(filtFs3[751:910], verbose=TRUE)
derepRs3.6 <- derepFastq(filtRs3[751:910], verbose=TRUE)
names(derepFs3.6) <- sample3.names[751:910]
names(derepRs3.6) <- sample3.names[751:910]
dadaFs3.6 <- dada(derepFs3.6, err=errF3, multithread=TRUE)
dadaRs3.6 <- dada(derepRs3.6, err=errR3, multithread=TRUE)
mergers3.6 <- mergePairs(dadaFs3.6, derepFs3.6, dadaRs3.6, derepRs3.6, verbose=TRUE)
head(mergers3.6[[1]])
seqtab3.6 <- makeSequenceTable(mergers3.6)
dim(seqtab3.6)
table(nchar(getSequences(seqtab3.6)))
seqtab3.6.nochim <- removeBimeraDenovo(seqtab3.6, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab3.6.nochim)
sum(seqtab3.6.nochim)/sum(seqtab3.6)
getN <- function(x) sum(getUniques(x))
track3.6 <- cbind(out3[751:910], sapply(dadaFs3.6, getN), sapply(dadaRs3.6, getN), sapply(mergers3.6, getN), rowSums(seqtab3.6.nochim))
rownames(track3.6) <- sample3.names[751:910]
head(track3.6)
rm(derepRs3.6)
rm(derepFs3.6)

taxa3.1 <- assignTaxonomy(seqtab3.1.nochim, "~/Downloads/CompGenomicsProject/Liu/Liu/taxa/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa3.1 <- addSpecies(taxa3.1, "~/Downloads/CompGenomicsProject/Liu/Liu/taxa/silva_species_assignment_v138.1.fa.gz")
ps3.1 <- phyloseq(otu_table(seqtab3.1.nochim, taxa_are_rows=FALSE), sample3.names[0:150], tax_table(taxa3.1))
taxa3.2 <- assignTaxonomy(seqtab3.2.nochim, "~/Downloads/CompGenomicsProject/Liu/Liu/taxa/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa3.2 <- addSpecies(taxa3.2, "~/Downloads/CompGenomicsProject/Liu/Liu/taxa/silva_species_assignment_v138.1.fa.gz")
ps3.2 <- phyloseq(otu_table(seqtab3.2.nochim, taxa_are_rows=FALSE), sample3.names[151:300], tax_table(taxa3.2))
taxa3.3 <- assignTaxonomy(seqtab3.3.nochim, "~/Downloads/CompGenomicsProject/Liu/Liu/taxa/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa3.3 <- addSpecies(taxa3.3, "~/Downloads/CompGenomicsProject/Liu/Liu/taxa/silva_species_assignment_v138.1.fa.gz")
ps3.3 <- phyloseq(otu_table(seqtab3.3.nochim, taxa_are_rows=FALSE), sample3.names[301:450], tax_table(taxa3.3))
taxa3.4 <- assignTaxonomy(seqtab3.4.nochim, "~/Downloads/CompGenomicsProject/Liu/Liu/taxa/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa3.4 <- addSpecies(taxa3.4, "~/Downloads/CompGenomicsProject/Liu/Liu/taxa/silva_species_assignment_v138.1.fa.gz")
ps3.4 <- phyloseq(otu_table(seqtab3.4.nochim, taxa_are_rows=FALSE), sample3.names[451:600], tax_table(taxa3.4))
taxa3.5 <- assignTaxonomy(seqtab3.5.nochim, "~/Downloads/CompGenomicsProject/Liu/Liu/taxa/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa3.5 <- addSpecies(taxa3.5, "~/Downloads/CompGenomicsProject/Liu/Liu/taxa/silva_species_assignment_v138.1.fa.gz")
ps3.5 <- phyloseq(otu_table(seqtab3.5.nochim, taxa_are_rows=FALSE), sample3.names[601:750], tax_table(taxa3.5))
taxa3.6 <- assignTaxonomy(seqtab3.6.nochim, "~/Downloads/CompGenomicsProject/Liu/Liu/taxa/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa3.6 <- addSpecies(taxa3.6, "~/Downloads/CompGenomicsProject/Liu/Liu/taxa/silva_species_assignment_v138.1.fa.gz")
ps3.6 <- phyloseq(otu_table(seqtab3.6.nochim, taxa_are_rows=FALSE), sample3.names[751:910], tax_table(taxa3.6))


taxa_names_physeq3.1 <- taxa_names(ps3.1)
taxa_names_physeq3.2 <- taxa_names(ps3.2)
overlapping_taxa3.1 <- intersect(taxa_names_physeq3.1, taxa_names_physeq3.1)
unique_taxa_physeq3.1 <- setdiff(taxa_names_physeq3.1, taxa_names_physeq3.2)
unique_taxa_physeq3.2 <- setdiff(taxa_names_physeq3.2, taxa_names_physeq3.1)
merged_physeq3.12 <- merge_phyloseq(ps3.1, ps3.2, taxa_are_rows = FALSE, common_by_taxa = TRUE)
if (length(unique_taxa_physeq3.1) > 0) {
  physeq1_unique_taxa3.1 <- prune_taxa(unique_taxa_physeq3.1, ps3.1)
  merged_physeq3.12 <- merge_phyloseq(merged_physeq3.12, physeq1_unique_taxa3.1, taxa_are_rows = FALSE, sample = FALSE)
}
if (length(unique_taxa_physeq3.2) > 0) {
  physeq2_unique_taxa3.2 <- prune_taxa(unique_taxa_physeq3.2, ps3.2)
  merged_physeq3.12 <- merge_phyloseq(merged_physeq3.12, physeq2_unique_taxa3.2, tax = FALSE, sample = FALSE)
}

taxa_names_physeq3.12 <- taxa_names(merged_physeq3.12)
taxa_names_physeq3.3 <- taxa_names(ps3.3)
overlapping_taxa3.123 <- intersect(taxa_names_physeq3.12, taxa_names_physeq3.3)
unique_taxa_physeq3.12 <- setdiff(taxa_names_physeq3.12, taxa_names_physeq3.3)
unique_taxa_physeq3.3 <- setdiff(taxa_names_physeq3.3, taxa_names_physeq3.12)
merged_physeq3.123 <- merge_phyloseq(merged_physeq3.12, ps3.3, taxa_are_rows = FALSE, common_by_taxa = TRUE)
if (length(unique_taxa_physeq3.12) > 0) {
  physeq1_unique_taxa3.12 <- prune_taxa(unique_taxa_physeq3.12, merged_physeq3.12)
  merged_physeq3.123 <- merge_phyloseq(merged_physeq3.123, physeq1_unique_taxa3.12, taxa_are_rows = FALSE, sample = FALSE)
}
if (length(unique_taxa_physeq3.3) > 0) {
  physeq2_unique_taxa3.3 <- prune_taxa(unique_taxa_physeq3.3, ps3.3)
  merged_physeq3.123 <- merge_phyloseq(merged_physeq3.123, physeq2_unique_taxa3.3, tax = FALSE, sample = FALSE)
}

taxa_names_physeq3.123 <- taxa_names(merged_physeq3.123)
taxa_names_physeq3.4 <- taxa_names(ps3.4)
overlapping_taxa3.1234 <- intersect(taxa_names_physeq3.123, taxa_names_physeq3.4)
unique_taxa_physeq3.123 <- setdiff(taxa_names_physeq3.123, taxa_names_physeq3.4)
unique_taxa_physeq3.4 <- setdiff(taxa_names_physeq3.4, taxa_names_physeq3.123)
merged_physeq3.1234 <- merge_phyloseq(merged_physeq3.123, ps3.4, taxa_are_rows = FALSE, common_by_taxa = TRUE)
if (length(unique_taxa_physeq3.123) > 0) {
  physeq1_unique_taxa3.123 <- prune_taxa(unique_taxa_physeq3.123, merged_physeq3.123)
  merged_physeq3.1234 <- merge_phyloseq(merged_physeq3.1234, physeq1_unique_taxa3.123, taxa_are_rows = FALSE, sample = FALSE)
}
if (length(unique_taxa_physeq3.4) > 0) {
  physeq2_unique_taxa3.4 <- prune_taxa(unique_taxa_physeq3.4, ps3.4)
  merged_physeq3.1234 <- merge_phyloseq(merged_physeq3.1234, physeq2_unique_taxa3.4, tax = FALSE, sample = FALSE)
}

taxa_names_physeq3.1234 <- taxa_names(merged_physeq3.1234)
taxa_names_physeq3.5 <- taxa_names(ps3.5)
overlapping_taxa3.12345 <- intersect(taxa_names_physeq3.1234, taxa_names_physeq3.5)
unique_taxa_physeq3.1234 <- setdiff(taxa_names_physeq3.1234, taxa_names_physeq3.5)
unique_taxa_physeq3.5 <- setdiff(taxa_names_physeq3.5, taxa_names_physeq3.1234)
merged_physeq3.12345 <- merge_phyloseq(merged_physeq3.1234, ps3.5, taxa_are_rows = FALSE, common_by_taxa = TRUE)
if (length(unique_taxa_physeq3.1234) > 0) {
  physeq1_unique_taxa3.1234 <- prune_taxa(unique_taxa_physeq3.1234, merged_physeq3.1234)
  merged_physeq3.12345 <- merge_phyloseq(merged_physeq3.12345, physeq1_unique_taxa3.1234, taxa_are_rows = FALSE, sample = FALSE)
}
if (length(unique_taxa_physeq3.5) > 0) {
  physeq2_unique_taxa3.5 <- prune_taxa(unique_taxa_physeq3.5, ps3.5)
  merged_physeq3.12345 <- merge_phyloseq(merged_physeq3.12345, physeq2_unique_taxa3.5, tax = FALSE, sample = FALSE)
}

taxa_names_physeq3.12345 <- taxa_names(merged_physeq3.12345)
taxa_names_physeq3.6 <- taxa_names(ps3.6)
overlapping_taxa3.123456 <- intersect(taxa_names_physeq3.12345, taxa_names_physeq3.6)
unique_taxa_physeq3.12345 <- setdiff(taxa_names_physeq3.12345, taxa_names_physeq3.6)
unique_taxa_physeq3.6 <- setdiff(taxa_names_physeq3.6, taxa_names_physeq3.12345)
merged_physeq3.123456 <- merge_phyloseq(merged_physeq3.12345, ps3.6, taxa_are_rows = FALSE, common_by_taxa = TRUE)
if (length(unique_taxa_physeq3.12345) > 0) {
  physeq1_unique_taxa3.12345 <- prune_taxa(unique_taxa_physeq3.12345, merged_physeq3.12345)
  merged_physeq3.123456 <- merge_phyloseq(merged_physeq3.123456, physeq1_unique_taxa3.12345, taxa_are_rows = FALSE, sample = FALSE)
}
if (length(unique_taxa_physeq3.6) > 0) {
  physeq2_unique_taxa3.6 <- prune_taxa(unique_taxa_physeq3.6, ps3.6)
  merged_physeq3.123456 <- merge_phyloseq(merged_physeq3.123456, physeq2_unique_taxa3.6, tax = FALSE, sample = FALSE)
}

abundance_matrix_merged3.123456 <- as.data.frame(otu_table(merged_physeq3.123456))
tax_table3.123456 <- tax_table(merged_physeq3.123456)
tax_info3.123456 <- apply(tax_table3.123456, 1, function(x) paste(na.omit(x), collapse = ";"))
colnames(abundance_matrix_merged3.123456) <- tax_info3.123456
write.csv(abundance_matrix_merged3.123456, file = "abundance_data_with_taxonomy_wspecies_Gevers1.csv", row.names = TRUE)

taxa_names_physeq12 <- taxa_names(merged_physeq_spec)
taxa_names_physeq3.123456 <- taxa_names(merged_physeq3.123456)
overlapping_taxa123.123456 <- intersect(taxa_names_physeq12, taxa_names_physeq3.123456)
unique_taxa_physeq12 <- setdiff(taxa_names_physeq12, taxa_names_physeq3.123456)
unique_taxa_physeq3.123456 <- setdiff(taxa_names_physeq3.123456, taxa_names_physeq12)
merged_physeq123.123456 <- merge_phyloseq(merged_physeq_spec, merged_physeq3.123456, taxa_are_rows = FALSE, common_by_taxa = TRUE)
if (length(unique_taxa_physeq12) > 0) {
  physeq1_unique_taxa12 <- prune_taxa(unique_taxa_physeq12, merged_physeq_spec)
  merged_physeq123.123456 <- merge_phyloseq(merged_physeq123.123456, physeq1_unique_taxa12, taxa_are_rows = FALSE, sample = FALSE)
}
if (length(unique_taxa_physeq3.123456) > 0) {
  physeq2_unique_taxa3.123456 <- prune_taxa(unique_taxa_physeq3.123456, merged_physeq3.123456)
  merged_physeq123.123456 <- merge_phyloseq(merged_physeq123.123456, physeq2_unique_taxa3.123456, tax = FALSE, sample = FALSE)
}

abundance_matrix_merged123.123456 <- as.data.frame(otu_table(merged_physeq123.123456))
tax_table123.123456 <- tax_table(merged_physeq123.123456)
tax_info123.123456 <- apply(tax_table123.123456, 1, function(x) paste(na.omit(x), collapse = ";"))
colnames(abundance_matrix_merged123.123456) <- tax_info123.123456
write.csv(abundance_matrix_merged123.123456, file = "abundance_data_with_taxonomy_wspecies_LiuLloydGevers1.csv", row.names = TRUE)
