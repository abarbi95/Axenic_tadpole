# 01 - DADA2: demultiplexed reads -> ASV table + taxonomy
# Axenic tadpole 16S rRNA gene amplicons (V3-V4, MiSeq v3 2x300, Macrogen Europe)
#
# Input  : raw_data/<sample>_1.fastq.gz and <sample>_2.fastq.gz  (primers still present)
#          raw_data/silva_nr99_v138.1_train_set.fa.gz
#          raw_data/silva_species_assignment_v138.1.fa.gz
# Output : 01.dada2_output/seqtab_nochim_16S_axenic_td.rds
#          01.dada2_output/tax_atd.rds
# Run from the project folder (all paths are relative).
# dada2 version recorded when this script was run: 1.36.0
# Packages: install beforehand (BiocManager::install("dada2")).

library(dada2)
packageVersion("dada2")

# --- files -------------------------------------------------------------------
fnFs <- sort(list.files("./raw_data", pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files("./raw_data", pattern = "_2.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Quality profiles (used to choose the truncation lengths below)
plotQualityProfile(fnFs[1:2]); plotQualityProfile(fnFs[23:24])
plotQualityProfile(fnRs[1:2]); plotQualityProfile(fnRs[23:24])

# --- filter and trim ---------------------------------------------------------
filtFs <- file.path("./raw_data/dada2_filtered", "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("./raw_data/dada2_filtered", "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# trimLeft removes the primers (forward 17 nt, reverse 21 nt);
# truncLen = 280 / 240 bp where read quality drops.
# On Windows set multithread = FALSE.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(280, 240),
                     maxN = 0, maxEE = c(3, 3), truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE, trimLeft = c(17, 21))
head(out)

# Samples with no reads passing the filter produce no file: drop them
filtFs <- filtFs[file.exists(filtFs)]
filtRs <- filtRs[file.exists(filtRs)]

# --- error model, inference, merging -----------------------------------------
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)

dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

# --- sequence table, length window, chimeras ---------------------------------
seqtab <- makeSequenceTable(mergers)
dim(seqtab)                                    # 42 samples x 4455 ASVs when run
table(nchar(getSequences(seqtab)))

# keep amplicons of the expected length (removes non-specific products)
seqtab2 <- seqtab[, nchar(colnames(seqtab)) %in% 390:430]

seqtab.nochim <- removeBimeraDenovo(seqtab2, method = "consensus",
                                    multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)                             # 42 x 1354 when run
sum(seqtab.nochim) / sum(seqtab2)              # 0.9699: fraction of reads kept after chimera removal

dir.create("./01.dada2_output", showWarnings = FALSE)
saveRDS(seqtab.nochim, "./01.dada2_output/seqtab_nochim_16S_axenic_td.rds")

# --- taxonomy (naive Bayesian classifier, SILVA v138.1) -----------------------
taxa <- assignTaxonomy(seqtab.nochim, "./raw_data/silva_nr99_v138.1_train_set.fa.gz",
                       multithread = TRUE)
taxa <- addSpecies(taxa, "./raw_data/silva_species_assignment_v138.1.fa.gz",
                   allowMultiple = TRUE)
saveRDS(taxa, "./01.dada2_output/tax_atd.rds")
