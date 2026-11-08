# 03 - phyloseq object, taxonomic filtering, decontamination (negative controls)
# Input : 01.dada2_output/{seqtab_nochim_16S_axenic_td.rds, tax_atd.rds, 16S_atd_phylogenetic-tree.rds}
#         axenic_tadpole_metadata1.csv, add_refseq.R
# Output: 02.phyloseq_output/physeq_atd.rds (+ ASV sequences as .rds/.csv)
# Versions recorded when run: phyloseq 1.42.0
# Run from the project folder (relative paths).

library(phyloseq)
library(Biostrings)
library(ape)
library(readr)
library(tibble)
library(dplyr)
library(ggplot2)
library(decontam)       # BiocManager::install("decontam")

# add_refseq() comes from the microbiomeutilities package (Shetty & Lahti); see README
source("add_refseq.R")

# --- load data ---------------------------------------------------------------
seqtab.nochim <- readRDS("./01.dada2_output/seqtab_nochim_16S_axenic_td.rds")
taxa          <- readRDS("./01.dada2_output/tax_atd.rds")
fitGTR        <- readRDS("./01.dada2_output/16S_atd_phylogenetic-tree.rds")

map <- read_csv("axenic_tadpole_metadata1.csv") %>% column_to_rownames("sample_ID")

# The metadata has 43 rows but sample Andrea8 (CA8) has no sequence data (42 samples
# in the sequence table); phyloseq() keeps only the samples present in both.
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               tax_table(taxa), phy_tree(fitGTR$tree), sample_data(map))
ps <- add_refseq(ps)                      # ASV sequences -> refseq slot; ASVs renamed ASV1, ASV2, ...

# Root the tree (needed for UniFrac); the outgroup is chosen at random, seed fixed
set.seed(120)
phy_tree(ps) <- root(phy_tree(ps), sample(taxa_names(ps), 1), resolve.root = TRUE)
is.rooted(phy_tree(ps))

# --- taxonomic filtering -----------------------------------------------------
print(table(tax_table(ps)[, "Kingdom"], useNA = "ifany"))   # check: only Bacteria expected
print(table(tax_table(ps)[, "Phylum"],  useNA = "ifany"))

# Remove ASVs without a phylum assignment (18 when run; probably artefacts)
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

# Remove chloroplasts and mitochondria.
# NB: subset_taxa() also drops ASVs whose Order / Family is NA (NA comparisons count as FALSE).
ps <- subset_taxa(ps, !Order == "Chloroplast")
ps <- subset_taxa(ps, !Family == "Mitochondria")

# --- contaminants identified from the negative controls (NC = NTCs) -----------
sample_data(ps)$is.neg <- sample_data(ps)$group == "NC"
contamdf.prev <- isContaminant(ps, method = "prevalence", neg = "is.neg")   # default threshold 0.1
table(contamdf.prev$contaminant)                                            # 8 contaminants when run

# Prevalence in true samples vs negative controls (QC plot)
ps.neg.pres <- transform_sample_counts(prune_samples(sample_data(ps)$group == "NC", ps),
                                       function(x) 1 * (x > 0))
ps.pos.pres <- transform_sample_counts(prune_samples(sample_data(ps)$group != "NC", ps),
                                       function(x) 1 * (x > 0))
df.pres <- data.frame(prevalence.pos = taxa_sums(ps.pos.pres),
                      prevalence.neg = taxa_sums(ps.neg.pres),
                      contaminant    = contamdf.prev$contaminant)
ggplot(df.pres, aes(x = prevalence.neg, y = prevalence.pos, color = contaminant)) +
  geom_point() + theme_bw()

ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps)
ps.noncontam <- subset_samples(ps.noncontam, group != "NC")               # drop NC samples
ps.noncontam <- filter_taxa(ps.noncontam, function(x) sum(x) > 0, TRUE)   # drop ASVs absent from all remaining samples

# Keep ASVs seen more than 3 times in more than 5% of the samples
ps.filtered <- filter_taxa(ps.noncontam, function(x) sum(x > 3) > (0.05 * length(x)), TRUE)
cat("Samples:", nsamples(ps.filtered), " ASVs:", ntaxa(ps.filtered), "\n")
# manuscript: 921 ASVs in 39 samples (check whether that figure is before or after rarefaction, script 04)

# --- save --------------------------------------------------------------------
dir.create("./02.phyloseq_output", showWarnings = FALSE)
ASV_sequence_df <- as.data.frame(refseq(ps.filtered))
saveRDS(ASV_sequence_df, "./02.phyloseq_output/ASV_atd_sequences.rds")
write.csv(ASV_sequence_df, "./02.phyloseq_output/ASV_atd_sequences.csv")
saveRDS(ps.filtered, "./02.phyloseq_output/physeq_atd.rds")
