# 02 - Phylogenetic tree of the ASVs (needed for weighted UniFrac)
# Input : 01.dada2_output/seqtab_nochim_16S_axenic_td.rds
# Output: 01.dada2_output/16S_atd_phylogenetic-tree.rds, 16S_atd_GTR_tree.phy
# Versions recorded when run: DECIPHER 2.30.0, phangorn 2.12.1
# Run from the project folder (relative paths).

library(DECIPHER)
library(dada2)
library(phangorn)

seqtab.nochim <- readRDS("./01.dada2_output/seqtab_nochim_16S_axenic_td.rds")
sequences <- getSequences(seqtab.nochim)
names(sequences) <- sequences

# Multiple sequence alignment (DECIPHER)
alignment <- AlignSeqs(DNAStringSet(sequences), anchor = NA)

# Neighbour-joining tree, then a GTR+G+I maximum-likelihood tree started from it
phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
dm     <- dist.ml(phang.align)
treeNJ <- NJ(dm)                      # note: tip order differs from sequence order
fit    <- pml(treeNJ, data = phang.align)
fitGTR <- update(fit, k = 4, inv = 0.2)
fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                    rearrangement = "NNI", control = pml.control(trace = 0))

saveRDS(fitGTR, "./01.dada2_output/16S_atd_phylogenetic-tree.rds")
write.tree(fitGTR$tree, file = "./01.dada2_output/16S_atd_GTR_tree.phy")
