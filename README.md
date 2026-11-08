# Analysis code - "An axenic tadpole model for the study of the amphibian microbiome"

Scripts for the *Discoglossus pictus* axenic tadpole study. Run each script from the
project folder (paths are relative). Scripts are numbered in the order they are run.

| Script | What it does | Main inputs | Main outputs |
|---|---|---|---|
| `01_dada2_16S.R` | 16S V3-V4 reads -> ASV table, SILVA 138.1 taxonomy | raw FASTQ (ENA), SILVA files | `01.dada2_output/*.rds` |
| `02_phylogenetic_tree.R` | Alignment (DECIPHER) and GTR+G+I tree (phangorn) | ASV table | `16S_atd_phylogenetic-tree.rds` |
| `03_phyloseq_filtering.R` | phyloseq object, taxonomic filtering, decontam with NTCs | `axenic_tadpole_metadata1.csv`, `add_refseq.R` | `02.phyloseq_output/physeq_atd.rds` |
| `04_alpha_beta_diversity.Rmd` | Rarefaction, alpha/beta diversity, composition, Table 2 | `physeq_atd.rds`, `Wd.R` | tables and figures |
| `05_bacterial_DNA_load.R` | qPCR bacterial DNA load, Kruskal-Wallis, pairwise Wilcoxon | `qPCR_bload_CpT_values.csv` | Figure 2 |
| `06_body_size.R` | Length and weight (Gamma GLM, emmeans) | `Axenic tadpoles measurements .csv` | size plots |
| `07_DE_RNAseq.R` | RNA-seq differential expression (edgeR) and PCA | Salmon quantifications, `AxenicTadpole_RNA_metadata.csv`, `tx2gene_17709.txt` | DE tables, PCA |

## Third-party code (included unchanged)
- `Wd.R`: Wd* test for heterogeneous dispersions, from Hamidi et al. (2019), *Microbiome* 7:51 (supplementary file).
- `add_refseq.R`: `add_refseq()` from the microbiomeutilities R package (Shetty & Lahti).
- The RNA-seq script was originally written by Kevin Mulder and adapted by Andrea Barbi.

## Notes
- Software versions are printed with `sessionInfo()` (Rmd and script 07); the versions recorded when the
  analyses were run were dada2 1.36.0, DECIPHER 2.30.0, phangorn 2.12.1, phyloseq 1.42.0, R 4.4.2.
- Permutation tests (script 04) use a fixed seed; p-values may differ marginally from the paper.
- Raw reads: ENA accession (to be added).
