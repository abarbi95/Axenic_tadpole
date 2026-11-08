# 05 - Bacterial DNA load (16S rRNA gene copies per tadpole, qPCR): comparison between groups
# Input : qPCR_bload_CpT_values.csv  (columns 1-4: ID, SQ = starting quantity, group, CpT)
# Output: results_bload/BLoad_boxplot.{tiff,pdf}, load_statistics.txt
# Run from the project folder (relative paths).
#
# CpT (copies per tadpole) = (SQ - mean SQ of the NTCs) x 50   [x 50: 100 uL eluate / 2 uL template]
# Group codes in the raw file: A = AA, AP = AP, X = CA, P = CP, NC = NTC.

library(readr)
library(dplyr)
library(ggplot2)
library(ggsci)
library(rstatix)
library(scales)

dir.create("results_bload", showWarnings = FALSE)

datab <- read_csv("qPCR_bload_CpT_values.csv", col_types = cols(group = col_character()))
# use read_csv2() if the columns come out merged (";" separator, decimal commas)
datab <- datab[, 1:4]
names(datab) <- c("ID", "SQ", "group", "CpT")
stopifnot(is.numeric(datab$SQ), is.numeric(datab$CpT))

# Recompute CpT from SQ at full precision and drop the NTCs
nc_mean <- mean(datab$SQ[datab$group == "NC"])
datab <- datab %>%
  filter(group != "NC") %>%
  mutate(CpT   = (SQ - nc_mean) * 50,
         group = recode(group, A = "AA", X = "CA", P = "CP"),
         group = factor(group, levels = c("AA", "AP", "CA", "CP")))
table(datab$group, useNA = "ifany")     # expected when run: AA 10, AP 10, CA 10, CP 9
# (the 16S metadata lists 10 CP tadpoles with a load: check n against the manuscript)
# NB: CpT <= 0 cannot be shown on the log axis

# --- statistics --------------------------------------------------------------
sink("results_bload/load_statistics.txt")
cat("Shapiro-Wilk per group (CpT)\n")
print(by(datab$CpT, datab$group, shapiro.test))
kw <- kruskal.test(CpT ~ group, data = datab)
print(kw)
stat.test <- datab %>% wilcox_test(CpT ~ group, p.adjust.method = "BH")   # pairwise Wilcoxon, BH
print(as.data.frame(stat.test))
sink()
kw_label <- paste0("Kruskal-Wallis, p = ", format.pval(kw$p.value, digits = 3))

# Significant pairs only; brackets stacked on the log axis, shortest span first
lv    <- levels(datab$group)
y_top <- max(datab$CpT)
stat.plot <- stat.test %>%
  filter(p.adj < 0.05) %>%
  mutate(span = match(group2, lv) - match(group1, lv)) %>%
  arrange(span) %>%
  mutate(y.position = y_top * 10^(0.3 * row_number()))
br <- stat.plot %>%
  mutate(x1 = match(group1, lv), x2 = match(group2, lv), tip = y.position / 10^0.08)

# --- figure ------------------------------------------------------------------
group_cols <- setNames(pal_jco()(4), lv)
brks <- 10^(floor(log10(min(datab$CpT))):ceiling(log10(max(datab$CpT))))

p <- ggplot(datab, aes(x = group, y = CpT, color = group, fill = group)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_point(size = 4, alpha = 0.5, position = position_jitter(width = 0.2, seed = 1)) +
  # mean +/- SEM, computed on log10 values because of the log scale
  stat_summary(fun.data = function(y) data.frame(y = mean(y),
                                                 ymin = mean(y) - sd(y) / sqrt(length(y)),
                                                 ymax = mean(y) + sd(y) / sqrt(length(y))),
               geom = "errorbar", width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5, shape = 18) +
  geom_segment(data = br, aes(x = x1, xend = x2, y = y.position, yend = y.position, color = group1),
               inherit.aes = FALSE, linewidth = 0.7) +
  geom_segment(data = br, aes(x = x1, xend = x1, y = y.position, yend = tip, color = group1),
               inherit.aes = FALSE, linewidth = 0.7) +
  geom_segment(data = br, aes(x = x2, xend = x2, y = y.position, yend = tip, color = group1),
               inherit.aes = FALSE, linewidth = 0.7) +
  geom_text(data = br, aes(x = (x1 + x2) / 2, y = y.position, label = p.adj.signif, color = group1),
            inherit.aes = FALSE, vjust = 0.1, size = 6) +
  annotate("text", x = -Inf, y = Inf, label = kw_label, hjust = -0.05, vjust = 1.5,
           size = 5.5, color = "black") +
  scale_y_log10(breaks = brks, labels = trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(min(datab$CpT) * 0.7, max(stat.plot$y.position) * 10^0.5)) +
  scale_color_manual(values = group_cols) +
  scale_fill_manual(values = group_cols) +
  xlab("Group") + ylab("Bacterial DNA Load (copies per tadpole)") +
  theme_bw() +
  theme(axis.line = element_line(color = "black"),
        text = element_text(size = 17),
        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(face = "bold", margin = margin(r = 10)))
p

ggsave("results_bload/BLoad_boxplot.tiff", p, width = 8, height = 6, units = "in",
       dpi = 600, compression = "lzw")
ggsave("results_bload/BLoad_boxplot.pdf", p, width = 8, height = 6, units = "in")

# --- Crawford-Howell test: is the outlying AA sample (AA3) different in DNA load from the other AA samples? ---
# (single case vs. small control sample, on log10 values; reported in the Results)
crawford_howell <- function(case, controls) {
  n <- length(controls)
  t <- (case - mean(controls)) / (sd(controls) * sqrt((n + 1) / n))
  p <- 2 * pt(-abs(t), df = n - 1)
  c(t = t, df = n - 1, p = p)
}
# Uses the 16S metadata file (columns sample_ID, at_ID, group, bacterial_load); AA3 = at_ID "AA3" (Andrea23).
meta <- read.csv("axenic_tadpole_metadata1.csv")
aa <- meta[meta$group == "AA" & !is.na(meta$bacterial_load), ]
print(crawford_howell(log10(aa$bacterial_load[aa$at_ID == "AA3"]),
                      log10(aa$bacterial_load[aa$at_ID != "AA3"])))
# expected from the manuscript: t(8) = -1.04, P = 0.33
