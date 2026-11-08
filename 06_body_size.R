# 06 - Body length and weight of the tadpoles: comparison between groups
# Input : Axenic tadpoles measurements.csv  (columns: Group, Length, Weight, Sample.ID;
#         groups AA, AP, CA, CP). 
# Output: results_size/size_models.txt, Weight_boxplot.{tiff,pdf}, Length_boxplot.{tiff,pdf}
# Run from the project folder (relative paths).

library(emmeans)
library(ggplot2)
library(ggsci)

dir.create("results_size", showWarnings = FALSE)
data <- read.csv("Axenic tadpoles measurements.csv")
data$Group <- factor(data$Group, levels = c("AA", "AP", "CA", "CP"))

# --- Gamma GLMs + pairwise comparisons (emmeans) -------------------------------
m1 <- glm(Length ~ Group, data = data, family = Gamma)
m2 <- glm(Weight ~ Group, data = data, family = Gamma)

sink("results_size/size_models.txt")
cat("Length\n"); print(summary(m1)); print(emmeans(m1, pairwise ~ Group))
cat("\nWeight\n"); print(summary(m2)); print(emmeans(m2, pairwise ~ Group))
sink()

# Model assumptions (manuscript: checked with the 'performance' package)
# performance::check_model(m1); performance::check_model(m2)   # added; not in the original script

# --- figures -----------------------------------------------------------------
size_plot <- function(var, ylab) {
  ggplot(data, aes(y = .data[[var]], x = Group, color = Group, fill = Group)) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA) +
    geom_point(alpha = 0.6, size = 4) +
    xlab("Group") + ylab(ylab) +
    theme_bw() + scale_color_jco() + scale_fill_jco() +
    theme(axis.line = element_line(color = "black"),
          text = element_text(size = 17),
          axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
          axis.title.y = element_text(face = "bold", margin = margin(r = 10)))
}
p_w <- size_plot("Weight", "Weight (g)")
# CHECK the unit: the original label says mm, but the values (about 1-1.6) look like cm
p_l <- size_plot("Length", "Length (mm)")

for (nm in c("Weight", "Length")) {
  p <- if (nm == "Weight") p_w else p_l
  ggsave(sprintf("results_size/%s_boxplot.tiff", nm), p, width = 8, height = 6, units = "in",
         dpi = 600, compression = "lzw")
  ggsave(sprintf("results_size/%s_boxplot.pdf", nm), p, width = 8, height = 6, units = "in")
}
