# header -----------------------------------------------------------------------
#
# TITLE:   4. PFAS correlation plot
#
# PURPOSE: visualize and analyze the relationship between exposure variables in a dataset
#
# DATE:    January 20, 2024
#
# CODE REVIEW:
#
# Reviewed by Ellie Holzhausen (EAH) on April 23, 2026
# by Hannan Li (HL) on May 25, 2026
#
# set up -----------------------------------------------------------------------

#clear the workspace
rm(list = ls())

# Load packages
library(dplyr)
library(ggtext)
library(ggcorrplot)
library(ggplot2)
library(here)
library(tidyverse)
library(reshape2)
library(gridExtra)

# ── Load data ──────────────────────────────────────────────────────────────────
df <- read.csv(here::here("out_files", "PFAS1m_micro1m_species.csv")) %>% 
  select(-"X")

# PFAS values in this file are already:
#   1. Below-MDL imputed with MDL/sqrt(2)
#   2. Complete case (participants with NA PFAS excluded)
#   3. Log2-transformed
# See Script 1 (Data Cleaning) for details

# Define PFAS columns
pfas_cols <- c("PFBS_pgmL", "PFHxS_pgmL", "PFNA_pgmL", "PFOA_pgmL", "PFOS_pgmL")

pfas_data <- df %>%
  select(all_of(pfas_cols)) %>%
  dplyr::rename(
    PFBS  = PFBS_pgmL,
    PFHxS = PFHxS_pgmL,
    PFNA  = PFNA_pgmL,
    PFOA  = PFOA_pgmL,
    PFOS  = PFOS_pgmL
  ) %>%
  mutate(across(everything(), as.numeric))

# Set the PFAS order 
pfas_order <- c("PFBS", "PFHxS", "PFOA", "PFOS", "PFNA")
n_pfas <- length(pfas_order) 
pfas_data  <- pfas_data[, pfas_order]

# ── Normality testing ──────────────────────────────────────────────────────────
# Shapiro-Wilk test for each PFAS (H0: data is normally distributed)
# If p < 0.05, normality is rejected and use Spearman
# If all p >= 0.05, normality holds and use Pearson

shapiro_results <- sapply(pfas_order, function(col) {
  shapiro.test(na.omit(pfas_data[[col]]))$p.value
})

shapiro_df <- data.frame(
  PFAS    = names(shapiro_results),
  p_value = round(shapiro_results, 4),
  Normal  = ifelse(shapiro_results >= 0.05, "Yes", "No")
)

print(shapiro_df) # Based on result, Spearman correlation

# ── Correlation matrix ─────────────────────────────────────────────────────────
corr_mat <- round(cor(na.omit(pfas_data), method = "spearman"), 2)
corr_mat[upper.tri(corr_mat)] <- NA
melted_corr_mat <- reshape2::melt(corr_mat, na.rm = TRUE)
colnames(melted_corr_mat) <- c("Var1", "Var2", "value")

# ── P-values ───────────────────────────────────────────────────────────────────
n <- ncol(pfas_data)
p_values <- matrix(NA, nrow = n, ncol = n)
for (i in 1:n) {
  for (j in 1:n) {
    p_values[i, j] <- cor.test(pfas_data[, i], pfas_data[, j],
                               method = "spearman",
                               na.action = na.omit)$p.value
  }
}
# Warnings here because of tied values from MDL imputation. When ties exist, 
# Spearman cannot compute an exact p-value and falls back to a normal approximation instead. So can be ignored.
p_values[upper.tri(p_values)] <- NA
p_vals_melted <- reshape2::melt(p_values, na.rm = TRUE)
colnames(p_vals_melted) <- c("Var1", "Var2", "value")
melted_corr_mat$p_vals <- p_vals_melted$value

# ── Significance labels ────────────────────────────────────────────────────────
add_star <- function(x) ifelse(x < 0.05 & x != 0, "*", "")
melted_corr_mat$label <- paste0(melted_corr_mat$value,
                                add_star(melted_corr_mat$p_vals))

# ── Color palette ──────────────────────────────────────────────────────────────
heatmap_colors <- c("indianred3", "indianred2", "indianred1", "white",
                    "#91BFDB", "#4575B4", "#457")
palette_continuous <- colorRampPalette(heatmap_colors)

# Label colors: PFBS = short-chain (navy), rest = long-chain (navy)
label_colors <- c("navy", "navy", "navy", "navy", "navy")

# ── Heatmap ────────────────────────────────────────────────────────────────────
p_heat <- ggplot(data = melted_corr_mat,
                 aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = label), color = "black", size = 5) +
  scale_fill_gradientn(colors = palette_continuous(100), limits = c(-1, 1)) +
  scale_x_discrete(limits = pfas_order) +
  scale_y_discrete(limits = pfas_order) +
  labs(fill = "Spearman\nCorrelation") +
  theme_bw() +
  theme(
    axis.text.x  = element_text(size = 13, hjust = 0.5,
                                face = "bold", color = label_colors),
    axis.text.y  = element_text(size = 13, face = "bold",
                                color = label_colors),
    axis.title   = element_blank(),
    legend.title = element_text(size = 13, margin = margin(b = 8)),
    legend.text  = element_text(size = 12),
    legend.key.size = unit(2, "lines"),
    panel.grid   = element_blank()
  )

# ── Density plots ──────────────────────────────────────────────────────────────
single_color <- "#F46D43"

density_plots <- lapply(pfas_order, function(col) {
  clean_vals <- na.omit(pfas_data[[col]])
  ggplot(data.frame(value = clean_vals), aes(x = value)) +
    geom_density(fill = single_color, color = NA) +
    ggtitle(col) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 10, face = "bold",
                                    color = label_colors[which(pfas_order == col)]),
          plot.margin = margin(0, 5, 0, 5))
})

# ── Embed density plots directly into heatmap ──────────────────────────────────

# Add top margin to heatmap to make room for density plots above
p_heat2 <- p_heat +
  theme(plot.margin = margin(t = 70, r = 10, b = 5, l = 5))  # increase t to push more space

# Inject each density plot as an annotation grob above its column tile
for (i in seq_along(pfas_order)) {
  grob_i <- ggplotGrob(density_plots[[i]])
  p_heat2 <- p_heat2 +
    annotation_custom(
      grob = grob_i,
      xmin = i - 0.5,   # left edge of tile i
      xmax = i + 0.5,   # right edge of tile i
      ymin = n_pfas + 0.6,   # just above the top tile row
      ymax = n_pfas + 1.3    # height of density; increase to make taller
    )
}

# clip = "off" allows drawing outside the plot panel
p_final <- p_heat2 + coord_cartesian(clip = "off")
print(p_final)

# ── Save ───────────────────────────────────────────────────────────────────────
ggsave("out_figures/PFAS_1m_Correlation_Density.pdf", p_final,
       width = 7, height = 5, units = "in", dpi = 600)

# END