
# ==========================================================================
# MANUSCRIPT FIGURES  —  single self-contained, portable generator.
# Unzip this folder, then run FROM INSIDE it (paths are resolved with here()):
#       Rscript flux_manuscript_figures.R
# Reads input from ./data/, writes figures to ./figures/ and tables to ./tables/.
# R packages: dplyr ggplot2 ggrepel lme4 emmeans car patchwork openxlsx here,
#   plus ComplexHeatmap + circlize (Bioconductor).
#
#   PART 0  Cohort matching & consistency checks
#   PART 1  Observed vs predicted growth, 3 media (origin fit, n = 66)
#   PART 2  Jaccard capability-similarity clustering heatmap (assayed species)
#   PART 3  Four capability box plots with capability x dose mixed-model test
#
#  Code review: Daniel Fassler on August 3, 2026
#               Rechecked by Devendra Paudel on August 3, 2026
#
# --------------------------------------------------------------------------
#
# set up -----------------------------------------------------------------------

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(ggrepel)
  library(lme4); library(emmeans); library(car)
  library(ComplexHeatmap); library(circlize)
  library(grid); library(openxlsx); library(patchwork)
  library(stringr); library(here)
})

# setwd() removed for portability: run this script with the working
# directory set to the project root (the folder containing this "code" dir)

# Input data folder
DAT_YUAN     <- "invitro_data"       # Yuan in-vitro data: formated_data.csv
DAT_JOHANNES <- "insilico_data"   # Johannes FVA files: FVA_*.csv

# ---- FVA file suffix control ---------------------------------------------
FVA_SUFFIX        <- "_73"   # used by heatmap (Part 2) + box plots (Part 3)
FVA_SUFFIX_GROWTH <- "_73"   # Part 1 growth correlation on the 73-species set (matches "73 models / 66 tested")
# --------------------------------------------------------------------------

# Output folders
FIG <- here::here("out_figures")
TAB <- here::here("out_files")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)


m_adapted <- read.csv(file.path(TAB, "formated_data_adapted.csv"), stringsAsFactors = FALSE)

# Phylum annotation

## ── 0b. Resolve species names: GTDB preferred, SILVA fallback ─────────────
# This matches the in-vitro script logic exactly and fixes 43 SILVA errors
# (most critically: JEB00298 = B. adolescentis mislabeled as B. caccae in SILVA)
gtdb_col <- grep("GTDB", names(m_adapted), value = TRUE)[1]


m_adapted <- m_adapted %>%
  mutate(
    phylum = if_else(
      grepl("p__", .data[[gtdb_col]], fixed = TRUE),
      sub(".*p__([^;]+);.*", "\\1", .data[[gtdb_col]]),
      NA_character_
    )
  )
m_adapted$phylum[m_adapted$Species %in%
               c("Butyricimonas synergistica",
                 "Parabacteroides distasonis")] <- "Bacteroidota"

cat("  Unique species after GTDB resolution:", length(unique(m_adapted$Species)), "\n")

## ── 0d. Cohort species matching ────────────────────────────────────────────
# Read cohort species list (same source as in-vitro script)
# Expects cohort_species.csv in out_files/ with column 'species_name'
# Generate from in-vitro script:
#   write.csv(cohort_species, here::here('out_files', 'cohort_species.csv'), row.names = FALSE)
cohort_file <- here::here("out_files", "cohort_species.csv")
if (file.exists(cohort_file)) {
  cohort_sp <- read.csv(cohort_file, stringsAsFactors = FALSE)$species_name
  
  # Match on resolved species names (GTDB-preferred)
  in_cohort <- m_adapted %>%
    filter(Drug == "Veh") %>%              # one row per strain for counting
    distinct(StrainID, Species) %>%
    mutate(in_cohort = Species %in% cohort_sp)
  
  cat("\n  === Cohort matching summary ===\n")
  cat("  Total strains in Yuan dataset:          ",
      n_distinct(m_adapted$StrainID), "\n")
  cat("  Total species in Yuan dataset:          ",
      n_distinct(m_adapted$Species), "\n")
  cat("  Strains matched to cohort species:      ",
      sum(in_cohort$in_cohort), "\n")
  cat("  Species matched to cohort species:      ",
      n_distinct(in_cohort$Species[in_cohort$in_cohort]), "\n")
  
  # Flag cohort-matched strains in main data
  cohort_strains <- in_cohort$StrainID[in_cohort$in_cohort]
  m_adapted$in_cohort <- m_adapted$StrainID %in% cohort_strains
  
  # Write matching summary
  write.csv(in_cohort,
            file.path(TAB, "cohort_strain_matching.csv"),
            row.names = FALSE)
  cat("  -- Written to tables/cohort_strain_matching.csv\n")
} else {
  cat("  NOTE: cohort_species.csv not found in out_files/ — skipping cohort match.\n")
  cat("        Export from in-vitro script:\n")
  cat("        write.csv(cohort_species, here::here('out_files','cohort_species.csv'), row.names=FALSE)\n")
  m_adapted$in_cohort <- NA
}

## ── 0e. Filter to PFOS conditions (drop PFOSNO) ───────────────────────────
m0 <- m_adapted %>%
  filter(Drug %in% c("Veh", "PFOS-low", "PFOS-high")) %>%
  mutate(Drug = factor(Drug, levels = c("Veh", "PFOS-low", "PFOS-high")))

cat("\n  Rows retained (Veh + PFOS-low + PFOS-high):", nrow(m0), "\n")
cat("  Unique strains:", n_distinct(m0$StrainID), "\n\n")

# Full (pre-cohort-filter) copy for the Part 1 growth correlation (n = 66).
# The cohort filter below must NOT affect the validation correlation.
m0_full <- m0

## ── 0e.1 Filter to cohort-matched strains only ────────────────────────────
# NOTE: this filter also affects the Part 1 growth correlation below. If your
# manuscript reports n = 66 for the growth correlation, that number comes from
# running WITHOUT cohort_species.csv present. Keep that in mind when verifying.
if (exists("cohort_strains") && length(cohort_strains) > 0) {
  m0 <- m0 %>% filter(in_cohort == TRUE)
  cat("  Filtered to cohort-matched strains only\n")
  cat("  Unique strains after cohort filter:", n_distinct(m0$StrainID), "\n")
  cat("  Unique species after cohort filter:", n_distinct(m0$Species), "\n\n")
} else {
  cat("  NOTE: cohort_species.csv not found — using all strains\n\n")
}

## ── 0f. Check FVA file species coverage (on the _73 files used downstream) ──
check_fva_coverage <- function(fva_file, growth_species, label) {
  if (!file.exists(file.path(DAT_JOHANNES, fva_file))) {
    cat("  WARNING: FVA file not found:", fva_file, "\n")
    return(invisible(NULL))
  }
  fva <- read.csv(file.path(DAT_JOHANNES, fva_file), check.names = FALSE)
  fva_sp  <- unique(fva$Species)
  matched <- intersect(growth_species, fva_sp)
  only_growth <- setdiff(growth_species, fva_sp)
  only_fva    <- setdiff(fva_sp, growth_species)
  cat(sprintf("  %-50s: %d growth species, %d FVA species, %d matched\n",
              label, length(growth_species), length(fva_sp), length(matched)))
  if (length(only_growth) > 0)
    cat("    In growth data but NOT in FVA:", paste(only_growth, collapse = ", "), "\n")
  if (length(only_fva) > 0)
    cat("    In FVA but NOT in growth data:", paste(only_fva, collapse = ", "), "\n")
}

growth_sp_veh <- unique(m0$Species[m0$Drug == "Veh"])

cat("PART 0f — FVA file coverage checks (", FVA_SUFFIX, " files):\n", sep = "")
fva_files <- c(
  sprintf("FVA_Secretion_Capacity_HMO_Complex%s.csv",   FVA_SUFFIX),
  sprintf("FVA_Uptake_Capacity_HMO_Complex%s.csv",      FVA_SUFFIX),
  sprintf("FVA_Secretion_Capacity_HMO_BM_1month%s.csv", FVA_SUFFIX),
  sprintf("FVA_Uptake_Capacity_HMO_BM_1month%s.csv",    FVA_SUFFIX)
)
for (f in fva_files) {
  check_fva_coverage(f, growth_sp_veh,
                     sub("FVA_|_Capacity_HMO_|\\.csv", "", f))
}
cat("\n")

# ==========================================================================
# PART 1 — Observed vs predicted growth (3 media)
#   Reads FVA_SUFFIX_GROWTH files ("" = original 76-sp, preserves n=66).
# ==========================================================================
cat("PART 1 — Growth correlations (observed vs predicted):\n")

growth_df <- function(tag) {
  # "Complex" medium files have no "_BM_" infix; only the "1month"
  # (breast-milk) medium files do -- matches the naming used in Part 0f above.
  medium_tag <- if (tag == "Complex") tag else paste0("BM_", tag)
  fname <- sprintf("FVA_Secretion_Capacity_HMO_%s%s.csv", medium_tag, FVA_SUFFIX_GROWTH)
  sec <- read.csv(file.path(DAT_JOHANNES, fname), check.names = FALSE)
  
  gd <- m0_full[m0_full$Drug == "Veh", ]
  
  # FVA files use SILVA species names, but m0 uses GTDB-resolved names.
  # Merging on GTDB drops every strain where SILVA != GTDB (~43 of them),
  # collapsing n to ~29. Pick whichever name key pairs MORE species so the
  # validation spans all tested species (~66). Display labels stay GTDB.
  n_gtdb  <- length(intersect(unique(gd$Species),       unique(sec$Species)))
  n_silva <- length(intersect(unique(gd$SILVA_Species), unique(sec$Species)))
  key <- if (n_silva >= n_gtdb) "SILVA_Species" else "Species"
  cat(sprintf("    [%s] FVA species = %d | name matches: GTDB = %d, SILVA = %d  -> merging on %s\n",
              tag, n_distinct(sec$Species), n_gtdb, n_silva, key))
  
  gd$.mergekey <- gd[[key]]
  
  merge(gd[, c("Species", ".mergekey", "median_value")],
        sec[, c("Species", "EX_biomass(e)")],
        by.x = ".mergekey", by.y = "Species") %>%
    group_by(.mergekey) %>%                       # group on the SAME key used to merge
    summarise(
      Species = dplyr::first(Species),            # keep a GTDB name for the plot label
      x       = median(median_value, na.rm = TRUE),
      y       = median(`EX_biomass(e)`, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(!is.na(x) & !is.na(y))
}

growth_plot <- function(df, title, outfile) {
  r  <- cor(df$x, df$y)
  ct <- cor.test(df$x, df$y)
  n  <- nrow(df)
  p <- ggplot(df, aes(x, y)) +
    geom_smooth(method = "lm", formula = y ~ 0 + x, se = TRUE,
                color = "#2166AC", fill = "grey85", linewidth = 1) +
    geom_point(size = 2.4, alpha = 0.75, color = "grey20") +
    geom_text_repel(aes(label = Species), size = 2.6, max.overlaps = 14,
                    seed = 1, segment.size = 0.2, color = "grey35") +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.3,
             size = 5, fontface = "bold",
             label = sprintf("r = %.2f   (n = %d)", r, n)) +
    scale_x_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(title    = title,
         x        = "Observed growth (median AUC, Vehicle)",
         y        = expression("In silico maximal growth rate (" * h^{-1} * ")")) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 13.5),
      axis.title = element_text(size = 13, colour = "black")
    )
  ggsave(file.path(FIG, outfile), p,
         width = 8.5, height = 5.4, dpi = 300, bg = "white")
  cat(sprintf("  wrote %-52s r=%.2f p=%.3g n=%d\n",
              outfile, ct$estimate, ct$p.value, n))
  list(result = data.frame(medium      = title,
                           n_species   = n,
                           pearson_r   = round(r, 3),
                           p_value     = signif(ct$p.value, 3)),
       plot = p)
}

gC  <- growth_plot(growth_df("Complex"),
                   "Defined complex medium",
                   "MedianGrowth_PredictedGrowth_ComplexMedium.png")
gB1 <- growth_plot(growth_df("1month"),
                   "Breast milk (1 month postpartum) medium",
                   "MedianGrowth_PredictedGrowth_BM_1month_medium.png")

growth_res <- bind_rows(gC$result, gB1$result)

# 2-panel figure (A = complex, B = 1-month milk)
g2 <- gC$plot + gB1$plot +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))
ggsave(file.path(FIG, "Growth_Complex_vs_BM1month_2panel.png"),
       g2, width = 14, height = 5.4, dpi = 300, bg = "white")
cat("  wrote Growth_Complex_vs_BM1month_2panel.png\n\n")

# ==========================================================================
# PART 2 — Jaccard capability-similarity clustering heatmap  (uses _73 files)
# ==========================================================================
cat("PART 2 — Jaccard clustering heatmap:\n")


HMED <- "1month"

upH  <- read.csv(file.path(DAT_JOHANNES,
                           sprintf("FVA_Uptake_Capacity_HMO_BM_%s%s.csv",    HMED, FVA_SUFFIX)),
                 check.names = FALSE)
secH <- read.csv(file.path(DAT_JOHANNES,
                           sprintf("FVA_Secretion_Capacity_HMO_BM_%s%s.csv", HMED, FVA_SUFFIX)),
                 check.names = FALSE)

exU <- grep("^EX_", names(upH),  value = TRUE)
exS <- grep("^EX_", names(secH), value = TRUE)
exU <- exU[!grepl("EX_biomass", exU, ignore.case = TRUE)]
exS <- exS[!grepl("EX_biomass", exS, ignore.case = TRUE)]

U  <- as.data.frame(lapply(exU, function(c)
  as.integer((-upH[[c]]) > 0)))
names(U)  <- paste0("U_", exU); U$Species  <- upH$Species

Sx <- as.data.frame(lapply(exS, function(c)
  as.integer(secH[[c]] > 0)))
names(Sx) <- paste0("S_", exS); Sx$Species <- secH$Species

Fm <- merge(U, Sx, by = "Species")
rownames(Fm) <- Fm$Species; Fm$Species <- NULL

# Variable capabilities only (10–90% prevalence)
prev <- colMeans(Fm)
info <- names(prev)[prev >= 0.1 & prev <= 0.9]
cat(sprintf("  Variable capabilities (10-90%% prevalence): %d of %d\n",
            length(info), ncol(Fm)))

modal <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

sus <- m0 %>%
  filter(Drug == "PFOS-high") %>%
  group_by(Species) %>%
  summarise(
    inhib  = mean(perc_inhibition, na.rm = TRUE),
    phylum = modal(phylum),
    .groups = "drop"
  )

inhib_all <- setNames(sus$inhib,  sus$Species)
phy_all   <- setNames(sus$phylum, sus$Species)

# Restrict to species with PFOS-high data (no grey annotation cells)
sp <- intersect(rownames(Fm), sus$Species[!is.na(sus$inhib)])
cat(sprintf("  Species with PFOS-high data for annotation: %d\n", length(sp)))

# Warn if FVA has species not in growth data (join key check)
fva_only <- setdiff(rownames(Fm), unique(m0$Species))
if (length(fva_only) > 0) {
  cat("  WARNING: these FVA species have no match in growth data",
      "(check GTDB vs SILVA join key):\n")
  cat("  ", paste(fva_only, collapse = "\n   "), "\n")
}

Fi  <- Fm[sp, info]
inh <- inhib_all[sp]
phy <- phy_all[sp]
phy[is.na(phy)] <- "unclassified"

D  <- dist(Fi, method = "binary")
hc <- hclust(D, method = "ward.D2")
SS <- 1 - as.matrix(D)

mds1 <- cmdscale(D, k = 1)[, 1]
o0   <- order(mds1)
k    <- max(3L, round(length(o0) / 3))
if (mean(SS[tail(o0, k), tail(o0, k)]) >
    mean(SS[o0[1:k], o0[1:k]])) mds1 <- -mds1

basecol <- c(
  Bacteroidota     = "#B2182B", Bacillota_A    = "#F4A582",
  Actinomycetota   = "#1B7837", Bacillota      = "#67A9CF",
  Bacillota_C      = "#D6604D", Pseudomonadota = "#2166AC",
  Fusobacteriota   = "#999999", unclassified   = "grey80"
)
extra   <- setdiff(unique(phy), names(basecol))
phy_pal <- c(basecol, setNames(rep("grey60", length(extra)), extra))
phy_used <- phy_pal[intersect(names(phy_pal), unique(phy))]

inh_col <- colorRamp2(c(min(inh), max(inh)), c("white", "#B2182B"))

ha <- HeatmapAnnotation(
  `PFOS inhibition (%)` = inh, Phylum = phy,
  col = list(`PFOS inhibition (%)` = inh_col, Phylum = phy_used),
  annotation_name_gp = gpar(fontsize = 8.5, fontface = "bold")
)
hl <- rowAnnotation(
  `PFOS inhibition (%)` = inh, Phylum = phy,
  col = list(`PFOS inhibition (%)` = inh_col, Phylum = phy_used),
  show_annotation_name = FALSE, show_legend = FALSE
)
ht <- Heatmap(
  SS, name = "Capability\nsimilarity",
  col = colorRamp2(c(min(SS), 1), c("#f7fbff", "#08306b")),
  cluster_rows = hc, cluster_columns = hc, border = TRUE,
  row_dend_reorder    = mds1, column_dend_reorder = mds1,
  row_dend_width      = unit(16, "mm"),
  column_dend_height  = unit(16, "mm"),
  show_row_names      = TRUE, show_column_names = TRUE,
  row_names_gp    = gpar(fontsize = 9, fontface = "italic"),
  column_names_gp = gpar(fontsize = 9, fontface = "italic"),
  top_annotation = ha, left_annotation = hl
)

png(file.path(FIG, "FVA_Capability_Jaccard_Clustering.png"),
    width = 10, height = 9, units = "in", res = 300)
draw(ht,
     heatmap_legend_side      = "right",
     annotation_legend_side   = "right",
     merge_legend             = TRUE)
invisible(dev.off())
cat(sprintf("  wrote FVA_Capability_Jaccard_Clustering.png  (n=%d species, %d capabilities)\n\n",
            nrow(SS), length(info)))

# ==========================================================================
# PART 3 — Capability box plots with mixed-model interaction test (uses _73)
# ==========================================================================
cat("PART 3 — Box plots with mixed-model stats:\n")

BMED         <- "BM_1month"
custom_order <- c("Veh", "PFOS-low", "PFOS-high")
DRUG_LABELS  <- c(Veh = "Vehicle", "PFOS-low" = "Low PFOS",
                  "PFOS-high" = "High PFOS")

sig_stars <- function(p) {
  if (is.na(p))   return("ns")
  if (p < .001)   return("***")
  if (p < .01)    return("**")
  if (p < .05)    return("*")
  return("ns")
}

load_bin <- function(file, flip) {
  d <- read.csv(file.path(DAT_JOHANNES, file), check.names = FALSE)
  for (c in names(d)[-1])
    d[[c]] <- as.integer((if (flip) -d[[c]] else d[[c]]) > 0)
  d
}

box_plot <- function(file, flip, vmh_id, kind, title, outfile) {
  dmat <- load_bin(file, flip)
  
  # Join on resolved Species (GTDB-preferred) — same key as m0
  data <- merge(m0, dmat, by = "Species")
  
  # Check join quality
  n_joined  <- n_distinct(data$Species)
  n_fva     <- n_distinct(dmat$Species)
  n_growth  <- n_distinct(m0$Species)
  cat(sprintf("  [%s] FVA species=%d, growth species=%d, joined=%d\n",
              outfile, n_fva, n_growth, n_joined))
  
  data <- data[data$Drug %in% custom_order, ]
  data$Drug <- factor(data$Drug, levels = custom_order)
  
  col <- paste0("EX_", vmh_id, "(e)")
  if (!col %in% names(data)) stop("Column not found in FVA file: ", col)
  
  data$statf <- factor(data[[col]], levels = c(0, 1))
  data <- data[!is.na(data$statf), ]
  
  base2 <- if (kind == "Uptake") {
    c("Uptake: no", "Uptake: yes")
  } else {
    c("Secretion: no", "Secretion: yes")
  }
  
  nsp_total <- length(unique(data$Species))
  nsp_stat  <- tapply(data$Species, data$statf,
                      function(x) length(unique(x)))
  
  # Guard: need both capability groups to have >1 species for mixed model
  if (any(nsp_stat < 2)) {
    cat(sprintf("  WARNING: %s — too few species in one capability group; skipping model.\n",
                vmh_id))
    return(NULL)
  }
  
  labs2 <- paste0(base2, "  (n=", nsp_stat[c("0","1")], " sp.)")
  data$status <- factor(
    labs2[as.integer(as.character(data$statf)) + 1],
    levels = labs2
  )
  
  fit <- lmer(log(median_value) ~ statf * Drug + (1 | Species),
              data = data, REML = FALSE)
  
  ix <- grep("^statf1:Drug", names(fixef(fit)), value = TRUE)
  lh <- tryCatch(
    linearHypothesis(fit, paste(ix, "= 0"), test = "Chisq"),
    error = function(e) NULL
  )
  ip   <- if (is.null(lh)) NA else lh$`Pr(>Chisq)`[2]
  ichi <- if (is.null(lh)) NA else lh$Chisq[2]
  idf  <- if (is.null(lh)) NA else lh$Df[2]
  
  emm <- emmeans(fit, ~ Drug | statf, lmer.df = "asymptotic")
  ct  <- as.data.frame(contrast(emm, method = "trt.vs.ctrl", ref = 1))
  ct$drug   <- gsub("[()]", "", sub(" - Veh$", "", ct$contrast))
  ct$x      <- c("PFOS-low" = 1, "PFOS-high" = 2)[ct$drug]
  ct$status <- factor(
    labs2[as.integer(as.character(ct$statf)) + 1],
    levels = labs2
  )
  ct$lab <- vapply(ct$p.value, sig_stars, character(1))
  
  plotd <- data[data$Drug %in% c("PFOS-low", "PFOS-high"), ]
  plotd$Drug <- factor(as.character(plotd$Drug),
                       levels = c("PFOS-low", "PFOS-high"))
  
  prng <- max(plotd$perc, na.rm = TRUE)
  gmax <- tapply(plotd$perc, list(plotd$status, plotd$Drug),
                 max, na.rm = TRUE)
  ct$y <- mapply(function(s, d) gmax[s, d],
                 as.character(ct$status), ct$drug) + 0.04 * prng
  
  p <- ggplot(plotd, aes(Drug, perc, fill = Drug)) +
    geom_hline(yintercept = 100, linetype = "dashed",
               color = "grey40", linewidth = 0.5) +
    geom_boxplot(alpha = 0.6, color = "black", linewidth = 0.8,
                 outlier.shape = NA) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.5,
                size = 2, color = "black", stroke = 0.2) +
    geom_text(data = subset(ct, p.value < 0.05),
              aes(x = x, y = y, label = lab),
              inherit.aes = FALSE, size = 6,
              fontface = "bold", vjust = 0.3) +
    facet_wrap(~ status) +
    scale_fill_manual(
      values = c("PFOS-low" = "#21908C", "PFOS-high" = "#FDE725")
    ) +
    scale_x_discrete(labels = DRUG_LABELS) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(title = title, x = NULL, y = "Growth (% of Vehicle)") +
    theme_classic(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold", size = 13),
      strip.text       = element_text(size = 12, face = "bold"),
      strip.background = element_blank(),
      axis.text        = element_text(color = "black", face = "bold"),
      axis.title.y     = element_text(face = "bold"),
      legend.position  = "none"
    )
  
  ggsave(file.path(FIG, outfile), p,
         width = 7.6, height = 4.5, dpi = 300, bg = "white")
  cat(sprintf("  wrote %-38s interaction p = %.2e  (n=%d species)\n",
              outfile, ip, nsp_total))
  
  list(
    interaction = data.frame(
      metabolite              = vmh_id,
      capability              = tolower(kind),
      n_species               = nsp_total,
      n_without_capability    = as.integer(nsp_stat["0"]),
      n_with_capability       = as.integer(nsp_stat["1"]),
      interaction_chisq       = round(ichi, 2),
      interaction_df          = idf,
      interaction_p           = signif(ip, 3)
    ),
    contrasts = data.frame(
      metabolite              = vmh_id,
      capability              = tolower(kind),
      has_capability          = as.integer(as.character(ct$statf)),
      dose                    = ct$drug,
      estimate_vs_vehicle     = round(ct$estimate, 2),
      SE                      = round(ct$SE, 2),
      p_value                 = signif(ct$p.value, 3)
    ),
    plot = p
  )
}

box_res <- list(
  box_plot(sprintf("FVA_Uptake_Capacity_HMO_%s%s.csv",    BMED, FVA_SUFFIX),
           TRUE,  "ala_D",  "Uptake",
           "D-alanine uptake",
           "Uptake_Plot_ala_D.png"),
  box_plot(sprintf("FVA_Uptake_Capacity_HMO_%s%s.csv",    BMED, FVA_SUFFIX),
           TRUE,  "lcts",   "Uptake",
           "Lactose uptake",
           "Uptake_Plot_Lactose.png"),
  box_plot(sprintf("FVA_Secretion_Capacity_HMO_%s%s.csv", BMED, FVA_SUFFIX),
           FALSE, "ppa",    "Secretion",
           "Propionate secretion",
           "Secretion_Plot_PropionicAcid.png"),
  box_plot(sprintf("FVA_Secretion_Capacity_HMO_%s%s.csv", BMED, FVA_SUFFIX),
           FALSE, "isoval", "Secretion",
           "Isovalerate secretion",
           "Secretion_Plot_IsovalericAcid.png")
)

# only 1 species lacks lactose uptake capability among your 24 cohort-matched species.
# The guard condition requires at least 2 species in each group, so the lactose uptake model is skipped
# entirely because we can't fit a mixed-effects model with species random intercepts when one group has n=1

# Remove NULLs (skipped plots)
box_res <- Filter(Negate(is.null), box_res)

# 4-panel figure
shorts <- c("D-alanine uptake", "Lactose uptake",
            "Propionate secretion", "Isovalerate secretion")
panels <- lapply(seq_along(box_res), function(i)
  box_res[[i]]$plot +
    labs(title = shorts[i]) +
    theme(plot.title = element_text(face = "bold", size = 12)))

if (length(panels) == 4) {
  fig4 <- (panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]]) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 16))
  ggsave(file.path(FIG, "Boxplots_capability_4panel.png"),
         fig4, width = 13, height = 9.5, dpi = 300, bg = "white")
  cat("  wrote Boxplots_capability_4panel.png\n")
}

# ── Write statistical results ───────────────────────────────────────────────
write.xlsx(
  list(
    GrowthCorrelations   = growth_res,
    BoxplotInteractions  = bind_rows(lapply(box_res, `[[`, "interaction")),
    BoxplotDoseContrasts = bind_rows(lapply(box_res, `[[`, "contrasts"))
  ),
  file.path(TAB, "flux_model_results.xlsx")
)
cat("  wrote tables/flux_model_results.xlsx\n")

cat("\nDone — figures -> figures/  |  tables -> tables/\n")