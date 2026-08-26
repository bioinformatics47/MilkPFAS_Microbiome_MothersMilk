# header -----------------------------------------------------------------------
#
# TITLE:   12. Pathway-Metabolite Integration
#
# PURPOSE: Two complementary approaches linking PFAS-associated microbial gene
#          functional pathways with PFAS-associated fecal metabolites at 6m:
#
#   PART A — Spearman Correlation (unadjusted)
#     Mirrors Script 11 taxa-metabolite approach. Computes per-sample Spearman
#     correlations between CLR-transformed pathway abundances and log2 metabolite
#     intensities. Aggregated to pathway category x metabolite superclass summary
#     heatmaps (count of FDR-significant pairs; mean Spearman rho).
#
#   PART B — MWAS-Style Covariate-Adjusted Regression
#     For each PFAS-associated pathway x PFAS-associated metabolite pair, runs
#     covariate-adjusted linear regression:
#       metabolite ~ CLR_pathway + SES + gestational_age +
#                    mode_of_delivery + breastmilk_per_day + birthweight
#     Same covariate set and dummy variable coding as MWAS (Script 10).
#     Aggregated to pathway category x metabolite superclass summary heatmaps
#     (count of BH-significant pairs; mean beta direction).
#
# RATIONALE FOR TWO APPROACHES:
#   Spearman correlation is bivariate and unadjusted — provides a simple
#   description of co-occurrence patterns but may reflect confounding.
#   MWAS-style regression adjusts for the same covariates used throughout
#   the pipeline, providing covariate-adjusted pathway-metabolite associations
#   directly comparable to PFAS-taxa and PFAS-metabolite analyses.
#
# NOTE ON BIOLOGICAL INTERPRETABILITY:
#   HUMAnN pathway scores reflect microbial gene content (genetic potential)
#   rather than actual metabolic output. Most MetaCyc pathways represent
#   intracellular biosynthetic processes. The fecal metabolome integrates
#   host, dietary, and microbial contributions. Near-zero Spearman correlations
#   at 6m are consistent with these mechanistic limitations. Results from both
#   approaches are presented as supplementary exploratory analyses. The
#   taxa-metabolite correlation (Script 11) remains the primary integrative
#   analysis.
#
# INPUTS (from out_files/):
#   - PATHWAY_CLR_continuous_1m_6m.csv    (pathway PFAS regression results)
#   - analysis_df_cont_1m_6m.csv          (per-sample CLR pathway abundances)
#   - pathway_category_lookup.csv         (col_name -> Pathway -> Category)
#   - MWAS_continuous_c18_6m.csv          (MWAS C18 results, 6m)
#   - MWAS_continuous_hilic_6m.csv        (MWAS HILIC results, 6m)
#   - PFAS1m_c18_6m.csv                   (per-sample C18 intensities + covariates)
#   - PFAS1m_hilic_6m.csv                 (per-sample HILIC intensities + covariates)
#   - metabolite_cols_c18.rds             (C18 column list after QC)
#   - metabolite_cols_hilic.rds           (HILIC column list after QC)
#   - cname_superclass_lookup.csv         (manually annotated CNAME -> Super Class)
#
# OUTPUTS (to out_figures/ and out_files/):
#   Part A: pathway_metabolite_spearman_6m.csv
#           heatmap_pathway_metabolite_spearman_count_6m.pdf
#           heatmap_pathway_metabolite_spearman_rho_6m.pdf
#           pathway_metabolite_spearman_count_matrix_6m.csv
#           pathway_metabolite_spearman_rho_matrix_6m.csv
#
# DATE: May 2026
#
# set up -----------------------------------------------------------------------
rm(list = ls())
options(scipen = 0)

library(tidyverse)
library(dplyr)
library(here)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(stringr)
library(broom)

# Load combined label names for co-eluters
display_c18   <- read.csv(here::here("out_files", "c18_display_names.csv"),
                          stringsAsFactors = FALSE)
display_hilic <- read.csv(here::here("out_files", "hilic_display_names.csv"),
                          stringsAsFactors = FALSE)

# ── PFAS variable lists ────────────────────────────────────────────────────────
pfas_vars_continuous <- c("PFBS", "PFHxS", "PFNA", "PFOA", "PFOS")

# ── Covariates — same as MWAS (Script 10) ─────────────────────────────────────
# gestational_age_cat and mode_of_delivery_cat converted to dummies below
covariates <- c("SES_index_final", "baby_birthweight_kg",
                "breastmilk_per_day", "gest_Early", "gest_Late",
                "mode_of_delivery_bin")

# ── Ordered display levels ─────────────────────────────────────────────────────
pwy_cats_ordered <- c("Amino Acids", "Carbohydrates", "Lipids",
                      "Nucleotides", "Vitamins & Cofactors",
                      "Cell Wall", "SCFAs", "Other")

met_classes_ordered <- c("Nucleosides & Bases",
                         "Amino Acids & Derivatives",
                         "Sugars & Glycans",
                         "Aromatic & Indole Compounds",
                         "Lipids & Fatty Acids",
                         "Organic Acids & SCFAs",
                         "Steroids",
                         "Vitamins & Cofactors",
                         "Other")

# ── Create out_figures if needed ───────────────────────────────────────────────
if (!dir.exists(here::here("out_figures"))) {
  dir.create(here::here("out_figures"))
  cat("Created: out_figures/\n")
}

# ── Load QC column lists ───────────────────────────────────────────────────────
keep_c18   <- readRDS(here::here("out_files", "metabolite_cols_c18.rds"))
keep_hilic <- readRDS(here::here("out_files", "metabolite_cols_hilic.rds"))

# ── Load superclass lookup ─────────────────────────────────────────────────────
cname_superclass_lookup <- read.csv(here::here("out_files",
                                               "cname_superclass_lookup.csv"))

# ── Load pathway short name lookup ─────────────────────────────────────────
pwy_short_lookup <- read.csv(here::here("out_files", "pathway_short_names.csv"))
pwy_short_vec    <- setNames(pwy_short_lookup$Short_Name,
                             pwy_short_lookup$Pathway)
pwy_short_vec["P4-PWY superpathway of L-lysine, L-threonine and L-methionine biosynthesis I"] <- "Lys/Thr/Met biosyn. I"
cat("Pathway short names loaded:", length(pwy_short_vec), "\n")

# Both dot-format and hyphen-format entries needed to cover all lookup paths
pwy_short_vec["GLCMANNANAUT-PWY superpathway of N-acetylglucosamine, N-acetylmannosamine and N-acetylneuraminate degradation"] <- "GlcNAc/ManNAc/Neu5Ac Degrad."
pwy_short_vec["GLCMANNANAUT.PWY.superpathway.of.N.acetylglucosamine..N.acetylmannosamine.and.N.acetylneuraminate.degradation"] <- "GlcNAc/ManNAc/Neu5Ac Degrad."

# SHARED SETUP: LOAD AND PREPARE DATA USED BY BOTH PARTS
# ── Load PFAS-associated pathways at 6m (p < 0.05) ────────────────────────────
pwy_results_6m <- read.csv(here::here("out_files",
                                      "PATHWAY_CLR_continuous_1m_6m.csv")) %>%
  filter(!is.na(p_val), p_val < 0.05,
         predictor %in% pfas_vars_continuous) %>%
  mutate(
    plot_beta = dplyr::coalesce(as.numeric(Beta_IQR), as.numeric(betas)),
    Direction = ifelse(plot_beta > 0, "Positive", "Negative")
  )

sig_pathways_6m <- unique(pwy_results_6m$name)
cat("Significant pathways at 6m (p < 0.05):", length(sig_pathways_6m), "\n")
cat("Pathway categories:\n")
print(table(pwy_results_6m %>%
              dplyr::select(name, Category) %>%
              distinct() %>% pull(Category)))

# VERIFICATION: category names match pwy_cats_ordered
cats_in_data <- sort(unique(pwy_results_6m$Category))
missing_cats <- setdiff(cats_in_data, pwy_cats_ordered)
if (length(missing_cats) > 0) {
  cat("WARNING: categories in data NOT in pwy_cats_ordered:\n")
  print(missing_cats)
} else {
  cat("OK: all", length(cats_in_data), "pathway categories covered.\n")
}

# ── Load CLR pathway data (6m samples only) ───────────────────────────────────
pwy_col_lookup <- read.csv(here::here("out_files",
                                      "pathway_category_lookup.csv"))

sig_pwy_cols <- pwy_col_lookup %>%
  filter(Pathway %in% sig_pathways_6m) %>%
  dplyr::select(col_name, Pathway, Category)

pwy_data_all <- read.csv(here::here("out_files", "analysis_df_cont_1m_6m.csv"),
                         check.names = FALSE)

pwy_clr_6m <- pwy_data_all %>%
  filter(str_detect(merge_id_dyad, "-06$")) %>%
  dplyr::select(merge_id_dyad, any_of(sig_pwy_cols$col_name))

sig_pwy_cols <- sig_pwy_cols %>%
  filter(col_name %in% colnames(pwy_clr_6m))

cat("CLR pathway columns matched in data:", nrow(sig_pwy_cols), "\n")

# ── Load PFAS-associated metabolites at 6m (p < 0.05 MWAS) ───────────────────
mwas_c18_6m   <- read.csv(here::here("out_files", "MWAS_continuous_c18_6m.csv"))
mwas_hilic_6m <- read.csv(here::here("out_files", "MWAS_continuous_hilic_6m.csv"))

sig_c18 <- mwas_c18_6m %>%
  filter(P_value < 0.05, Metabolite %in% keep_c18,
         PFAS != "PFBS_pgmL") %>%
  pull(Metabolite) %>% unique()

sig_hilic <- mwas_hilic_6m %>%
  filter(P_value < 0.05, Metabolite %in% keep_hilic,
         PFAS != "PFBS_pgmL") %>%
  pull(Metabolite) %>% unique()

both       <- intersect(sig_c18, sig_hilic)
c18_only   <- setdiff(sig_c18,   both)
hilic_only <- setdiff(sig_hilic, both)

cat("Sig metabolites 6m — C18:", length(sig_c18),
    "| HILIC:", length(sig_hilic),
    "| Both:", length(both), "\n")

# ── Load per-sample metabolite intensities (6m) ───────────────────────────────
# These files also contain covariates needed for regression (Part B)
metab_c18_raw <- read.csv(here::here("out_files", "PFAS1m_c18_6m.csv"),
                          check.names = FALSE) %>%
  mutate(merge_id_dyad = str_replace(merge_id_dyad, "-01$", "-06"))

metab_hilic_raw <- read.csv(here::here("out_files", "PFAS1m_hilic_6m.csv"),
                            check.names = FALSE) %>%
  mutate(merge_id_dyad = str_replace(merge_id_dyad, "-01$", "-06"))

# Prepare covariate-adjusted data frame (same as MWAS Script 10)
prepare_data <- function(df) {
  df %>%
    mutate(
      gest_Early           = ifelse(gestational_age_cat == "Early",      1L, 0L),
      gest_Late            = ifelse(gestational_age_cat == "Late",       1L, 0L),
      mode_of_delivery_bin = ifelse(mode_of_delivery_cat == "C-Section", 1L, 0L),
      breastmilk_per_day   = as.numeric(scale(breastmilk_per_day)),
      SES_index_final      = as.numeric(scale(SES_index_final)),
      baby_birthweight_kg  = as.numeric(scale(baby_birthweight_kg))
    )
}

metab_c18_prep   <- prepare_data(metab_c18_raw)
metab_hilic_prep <- prepare_data(metab_hilic_raw)

# Build metabolite intensity subsets (with _C18/_HILIC suffixes for overlapping)
metab_c18_sub <- metab_c18_prep %>%
  dplyr::select(merge_id_dyad, all_of(covariates),
                any_of(c(c18_only, both)))

metab_hilic_sub <- metab_hilic_prep %>%
  dplyr::select(merge_id_dyad, any_of(c(hilic_only, both)))

if (length(both) > 0) {
  metab_c18_sub   <- metab_c18_sub %>%
    rename_with(~ paste0(., "_C18"),   any_of(both))
  metab_hilic_sub <- metab_hilic_sub %>%
    rename_with(~ paste0(., "_HILIC"), any_of(both))
}

all_met_labels <- c(c18_only, hilic_only,
                    paste0(both, "_C18"), paste0(both, "_HILIC"))

# ── Merge pathway CLR + metabolite data ───────────────────────────────────────
# C18 side carries the covariates; HILIC side adds only additional metabolites
merged_6m <- pwy_clr_6m %>%
  inner_join(metab_c18_sub,   by = "merge_id_dyad") %>%
  inner_join(metab_hilic_sub %>%
               dplyr::select(merge_id_dyad, any_of(all_met_labels)),
             by = "merge_id_dyad")

cat("Merged samples for analysis:", nrow(merged_6m), "\n")

pwy_cols_use <- sig_pwy_cols$col_name[sig_pwy_cols$col_name %in%
                                        colnames(merged_6m)]
met_cols_use <- all_met_labels[all_met_labels %in% colnames(merged_6m)]

cat("Pathway columns:", length(pwy_cols_use),
    "| Metabolite columns:", length(met_cols_use),
    "| Total pairs:", length(pwy_cols_use) * length(met_cols_use), "\n")

# ── Superclass annotation function ────────────────────────────────────────────
add_superclass <- function(df) {
  df %>%
    left_join(cname_superclass_lookup %>%
                dplyr::select(CNAME, Super_Class, Class),
              by = c("met_base" = "CNAME")) %>%
    mutate(
      Super_Class = replace_na(Super_Class, "Other"),
      Met_SuperClass = case_when(
        met_base %in% c("NICOTINAMIDE", "NICOTINAMIDE(B3)", "THIAMINE",
                        "PYRIDOXAL", "PYRIDOXAL(B6)", "PYRIDOXAMINE",
                        "PYRIDOXINE", "ALPHA-TOCOPHEROL", "ASCORBATE",
                        "D-PANTOTHENIC ACID_c18", "D-PANTOTHENIC ACID_hilic",
                        "PANTOTHENIC ACID(B5)", "4-PYRIDOXATE",
                        "DETHIOBIOTIN_c18", "DETHIOBIOTIN_hilic") ~ "Vitamins & Cofactors",
        met_base %in% c("FA4:0(BUTYRATE/ISOBUTYRATE)",
                        "FA5:0(VALERATE, ISOVALERATE, OTHERS)",
                        "SUCCINATE_c18", "SUCCINATE_hilic", "SUCCINIC ACID",
                        "PYRUVATE", "OXALIC ACID", "GLUTARATE",
                        "METHYLMALONIC ACID", "ETHYLMALONIC ACID",
                        "ITACONATE", "GLYCERATE", "GLYCERIC ACID",
                        "6-CARBOXYHEXANOATE", "TARTARIC ACID",
                        "D-SACCHARIC ACID", "GALACTARATE",
                        "ALPHA-KETOGLUTARIC ACID", "3-METHYGLUTARIC ACID",
                        "MONO-METHYL GLUTARATE") ~ "Organic Acids & SCFAs",
        str_detect(Super_Class, regex("nucleoside|nucleotide|purine|pyrimidine",
                                      ignore_case = TRUE)) ~ "Nucleosides & Bases",
        str_detect(Super_Class, regex(
          "amino acid|organic acid|carboxylic|peptide|nitrogen compound|organonitrogen",
          ignore_case = TRUE)) ~ "Amino Acids & Derivatives",
        str_detect(Super_Class, regex(
          "carbohydrate|sugar|glycan|saccharide|oxygen compound",
          ignore_case = TRUE)) ~ "Sugars & Glycans",
        str_detect(Super_Class, regex(
          "organoheterocyclic|alkaloid|indole|benzene|phenylpropan|polyketide|benzenoid",
          ignore_case = TRUE)) ~ "Aromatic & Indole Compounds",
        str_detect(Super_Class, regex(
          "lipid|fatty|sterol|sphingo|glycerophos|glycerolipid|prenol",
          ignore_case = TRUE)) ~ "Lipids & Fatty Acids",
        str_detect(Super_Class, regex("steroid|terpenoid",
                                      ignore_case = TRUE)) ~ "Steroids",
        TRUE ~ "Other"
      )
    )
}

# PART A: SPEARMAN CORRELATION (UNADJUSTED)-------------------------------------
# Mirrors Script 11 taxa-metabolite approach
spearman_results <- list()
idx <- 1L

for (pwy_col in pwy_cols_use) {
  for (met_col in met_cols_use) {
    x  <- merged_6m[[pwy_col]]
    y  <- merged_6m[[met_col]]
    ok <- sum(is.finite(x) & is.finite(y))
    if (ok < 3) {
      rho_val <- NA_real_; p_val <- NA_real_
    } else {
      ct      <- cor.test(x, y, method = "spearman", exact = FALSE)
      rho_val <- as.numeric(ct$estimate)
      p_val   <- as.numeric(ct$p.value)
    }
    spearman_results[[idx]] <- data.frame(
      pwy_col = pwy_col, met_label = met_col,
      rho = rho_val, p_val = p_val,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
}

spearman_df <- bind_rows(spearman_results) %>%
  mutate(FDR = p.adjust(p_val, method = "BH")) %>%
  left_join(sig_pwy_cols %>% dplyr::select(col_name, Pathway, Category),
            by = c("pwy_col" = "col_name")) %>%
  mutate(met_base = gsub("_C18$|_HILIC$", "", met_label))

cat("Significant pairs (FDR < 0.05):", sum(spearman_df$FDR < 0.05, na.rm = TRUE), "\n")
cat("Significant pairs (p < 0.05):",   sum(spearman_df$p_val < 0.05, na.rm = TRUE), "\n")
cat("Rho range:", round(min(spearman_df$rho, na.rm = TRUE), 3),
    "to", round(max(spearman_df$rho, na.rm = TRUE), 3), "\n")

# Add superclass annotation
spearman_df <- add_superclass(spearman_df)

cat("\nMetabolite superclass distribution (Spearman):\n")
print(table(spearman_df$Met_SuperClass))

n_other <- sum(spearman_df$Met_SuperClass == "Other", na.rm = TRUE)
if (n_other > 0) {
  cat("\nCNAMEs still in 'Other':\n")
  print(spearman_df %>% filter(Met_SuperClass == "Other") %>%
          distinct(met_base, Super_Class) %>% arrange(met_base))
} else {
  cat("OK: no metabolites in 'Other'.\n")
}

# Save full results
write.csv(spearman_df,
          here::here("out_files", "pathway_metabolite_spearman_6m.csv"),
          row.names = FALSE)


# PART A2: SPEARMAN CORRELATION — 1m-1m DATA
# Same approach as Part A, using PFAS-associated pathways and metabolites at 1m

# ── Load PFAS-associated pathways at 1m (p < 0.05) ────────────────────────────
pwy_results_1m <- read.csv(here::here("out_files",
                                      "PATHWAY_CLR_continuous_1m_1m.csv")) %>%
  filter(!is.na(p_val), p_val < 0.05,
         predictor %in% pfas_vars_continuous) %>%
  mutate(
    plot_beta = dplyr::coalesce(as.numeric(Beta_IQR), as.numeric(betas)),
    Direction = ifelse(plot_beta > 0, "Positive", "Negative")
  )

sig_pathways_1m <- unique(pwy_results_1m$name)
cat("Significant pathways at 1m (p < 0.05):", length(sig_pathways_1m), "\n")

# ── Pathway columns for 1m ─────────────────────────────────────────────────────
sig_pwy_cols_1m <- pwy_col_lookup %>%
  filter(Pathway %in% sig_pathways_1m) %>%
  dplyr::select(col_name, Pathway, Category)

pwy_data_1m <- read.csv(here::here("out_files", "analysis_df_cont_1m_1m.csv"),
                        check.names = FALSE)

pwy_clr_1m <- pwy_data_1m %>%
  dplyr::select(merge_id_dyad, any_of(sig_pwy_cols_1m$col_name))

sig_pwy_cols_1m <- sig_pwy_cols_1m %>%
  filter(col_name %in% colnames(pwy_clr_1m))

cat("CLR pathway columns matched in 1m data:", nrow(sig_pwy_cols_1m), "\n")

# ── Load PFAS-associated metabolites at 1m (p < 0.05 MWAS) ───────────────────
mwas_c18_1m   <- read.csv(here::here("out_files", "MWAS_continuous_c18_1m.csv"))
mwas_hilic_1m <- read.csv(here::here("out_files", "MWAS_continuous_hilic_1m.csv"))

sig_c18_1m <- mwas_c18_1m %>%
  filter(P_value < 0.05, Metabolite %in% keep_c18,
         PFAS != "PFBS_pgmL") %>%
  pull(Metabolite) %>% unique()

sig_hilic_1m <- mwas_hilic_1m %>%
  filter(P_value < 0.05, Metabolite %in% keep_hilic,
         PFAS != "PFBS_pgmL") %>%
  pull(Metabolite) %>% unique()

both_1m       <- intersect(sig_c18_1m, sig_hilic_1m)
c18_only_1m   <- setdiff(sig_c18_1m, both_1m)
hilic_only_1m <- setdiff(sig_hilic_1m, both_1m)

cat("Sig metabolites 1m — C18:", length(sig_c18_1m),
    "| HILIC:", length(sig_hilic_1m),
    "| Both:", length(both_1m), "\n")

# ── Load per-sample metabolite intensities (1m) ───────────────────────────────
metab_c18_1m_raw   <- read.csv(here::here("out_files", "PFAS1m_c18_1m.csv"),
                                check.names = FALSE)
metab_hilic_1m_raw <- read.csv(here::here("out_files", "PFAS1m_hilic_1m.csv"),
                                check.names = FALSE)

metab_c18_1m_prep   <- prepare_data(metab_c18_1m_raw)
metab_hilic_1m_prep <- prepare_data(metab_hilic_1m_raw)

metab_c18_1m_sub <- metab_c18_1m_prep %>%
  dplyr::select(merge_id_dyad, all_of(intersect(covariates, colnames(metab_c18_1m_prep))),
                any_of(c(c18_only_1m, both_1m)))

metab_hilic_1m_sub <- metab_hilic_1m_prep %>%
  dplyr::select(merge_id_dyad, any_of(c(hilic_only_1m, both_1m)))

if (length(both_1m) > 0) {
  metab_c18_1m_sub   <- metab_c18_1m_sub %>%
    rename_with(~ paste0(., "_C18"),   any_of(both_1m))
  metab_hilic_1m_sub <- metab_hilic_1m_sub %>%
    rename_with(~ paste0(., "_HILIC"), any_of(both_1m))
}

all_met_labels_1m <- c(c18_only_1m, hilic_only_1m,
                       paste0(both_1m, "_C18"), paste0(both_1m, "_HILIC"))

# ── Merge pathway CLR + metabolite data (1m) ──────────────────────────────────
merged_1m <- pwy_clr_1m %>%
  inner_join(metab_c18_1m_sub, by = "merge_id_dyad") %>%
  inner_join(metab_hilic_1m_sub %>%
               dplyr::select(merge_id_dyad, any_of(all_met_labels_1m)),
             by = "merge_id_dyad")

cat("Merged samples for 1m analysis:", nrow(merged_1m), "\n")

pwy_cols_use_1m <- sig_pwy_cols_1m$col_name[sig_pwy_cols_1m$col_name %in%
                                              colnames(merged_1m)]
met_cols_use_1m <- all_met_labels_1m[all_met_labels_1m %in% colnames(merged_1m)]

cat("Pathway columns:", length(pwy_cols_use_1m),
    "| Metabolite columns:", length(met_cols_use_1m),
    "| Total pairs:", length(pwy_cols_use_1m) * length(met_cols_use_1m), "\n")

# ── Run Spearman correlations (1m) ─────────────────────────────────────────────
spearman_results_1m <- list()
idx <- 1L

for (pwy_col in pwy_cols_use_1m) {
  for (met_col in met_cols_use_1m) {
    x  <- merged_1m[[pwy_col]]
    y  <- merged_1m[[met_col]]
    ok <- sum(is.finite(x) & is.finite(y))
    if (ok < 3) {
      rho_val <- NA_real_; p_val <- NA_real_
    } else {
      ct      <- cor.test(x, y, method = "spearman", exact = FALSE)
      rho_val <- as.numeric(ct$estimate)
      p_val   <- as.numeric(ct$p.value)
    }
    spearman_results_1m[[idx]] <- data.frame(
      pwy_col = pwy_col, met_label = met_col,
      rho = rho_val, p_val = p_val,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
}

spearman_df_1m <- bind_rows(spearman_results_1m) %>%
  mutate(FDR = p.adjust(p_val, method = "BH")) %>%
  left_join(sig_pwy_cols_1m %>% dplyr::select(col_name, Pathway, Category),
            by = c("pwy_col" = "col_name")) %>%
  mutate(met_base = gsub("_C18$|_HILIC$", "", met_label))

spearman_df_1m <- add_superclass(spearman_df_1m)

cat("1m Significant pairs (FDR < 0.05):", sum(spearman_df_1m$FDR < 0.05, na.rm = TRUE), "\n")
cat("1m Significant pairs (p < 0.05):",   sum(spearman_df_1m$p_val < 0.05, na.rm = TRUE), "\n")
cat("1m Rho range:", round(min(spearman_df_1m$rho, na.rm = TRUE), 3),
    "to", round(max(spearman_df_1m$rho, na.rm = TRUE), 3), "\n")

# Save full results
write.csv(spearman_df_1m,
          here::here("out_files", "pathway_metabolite_spearman_1m.csv"),
          row.names = FALSE)

# ==============================================================================
# PART C: INDIVIDUAL PATHWAY x METABOLITE HEATMAPS
# Layout mirrors Script 11 taxa-metabolite correlation:
#   - Pathways on x-axis (columns), metabolites on y-axis (rows)
#   - Top annotation: pathway category bar + PFAS direction bars per pathway
#   - Left annotation: PFAS direction bars per metabolite
#   - purple3 = positively PFAS-associated, orange3 = negatively PFAS-associated
#   - Asterisk (*) = FDR < 0.05
#   - Iterative mutual filter applied
#   - Within each pathway category: negative pathways first then positive
# ==============================================================================

# ── build_pwy_met_heatmap() ───────────────────────────────────────────────────
build_pwy_met_heatmap <- function(result_df,
                                  effect_col   = "rho",
                                  mwas_c18     = mwas_c18_6m,
                                  mwas_hilic   = mwas_hilic_6m,
                                  pwy_results  = pwy_results_6m,
                                  fdr_thresh   = 0.05,
                                  min_assoc    = 5,
                                  cellwidth    = 10,
                                  cellheight   = 10,
                                  fontsize_col = 7,
                                  fontsize_row = 9,
                                  bar_left     = 0.4,
                                  bar_y        = 0.95) {
  
  # ── Filter to FDR-significant pairs ───────────────────────────────────────
  sig_df <- result_df %>%
    filter(!is.na(FDR), FDR < fdr_thresh,
           Met_SuperClass != "Other",
           Category != "Other")
  
  if (nrow(sig_df) == 0) {
    cat("No significant pairs at FDR <", fdr_thresh, "- skipping.\n")
    return(NULL)
  }
  
  sig_pathways_plot <- unique(sig_df$Pathway)
  sig_mets_plot     <- unique(sig_df$met_label)
  
  cat("Significant pairs before filter:", nrow(sig_df), "\n")
  cat("Pathways:", length(sig_pathways_plot),
      "| Metabolites:", length(sig_mets_plot), "\n")
  
  # ── Iterative mutual filter ────────────────────────────────────────────────
  repeat {
    prev_pwy <- sig_pathways_plot
    prev_met <- sig_mets_plot
    
    sig_mets_plot <- sig_df %>%
      filter(Pathway %in% sig_pathways_plot) %>%
      count(met_label) %>%
      filter(n >= min_assoc) %>%
      pull(met_label)
    
    sig_pathways_plot <- sig_df %>%
      filter(met_label %in% sig_mets_plot) %>%
      count(Pathway) %>%
      filter(n >= min_assoc) %>%
      pull(Pathway)
    
    if (setequal(sig_pathways_plot, prev_pwy) &&
        setequal(sig_mets_plot,     prev_met)) break
  }
  
  cat("After iterative mutual filter (>=", min_assoc, "associations):",
      length(sig_pathways_plot), "pathways x", length(sig_mets_plot), "metabolites\n")
  
  # Second iterative filter: require direct PFAS association for metabolites
  # Consistent with Script 11 taxa-metabolite pipeline
  mwas_combined_local <- bind_rows(
    mwas_c18   %>% dplyr::select(Metabolite, PFAS, Estimate, P_value),
    mwas_hilic %>% dplyr::select(Metabolite, PFAS, Estimate, P_value)
  ) %>%
    filter(P_value < 0.05, !is.na(Estimate)) %>%
    mutate(PFAS_clean = gsub("_pgmL$", "", PFAS)) %>%
    filter(PFAS_clean %in% pfas_vars_continuous) %>%
    pull(Metabolite) %>%
    unique()
  
  repeat {
    prev_pwy <- sig_pathways_plot
    prev_met <- sig_mets_plot
    
    # Restrict to metabolites with direct PFAS association
    sig_mets_plot <- sig_mets_plot[sig_mets_plot %in% mwas_combined_local]
    
    # Re-apply metabolite threshold
    sig_mets_plot <- sig_df %>%
      filter(Pathway %in% sig_pathways_plot,
             met_label %in% sig_mets_plot) %>%
      count(met_label) %>%
      filter(n >= min_assoc) %>%
      pull(met_label)
    
    # Keep only those with direct PFAS association
    sig_mets_plot <- sig_mets_plot[sig_mets_plot %in% mwas_combined_local]
    
    # Re-apply pathway threshold
    sig_pathways_plot <- sig_df %>%
      filter(met_label %in% sig_mets_plot) %>%
      count(Pathway) %>%
      filter(n >= min_assoc) %>%
      pull(Pathway)
    
    if (setequal(sig_pathways_plot, prev_pwy) &&
        setequal(sig_mets_plot,     prev_met)) break
  }
  
  cat("After final filter (direct PFAS + thresholds):",
      length(sig_pathways_plot), "pathways x", length(sig_mets_plot), "metabolites\n")
  
  if (length(sig_pathways_plot) == 0 || length(sig_mets_plot) == 0) {
    cat("Nothing survived final filter - skipping.\n")
    return(NULL)
  }
  
  # ── Build effect matrix: rows = metabolites, cols = pathways ──────────────
  plot_df         <- result_df %>%
    filter(Pathway %in% sig_pathways_plot,
           met_label %in% sig_mets_plot)
  plot_df$eff_val <- plot_df[[effect_col]]
  
  eff_wide <- plot_df %>%
    group_by(met_label, Pathway) %>%
    summarise(eff = mean(eff_val, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Pathway, values_from = eff) %>%
    column_to_rownames("met_label") %>%
    as.matrix()
  
  # ── PFAS direction for pathways (top annotation) ───────────────────────────
  pwy_dir <- pwy_results %>%
    filter(name %in% sig_pathways_plot) %>%
    mutate(
      PFAS_clean = predictor,
      direction  = ifelse(plot_beta > 0, "Positive", "Negative")
    ) %>%
    dplyr::select(Pathway = name, PFAS_clean, direction) %>%
    filter(PFAS_clean %in% pfas_vars_continuous) %>%
    distinct()
  
  annot_col <- pwy_dir %>%
    pivot_wider(names_from  = PFAS_clean,
                values_from = direction,
                values_fill = NA_character_) %>%
    column_to_rownames("Pathway")
  
  for (p in pfas_vars_continuous) {
    if (!p %in% colnames(annot_col)) annot_col[[p]] <- NA_character_
  }
  annot_col <- annot_col[, colSums(!is.na(annot_col)) > 0, drop = FALSE]
  
  # Pathway category annotation
  pwy_cat_annot <- result_df %>%
    distinct(Pathway, Category) %>%
    filter(Pathway %in% sig_pathways_plot) %>%
    distinct(Pathway, .keep_all = TRUE) %>%
    column_to_rownames("Pathway") %>%
    dplyr::rename(`Pathway Category` = Category)
  
  # Dominant direction per pathway for within-category ordering
  pwy_dom_dir <- pwy_dir %>%
    count(Pathway, direction) %>%
    group_by(Pathway) %>%
    slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::select(Pathway, dominant_dir = direction)
  
  # Combine: category bar first, then PFAS direction bars
  annot_col_combined <- cbind(
    pwy_cat_annot[rownames(annot_col), , drop = FALSE],
    annot_col
  )
  
  # ── PFAS direction for metabolites (left annotation) ───────────────────────
  mwas_dir_met <- bind_rows(
    mwas_c18   %>% dplyr::select(Metabolite, PFAS, Estimate),
    mwas_hilic %>% dplyr::select(Metabolite, PFAS, Estimate)
  ) %>%
    filter(!is.na(Estimate)) %>%
    mutate(
      direction  = ifelse(Estimate > 0, "Positive", "Negative"),
      PFAS_clean = gsub("_pgmL$", "", PFAS)
    ) %>%
    filter(PFAS_clean %in% pfas_vars_continuous) %>%
    dplyr::select(Metabolite, PFAS_clean, direction) %>%
    distinct()
  
  met_base_map <- result_df %>%
    distinct(met_label, met_base) %>%
    filter(met_label %in% rownames(eff_wide))
  
  annot_row <- met_base_map %>%
    left_join(mwas_dir_met, by = c("met_base" = "Metabolite"),
              relationship = "many-to-many") %>%
    filter(!is.na(PFAS_clean)) %>%
    pivot_wider(names_from  = PFAS_clean,
                values_from = direction,
                values_fill = NA_character_) %>%
    dplyr::select(-met_base) %>%
    distinct(met_label, .keep_all = TRUE) %>%
    column_to_rownames("met_label")
  
  for (p in pfas_vars_continuous) {
    if (!p %in% colnames(annot_row)) annot_row[[p]] <- NA_character_
  }
  annot_row <- annot_row[, pfas_vars_continuous, drop = FALSE]
  annot_row <- annot_row[, colSums(!is.na(annot_row)) > 0, drop = FALSE]
  
  cat("PFAS in pathway annotation:", paste(colnames(annot_col), collapse = ", "), "\n")
  cat("PFAS in metabolite annotation:", paste(colnames(annot_row), collapse = ", "), "\n")
  
  # ── Pathway column order: category → negative first → alphabetical ─────────
  pwy_cat_order_df <- result_df %>%
    distinct(Pathway, Category) %>%
    filter(Pathway %in% colnames(eff_wide)) %>%
    left_join(pwy_dom_dir, by = "Pathway") %>%
    mutate(dominant_dir = replace_na(dominant_dir, "Positive"))
  
  cats_present <- pwy_cats_ordered[pwy_cats_ordered %in%
                                     unique(pwy_cat_order_df$Category) &
                                     pwy_cats_ordered != "Other"]
  
  pwy_ordered <- pwy_cat_order_df %>%
    mutate(
      Category  = factor(Category, levels = cats_present),
      dir_order = ifelse(dominant_dir == "Negative", 1, 2)
    ) %>%
    arrange(Category, dir_order, Pathway) %>%
    pull(Pathway)
  pwy_ordered <- pwy_ordered[pwy_ordered %in% colnames(eff_wide)]
  
  # ── Metabolite row order: negative first → by superclass → alphabetical ────
  # Dominant direction — most significant PFAS association per metabolite
  met_dom_dir <- bind_rows(
    mwas_c18   %>% dplyr::select(Metabolite, PFAS, Estimate, P_value),
    mwas_hilic %>% dplyr::select(Metabolite, PFAS, Estimate, P_value)
  ) %>%
    filter(!is.na(Estimate), !is.na(P_value)) %>%
    mutate(PFAS_clean = gsub("_pgmL$", "", PFAS)) %>%
    filter(PFAS_clean %in% pfas_vars_continuous) %>%
    filter(Metabolite %in% (met_base_map %>% pull(met_base))) %>%
    group_by(Metabolite) %>%
    slice_min(order_by = P_value, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(dominant_dir = ifelse(Estimate > 0, "Positive", "Negative")) %>%
    # Manual overrides for metabolites with mixed directions where
    # biological context supports negative classification
    mutate(dominant_dir = case_when(
      Metabolite == "CHOLINE"         ~ "Negative",
      Metabolite == "L-ALLOTHREONINE" ~ "Negative",
      Metabolite == "L-THREONINE"     ~ "Negative",
      TRUE                            ~ dominant_dir
    )) %>%
    left_join(met_base_map %>% rename(Metabolite = met_base),
              by = "Metabolite") %>%
    dplyr::select(met_label, dominant_dir)
  
  met_class_order <- result_df %>%
    distinct(met_label, met_base, Met_SuperClass) %>%
    distinct(met_label, .keep_all = TRUE)
  
  met_max_eff <- plot_df %>%
    group_by(met_label) %>%
    summarise(max_eff = max(abs(eff_val), na.rm = TRUE), .groups = "drop")
  
  met_ordered <- met_max_eff %>%
    left_join(met_dom_dir, by = "met_label") %>%
    left_join(met_class_order %>%
                dplyr::select(met_label, met_base, Met_SuperClass),
              by = "met_label") %>%
    mutate(
      dominant_dir   = replace_na(dominant_dir, "Positive"),
      Met_SuperClass = replace_na(Met_SuperClass, "Other"),
      met_base       = ifelse(is.na(met_base), met_label, met_base),
      dir_order      = ifelse(dominant_dir == "Negative", 1, 2),
      class_order    = match(Met_SuperClass, met_classes_ordered),
      class_order    = ifelse(is.na(class_order), 99L, class_order)
    ) %>%
    arrange(dir_order, class_order, met_base) %>%
    pull(met_label)
  met_ordered <- met_ordered[met_ordered %in% rownames(eff_wide)]
  
  # Reorder matrix
  eff_wide <- eff_wide[met_ordered, pwy_ordered, drop = FALSE]
  
  # Align annotations
  annot_col_combined <- annot_col_combined[colnames(eff_wide), , drop = FALSE]
  annot_row          <- annot_row[rownames(eff_wide), , drop = FALSE]
  
  # ── Significance labels ────────────────────────────────────────────────────
  sig_wide <- result_df %>%
    filter(met_label %in% rownames(eff_wide),
           Pathway   %in% colnames(eff_wide),
           !is.na(FDR)) %>%
    group_by(met_label, Pathway) %>%
    summarise(min_fdr = min(FDR, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Pathway, values_from = min_fdr) %>%
    column_to_rownames("met_label") %>%
    as.matrix()
  
  for (pwy in colnames(eff_wide)) {
    if (!pwy %in% colnames(sig_wide)) sig_wide <- cbind(sig_wide, NA_real_)
  }
  sig_wide   <- sig_wide[rownames(eff_wide), colnames(eff_wide), drop = FALSE]
  sig_labels <- matrix(
    ifelse(!is.na(sig_wide) & sig_wide < fdr_thresh, "*", ""),
    nrow = nrow(sig_wide), ncol = ncol(sig_wide),
    dimnames = dimnames(sig_wide)
  )
  
  # ── Strip met_label suffixes and apply display names (LAST step) ──────────
  # Manual display name overrides — same as Script 11 taxa-metabolite heatmap
  manual_display_overrides <- c(
    "6-DEOXY-L-GALACTOSE (FUCOSE)"              = "Fucose; L-Rhamnose",
    "BETA-ALANINE"                               = "Beta-Alanine; Sarcosine; D-Alanine",
    "SUCCINATE"                                  = "Succinate; Methylmalonic Acid",
    "L-THREONINE"                                = "Threonine; L-Allothreonine",
    "5-OXO-D-PROLINE"                            = "5-Oxoproline",
    "GLUTARATE"                                  = "Glutarate; Ethylmalonic Acid",
    "L-GLUTAMIC ACID"                            = "Glutamic Acid; N-Methyl-D-Aspartate",
    "(RS)-MEVALONIC ACID LITHIUM SALT"           = "Mevalonic Acid",
    "MANDELIC ACID"                              = "Mandelic Acid; Hydroxyphenylacetate",
    "ALPHA-AMINOADIPATE"                         = "Aminoadipate; N-Methyl-Glutamate",
    "METHYL VANILLATE"                           = "Methyl Vanillate; Homovanillate",
    "XYLITOL"                                    = "Xylitol; Ribitol",
    "D-FRUCTOSE"                                 = "Hexose",
    "PETROSELINIC ACID"                          = "Petroselinic Acid; Elaidic Acid",
    "HEPTADECANOATE"                             = "Heptadecanoate",
    "LAURIC ACID"                                = "Lauric Acid",
    "NICOTINAMIDE"                               = "Nicotinamide",
    "L-ASPARAGINE"                               = "Asparagine",
    "ASPARAGINE"                                 = "Asparagine",
    "2-QUINOLINECARBOXYLIC ACID"                 = "Quinolinecarboxylic Acid",
    "CYTIDINE 2',3'-CYCLIC MONO-PHOS-PHATE"     = "Cytidine",
    "N-ALPHA-ACETYL-L-ASPARAGINE"               = "Na-Acetyl-Asparagine",
    "NALPHA-ACETYL-L-LYSINE"                     = "Na-Acetyl-Lysine",
    "FA5:0(VALERATE, ISOVALERATE, OTHERS)"       = "Valerate/Isovalerate",
    "FA4:0(BUTYRATE/ISOBUTYRATE)"               = "Butyrate/Isobutyrate",
    "MG(14:0/0:0/0:0)"                          = "MG(14:0)",
    "4-AMINOBUTANOATE"                           = "GABA; Aminobutyrate",
    "2-METHYLMALEATE"                            = "Methylmaleate; Itaconate",
    "D-ASPARTATE"                               = "Aspartate",
    "TYRAMINE"                                   = "Tyramine; Phenylethanolamine",
    "4-HYDROXY-L-PHENYLGLYCINE"                 = "Pyridoxal; Hydroxyphenylglycine",
    "PYRIDOXINE"                                 = "Pyridoxine; Noradrenaline",
    "NORMETANEPHRINE"                            = "Normetanephrine; Epinephrine",
    "2-METHYLHIPPURATE"                          = "Methylhippurate",
    "D-(+)-GLUCOSAMINE"                          = "Glucosamine; Mannosamine",
    "D-GALACTOSE"                               = "Hexose (HILIC)",
    "GAMMA-LINOLENIC ACID"                       = "Linolenic Acid"
  )
  
  basic_clean <- function(nm) {
    nm %>%
      gsub("^L-", "", .) %>%
      gsub("^D-", "", .) %>%
      gsub("^DL-", "", .) %>%
      gsub("\\(B[0-9]+\\)", "", .) %>%
      stringr::str_to_title(.) %>%
      gsub("N-Acetyl([- ])", "N-acetyl\\1", .) %>%
      trimws(.)
  }
  
  # Build display name map from current rownames
  row_display_map <- setNames(rownames(eff_wide), rownames(eff_wide))
  
  for (nm in rownames(eff_wide)) {
    if (grepl("_[Cc]18$", nm)) {
      base   <- sub("_[Cc]18$", "", nm)
      suffix <- " (C18)"
    } else if (grepl("_[Hh][Ii][Ll][Ii][Cc]$", nm)) {
      base   <- sub("_[Hh][Ii][Ll][Ii][Cc]$", "", nm)
      suffix <- " (HILIC)"
    } else {
      base   <- nm
      suffix <- ""
    }
    
    if (base %in% names(manual_display_overrides)) {
      row_display_map[nm] <- paste0(manual_display_overrides[base], suffix)
    } else {
      match_c18   <- display_c18[display_c18$CNAME == base, "display_name"]
      match_hilic <- display_hilic[display_hilic$CNAME == base, "display_name"]
      if (length(match_c18) > 0 && match_c18[1] != base) {
        row_display_map[nm] <- paste0(basic_clean(match_c18[1]), suffix)
      } else if (length(match_hilic) > 0 && match_hilic[1] != base) {
        row_display_map[nm] <- paste0(basic_clean(match_hilic[1]), suffix)
      } else {
        row_display_map[nm] <- paste0(basic_clean(base), suffix)
      }
    }
  }
  
  # Deduplicate after applying display names
  clean_row_names <- row_display_map[rownames(eff_wide)]
  dup_row_mask    <- duplicated(clean_row_names)
  if (any(dup_row_mask)) {
    eff_wide        <- eff_wide[  !dup_row_mask, , drop = FALSE]
    sig_labels      <- sig_labels[!dup_row_mask, , drop = FALSE]
    annot_row       <- annot_row[ !dup_row_mask, , drop = FALSE]
    clean_row_names <- clean_row_names[!dup_row_mask]
  }
  
  rownames(eff_wide)   <- clean_row_names
  rownames(sig_labels) <- clean_row_names
  rownames(annot_row)  <- clean_row_names
  
  cat("  Display names applied:", sum(row_display_map != names(row_display_map)), "\n")
  
  # ── Apply short pathway names for display ─────────────────────────────────
  short_col <- ifelse(colnames(eff_wide) %in% names(pwy_short_vec),
                      pwy_short_vec[colnames(eff_wide)],
                      colnames(eff_wide))
  n_unmatched <- sum(!colnames(eff_wide) %in% names(pwy_short_vec))
  if (n_unmatched > 0) {
    cat("WARNING:", n_unmatched, "pathways not in short name lookup:\n")
    print(colnames(eff_wide)[!colnames(eff_wide) %in% names(pwy_short_vec)])
  }
  colnames(eff_wide)         <- short_col
  colnames(sig_labels)       <- short_col
  rownames(annot_col_combined) <- short_col
  
  # ── Annotation colors ──────────────────────────────────────────────────────
  dir_colors <- c("Positive" = "purple3", "Negative" = "orange3")
  
  cats_present_final <- cats_present[cats_present %in%
                                       unique(pwy_cat_annot$`Pathway Category`)]
  cat_colors <- setNames(
    colorRampPalette(brewer.pal(8, "Set2"))(length(cats_present_final)),
    cats_present_final
  )
  
  all_pfas_annot <- unique(c(
    colnames(annot_col_combined)[colnames(annot_col_combined) != "Pathway Category"],
    colnames(annot_row)
  ))
  pfas_color_list <- setNames(
    lapply(all_pfas_annot, function(p) dir_colors),
    all_pfas_annot
  )
  
  annot_colors <- c(
    list(`Pathway Category` = cat_colors),
    pfas_color_list
  )
  
  # ── Build pheatmap silently ────────────────────────────────────────────────
  pal       <- colorRampPalette(brewer.pal(11, "RdBu"))(100)
  eff_range <- max(abs(eff_wide), na.rm = TRUE)
  if (!is.finite(eff_range) || eff_range == 0) eff_range <- 0.7
  breaks    <- seq(-eff_range, eff_range, length.out = 101)
  
  p <- pheatmap::pheatmap(
    mat                  = eff_wide,
    color                = pal,
    breaks               = breaks,
    display_numbers      = sig_labels,
    fontsize_number      = 7,
    number_color         = "black",
    annotation_col       = annot_col_combined,
    annotation_row       = annot_row,
    annotation_colors    = annot_colors,
    annotation_names_col = TRUE,
    annotation_names_row = TRUE,
    annotation_legend    = FALSE,
    legend               = FALSE,
    cluster_rows         = FALSE,
    cluster_cols         = FALSE,
    show_rownames        = TRUE,
    show_colnames        = TRUE,
    cellwidth            = cellwidth,
    cellheight           = cellheight,
    fontsize_row         = fontsize_row,
    fontsize_col         = fontsize_col,
    angle_col            = 45,
    border_color         = NA,
    na_col               = "white",
    silent               = TRUE
  )
  
  # ── Figure sizing ──────────────────────────────────────────────────────────
  n_rows        <- nrow(eff_wide)
  n_cols        <- ncol(eff_wide)
  max_col_chars <- max(nchar(colnames(eff_wide)))
  max_row_chars <- max(nchar(rownames(eff_wide)))
  
  col_label_h   <- max_col_chars * fontsize_col * 0.55 * sin(pi / 4) / 72 + 0.2
  
  # legend_h_in: inches reserved at bottom for horizontal legend strip
  legend_h_in   <- 1.8
  
  fig_h         <- max(5.0, n_rows * cellheight / 72 + col_label_h + 2.0 + legend_h_in)
  
  annot_row_w   <- ncol(annot_row) * cellwidth / 72 + 0.3
  row_label_w   <- max_row_chars * fontsize_row * 0.60 / 72 + 0.8
  fig_w         <- max(8.0, n_cols * cellwidth / 72 + annot_row_w + row_label_w + 1.2)
  
  return(list(p             = p,
              fig_w         = fig_w,
              fig_h         = fig_h,
              legend_h_in   = legend_h_in,
              eff_range     = eff_range,
              effect_col    = effect_col,
              cat_colors    = cat_colors,
              cats_present  = cats_present_final,
              bar_left      = bar_left,
              bar_y         = bar_y))
}

# ── draw_pwy_met_heatmap() ────────────────────────────────────────────────────
  draw_pwy_met_heatmap <- function(hm, legend_y = NULL,
                                   effect_label = "Spearman Correlation") {
    
    if (is.null(hm)) return(invisible(NULL))
    
    legend_h_in <- hm$legend_h_in
    
    grid.newpage()
    
    # Heatmap: full width, everything above the bottom legend strip
    pushViewport(viewport(
      x      = unit(0, "npc"),
      y      = unit(legend_h_in, "inches"),
      width  = unit(1, "npc"),
      height = unit(1, "npc") - unit(legend_h_in, "inches"),
      just   = c("left", "bottom")
    ))
    grid.draw(hm$p$gtable)
    popViewport()
    
    # ── Bottom legend strip ───────────────────────────────────────────────────
    pushViewport(viewport(
      x      = unit(0, "npc"),
      y      = unit(0, "npc"),
      width  = unit(1, "npc"),
      height = unit(legend_h_in, "inches"),
      just   = c("left", "bottom")
    ))
    
    pal <- colorRampPalette(brewer.pal(11, "RdBu"))(100)
    n   <- length(pal)
    
    # ── Horizontal colour bar ─────────────────────────────────────────────────
    bar_left  <- hm$bar_left
    bar_right <- bar_left + 2.8
    bar_y     <- hm$bar_y
    bar_h     <- 0.20
    each_w    <- (bar_right - bar_left) / n
    
    # Title above bar
    grid.text(effect_label,
              x    = unit(bar_left, "inches"),
              y    = unit(bar_y + bar_h / 2 + 0.22, "inches"),
              gp   = gpar(fontsize = 10, fontface = "bold"),
              just = c("left", "bottom"))
    
    for (i in seq_len(n)) {
      grid.rect(x      = unit(bar_left + (i - 1) * each_w, "inches"),
                y      = unit(bar_y, "inches"),
                width  = unit(each_w, "inches"),
                height = unit(bar_h, "inches"),
                just   = c("left", "center"),
                gp     = gpar(fill = pal[i], col = NA))
    }
    
    # Tick labels below bar
    eff_range  <- hm$eff_range
    tick_vals  <- c(-eff_range, -eff_range / 2, 0, eff_range / 2, eff_range)
    tick_fracs <- (tick_vals + eff_range) / (2 * eff_range)
    
    for (i in seq_along(tick_vals)) {
      grid.text(sprintf("%.2f", tick_vals[i]),
                x    = unit(bar_left + tick_fracs[i] * (bar_right - bar_left), "inches"),
                y    = unit(bar_y - bar_h / 2 - 0.08, "inches"),
                just = c("center", "top"),
                gp   = gpar(fontsize = 8))
    }
    
    # ── PFAS Association legend (right of colour bar) ─────────────────────────
    dir_left <- bar_right + 0.6
    
    grid.text("PFAS Association:",
              x    = unit(dir_left, "inches"),
              y    = unit(bar_y + bar_h / 2 + 0.22, "inches"),
              gp   = gpar(fontsize = 9, fontface = "bold"),
              just = c("left", "bottom"))
    
    bw <- unit(0.22, "inches")
    bh <- unit(0.14, "inches")
    
    grid.rect(x      = unit(dir_left, "inches"),
              y      = unit(bar_y + 0.02, "inches"),
              width  = bw, height = bh,
              just   = c("left", "center"),
              gp     = gpar(fill = "purple3", col = NA))
    grid.text("Positive",
              x    = unit(dir_left + 0.30, "inches"),
              y    = unit(bar_y + 0.02, "inches"),
              gp   = gpar(fontsize = 9),
              just = c("left", "center"))
    
    grid.rect(x      = unit(dir_left + 1.2, "inches"),
              y      = unit(bar_y + 0.02, "inches"),
              width  = bw, height = bh,
              just   = c("left", "center"),
              gp     = gpar(fill = "orange3", col = NA))
    grid.text("Negative",
              x    = unit(dir_left + 1.50, "inches"),
              y    = unit(bar_y + 0.02, "inches"),
              gp   = gpar(fontsize = 9),
              just = c("left", "center"))
    
    # ── Pathway Category legend (right of PFAS, two horizontal rows) ──────────
    cats       <- hm$cats_present
    cat_colors <- hm$cat_colors
    pwy_left   <- dir_left + 2.5   # ← adjust to shift whole block left/right
    
    grid.text("Pathway Category:",
              x    = unit(pwy_left, "inches"),
              y    = unit(bar_y + bar_h / 2 + 0.22, "inches"),
              gp   = gpar(fontsize = 9, fontface = "bold"),
              just = c("left", "bottom"))
    
    # Split categories into two rows
    n_top <- ceiling(length(cats) / 2)
    row1  <- cats[seq_len(n_top)]
    row2  <- cats[seq(n_top + 1, length(cats))]
    
    # Row 1 (upper, aligned with PFAS boxes)
    cat_x <- pwy_left
    for (cat_name in row1) {
      grid.rect(x      = unit(cat_x, "inches"),
                y      = unit(bar_y + 0.02, "inches"),
                width  = bw, height = bh,
                just   = c("left", "center"),
                gp     = gpar(fill = cat_colors[cat_name], col = NA))
      grid.text(cat_name,
                x    = unit(cat_x + 0.28, "inches"),
                y    = unit(bar_y + 0.02, "inches"),
                gp   = gpar(fontsize = 7),
                just = c("left", "center"))
      cat_x <- cat_x + 0.28 + nchar(cat_name) * 0.055 + 0.12
    }
    
    # Row 2 (lower)
    cat_x <- pwy_left
    for (cat_name in row2) {
      grid.rect(x      = unit(cat_x, "inches"),
                y      = unit(bar_y - 0.20, "inches"),
                width  = bw, height = bh,
                just   = c("left", "center"),
                gp     = gpar(fill = cat_colors[cat_name], col = NA))
      grid.text(cat_name,
                x    = unit(cat_x + 0.28, "inches"),
                y    = unit(bar_y - 0.20, "inches"),
                gp   = gpar(fontsize = 7),
                just = c("left", "center"))
      cat_x <- cat_x + 0.28 + nchar(cat_name) * 0.055 + 0.12
    }
    
    popViewport()
  }


# ── Part C.1: Individual heatmap - Spearman ───────────────────────────────────
cat("--- Part C.1: Individual Spearman heatmap ---\n")

hm_spearman <- build_pwy_met_heatmap(
  result_df    = spearman_df,
  effect_col   = "rho",
  mwas_c18     = mwas_c18_6m,
  mwas_hilic   = mwas_hilic_6m,
  pwy_results  = pwy_results_6m,
  fdr_thresh   = 0.05,
  min_assoc    = 5,
  cellwidth    = 10,
  cellheight   = 10,
  fontsize_col = 7,
  fontsize_row = 9,
  bar_left     = 1.0,
  bar_y        = 1.9
)

if (!is.null(hm_spearman)) {
  pdf(here::here("out_figures",
                 "heatmap_pathway_metabolite_individual_spearman_6m.pdf"),
      width  = hm_spearman$fig_w,
      height = hm_spearman$fig_h)
  draw_pwy_met_heatmap(hm_spearman,
                       legend_y     = 0.90,
                       effect_label = "Spearman Correlation")
  dev.off()
  cat("Saved: heatmap_pathway_metabolite_individual_spearman_6m.pdf\n")
}


# ── Part C.2: Individual heatmap - Spearman 1m ────────────────────────────────
cat("\n--- Part C.2: Individual Spearman heatmap (1m) ---\n")

if (length(pwy_cols_use_1m) > 0 && length(met_cols_use_1m) > 0) {
  hm_spearman_1m <- build_pwy_met_heatmap(
    result_df    = spearman_df_1m,
    effect_col   = "rho",
    mwas_c18     = mwas_c18_1m,
    mwas_hilic   = mwas_hilic_1m,
    pwy_results  = pwy_results_1m,
    fdr_thresh   = 0.05,
    min_assoc    = 2,
    cellwidth    = 10,
    cellheight   = 10,
    fontsize_col = 7,
    fontsize_row = 9,
    bar_left     = 1.0,
    bar_y        = 1.9
  )

  if (!is.null(hm_spearman_1m)) {
    pdf(here::here("out_figures",
                   "heatmap_pathway_metabolite_individual_spearman_1m.pdf"),
        width  = hm_spearman_1m$fig_w,
        height = hm_spearman_1m$fig_h)
    draw_pwy_met_heatmap(hm_spearman_1m,
                         legend_y     = 0.90,
                         effect_label = "Spearman Correlation")
    dev.off()
    cat("Saved: heatmap_pathway_metabolite_individual_spearman_1m.pdf\n")
  } else {
    cat("1m heatmap: no pairs survived filter.\n")
  }
} else {
  cat("Insufficient 1m pathway/metabolite data for individual heatmap.\n")
}

