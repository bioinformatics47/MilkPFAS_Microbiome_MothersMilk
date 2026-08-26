# header -----------------------------------------------------------------------
#
# TITLE:   4. Alpha Diversity Analyses.R
#
# PURPOSE: Single-exposure and mixture (qgcomp) associations between
#          1-month PFAS and alpha diversity at 1 month and 6 months
#
# DATE:    February 2026
#
#   Reviewed by Ellie Holzhausen (EAH) on April 23, 2026
# by Haonan Li (HL) on May 25, 2026
#
# set up -----------------------------------------------------------------------
rm(list = ls())
options(scipen = 0)

library(tidyr); library(reshape); library(dplyr)
library(purrr); library(stringr); library(lme4); library(writexl)
library(pscl); library(tibble); library(readxl)
library(here); library(ggplot2); library(vegan)
library(tidyverse); library(qgcomp)
library(patchwork)
library(flextable)
library(officer)


# Read alpha diversity ---------------------------------------------------------
alpha_div <- read.table(
  here::here("input", "alphaDiv_repeated_rarefied_100_data_readDepth_1000000_mothersMilk_replacementFALSE.tsv"),
  sep = '\t', header = TRUE)

# Load IDs manifest to map MG sample IDs to merge_id_dyad
IDs <- read_excel(here::here("input", "GORAN_MICROBIOME_MANIFEST 10-2022.xlsx"))
colnames(IDs) <- IDs[1,]
IDs <- IDs[-1,]

# Fix duplicate column names
colnames(IDs) <- make.unique(colnames(IDs), sep = "_")

IDs$dyad_id2      <- str_pad(substr(IDs$`Old Together`, 9, 11), width = 4, pad = "0")
IDs$timepoint_str <- substr(IDs$`Old Together`, 1, 2)
IDs$merge_id_dyad <- paste0("MM-", IDs$dyad_id2, "-", IDs$timepoint_str)
IDs$core_id       <- paste0("MG", IDs$`CORE ID`)

# Keep only 1m and 6m
IDs_1m6m <- IDs %>%
  filter(timepoint_str %in% c("01", "06")) %>%
  dplyr::select(core_id, merge_id_dyad, timepoint_str) %>%
  mutate(timepoint = as.integer(timepoint_str))

# Merge alpha diversity with IDs
# MG prefix already present in SampleID — no paste0 needed
# Lab controls (L, S, W samples) dropped via inner_join (no match in IDs_1m6m)
alpha_div <- alpha_div %>%
  dplyr::rename(core_id = SampleID) %>%
  inner_join(IDs_1m6m, by = "core_id") %>%   # inner_join drops unmatched (controls)
  dplyr::rename(
    linkerID                 = merge_id_dyad,
    ShannonDiv_repeated_rare = Shannon,
    Richness_repeated_rare   = Richness,
    Simpson_repeated_rare    = Simpson,
    Evenness_repeated_rare   = Evenness
  )

cat("Alpha div total rows:", nrow(alpha_div), "\n")
cat("Timepoint distribution:\n")
print(table(alpha_div$timepoint_str))
cat("Sample linkerID format:", head(alpha_div$linkerID), "\n")

alpha_div <- alpha_div %>%
  dplyr::rename(merge_id_dyad = linkerID) %>%
  dplyr::select(merge_id_dyad, timepoint, ShannonDiv_repeated_rare,
                Richness_repeated_rare, Simpson_repeated_rare, Evenness_repeated_rare)

hist(alpha_div$ShannonDiv_repeated_rare)
hist(alpha_div$Richness_repeated_rare)
hist(alpha_div$Simpson_repeated_rare)
hist(alpha_div$Evenness_repeated_rare)

# Read the 4 data files --------------------------------------------------------
PFAS1m_micro1m       <- read.csv(here::here("out_files", "PFAS1m_micro1m_species.csv")) %>% dplyr::select(-X)
PFAS1m_micro6m       <- read.csv(here::here("out_files", "PFAS1m_micro6m_species.csv")) %>% dplyr::select(-X)
PFAS1mDetect_micro1m <- read.csv(here::here("out_files", "PFAS1mDetect_micro1m_species.csv")) %>% dplyr::select(-X)
PFAS1mDetect_micro6m <- read.csv(here::here("out_files", "PFAS1mDetect_micro6m_species.csv")) %>% dplyr::select(-X)

# Variable lists ---------------------------------------------------------------
# Continuously measured PFAS
pfas_vars_continuous <- c("PFBS", "PFHxS", "PFNA", "PFOA", "PFOS")

# Binary classified PFAS
pfas_vars_binary <- c("N.MeFOSAA", "PFBA", "PFDA", "PFDoA", "PFHpA")

# Rename dictionary for continuous files:
pfas_rename_continuous <- c(
  "PFBS"  = "PFBS_pgmL",
  "PFHxS" = "PFHxS_pgmL",
  "PFNA"  = "PFNA_pgmL",
  "PFOA"  = "PFOA_pgmL",
  "PFOS"  = "PFOS_pgmL"
)

# Covariates used in all models
# mode_of_delivery coded as binary (0 = Vaginal, 1 = C-Section)
# gestational_age coded as two dummy variables with Ontime as reference (0):
#   gest_Early = 1 if Early, 0 otherwise
#   gest_Late  = 1 if Late,  0 otherwise
covariates <- c(
  "SES_index_final", "gest_Early", "gest_Late", "breastmilk_per_day",
  "mode_of_delivery_bin", "baby_birthweight_kg"
)

# Alpha diversity metrics
alpha_metrics <- c(
  "ShannonDiv_repeated_rare",
  "Simpson_repeated_rare",
  "Richness_repeated_rare",
  "Evenness_repeated_rare"
)

# PFAS ordering for output tables
# Continuous: 5 individual PFAS + Mixture at the end
# Binary:     5 individual PFAS + n_detect at the end
pfas_order_continuous <- c("PFBS", "PFHxS", "PFNA", "PFOA", "PFOS", "Mixture")
pfas_order_binary     <- c("N.MeFOSAA", "PFBA", "PFDA", "PFDoA", "PFHpA", "n_detect")

# Compute IQR of each log2-transformed continuous PFAS-------------------------
# Used for post-hoc IQR scaling of individual betas and for qgcomp IQR inputs
compute_pfas_iqr <- function(data) {
  sapply(pfas_vars_continuous, function(pfas) {
    IQR(data[[pfas]], na.rm = TRUE)
  })
}

# Prepare_data------------------------------------------------------------------
prepare_data <- function(pfas_data, alpha_timepoint, scenario_name,
                         pfas_type = "continuous") {
  
  # Filter alpha diversity for the specified timepoint
  alpha_div_filtered <- alpha_div %>%
    dplyr::filter(timepoint == alpha_timepoint) %>%
    dplyr::select(merge_id_dyad, all_of(alpha_metrics))
  
  # Join with PFAS data
  comb <- left_join(pfas_data, alpha_div_filtered, by = "merge_id_dyad")
  
  if (pfas_type == "continuous") {
    comb <- comb %>%
      dplyr::rename(!!!pfas_rename_continuous)
    pfas_vars <- pfas_vars_continuous
    
  } else {
    detect_rename <- setNames(
      paste0(pfas_vars_binary, "_Detect"),
      pfas_vars_binary
    )
    comb <- comb %>%
      dplyr::rename(!!!detect_rename)
    pfas_vars <- pfas_vars_binary
    comb <- comb %>%
      mutate(across(all_of(pfas_vars),
                    ~ factor(., levels = c("non-detect", "detect"))))
  }
  
  comb <- comb %>%
    mutate(
      gest_Early           = ifelse(gestational_age_cat  == "Early",     1, 0),
      gest_Late            = ifelse(gestational_age_cat  == "Late",      1, 0),
      mode_of_delivery_bin = ifelse(mode_of_delivery_cat == "C-Section", 1, 0)
    )
  
  # Remove rows with NAs in PFAS variables
  n_before <- nrow(comb)
  comb <- comb[complete.cases(comb[, pfas_vars]), ]
  cat("Rows removed due to NA in PFAS:", n_before - nrow(comb), "\n")
  
  # Remove rows with NAs in alpha diversity metrics
  n_before <- nrow(comb)
  comb <- comb[complete.cases(comb[, alpha_metrics]), ]
  cat("Rows removed due to NA in alpha diversity:", n_before - nrow(comb), "\n")
  cat("Final N:", nrow(comb), "\n")
  
  attr(comb, "pfas_vars") <- pfas_vars
  return(comb)
}

# Run_analysis for single-exposure (lm models)----------------------------------
run_analysis <- function(data, scenario_name) {
  pfas_vars    <- attr(data, "pfas_vars")
  results_list <- list()
  
  for (pfas in pfas_vars) {
    cat("  Running models for", pfas, "...\n")
    
    results_list[[pfas]] <- list(
      shannon  = lm(as.formula(paste("ShannonDiv_repeated_rare ~", pfas, "+", paste(covariates, collapse = "+"))), data = data),
      simpson  = lm(as.formula(paste("Simpson_repeated_rare ~",    pfas, "+", paste(covariates, collapse = "+"))), data = data),
      richness = lm(as.formula(paste("Richness_repeated_rare ~",   pfas, "+", paste(covariates, collapse = "+"))), data = data),
      evenness = lm(as.formula(paste("Evenness_repeated_rare ~",   pfas, "+", paste(covariates, collapse = "+"))), data = data)
    )
  }
  return(results_list)
}


# Run_ndetect  (n_detect mixture score models)----------------------------------
# n_detect is a continuous count (0-5) of the number of PFAS detected per
# participant. It is not in pfas_vars_binary and is passed through prepare_data
# unchanged. Models use the same covariates as run_analysis
run_ndetect <- function(data, scenario_name) {
  results_list <- list()
  results_list[["n_detect"]] <- list(
    shannon  = lm(ShannonDiv_repeated_rare ~ n_detect + SES_index_final + gest_Early + gest_Late +
                    breastmilk_per_day + mode_of_delivery_bin + baby_birthweight_kg, data = data),
    simpson  = lm(Simpson_repeated_rare   ~ n_detect + SES_index_final + gest_Early + gest_Late +
                    breastmilk_per_day + mode_of_delivery_bin + baby_birthweight_kg, data = data),
    richness = lm(Richness_repeated_rare  ~ n_detect + SES_index_final + gest_Early + gest_Late +
                    breastmilk_per_day + mode_of_delivery_bin + baby_birthweight_kg, data = data),
    evenness = lm(Evenness_repeated_rare  ~ n_detect + SES_index_final + gest_Early + gest_Late +
                    breastmilk_per_day + mode_of_delivery_bin + baby_birthweight_kg, data = data)
  )
  return(results_list)
}

# Create_output_table  (single-exposure results)--------------------------------
create_output_table <- function(results_list, pfas_order, pfas_iqr = NULL) {
  
  output_table <- data.frame()
  
  for (pfas in names(results_list)) {
    models <- results_list[[pfas]]
    
    for (metric in c("Richness", "Evenness", "Shannon", "Simpson")) {
      model_obj  <- models[[tolower(metric)]]
      coef_table <- summary(model_obj)$coefficients
      
      coef_row <- grep(paste0("^", pfas), rownames(coef_table), value = TRUE)
      output_table[paste0(metric, "_", pfas), "Beta"] <- coef_table[coef_row, "Estimate"]
      output_table[paste0(metric, "_", pfas), "SE"]   <- coef_table[coef_row, "Std. Error"]
      output_table[paste0(metric, "_", pfas), "P"]    <- coef_table[coef_row, "Pr(>|t|)"]
      # IQR-scaled beta — only for continuous PFAS (not binary, not n_detect)
      if (!is.null(pfas_iqr) && pfas %in% names(pfas_iqr)) {
        output_table[paste0(metric, "_", pfas), "Beta_IQR"] <-
          coef_table[coef_row, "Estimate"] * pfas_iqr[pfas]
        output_table[paste0(metric, "_", pfas), "SE_IQR"] <-
          coef_table[coef_row, "Std. Error"] * pfas_iqr[pfas]
      }
    }
  }
  
  output_table <- rownames_to_column(output_table, var = "Measure") %>%
    mutate(Measure = trimws(Measure)) %>%
    separate(Measure, into = c("Measure", "PFAS"),
             sep = "_", remove = TRUE, extra = "merge", fill = "right")
  
  metric_order <- c("Richness", "Evenness", "Shannon", "Simpson")
  
  output_ordered <- output_table %>%
    mutate(
      Measure      = str_trim(Measure),
      PFAS         = str_trim(PFAS),
      measure_rank = match(Measure, metric_order),
      pfas_rank    = match(PFAS, pfas_order)   # proper ordering by PFAS name
    ) %>%
    arrange(measure_rank, pfas_rank) %>%
    select(Measure, PFAS, Beta, SE, P, any_of(c("Beta_IQR", "SE_IQR"))) %>%
    mutate(Measure = factor(Measure, levels = metric_order, ordered = TRUE))
  
  return(output_ordered)
}

# Run_qgcomp_combined: runs both raw (q=4) and IQR-scaled (q=NULL) qgcomp-------
# Returns Beta, SE, P (raw q=4 scale) and Beta_IQR, SE_IQR, P_IQR (IQR scale)
# in a single data frame per metric — no need to run twice separately
# PFAS divided by their IQR before modeling so psi represents change per
# simultaneous 1-IQR increase across all PFAS — comparable to IQR-scaled
# individual betas
# NOTE: qgcomp is NOT run for binary PFAS scenarios. Binary exposures violate
# qgcomp's assumption of continuous, orderable exposures.

run_qgcomp_combined <- function(data, pfas_iqr, scenario_name) {
  
  pfas_vars <- attr(data, "pfas_vars")
  
  # IQR-scaled copy of data
  data_iqr <- data
  for (pfas in pfas_vars_continuous) {
    data_iqr[[pfas]] <- data_iqr[[pfas]] / pfas_iqr[pfas]
  }
  
  qgcomp_results <- data.frame()
  
  for (metric in c("Richness", "Evenness", "Shannon", "Simpson")) {
    
    outcome_col <- c(
      Richness = "Richness_repeated_rare",
      Evenness = "Evenness_repeated_rare",
      Shannon  = "ShannonDiv_repeated_rare",
      Simpson  = "Simpson_repeated_rare"
    )[[metric]]
    
    formula_qg <- as.formula(
      paste0(outcome_col, " ~ ",
             paste(c(pfas_vars_continuous, covariates), collapse = " + "))
    )
    
    # Raw model (q=4)
    qgmod_raw <- qgcomp.boot(
      formula_qg, expnms = pfas_vars_continuous,
      data = data, family = gaussian(),
      q = 4, B = 5000, seed = 2024
    )
    psi_se_raw <- summary(qgmod_raw)$coefficients["psi1", "Std. Error"]
    
    # IQR-scaled model (q=NULL)
    qgmod_iqr <- qgcomp.boot(
      formula_qg, expnms = pfas_vars_continuous,
      data = data_iqr, family = gaussian(),
      q = NULL, B = 5000, seed = 2024
    )
    psi_se_iqr <- summary(qgmod_iqr)$coefficients["psi1", "Std. Error"]
    
    qgcomp_results <- rbind(qgcomp_results, data.frame(
      Measure  = metric,
      PFAS     = "Mixture",
      Beta     = qgmod_raw$psi,
      SE       = psi_se_raw,
      P        = qgmod_raw$pval[2],
      Beta_IQR = qgmod_iqr$psi,
      SE_IQR   = psi_se_iqr,
      stringsAsFactors = FALSE
    ))
  }
  return(qgcomp_results)
}

# MAIN ANALYSIS-----------------------------------------------------------------
# Scenario 1: Continuous 1m PFAS + 1m alpha diversity
data_cont_1m_1m <- prepare_data(PFAS1m_micro1m, alpha_timepoint = 1,
                                scenario_name = "Continuous 1m PFAS + 1m Microbiome",
                                pfas_type = "continuous")
pfas_iqr_1m     <- compute_pfas_iqr(data_cont_1m_1m)
cat("PFAS IQRs (1m data):\n"); print(round(pfas_iqr_1m, 4))

results_cont_1m_1m  <- run_analysis(data_cont_1m_1m, "Continuous 1m PFAS + 1m Microbiome")
output_cont_1m_1m   <- create_output_table(results_cont_1m_1m, pfas_order_continuous,
                                           pfas_iqr = pfas_iqr_1m)
qgcomp_cont_1m_1m <- run_qgcomp_combined(data_cont_1m_1m, pfas_iqr_1m,
                                         "Continuous 1m PFAS + 1m Microbiome")

# Combine: individual betas + mixture
combined_cont_1m_1m <- bind_rows(output_cont_1m_1m, qgcomp_cont_1m_1m)

# Save file
write.csv(combined_cont_1m_1m,
          "out_files/alphaDiv_continuous_1m_1m.csv", row.names = FALSE)


# Scenario 2: Continuous 1m PFAS + 6m alpha diversity
data_cont_1m_6m <- prepare_data(PFAS1m_micro6m, alpha_timepoint = 6,
                                scenario_name = "Continuous 1m PFAS + 6m Microbiome",
                                pfas_type = "continuous")
pfas_iqr_6m     <- compute_pfas_iqr(data_cont_1m_6m)
cat("PFAS IQRs (6m data):\n"); print(round(pfas_iqr_6m, 4))

results_cont_1m_6m  <- run_analysis(data_cont_1m_6m, "Continuous 1m PFAS + 6m Microbiome")
output_cont_1m_6m   <- create_output_table(results_cont_1m_6m, pfas_order_continuous,
                                           pfas_iqr = pfas_iqr_6m)
qgcomp_cont_1m_6m <- run_qgcomp_combined(data_cont_1m_6m, pfas_iqr_6m,
                                         "Continuous 1m PFAS + 6m Microbiome")

#Combine: individual betas + mixture
combined_cont_1m_6m <- bind_rows(output_cont_1m_6m, qgcomp_cont_1m_6m)

# Save file
write.csv(combined_cont_1m_6m,
          "out_files/alphaDiv_continuous_1m_6m.csv", row.names = FALSE)

# Scenario 3: Binary 1m PFAS + 1m alpha diversity
data_bin_1m_1m    <- prepare_data(PFAS1mDetect_micro1m, alpha_timepoint = 1,
                                  scenario_name = "Binary 1m PFAS + 1m Microbiome",
                                  pfas_type = "binary")
results_bin_1m_1m <- run_analysis(data_bin_1m_1m, "Binary 1m PFAS + 1m Microbiome")
output_bin_1m_1m  <- create_output_table(results_bin_1m_1m, pfas_order_binary)
results_ndetect_1m_1m <- run_ndetect(data_bin_1m_1m, "n_detect + 1m Microbiome")
output_ndetect_1m_1m  <- create_output_table(results_ndetect_1m_1m, pfas_order_binary)
combined_bin_1m_1m <- rbind(output_bin_1m_1m, output_ndetect_1m_1m)
write.csv(combined_bin_1m_1m, "out_files/alphaDiv_binary_1m_1m.csv", row.names = FALSE)

# Scenario 4: Binary 1m PFAS + 6m alpha diversity
data_bin_1m_6m    <- prepare_data(PFAS1mDetect_micro6m, alpha_timepoint = 6,
                                  scenario_name = "Binary 1m PFAS + 6m Microbiome",
                                  pfas_type = "binary")
results_bin_1m_6m <- run_analysis(data_bin_1m_6m, "Binary 1m PFAS + 6m Microbiome")
output_bin_1m_6m  <- create_output_table(results_bin_1m_6m, pfas_order_binary)
results_ndetect_1m_6m <- run_ndetect(data_bin_1m_6m, "n_detect + 6m Microbiome")
output_ndetect_1m_6m  <- create_output_table(results_ndetect_1m_6m, pfas_order_binary)
combined_bin_1m_6m <- rbind(output_bin_1m_6m, output_ndetect_1m_6m)
write.csv(combined_bin_1m_6m, "out_files/alphaDiv_binary_1m_6m.csv", row.names = FALSE)


# Visualize PFAS-Alpha Diversity Results----------------------------------------

# ── Load saved results ─────────────────────────────────────────────────────────
combined_cont_1m_1m <- read.csv("out_files/alphaDiv_continuous_1m_1m.csv")
combined_cont_1m_6m <- read.csv("out_files/alphaDiv_continuous_1m_6m.csv")
combined_bin_1m_1m  <- read.csv("out_files/alphaDiv_binary_1m_1m.csv")
combined_bin_1m_6m  <- read.csv("out_files/alphaDiv_binary_1m_6m.csv")

# ── Add scenario labels ────────────────────────────────────────────────────────
combined_cont_1m_1m$Scenario <- "1m PFAS → 1m Microbiome"
combined_cont_1m_6m$Scenario <- "1m PFAS → 6m Microbiome"
combined_bin_1m_1m$Scenario  <- "1m PFAS → 1m Microbiome"
combined_bin_1m_6m$Scenario  <- "1m PFAS → 6m Microbiome"

# Combine continuous and binary separately
all_continuous <- bind_rows(combined_cont_1m_1m, combined_cont_1m_6m)

all_binary <- bind_rows(combined_bin_1m_1m, combined_bin_1m_6m) %>%
  mutate(PFAS = recode(PFAS,
                       "N.MeFOSAA" = "N-MeFOSAA",
                       "n_detect"  = "N-detect"))

# Format cell labels and set factor levels
# Cell shows: Beta value with * if p < 0.05
# Richness has large Betas (integer-like), others are small decimals
# so we use different decimal formatting per metric

format_beta <- function(beta, measure) {
  ifelse(measure == "Richness",
         formatC(beta, format = "f", digits = 1),   # e.g. 210.3
         formatC(beta, format = "f", digits = 2))   # e.g. 0.052
}

all_continuous <- all_continuous %>%
  mutate(
    cell_label = paste0(format_beta(dplyr::coalesce(Beta_IQR, Beta), Measure),
                        ifelse(P < 0.05, "*", "")),
    Significant = ifelse(P < 0.05, "Yes", "No"),
    PFAS     = factor(PFAS,    levels = c("PFBS", "PFHxS", "PFNA", "PFOA", "PFOS", "Mixture")),
    Measure  = factor(Measure, levels = c("Richness", "Evenness", "Shannon", "Simpson")),
    Scenario = factor(Scenario, levels = c("1m PFAS → 1m Microbiome",
                                           "1m PFAS → 6m Microbiome"))
  )

all_binary <- all_binary %>%
  mutate(
    cell_label = paste0(format_beta(Beta, Measure),
                        ifelse(P < 0.05, "*", "")),
    Significant = ifelse(P < 0.05, "Yes", "No"),
    PFAS     = factor(PFAS,    levels = c("N-MeFOSAA", "PFBA", "PFDA", "PFDoA", "PFHpA", "N-detect")),
    Measure  = factor(Measure, levels = c("Richness", "Evenness", "Shannon", "Simpson")),
    Scenario = factor(Scenario, levels = c("1m PFAS → 1m Microbiome",
                                           "1m PFAS → 6m Microbiome"))
  )

# Plot function
# make_tablefig <- function(data) {
#   ggplot(data, aes(x = PFAS, y = Measure)) +
#     geom_tile(aes(fill = Significant), color = "white", linewidth = 0.5) +
#     scale_fill_manual(
#       values = c("Yes" = "#F8D7D7", "No" = "white"),
#       guide  = "none"
#     ) +
make_tablefig <- function(data, sig_color = "#F8D7D7") {
  ggplot(data, aes(x = PFAS, y = Measure)) +
    geom_tile(aes(fill = Significant), color = "white", linewidth = 0.5) +
    scale_fill_manual(
      values = c("Yes" = sig_color, "No" = "white"),
      guide  = "none"
    ) +
    geom_text(aes(label = cell_label), size = 3.2, color = "black") +
    facet_wrap(~Scenario, nrow = 1) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 10,
                                      face = "bold", color = "black"),
      axis.text.y      = element_text(size = 12, face = "bold", color = "black"),
      axis.title       = element_blank(),
      panel.grid       = element_blank(),
      panel.spacing    = unit(1.5, "lines"),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
      strip.text       = element_text(face = "bold", size = 10, color = "white"),
      strip.background = element_rect(color = "black", fill = "black", linewidth = 1),
      legend.position  = "none",
      plot.margin      = margin(5, 15, 5, 5),
      plot.caption     = element_text(size = 7, hjust = 0, color = "grey40")
    ) +
    labs(caption = NULL)
}

# Create plots
p_continuous <- make_tablefig(all_continuous, sig_color = "#F8D7D7")
print(p_continuous)

p_binary     <- make_tablefig(all_binary, sig_color = "#D7EBF8")
print(p_binary)

# Save plots
ggsave(
  filename = here::here("out_figures", "PFAS_alpha_tablefig_continuous.pdf"),
  plot     = p_continuous,
  width    = 9,
  height   = 2.5,
  dpi      = 600,
  bg       = "white"
)

ggsave(
  filename = here::here("out_figures", "PFAS_alpha_tablefig_binary.pdf"),
  plot     = p_binary,
  width    = 9,
  height   = 2.5,
  dpi      = 600,
  bg       = "white"
)

# END



# Alpha Diversity Table — flextable output
# Layout: rows = PFAS, columns = 4 alpha metrics × 2 scenarios (1m, 6m)
#         grouped by timepoint block, then continuous / binary sub-blocks
#
# Estimate shown: IQR-scaled Beta (Beta_IQR) for continuous,
#                 raw Beta for binary/n_detect
# Significant p-values (p < 0.05) → italic bold

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(here)


# ── 1. Load saved CSVs ────────────────────────────────────────────────────────
combined_cont_1m_1m <- read_csv("out_files/alphaDiv_continuous_1m_1m.csv")
combined_cont_1m_6m <- read_csv("out_files/alphaDiv_continuous_1m_6m.csv")
combined_bin_1m_1m  <- read_csv("out_files/alphaDiv_binary_1m_1m.csv")
combined_bin_1m_6m  <- read_csv("out_files/alphaDiv_binary_1m_6m.csv")

# ── 2. Choose estimate column ─────────────────────────────────────────────────
# Continuous: use Beta_IQR; binary/n_detect: use Beta
prep_continuous <- function(df, timepoint_label) {
  df %>%
    mutate(
      est  = coalesce(Beta_IQR, Beta),  # Beta_IQR for PFAS, Beta for Mixture fallback
      pfas_type = "Continuous",
      timepoint = timepoint_label,
      PFAS = recode(PFAS, "Mixture" = "Mixture (qgcomp)")
    ) %>%
    dplyr::select(PFAS, Measure, est, P, pfas_type, timepoint)
}

prep_binary <- function(df, timepoint_label) {
  df %>%
    mutate(
      est  = Beta,
      pfas_type = "Binary",
      timepoint = timepoint_label,
      PFAS = recode(PFAS,
                    "N.MeFOSAA" = "N-MeFOSAA",
                    "n_detect"  = "N-detect")
    ) %>%
    dplyr::select(PFAS, Measure, est, P, pfas_type, timepoint)
}

all_data <- bind_rows(
  prep_continuous(combined_cont_1m_1m, "1m PFAS → 1m Microbiome"),
  prep_continuous(combined_cont_1m_6m, "1m PFAS → 6m Microbiome"),
  prep_binary(combined_bin_1m_1m,      "1m PFAS → 1m Microbiome"),
  prep_binary(combined_bin_1m_6m,      "1m PFAS → 6m Microbiome")
)

# ── 3. Format cell content ────────────────────────────────────────────────────
# Richness has large betas → 1 decimal; others → 3 decimals
# Cell text: "β  (p = x.xxx)"
all_data <- all_data %>%
  mutate(
    est_fmt = case_when(
      Measure == "Richness" ~ formatC(est, format = "f", digits = 1),
      TRUE                  ~ formatC(est, format = "f", digits = 3)
    ),
    p_fmt = formatC(P, format = "f", digits = 3),
    cell_text = paste0(est_fmt, "\n(p = ", p_fmt, ")"),
    sig = P < 0.05
  )

# ── 4. Pivot to wide ──────────────────────────────────────────────────────────
# Columns: metric × timepoint × pfas_type combinations
# We'll build one combined table with a hierarchical column structure

metrics <- c("Richness", "Evenness", "Shannon", "Simpson")
timepoints <- c("1m PFAS → 1m Microbiome", "1m PFAS → 6m Microbiome")
pfas_types <- c("Continuous", "Binary")

# PFAS row ordering
pfas_order_cont <- c("PFBS", "PFHxS", "PFNA", "PFOA", "PFOS", "Mixture (qgcomp)")
pfas_order_bin  <- c("N-MeFOSAA", "PFBA", "PFDA", "PFDoA", "PFHpA", "N-detect")

# Build separate wide tables for continuous and binary, then row-bind
make_wide <- function(data_sub, pfas_order) {
  base <- data_sub %>%
    mutate(col_id = paste(timepoint, Measure, sep = "|||"))
  
  wide_text <- base %>%
    dplyr::select(PFAS, col_id, cell_text) %>%
    pivot_wider(id_cols = PFAS, names_from = col_id,
                values_from = cell_text) %>%
    rename_with(~ paste0("cell_text~~", .), -PFAS)
  
  wide_sig <- base %>%
    dplyr::select(PFAS, col_id, sig) %>%
    pivot_wider(id_cols = PFAS, names_from = col_id,
                values_from = sig) %>%
    rename_with(~ paste0("sig~~", .), -PFAS)
  
  left_join(wide_text, wide_sig, by = "PFAS") %>%
    mutate(PFAS = factor(PFAS, levels = pfas_order)) %>%
    arrange(PFAS) %>%
    mutate(PFAS = as.character(PFAS))
}

wide_cont <- make_wide(
  filter(all_data, pfas_type == "Continuous"),
  pfas_order_cont
) %>% mutate(pfas_type = "Continuous")

wide_bin <- make_wide(
  filter(all_data, pfas_type == "Binary"),
  pfas_order_bin
) %>% mutate(pfas_type = "Binary")

# ── 5. Define column order ─────────────────────────────────────────────────────
# Want: PFAS | [1m→1m: R, E, S, Si] | [1m→6m: R, E, S, Si]  (for cell_text cols)
tp1 <- "1m PFAS → 1m Microbiome"
tp2 <- "1m PFAS → 6m Microbiome"

cell_cols <- c(
  outer(
    c(tp1, tp2),
    metrics,
    FUN = function(tp, m) paste0("cell_text~~", tp, "|||", m)
  )
)
# Flatten in the right order: tp1+all metrics, then tp2+all metrics
wide_all <- bind_rows(wide_cont, wide_bin)

# Derive column order directly from actual column names present
cell_cols_ordered <- grep("^cell_text", names(wide_all), value = TRUE)
sig_cols_ordered  <- grep("^sig",       names(wide_all), value = TRUE)

wide_all <- wide_all %>%
  dplyr::select(pfas_type, PFAS, all_of(cell_cols_ordered), all_of(sig_cols_ordered))

# Pretty display column names (will be overridden by header spanning)
display_names <- c("PFAS Type", "PFAS", rep(metrics, 2))
names(wide_all)[1:ncol(wide_all)] <- make.names(names(wide_all))  # safety

# ── 7. Build flextable ────────────────────────────────────────────────────────
# Use actual column names as they appear in wide_all (dots, not special chars)
val_cols_actual <- grep("^cell_text", names(wide_all), value = TRUE)
sig_cols_actual <- grep("^sig",       names(wide_all), value = TRUE)

wide_display <- wide_all
names(wide_display)[match(val_cols_actual, names(wide_display))] <- paste0("val_", seq_along(val_cols_actual))
names(wide_display)[match(sig_cols_actual, names(wide_display))] <- paste0("sig_", seq_along(sig_cols_actual))

# Columns to show in table (hide sig_ columns — used only for formatting)
val_cols <- paste0("val_", 1:8)
show_cols <- c("pfas_type", "PFAS", val_cols)

# # Add a spacer blank row between Continuous and Binary
# spacer_row <- wide_display[1, ]
# spacer_row[, ] <- ""
# spacer_row$pfas_type <- "___spacer___"
# # Fix sig columns back to logical NA (avoid type conflict in bind_rows)
# sig_cols_in_display <- grep("^sig", names(spacer_row), value = TRUE)
# spacer_row[, sig_cols_in_display] <- NA
# 
# wide_final <- bind_rows(
#   filter(wide_display, pfas_type == "Continuous"),
#   spacer_row,
#   filter(wide_display, pfas_type == "Binary")
# )
wide_final <- bind_rows(
  filter(wide_display, pfas_type == "Continuous"),
  filter(wide_display, pfas_type == "Binary")
)

# Build flextable on display columns only
ft_data <- wide_final[, show_cols]

# ft <- flextable(ft_data) %>%
ft <- flextable(ft_data) %>%
  merge_v(j = "pfas_type") %>%
  rotate(j = "pfas_type", rotation = "btlr", part = "body") %>%
  # ── Column headers ──────────────────────────────────────────────────────────
  set_header_labels(
    pfas_type = "",
    PFAS      = "PFAS",
    val_1 = "Richness", val_2 = "Evenness", val_3 = "Shannon", val_4 = "Simpson",
    val_5 = "Richness", val_6 = "Evenness", val_7 = "Shannon", val_8 = "Simpson"
  ) %>%
  # Spanning header row for timepoints
  add_header_row(
    values = c("", "", "1m PFAS → 1m Microbiome", "1m PFAS → 6m Microbiome"),
    colwidths = c(1, 1, 4, 4)
  ) %>%
  # ── Basic styling ───────────────────────────────────────────────────────────
  theme_booktabs() %>%
  fontsize(size = 8, part = "all") %>%
  font(fontname = "Arial", part = "all") %>%
  align(align = "center", part = "all") %>%
  align(j = c("pfas_type", "PFAS"), align = "left", part = "body") %>%
  bold(part = "header") %>%
  # ── Column widths ───────────────────────────────────────────────────────────
  width(j = "pfas_type", width = 0.45) %>%
  fontsize(j = "pfas_type", size = 12, part = "body") %>%
  width(j = "PFAS",      width = 0.95) %>%
  width(j = val_cols,    width = 0.85) %>%
  # ── Row heights ─────────────────────────────────────────────────────────────
  height_all(height = 0.45) %>%
  # ── Header shading ──────────────────────────────────────────────────────────
  bg(i = 1, bg = "black", part = "header") %>%  # top spanning row — dark
  bg(i = 2, bg = "grey", part = "header") %>%  # metric row — medium blue
  color(i = 1, color = "white", part = "header") %>%
  color(i = 2, color = "black", part = "header") %>%
  fontsize(i = 1, size = 11, part = "header") %>%  # top spanning row text size
  fontsize(i = 2, size = 9,  part = "header") %>%  # metric row text size
  # ── Section label rows (pfas_type col) ─────────────────────────────────────
  bold(j = "pfas_type", part = "body") %>%
  color(j = "pfas_type", color = "black", part = "body") %>%
  # ── Borders ─────────────────────────────────────────────────────────────────
  # border_outer(part = "all", border = fp_border(color = "black", width = 1.5)) %>%
  border_outer(part = "all", border = fp_border(color = "black", width = 1.5)) %>%
  border_inner_h(part = "body", border = fp_border(color = "#CCCCCC", width = 0.5)) %>%
  border_inner_v(part = "body", border = fp_border(color = "#CCCCCC", width = 0.5)) %>%
  hline(i = nrow(filter(wide_display, pfas_type == "Continuous")),
        border = fp_border(color = "black", width = 2), part = "body") %>%
  # Thick border between metric groups (after col 5 = val_4 / val_5 boundary)
  vline(j = "val_4", border = fp_border(color = "black", width = 1.5), part = "all") %>%
  vline(j = "PFAS",  border = fp_border(color = "black", width = 1.0), part = "all")

# ── 8. Apply bold to significant p-values ──────────────────────────────
# We iterate over val columns (1..8) and check corresponding sig_ column
  for (col_i in 1:8) {
    val_col <- paste0("val_", col_i)
    sig_col <- paste0("sig_", col_i)
  
  sig_vec <- wide_final[[sig_col]]
  
  # Convert to logical safely
  sig_logical <- suppressWarnings(as.logical(sig_vec))
  sig_logical[is.na(sig_logical)] <- FALSE
  
  sig_rows <- which(sig_logical & wide_final$pfas_type != "___spacer___")
  
  if (length(sig_rows) > 0) {
    ft <- bold(ft,   i = sig_rows, j = val_col, part = "body")
    ft <- color(ft,  i = sig_rows, j = val_col, color = "black", part = "body")
  }
}

# ── 10. Save to Word ──────────────────────────────────────────────────────────
doc <- read_docx() %>%
  body_end_section_landscape() %>%
  body_add_flextable(ft)

print(doc, target = "out_files/Table_AlphaDiversity_PFAS.docx")
cat("✓ Table saved to: out_files/Table_AlphaDiversity_PFAS.docx\n")



# Alpha Diversity Scatter Plots — Continuous PFAS
# Two figures: 1m→1m and 1m→6m
# Each: 4 rows (metrics) × 5 cols (PFAS), with regression line + 95% CI
#
# Reads the 4 merged analysis datasets produced in Script 4
# (data_cont_1m_1m and data_cont_1m_6m must be in environment, or re-run
#  prepare_data() calls below)

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(here)

# ── If running standalone, re-run prepare_data to get the merged datasets ─────
# source("4. Alpha Diversity Analyses.R")   # uncomment if needed
# Otherwise assumes data_cont_1m_1m / data_cont_1m_6m are already in environment

# ── PFAS and metric labels ────────────────────────────────────────────────────
pfas_labels <- c(
  PFBS  = "PFBS",
  PFHxS = "PFHxS",
  PFNA  = "PFNA",
  PFOA  = "PFOA",
  PFOS  = "PFOS"
)

metric_labels <- c(
  ShannonDiv_repeated_rare = "Shannon Diversity",
  Richness_repeated_rare   = "Richness",
  Simpson_repeated_rare    = "Simpson Index",
  Evenness_repeated_rare   = "Evenness"
)

pfas_vars    <- names(pfas_labels)
metric_vars  <- names(metric_labels)

# ── Helper: reshape one dataset to long format ────────────────────────────────
make_long <- function(data, scenario_label) {
  data %>%
    dplyr::select(all_of(pfas_vars), all_of(metric_vars)) %>%
    pivot_longer(cols = all_of(pfas_vars),
                 names_to  = "PFAS",
                 values_to = "pfas_val") %>%
    pivot_longer(cols = all_of(metric_vars),
                 names_to  = "metric",
                 values_to = "metric_val") %>%
    mutate(
      PFAS     = factor(PFAS,   levels = pfas_vars,   labels = pfas_labels),
      metric   = factor(metric, levels = metric_vars, labels = metric_labels),
      scenario = scenario_label
    )
}

# ── Re-create merged datasets from saved CSVs ─────────────────────────────────
PFAS1m_micro1m <- read.csv(here::here("out_files", "PFAS1m_micro1m_species.csv")) %>% dplyr::select(-X)
PFAS1m_micro6m <- read.csv(here::here("out_files", "PFAS1m_micro6m_species.csv")) %>% dplyr::select(-X)

pfas_rename_continuous <- c(
  "PFBS"  = "PFBS_pgmL",
  "PFHxS" = "PFHxS_pgmL",
  "PFNA"  = "PFNA_pgmL",
  "PFOA"  = "PFOA_pgmL",
  "PFOS"  = "PFOS_pgmL"
)

alpha_1m <- alpha_div %>% filter(timepoint == 1)
alpha_6m <- alpha_div %>% filter(timepoint == 6)

data_cont_1m_1m <- PFAS1m_micro1m %>%
  dplyr::rename(any_of(pfas_rename_continuous)) %>%
  left_join(alpha_1m, by = "merge_id_dyad")

data_cont_1m_6m <- PFAS1m_micro6m %>%
  dplyr::rename(any_of(pfas_rename_continuous)) %>%
  left_join(alpha_6m, by = "merge_id_dyad")

long_1m_1m <- make_long(data_cont_1m_1m, "1m PFAS \u2192 1m Microbiome")
long_1m_6m <- make_long(data_cont_1m_6m, "1m PFAS \u2192 6m Microbiome")

# ── Plot theme ────────────────────────────────────────────────────────────────
scatter_theme <- theme_classic(base_size = 10) +
  theme(
    strip.text.x       = element_text(face = "bold", size = 9, color = "white"),
    strip.text.y       = element_text(face = "bold", size = 9, color = "white"),
    strip.background   = element_rect(fill = "#2E4057", color = NA),
    axis.title         = element_text(size = 9),
    axis.text          = element_text(size = 7, color = "black"),
    panel.border       = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
    panel.spacing      = unit(0.6, "lines"),
    plot.title         = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle      = element_text(size = 8, hjust = 0.5, color = "grey40"),
    legend.position    = "none"
  )

# ── Plot function ─────────────────────────────────────────────────────────────
make_scatter_metric <- function(long_data, metric_label, title_label) {
  long_data %>%
    filter(metric == metric_label) %>%
    ggplot(aes(x = pfas_val, y = metric_val)) +
    geom_point(alpha = 0.35, size = 1.2, color = "#4A6FA5") +
    geom_smooth(method = "lm", se = TRUE,
                color  = "#C0392B", fill = "#F8D7D7",
                linewidth = 0.7, alpha = 0.25) +
    facet_wrap(~ PFAS, nrow = 1, scales = "free_x") +
    labs(
      title = paste0(title_label, " | ", metric_label),
      x     = "PFAS concentration (log2 pg/mL)",
      y     = metric_label
    ) +
    scatter_theme
}

metrics_to_plot <- c("Shannon Diversity", "Richness", "Simpson Index", "Evenness")

for (metric_label in metrics_to_plot) {
  p1 <- make_scatter_metric(long_1m_1m, metric_label, "1-Month PFAS -> 1-Month Microbiome")
  p2 <- make_scatter_metric(long_1m_6m, metric_label, "1-Month PFAS -> 6-Month Microbiome")
  
  combined <- p1 / p2   # stacked: 1m on top, 6m on bottom
  
  metric_slug <- gsub(" ", "_", tolower(metric_label))
  ggsave(
    filename = here::here("out_figures", paste0("scatter_alphaDiv_", metric_slug, ".pdf")),
    plot     = combined,
    width    = 12, height = 8, dpi = 600, bg = "white"
  )
  cat("✓ Saved:", paste0("scatter_alphaDiv_", metric_slug, ".pdf"), "\n")
}
