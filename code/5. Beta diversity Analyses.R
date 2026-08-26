# header -----------------------------------------------------------------------
#
# TITLE:   4. Beta Diversity Analyses.R
#
# PURPOSE: Single-exposure (unadjusted) PERMANOVA and mixture analyses for
#          1-month PFAS associations with beta diversity at 1m and 6m
#
# DATE:    February 2026
#
#
# ANALYSIS STRUCTURE:
#    Continuously measured PFAS
#     (1) Single-exposure PERMANOVA per compound → R2 and p-value (unadjusted)
#     (2) qgcomp.boot on PCO1 and PCO2 from Bray-Curtis PCoA → Psi [95% CI]
#         PCoA converts the distance matrix into per-sample scalars so that
#         qgcomp.boot (gaussian) applies directly
#         Two mixture models run per outcome:
#           - Raw (q=4): PFAS log2-transformed, quantile-binned into quartiles
#             Psi = change per 1-quartile increase across all PFAS simultaneously
#           - IQR-scaled (q=NULL): each log2-PFAS divided by its IQR prior to
#             modeling so Psi = change per simultaneous 1-IQR increase across
#             all PFAS — directly comparable to IQR-scaled individual betas
#         B = 5000 bootstrap iterations, seed = 2024
#     Output: PERMANOVA CSV + qgcomp CSV per scenario; combined Excel workbook
#
#   Binary classified PFAS
#     (1) Single-exposure PERMANOVA per compound → R2 and p-value (unadjusted)
#     (2) n_detect PERMANOVA → R2 and p-value
#         n_detect is a continuous count (0-5), single predictor, so plain
#         PERMANOVA is sufficient — no qgcomp needed for a single variable.
#         IQR scaling not applied to binary PFAS — no logical IQR for 0/1 vars
#     Output: PERMANOVA CSV per scenario; combined Excel workbook
#
# All PERMANOVA models are unadjusted (standard practice for beta diversity).
# qgcomp models are also unadjusted — covariates not included for beta diversity.
#
# # CODE REVIEW:
# Reviewed by Ellie Holzhausen (EAH) on April 23, 2026
# Haonan Li (HL) on May 25, 2026
#
# set up -----------------------------------------------------------------------
rm(list = ls())
options(scipen = 0)

library(tidyr); library(reshape); library(dplyr)
library(purrr); library(stringr); library(readxl)
library(here); library(vegan); library(writexl)
library(tibble); library(qgcomp)
library(openxlsx)

# Read data files and process---------------------------------------------------
# Read beta diversity matrix
beta_div <- read.table(
  here::here("input", "betaDiv_bray_distance_repeated_rarefied_100_data_readDepth_1000000_mothersMilk_replacementFALSE.tsv"),
  sep = '\t', row.names = 1
)
colnames(beta_div) <- t(beta_div[1, ])
beta_div           <- beta_div[-1, ]
for (i in 1:ncol(beta_div)) beta_div[, i] <- as.numeric(beta_div[, i])
# Safety net: coerce any failed numeric conversions to 0 (none expected)
beta_div[is.na(beta_div)] <- 0
beta_div <- beta_div + t(beta_div)

# Read and process IDs
IDs <- read_excel(here::here("input", "GORAN_MICROBIOME_MANIFEST 10-2022.xlsx"))
colnames(IDs)     <- IDs[1, ]
IDs               <- IDs[-1, ]
IDs$dyad_id       <- substr(IDs$`Old Together`, 9, 11)
IDs$timepoint     <- substr(IDs$`Old Together`, 1, 2)
IDs$merge_id_dyad <- paste0("MM-", str_pad(IDs$dyad_id, width = 4, pad = "0"),
                            "-", IDs$timepoint)
IDs$`CORE ID`     <- paste0("MG", IDs$`CORE ID`)
colnames(IDs)[1]  <- "merge_id"
IDs               <- dplyr::select(IDs, c("merge_id", "merge_id_dyad"))

# Replace MG IDs with merge_id_dyad
id_lookup    <- setNames(IDs$merge_id_dyad, IDs$merge_id)
rownames_new <- id_lookup[rownames(beta_div)]
colnames_new <- id_lookup[colnames(beta_div)]
rownames_new[is.na(rownames_new)] <- rownames(beta_div)[is.na(rownames_new)]
colnames_new[is.na(colnames_new)] <- colnames(beta_div)[is.na(colnames_new)]
rownames(beta_div) <- rownames_new
colnames(beta_div) <- colnames_new
valid_samples      <- IDs$merge_id_dyad
beta_div <- beta_div[rownames(beta_div) %in% valid_samples,
                     colnames(beta_div) %in% valid_samples]

# Read PFAS metadata files
PFAS1m_micro1m       <- read.csv(here::here("out_files", "PFAS1m_micro1m_species.csv"))       %>% dplyr::select(-X)
PFAS1m_micro6m       <- read.csv(here::here("out_files", "PFAS1m_micro6m_species.csv"))       %>% dplyr::select(-X)
PFAS1mDetect_micro1m <- read.csv(here::here("out_files", "PFAS1mDetect_micro1m_species.csv")) %>% dplyr::select(-X)
PFAS1mDetect_micro6m <- read.csv(here::here("out_files", "PFAS1mDetect_micro6m_species.csv")) %>% dplyr::select(-X)

# Variable lists
pfas_rename_continuous <- c(
  "PFBS_pgmL"  = "PFBS",
  "PFHxS_pgmL" = "PFHxS",
  "PFNA_pgmL"  = "PFNA",
  "PFOA_pgmL"  = "PFOA",
  "PFOS_pgmL"  = "PFOS"
)

pfas_vars_continuous <- c("PFBS", "PFHxS", "PFNA", "PFOA", "PFOS")
pfas_vars_binary     <- c("N.MeFOSAA", "PFBA", "PFDA", "PFDoA", "PFHpA")

# Rename continuous PFAS columns
PFAS1m_micro1m <- PFAS1m_micro1m %>%
  dplyr::rename(!!!setNames(names(pfas_rename_continuous), pfas_rename_continuous))
PFAS1m_micro6m <- PFAS1m_micro6m %>%
  dplyr::rename(!!!setNames(names(pfas_rename_continuous), pfas_rename_continuous))

# Rename binary PFAS columns and set reference level
detect_rename <- setNames(paste0(pfas_vars_binary, "_Detect"), pfas_vars_binary)

PFAS1mDetect_micro1m <- PFAS1mDetect_micro1m %>%
  dplyr::rename(!!!detect_rename) %>%
  mutate(across(all_of(pfas_vars_binary),
                ~ factor(., levels = c("non-detect", "detect"))))

PFAS1mDetect_micro6m <- PFAS1mDetect_micro6m %>%
  dplyr::rename(!!!detect_rename) %>%
  mutate(across(all_of(pfas_vars_binary),
                ~ factor(., levels = c("non-detect", "detect"))))

# Compute IQR of each log2-transformed continuous PFAS-------------------------
compute_pfas_iqr <- function(data) {
  sapply(pfas_vars_continuous, function(pfas) {
    IQR(data[[pfas]], na.rm = TRUE)
  })
}

# FUNCTION: prepare_beta_data---------------------------------------------------
# Subsets beta matrix and metadata to shared IDs, removes rows with missing
# PFAS values, and extracts PCO1/PCO2 via PCoA.
# extract_pco = TRUE  for continuous scenarios (PCO scores needed for qgcomp)

prepare_beta_data <- function(metadata, beta_matrix, pfas_vars, scenario_name,
                              extract_pco = FALSE) {
  
  shared_ids <- intersect(metadata$merge_id_dyad, rownames(beta_matrix))
  cat("Shared samples:", length(shared_ids), "\n")
  if (length(shared_ids) == 0) stop("No shared samples!")
  
  beta_sub <- beta_matrix[shared_ids, shared_ids]
  meta_sub <- metadata[metadata$merge_id_dyad %in% shared_ids, ]
  meta_sub <- meta_sub[match(rownames(beta_sub), meta_sub$merge_id_dyad), ]
  
  if (!all(rownames(beta_sub) == meta_sub$merge_id_dyad)) {
    stop("Order mismatch between beta matrix and metadata!")
  }
  
  # Remove rows with missing PFAS values
  complete_idx <- complete.cases(meta_sub[, pfas_vars])
  n_before     <- nrow(meta_sub)
  meta_sub     <- meta_sub[complete_idx, ]
  beta_sub     <- beta_sub[complete_idx, complete_idx]
  cat("Samples removed due to missing PFAS:", n_before - nrow(meta_sub), "\n")
  cat("Final N:", nrow(meta_sub), "\n")
  
  var_pco1 <- NULL
  var_pco2 <- NULL
  
  # PCoA only for continuous scenarios — converts distance matrix to per-sample
  # scalars (PCO1, PCO2) so qgcomp.boot can be applied as an outcome model
  if (extract_pco) {
    pcoa_result <- cmdscale(as.dist(beta_sub), k = 2, eig = TRUE)
    pos_eig     <- pcoa_result$eig[pcoa_result$eig > 0]
    var_pco1    <- round(pcoa_result$eig[1] / sum(pos_eig) * 100, 1)
    var_pco2    <- round(pcoa_result$eig[2] / sum(pos_eig) * 100, 1)
    cat("  PCO1 explains:", var_pco1, "% of Bray-Curtis variance\n")
    cat("  PCO2 explains:", var_pco2, "% of Bray-Curtis variance\n")
    meta_sub$PCO1 <- pcoa_result$points[, 1]
    meta_sub$PCO2 <- pcoa_result$points[, 2]
  }
  
  return(list(
    beta     = beta_sub,
    metadata = meta_sub,
    var_pco1 = var_pco1,
    var_pco2 = var_pco2
  ))
}

# FUNCTION: run_permanova_single------------------------------------------------
# Unadjusted PERMANOVA, one PFAS at a time. Reports R2 and p-value.
# Used for both continuous and binary PFAS.
run_permanova_single <- function(beta_matrix, metadata, pfas_list,
                                 scenario_name) {
  results <- list()
  set.seed(2024)
  for (pfas in pfas_list) {
    if (!pfas %in% names(metadata)) {
      cat("  Warning:", pfas, "not found. Skipping.\n"); next
    }
    result <- adonis2(
      as.formula(paste0("beta_matrix ~ ", pfas)),
      data         = metadata,
      method       = "bray",
      by           = "margin",
      permutations = 9999
    )
    results[[pfas]] <- result
    cat("  Completed:", pfas, "\n")
  }
  return(results)
}


# Check the distirbution of n_detect values
print(table(PFAS1mDetect_micro1m$n_detect))
print(table(PFAS1mDetect_micro6m$n_detect))
# Good spread  across the range. Distribution looks fair and reasonable for PERMANOVA

# FUNCTION: run_permanova_ndetect-----------------------------------------------
# n_detect (continuous count 0-5) as a single numeric predictor.
# Reports R2 and p-value. No qgcomp needed — single variable.

run_permanova_ndetect <- function(beta_matrix, metadata, scenario_name) {
  
  cat("\n===== n_detect PERMANOVA:", scenario_name, "=====\n")
  if (!"n_detect" %in% names(metadata)) stop("n_detect not found!")
  
  metadata$n_detect <- as.numeric(metadata$n_detect)
  cat("  n_detect range:", min(metadata$n_detect, na.rm = TRUE),
      "to", max(metadata$n_detect, na.rm = TRUE), "\n")
  print(table(metadata$n_detect))
  
  set.seed(2024)
  result <- adonis2(
    beta_matrix ~ n_detect,
    data         = metadata,
    method       = "bray",
    by           = "margin",
    permutations = 9999
  )
  return(result)
}

# FUNCTION: run_qgcomp_pco_combined---------------------------------------------
# qgcomp.boot with gaussian() on PCO1 and PCO2 as outcomes.
# Approach: PCoA converts the Bray-Curtis distance matrix into per-sample
# coordinate scores (PCO1, PCO2) — continuous scalars
# B = 5000 : bootstrapped CIs, consistent with alpha diversity qgcomp
# Runs both log2 transformed (q=4) and IQR-scaled (q=NULL) qgcomp on PCO1 and PCO2
# Returns transformed Psi (Beta), IQR-scaled Psi (Beta_IQR)

run_qgcomp_pco_combined <- function(metadata, pfas_vars, pfas_iqr,
                                    scenario_name, B = 5000) {
  
  complete_idx <- complete.cases(metadata[, c(pfas_vars, "PCO1", "PCO2")])
  meta         <- metadata[complete_idx, ]
  cat("  N:", nrow(meta), "\n")
  
  # IQR-scaled copy
  meta_iqr <- meta
  for (pfas in pfas_vars_continuous) {
    meta_iqr[[pfas]] <- meta_iqr[[pfas]] / pfas_iqr[pfas]
  }
  
  results <- list()
  set.seed(2024)
  for (outcome_var in c("PCO1", "PCO2")) {
    
    cat("  Running outcome:", outcome_var, "...\n")
    
    formula_qg <- as.formula(paste(outcome_var, "~",
                                   paste(pfas_vars, collapse = " + ")))
    
    # Raw model (q=4)
    fit_raw <- qgcomp.boot(
      f = formula_qg, data = meta, expnms = pfas_vars,
      family = gaussian(), q = 4, B = B, seed = 2024, rr = FALSE
    )
    
    # IQR-scaled model (q=NULL)
    fit_iqr <- qgcomp.boot(
      f = formula_qg, data = meta_iqr, expnms = pfas_vars,
      family = gaussian(), q = NULL, B = B, seed = 2024, rr = FALSE
    )
    
    results[[outcome_var]] <- list(
      fit_raw = fit_raw,
      fit_iqr = fit_iqr
    )
    
    cat("  Raw Psi:", round(fit_raw$coef["psi1"], 5),
        "| IQR Psi:", round(fit_iqr$coef["psi1"], 5),
        "| p:", round(fit_raw$pval[2], 4), "\n")
  }
  return(results)
}


# EXTRACTION: single-exposure PERMANOVA-----------------------------------------
extract_single_results <- function(results_list, scenario_name) {
  do.call(rbind, lapply(names(results_list), function(pfas) {
    res <- as.data.frame(results_list[[pfas]])
    data.frame(
      PFAS        = pfas,
      R2          = round(res$R2[1], 4),
      F_statistic = round(res$F[1], 4),
      P           = round(res$`Pr(>F)`[1], 4),
      Scenario    = scenario_name,
      Type        = "Single-exposure (PERMANOVA)",
      stringsAsFactors = FALSE
    )
  }))
}


# EXTRACTION: n_detect PERMANOVA------------------------------------------------
extract_ndetect_permanova_result <- function(result, scenario_name) {
  res <- as.data.frame(result)
  data.frame(
    PFAS        = "N-detect",
    R2          = round(res$R2[1], 4),
    F_statistic = round(res$F[1], 4),
    P           = round(res$`Pr(>F)`[1], 4),
    Scenario    = scenario_name,
    Type        = "Mixture (N-detect PERMANOVA)",
    stringsAsFactors = FALSE
  )
}


# EXTRACTION: qgcomp on PCO scores----------------------------------------------
# One mixture row per outcome (PCO1 and PCO2)
extract_qgcomp_pco_results <- function(qgcomp_output,
                                       scenario_name, var_pco1, var_pco2) {
  
  pco_var_explained <- c(PCO1 = var_pco1, PCO2 = var_pco2)
  rows <- list()
  
  for (outcome_var in names(qgcomp_output)) {
    
    fit_raw       <- qgcomp_output[[outcome_var]]$fit_raw
    fit_iqr       <- qgcomp_output[[outcome_var]]$fit_iqr
    outcome_label <- paste0(outcome_var, " (", pco_var_explained[outcome_var], "% var)")
    
    rows[[paste(outcome_var, "Mixture")]] <- data.frame(
      Outcome  = outcome_label,
      Predictor = "Mixture (Psi)",
      Estimate  = round(fit_raw$coef["psi1"], 5),
      CI_lo     = round(fit_raw$ci[1], 5),
      CI_hi     = round(fit_raw$ci[2], 5),
      P         = round(fit_raw$pval[2], 4),
      Estimate_IQR = round(fit_iqr$coef["psi1"], 5),
      CI_lo_IQR    = round(fit_iqr$ci[1], 5),
      CI_hi_IQR    = round(fit_iqr$ci[2], 5),
      N         = nrow(fit_raw$fit$model),
      Scenario  = scenario_name,
      Type      = "Mixture (qgcomp.boot on PCO)",
      stringsAsFactors = FALSE
    )
  }
  return(do.call(rbind, rows))
}

# add significance flag
add_sig <- function(df) mutate(df, Significant = ifelse(P < 0.05, "Yes", "No"))

# MAIN ANALYSIS-----------------------------------------------------------------
# Scenario 1: Continuous 1m PFAS + 1m beta diversity
data_cont_1m_1m <- prepare_beta_data(
  PFAS1m_micro1m, beta_div, pfas_vars_continuous,
  "Continuous 1m PFAS + 1m Microbiome",
  extract_pco = TRUE
)
single_cont_1m_1m     <- run_permanova_single(
  data_cont_1m_1m$beta, data_cont_1m_1m$metadata,
  pfas_vars_continuous, "Continuous 1m PFAS + 1m Microbiome"
)
pfas_iqr_1m           <- compute_pfas_iqr(data_cont_1m_1m$metadata)
cat("PFAS IQRs (1m data):\n"); print(round(pfas_iqr_1m, 4))

qgcomp_pco_cont_1m_1m <- run_qgcomp_pco_combined(
  data_cont_1m_1m$metadata, pfas_vars_continuous, pfas_iqr_1m,
  "Continuous 1m PFAS + 1m Microbiome"
)

output_cont_1m_1m_permanova <- add_sig(
  extract_single_results(single_cont_1m_1m, "1m PFAS + 1m Microbiome")
)
output_cont_1m_1m_qgcomp <- extract_qgcomp_pco_results(
  qgcomp_pco_cont_1m_1m,
  "1m PFAS + 1m Microbiome",
  data_cont_1m_1m$var_pco1, data_cont_1m_1m$var_pco2
)


# Scenario 2: Continuous 1m PFAS + 6m beta diversity
data_cont_1m_6m <- prepare_beta_data(
  PFAS1m_micro6m, beta_div, pfas_vars_continuous,
  "Continuous 1m PFAS + 6m Microbiome",
  extract_pco = TRUE
)
single_cont_1m_6m     <- run_permanova_single(
  data_cont_1m_6m$beta, data_cont_1m_6m$metadata,
  pfas_vars_continuous, "Continuous 1m PFAS + 6m Microbiome"
)
pfas_iqr_6m           <- compute_pfas_iqr(data_cont_1m_6m$metadata)
cat("PFAS IQRs (6m data):\n"); print(round(pfas_iqr_6m, 4))

qgcomp_pco_cont_1m_6m <- run_qgcomp_pco_combined(
  data_cont_1m_6m$metadata, pfas_vars_continuous, pfas_iqr_6m,
  "Continuous 1m PFAS + 6m Microbiome"
)
output_cont_1m_6m_permanova <- add_sig(
  extract_single_results(single_cont_1m_6m, "1m PFAS + 6m Microbiome")
)
output_cont_1m_6m_qgcomp <- extract_qgcomp_pco_results(
  qgcomp_pco_cont_1m_6m,
  "1m PFAS + 6m Microbiome",
  data_cont_1m_6m$var_pco1, data_cont_1m_6m$var_pco2
)


# Scenario 3: Binary 1m PFAS + 1m beta diversity
data_bin_1m_1m <- prepare_beta_data(
  PFAS1mDetect_micro1m, beta_div, pfas_vars_binary,
  "Binary 1m PFAS + 1m Microbiome",
  extract_pco = FALSE
)
single_bin_1m_1m   <- run_permanova_single(
  data_bin_1m_1m$beta, data_bin_1m_1m$metadata,
  pfas_vars_binary, "Binary 1m PFAS + 1m Microbiome"
)
ndetect_perm_1m_1m <- run_permanova_ndetect(
  data_bin_1m_1m$beta, data_bin_1m_1m$metadata,
  "Binary 1m PFAS + 1m Microbiome"
)
output_bin_1m_1m_permanova <- add_sig(rbind(
  extract_single_results(single_bin_1m_1m, "1m PFAS + 1m Microbiome"),
  extract_ndetect_permanova_result(ndetect_perm_1m_1m, "1m PFAS + 1m Microbiome")
))


# Scenario 4: Binary 1m PFAS + 6m beta diversity
data_bin_1m_6m <- prepare_beta_data(
  PFAS1mDetect_micro6m, beta_div, pfas_vars_binary,
  "Binary 1m PFAS + 6m Microbiome",
  extract_pco = FALSE
)
single_bin_1m_6m   <- run_permanova_single(
  data_bin_1m_6m$beta, data_bin_1m_6m$metadata,
  pfas_vars_binary, "Binary 1m PFAS + 6m Microbiome"
)
ndetect_perm_1m_6m <- run_permanova_ndetect(
  data_bin_1m_6m$beta, data_bin_1m_6m$metadata,
  "Binary 1m PFAS + 6m Microbiome"
)
output_bin_1m_6m_permanova <- add_sig(rbind(
  extract_single_results(single_bin_1m_6m, "1m PFAS + 6m Microbiome"),
  extract_ndetect_permanova_result(ndetect_perm_1m_6m, "1m PFAS + 6m Microbiome")
))


# SAVE RESULTS------------------------------------------------------------------
write.csv(output_cont_1m_1m_permanova, here::here("out_files", "betaDiv_continuous_1m_1m_PERMANOVA.csv"),  row.names = FALSE)
write.csv(output_cont_1m_1m_qgcomp,   here::here("out_files", "betaDiv_continuous_1m_1m_qgcompPCO.csv"),  row.names = FALSE)

write.csv(output_cont_1m_6m_permanova, here::here("out_files", "betaDiv_continuous_1m_6m_PERMANOVA.csv"),  row.names = FALSE)
write.csv(output_cont_1m_6m_qgcomp,   here::here("out_files", "betaDiv_continuous_1m_6m_qgcompPCO.csv"),  row.names = FALSE)

write.csv(output_bin_1m_1m_permanova, here::here("out_files", "betaDiv_binary_1m_1m_PERMANOVA.csv"), row.names = FALSE)
write.csv(output_bin_1m_6m_permanova, here::here("out_files", "betaDiv_binary_1m_6m_PERMANOVA.csv"), row.names = FALSE)



# EXPORT BETA DIVERSITY RESULTS TO EXCEL TOGETHER-------------------------------
# Format helpers
fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

fmt_r2 <- function(r2, p) {
  ifelse(is.na(r2), NA_character_,
         ifelse(p < 0.05,
                paste0(sprintf("%.4f", r2), "*"),
                sprintf("%.4f", r2)))
}

fmt_psi <- function(est, lo, hi, p) {
  ifelse(is.na(est), NA_character_,
         ifelse(p < 0.05,
                paste0(sprintf("%.4f", est), " [",
                       sprintf("%.4f", lo), ", ",
                       sprintf("%.4f", hi), "]*"),
                paste0(sprintf("%.4f", est), " [",
                       sprintf("%.4f", lo), ", ",
                       sprintf("%.4f", hi), "]")))
}

# Section 1: Continuous PFAS — Single-exposure PERMANOVA
pfas_order_cont <- c("PFBS", "PFHxS", "PFNA", "PFOA", "PFOS")

perm_11 <- output_cont_1m_1m_permanova %>% 
  filter(PFAS %in% pfas_order_cont) %>%
  mutate(PFAS = factor(PFAS, levels = pfas_order_cont)) %>%
  arrange(PFAS)

perm_16 <- output_cont_1m_6m_permanova %>%
  filter(PFAS %in% pfas_order_cont) %>%
  mutate(PFAS = factor(PFAS, levels = pfas_order_cont)) %>%
  arrange(PFAS)

cont_perm_section <- data.frame(
  Exposure             = pfas_order_cont,
  `R2 (1m Microbiome)` = fmt_r2(perm_11$R2, perm_11$P),
  `P (1m Microbiome)`  = fmt_p(perm_11$P),
  `R2 (6m Microbiome)` = fmt_r2(perm_16$R2, perm_16$P),
  `P (6m Microbiome)`  = fmt_p(perm_16$P),
  check.names = FALSE, stringsAsFactors = FALSE
)

# Section 2: Continuous PFAS — qgcomp PCO mixture
qpco_11 <- output_cont_1m_1m_qgcomp
qpco_16 <- output_cont_1m_6m_qgcomp

cont_qgcomp_section <- data.frame(
  Exposure                       = gsub(" \\(.*\\)", "", qpco_11$Outcome),
  `Psi [95CI] (1m Microbiome)`   = fmt_psi(qpco_11$Estimate, qpco_11$CI_lo, qpco_11$CI_hi, qpco_11$P),
  `Psi_IQR [95CI] (1m Microbiome)` = fmt_psi(qpco_11$Estimate_IQR, qpco_11$CI_lo_IQR, qpco_11$CI_hi_IQR, qpco_11$P),
  `P (1m Microbiome)`            = fmt_p(qpco_11$P),
  `Psi [95CI] (6m Microbiome)`   = fmt_psi(qpco_16$Estimate, qpco_16$CI_lo, qpco_16$CI_hi, qpco_16$P),
  `Psi_IQR [95CI] (6m Microbiome)` = fmt_psi(qpco_16$Estimate_IQR, qpco_16$CI_lo_IQR, qpco_16$CI_hi_IQR, qpco_16$P),
  `P (6m Microbiome)`            = fmt_p(qpco_16$P),
  check.names = FALSE, stringsAsFactors = FALSE
)

# Rename qgcomp columns to match the rest before rbind
names(cont_qgcomp_section) <- names(cont_perm_section)

# Section 3: Binary PFAS — Single-exposure PERMANOVA
pfas_order_bin <- c("N-MeFOSAA", "PFBA", "PFDA", "PFDoA", "PFHpA")

bin_single_11 <- output_bin_1m_1m_permanova %>%
  filter(Type == "Single-exposure (PERMANOVA)") %>%
  mutate(PFAS = recode(PFAS, "N.MeFOSAA" = "N-MeFOSAA")) %>%
  mutate(PFAS = factor(PFAS, levels = pfas_order_bin)) %>%
  arrange(PFAS)

bin_single_16 <- output_bin_1m_6m_permanova %>%
  filter(Type == "Single-exposure (PERMANOVA)") %>%
  mutate(PFAS = recode(PFAS, "N.MeFOSAA" = "N-MeFOSAA")) %>%
  mutate(PFAS = factor(PFAS, levels = pfas_order_bin)) %>%
  arrange(PFAS)

bin_single_section <- data.frame(
  Exposure             = as.character(bin_single_11$PFAS),
  `R2 (1m Microbiome)` = fmt_r2(bin_single_11$R2, bin_single_11$P),
  `P (1m Microbiome)`  = fmt_p(bin_single_11$P),
  `R2 (6m Microbiome)` = fmt_r2(bin_single_16$R2, bin_single_16$P),
  `P (6m Microbiome)`  = fmt_p(bin_single_16$P),
  check.names = FALSE, stringsAsFactors = FALSE
)

# Section 4: Binary PFAS — N-detect PERMANOVA
bin_ndet_11 <- output_bin_1m_1m_permanova %>% filter(Type == "Mixture (N-detect PERMANOVA)")
bin_ndet_16 <- output_bin_1m_6m_permanova %>% filter(Type == "Mixture (N-detect PERMANOVA)")

bin_ndet_section <- data.frame(
  Exposure             = "N-detect",
  `R2 (1m Microbiome)` = fmt_r2(bin_ndet_11$R2, bin_ndet_11$P),
  `P (1m Microbiome)`  = fmt_p(bin_ndet_11$P),
  `R2 (6m Microbiome)` = fmt_r2(bin_ndet_16$R2, bin_ndet_16$P),
  `P (6m Microbiome)`  = fmt_p(bin_ndet_16$P),
  check.names = FALSE, stringsAsFactors = FALSE
)

# Combine all sections with labels
blank_row <- data.frame(
  Exposure = "", `R2 (1m Microbiome)` = "", `P (1m Microbiome)` = "",
  `R2 (6m Microbiome)` = "", `P (6m Microbiome)` = "",
  check.names = FALSE, stringsAsFactors = FALSE
)

label_row <- function(label) {
  data.frame(
    Exposure = label, `R2 (1m Microbiome)` = "", `P (1m Microbiome)` = "",
    `R2 (6m Microbiome)` = "", `P (6m Microbiome)` = "",
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

# Restore correct column names for qgcomp section
names(cont_qgcomp_section) <- c(
  "Exposure",
  "Psi [95% CI] (1m Microbiome)",
  "Psi_IQR [95% CI] (1m Microbiome)",
  "P (1m Microbiome)",
  "Psi [95% CI] (6m Microbiome)",
  "Psi_IQR [95% CI] (6m Microbiome)",
  "P (6m Microbiome)"
)

# Round all numeric-looking values to 2 decimal places in each section
round_section <- function(df) {
  df[] <- lapply(df, function(x) {
    if (is.character(x)) {
      sapply(x, function(val) {
        if (is.na(val) || val == "") return(val)
        # Find all numbers with 3+ decimal places and round them
        matches <- gregexpr("-?[0-9]+\\.[0-9]{3,}", val, perl = TRUE)
        regmatches(val, matches) <- lapply(regmatches(val, matches), function(nums) {
          sprintf("%.2f", as.numeric(nums))
        })
        val
      }, USE.NAMES = FALSE)
    } else {
      round(x, 2)
    }
  })
  return(df)
}

cont_perm_section    <- round_section(cont_perm_section)
cont_qgcomp_section  <- round_section(cont_qgcomp_section)
bin_single_section   <- round_section(bin_single_section)
bin_ndet_section     <- round_section(bin_ndet_section)

# Write two separate blocks with correct headers
wb <- createWorkbook()
addWorksheet(wb, "Beta Diversity Results")

# Block 1: Continuous PERMANOVA
writeData(wb, "Beta Diversity Results", 
          data.frame(Section = "--- Continuous PFAS: Single-Exposure PERMANOVA ---"),
          startRow = 1, startCol = 1, colNames = FALSE)
writeData(wb, "Beta Diversity Results", cont_perm_section, startRow = 2, startCol = 1)

# Block 2: Continuous qgcomp
writeData(wb, "Beta Diversity Results",
          data.frame(Section = "--- Continuous PFAS: Mixture (qgcomp.boot on PCoA Scores) ---"),
          startRow = 10, startCol = 1, colNames = FALSE)
writeData(wb, "Beta Diversity Results", cont_qgcomp_section, startRow = 11, startCol = 1)

# Block 3: Binary PERMANOVA
writeData(wb, "Beta Diversity Results",
          data.frame(Section = "--- Binary PFAS: Single-Exposure PERMANOVA ---"),
          startRow = 15, startCol = 1, colNames = FALSE)
writeData(wb, "Beta Diversity Results", bin_single_section, startRow = 16, startCol = 1)

# Block 4: Binary N-detect
writeData(wb, "Beta Diversity Results",
          data.frame(Section = "--- Binary PFAS: Mixture (N-detect PERMANOVA) ---"),
          startRow = 23, startCol = 1, colNames = FALSE)
writeData(wb, "Beta Diversity Results", bin_ndet_section, startRow = 24, startCol = 1)

setColWidths(wb, "Beta Diversity Results", cols = 1:7, widths = c(25, 25, 25, 12, 25, 25, 12))
saveWorkbook(wb, here::here("out_files", "betaDiv_results_table.xlsx"), overwrite = TRUE)

# END
