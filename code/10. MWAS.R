# TITLE:   2. MWAS.R
#
# PURPOSE: Metabolome-wide association study (MWAS) — 1-month human milk PFAS
#          and fecal metabolomics at 1m and 6m. Individual PFAS via lm(),
#          continuous mixture via qgcomp.glm.boot(), binary mixture via
#          n_detect lm(). Mirrors PFAS-microbiome pipeline.
#
# DATE:    May 2025
#
# NOTES:
#   - 8 scenarios: continuous/binary × 1m/6m × C18/HILIC
#   - Metabolites already log2 transformed from preprocessing
#   - PFAS already log2 transformed and MDL-imputed from cleaning script
#   - qgcomp: q=4, B=5000, seed=2024, gaussian family
#
# Code Review
# Ellie Holzhausen (EAH) on April 27, 2026
# Haonan Li (HL) on May 25, 2026
#
#-------------------------------------------------------------------------------

rm(list = ls())
options(scipen = 100)

library(tidyverse)
library(dplyr)
library(broom)
library(qgcomp)
library(parallel)
library(pheatmap)
library(ggplot2)
library(here)

# Parallelization
n_cores <- detectCores() - 2
N_BOOT  <- 5000

# Load saved files from Ellie annotated sensitivity analysis --------------------

# input_dir  <- "out_files/ellie_annotated_sensitivity"
# output_dir <- "out_files/ellie_annotated_sensitivity_mwas"
# 
# dir.create(here::here(output_dir), showWarnings = FALSE, recursive = TRUE)
# 
# PFAS1m_c18_1m          <- read.csv(here::here(input_dir, "PFAS1m_c18_1m.csv"),         check.names = FALSE)
# PFAS1m_c18_6m          <- read.csv(here::here(input_dir, "PFAS1m_c18_6m.csv"),         check.names = FALSE)
# PFAS1m_hilic_1m        <- read.csv(here::here(input_dir, "PFAS1m_hilic_1m.csv"),       check.names = FALSE)
# PFAS1m_hilic_6m        <- read.csv(here::here(input_dir, "PFAS1m_hilic_6m.csv"),       check.names = FALSE)
# 
# PFAS1mDetect_c18_1m    <- read.csv(here::here(input_dir, "PFAS1mDetect_c18_1m.csv"),   check.names = FALSE)
# PFAS1mDetect_c18_6m    <- read.csv(here::here(input_dir, "PFAS1mDetect_c18_6m.csv"),   check.names = FALSE)
# PFAS1mDetect_hilic_1m  <- read.csv(here::here(input_dir, "PFAS1mDetect_hilic_1m.csv"), check.names = FALSE)
# PFAS1mDetect_hilic_6m  <- read.csv(here::here(input_dir, "PFAS1mDetect_hilic_6m.csv"), check.names = FALSE)
# 
# keep_c18   <- readRDS(here::here(input_dir, "metabolite_cols_c18.rds"))
# keep_hilic <- readRDS(here::here(input_dir, "metabolite_cols_hilic.rds"))

# Load saved files from cleaning script-----------------------------------------
PFAS1m_c18_1m          <- read.csv(here::here("out_files", "PFAS1m_c18_1m.csv"),         check.names = FALSE)
PFAS1m_c18_6m          <- read.csv(here::here("out_files", "PFAS1m_c18_6m.csv"),         check.names = FALSE)
PFAS1m_hilic_1m        <- read.csv(here::here("out_files", "PFAS1m_hilic_1m.csv"),       check.names = FALSE)
PFAS1m_hilic_6m        <- read.csv(here::here("out_files", "PFAS1m_hilic_6m.csv"),       check.names = FALSE)
PFAS1mDetect_c18_1m    <- read.csv(here::here("out_files", "PFAS1mDetect_c18_1m.csv"),   check.names = FALSE)
PFAS1mDetect_c18_6m    <- read.csv(here::here("out_files", "PFAS1mDetect_c18_6m.csv"),   check.names = FALSE)
PFAS1mDetect_hilic_1m  <- read.csv(here::here("out_files", "PFAS1mDetect_hilic_1m.csv"), check.names = FALSE)
PFAS1mDetect_hilic_6m  <- read.csv(here::here("out_files", "PFAS1mDetect_hilic_6m.csv"), check.names = FALSE)

keep_c18   <- readRDS(here::here("out_files", "metabolite_cols_c18.rds"))
keep_hilic <- readRDS(here::here("out_files", "metabolite_cols_hilic.rds"))

# Convert detect columns separately — more reliable than inside mutate/across
convert_detect_cols <- function(df) {
  detect_cols <- grep("_Detect$", colnames(df), value = TRUE)
  for (col in detect_cols) {
    if (is.character(df[[col]])) {
      df[[col]] <- ifelse(trimws(df[[col]]) == "detect", 1L,
                          ifelse(trimws(df[[col]]) == "non-detect", 0L, NA_integer_))
    }
  }
  return(df)
}

# Apply immediately after read.csv in Section 1 — before prepare_mwas_data
PFAS1mDetect_c18_1m   <- convert_detect_cols(PFAS1mDetect_c18_1m)
PFAS1mDetect_c18_6m   <- convert_detect_cols(PFAS1mDetect_c18_6m)
PFAS1mDetect_hilic_1m <- convert_detect_cols(PFAS1mDetect_hilic_1m)
PFAS1mDetect_hilic_6m <- convert_detect_cols(PFAS1mDetect_hilic_6m)


# Variable lists----------------------------------------------------------------
pfas_vars_continuous <- c("PFBS_pgmL", "PFHxS_pgmL", "PFNA_pgmL",
                          "PFOA_pgmL", "PFOS_pgmL")

pfas_vars_binary     <- c("N.MeFOSAA_Detect", "PFBA_Detect", "PFDA_Detect",
                          "PFDoA_Detect", "PFHpA_Detect")

# Base covariates — gestational age dummies created inside prepare_mwas_data
covariates <- c("breastmilk_per_day", "SES_index_final", "baby_birthweight_kg",
                "gest_Early", "gest_Late", "mode_of_delivery_bin")


# Prepare_mwas_data-------------------------------------------------------------
# Creates gestational age dummies and mode of delivery binary
prepare_mwas_data <- function(df) {
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

PFAS1m_c18_1m         <- prepare_mwas_data(PFAS1m_c18_1m)
PFAS1m_c18_6m         <- prepare_mwas_data(PFAS1m_c18_6m)
PFAS1m_hilic_1m       <- prepare_mwas_data(PFAS1m_hilic_1m)
PFAS1m_hilic_6m       <- prepare_mwas_data(PFAS1m_hilic_6m)
PFAS1mDetect_c18_1m   <- prepare_mwas_data(PFAS1mDetect_c18_1m)
PFAS1mDetect_c18_6m   <- prepare_mwas_data(PFAS1mDetect_c18_6m)
PFAS1mDetect_hilic_1m <- prepare_mwas_data(PFAS1mDetect_hilic_1m)
PFAS1mDetect_hilic_6m <- prepare_mwas_data(PFAS1mDetect_hilic_6m)


# Compute PFAS IQRs for post-hoc beta scaling (continuous only)
# IQRs computed from 1m C18 dataset — same PFAS values across all scenarios
# since PFAS come from same 1m measurements
pfas_iqr_cont <- sapply(pfas_vars_continuous, function(p)
  IQR(PFAS1m_c18_1m[[p]], na.rm = TRUE))
cat("PFAS IQRs (continuous):\n")
print(round(pfas_iqr_cont, 4))

# Run_mwas: individual lm per PFAS × metabolite---------------------------------
run_mwas <- function(data, pfas_vars, metabolite_cols,
                     scenario_name, covs = covariates,
                     pfas_iqr = NULL) {
  
  cat("\n===== MWAS:", scenario_name, "=====\n")
  cat("  PFAS:", length(pfas_vars),
      "| Metabolites:", length(metabolite_cols), "\n")
  
  results_all <- list()
  
  for (pfas_var in pfas_vars) {
    
    cat("  Running:", pfas_var, "\n")
    results <- vector("list", length(metabolite_cols))
    
    for (i in seq_along(metabolite_cols)) {
      met          <- metabolite_cols[i]
      vars_needed  <- c(covs, pfas_var, met)
      ithData      <- na.omit(data[, vars_needed])
      
      tryCatch({
        formula  <- as.formula(
          paste0("`", met, "` ~ ", pfas_var, " + ",
                 paste(covs, collapse = " + "))
        )
        model    <- lm(formula, data = ithData)
        tidy_res <- broom::tidy(model)
        ci       <- confint(model)
        pfas_row <- tidy_res %>% filter(term == pfas_var)
        
        results[[i]] <- data.frame(
          Metabolite = met,
          Estimate   = pfas_row$estimate,
          CI_lo      = ci[pfas_var, 1],
          CI_hi      = ci[pfas_var, 2],
          P_value    = pfas_row$p.value,
          N          = nrow(ithData),
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        results[[i]] <<- data.frame(
          Metabolite = met, Estimate = NA_real_,
          CI_lo = NA_real_, CI_hi = NA_real_,
          P_value = NA_real_, N = NA_integer_,
          stringsAsFactors = FALSE
        )
      })
    }
    
    combined     <- bind_rows(results)
    combined$FDR <- p.adjust(combined$P_value, method = "BH")
    combined$PFAS <- pfas_var
    
    # Post-hoc IQR scaling — only if pfas_iqr provided (continuous scenarios)
    if (!is.null(pfas_iqr) && pfas_var %in% names(pfas_iqr)) {
      combined$Beta_IQR <- combined$Estimate * pfas_iqr[pfas_var]
      combined$CI_lo_IQR <- combined$CI_lo  * pfas_iqr[pfas_var]
      combined$CI_hi_IQR <- combined$CI_hi  * pfas_iqr[pfas_var]
    } else {
      combined$Beta_IQR  <- NA_real_
      combined$CI_lo_IQR <- NA_real_
      combined$CI_hi_IQR <- NA_real_
    }
    
    results_all[[pfas_var]] <- combined
  }
  
  bind_rows(results_all)
}


# Run_ndetect_mwas: n_detect mixture for binary PFAS
run_ndetect_mwas <- function(data, metabolite_cols,
                             scenario_name, covs = covariates) {
  
  results <- vector("list", length(metabolite_cols))
  
  for (i in seq_along(metabolite_cols)) {
    met         <- metabolite_cols[i]
    vars_needed <- c(covs, "n_detect", met)
    ithData     <- na.omit(data[, vars_needed])
    
    tryCatch({
      formula  <- as.formula(
        paste0("`", met, "` ~ n_detect + ",
               paste(covs, collapse = " + "))
      )
      model    <- lm(formula, data = ithData)
      tidy_res <- broom::tidy(model)
      ci       <- confint(model)
      nd_row   <- tidy_res %>% filter(term == "n_detect")
      
      results[[i]] <- data.frame(
        Metabolite = met,
        Estimate   = nd_row$estimate,
        CI_lo      = ci["n_detect", 1],
        CI_hi      = ci["n_detect", 2],
        P_value    = nd_row$p.value,
        N          = nrow(ithData),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      results[[i]] <<- data.frame(
        Metabolite = met, Estimate = NA_real_,
        CI_lo = NA_real_, CI_hi = NA_real_,
        P_value = NA_real_, N = NA_integer_,
        stringsAsFactors = FALSE
      )
    })
  }
  
  combined      <- bind_rows(results)
  combined$FDR  <- p.adjust(combined$P_value, method = "BH")
  combined$PFAS <- "N-detect"
  return(combined)
}


# Run_qgcomp_mwas: mixture model for continuous PFAS
# Parallelized across metabolites
# q=4, B=5000, seed=2024, gaussian family
run_qgcomp_mwas <- function(data, pfas_vars, metabolite_cols,
                            scenario_name, covs = covariates,
                            pfas_iqr = NULL,
                            B = N_BOOT, seed = 2024) {
  
  cat("\n===== qgcomp MWAS:", scenario_name, "=====\n")
  cat("  Metabolites:", length(metabolite_cols),
      "| B:", B, "| Cores:", n_cores, "\n")
  
  # IQR-scaled copy of data for q=NULL model
  data_iqr <- data
  if (!is.null(pfas_iqr)) {
    for (pfas in pfas_vars) {
      if (pfas %in% names(pfas_iqr))
        data_iqr[[pfas]] <- data_iqr[[pfas]] / pfas_iqr[pfas]
    }
  }
  
  results_list <- mclapply(seq_along(metabolite_cols), function(i) {
    
    met              <- metabolite_cols[i]
    vars_needed      <- c(covs, pfas_vars, met)
    ithData          <- na.omit(data[,     vars_needed])
    ithData_iqr      <- na.omit(data_iqr[, vars_needed])
    
    out <- list(
      Metabolite = met, Estimate = NA_real_,
      CI_lo = NA_real_, CI_hi = NA_real_,
      P_value = NA_real_, Beta_IQR = NA_real_,
      CI_lo_IQR = NA_real_, CI_hi_IQR = NA_real_,
      N = nrow(ithData)
    )
    
    if (nrow(ithData) < 20) return(out)
    
    formula <- as.formula(
      paste0("`", met, "` ~ ",
             paste(c(pfas_vars, covs), collapse = " + "))
    )
    
    tryCatch({
      # Raw model (q=4)
      qgmod_raw <- qgcomp.glm.boot(
        f      = formula,
        data   = ithData,
        expnms = pfas_vars,
        q      = 4,
        family = gaussian(),
        B      = B,
        seed   = seed,
        rr     = FALSE
      )
      # IQR-scaled model (q=NULL)
      qgmod_iqr <- qgcomp.glm.boot(
        f      = formula,
        data   = ithData_iqr,
        expnms = pfas_vars,
        q      = NULL,
        family = gaussian(),
        B      = B,
        seed   = seed,
        rr     = FALSE
      )
      out <- list(
        Metabolite = met,
        Estimate   = as.numeric(qgmod_raw$psi[1]),
        CI_lo      = as.numeric(qgmod_raw$ci[1]),
        CI_hi      = as.numeric(qgmod_raw$ci[2]),
        P_value    = as.numeric(qgmod_raw$pval[2]),
        Beta_IQR   = as.numeric(qgmod_iqr$psi[1]),
        CI_lo_IQR  = as.numeric(qgmod_iqr$ci[1]),
        CI_hi_IQR  = as.numeric(qgmod_iqr$ci[2]),
        N          = nrow(ithData)
      )
    }, error = function(e) {
      warning(paste("qgcomp failed for", met, ":", e$message))
    })
    
    return(out)
    
  }, mc.cores = n_cores)
  
  combined      <- bind_rows(lapply(results_list, as.data.frame))
  combined$FDR  <- p.adjust(combined$P_value, method = "BH")
  combined$PFAS <- "Mixture"
  
  cat("  Converged:", sum(!is.na(combined$Estimate)),
      "| Failed:", sum(is.na(combined$Estimate)), "\n")
  
  return(combined)
}


# MAIN ANALYSIS-----------------------------------------------------------------
# Scenario 1: Continuous 1m PFAS + 1m C18
mwas_cont_c18_1m   <- run_mwas(PFAS1m_c18_1m, pfas_vars_continuous,
                               keep_c18, "Cont 1m PFAS + C18 1m",
                               pfas_iqr = pfas_iqr_cont)
qgcomp_cont_c18_1m <- run_qgcomp_mwas(PFAS1m_c18_1m, pfas_vars_continuous,
                                      keep_c18, "Cont 1m PFAS + C18 1m",
                                      pfas_iqr = pfas_iqr_cont)
write.csv(bind_rows(mwas_cont_c18_1m, qgcomp_cont_c18_1m),
          here::here("out_files", "MWAS_continuous_c18_1m.csv"), row.names = FALSE)

# Scenario 2: Continuous 1m PFAS + 6m C18
mwas_cont_c18_6m   <- run_mwas(PFAS1m_c18_6m, pfas_vars_continuous,
                               keep_c18, "Cont 1m PFAS + C18 6m",
                               pfas_iqr = pfas_iqr_cont)
qgcomp_cont_c18_6m <- run_qgcomp_mwas(PFAS1m_c18_6m, pfas_vars_continuous,
                                      keep_c18, "Cont 1m PFAS + C18 6m",
                                      pfas_iqr = pfas_iqr_cont)
write.csv(bind_rows(mwas_cont_c18_6m, qgcomp_cont_c18_6m),
          here::here("out_files", "MWAS_continuous_c18_6m.csv"), row.names = FALSE)

# Scenario 3: Continuous 1m PFAS + 1m HILIC
mwas_cont_hilic_1m   <- run_mwas(PFAS1m_hilic_1m, pfas_vars_continuous,
                                 keep_hilic, "Cont 1m PFAS + HILIC 1m",
                                 pfas_iqr = pfas_iqr_cont)
qgcomp_cont_hilic_1m <- run_qgcomp_mwas(PFAS1m_hilic_1m, pfas_vars_continuous,
                                        keep_hilic, "Cont 1m PFAS + HILIC 1m",
                                        pfas_iqr = pfas_iqr_cont)
write.csv(bind_rows(mwas_cont_hilic_1m, qgcomp_cont_hilic_1m),
          here::here("out_files", "MWAS_continuous_hilic_1m.csv"), row.names = FALSE)

# Scenario 4: Continuous 1m PFAS + 6m HILIC
mwas_cont_hilic_6m   <- run_mwas(PFAS1m_hilic_6m, pfas_vars_continuous,
                                 keep_hilic, "Cont 1m PFAS + HILIC 6m",
                                 pfas_iqr = pfas_iqr_cont)
qgcomp_cont_hilic_6m <- run_qgcomp_mwas(PFAS1m_hilic_6m, pfas_vars_continuous,
                                        keep_hilic, "Cont 1m PFAS + HILIC 6m",
                                        pfas_iqr = pfas_iqr_cont)
write.csv(bind_rows(mwas_cont_hilic_6m, qgcomp_cont_hilic_6m),
          here::here("out_files", "MWAS_continuous_hilic_6m.csv"), row.names = FALSE)

# Scenario 5: Binary 1m PFAS + 1m C18
mwas_bin_c18_1m    <- run_mwas(PFAS1mDetect_c18_1m, pfas_vars_binary,
                               keep_c18, "Bin 1m PFAS + C18 1m")
ndetect_c18_1m     <- run_ndetect_mwas(PFAS1mDetect_c18_1m,
                                       keep_c18, "Bin 1m PFAS + C18 1m")
write.csv(bind_rows(mwas_bin_c18_1m, ndetect_c18_1m),
          here::here("out_files", "MWAS_binary_c18_1m.csv"), row.names = FALSE)

# Scenario 6: Binary 1m PFAS + 6m C18
mwas_bin_c18_6m    <- run_mwas(PFAS1mDetect_c18_6m, pfas_vars_binary,
                               keep_c18, "Bin 1m PFAS + C18 6m")
ndetect_c18_6m     <- run_ndetect_mwas(PFAS1mDetect_c18_6m,
                                       keep_c18, "Bin 1m PFAS + C18 6m")
write.csv(bind_rows(mwas_bin_c18_6m, ndetect_c18_6m),
          here::here("out_files", "MWAS_binary_c18_6m.csv"), row.names = FALSE)

# Scenario 7: Binary 1m PFAS + 1m HILIC
mwas_bin_hilic_1m  <- run_mwas(PFAS1mDetect_hilic_1m, pfas_vars_binary,
                               keep_hilic, "Bin 1m PFAS + HILIC 1m")
ndetect_hilic_1m   <- run_ndetect_mwas(PFAS1mDetect_hilic_1m,
                                       keep_hilic, "Bin 1m PFAS + HILIC 1m")
write.csv(bind_rows(mwas_bin_hilic_1m, ndetect_hilic_1m),
          here::here("out_files", "MWAS_binary_hilic_1m.csv"), row.names = FALSE)

# Scenario 8: Binary 1m PFAS + 6m HILIC
mwas_bin_hilic_6m  <- run_mwas(PFAS1mDetect_hilic_6m, pfas_vars_binary,
                               keep_hilic, "Bin 1m PFAS + HILIC 6m")
ndetect_hilic_6m   <- run_ndetect_mwas(PFAS1mDetect_hilic_6m,
                                       keep_hilic, "Bin 1m PFAS + HILIC 6m")
write.csv(bind_rows(mwas_bin_hilic_6m, ndetect_hilic_6m),
          here::here("out_files", "MWAS_binary_hilic_6m.csv"), row.names = FALSE)



# Quick summary of significant hits---------------------------------------------
result_files <- list(
  "Continuous C18 1m"   = "MWAS_continuous_c18_1m.csv",
  "Continuous C18 6m"   = "MWAS_continuous_c18_6m.csv",
  "Continuous HILIC 1m" = "MWAS_continuous_hilic_1m.csv",
  "Continuous HILIC 6m" = "MWAS_continuous_hilic_6m.csv",
  "Binary C18 1m"       = "MWAS_binary_c18_1m.csv",
  "Binary C18 6m"       = "MWAS_binary_c18_6m.csv",
  "Binary HILIC 1m"     = "MWAS_binary_hilic_1m.csv",
  "Binary HILIC 6m"     = "MWAS_binary_hilic_6m.csv"
)

for (nm in names(result_files)) {
  df  <- read.csv(here::here("out_files", result_files[[nm]]))
  sig <- sum(df$P_value < 0.05, na.rm = TRUE)
  fdr <- sum(df$FDR < 0.05, na.rm = TRUE)
  cat(nm, "— p<0.05:", sig, "| FDR<0.05:", fdr, "\n")
}


