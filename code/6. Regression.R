# header -----------------------------------------------------------------------
#
# TITLE:   6. Taxa Associations - CLR Transformation (All Taxonomic Levels)
#
# PURPOSE: CLR-transformed taxa associations with 1-month PFAS at 1m and 6m
#          microbiome. Uses lm() on CLR-transformed counts for individual PFAS
#          and qgcomp.glm.boot() with gaussian family for mixture analysis.
#          Runs across all taxonomic levels (species through phylum) and saves
#          combined results files for dendrogram plotting.
#
# Code Review:
#   Reviewed by Ellie Holzhausen (EAH) on April 27,2026
# by Haonan Li (HL) on May 25, 2026
#
# IQR SCALING:
#          To facilitate comparison of effect sizes across individual PFAS
#          compounds and with the mixture estimate, all log2-transformed PFAS
#          concentrations are divided by their respective interquartile ranges
#          (IQRs) prior to analysis or post-hoc scaled:
#            - Individual PFAS: raw betas multiplied post-hoc by each PFAS IQR
#              (Beta_IQR = beta x IQR). P-values unchanged.
#            - Mixture (qgcomp): re-run with IQR-scaled inputs (q=NULL) so
#              psi represents change per simultaneous 1-IQR increase across
#              all PFAS. Raw q=4 model also retained (Beta, P).
#          IQR scaling applied to continuous (Scenarios 1-2) and semi-
#          quantitative (Scenarios 5-6) scenarios only. Binary scenarios
#          (3-4, 7-8) are not scaled — no logical IQR for detect/non-detect.
#
# NOTE:    CLR transformation converts compositional count data to continuous
#          symmetric data, enabling linear models and qgcomp with gaussian
#          family. No offset needed — CLR accounts for library size
#          compositionally.
#          Each Model run time is approximately 25 minutes (especially the continuous with qgcomp)
#          Binary models are less than a minute
#
# set up -----------------------------------------------------------------------
rm(list = ls())
options(scipen = 0)

library(tidyverse)
library(dplyr)
library(readr)
library(here)
library(qgcomp)
library(tibble)
library(zCompositions)
library(compositions)
library(cowplot)
library(parallel)


n_cores <- detectCores() - 2
N_BOOT  <- 5000


# Load Taxonomy tables----------------------------------------------------------
tax_files <- list(
  species = "taxonomyDictionary_brack_jan_species_withLineage_bacteriaOnly.tsv",
  genus   = "taxonomyDictionary_brack_jan_genus_withLineage_bacteriaOnly.tsv",
  family  = "taxonomyDictionary_brack_jan_family_withLineage_bacteriaOnly.tsv",
  order   = "taxonomyDictionary_brack_jan_order_withLineage_bacteriaOnly.tsv",
  class   = "taxonomyDictionary_brack_jan_class_withLineage_bacteriaOnly.tsv",
  phylum  = "taxonomyDictionary_brack_jan_phylum_withLineage_bacteriaOnly.tsv"
)

taxTables_full <- lapply(tax_files, function(f) {
  read_tsv(here::here('input', f)) %>%
    mutate(taxonomy_id = paste0("X", taxonomy_id))
})

taxTables <- lapply(taxTables_full, function(t) {
  t %>% dplyr::select(name, taxonomy_id)
})


# Variable lists----------------------------------------------------------------
pfas_rename_continuous <- c(
  "PFBS"  = "PFBS_pgmL",
  "PFHxS" = "PFHxS_pgmL",
  "PFNA"  = "PFNA_pgmL",
  "PFOA"  = "PFOA_pgmL",
  "PFOS"  = "PFOS_pgmL"
)

pfas_vars_continuous <- c("PFBS", "PFHxS", "PFNA", "PFOA", "PFOS")
pfas_vars_binary     <- c("N.MeFOSAA", "PFBA", "PFDA", "PFDoA", "PFHpA")

pfas_rename_semiquant <- c(
  "N.MeFOSAA" = "N.MeFOSAA_pgmL",
  "PFBA"      = "PFBA_pgmL",
  "PFDA"      = "PFDA_pgmL",
  "PFDoA"     = "PFDoA_pgmL",
  "PFHpA"     = "PFHpA_pgmL"
)

pfas_vars_semiquant <- c("N.MeFOSAA", "PFBA", "PFDA", "PFDoA", "PFHpA")

# gestational_age coded as two dummies with Ontime as reference (0):
#   gest_Early = 1 if Early, 0 otherwise
#   gest_Late  = 1 if Late,  0 otherwise
# mode_of_delivery_bin: binary with Vaginal as reference (0 = Vaginal, 1 = C-Section)
# when rare categories are absent from a resample
covariates <- c("breastmilk_per_day", "SES_index_final", "baby_birthweight_kg",
                "gest_Early", "gest_Late", "mode_of_delivery_bin")

covariates_6m <- c(covariates, "age_of_solid_foods")

selectCovariates <- c("merge_id_dyad", "breastmilk_per_day", "SES_index_final",
                      "baby_birthweight_kg", "gest_Early", "gest_Late",
                      "mode_of_delivery_bin")

taxa_levels <- c("species", "genus", "family", "order", "class", "phylum")

level_prefix <- c(
  species = "s", genus = "g", family = "f",
  order = "o", class = "c", phylum = "p"
)

# Compute IQR of each log2-transformed continuous PFAS-------------------------
compute_pfas_iqr <- function(data) {
  sapply(pfas_vars_continuous, function(pfas) {
    IQR(data[[pfas]], na.rm = TRUE)
  })
}

# General IQR function for any PFAS variable list (used for semi-quant scenarios)
compute_pfas_iqr_custom <- function(data, pfas_vars) {
  sapply(pfas_vars, function(pfas) {
    IQR(data[[pfas]], na.rm = TRUE)
  })
}

# Helper functions--------------------------------------------------------------

# add_level_lineage: adds taxa_level, taxa_name, taxa_full
add_level_lineage <- function(df, level, taxTable_full) {
  prefix <- level_prefix[[level]]
  df %>%
    mutate(taxa_level = tools::toTitleCase(level)) %>%
    left_join(
      taxTable_full %>% dplyr::select(taxonomy_id, taxonomic_lineage),
      by = "taxonomy_id"
    ) %>%
    mutate(
      taxa_name = name,
      taxa_full = paste0(taxonomic_lineage, "_", prefix, "__", gsub(" ", "_", name))
    ) %>%
    dplyr::select(-taxonomic_lineage)
}


# prepare_taxa_data_clr: Re-close -> CZM imputation -> CLR transformation
prepare_taxa_data_clr <- function(data, pfas_vars, scenario_name,
                                  covs = covariates) {
  data$gest_Early           <- ifelse(data$gestational_age_cat  == "Early",     1, 0)
  data$gest_Late            <- ifelse(data$gestational_age_cat  == "Late",      1, 0)
  data$mode_of_delivery_bin <- ifelse(data$mode_of_delivery_cat == "C-Section", 1, 0)
  
  # Scale continuous covariates
  data$breastmilk_per_day  <- as.numeric(scale(data$breastmilk_per_day))
  data$SES_index_final     <- as.numeric(scale(data$SES_index_final))
  data$baby_birthweight_kg <- as.numeric(scale(data$baby_birthweight_kg))
  
  # Identify taxa columns
  taxa_cols <- grep("^X[0-9]", colnames(data), value = TRUE)
  cat("Taxa columns found:", length(taxa_cols), "\n")
  
  # Step 1: Input is Bracken relative abundance
  # After taxa filtering in Script 1, rows no longer sum to 1
  MB_counts <- as.matrix(data[, taxa_cols])
  cat("Input is Bracken relative abundance. Range:", round(range(MB_counts), 8), "\n")
  cat("Zeros in RA matrix:", sum(MB_counts == 0),
      "(", round(mean(MB_counts == 0) * 100, 1), "% )\n")
  
  # Step 1b: Re-close compositions after taxa subsetting
  # Taxa filtering removed columns, breaking compositional closure (rows != 1).
  # Re-normalization restores closure required by CZM and CLR.
  row_sums <- rowSums(MB_counts)
  cat("Row sums before re-closure — Range:", round(range(row_sums), 4), "\n")
  MB_counts <- MB_counts / row_sums
  cat("Row sums after re-closure — Range:", round(range(rowSums(MB_counts)), 6), "\n")
  
  # Step 2: CZM zero imputation
  if (sum(MB_counts == 0) > 0) {
    cat("Running CZM zero imputation...\n")
    MB_tss_imputed <- cmultRepl(
      MB_counts, label = 0, method = "CZM",
      z.warning = 1, z.delete = FALSE
    )
  } else {
    cat("No zeros found — skipping CZM.\n")
    MB_tss_imputed <- MB_counts           
  }
  
  # Step 3: CLR transformation
  MB_clr        <- t(apply(MB_tss_imputed, 1, clr))
  MB_clr_df     <- as.data.frame(MB_clr)
  colnames(MB_clr_df) <- paste0(taxa_cols, "_CLR")
  taxa_cols_clr <- colnames(MB_clr_df)
  cat("CLR done. Range:", round(range(MB_clr_df), 4), "\n")
  
  extra_cols <- intersect("n_detect", colnames(data))
  extra_covs <- setdiff(covs, covariates)
  
  data_clr <- data[, c(selectCovariates, extra_covs, extra_cols, pfas_vars)] %>%
    bind_cols(MB_clr_df)
  
  cat("Final N:", nrow(data_clr), "\n")
  
  attr(data_clr, "taxa_cols")     <- taxa_cols
  attr(data_clr, "taxa_cols_clr") <- taxa_cols_clr
  return(data_clr)
}


# save_clr_matrix: saves merge_id_dyad + CLR columns to out_files/
save_clr_matrix <- function(data_clr, filename) {
  taxa_cols_clr <- attr(data_clr, "taxa_cols_clr")
  data_clr %>%
    dplyr::select(merge_id_dyad, all_of(taxa_cols_clr)) %>%
    write.csv(here::here("out_files", filename), row.names = FALSE)
  cat("Saved:", filename, "\n")
}


# SAVE PER-SAMPLE CLR MATRICES — SPECIES LEVEL, 1m AND 6m-----------------------
# Re-run just species CLR prep for 1m and 6m
PFAS1m_micro1m_sp <- read.csv(here::here("out_files", "PFAS1m_micro1m_species.csv")) %>%
  dplyr::rename(!!!pfas_rename_continuous)
data_1m_sp <- prepare_taxa_data_clr(PFAS1m_micro1m_sp, pfas_vars_continuous,
                                    "species - 1m")
save_clr_matrix(data_1m_sp, "CLR_perSample_species_1m.csv")

PFAS1m_micro6m_sp <- read.csv(here::here("out_files", "PFAS1m_micro6m_species.csv")) %>%
  dplyr::rename(!!!pfas_rename_continuous)
data_6m_sp <- prepare_taxa_data_clr(PFAS1m_micro6m_sp, pfas_vars_continuous,
                                    "species - 6m")
save_clr_matrix(data_6m_sp, "CLR_perSample_species_6m.csv")


# run_clr_models: individual lm per PFAS per CLR taxon--------------------------
# pfas_iqr: named vector of IQRs for post-hoc scaling (continuous only)
run_clr_models <- function(data, pfas_vars, scenario_name, covs = covariates,
                           pfas_iqr = NULL) {
  
  taxa_cols_clr <- attr(data, "taxa_cols_clr")
  predictor_list <- character(0); outcome_list  <- character(0)
  betas_list     <- numeric(0);   se_list       <- numeric(0)
  ci_lo_list     <- numeric(0);   ci_hi_list    <- numeric(0)
  p_values_list  <- numeric(0);   n_list        <- numeric(0)
  
  for (i in taxa_cols_clr) {
    for (j in pfas_vars) {
      
      extra_covs   <- setdiff(covs, covariates)
      ithDataframe <- na.omit(data[, c(selectCovariates, extra_covs, i, j)])
      taxa_id      <- sub("_CLR$", "", i)
      
      formula <- paste(i, "~", j, "+",
                       paste(covs, collapse = " + ")) %>% as.formula()
      
      tryCatch({
        model    <- lm(formula, data = ithDataframe)
        ct       <- summary(model)$coefficients
        ci       <- confint(model)
        coef_row <- grep(paste0("^", j), rownames(ct), value = TRUE)[1]
        
        predictor_list <- c(predictor_list, j)
        outcome_list   <- c(outcome_list,   taxa_id)
        betas_list     <- c(betas_list,     ct[coef_row, "Estimate"])
        se_list        <- c(se_list,        ct[coef_row, "Std. Error"])
        ci_lo_list     <- c(ci_lo_list,     ci[coef_row, 1])
        ci_hi_list     <- c(ci_hi_list,     ci[coef_row, 2])
        p_values_list  <- c(p_values_list,  ct[coef_row, "Pr(>|t|)"])
        n_list         <- c(n_list,         nrow(ithDataframe))
        
      }, error = function(e) {
        warning(paste("Model failed for", taxa_id, "~", j, ":", e$message))
        predictor_list <<- c(predictor_list, j)
        outcome_list   <<- c(outcome_list,   taxa_id)
        betas_list     <<- c(betas_list,     NA)
        se_list        <<- c(se_list,        NA)
        ci_lo_list     <<- c(ci_lo_list,     NA)
        ci_hi_list     <<- c(ci_hi_list,     NA)
        p_values_list  <<- c(p_values_list,  NA)
        n_list         <<- c(n_list,         NA)
      })
    }
  }
  
  out <- data.frame(
    predictor = predictor_list, outcome = outcome_list,
    betas     = betas_list,     SE      = se_list,
    ci_lo     = ci_lo_list,     ci_hi   = ci_hi_list,
    p_val     = p_values_list,  n       = n_list,
    stringsAsFactors = FALSE
  )
  
  # Post-hoc IQR scaling — only applied when pfas_iqr is provided (continuous)
  if (!is.null(pfas_iqr)) {
    out <- out %>%
      mutate(
        Beta_IQR = ifelse(predictor %in% names(pfas_iqr),
                          betas * pfas_iqr[predictor], NA_real_),
        SE_IQR   = ifelse(predictor %in% names(pfas_iqr),
                          SE    * pfas_iqr[predictor], NA_real_)
      )
  }
  
  return(out)
}


# Apply FDR (BH correction per PFAS, merge taxonomy name)
apply_fdr <- function(model_outputs, pfas_vars, taxTable) {
  store_fdrs <- NULL
  for (pfas in pfas_vars) {
    ith_pfas      <- subset(model_outputs, predictor == pfas)
    ith_pfas$FDR  <- p.adjust(ith_pfas$p_val, method = "BH")
    store_fdrs    <- rbind(store_fdrs, ith_pfas)
  }
  store_fdrs %>%
    dplyr::rename(taxonomy_id = outcome) %>%
    merge(taxTable, all.x = TRUE, by = "taxonomy_id")
}

# run_ndetect_clr: lm with n_detect
run_ndetect_clr <- function(data, level, scenario_name, covs = covariates) {
  
  taxa_cols_clr <- attr(data, "taxa_cols_clr")
  outcome_list  <- character(0); betas_list    <- numeric(0)
  ci_lo_list    <- numeric(0);   ci_hi_list    <- numeric(0)
  p_values_list <- numeric(0);   n_list        <- numeric(0)
  
  for (i in taxa_cols_clr) {
    
    extra_covs   <- setdiff(covs, covariates)
    ithDataframe <- na.omit(data[, c(selectCovariates, extra_covs, "n_detect", i)])
    taxa_id      <- sub("_CLR$", "", i)
    
    formula <- paste(i, "~ n_detect +",
                     paste(covs, collapse = " + ")) %>% as.formula()
    
    print(paste(taxa_id, "n_detect"))
    
    tryCatch({
      model    <- lm(formula, data = ithDataframe)
      ct       <- summary(model)$coefficients
      ci       <- confint(model)
      
      outcome_list  <- c(outcome_list,  taxa_id)
      betas_list    <- c(betas_list,    ct["n_detect", "Estimate"])
      ci_lo_list    <- c(ci_lo_list,    ci["n_detect", 1])
      ci_hi_list    <- c(ci_hi_list,    ci["n_detect", 2])
      p_values_list <- c(p_values_list, ct["n_detect", "Pr(>|t|)"])
      n_list        <- c(n_list,        nrow(ithDataframe))
      
    }, error = function(e) {
      warning(paste("n_detect model failed for", taxa_id, ":", e$message))
      outcome_list  <<- c(outcome_list,  taxa_id)
      betas_list    <<- c(betas_list,    NA)
      ci_lo_list    <<- c(ci_lo_list,    NA)
      ci_hi_list    <<- c(ci_hi_list,    NA)
      p_values_list <<- c(p_values_list, NA)
      n_list        <<- c(n_list,        NA)
    })
  }
  
  data.frame(
    predictor = "N-detect", outcome = outcome_list,
    betas = betas_list, ci_lo = ci_lo_list, ci_hi = ci_hi_list,
    p_val = p_values_list, n = n_list,
    stringsAsFactors = FALSE
  ) %>%
    mutate(FDR = p.adjust(p_val, method = "BH")) %>%
    dplyr::rename(taxonomy_id = outcome) %>%
    merge(taxTables[[level]], all.x = TRUE, by = "taxonomy_id")
}


# run_qgcomp_clr_combined: runs both raw (q=4) and IQR-scaled (q=NULL) qgcomp--
# Returns Beta (raw), Beta_IQR (IQR-scaled), CIs, and single P from raw model
# NOTE: only called for continuous PFAS scenarios
run_qgcomp_clr_combined <- function(data, pfas_vars, pfas_iqr, level,
                                    scenario_name, covs = covariates) {
  
  taxa_cols_clr <- attr(data, "taxa_cols_clr")
  
  # IQR-scaled copy of data
  data_iqr <- data
  for (pfas in pfas_vars) {
    data_iqr[[pfas]] <- data_iqr[[pfas]] / pfas_iqr[pfas]
  }
  
  results_list <- mclapply(seq_along(taxa_cols_clr), function(idx) {
    
    i                <- taxa_cols_clr[idx]
    extra_covs       <- setdiff(covs, covariates)
    ithDataframe     <- na.omit(data[,     c(selectCovariates, extra_covs, pfas_vars, i)])
    ithDataframe_iqr <- na.omit(data_iqr[, c(selectCovariates, extra_covs, pfas_vars, i)])
    taxa_id          <- sub("_CLR$", "", i)
    
    out <- list(
      outcome   = taxa_id,  betas     = NA_real_,
      ci_lo     = NA_real_, ci_hi     = NA_real_,
      p_val     = NA_real_, Beta_IQR  = NA_real_,
      CI_lo_IQR = NA_real_, CI_hi_IQR = NA_real_,
      n         = nrow(ithDataframe)
    )
    
    if (nrow(ithDataframe) < 20) return(out)
    
    formula <- as.formula(
      paste0(i, " ~ ", paste(c(pfas_vars, covs), collapse = " + "))
    )
    
    tryCatch({
      # Raw model (q=4)
      qgmod_raw <- qgcomp.glm.boot(
        f = formula, data = ithDataframe, expnms = pfas_vars,
        q = 4, family = gaussian(), B = N_BOOT, seed = 2024, rr = FALSE
      )
      # IQR-scaled model (q=NULL)
      qgmod_iqr <- qgcomp.glm.boot(
        f = formula, data = ithDataframe_iqr, expnms = pfas_vars,
        q = NULL, family = gaussian(), B = N_BOOT, seed = 2024, rr = FALSE
      )
      out <- list(
        outcome   = taxa_id,
        betas     = as.numeric(qgmod_raw$psi[1]),
        ci_lo     = as.numeric(qgmod_raw$ci[1]),
        ci_hi     = as.numeric(qgmod_raw$ci[2]),
        p_val     = as.numeric(qgmod_raw$pval[2]),
        Beta_IQR  = as.numeric(qgmod_iqr$psi[1]),
        CI_lo_IQR = as.numeric(qgmod_iqr$ci[1]),
        CI_hi_IQR = as.numeric(qgmod_iqr$ci[2]),
        n         = nrow(ithDataframe)
      )
    }, error = function(e) {
      warning(paste("qgcomp combined failed for", taxa_id, ":", e$message))
    })
    
    return(out)
    
  }, mc.cores = n_cores)
  
  model_outputs <- bind_rows(lapply(results_list, as.data.frame)) %>%
    mutate(predictor = "Mixture",
           FDR       = p.adjust(p_val, method = "BH")) %>%
    dplyr::rename(taxonomy_id = outcome) %>%
    merge(taxTables[[level]], all.x = TRUE, by = "taxonomy_id")
  
  cat("  Converged:", sum(!is.na(model_outputs$betas)),
      "| Failed:", sum(is.na(model_outputs$betas)),
      "| Total:", nrow(model_outputs), "\n")
  
  return(model_outputs)
}


# Combined result collectors----------------------------------------------------
combined_cont_1m_1m      <- list()
combined_cont_1m_6m      <- list()
combined_bin_1m_1m       <- list()
combined_bin_1m_6m       <- list()
combined_semiquant_1m_1m <- list()
combined_semiquant_1m_6m <- list()
combined_sens_cont_1m_6m <- list()
combined_sens_bin_1m_6m  <- list()
combined_cont_1m_6m_bf   <- list()
combined_bin_1m_6m_bf    <- list() 

# SCENARIO 1: Continuous 1m PFAS + 1m microbiome--------------------------------
data_cont_1m_1m       <- list()
output_cont_1m_1m     <- list()
qgcomp_out_cont_1m_1m <- list()

# Compute IQR from species-level data (same PFAS values across all levels)
PFAS1m_micro1m_iqr <- read.csv(here::here("out_files", "PFAS1m_micro1m_species.csv")) %>%
  dplyr::rename(!!!pfas_rename_continuous)
pfas_iqr_1m <- compute_pfas_iqr(PFAS1m_micro1m_iqr)
cat("PFAS IQRs (1m data):\n"); print(round(pfas_iqr_1m, 4))

for (level in taxa_levels) {
  PFAS1m_micro1m <- read.csv(here::here("out_files", paste0("PFAS1m_micro1m_", level, ".csv"))) %>%
    dplyr::rename(!!!pfas_rename_continuous)
  data_cont_1m_1m[[level]] <- prepare_taxa_data_clr(
    PFAS1m_micro1m, pfas_vars_continuous,
    paste(level, "- Continuous 1m PFAS + 1m Microbiome")
  )
  clr_res <- run_clr_models(
    data_cont_1m_1m[[level]], pfas_vars_continuous,
    paste(level, "- Continuous 1m PFAS + 1m Microbiome"),
    pfas_iqr = pfas_iqr_1m
  )
  output_cont_1m_1m[[level]] <- apply_fdr(clr_res, pfas_vars_continuous, taxTables[[level]])
  qgcomp_out_cont_1m_1m[[level]] <- run_qgcomp_clr_combined(
    data_cont_1m_1m[[level]], pfas_vars_continuous, pfas_iqr_1m, level,
    paste(level, "- Continuous 1m PFAS + 1m Microbiome")
  )
  combined_cont_1m_1m[[level]] <- bind_rows(
    add_level_lineage(output_cont_1m_1m[[level]],        level, taxTables_full[[level]]),
    add_level_lineage(qgcomp_out_cont_1m_1m[[level]],    level, taxTables_full[[level]])
  )
}
write.csv(bind_rows(combined_cont_1m_1m),
          here::here("out_files", "COMBINED_CLR_continuous_1m_1m.csv"), row.names = FALSE)


# SCENARIO 2: Continuous 1m PFAS + 6m microbiome--------------------------------
data_cont_1m_6m       <- list()
output_cont_1m_6m     <- list()
qgcomp_out_cont_1m_6m <- list()

# Compute IQR from species-level data
PFAS1m_micro6m_iqr <- read.csv(here::here("out_files", "PFAS1m_micro6m_species.csv")) %>%
  dplyr::rename(!!!pfas_rename_continuous)
pfas_iqr_6m <- compute_pfas_iqr(PFAS1m_micro6m_iqr)
cat("PFAS IQRs (6m data):\n"); print(round(pfas_iqr_6m, 4))

for (level in taxa_levels) {
  PFAS1m_micro6m <- read.csv(here::here("out_files", paste0("PFAS1m_micro6m_", level, ".csv"))) %>%
    dplyr::rename(!!!pfas_rename_continuous)
  data_cont_1m_6m[[level]] <- prepare_taxa_data_clr(
    PFAS1m_micro6m, pfas_vars_continuous,
    paste(level, "- Continuous 1m PFAS + 6m Microbiome")
  )
  clr_res <- run_clr_models(
    data_cont_1m_6m[[level]], pfas_vars_continuous,
    paste(level, "- Continuous 1m PFAS + 6m Microbiome"),
    pfas_iqr = pfas_iqr_6m
  )
  output_cont_1m_6m[[level]] <- apply_fdr(clr_res, pfas_vars_continuous, taxTables[[level]])
  qgcomp_out_cont_1m_6m[[level]] <- run_qgcomp_clr_combined(
    data_cont_1m_6m[[level]], pfas_vars_continuous, pfas_iqr_6m, level,
    paste(level, "- Continuous 1m PFAS + 6m Microbiome")
  )
  combined_cont_1m_6m[[level]] <- bind_rows(
    add_level_lineage(output_cont_1m_6m[[level]],        level, taxTables_full[[level]]),
    add_level_lineage(qgcomp_out_cont_1m_6m[[level]],    level, taxTables_full[[level]])
  )
}
write.csv(bind_rows(combined_cont_1m_6m),
          here::here("out_files", "COMBINED_CLR_continuous_1m_6m.csv"), row.names = FALSE)



# SCENARIO 3: Binary 1m PFAS + 1m microbiome------------------------------------
data_bin_1m_1m    <- list()
output_bin_1m_1m  <- list()
ndetect_1m_1m     <- list()

for (level in taxa_levels) {
  PFAS1mDetect_micro1m <- read.csv(here::here("out_files", paste0("PFAS1mDetect_micro1m_", level, ".csv")))
  detect_cols <- grep("_Detect$", colnames(PFAS1mDetect_micro1m), value = TRUE)
  PFAS1mDetect_micro1m$n_detect <- rowSums(PFAS1mDetect_micro1m[detect_cols] == "detect", na.rm = TRUE)
  PFAS1mDetect_micro1m <- PFAS1mDetect_micro1m %>%
    dplyr::rename(!!!setNames(paste0(pfas_vars_binary, "_Detect"), pfas_vars_binary)) %>%
    mutate(across(all_of(pfas_vars_binary), ~ factor(., levels = c("non-detect", "detect"))))
  data_bin_1m_1m[[level]] <- prepare_taxa_data_clr(
    PFAS1mDetect_micro1m, pfas_vars_binary,
    paste(level, "- Binary 1m PFAS + 1m Microbiome")
  )
  clr_res <- run_clr_models(
    data_bin_1m_1m[[level]], pfas_vars_binary,
    paste(level, "- Binary 1m PFAS + 1m Microbiome")
  )
  output_bin_1m_1m[[level]] <- apply_fdr(clr_res, pfas_vars_binary, taxTables[[level]])
  ndetect_1m_1m[[level]]    <- run_ndetect_clr(
    data_bin_1m_1m[[level]], level,
    paste(level, "- Binary 1m PFAS + 1m Microbiome")
  )
  combined_bin_1m_1m[[level]] <- bind_rows(
    add_level_lineage(output_bin_1m_1m[[level]], level, taxTables_full[[level]]),
    add_level_lineage(ndetect_1m_1m[[level]],    level, taxTables_full[[level]])
  )
}
write.csv(bind_rows(combined_bin_1m_1m),
          here::here("out_files", "COMBINED_CLR_binary_1m_1m.csv"), row.names = FALSE)

# SCENARIO 4: Binary 1m PFAS + 6m microbiome------------------------------------
data_bin_1m_6m    <- list()
output_bin_1m_6m  <- list()
ndetect_1m_6m     <- list()

for (level in taxa_levels) {
  PFAS1mDetect_micro6m <- read.csv(here::here("out_files", paste0("PFAS1mDetect_micro6m_", level, ".csv")))
  detect_cols <- grep("_Detect$", colnames(PFAS1mDetect_micro6m), value = TRUE)
  PFAS1mDetect_micro6m$n_detect <- rowSums(PFAS1mDetect_micro6m[detect_cols] == "detect", na.rm = TRUE)
  PFAS1mDetect_micro6m <- PFAS1mDetect_micro6m %>%
    dplyr::rename(!!!setNames(paste0(pfas_vars_binary, "_Detect"), pfas_vars_binary)) %>%
    mutate(across(all_of(pfas_vars_binary), ~ factor(., levels = c("non-detect", "detect"))))
  data_bin_1m_6m[[level]] <- prepare_taxa_data_clr(
    PFAS1mDetect_micro6m, pfas_vars_binary,
    paste(level, "- Binary 1m PFAS + 6m Microbiome")
  )
  clr_res <- run_clr_models(
    data_bin_1m_6m[[level]], pfas_vars_binary,
    paste(level, "- Binary 1m PFAS + 6m Microbiome")
  )
  output_bin_1m_6m[[level]] <- apply_fdr(clr_res, pfas_vars_binary, taxTables[[level]])
  ndetect_1m_6m[[level]]    <- run_ndetect_clr(
    data_bin_1m_6m[[level]], level,
    paste(level, "- Binary 1m PFAS + 6m Microbiome")
  )
  combined_bin_1m_6m[[level]] <- bind_rows(
    add_level_lineage(output_bin_1m_6m[[level]], level, taxTables_full[[level]]),
    add_level_lineage(ndetect_1m_6m[[level]],    level, taxTables_full[[level]])
  )
}
write.csv(bind_rows(combined_bin_1m_6m),
          here::here("out_files", "COMBINED_CLR_binary_1m_6m.csv"), row.names = FALSE)

# SCENARIO 5 (SENSITIVITY): Semi-quantitative 1m PFAS + 1m microbiome----------
data_semiquant_1m_1m       <- list()
output_semiquant_1m_1m     <- list()
qgcomp_out_semiquant_1m_1m <- list()

# Compute IQR from semi-quantitative PFAS values
PFAS1mDetect_micro1m_iqr <- read.csv(here::here("out_files", "PFAS1mDetect_micro1m_species.csv")) %>%
  dplyr::select(-any_of(pfas_vars_binary), -any_of("n_detect")) %>%
  dplyr::rename(!!!pfas_rename_semiquant)
pfas_iqr_semiquant_1m <- compute_pfas_iqr_custom(PFAS1mDetect_micro1m_iqr, pfas_vars_semiquant)
cat("PFAS IQRs (semi-quant 1m data):\n"); print(round(pfas_iqr_semiquant_1m, 4))

for (level in taxa_levels) {
  PFAS1mDetect_micro1m <- read.csv(here::here("out_files", paste0("PFAS1mDetect_micro1m_", level, ".csv"))) %>%
    dplyr::select(-any_of(pfas_vars_binary), -any_of("n_detect")) %>%
    dplyr::rename(!!!pfas_rename_semiquant)
  data_semiquant_1m_1m[[level]] <- prepare_taxa_data_clr(
    PFAS1mDetect_micro1m, pfas_vars_semiquant,
    paste(level, "- Semi-quantitative 1m PFAS + 1m Microbiome")
  )
  clr_res <- run_clr_models(
    data_semiquant_1m_1m[[level]], pfas_vars_semiquant,
    paste(level, "- Semi-quantitative 1m PFAS + 1m Microbiome"),
    pfas_iqr = pfas_iqr_semiquant_1m
  )
  output_semiquant_1m_1m[[level]] <- apply_fdr(clr_res, pfas_vars_semiquant, taxTables[[level]])
  qgcomp_out_semiquant_1m_1m[[level]] <- run_qgcomp_clr_combined(
    data_semiquant_1m_1m[[level]], pfas_vars_semiquant, pfas_iqr_semiquant_1m, level,
    paste(level, "- Semi-quantitative 1m PFAS + 1m Microbiome")
  )
  combined_semiquant_1m_1m[[level]] <- bind_rows(
    add_level_lineage(output_semiquant_1m_1m[[level]],     level, taxTables_full[[level]]),
    add_level_lineage(qgcomp_out_semiquant_1m_1m[[level]], level, taxTables_full[[level]])
  )
}
write.csv(bind_rows(combined_semiquant_1m_1m),
          here::here("out_files", "COMBINED_CLR_semiquant_1m_1m.csv"), row.names = FALSE)


# SCENARIO 6 (SENSITIVITY): Semi-quantitative 1m PFAS + 6m microbiome-----------
data_semiquant_1m_6m       <- list()
output_semiquant_1m_6m     <- list()
qgcomp_out_semiquant_1m_6m <- list()

# Compute IQR from semi-quantitative PFAS values (species level)
PFAS1mDetect_micro6m_iqr <- read.csv(here::here("out_files", "PFAS1mDetect_micro6m_species.csv")) %>%
  dplyr::select(-any_of(pfas_vars_binary), -any_of("n_detect")) %>%
  dplyr::rename(!!!pfas_rename_semiquant)
pfas_iqr_semiquant_6m <- compute_pfas_iqr_custom(PFAS1mDetect_micro6m_iqr, pfas_vars_semiquant)
cat("PFAS IQRs (semi-quant 6m data):\n"); print(round(pfas_iqr_semiquant_6m, 4))

for (level in taxa_levels) {
  PFAS1mDetect_micro6m <- read.csv(here::here("out_files", paste0("PFAS1mDetect_micro6m_", level, ".csv"))) %>%
    dplyr::select(-any_of(pfas_vars_binary), -any_of("n_detect")) %>%
    dplyr::rename(!!!pfas_rename_semiquant)
  data_semiquant_1m_6m[[level]] <- prepare_taxa_data_clr(
    PFAS1mDetect_micro6m, pfas_vars_semiquant,
    paste(level, "- Semi-quantitative 1m PFAS + 6m Microbiome")
  )
  clr_res <- run_clr_models(
    data_semiquant_1m_6m[[level]], pfas_vars_semiquant,
    paste(level, "- Semi-quantitative 1m PFAS + 6m Microbiome"),
    pfas_iqr = pfas_iqr_semiquant_6m
  )
  output_semiquant_1m_6m[[level]] <- apply_fdr(clr_res, pfas_vars_semiquant, taxTables[[level]])
  qgcomp_out_semiquant_1m_6m[[level]] <- run_qgcomp_clr_combined(
    data_semiquant_1m_6m[[level]], pfas_vars_semiquant, pfas_iqr_semiquant_6m, level,
    paste(level, "- Semi-quantitative 1m PFAS + 6m Microbiome")
  )
  combined_semiquant_1m_6m[[level]] <- bind_rows(
    add_level_lineage(output_semiquant_1m_6m[[level]],     level, taxTables_full[[level]]),
    add_level_lineage(qgcomp_out_semiquant_1m_6m[[level]], level, taxTables_full[[level]])
  )
}
write.csv(bind_rows(combined_semiquant_1m_6m),
          here::here("out_files", "COMBINED_CLR_semiquant_1m_6m.csv"), row.names = FALSE)


# SCENARIO 7 (SENSITIVITY): Continuous 1m PFAS + 6m microbiome + age_of_solid_foods----
data_cont_1m_6m_sens       <- list()
output_cont_1m_6m_sens     <- list()
qgcomp_out_cont_1m_6m_sens <- list()

for (level in taxa_levels) {
  PFAS1m_micro6m <- read.csv(here::here("out_files", paste0("PFAS1m_micro6m_", level, ".csv"))) %>%
    dplyr::rename(!!!pfas_rename_continuous)
  data_cont_1m_6m_sens[[level]] <- prepare_taxa_data_clr(
    PFAS1m_micro6m, pfas_vars_continuous,
    paste(level, "- Continuous 1m PFAS + 6m Sensitivity"),
    covs = covariates_6m
  )
  clr_res <- run_clr_models(
    data_cont_1m_6m_sens[[level]], pfas_vars_continuous,
    paste(level, "- Continuous 1m PFAS + 6m Sensitivity"),
    covs = covariates_6m,
    pfas_iqr = pfas_iqr_6m
  )
  output_cont_1m_6m_sens[[level]] <- apply_fdr(clr_res, pfas_vars_continuous, taxTables[[level]])
  qgcomp_out_cont_1m_6m_sens[[level]] <- run_qgcomp_clr_combined(
    data_cont_1m_6m_sens[[level]], pfas_vars_continuous, pfas_iqr_6m, level,
    paste(level, "- Continuous 1m PFAS + 6m Sensitivity"),
    covs = covariates_6m
  )
  combined_sens_cont_1m_6m[[level]] <- bind_rows(
    add_level_lineage(output_cont_1m_6m_sens[[level]],          level, taxTables_full[[level]]),
    add_level_lineage(qgcomp_out_cont_1m_6m_sens[[level]],      level, taxTables_full[[level]])
  )
}
write.csv(bind_rows(combined_sens_cont_1m_6m),
          here::here("out_files", "COMBINED_CLR_sensitivity_continuous_1m_6m.csv"), row.names = FALSE)

# SCENARIO 8 (SENSITIVITY): Binary 1m PFAS + 6m microbiome + age_of_solid_foods----
data_bin_1m_6m_sens    <- list()
output_bin_1m_6m_sens  <- list()
ndetect_1m_6m_sens     <- list()

for (level in taxa_levels) {
  PFAS1mDetect_micro6m <- read.csv(here::here("out_files", paste0("PFAS1mDetect_micro6m_", level, ".csv")))
  detect_cols <- grep("_Detect$", colnames(PFAS1mDetect_micro6m), value = TRUE)
  PFAS1mDetect_micro6m$n_detect <- rowSums(PFAS1mDetect_micro6m[detect_cols] == "detect", na.rm = TRUE)
  PFAS1mDetect_micro6m <- PFAS1mDetect_micro6m %>%
    dplyr::rename(!!!setNames(paste0(pfas_vars_binary, "_Detect"), pfas_vars_binary)) %>%
    mutate(across(all_of(pfas_vars_binary), ~ factor(., levels = c("non-detect", "detect"))))
  data_bin_1m_6m_sens[[level]] <- prepare_taxa_data_clr(
    PFAS1mDetect_micro6m, pfas_vars_binary,
    paste(level, "- Binary 1m PFAS + 6m Sensitivity"),
    covs = covariates_6m
  )
  clr_res <- run_clr_models(
    data_bin_1m_6m_sens[[level]], pfas_vars_binary,
    paste(level, "- Binary 1m PFAS + 6m Sensitivity"),
    covs = covariates_6m
  )
  output_bin_1m_6m_sens[[level]] <- apply_fdr(clr_res, pfas_vars_binary, taxTables[[level]])
  ndetect_1m_6m_sens[[level]]    <- run_ndetect_clr(
    data_bin_1m_6m_sens[[level]], level,
    paste(level, "- Binary 1m PFAS + 6m Sensitivity"),
    covs = covariates_6m
  )
  combined_sens_bin_1m_6m[[level]] <- bind_rows(
    add_level_lineage(output_bin_1m_6m_sens[[level]], level, taxTables_full[[level]]),
    add_level_lineage(ndetect_1m_6m_sens[[level]],    level, taxTables_full[[level]])
  )
}
write.csv(bind_rows(combined_sens_bin_1m_6m),
          here::here("out_files", "COMBINED_CLR_sensitivity_binary_1m_6m.csv"), row.names = FALSE)


# SCENARIO 9 (SENSITIVITY): Continuous 1m PFAS + 6m microbiome — BF ≥1/day at 6m ----
# Filter to participants with breastmilk_per_day >= 1 at 6 months (from meta_trim)
meta_trim <- read.csv(here::here("out_files", "meta_trim.csv"))
ids_bf_6m <- meta_trim %>%
  filter(timepoint == 6, breastmilk_per_day >= 0.5) %>%
  pull(merge_id_dyad)
cat("Participants with BF >= 0.5/day at 6m:", length(ids_bf_6m), "\n")

# Recompute pfas_iqr_6m in case session was restarted
if (!exists("pfas_iqr_6m")) {
  PFAS1m_micro6m_iqr <- read.csv(here::here("out_files", "PFAS1m_micro6m_species.csv")) %>%
    dplyr::rename(!!!pfas_rename_continuous)
  pfas_iqr_6m <- compute_pfas_iqr(PFAS1m_micro6m_iqr)
  cat("pfas_iqr_6m recomputed:\n"); print(round(pfas_iqr_6m, 4))
}

data_cont_1m_6m_bf       <- list()
output_cont_1m_6m_bf     <- list()
qgcomp_out_cont_1m_6m_bf <- list()
combined_cont_1m_6m_bf   <- list()

for (level in taxa_levels) {
  PFAS1m_micro6m_bf <- read.csv(here::here("out_files", paste0("PFAS1m_micro6m_", level, ".csv"))) %>%
    dplyr::rename(!!!pfas_rename_continuous) %>%
    filter(merge_id_dyad %in% ids_bf_6m)
  cat("Scenario 9 |", level, "| N after BF filter:", nrow(PFAS1m_micro6m_bf), "\n")
  
  data_cont_1m_6m_bf[[level]] <- prepare_taxa_data_clr(
    PFAS1m_micro6m_bf, pfas_vars_continuous,
    paste(level, "- Continuous 1m PFAS + 6m BF Sensitivity")
  )
  clr_res <- run_clr_models(
    data_cont_1m_6m_bf[[level]], pfas_vars_continuous,
    paste(level, "- Continuous 1m PFAS + 6m BF Sensitivity"),
    pfas_iqr = pfas_iqr_6m
  )
  output_cont_1m_6m_bf[[level]]     <- apply_fdr(clr_res, pfas_vars_continuous, taxTables[[level]])
  qgcomp_out_cont_1m_6m_bf[[level]] <- run_qgcomp_clr_combined(
    data_cont_1m_6m_bf[[level]], pfas_vars_continuous, pfas_iqr_6m, level,
    paste(level, "- Continuous 1m PFAS + 6m BF Sensitivity")
  )
  combined_cont_1m_6m_bf[[level]] <- bind_rows(
    add_level_lineage(output_cont_1m_6m_bf[[level]],     level, taxTables_full[[level]]),
    add_level_lineage(qgcomp_out_cont_1m_6m_bf[[level]], level, taxTables_full[[level]])
  )
}
write.csv(bind_rows(combined_cont_1m_6m_bf),
          here::here("out_files", "COMBINED_CLR_sensitivity_cont_1m_6m_bf.csv"), row.names = FALSE)


# SCENARIO 10 (SENSITIVITY): Binary 1m PFAS + 6m microbiome — BF ≥1/day at 6m ----
data_bin_1m_6m_bf    <- list()
output_bin_1m_6m_bf  <- list()
ndetect_1m_6m_bf     <- list()
combined_bin_1m_6m_bf <- list()

for (level in taxa_levels) {
  PFAS1mDetect_micro6m_bf <- read.csv(here::here("out_files", paste0("PFAS1mDetect_micro6m_", level, ".csv"))) %>%
    filter(merge_id_dyad %in% ids_bf_6m)
  cat("Scenario 10 |", level, "| N after BF filter:", nrow(PFAS1mDetect_micro6m_bf), "\n")
  
  detect_cols <- grep("_Detect$", colnames(PFAS1mDetect_micro6m_bf), value = TRUE)
  PFAS1mDetect_micro6m_bf$n_detect <- rowSums(PFAS1mDetect_micro6m_bf[detect_cols] == "detect", na.rm = TRUE)
  PFAS1mDetect_micro6m_bf <- PFAS1mDetect_micro6m_bf %>%
    dplyr::rename(!!!setNames(paste0(pfas_vars_binary, "_Detect"), pfas_vars_binary)) %>%
    mutate(across(all_of(pfas_vars_binary), ~ factor(., levels = c("non-detect", "detect"))))
  
  data_bin_1m_6m_bf[[level]] <- prepare_taxa_data_clr(
    PFAS1mDetect_micro6m_bf, pfas_vars_binary,
    paste(level, "- Binary 1m PFAS + 6m BF Sensitivity")
  )
  clr_res <- run_clr_models(
    data_bin_1m_6m_bf[[level]], pfas_vars_binary,
    paste(level, "- Binary 1m PFAS + 6m BF Sensitivity")
  )
  output_bin_1m_6m_bf[[level]]  <- apply_fdr(clr_res, pfas_vars_binary, taxTables[[level]])
  ndetect_1m_6m_bf[[level]]     <- run_ndetect_clr(
    data_bin_1m_6m_bf[[level]], level,
    paste(level, "- Binary 1m PFAS + 6m BF Sensitivity")
  )
  combined_bin_1m_6m_bf[[level]] <- bind_rows(
    add_level_lineage(output_bin_1m_6m_bf[[level]], level, taxTables_full[[level]]),
    add_level_lineage(ndetect_1m_6m_bf[[level]],    level, taxTables_full[[level]])
  )
}
write.csv(bind_rows(combined_bin_1m_6m_bf),
          here::here("out_files", "COMBINED_CLR_sensitivity_bin_1m_6m_bf.csv"), row.names = FALSE)


# Taxa Heatmap Visualization----------------------------------------------------
# set up
rm(list = ls())

# PFAS display order
pfas_order_continuous <- c("PFBS", "PFHxS", "PFNA", "PFOA", "PFOS")
pfas_order_binary     <- c("N-MeFOSAA", "PFBA", "PFDA", "PFDoA", "PFHpA")
binary_display_rename <- c("N.MeFOSAA" = "N-MeFOSAA")


# SECTION 1 — Load all combined results

combined_files <- list(
  cont_1m_1m        = here::here("out_files", "COMBINED_CLR_continuous_1m_1m.csv"),
  cont_1m_6m        = here::here("out_files", "COMBINED_CLR_continuous_1m_6m.csv"),
  bin_1m_1m         = here::here("out_files", "COMBINED_CLR_binary_1m_1m.csv"),
  bin_1m_6m         = here::here("out_files", "COMBINED_CLR_binary_1m_6m.csv"),
  semiquant_1m_1m   = here::here("out_files", "COMBINED_CLR_semiquant_1m_1m.csv"),
  semiquant_1m_6m   = here::here("out_files", "COMBINED_CLR_semiquant_1m_6m.csv"),
  sens_cont_1m_6m   = here::here("out_files", "COMBINED_CLR_sensitivity_continuous_1m_6m.csv"),
  sens_bin_1m_6m    = here::here("out_files", "COMBINED_CLR_sensitivity_binary_1m_6m.csv"),
  sens_cont_1m_6m_bf = here::here("out_files", "COMBINED_CLR_sensitivity_cont_1m_6m_bf.csv"),
  sens_bin_1m_6m_bf  = here::here("out_files", "COMBINED_CLR_sensitivity_bin_1m_6m_bf.csv")
)

scenario_labels <- c(
  cont_1m_1m      = "Continuous 1m PFAS + 1m Microbiome",
  cont_1m_6m      = "Continuous 1m PFAS + 6m Microbiome",
  bin_1m_1m       = "Binary 1m PFAS + 1m Microbiome",
  bin_1m_6m       = "Binary 1m PFAS + 6m Microbiome",
  semiquant_1m_1m = "Semi-quantitative 1m PFAS + 1m Microbiome",
  semiquant_1m_6m = "Semi-quantitative 1m PFAS + 6m Microbiome",
  sens_cont_1m_6m = "Sensitivity: Continuous 1m PFAS + 6m Microbiome",
  sens_bin_1m_6m  = "Sensitivity: Binary 1m PFAS + 6m Microbiome",
  sens_cont_1m_6m_bf = "Sensitivity BF: Continuous 1m PFAS + 6m Microbiome",
  sens_bin_1m_6m_bf  = "Sensitivity BF: Binary 1m PFAS + 6m Microbiome"
)

# Load main results
all_data <- lapply(names(combined_files), function(scen) {
  read.csv(combined_files[[scen]]) %>%
    mutate(scenario = scenario_labels[[scen]])
}) %>% bind_rows()

# FUNCTION: make_taxa_heatmap
make_taxa_heatmap <- function(df,
                              output_file,
                              pfas_order,
                              top_n          = NULL,
                              display_rename = NULL) {
  
  plot_title <- unique(df$scenario)
  
  df <- df %>%
    filter(!is.na(p_val), p_val < 0.05) %>%
    mutate(
      betas     = as.numeric(betas),
      Beta_IQR  = as.numeric(Beta_IQR),
      FDR       = as.numeric(FDR),
      p_val     = as.numeric(p_val),
      taxa_name = gsub("\\[|\\]", "", taxa_name),
      # Use IQR-scaled beta where available, fall back to raw beta
      plot_beta = dplyr::coalesce(Beta_IQR, betas)
    )
  
  cat("p < 0.05 associations:", nrow(df), "\n")
  cat("of which FDR < 0.05:  ", sum(df$FDR < 0.05, na.rm = TRUE), "\n")
  
  if (nrow(df) == 0) {
    cat("  -> No significant associations. Skipping.\n")
    return(invisible(NULL))
  }
  
  if (!is.null(display_rename)) {
    df <- df %>% mutate(predictor = recode(predictor, !!!display_rename))
  }
  
  # Rank taxa by best p_val, cap if top_n specified, then sort alphabetically
  all_taxa_ranked <- df %>%
    group_by(taxa_name) %>%
    slice_min(p_val, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(p_val) %>%
    pull(taxa_name)
  
  if (!is.null(top_n)) {
    all_taxa_ranked <- head(all_taxa_ranked, top_n)
    cat("Capped to top", top_n, "taxa by best p-value.\n")
  }
  
  all_taxa <- sort(all_taxa_ranked, decreasing = TRUE)
  
  df_top <- df %>%
    filter(taxa_name %in% all_taxa) %>%
    mutate(
      taxa_name = factor(taxa_name, levels = all_taxa),
      predictor = factor(predictor, levels = pfas_order)
    ) %>%
    tidyr::complete(predictor, taxa_name) %>%
    mutate(
      plot_beta = dplyr::coalesce(Beta_IQR, betas),
      star = case_when(
        !is.na(FDR) & FDR < 0.05 ~ "*",
        TRUE                      ~ ""
      )
    )
  
  # Figure dimensions
  n_taxa  <- length(all_taxa)
  n_pfas  <- length(pfas_order)
  panel_h <- max(3.0, n_taxa * 0.30)
  panel_w <- max(2.0, n_pfas * 0.55)
  leg_w   <- 1.40
  fig_w   <- panel_w + leg_w
  fig_h   <- panel_h + 1.0
  
  p_heat <- ggplot(df_top, aes(x = predictor, y = taxa_name, fill = plot_beta)) +
    geom_tile(color = "grey50", linewidth = 0.15, linetype = "dotted", na.rm = FALSE) +
    geom_text(aes(label = star), size = 4, color = "black", vjust = 0.75, na.rm = TRUE) +
    scale_fill_gradient2(
      low      = "#d7191c",
      mid      = "white",
      high     = "#2c7bb6",
      midpoint = 0,
      na.value = "white",
      name     = "Effect Size\n(IQR-scaled\nCLR units)"
    ) +
    scale_x_discrete(position = "bottom") +
    coord_fixed(ratio = 0.5) +
    labs(title = NULL, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid       = element_blank(),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 9,
                                      face = "bold", color = "black"),
      axis.text.y      = element_text(size = 8, face = "italic", color = "black"),
      legend.position  = "none",
      plot.margin      = margin(5, 5, 5, 5)
    )
  
  p_for_legend <- p_heat +
    theme(legend.position = "right",
          legend.title    = element_text(size = 9),
          legend.text     = element_text(size = 8)) +
    guides(fill = guide_colorbar(
      barheight      = unit(1.2, "in"),
      barwidth       = unit(0.15, "in"),
      ticks          = TRUE,
      title.position = "top"
    ))
  
  leg_gt   <- cowplot::get_legend(p_for_legend)
  leg_plot <- cowplot::ggdraw(leg_gt)
  
  note_plot <- ggplot() +
    theme_void() +
    annotate("text", x = 0.05, y = 0.70, hjust = 0, vjust = 1,
             label = "Colored: p < 0.05", size = 2.8, color = "black") +
    annotate("text", x = 0.05, y = 0.30, hjust = 0, vjust = 1,
             label = "* FDR < 0.05",     size = 2.8, color = "black") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE)
  
  bar_h_in   <- 1.2
  bar_h_frac <- bar_h_in / fig_h
  annot_h    <- 0.35 / fig_h
  bar_y_mid  <- 0.50
  bar_ymin   <- bar_y_mid - bar_h_frac / 2
  bar_ymax   <- bar_y_mid + bar_h_frac / 2
  annot_ymax <- bar_ymin + 0.04
  annot_ymin <- annot_ymax - annot_h
  
  right_col <- cowplot::ggdraw() +
    cowplot::draw_plot(leg_plot,
                       x = -0.20, y = bar_ymin,
                       width = 1, height = bar_h_frac + 0.12) +
    cowplot::draw_plot(note_plot,
                       x = 0, y = annot_ymin,
                       width = 1, height = annot_h)
  
  p_final <- cowplot::plot_grid(p_heat, right_col, ncol = 2,
                                rel_widths = c(panel_w, leg_w))
  
  print(p_final)
  ggsave(here::here("out_figures", output_file),
         plot = p_final, width = fig_w, height = fig_h,
         device = cairo_pdf, units = "in")
  cat("Saved:", output_file, "(", round(fig_w, 1), "x", round(fig_h, 1), "in )\n")
  return(invisible(p_final))
}


# Species-only taxa association heatmaps (top 40, p < 0.05)

get_species <- function(scenario_name) {
  all_data %>% filter(scenario == scenario_labels[[scenario_name]],
                      taxa_level == "Species")
}

# Scenario 1: Continuous 1m PFAS + 1m Microbiome
make_taxa_heatmap(
  df          = get_species("cont_1m_1m"),
  output_file = "heatmap_CLR_cont_1m_1m_species.pdf",
  pfas_order  = c(pfas_order_continuous, "Mixture"),
  top_n       = 40
)

# Scenario 2: Continuous 1m PFAS + 6m Microbiome
make_taxa_heatmap(
  df          = get_species("cont_1m_6m"),
  output_file = "heatmap_CLR_cont_1m_6m_species.pdf",
  pfas_order  = c(pfas_order_continuous, "Mixture"),
  top_n       = 40
)

# Scenario 3: Binary 1m PFAS + 1m Microbiome
make_taxa_heatmap(
  df             = get_species("bin_1m_1m"),
  output_file    = "heatmap_CLR_bin_1m_1m_species.pdf",
  pfas_order     = c(pfas_order_binary, "N-detect"),
  top_n          = 40,
  display_rename = binary_display_rename
)

# Scenario 4: Binary 1m PFAS + 6m Microbiome
make_taxa_heatmap(
  df             = get_species("bin_1m_6m"),
  output_file    = "heatmap_CLR_bin_1m_6m_species.pdf",
  pfas_order     = c(pfas_order_binary, "N-detect"),
  top_n          = 40,
  display_rename = binary_display_rename
)

# Scenario 5: Semi-quantitative 1m PFAS + 1m Microbiome
make_taxa_heatmap(
  df             = get_species("semiquant_1m_1m"),
  output_file    = "heatmap_CLR_semiquant_1m_1m_species.pdf",
  pfas_order     = c(pfas_order_binary, "Mixture"),
  top_n          = 40,
  display_rename = binary_display_rename
)

# Scenario 6: Semi-quantitative 1m PFAS + 6m Microbiome
make_taxa_heatmap(
  df             = get_species("semiquant_1m_6m"),
  output_file    = "heatmap_CLR_semiquant_1m_6m_species.pdf",
  pfas_order     = c(pfas_order_binary, "Mixture"),
  top_n          = 40,
  display_rename = binary_display_rename
)

# Scenario 7: Sensitivity Continuous 1m PFAS + 6m Microbiome
make_taxa_heatmap(
  df          = get_species("sens_cont_1m_6m"),
  output_file = "heatmap_CLR_cont_1m_6m_sensitivity_species.pdf",
  pfas_order  = c(pfas_order_continuous, "Mixture"),
  top_n       = 40
)

# Scenario 8: Sensitivity Binary 1m PFAS + 6m Microbiome
make_taxa_heatmap(
  df             = get_species("sens_bin_1m_6m"),
  output_file    = "heatmap_CLR_bin_1m_6m_sensitivity_species.pdf",
  pfas_order     = c(pfas_order_binary, "N-detect"),
  top_n          = 40,
  display_rename = binary_display_rename
)

# Scenario 9: Sensitivity BF Continuous 1m PFAS + 6m Microbiome
make_taxa_heatmap(
  df          = get_species("sens_cont_1m_6m_bf"),
  output_file = "heatmap_CLR_cont_1m_6m_sensitivity_bf_species.pdf",
  pfas_order  = c(pfas_order_continuous, "Mixture"),
  top_n       = 40
)

# Scenario 10: Sensitivity BF Binary 1m PFAS + 6m Microbiome
make_taxa_heatmap(
  df             = get_species("sens_bin_1m_6m_bf"),
  output_file    = "heatmap_CLR_bin_1m_6m_sensitivity_bf_species.pdf",
  pfas_order     = c(pfas_order_binary, "N-detect"),
  top_n          = 40,
  display_rename = binary_display_rename
)


# Beta correlation: Scenario 2 (full) vs Scenario 9 (BF ≥1/day) ---------------
# Load the two combined result files (species level, continuous PFAS only)
cor_full <- read.csv(here::here("out_files", "COMBINED_CLR_continuous_1m_6m.csv")) %>%
  filter(taxa_level == "Species",
         predictor %in% c(pfas_order_continuous, "Mixture")) %>%
  dplyr::select(taxonomy_id, predictor, beta_full = Beta_IQR, p_full = p_val)

cor_bf <- read.csv(here::here("out_files", "COMBINED_CLR_sensitivity_cont_1m_6m_bf.csv")) %>%
  filter(taxa_level == "Species",
         predictor %in% c(pfas_order_continuous, "Mixture")) %>%
  dplyr::select(taxonomy_id, predictor, beta_bf = Beta_IQR, p_bf = p_val)

cor_dat <- inner_join(cor_full, cor_bf, by = c("taxonomy_id", "predictor")) %>%
  filter(!is.na(beta_full), !is.na(beta_bf)) %>%
  mutate(
    predictor = factor(predictor, levels = c(pfas_order_continuous, "Mixture")),
    sig_group = case_when(
      p_full < 0.05 & p_bf < 0.05  ~ "Significant in both",
      p_full < 0.05 & p_bf >= 0.05 ~ "Full only",
      p_full >= 0.05 & p_bf < 0.05 ~ "BF only",
      TRUE                          ~ "Neither"
    ),
    sig_group = factor(sig_group, levels = c("Significant in both",
                                             "Full only", "BF only", "Neither"))
  )

# Compute per PFAS r
r_by_pfas <- cor_dat %>%
  group_by(predictor) %>%
  summarise(r = round(cor(beta_full, beta_bf, use = "complete.obs"), 3),
            .groups = "drop") %>%
  mutate(label = paste0("r = ", r))

# Overall Pearson correlation
r_val <- cor(cor_dat$beta_full, cor_dat$beta_bf, use = "complete.obs")
cat("Overall beta correlation (Scenario 2 vs 9):", round(r_val, 3), "\n")

# Per-PFAS correlations
cor_dat %>%
  group_by(predictor) %>%
  summarise(r = round(cor(beta_full, beta_bf, use = "complete.obs"), 3),
            n = n()) %>%
  print()

# Plot
p_cor <- ggplot(cor_dat, aes(x = beta_full, y = beta_bf, color = sig_group)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey50", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_color_manual(
    values = c("Significant in both" = "#2c7bb6",
               "Full only"           = "#f4a582",
               "BF only"             = "#92c5de",
               "Neither"             = "grey75"),
    name = "Significance"
  ) +
  facet_wrap(~ predictor, ncol = 3) +        # changed from nrow=1 to ncol=3
  geom_text(data = r_by_pfas,
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.2, vjust = 1.5, size = 4.5,
            color = "black", inherit.aes = FALSE) +
  labs(
    x     = "Estimate (IQR-scaled, Full sample)",
    y     = "Estimate (IQR-scaled, Limited to breastfeeding)",
    title = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold", size = 14),
    legend.position  = "bottom",
    legend.text      = element_text(size = 13),
    legend.title     = element_text(size = 14, face = "bold"),
    legend.key.size  = unit(0.5, "cm"),
    axis.title.x     = element_text(face = "bold", color = "black", size = 14),
    axis.title.y     = element_text(face = "bold", color = "black", size = 14),
    axis.text.x      = element_text(color = "black", size = 12),
    axis.text.y      = element_text(color = "black", size = 12)
  )

print(p_cor)
ggsave(here::here("out_figures", "beta_correlation_cont_1m_6m_vs_bf_sensitivity.pdf"),
       plot = p_cor, width = 10, height = 8, device = cairo_pdf, units = "in")
ggsave(here::here("out_figures", "beta_correlation_cont_1m_6m_vs_bf_sensitivity.png"),
       plot = p_cor, width = 10, height = 8, dpi = 600, units = "in")


# Check taxa significant in full but not BF, and vice versa
# Check PFBS significant taxa: full vs BF only
pfbs_check <- cor_dat %>%
  filter(predictor == "PFBS") %>%
  mutate(sig_group = case_when(
    p_full < 0.05 & p_bf < 0.05  ~ "Significant in both",
    p_full < 0.05 & p_bf >= 0.05 ~ "Full only",
    p_full >= 0.05 & p_bf < 0.05 ~ "BF only",
    TRUE                          ~ "Neither"
  )) %>%
  filter(sig_group != "Neither") %>%
  left_join(
    read.csv(here::here("out_files", "COMBINED_CLR_continuous_1m_6m.csv")) %>%
      filter(taxa_level == "Species") %>%
      dplyr::select(taxonomy_id, taxa_name) %>%
      distinct(),
    by = "taxonomy_id"
  ) %>%
  dplyr::select(taxa_name, beta_full, p_full, beta_bf, p_bf, sig_group) %>%
  arrange(sig_group, p_full)

cat("PFBS significant taxa (full vs BF subset):\n")
print(pfbs_check)
cat("Summary:\n")
print(table(pfbs_check$sig_group))

# Beta correlation: Scenario 4 (full binary) vs Scenario 10 (BF ≥0.5/day) -------
cor_full_bin <- read.csv(here::here("out_files", "COMBINED_CLR_binary_1m_6m.csv")) %>%
  mutate(predictor = recode(predictor, "N.MeFOSAA" = "N-MeFOSAA")) %>%
  filter(taxa_level == "Species",
         predictor %in% c(pfas_order_binary, "N-detect")) %>%
  dplyr::select(taxonomy_id, predictor, beta_full = betas, p_full = p_val)

cor_bf_bin <- read.csv(here::here("out_files", "COMBINED_CLR_sensitivity_bin_1m_6m_bf.csv")) %>%
  mutate(predictor = recode(predictor, "N.MeFOSAA" = "N-MeFOSAA")) %>%
  filter(taxa_level == "Species",
         predictor %in% c(pfas_order_binary, "N-detect")) %>%
  dplyr::select(taxonomy_id, predictor, beta_bf = betas, p_bf = p_val)

cor_dat_bin <- inner_join(cor_full_bin, cor_bf_bin, by = c("taxonomy_id", "predictor")) %>%
  filter(!is.na(beta_full), !is.na(beta_bf)) %>%
  mutate(
    predictor = factor(predictor, levels = c(pfas_order_binary, "N-detect")),
    sig_group = case_when(
      p_full < 0.05 & p_bf < 0.05  ~ "Significant in both",
      p_full < 0.05 & p_bf >= 0.05 ~ "Full only",
      p_full >= 0.05 & p_bf < 0.05 ~ "BF only",
      TRUE                          ~ "Neither"
    ),
    sig_group = factor(sig_group, levels = c("Significant in both",
                                             "Full only", "BF only", "Neither"))
  )

# Per-PFAS correlations
r_by_pfas_bin <- cor_dat_bin %>%
  group_by(predictor) %>%
  summarise(r = round(cor(beta_full, beta_bf, use = "complete.obs"), 3),
            .groups = "drop") %>%
  mutate(label = paste0("r = ", r))

# Overall correlation
r_val_bin <- cor(cor_dat_bin$beta_full, cor_dat_bin$beta_bf, use = "complete.obs")
cat("Overall beta correlation binary (Scenario 4 vs 10):", round(r_val_bin, 3), "\n")

# Per-PFAS correlations printed
cor_dat_bin %>%
  group_by(predictor) %>%
  summarise(r = round(cor(beta_full, beta_bf, use = "complete.obs"), 3),
            n = n()) %>%
  print()

# Plot
p_cor_bin <- ggplot(cor_dat_bin, aes(x = beta_full, y = beta_bf, color = sig_group)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey50", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_color_manual(
    values = c("Significant in both" = "#2c7bb6",
               "Full only"           = "#f4a582",
               "BF only"             = "#92c5de",
               "Neither"             = "grey75"),
    name = "Significance"
  ) +
  facet_wrap(~ predictor, ncol = 3) +
  geom_text(data        = r_by_pfas_bin,
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.2, vjust = 1.5, size = 4.5,
            color = "black", inherit.aes = FALSE) +
  labs(
    x     = "Estimate (Full sample)",
    y     = "Estimate (Limited to breastfeeding)",
    title = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold", size = 14),
    legend.position  = "bottom",
    legend.text      = element_text(size = 13),
    legend.title     = element_text(size = 14, face = "bold"),
    legend.key.size  = unit(0.5, "cm"),
    axis.title.x     = element_text(face = "bold", color = "black", size = 14),
    axis.title.y     = element_text(face = "bold", color = "black", size = 14),
    axis.text.x      = element_text(color = "black", size = 12),
    axis.text.y      = element_text(color = "black", size = 12)
  )

print(p_cor_bin)
ggsave(here::here("out_figures", "beta_correlation_bin_1m_6m_vs_bf_sensitivity.pdf"),
       plot = p_cor_bin, width = 10, height = 8, device = cairo_pdf, units = "in")
ggsave(here::here("out_figures", "beta_correlation_bin_1m_6m_vs_bf_sensitivity.png"),
       plot = p_cor_bin, width = 10, height = 8, dpi = 600, units = "in")



# ==============================================================================
# SCATTER PLOTS: FDR-significant Bifidobacterium, Lachnospiraceae, Enterobacter
# Continuous PFAS x 6m Microbiome
# Saves all plots as a single multi-page PDF
# ==============================================================================

library(ggplot2)
library(dplyr)
library(stringr)
library(Cairo)
library(officer)

# Load regression results and analysis dataframe --------------------------------

results_6m <- read.csv(here::here("out_files", "COMBINED_CLR_continuous_1m_6m.csv")) %>%
  filter(taxa_level == "Species",
         FDR < 0.05) %>%
  mutate(
    taxa_name  = gsub("\\[|\\]", "", taxa_name),
    plot_beta  = dplyr::coalesce(as.numeric(Beta_IQR), as.numeric(betas)),
    Genus      = word(taxa_name, 1)
  ) %>%
  filter(Genus %in% c("Bifidobacterium", "Enterobacter",
                      "Blautia", "Dorea", "Anaerostipes", "Roseburia",
                      "Coprococcus", "Ruminococcus", "Lachnoclostridium",
                      "Anaerobutyricum", "Mediterraneibacter", "Simiaoa",
                      "[Clostridium]"))

cat("FDR-significant pairs to plot:", nrow(results_6m), "\n")
cat("Unique taxa:", n_distinct(results_6m$taxa_name), "\n")
cat("Unique PFAS:", n_distinct(results_6m$predictor), "\n")
print(results_6m %>% count(Genus, taxa_name, predictor) %>% arrange(Genus))

# Load analysis dataframe (per-participant PFAS + CLR taxa) --------------------

pfas_vars  <- c("PFBS", "PFHxS", "PFNA", "PFOA", "PFOS")

# Load PFAS data for x-axis (log2-transformed, IQR-scaled)
pfas_df <- read.csv(here::here("out_files", "PFAS1m_micro6m_species.csv")) %>%
  dplyr::rename(
    PFBS  = PFBS_pgmL,
    PFHxS = PFHxS_pgmL,
    PFNA  = PFNA_pgmL,
    PFOA  = PFOA_pgmL,
    PFOS  = PFOS_pgmL
  ) %>%
  dplyr::select(merge_id_dyad, PFBS, PFHxS, PFNA, PFOA, PFOS,
                mode_of_delivery_cat, baby_birthweight_kg, gestational_age_cat,
                breastmilk_per_day, SES_index_final, age_of_solid_foods)

# Load CLR-transformed taxa for y-axis
clr_df <- read.csv(here::here("out_files", "CLR_perSample_species_6m.csv"))
# CLR columns are named X{taxonomy_id}_CLR — strip the _CLR suffix to match taxonomy_id
colnames(clr_df) <- sub("_CLR$", "", colnames(clr_df))

# Merge PFAS + CLR by merge_id_dyad
analysis_df <- inner_join(pfas_df, clr_df, by = "merge_id_dyad")
cat("Merged analysis_df rows:", nrow(analysis_df), "\n")
cat("Sample CLR values for check:\n")
print(summary(analysis_df[, grep("^X", colnames(analysis_df))[1:3]]))

# PFAS IQR values for x-axis scaling ------------------------------------------

pfas_iqr <- sapply(pfas_vars, function(p) IQR(analysis_df[[p]], na.rm = TRUE))
cat("\nPFAS IQRs:\n"); print(round(pfas_iqr, 3))

# Genus group label for plot strip ---------------------------------------------

get_genus_group <- function(taxa_name) {
  genus <- word(taxa_name, 1)
  lachn_genera <- c("Blautia", "Dorea", "Anaerostipes", "Roseburia",
                    "Coprococcus", "Ruminococcus", "Lachnoclostridium",
                    "Anaerobutyricum", "Mediterraneibacter", "Simiaoa",
                    "[Clostridium]")
  if (genus == "Bifidobacterium") return("Bifidobacterium")
  if (genus == "Enterobacter")    return("Enterobacter")
  if (genus %in% lachn_genera)    return("Lachnospiraceae")
  return("Other")
}

# Group colors
group_cols <- c(
  "Bifidobacterium" = "#D55E00",
  "Lachnospiraceae" = "#0072B2",
  "Enterobacter"    = "#009E73"
)

# Scatter plot function --------------------------------------------------------

make_scatter <- function(taxa_col, pfas_predictor, taxa_label,
                         beta_iqr, fdr_val, genus_group, df) {
  
  # Get CLR column name from taxonomy_id
  tax_row <- results_6m %>%
    filter(taxa_name == taxa_label, predictor == pfas_predictor) %>%
    slice(1)
  
  clr_col <- tax_row$taxonomy_id
  
  if (!clr_col %in% colnames(df)) {
    cat("  Column not found:", clr_col, "for", taxa_label, "\n")
    return(NULL)
  }
  
  plot_df <- df %>%
    dplyr::select(all_of(c(pfas_predictor, clr_col))) %>%
    dplyr::rename(PFAS = 1, CLR = 2) %>%
    filter(!is.na(PFAS), !is.na(CLR))
  
  # IQR-scale x axis for display
  iqr_val   <- pfas_iqr[pfas_predictor]
  plot_df   <- plot_df %>% mutate(PFAS_scaled = PFAS / iqr_val)
  
  pt_color  <- group_cols[genus_group]
  fdr_label <- formatC(fdr_val, format = "f", digits = 3)
  b_label   <- formatC(beta_iqr, format = "f", digits = 3)
  
  ggplot(plot_df, aes(x = PFAS_scaled, y = CLR)) +
    geom_point(color = pt_color, alpha = 0.55, size = 1.8) +
    geom_smooth(method = "lm", se = TRUE, color = "black",
                fill = "grey80", linewidth = 0.8) +
    annotate("text",
             x = -Inf, y = Inf,
             label = paste0("beta = ", b_label, "\nFDR = ", fdr_label),
             hjust = -0.15, vjust = 1.3,
             size = 4.5, color = "black") +
    labs(
      x     = paste0("log2-", pfas_predictor, " (IQR-scaled)"),
      y = bquote("CLR(" ~ italic(.(taxa_label)) ~ ")"),
      title = paste0(taxa_label, " - ", pfas_predictor),
      subtitle = genus_group
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      axis.line.x       = element_line(color = "black", linewidth = 0.4),
      axis.line.y       = element_line(color = "black", linewidth = 0.4),
      axis.text         = element_text(color = "black", size = 12),
      axis.title        = element_text(color = "black", size = 13, face = "bold"),
      legend.text       = element_text(size = 12),
      legend.title      = element_text(size = 13, face = "bold"),
      legend.key.size   = unit(0.5, "cm"),
      plot.title        = element_text(size = 11, face = "italic"),
      plot.subtitle     = element_text(size = 10, color = pt_color, face = "bold"),
      plot.margin       = margin(10, 15, 10, 10)
    )
}

# Helper: italicize species name in title (plain text approximation)
italic_name <- function(x) x   # ggplot titles don't support mixed formatting
# use as-is; italics handled by theme(plot.title)

# Generate all plots -----------------------------------------------------------
# Title page and separator functions ------------------------------------------
make_title_page <- function() {
  ggplot() +
    theme_void() +
    annotate("text", x = 0.5, y = 0.65, hjust = 0.5,
             label = "Figure S15. FDR-significant associations between 1-month human milk PFAS\nand infant gut microbiome taxa at 6 months of age",
             size = 5, fontface = "bold") +
    annotate("text", x = 0.5, y = 0.40, hjust = 0.5,
             label = "Scatter plots show CLR-transformed species abundance vs IQR-scaled log2-PFAS.\nOrdered by taxon group: Bifidobacterium (negative), Lachnospiraceae (positive), Enterobacter (positive).\nbeta = IQR-scaled regression coefficient; FDR = Benjamini-Hochberg adjusted p-value.",
             size = 3.5, color = "grey30") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
}

make_separator <- function(label, color) {
  ggplot() +
    theme_void() +
    annotate("text", x = 0.5, y = 0.5, hjust = 0.5,
             label = label, size = 10, fontface = "bold.italic",
             color = color) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
}

plot_list <- list()
if (exists("last_group")) rm(last_group)

# Sort results by genus group order before plotting
results_6m <- results_6m %>%
  mutate(
    genus_group_order = case_when(
      Genus == "Bifidobacterium" ~ 1,
      Genus %in% c("Blautia", "Dorea", "Anaerostipes", "Roseburia",
                   "Coprococcus", "Ruminococcus", "Lachnoclostridium",
                   "Anaerobutyricum", "Mediterraneibacter", "Simiaoa",
                   "[Clostridium]") ~ 2,
      Genus == "Enterobacter" ~ 3,
      TRUE ~ 4
    )
  ) %>%
  arrange(genus_group_order, taxa_name, predictor)

for (i in seq_len(nrow(results_6m))) {
  
  row         <- results_6m[i, ]
  taxa_label  <- row$taxa_name
  pfas_pred   <- row$predictor
  beta_iqr    <- as.numeric(row$plot_beta)
  fdr_val     <- as.numeric(row$FDR)
  genus_group <- get_genus_group(taxa_label)
  
  if (pfas_pred == "Mixture") next   # skip mixture — no single column to plot
  # Track group changes (no separator page added)
  current_group <- get_genus_group(taxa_label)
  last_group <- current_group
  cat(sprintf("  Plotting: %s ~ %s (FDR=%.3f)\n", taxa_label, pfas_pred, fdr_val))
  
  p <- make_scatter(
    taxa_col      = taxa_label,
    pfas_predictor = pfas_pred,
    taxa_label    = taxa_label,
    beta_iqr      = beta_iqr,
    fdr_val       = fdr_val,
    genus_group   = genus_group,
    df            = analysis_df
  )
  
  if (!is.null(p)) plot_list[[length(plot_list) + 1]] <- p
}

cat("\nTotal plots generated:", length(plot_list), "\n")

# Save as multi-page PDF + Word document --------------------------------------
# 4 plots per page (2 columns x 2 rows)

out_pdf  <- here::here("out_figures", "scatter_FDR_bifido_lach_enterobacter_6m.pdf")
out_docx <- here::here("out_figures", "scatter_FDR_bifido_lach_enterobacter_6m.docx")

plots_per_page <- 4
n_pages <- ceiling(length(plot_list) / plots_per_page)

# ---- Save PDF ---------------------------------------------------------------

cairo_pdf(out_pdf, width = 10, height = 8, onefile = TRUE, family = "Arial")

print(make_title_page())

for (pg in seq_len(n_pages)) {
  idx_start <- (pg - 1) * plots_per_page + 1
  idx_end   <- min(pg * plots_per_page, length(plot_list))
  page_plots <- plot_list[idx_start:idx_end]

  while (length(page_plots) < plots_per_page) {
    page_plots[[length(page_plots) + 1]] <- ggplot() + theme_void()
  }

  grid_plot <- cowplot::plot_grid(
    plotlist = page_plots,
    ncol = 2,
    nrow = 2
  )

  print(grid_plot)
}

dev.off()
cat("Saved PDF:", out_pdf, "\n")

# ---- Save Word document directly -------------------------------------------

doc <- read_docx() %>%
  body_set_default_section(
    prop_section(
      page_size    = page_size(orient = "landscape", width = 11, height = 8.5),
      page_margins = page_mar(top = 0.5, bottom = 0.5, left = 0.5, right = 0.5)
    )
  )

doc <- body_add_gg(doc, value = make_title_page(), width = 10, height = 7.5)
doc <- body_add_break(doc)

# Strip separator pages — keep only actual scatter plots (those with ggplot layers)
scatter_only <- plot_list[vapply(plot_list, function(p) length(p$layers) > 0, logical(1))]

n_pages_word <- ceiling(length(scatter_only) / plots_per_page)

for (pg in seq_len(n_pages_word)) {
  idx_start  <- (pg - 1) * plots_per_page + 1
  idx_end    <- min(pg * plots_per_page, length(scatter_only))
  page_plots <- scatter_only[idx_start:idx_end]
  
  while (length(page_plots) < plots_per_page) {
    page_plots[[length(page_plots) + 1]] <- ggplot() + theme_void()
  }
  
  grid_plot <- cowplot::plot_grid(
    plotlist    = page_plots,
    ncol        = 2,
    nrow        = 2,
    rel_widths  = c(1, 1),
    rel_heights = c(1, 1)
  )
  
  doc <- body_add_gg(doc, value = grid_plot, width = 10, height = 7.3)
  
  if (pg < n_pages_word) doc <- body_add_break(doc)
}

print(doc, target = out_docx)
cat("Saved Word document:", out_docx, "\n")






# Save files for supplementary file---------------------------------------------
library(readr)
library(openxlsx)
library(here)

# Load the four files
cont_1m1m <- read_csv(here("out_files", "COMBINED_CLR_continuous_1m_1m.csv"))
cont_1m6m <- read_csv(here("out_files", "COMBINED_CLR_continuous_1m_6m.csv"))
bin_1m1m  <- read_csv(here("out_files", "COMBINED_CLR_binary_1m_1m.csv"))
bin_1m6m  <- read_csv(here("out_files", "COMBINED_CLR_binary_1m_6m.csv"))

# Create workbook
wb <- createWorkbook()

# Add worksheets
addWorksheet(wb, "PFAS1m_Micro1m_Continuous")
writeData(wb, "PFAS1m_Micro1m_Continuous", cont_1m1m)

addWorksheet(wb, "PFAS1m_Micro6m_Continuous")
writeData(wb, "PFAS1m_Micro6m_Continuous", cont_1m6m)

addWorksheet(wb, "PFAS1m_Micro1m_Binary")
writeData(wb, "PFAS1m_Micro1m_Binary", bin_1m1m)

addWorksheet(wb, "PFAS1m_Micro6m_Binary")
writeData(wb, "PFAS1m_Micro6m_Binary", bin_1m6m)

# Save workbook
saveWorkbook(
  wb,
  here("out_files", "Supplemental_Table_T1.xlsx"),
  overwrite = TRUE
)

cat("Saved: out_files/Supplemental_Table_T1.xlsx\n")
