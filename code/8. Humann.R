# header -----------------------------------------------------------------------
#
# TITLE:   8. Pathway Associations - CLR Transformation
#
# PURPOSE: CLR-transformed HUMAnN pathway associations with 1-month PFAS
#          at 1m and 6m microbiome. Uses lm() for individual PFAS,
#          qgcomp.glm.boot() (B=5000) for continuous mixture, and n_detect
#          as binary mixture proxy. Cross-sectional design — no random
#          effects. Mirrors taxa CLR pipeline.
#
# DATE:    March 2026
#
# NOTES:
#   - CPM values are already depth-normalized by HUMAnN → TSS step skipped
#   - CZM zero imputation applied before CLR (same as taxa pipeline)
#   - Prevalence filter (>=25%) applied across all 1m+6m MG samples combined
#     before splitting by timepoint (consistent with taxa pipeline)
#   - Bootstrap qgcomp used with q=4, B=5000, seed=2024, gaussian family
#   - Pathway names converted to valid R column names; original names
#     preserved in output via lookup table
#
# IQR SCALING:
#   - Individual PFAS: raw betas multiplied post-hoc by each PFAS IQR
#     (Beta_IQR = beta x IQR). P-values unchanged.
#   - Mixture (qgcomp): re-run with IQR-scaled inputs (q=NULL) so psi
#     represents change per simultaneous 1-IQR increase across all PFAS.
#     Raw q=4 model also retained (Beta, P).
#   - IQR scaling applied to continuous PFAS only (Scenarios 1-2).
#     Binary scenarios (3-4) are not scaled.
#   - Direction for alluvial plots and pattern analysis uses IQR-scaled
#     beta (plot_beta) for continuous; raw beta for binary.
#          Each Model run time is approximately 12 minutes (especially the continuous with qgcomp)
#          Binary models are less than a minute
#
# Code Review
# Ellie Holzhausen (EAH) on April 27, 2026
# 
#
# set up -----------------------------------------------------------------------
rm(list = ls())
options(scipen = 0)

library(dplyr)
library(readr)
library(readxl)
library(qgcomp)
library(tibble)
library(zCompositions)
library(compositions)
library(parallel)
library(ggalluvial)
library(patchwork)
library(here)
library(circlize)
library(stringr)
library(RColorBrewer)
library(tidyverse)

# Parallelization setup
n_cores <- detectCores() - 2   # leaving 2 cores out of 15

# Bootstrap iterations for qgcomp
N_BOOT <- 5000


# Variable lists
pfas_rename_continuous <- c(
  "PFBS"  = "PFBS_pgmL",
  "PFHxS" = "PFHxS_pgmL",
  "PFNA"  = "PFNA_pgmL",
  "PFOA"  = "PFOA_pgmL",
  "PFOS"  = "PFOS_pgmL"
)

pfas_vars_continuous <- c("PFBS", "PFHxS", "PFNA", "PFOA", "PFOS")
pfas_vars_binary     <- c("N.MeFOSAA", "PFBA", "PFDA", "PFDoA", "PFHpA")

covariates <- c("breastmilk_per_day", "SES_index_final", "baby_birthweight_kg",
                "gest_Early", "gest_Late", "mode_of_delivery_bin")

selectCovariates <- c("merge_id_dyad", covariates)

# Compute IQR of each log2-transformed continuous PFAS-------------------------
compute_pfas_iqr <- function(data) {
  sapply(pfas_vars_continuous, function(pfas) {
    IQR(data[[pfas]], na.rm = TRUE)
  })
}

# Load and process HUMAnN collapsed pathway data--------------------------------

# ── Load IDs manifest to link MG IDs → dyad_id + timepoint
IDs <- read_excel(here::here("input", "GORAN_MICROBIOME_MANIFEST 10-2022.xlsx"))
colnames(IDs) <- IDs[1, ]
IDs <- IDs[-1, ]

# Fix duplicate column names before any operations
colnames(IDs) <- make.unique(colnames(IDs), sep = "_")

IDs$dyad_id       <- substr(IDs$`Old Together`, 9, 11)
IDs$timepoint     <- substr(IDs$`Old Together`, 1, 2)
IDs$dyad_id2      <- str_pad(IDs$dyad_id, width = 4, pad = "0")
IDs$merge_id_dyad <- paste0("MM-", IDs$dyad_id2, "-", IDs$timepoint)
IDs$core_id       <- paste0("MG", IDs$`CORE ID`)

IDs_1m6m <- IDs %>%
  filter(timepoint %in% c("01", "06")) %>%
  dplyr::select(core_id, merge_id_dyad, timepoint)

# Load collapsed HUMAnN file
humann_raw <- read_tsv(
  here::here("input", "mothersmilk_humann_pathabundances_cpm_collapsed.tsv"),
  show_col_types = FALSE
)
colnames(humann_raw)[1] <- "Pathway"

# Keep only MG columns + Pathway column
humann_raw <- humann_raw %>%
  dplyr::select(Pathway, starts_with("MG"))

cat("  Total rows (raw):", nrow(humann_raw), "\n")
cat("  MG samples total:", ncol(humann_raw) - 1, "\n")

# Remove UNMAPPED, UNINTEGRATED, and species-stratified rows
# Species-stratified rows (e.g. "PWY-123|g__Genus.s__Species") are removed
# from the CLR regression input because:
#   1. Analysis targets community-level pathway abundance, not per-species contributions
#   2. Species contributions are analyzed separately via the genus arc ring
#      visualization (avg_contrib / top5_contrib objects below)
#   3. Including stratified rows alongside community totals would introduce
#      collinearity (community total = sum of species contributions)
humann_clean <- humann_raw %>%
  filter(!Pathway %in% c("UNMAPPED", "UNINTEGRATED")) %>%
  filter(!str_detect(Pathway, "\\|"))

# Clean pathway names (keep descriptive part after ":" to avoid collision or being appear as duplicate)
humann_clean <- humann_clean %>%
  mutate(
    Pathway_orig  = Pathway,
    # Keep MetaCyc ID as prefix: "PWY-5659 valine biosynthesis"
    Pathway_clean = sub("^([^:]+):\\s*(.+)$", "\\1 \\2", Pathway),
    Pathway_clean = str_squish(Pathway_clean)
  )

# Verify duplicates
humann_clean %>%
  count(Pathway_clean) %>%
  filter(n > 1)
# Should return 0 rows

# Keep only MG samples from 1m and 6m timepoints
mg_1m6m <- IDs_1m6m$core_id
humann_clean <- humann_clean %>%
  dplyr::select(Pathway_orig, Pathway_clean,
                all_of(intersect(mg_1m6m, colnames(.))))



# Prevalence filter (>=25% across all 1m+6m samples combined)-------------------
# Consistent with taxa pipeline: filter before splitting by timepoint

# Build numeric matrix: rows = pathways, cols = MG samples
pwy_matrix <- humann_clean %>%
  dplyr::select(-Pathway_orig, -Pathway_clean) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

rownames(pwy_matrix) <- humann_clean$Pathway_clean

# Prevalence = proportion of samples with value > 0
prevalence  <- rowMeans(pwy_matrix > 0, na.rm = TRUE)
keep_paths  <- names(prevalence[prevalence >= 0.25])

cat("  Pathways before filter:", nrow(pwy_matrix), "\n")
cat("  Pathways after >=25% filter:", length(keep_paths), "\n")
cat("  Median zero rate after filter:",
    round(median(rowMeans(pwy_matrix[keep_paths, ] == 0)), 3), "\n")

pwy_matrix_filtered <- pwy_matrix[keep_paths, ]

# Pathway name lookup: clean name → valid R column name
pwy_name_map <- data.frame(
  col_name  = make.unique(make.names(keep_paths), sep = "_"),
  orig_name = keep_paths,
  stringsAsFactors = FALSE
)

rownames(pwy_matrix_filtered) <- pwy_name_map$col_name
pwy_cols_clr <- pwy_name_map$col_name

print(head(pwy_name_map, 3))

# Assign_pathway_category
assign_pathway_category <- function(df) {
  df %>%
    mutate(
      Category = case_when(
        
        # Cell Wall
        str_detect(Pathway, regex(
          paste0("peptidoglycan|UDP-N-acetyl|CMP-3-deoxy|ADP-L-glycero|",
                 "colanic acid|O-antigen|anhydromuropeptides|",
                 "lipopolysaccharide|lipid A|lipid IVA|teichoic|teichuronic|murein|",
                 "enterobacterial common antigen|chitin|dTDP|",
                 "CMP-legionaminate|CMP-pseudaminate|polymyxin resistance|chondroitin"),
          ignore_case = TRUE)) ~ "Cell Wall",
        
        # SCFAs
        str_detect(Pathway, regex(
          paste0("fermentation|butanediol|",
                 "butanol biosynthesis|isopropanol biosynthesis|propanediol"),
          ignore_case = TRUE)) ~ "SCFAs",
        
        # Lipids
        str_detect(Pathway, regex(
          paste0("CDP-diacylglycerol|fatty acid|phosphatidylglycerol|",
                 "dodecenoate|cis-vaccenate|gondoate|stearate|palmitoleate|palmitate|",
                 "isoprene|polyisoprenoid|phospholipid|methylerythritol|oleate|",
                 "phospholipases|petroselinate|arachidonate|docosahexaenoate|",
                 "icosapentaenoate|ceramide|sphingolipid|diacylglycerol|triacylglycerol|",
                 "phosphatidate|myristate|plasmalogen|phosphatidylcholine|phosphoinositide|",
                 "monoacylglycerol|enoyl-CoA|ergosterol|acyl-carrier protein|mevalonate|",
                 "terpenoid|carotenoid|hopanoid|cardiolipin|palmitoyl|farnesol|geranyl|",
                 "phytol|ketogenesis|sciadonate|taxadiene|juniperonate|tricosene|",
                 "very long.chain fatty acid|ultra.long.chain fatty acid"),
          ignore_case = TRUE)) ~ "Lipids",
        
        str_detect(Pathway, regex(
          paste0("pentose phosphate|rhamnose|ascorbate|fructuronate|sucrose|galactose|",
                 "glycogen|glucosamine|mannose|gluconeogenesis|glycolysis|stachyose|",
                 "glucose|fucose|D-galacturonate|arabinose|xylose|fructose|fructan|",
                 "anhydrofructose|galactarate|glucarate|galactitol|mannitol|inositol|",
                 "gluconate|Entner-Doudoroff|glucuronoside|glucuronate|hexitol|starch|",
                 "lactose|Bifidobacterium shunt|formaldehyde|methanol oxidation|",
                 "sulfoquinovose|mannan|threonate|enopyranuronate|glycol metabolism|",
                 "C4 photosynthetic carbon assimilation|TCA cycle|glyoxylate cycle|",
                 "glyoxylate bypass|methylglyoxal|hexuronide|hexuronate|",
                 "neuraminate degradation|Calvin-Benson-Bassham|mannosylglycerate|",
                 "erythronate|methylcitrate"),
          ignore_case = TRUE)) ~ "Carbohydrates",
        
        # Note: tetrapyrrole was already in this pattern — reordering alone fixes it
        str_detect(Pathway, regex(
          paste0("coenzyme A|NAD|flavin|pyridoxal|thiamine|folate|cobalamin|menaquinol|",
                 "naphthoquinol|hydroxymethyl dihydropterin|molybdopterin|heme b|",
                 "phosphopantothenate|biotin|oxononanoate|tetrapyrrole|ubiquinol|",
                 "cobyrinate|coenzyme M|factor 420|nicotinate|adenosylcobalamin|",
                 "demethylmenaquinol|naphthoate|dihydropterin"),
          ignore_case = TRUE)) ~ "Vitamins & Cofactors",
        
        # Amino Acids
        str_detect(Pathway, regex(
          paste0("valine|methionine|arginine|cysteine|alanine|serine|glycine|",
                 "tryptophan|threonine|isoleucine|lysine|histidine|ornithine|",
                 "aspartate|asparagine|seleno-amino acid|tyrosine|aromatic amino acid|",
                 "glutamyl cycle|branched chain amino acid|chorismate|phenylalanine|",
                 "leucine|glutamate|glutamine|proline|citrulline|sulfur amino acid|",
                 "urea cycle|polyamine|putrescine|spermidine|norspermidine|homocysteine|",
                 "carnitine|homoserine|aminobutanoate|GABA|creatinine|",
                 "ethanolamine utilization|oxobutanoate|5-oxo-L-proline|shikimate|",
                 "anthranilate|dehydroquinate|acetylornithine"),
          ignore_case = TRUE)) ~ "Amino Acids",
        
        # Nucleotides
        str_detect(Pathway, regex(
          paste0("purine|pyrimidine|inosine|nucleotide|adenine|adenosine|guanosine|",
                 "queuosine|5-aminoimidazole|UMP|ppGpp|preQ|wybutosine|CTP|UTP|",
                 "allantoin|ureide|deoxyribonucleoside|ribonucleoside|nucleobase|",
                 "ribonucleotide"),
          ignore_case = TRUE)) ~ "Nucleotides",
        
        TRUE ~ "Other"
      )
    )
}

# Save pathway-category lookup
pwy_category_lookup <- pwy_name_map %>%
  dplyr::rename(Pathway = orig_name) %>%
  assign_pathway_category() %>%
  dplyr::select(Pathway, Category, col_name)

write.csv(pwy_category_lookup,
          here::here("out_files", "pathway_category_lookup.csv"),
          row.names = FALSE)

# SECTION 4 — CZM imputation + CLR transformation-------------------------------
# No TSS step — CPM values are already depth-normalized by HUMAnN

# Transpose: rows = samples, cols = pathways
pwy_t <- t(pwy_matrix_filtered)
colnames(pwy_t) <- pwy_cols_clr

# CZM zero imputation (same settings as taxa pipeline)
n_zeros <- sum(pwy_t == 0)
cat("  Zeros before CZM:", n_zeros,
    "(", round(mean(pwy_t == 0) * 100, 1), "% )\n")

if (n_zeros > 0) {
  cat("  Running CZM...\n")
  pwy_imputed <- cmultRepl(
    pwy_t, label = 0, method = "CZM",
    z.warning = 1, z.delete = FALSE
  )
  cat("  CZM done. Zeros remaining:", sum(pwy_imputed == 0), "\n")
} else {
  cat("  No zeros — skipping CZM.\n")
  pwy_imputed <- pwy_t
}

# CLR transformation (per-sample geometric mean reference)
pwy_clr        <- t(apply(pwy_imputed, 1, clr))
pwy_clr_df     <- as.data.frame(pwy_clr)
colnames(pwy_clr_df) <- pwy_cols_clr
pwy_clr_df$core_id   <- rownames(pwy_clr_df)



# Link CLR data to merge_id_dyad and split by timepoint-------------------------

pwy_clr_linked <- pwy_clr_df %>%
  left_join(IDs_1m6m, by = "core_id") %>%
  filter(!is.na(merge_id_dyad)) %>%
  dplyr::select(-core_id, -timepoint)

# Note: merge_id_dyad for 1m ends in "-01", for 6m ends in "-06"
# PFAS saved files use the same convention:
#   PFAS1m_micro1m_species.csv → merge_id_dyad ends in "-01"
#   PFAS1m_micro6m_species.csv → merge_id_dyad ends in "-06"
# So no substitution needed — merge directly by merge_id_dyad

pwy_clr_1m <- pwy_clr_linked %>%
  filter(str_detect(merge_id_dyad, "-01$"))

pwy_clr_6m <- pwy_clr_linked %>%
  filter(str_detect(merge_id_dyad, "-06$"))


# Load PFAS + covariate metadata from saved files-------------------------------
# These files already have: log2 PFAS, imputed covariates, n_detect
# Just strip the taxa count columns

load_pfas_meta <- function(filename) {
  read.csv(here::here("out_files", filename)) %>%
    dplyr::select(-matches("^X\\d+")) %>%
    { if ("counts" %in% colnames(.)) dplyr::select(., -counts) else . }
}

# Add binary covariates (gestational_age_bin, mode_of_delivery_bin)
# Mirrors prepare_taxa_data_clr in pfas-taxa analyses
add_binary_covariates <- function(df) {
  df %>% mutate(
    gest_Early           = ifelse(gestational_age_cat  == "Early",     1L, 0L),
    gest_Late            = ifelse(gestational_age_cat  == "Late",      1L, 0L),
    mode_of_delivery_bin = ifelse(mode_of_delivery_cat == "C-Section", 1L, 0L)
  )
}

# Helper to scale continuous covariates — mirrors prepare_taxa_data_clr() in script 5b
add_scale_covariates <- function(df) {
  df %>% mutate(
    breastmilk_per_day  = as.numeric(scale(breastmilk_per_day)),
    SES_index_final     = as.numeric(scale(SES_index_final)),
    baby_birthweight_kg = as.numeric(scale(baby_birthweight_kg))
  )
}

# Continuous PFAS files
pfas_meta_cont_1m <- load_pfas_meta("PFAS1m_micro1m_species.csv") %>%
  dplyr::rename(!!!pfas_rename_continuous) %>%
  add_binary_covariates() %>%
  add_scale_covariates()

pfas_meta_cont_6m <- load_pfas_meta("PFAS1m_micro6m_species.csv") %>%
  dplyr::rename(!!!pfas_rename_continuous) %>%
  add_binary_covariates() %>%
  add_scale_covariates()

# Binary PFAS files
# Rename to N.MeFOSAA etc. to match pfas_vars_binary
pfas_meta_bin_1m <- load_pfas_meta("PFAS1mDetect_micro1m_species.csv") %>%
  dplyr::rename(!!!setNames(paste0(pfas_vars_binary, "_Detect"), pfas_vars_binary)) %>%
  mutate(across(all_of(pfas_vars_binary),
                ~ factor(., levels = c("non-detect", "detect")))) %>%
  add_binary_covariates() %>%
  add_scale_covariates()

pfas_meta_bin_6m <- load_pfas_meta("PFAS1mDetect_micro6m_species.csv") %>%
  dplyr::rename(!!!setNames(paste0(pfas_vars_binary, "_Detect"), pfas_vars_binary)) %>%
  mutate(across(all_of(pfas_vars_binary),
                ~ factor(., levels = c("non-detect", "detect")))) %>%
  add_binary_covariates() %>%
  add_scale_covariates()



# Merge PFAS metadata with CLR pathway data-------------------------------------
data_cont_1m_1m <- inner_join(pfas_meta_cont_1m, pwy_clr_1m, by = "merge_id_dyad")
data_cont_1m_6m <- inner_join(pfas_meta_cont_6m, pwy_clr_6m, by = "merge_id_dyad")
data_bin_1m_1m  <- inner_join(pfas_meta_bin_1m,  pwy_clr_1m, by = "merge_id_dyad")
data_bin_1m_6m  <- inner_join(pfas_meta_bin_6m,  pwy_clr_6m, by = "merge_id_dyad")

write.csv(data_cont_1m_1m, here::here("out_files", "analysis_df_cont_1m_1m.csv"), row.names = FALSE)
write.csv(data_cont_1m_6m, here::here("out_files", "analysis_df_cont_1m_6m.csv"), row.names = FALSE)
write.csv(data_bin_1m_1m,  here::here("out_files", "analysis_df_bin_1m_1m.csv"),  row.names = FALSE)
write.csv(data_bin_1m_6m,  here::here("out_files", "analysis_df_bin_1m_6m.csv"),  row.names = FALSE)

# Prepare Functions-------------------------------------------------------------
# Individual lm per PFAS × pathway
run_lm_pathways <- function(data, pfas_vars, pwy_cols, scenario_name,
                            covs = covariates, pfas_iqr = NULL) {
  
  predictor_list <- character(0); name_list      <- character(0)
  betas_list     <- numeric(0);   ci_lo_list     <- numeric(0)
  ci_hi_list     <- numeric(0);   p_values_list  <- numeric(0)
  n_list         <- numeric(0)
  
  for (pwy in pwy_cols) {
    for (pfas in pfas_vars) {
      
      vars_needed  <- c("merge_id_dyad", covs, pfas, pwy)
      ithDataframe <- na.omit(data[, vars_needed])
      
      formula <- as.formula(
        paste0("`", pwy, "` ~ ", pfas, " + ",
               paste(covs, collapse = " + "))
      )
      
      tryCatch({
        model    <- lm(formula, data = ithDataframe)
        ct       <- summary(model)$coefficients
        ci       <- confint(model)
        coef_row <- grep(paste0("^", gsub(".", "\\.", pfas, fixed = TRUE)),
                         rownames(ct), value = TRUE)[1]
        
        predictor_list <- c(predictor_list, pfas)
        name_list      <- c(name_list,      pwy)
        betas_list     <- c(betas_list,     ct[coef_row, "Estimate"])
        ci_lo_list     <- c(ci_lo_list,     ci[coef_row, 1])
        ci_hi_list     <- c(ci_hi_list,     ci[coef_row, 2])
        p_values_list  <- c(p_values_list,  ct[coef_row, "Pr(>|t|)"])
        n_list         <- c(n_list,         nrow(ithDataframe))
        
      }, error = function(e) {
        warning(paste("lm failed:", pwy, "~", pfas, ":", e$message))
        predictor_list <<- c(predictor_list, pfas)
        name_list      <<- c(name_list,      pwy)
        betas_list     <<- c(betas_list,     NA_real_)
        ci_lo_list     <<- c(ci_lo_list,     NA_real_)
        ci_hi_list     <<- c(ci_hi_list,     NA_real_)
        p_values_list  <<- c(p_values_list,  NA_real_)
        n_list         <<- c(n_list,         NA_integer_)
      })
    }
  }
  
  results <- data.frame(
    predictor = predictor_list, col_name = name_list,
    betas = betas_list, ci_lo = ci_lo_list, ci_hi = ci_hi_list,
    p_val = p_values_list, n = n_list,
    stringsAsFactors = FALSE
  ) %>%
    left_join(pwy_name_map, by = "col_name") %>%
    dplyr::rename(name = orig_name) %>%
    dplyr::select(-col_name) %>%
    group_by(predictor) %>%
    mutate(FDR = p.adjust(p_val, method = "BH")) %>%
    ungroup()
  
  # Post-hoc IQR scaling — only applied when pfas_iqr is provided (continuous)
  if (!is.null(pfas_iqr)) {
    results <- results %>%
      mutate(
        Beta_IQR = ifelse(predictor %in% names(pfas_iqr),
                          betas * pfas_iqr[predictor], NA_real_),
        SE_IQR   = ifelse(predictor %in% names(pfas_iqr),
                          (ci_hi - ci_lo) / (2 * 1.96) * pfas_iqr[predictor],
                          NA_real_)
      )
  }
  
  cat("  Ran:", nrow(results), "| Converged:", sum(!is.na(results$betas)), "\n")
  return(results)
}


# lm with n_detect--------------------------------------------------------------
run_ndetect_pathways <- function(data, pwy_cols, scenario_name,
                                 covs = covariates) {
  
  name_list     <- character(0); betas_list    <- numeric(0)
  ci_lo_list    <- numeric(0);   ci_hi_list    <- numeric(0)
  p_values_list <- numeric(0);   n_list        <- numeric(0)
  
  for (pwy in pwy_cols) {
    
    vars_needed  <- c("merge_id_dyad", covs, "n_detect", pwy)
    ithDataframe <- na.omit(data[, vars_needed])
    
    formula <- as.formula(
      paste0("`", pwy, "` ~ n_detect + ",
             paste(covs, collapse = " + "))
    )
    
    tryCatch({
      model    <- lm(formula, data = ithDataframe)
      ct       <- summary(model)$coefficients
      ci       <- confint(model)
      
      name_list     <- c(name_list,     pwy)
      betas_list    <- c(betas_list,    ct["n_detect", "Estimate"])
      ci_lo_list    <- c(ci_lo_list,    ci["n_detect", 1])
      ci_hi_list    <- c(ci_hi_list,    ci["n_detect", 2])
      p_values_list <- c(p_values_list, ct["n_detect", "Pr(>|t|)"])
      n_list        <- c(n_list,        nrow(ithDataframe))
      
    }, error = function(e) {
      warning(paste("n_detect lm failed:", pwy, ":", e$message))
      name_list     <<- c(name_list,     pwy)
      betas_list    <<- c(betas_list,    NA_real_)
      ci_lo_list    <<- c(ci_lo_list,    NA_real_)
      ci_hi_list    <<- c(ci_hi_list,    NA_real_)
      p_values_list <<- c(p_values_list, NA_real_)
      n_list        <<- c(n_list,        NA_integer_)
    })
  }
  
  data.frame(
    predictor = "N-detect", col_name = name_list,
    betas = betas_list, ci_lo = ci_lo_list, ci_hi = ci_hi_list,
    p_val = p_values_list, n = n_list,
    stringsAsFactors = FALSE
  ) %>%
    left_join(pwy_name_map, by = "col_name") %>%
    dplyr::rename(name = orig_name) %>%
    dplyr::select(-col_name) %>%
    mutate(FDR = p.adjust(p_val, method = "BH"))
}

# Mixture analysis--------------------------------------------------------------
# Runs both raw (q=4) and IQR-scaled (q=NULL) qgcomp in one pass
# NOTE: only called for continuous PFAS scenarios
run_qgcomp_pathways_combined <- function(data, pfas_vars, pfas_iqr, pwy_cols,
                                         scenario_name, covs = covariates,
                                         B = N_BOOT, seed = 2024) {
  cat("  Pathways:", length(pwy_cols), "| B:", B,
      "| Cores:", n_cores, "\n")
  
  # IQR-scaled copy of data
  data_iqr <- data
  for (pfas in pfas_vars) {
    data_iqr[[pfas]] <- data_iqr[[pfas]] / pfas_iqr[pfas]
  }
  
  results_list <- mclapply(seq_along(pwy_cols), function(i) {
    
    pwy              <- pwy_cols[i]
    vars_needed      <- c("merge_id_dyad", covs, pfas_vars, pwy)
    ithDataframe     <- na.omit(data[,     vars_needed])
    ithDataframe_iqr <- na.omit(data_iqr[, vars_needed])
    
    out <- list(
      col_name  = pwy,   betas     = NA_real_,
      ci_lo     = NA_real_, ci_hi  = NA_real_,
      p_val     = NA_real_, Beta_IQR  = NA_real_,
      CI_lo_IQR = NA_real_, CI_hi_IQR = NA_real_,
      n         = nrow(ithDataframe)
    )
    
    if (nrow(ithDataframe) < 20) return(out)
    
    formula <- as.formula(
      paste0("`", pwy, "` ~ ",
             paste(c(pfas_vars, covs), collapse = " + "))
    )
    
    tryCatch({
      # Raw model (q=4)
      qgmod_raw <- qgcomp.glm.boot(
        f = formula, data = ithDataframe, expnms = pfas_vars,
        q = 4, family = gaussian(), B = B, seed = seed, rr = FALSE
      )
      # IQR-scaled model (q=NULL)
      qgmod_iqr <- qgcomp.glm.boot(
        f = formula, data = ithDataframe_iqr, expnms = pfas_vars,
        q = NULL, family = gaussian(), B = B, seed = seed, rr = FALSE
      )
      out <- list(
        col_name  = pwy,
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
      warning(paste("qgcomp combined failed:", pwy, ":", e$message))
    })
    
    return(out)
    
  }, mc.cores = n_cores)
  
  bind_rows(lapply(results_list, as.data.frame)) %>%
    left_join(pwy_name_map, by = "col_name") %>%
    dplyr::rename(name = orig_name) %>%
    dplyr::select(-col_name) %>%
    mutate(
      predictor = "Mixture",
      FDR       = p.adjust(p_val, method = "BH")
    )
}



# SCENARIO 1: Continuous PFAS × 1m Pathway--------------------------------------
pfas_iqr_1m <- compute_pfas_iqr(pfas_meta_cont_1m)
cat("PFAS IQRs (1m data):\n"); print(round(pfas_iqr_1m, 4))

lm_cont_1m_1m     <- run_lm_pathways(
  data_cont_1m_1m, pfas_vars_continuous, pwy_cols_clr,
  "Continuous PFAS x 1m Pathway",
  pfas_iqr = pfas_iqr_1m
)
qgcomp_cont_1m_1m <- run_qgcomp_pathways_combined(
  data_cont_1m_1m, pfas_vars_continuous, pfas_iqr_1m, pwy_cols_clr,
  "Continuous PFAS x 1m Pathway"
)

write.csv(
  bind_rows(lm_cont_1m_1m, qgcomp_cont_1m_1m) %>%
    dplyr::rename(Pathway = name) %>%
    assign_pathway_category() %>%
    dplyr::rename(name = Pathway),
  here::here("out_files", "PATHWAY_CLR_continuous_1m_1m.csv"),
  row.names = FALSE
)

# SCENARIO 2: Continuous PFAS × 6m Pathway--------------------------------------
pfas_iqr_6m <- compute_pfas_iqr(pfas_meta_cont_6m)
cat("PFAS IQRs (6m data):\n"); print(round(pfas_iqr_6m, 4))

lm_cont_1m_6m     <- run_lm_pathways(
  data_cont_1m_6m, pfas_vars_continuous, pwy_cols_clr,
  "Continuous PFAS x 6m Pathway",
  pfas_iqr = pfas_iqr_6m
)
qgcomp_cont_1m_6m <- run_qgcomp_pathways_combined(
  data_cont_1m_6m, pfas_vars_continuous, pfas_iqr_6m, pwy_cols_clr,
  "Continuous PFAS x 6m Pathway"
)

write.csv(
  bind_rows(lm_cont_1m_6m, qgcomp_cont_1m_6m) %>%
    dplyr::rename(Pathway = name) %>%
    assign_pathway_category() %>%
    dplyr::rename(name = Pathway),
  here::here("out_files", "PATHWAY_CLR_continuous_1m_6m.csv"),
  row.names = FALSE
)

# SCENARIO 3: Binary PFAS × 1m Pathway------------------------------------------
lm_bin_1m_1m  <- run_lm_pathways(
  data_bin_1m_1m, pfas_vars_binary, pwy_cols_clr,
  "Binary PFAS x 1m Pathway"
)
ndetect_1m_1m <- run_ndetect_pathways(
  data_bin_1m_1m, pwy_cols_clr,
  "Binary PFAS x 1m Pathway"
)

write.csv(
  bind_rows(lm_bin_1m_1m, ndetect_1m_1m) %>%
    dplyr::rename(Pathway = name) %>%
    assign_pathway_category() %>%
    dplyr::rename(name = Pathway),
  here::here("out_files", "PATHWAY_CLR_binary_1m_1m.csv"),
  row.names = FALSE
)

# SCENARIO 4: Binary PFAS × 6m Pathway------------------------------------------
lm_bin_1m_6m  <- run_lm_pathways(
  data_bin_1m_6m, pfas_vars_binary, pwy_cols_clr,
  "Binary PFAS x 6m Pathway"
)
ndetect_1m_6m <- run_ndetect_pathways(
  data_bin_1m_6m, pwy_cols_clr,
  "Binary PFAS x 6m Pathway"
)

write.csv(
  bind_rows(lm_bin_1m_6m, ndetect_1m_6m) %>%
    dplyr::rename(Pathway = name) %>%
    assign_pathway_category() %>%
    dplyr::rename(name = Pathway),
  here::here("out_files", "PATHWAY_CLR_binary_1m_6m.csv"),
  row.names = FALSE
)

# Pathway Visualization---------------------------------------------------------
# PURPOSE: Circos plots with embedded genus arc ring, direction barchart,
#          species contribution analysis. Replaces separate donut plots —
#          genus contribution is now shown as arc ring inside each circos plot.

rm(list = ls())

select <- dplyr::select

# LOAD SAVED PATHWAY RESULTS (p < 0.05 filtered)-------------------------------
pwy_cont_1m_1m <- read.csv(here::here("out_files",
                                      "PATHWAY_CLR_continuous_1m_1m.csv")) %>%
  dplyr::rename(Pathway = name) %>%
  filter(!is.na(p_val), p_val < 0.05) %>%
  mutate(
    betas     = as.numeric(betas),
    Beta_IQR  = as.numeric(Beta_IQR),
    plot_beta = dplyr::coalesce(Beta_IQR, betas),
    Direction = ifelse(plot_beta > 0, "Positive", "Negative"),
    predictor = recode(predictor, "N.MeFOSAA" = "N-MeFOSAA"),
    Category  = recode(Category,
                       "Vitamins & Cofactors" = "Vit. & Cofactor",
                       "Amino Acids"          = "Amino Acid",
                       "Carbohydrates"        = "Carbohydrate",
                       "Lipids"               = "Lipid",
                       "Nucleotides"          = "Nucleotide",
                       "SCFAs"                = "SCFA"),
    scenario  = "Continuous 1m PFAS + 1m Pathway"
  )

pwy_cont_1m_6m <- read.csv(here::here("out_files",
                                      "PATHWAY_CLR_continuous_1m_6m.csv")) %>%
  dplyr::rename(Pathway = name) %>%
  filter(!is.na(p_val), p_val < 0.05) %>%
  mutate(
    betas     = as.numeric(betas),
    Beta_IQR  = as.numeric(Beta_IQR),
    plot_beta = dplyr::coalesce(Beta_IQR, betas),
    Direction = ifelse(plot_beta > 0, "Positive", "Negative"),
    predictor = recode(predictor, "N.MeFOSAA" = "N-MeFOSAA"),
    Category  = recode(Category,
                       "Vitamins & Cofactors" = "Vit. & Cofactor",
                       "Amino Acids"          = "Amino Acid",
                       "Carbohydrates"        = "Carbohydrate",
                       "Lipids"               = "Lipid",
                       "Nucleotides"          = "Nucleotide",
                       "SCFAs"                = "SCFA"),
    scenario  = "Continuous 1m PFAS + 6m Pathway"
  )

pwy_bin_1m_1m <- read.csv(here::here("out_files",
                                     "PATHWAY_CLR_binary_1m_1m.csv")) %>%
  dplyr::rename(Pathway = name) %>%
  filter(!is.na(p_val), p_val < 0.05) %>%
  mutate(
    Direction = ifelse(betas > 0, "Positive", "Negative"),
    predictor = recode(predictor, "N.MeFOSAA" = "N-MeFOSAA"),
    Category  = recode(Category,
                       "Vitamins & Cofactors" = "Vit. & Cofactor",
                       "Amino Acids"          = "Amino Acid",
                       "Carbohydrates"        = "Carbohydrate",
                       "Lipids"               = "Lipid",
                       "Nucleotides"          = "Nucleotide",
                       "SCFAs"                = "SCFA"),
    scenario  = "Binary 1m PFAS + 1m Pathway"
  )

pwy_bin_1m_6m <- read.csv(here::here("out_files",
                                     "PATHWAY_CLR_binary_1m_6m.csv")) %>%
  dplyr::rename(Pathway = name) %>%
  filter(!is.na(p_val), p_val < 0.05) %>%
  mutate(
    Direction = ifelse(betas > 0, "Positive", "Negative"),
    predictor = recode(predictor, "N.MeFOSAA" = "N-MeFOSAA"),
    Category  = recode(Category,
                       "Vitamins & Cofactors" = "Vit. & Cofactor",
                       "Amino Acids"          = "Amino Acid",
                       "Carbohydrates"        = "Carbohydrate",
                       "Lipids"               = "Lipid",
                       "Nucleotides"          = "Nucleotide",
                       "SCFAs"                = "SCFA"),
    scenario  = "Binary 1m PFAS + 6m Pathway"
  )


# VERIFICATION: Loaded data summary---------------------------------------------
verify_loaded <- function(df, label) {
  cat("\n---", label, "---\n")
  cat("  Total significant associations:", nrow(df), "\n")
  cat("  Unique pathways:               ", n_distinct(df$Pathway), "\n")
  cat("  Unique PFAS predictors:        ", n_distinct(df$predictor), "\n")
  cat("  Direction breakdown:\n")
  print(df %>% count(Direction) %>% as.data.frame())
  cat("  Per-predictor counts:\n")
  print(df %>% count(predictor, Direction) %>%
          pivot_wider(names_from = Direction, values_from = n,
                      values_fill = 0) %>% as.data.frame())
  cat("  Unique pathways per category:\n")
  print(df %>% distinct(Pathway, Category) %>%
          count(Category) %>% arrange(desc(n)) %>% as.data.frame())
}

verify_loaded(pwy_cont_1m_1m, "Continuous 1m PFAS + 1m Microbiome")
verify_loaded(pwy_cont_1m_6m, "Continuous 1m PFAS + 6m Microbiome")
verify_loaded(pwy_bin_1m_1m,  "Binary 1m PFAS + 1m Microbiome")
verify_loaded(pwy_bin_1m_6m,  "Binary 1m PFAS + 6m Microbiome")


# SHARED COLORS AND THEMES------------------------------------------------------
cat_cols <- c(
  "Amino Acid"          = "#0072B2",
  "Carbohydrate"        = "#E69F00",
  "Lipid"               = "#009E73",
  "Nucleotide"          = "#CC79A7",
  "Vit. & Cofactor" = "#56B4E9",
  "Cell Wall"            = "#D55E00",
  "SCFA"                = "#F0E442",
  "Other"                = "#999999"
)

# Direction colors
col_pos <- "#0072B2"   # blue  — positive associations
col_neg <- "#D55E00"   # vermillion — negative associations

pathway_barchart_theme <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white", color = NA),
      axis.line.x        = element_line(color = "grey30", linewidth = 0.5),
      axis.line.y        = element_line(color = "grey30", linewidth = 0.5),
      axis.text.x        = element_text(size = 20, angle = 30, hjust = 1,
                                        color = "black", face = "bold"),
      axis.text.y        = element_text(size = 20, colour = "black"),
      axis.title.y       = element_text(size = 22, face = "bold",
                                        margin = margin(r = 15)),
      axis.ticks.x       = element_line(color = "black", linewidth = 0.5),
      axis.ticks.y       = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length  = unit(0.2, "cm"),
      strip.text         = element_text(size = 12, face = "bold", color = "white"),
      strip.background   = element_rect(fill = "black", color = NA),
      legend.position    = "right",
      legend.title       = element_text(size = 12, face = "bold"),
      legend.text        = element_text(size = 10),
      legend.key.size    = unit(0.4, "cm"),
      plot.margin        = margin(10, 10, 10, 10)
    )
}

suspect_genera <- c(
  "Firmicutes", "Bacteroidetes", "Proteobacteria", "Actinobacteria",
  "Verrucomicrobia", "Fusobacteria", "Spirochaetes", "Tenericutes",
  "Clostridia", "Bacilli", "Gammaproteobacteria", "Alphaproteobacteria",
  "Betaproteobacteria", "Deltaproteobacteria", "Epsilonproteobacteria"
)


# Species contribution: load stratified HUMAnN file----------------------------
# Used for genus arc ring inside circos plots
# Weighted score = sum(mean_rel_contrib) per genus across direction-filtered
humann_strat_full <- read_tsv(
  here::here("input", "mothersmilk_humann_pathabundances_cpm.tsv"),
  show_col_types = FALSE
)
colnames(humann_strat_full)[1] <- "Pathway_raw"

community_rows <- humann_strat_full %>%
  filter(!str_detect(Pathway_raw, "\\|")) %>%
  filter(!Pathway_raw %in% c("UNMAPPED", "UNINTEGRATED"))

stratified_rows <- humann_strat_full %>%
  filter(str_detect(Pathway_raw, "\\|")) %>%
  filter(!str_detect(Pathway_raw, "\\|unclassified"))

cat("  Community rows:", nrow(community_rows), "\n")
cat("  Stratified rows:", nrow(stratified_rows), "\n")

clean_pathway_strat <- function(s) {
  s <- sub("^([^:]+):\\s*(.+)$", "\\1 \\2", s)
  str_squish(s)
}

meta_cols <- c("Pathway_raw", "pathway_part", "species_part", "pathway_clean",
               "species_clean", "genus_part", "species_epithet",
               "species_epithet_clean", "species_name", "genus")

stratified_parsed <- stratified_rows %>%
  mutate(
    pathway_part          = str_split_fixed(Pathway_raw, "\\|", 2)[, 1],
    species_part          = str_split_fixed(Pathway_raw, "\\|", 2)[, 2],
    pathway_clean         = clean_pathway_strat(pathway_part),
    species_clean         = str_remove(species_part, "^g__"),
    genus_part            = str_split_fixed(species_clean, "\\.", 2)[, 1],
    species_epithet       = str_remove(
      str_split_fixed(species_clean, "\\.s__", 2)[, 2], "^s__"
    ),
    species_epithet_clean = ifelse(
      str_detect(species_epithet, paste0("^", genus_part, "_")),
      str_remove(species_epithet, paste0("^", genus_part, "_")),
      species_epithet
    ),
    species_name = paste(
      gsub("_", " ", genus_part),
      gsub("_", " ", species_epithet_clean)
    ) %>% str_squish(),
    genus = gsub("_", " ", genus_part)
  ) %>%
  filter(!genus %in% suspect_genera)

sample_cols      <- setdiff(colnames(stratified_parsed), meta_cols)
comm_sample_cols <- intersect(sample_cols, colnames(community_rows))

# Build species CPM matrix (rows = species-pathway pairs, cols = samples)
species_mat <- stratified_parsed %>%
  dplyr::select(all_of(sample_cols)) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

# Build community total CPM matrix (rows = pathways, cols = samples)
community_clean <- community_rows %>%
  mutate(pathway_clean = clean_pathway_strat(Pathway_raw))

total_mat <- community_clean %>%
  dplyr::select(pathway_clean, all_of(comm_sample_cols)) %>%
  mutate(across(-pathway_clean, as.numeric)) %>%
  { m <- as.matrix(.[, -1]); rownames(m) <- .[["pathway_clean"]]; m }

# Align community totals to each species row and compute mean rel contrib
path_idx      <- match(stratified_parsed$pathway_clean, rownames(total_mat))
total_aligned <- total_mat[path_idx, , drop = FALSE]

# Element-wise: species_cpm / total_cpm (NA where total = 0)
valid_mask <- total_aligned > 0 & !is.na(total_aligned)
ratios     <- matrix(NA_real_, nrow = nrow(species_mat), ncol = ncol(species_mat))
ratios[valid_mask] <- species_mat[valid_mask] / total_aligned[valid_mask]

# mean_rel_contrib = mean(rel_contrib) across samples
mean_rc <- rowMeans(ratios, na.rm = TRUE)

# avg_contrib: mean relative contribution per species-pathway pair
avg_contrib <- stratified_parsed %>%
  dplyr::select(pathway_clean, species_name, genus) %>%
  mutate(mean_rel_contrib = mean_rc) %>%
  filter(!is.na(mean_rel_contrib), mean_rel_contrib > 0) %>%
  group_by(pathway_clean, species_name, genus) %>%
  summarise(mean_rel_contrib = mean(mean_rel_contrib, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(pathway_clean, desc(mean_rel_contrib))

cat("  Species-pathway pairs:", nrow(avg_contrib), "\n")

# Top 5 species per pathway by mean_rel_contrib
top5_contrib <- avg_contrib %>%
  group_by(pathway_clean) %>%
  slice_max(mean_rel_contrib, n = 5, with_ties = FALSE) %>%
  ungroup()

cat("Pathway name match check:\n")
cat("  Pathways in regression results:",
    n_distinct(pwy_cont_1m_1m$Pathway), "\n")
cat("  Pathways in stratified HUMAnN:",
    n_distinct(top5_contrib$pathway_clean), "\n")
cat("  Pathways that match:",
    sum(unique(pwy_cont_1m_1m$Pathway) %in%
          unique(top5_contrib$pathway_clean)), "\n")


# Genus weighted score function-------------------------------------------------
# Computes global genus scores and per-category proportions for arc ring
# direction_filter: "Negative" for 1m (dominant), "Positive" for 6m (dominant)
# Score = sum(mean_rel_contrib) per genus across direction-filtered sig pathways
genus_global_top5 <- function(sig_df, direction_filter) {
  sig_filtered <- sig_df %>% filter(Direction == direction_filter)
  
  joined <- sig_filtered %>%
    left_join(top5_contrib, by = c("Pathway" = "pathway_clean"),
              relationship = "many-to-many") %>%
    filter(!is.na(species_name))
  
  # Global genus scores
  genus_global <- joined %>%
    group_by(genus) %>%
    summarise(score = sum(mean_rel_contrib, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(score))
  
  top5 <- genus_global %>%
    slice_max(score, n = 5, with_ties = FALSE) %>%
    pull(genus)
  
  # Per-category proportions for arc ring coloring
  per_cat <- joined %>%
    mutate(genus_grp = ifelse(genus %in% top5, genus, "Other")) %>%
    group_by(Category, genus_grp) %>%
    summarise(score = sum(mean_rel_contrib, na.rm = TRUE), .groups = "drop") %>%
    group_by(Category) %>%
    mutate(total = sum(score), pct = score / total) %>%
    ungroup()
  
  list(top5 = top5, per_cat = per_cat, global = genus_global)
}

# Compute genus data for all 4 scenarios
# 1m: dominant direction = Negative; 6m: dominant direction = Positive
genus_cont_1m <- genus_global_top5(pwy_cont_1m_1m, "Negative")
genus_cont_6m <- genus_global_top5(pwy_cont_1m_6m, "Positive")
genus_bin_1m  <- genus_global_top5(pwy_bin_1m_1m,  "Negative")
genus_bin_6m  <- genus_global_top5(pwy_bin_1m_6m,  "Positive")

cat("  Continuous 1m top 5 genera (negative):",
    paste(genus_cont_1m$top5, collapse = ", "), "\n")
cat("  Continuous 6m top 5 genera (positive):",
    paste(genus_cont_6m$top5, collapse = ", "), "\n")
cat("  Binary 1m top 5 genera (negative):",
    paste(genus_bin_1m$top5,  collapse = ", "), "\n")
cat("  Binary 6m top 5 genera (positive):",
    paste(genus_bin_6m$top5,  collapse = ", "), "\n")

# Genus color assignments
warm_pal <- colorRampPalette(c("#8B0000", "#C0392B", "#E87474", "#F5AAAA", "#FDDCDC"))
cool_pal <- colorRampPalette(c("#003580", "#1565C0", "#4A90D9", "#90C4F0", "#C8E6FF"))

make_genus_col_map <- function(top5_genera, palette_fn) {
  cols <- palette_fn(length(top5_genera))
  setNames(cols, top5_genera)
}

# Overlap % helper for circos subtitle
overlap_pct <- function(sig_this, sig_other) {
  keys_this  <- sig_this  %>% tidyr::unite("k", predictor, Pathway) %>% pull(k) %>% unique()
  keys_other <- sig_other %>% tidyr::unite("k", predictor, Pathway) %>% pull(k) %>% unique()
  pct_unique <- round(sum(!keys_this %in% keys_other) / length(keys_this) * 100)
  pct_shared <- 100 - pct_unique
  list(unique = pct_unique, shared = pct_shared)
}


# CIRCOS PLOT with genus arc ring-----------------------------------------------

make_circos_plot <- function(sig_this, sig_other, genus_data, genus_col_map,
                             pfas_order, tp_label, title_str, output_file,
                             genus_legend_title) {
  
  # Edge list: PFAS -> Category, by direction
  edge_df <- sig_this %>%
    group_by(predictor, Category, Direction) %>%
    summarise(value = n(), .groups = "drop") %>%
    filter(!is.na(Category))
  
  if (nrow(edge_df) == 0) {
    cat("  No data for circos:", title_str, "\n")
    return(invisible(NULL))
  }
  
  pfas_in   <- intersect(pfas_order, unique(edge_df$predictor))
  cats_in   <- sort(unique(edge_df$Category))
  sec_ord   <- c(pfas_in, cats_in)
  n_pfas    <- length(pfas_in)
  n_cat     <- length(cats_in)
  
  # PFAS sectors: black for individual PFAS, grey for Mixture/N-detect
  # Pathway category sectors: orange shades
  mixture_predictors <- c("Mixture", "N-detect")
  
  pfas_cols <- setNames(
    ifelse(pfas_in %in% mixture_predictors, "grey60", "black"),
    pfas_in
  )
  
  # Orange shades for pathway categories
  n_cats <- length(cats_in)
  cat_sector_cols <- setNames(rep("#1B5E20", n_cats), cats_in)
  
  grid_cols <- c(pfas_cols, cat_sector_cols)
  names(grid_cols) <- c(pfas_in, cats_in)
  
  # Ribbon colors by direction
  link_cols <- ifelse(edge_df$Direction == "Positive",
                      adjustcolor(col_pos, 0.50),
                      adjustcolor(col_neg, 0.50))
  
  # Subtitle: overlap context
  ov       <- overlap_pct(sig_this, sig_other)
  other_tp <- if (tp_label == "1m") "6m" else "1m"
  subtitle <- sprintf(
    "%d%% of PFAS-pathway pairs unique to %s microbiome  |  %d%% also significant at %s",
    ov$unique, tp_label, ov$shared, other_tp
  )
  
  dir_label <- if (tp_label == "1m") "negative" else "positive"
  top5      <- genus_data$top5
  per_cat   <- genus_data$per_cat
  
  pdf(here::here("out_figures", output_file), width = 14, height = 13)
  
  circos.clear()
  circos.par(
    canvas.xlim = c(-1.1, 2.3),
    canvas.ylim = c(-1.6, 1.6),
    gap.after   = c(rep(3, n_pfas - 1), 16,
                    rep(3, n_cat  - 1), 16),
    start.degree = 90,
    points.overflow.warning = FALSE
  )
  
  chordDiagram(
    edge_df %>% dplyr::select(from = predictor, to = Category, value),
    order             = sec_ord,
    grid.col          = grid_cols,
    col               = link_cols,
    transparency      = 0,
    annotationTrack   = "grid",
    preAllocateTracks = list(
      list(track.height = 0.08),   # track 1: sector labels
      list(track.height = 0.11)    # track 2: genus arc ring
    )
  )
  
  # Track 2 — genus arc ring (global top 5 genera, proportioned per category)
  circos.trackPlotRegion(track.index = 2, bg.border = NA,
                         panel.fun = function(x, y) {
                           sector <- get.cell.meta.data("sector.index")
                           xlim   <- get.cell.meta.data("xlim")
                           ylim   <- get.cell.meta.data("ylim")
                           if (!sector %in% cats_in) return(invisible(NULL))
                           gd <- per_cat %>% filter(Category == sector)
                           if (nrow(gd) == 0) return(invisible(NULL))
                           # Top5 genera first, Other last
                           gd_top <- gd %>% filter(genus_grp != "Other") %>%
                             mutate(ord = match(genus_grp, top5)) %>% arrange(ord)
                           gd_oth  <- gd %>% filter(genus_grp == "Other")
                           gd      <- bind_rows(gd_top, gd_oth)
                           total_w <- xlim[2] - xlim[1]
                           x_cur   <- xlim[1]
                           for (i in seq_len(nrow(gd))) {
                             g     <- gd$genus_grp[i]
                             x_end <- x_cur + gd$pct[i] * total_w
                             col   <- if (g == "Other") "grey88" else genus_col_map[g]
                             circos.rect(x_cur, ylim[1], x_end, ylim[2], col = col, border = NA)
                             x_cur <- x_end
                           }
                         }
  )
  
  # Track 1 — sector labels
  circos.trackPlotRegion(track.index = 1, bg.border = NA,
                         panel.fun = function(x, y) {
                           xlim <- get.cell.meta.data("xlim")
                           ylim <- get.cell.meta.data("ylim")
                           circos.text(mean(xlim), ylim[1],
                                       get.cell.meta.data("sector.index"),
                                       facing = "clockwise", niceFacing = TRUE,
                                       adj = c(0, 0.5), cex = 1.4, font = 2)
                         }
  )
  
  mtext(subtitle, side = 1, line = 2.2, cex = 0.98, col = "grey30")
  title(title_str, cex.main = 1.5, font.main = 2, line = 0.5)
  
  # Legends (direction + genus) stacked on right side
  leg_x <- 1.20
  
  leg1 <- legend(
    x = leg_x + 0.15, y = 0.72, xjust = 0, yjust = 1,
    legend = c("Positive", "Negative"),
    pch    = 15,
    col    = c(adjustcolor(col_pos, 0.70), adjustcolor(col_neg, 0.70)),
    border = NA, bty = "n", cex = 1.70, pt.cex = 3.0,
    title      = "Association direction",
    title.font = 2, title.cex = 1.70,
    x.intersp  = 0.6, y.intersp = 1.3,
    text.width = 0.54
  )
  
  leg2_y <- leg1$rect$top - leg1$rect$h - 0.05
  legend(
    x = leg_x, y = leg2_y, xjust = 0, yjust = 1,
    legend = c(sapply(top5, function(g) as.expression(bquote(italic(.(g))))),
               list("Other genera")),
    pch    = 15,
    col    = c(genus_col_map[top5], "grey88"),
    border = NA, bty = "n", cex = 1.70, pt.cex = 3.0,
    title      = genus_legend_title,
    title.font = 2, title.cex = 1.70,
    x.intersp  = 0.6, y.intersp = 1.3,
    text.width = 0.63
  )
  
  mtext(
    "* Weighted score: sum of mean relative species contribution\n  across direction-filtered significant pathways",
    side = 1, line = 4.5, cex = 0.80, col = "grey45", adj = 0
  )
  
  circos.clear()
  dev.off()
  cat("  Saved:", output_file, "\n")
}


# Generate 4 circos plots-------------------------------------------------------

# PLOT 1: Continuous PFAS x 1m Pathway
make_circos_plot(
  sig_this      = pwy_cont_1m_1m,
  sig_other     = pwy_cont_1m_6m,
  genus_data    = genus_cont_1m,
  genus_col_map = make_genus_col_map(genus_cont_1m$top5, warm_pal),
  pfas_order = c("PFBS", "PFHxS", "PFNA", "PFOA", "PFOS", "Mixture"),
  genus_legend_title = "Top 5 genera by weighted score\n(1-month negative pathways)",
  tp_label      = "1m",
  title_str     = NULL,
  output_file   = "circos_pwy_continuous_1m_1m.pdf"
)

# PLOT 2: Continuous PFAS x 6m Pathway
make_circos_plot(
  sig_this      = pwy_cont_1m_6m,
  sig_other     = pwy_cont_1m_1m,
  genus_data    = genus_cont_6m,
  genus_col_map = make_genus_col_map(genus_cont_6m$top5, cool_pal),
  pfas_order = c("PFBS", "PFHxS", "PFNA", "PFOA", "PFOS", "Mixture"),
  genus_legend_title = "Top 5 genera by weighted score\n(6-month positive pathways)",
  tp_label      = "6m",
  title_str     = NULL,
  output_file   = "circos_pwy_continuous_1m_6m.pdf"
)

# PLOT 3: Binary PFAS x 1m Pathway
make_circos_plot(
  sig_this      = pwy_bin_1m_1m,
  sig_other     = pwy_bin_1m_6m,
  genus_data    = genus_bin_1m,
  genus_col_map = make_genus_col_map(genus_bin_1m$top5, warm_pal),
  pfas_order = c("N-MeFOSAA", "PFBA", "PFDA", "PFDoA", "PFHpA", "N-detect"),
  genus_legend_title = "Top 5 genera by weighted score\n(1-month negative pathways)",
  tp_label      = "1m",
  title_str     = NULL,
  output_file   = "circos_pwy_binary_1m_1m.pdf"
)

# PLOT 4: Binary PFAS x 6m Pathway
make_circos_plot(
  sig_this      = pwy_bin_1m_6m,
  sig_other     = pwy_bin_1m_1m,
  genus_data    = genus_bin_6m,
  genus_col_map = make_genus_col_map(genus_bin_6m$top5, cool_pal),
  pfas_order = c("N-MeFOSAA", "PFBA", "PFDA", "PFDoA", "PFHpA", "N-detect"),
  genus_legend_title = "Top 5 genera by weighted score\n(6-month positive pathways)",
  tp_label      = "6m",
  title_str     = NULL,
  output_file   = "circos_pwy_binary_1m_6m.pdf"
)


# DIRECTION PATTERN ANALYSIS----------------------------------------------------
# Classifies pathways that reached significance at each timepoint into patterns:
# 1m-only, 6m-only, consistent (same direction both), flipped (opposite direction)

build_pattern_bars <- function(sig_1m, sig_6m, model_label) {
  
  keys_1m <- sig_1m %>% tidyr::unite("key", predictor, Pathway) %>% pull(key)
  keys_6m <- sig_6m %>% tidyr::unite("key", predictor, Pathway) %>% pull(key)
  
  only_1m_keys <- setdiff(keys_1m, keys_6m)
  only_6m_keys <- setdiff(keys_6m, keys_1m)
  shared_keys  <- intersect(keys_1m, keys_6m)
  
  s1m_only <- sig_1m %>%
    tidyr::unite("key", predictor, Pathway, remove = FALSE) %>%
    filter(key %in% only_1m_keys)
  s6m_only <- sig_6m %>%
    tidyr::unite("key", predictor, Pathway, remove = FALSE) %>%
    filter(key %in% only_6m_keys)
  
  shared_df <- sig_1m %>%
    tidyr::unite("key", predictor, Pathway, remove = FALSE) %>%
    filter(key %in% shared_keys) %>%
    dplyr::select(key, Direction) %>%
    dplyr::rename(dir_1m = Direction) %>%
    inner_join(
      sig_6m %>%
        tidyr::unite("key", predictor, Pathway, remove = FALSE) %>%
        filter(key %in% shared_keys) %>%
        dplyr::select(key, Direction) %>%
        dplyr::rename(dir_6m = Direction),
      by = "key"
    )
  
  n_same <- sum(shared_df$dir_1m == shared_df$dir_6m)
  n_flip <- sum(shared_df$dir_1m != shared_df$dir_6m)
  
  tibble(
    category = factor(
      c("1m-only\nnegative", "1m-only\npositive",
        "Both\n(consistent)", "Both\n(flipped)",
        "6m-only\nnegative", "6m-only\npositive"),
      levels = c("1m-only\nnegative", "1m-only\npositive",
                 "Both\n(consistent)", "Both\n(flipped)",
                 "6m-only\nnegative", "6m-only\npositive")
    ),
    count = c(
      sum(s1m_only$Direction == "Negative"),
      sum(s1m_only$Direction == "Positive"),
      n_same, n_flip,
      sum(s6m_only$Direction == "Negative"),
      sum(s6m_only$Direction == "Positive")
    ),
    fill_col = c("#c62828", "#ef9a9a", "#4caf50", "#ff9800",
                 "#90caf9", "#1565c0"),
    model = model_label
  )
}

pattern_cont <- build_pattern_bars(pwy_cont_1m_1m, pwy_cont_1m_6m,
                                   "Continuous PFAS")
pattern_bin  <- build_pattern_bars(pwy_bin_1m_1m,  pwy_bin_1m_6m,
                                   "Binary PFAS")


# VERIFICATION: Direction pattern counts----------------------------------------
verify_patterns <- function(pattern_df, sig_1m, sig_6m, label) {
  cat("\n---", label, "---\n")
  cat("  Pattern counts:\n")
  print(pattern_df %>% dplyr::select(category, count) %>% as.data.frame())
  
  total_1m_only <- sum(pattern_df$count[pattern_df$category %in%
                                          c("1m-only\nnegative", "1m-only\npositive")])
  total_6m_only <- sum(pattern_df$count[pattern_df$category %in%
                                          c("6m-only\nnegative", "6m-only\npositive")])
  total_both    <- sum(pattern_df$count[pattern_df$category %in%
                                          c("Both\n(consistent)", "Both\n(flipped)")])
  
  keys_1m  <- sig_1m %>% tidyr::unite("key", predictor, Pathway) %>% pull(key)
  keys_6m  <- sig_6m %>% tidyr::unite("key", predictor, Pathway) %>% pull(key)
  n_shared <- length(intersect(keys_1m, keys_6m))
  
  cat("  1m-only total:              ", total_1m_only, "\n")
  cat("  6m-only total:              ", total_6m_only, "\n")
  cat("  Shared (both timepoints):   ", total_both,
      "| Expected from key overlap:", n_shared, "\n")
  cat("  Total 1m sig associations:  ", nrow(sig_1m), "\n")
  cat("  Total 6m sig associations:  ", nrow(sig_6m), "\n")
  expected_1m <- total_1m_only + total_both
  expected_6m <- total_6m_only + total_both
  cat("  Reconstructed 1m count:     ", expected_1m,
      ifelse(expected_1m == nrow(sig_1m), " OK", " WARNING: mismatch"), "\n")
  cat("  Reconstructed 6m count:     ", expected_6m,
      ifelse(expected_6m == nrow(sig_6m), " OK", " WARNING: mismatch"), "\n")
}

verify_patterns(pattern_cont, pwy_cont_1m_1m, pwy_cont_1m_6m, "Continuous PFAS")
verify_patterns(pattern_bin,  pwy_bin_1m_1m,  pwy_bin_1m_6m,  "Binary PFAS")


# Direction pattern bar charts--------------------------------------------------

# Continuous
p_pattern_cont <- ggplot(pattern_cont,
                         aes(x = category, y = count, fill = fill_col)) +
  geom_col(color = "black", linewidth = 0.4, width = 0.7) +
  geom_text(aes(label = count), vjust = -0.5, fontface = "bold", size = 9) +
  scale_fill_identity() +
  labs(x = NULL,
       y = "Number of\nPFAS-Pathway Associations",
       title = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  pathway_barchart_theme()
print(p_pattern_cont)
ggsave(here::here("out_figures", "direction_pattern_continuous.pdf"),
       plot = p_pattern_cont, width = 9.5, height = 6,
       device = cairo_pdf, units = "in")

# Binary
p_pattern_bin <- ggplot(pattern_bin,
                        aes(x = category, y = count, fill = fill_col)) +
  geom_col(color = "black", linewidth = 0.4, width = 0.7) +
  geom_text(aes(label = count), vjust = -0.5, fontface = "bold", size = 9) +
  scale_fill_identity() +
  labs(x = NULL,
       y = "Number of\nPFAS-Pathway Associations",
       title = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  pathway_barchart_theme()
print(p_pattern_bin)
ggsave(here::here("out_figures", "direction_pattern_binary.pdf"),
       plot = p_pattern_bin, width = 9.5, height = 6,
       device = cairo_pdf, units = "in")


# VERIFICATION: Species contribution summary------------------------------------
verify_species <- function(sig_1m, sig_6m, label) {
  cat("\n---", label, "---\n")
  all_sig_pathways <- unique(c(sig_1m$Pathway, sig_6m$Pathway))
  matched <- sum(all_sig_pathways %in% unique(top5_contrib$pathway_clean))
  cat("  Total unique sig pathways:   ", length(all_sig_pathways), "\n")
  cat("  Pathways matched to species: ", matched,
      "(", round(matched / length(all_sig_pathways) * 100, 1), "% )\n")
  unmatched <- all_sig_pathways[!all_sig_pathways %in%
                                  unique(top5_contrib$pathway_clean)]
  if (length(unmatched) > 0) {
    cat("  WARNING:", length(unmatched), "pathways have no species match:\n")
    print(head(unmatched, 5))
  } else {
    cat("  OK: all significant pathways matched to species\n")
  }
}

verify_species(pwy_cont_1m_1m, pwy_cont_1m_6m, "Continuous PFAS species contribution")
verify_species(pwy_bin_1m_1m,  pwy_bin_1m_6m,  "Binary PFAS species contribution")

# END

# TOP PATHWAYS SUMMARY FOR GRAPHICAL ABSTRACT ----------------------------------
# Extracts top 5 pathways by direction × timepoint × model type
# Uses |plot_beta| for continuous (IQR-scaled), |betas| for binary

get_top_pathways <- function(df, model_type, timepoint, n = 5) {
  
  beta_col <- if (model_type == "Continuous") "plot_beta" else "betas"
  
  df %>%
    filter(predictor != "Mixture") %>%          # exclude mixture for readability
    mutate(abs_beta = abs(.data[[beta_col]])) %>%
    group_by(Direction) %>%
    slice_max(abs_beta, n = n, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      Model     = model_type,
      Timepoint = timepoint,
      # Shorten pathway names: keep MetaCyc ID + first ~40 chars of description
      Pathway_short = str_trunc(Pathway, 55, side = "right"),
      beta_display  = round(.data[[beta_col]], 3)
    ) %>%
    dplyr::select(Model, Timepoint, Direction, Pathway_short,
                  Pathway, predictor, beta_display, p_val, FDR, Category)
}

top_cont_1m <- get_top_pathways(pwy_cont_1m_1m, "Continuous", "1m")
top_cont_6m <- get_top_pathways(pwy_cont_1m_6m, "Continuous", "6m")
top_bin_1m  <- get_top_pathways(pwy_bin_1m_1m,  "Binary",     "1m")
top_bin_6m  <- get_top_pathways(pwy_bin_1m_6m,  "Binary",     "6m")

top_all <- bind_rows(top_cont_1m, top_cont_6m, top_bin_1m, top_bin_6m)

# Print a clean summary to console --------------------------------------------
cat("\n========== TOP PATHWAYS FOR GRAPHICAL ABSTRACT ==========\n")

for (mod in c("Continuous", "Binary")) {
  for (tp in c("1m", "6m")) {
    for (dir in c("Negative", "Positive")) {
      
      sub_df <- top_all %>%
        filter(Model == mod, Timepoint == tp, Direction == dir)
      
      if (nrow(sub_df) == 0) next
      
      cat(sprintf("\n--- %s PFAS | %s Microbiome | %s Associations ---\n",
                  mod, tp, dir))
      
      sub_df %>%
        dplyr::select(Pathway_short, Category, beta_display, predictor) %>%
        arrange(desc(abs(beta_display))) %>%
        print(n = 5)
    }
  }
}

# Export to CSV ---------------------------------------------------------------
write.csv(
  top_all %>% arrange(Model, Timepoint, Direction, desc(abs(beta_display))),
  here::here("out_files", "top_pathways_graphical_abstract.csv"),
  row.names = FALSE
)
cat("\nSaved: out_files/top_pathways_graphical_abstract.csv\n")




# ==============================================================================
# CONSISTENT PATHWAYS ANALYSIS
# PURPOSE: Extract PFAS-pathway pairs significant at BOTH 1m and 6m in the
#          SAME direction ("Both consistent" bar in direction pattern charts).
#          Visualize trends across categories, PFAS predictors, and individual
#          pathways to identify storytelling candidates.
#
# DEPENDS ON: Objects already in environment from Script 8 visualization block:
#   pwy_cont_1m_1m, pwy_cont_1m_6m, pwy_bin_1m_1m, pwy_bin_1m_6m
#   cat_cols, col_pos, col_neg, pathway_barchart_theme()
# ==============================================================================


# ------------------------------------------------------------------------------
# HELPER: Extract consistent pairs from two timepoint data frames
# Returns the 1m rows with 6m beta/direction appended for side-by-side view
# ------------------------------------------------------------------------------

extract_consistent <- function(sig_1m, sig_6m, model_label,
                               beta_col_1m = "plot_beta",
                               beta_col_6m = "plot_beta") {
  
  keys_1m <- sig_1m %>% tidyr::unite("key", predictor, Pathway, remove = FALSE)
  keys_6m <- sig_6m %>% tidyr::unite("key", predictor, Pathway, remove = FALSE)
  
  shared_keys <- intersect(keys_1m$key, keys_6m$key)
  
  df_1m <- keys_1m %>%
    filter(key %in% shared_keys) %>%
    dplyr::select(key, predictor, Pathway, Category,
                  beta_1m  = all_of(beta_col_1m),
                  dir_1m   = Direction,
                  p_1m     = p_val,
                  fdr_1m   = FDR)
  
  df_6m <- keys_6m %>%
    filter(key %in% shared_keys) %>%
    dplyr::select(key,
                  beta_6m  = all_of(beta_col_6m),
                  dir_6m   = Direction,
                  p_6m     = p_val,
                  fdr_6m   = FDR)
  
  joined <- inner_join(df_1m, df_6m, by = "key") %>%
    filter(dir_1m == dir_6m) %>%   # keep only CONSISTENT (same direction)
    mutate(
      Model     = model_label,
      Direction = dir_1m,          # single direction label
      # Shorten pathway name for plotting
      Pathway_short = str_trunc(
        sub("^[A-Z0-9-]+ ", "", Pathway),   # strip MetaCyc ID prefix
        50, side = "right"
      )
    ) %>%
    dplyr::select(-key)
  
  cat(sprintf("  [%s] Consistent pairs: %d (%d negative, %d positive)\n",
              model_label,
              nrow(joined),
              sum(joined$Direction == "Negative"),
              sum(joined$Direction == "Positive")))
  
  return(joined)
}


# ------------------------------------------------------------------------------
# Extract consistent pairs
# NOTE: binary PFAS data frames use raw "betas" (no IQR scaling)
# ------------------------------------------------------------------------------

consist_cont <- extract_consistent(
  pwy_cont_1m_1m, pwy_cont_1m_6m,
  model_label = "Continuous PFAS",
  beta_col_1m = "plot_beta", beta_col_6m = "plot_beta"
)

consist_bin <- extract_consistent(
  pwy_bin_1m_1m, pwy_bin_1m_6m,
  model_label = "Binary PFAS",
  beta_col_1m = "betas", beta_col_6m = "betas"
)

# Save for reference
write.csv(consist_cont,
          here::here("out_files", "consistent_pathways_continuous.csv"),
          row.names = FALSE)
write.csv(consist_bin,
          here::here("out_files", "consistent_pathways_binary.csv"),
          row.names = FALSE)

cat("\nSaved consistent pathway tables.\n")

# ==============================================================================
# PLOT 3: Dumbbell / connected dot plot for top consistent pathways
# Shows beta at 1m (circle) vs 6m (square) for each pathway × predictor
# Faceted by direction; limited to pathways significant for ≥2 predictors
# ==============================================================================

make_dumbbell <- function(consist_df, model_label,
                          min_predictors = 2,
                          max_pathways   = 20) {
  
  pwy_counts <- consist_df %>%
    count(Pathway, Direction, name = "n_pred") %>%
    filter(n_pred >= min_predictors)
  
  plot_df <- consist_df %>%
    semi_join(pwy_counts, by = c("Pathway", "Direction")) %>%
    group_by(Direction) %>%
    mutate(pwy_rank = dense_rank(Pathway)) %>%
    filter(pwy_rank <= max_pathways) %>%
    ungroup() %>%
    pivot_longer(cols = c(beta_1m, beta_6m),
                 names_to = "Timepoint",
                 values_to = "Beta") %>%
    mutate(
      Timepoint = recode(Timepoint, beta_1m = "1 month", beta_6m = "6 months"),
      Timepoint = factor(Timepoint, levels = c("1 month", "6 months")),
      # Order facets: Negative on top, Positive below
      Direction = factor(Direction, levels = c("Negative", "Positive"))
    )
  
  # Manual shortening of specific pathway labels
  pathway_labels <- c(
    "NAD salvage pathway V (PNC V cycle)"          = "NAD salvage pathway",
    "L-N&delta;-acetylornithine biosynthesis"      = "acetylornithine biosynthesis",
    "flavin biosynthesis I (bacteria and plants)"  = "flavin biosynthesis",
    "colanic acid building blocks biosynthesis"    = "colanic acid biosynthesis",
    "acetylene degradation"                        = "acetylene degradation"
  )
  
  plot_df <- plot_df %>%
    mutate(Pathway_short = dplyr::recode(Pathway_short, !!!pathway_labels))
  
  if (nrow(plot_df) == 0) {
    cat(sprintf("  [%s] No pathways meet min_predictors >= %d threshold\n",
                model_label, min_predictors))
    return(invisible(NULL))
  }
  
  p <- ggplot(plot_df,
              aes(x = Beta, y = reorder(Pathway_short, Beta),
                  color = predictor, shape = Timepoint,
                  group = interaction(Pathway_short, predictor))) +
    geom_line(aes(group = interaction(Pathway_short, predictor)),
              color = "grey50", linewidth = 1.2) +
    geom_point(size = 5.5, stroke = 1.2) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey30", linewidth = 0.8) +
    scale_shape_manual(values = c("1 month" = 16, "6 months" = 15)) +
    scale_color_brewer(palette = "Dark2") +
    # ncol = 1 stacks Negative above Positive
    # scales = "free" keeps each panel's x and y independent
    # space = "free_y" makes each panel height proportional to its row count
    facet_wrap(~ Direction, ncol = 1, scales = "free", space = "free_y") +
    labs(
      title    = NULL,
      subtitle = NULL,
      x = "IQR-scaled β (change in CLR-transformed pathway abundance\nper 1-IQR increase in log\u2082-transformed PFAS concentration",
      y        = NULL,
      color    = "PFAS predictor",
      shape    = "Timepoint"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      strip.text         = element_text(size = 15, face = "bold", color = "white"),
      strip.background   = element_rect(fill = "black", color = NA),
      axis.text.y        = element_text(size = 13, color = "black"),
      axis.text.x        = element_text(size = 13, color = "black"),
      axis.title.x       = element_text(size = 13, face = "bold",
                                        color = "black", margin = margin(t = 10)),
      axis.line.x        = element_line(color = "black", linewidth = 0.5),
      axis.line.y        = element_line(color = "black", linewidth = 0.5),
      axis.ticks         = element_line(color = "black"),
      legend.position    = "right",
      legend.title       = element_text(size = 12, face = "bold"),
      legend.text        = element_text(size = 11),
      legend.key.size    = unit(0.5, "cm"),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white", color = NA),
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.spacing.y    = unit(0.8, "lines"),   # gap between Negative and Positive panels
      plot.title         = element_text(size = 15, face = "bold"),
      plot.subtitle      = element_text(size = 12, color = "grey30"),
      plot.margin        = margin(10, 15, 10, 10)
    ) +
    scale_y_discrete(expand = expansion(add = 0.4))  # tighten row padding
  
  return(p)
}

# --- Calling code with dynamic height -----------------------------------------

p_dumbbell_cont <- make_dumbbell(consist_cont, "Continuous PFAS",
                                 min_predictors = 2, max_pathways = 25)
if (!is.null(p_dumbbell_cont)) {
  n_pwy <- consist_cont %>%
    count(Pathway, Direction) %>%
    filter(n >= 2) %>%
    nrow()
  plot_height <- max(4, n_pwy * 0.65)  # 0.55 inches per pathway row; tune as needed
  print(p_dumbbell_cont)
  ggsave(here::here("out_figures", "consistent_dumbbell_continuous.pdf"),
         plot = p_dumbbell_cont, width = 10, height = plot_height, device = cairo_pdf)
  cat("  Saved: consistent_dumbbell_continuous.pdf\n")
}

p_dumbbell_bin <- make_dumbbell(consist_bin, "Binary PFAS",
                                min_predictors = 2, max_pathways = 25)
if (!is.null(p_dumbbell_bin)) {
  n_pwy <- consist_bin %>%
    count(Pathway, Direction) %>%
    filter(n >= 2) %>%
    nrow()
  plot_height <- max(4, n_pwy * 0.65)
  print(p_dumbbell_bin)
  ggsave(here::here("out_figures", "consistent_dumbbell_binary.pdf"),
         plot = p_dumbbell_bin, width = 10, height = plot_height, device = cairo_pdf)
  cat("  Saved: consistent_dumbbell_binary.pdf\n")
}









make_dumbbell <- function(consist_df, model_label,
                          min_predictors = 2,
                          max_pathways   = 20) {
  
  pwy_counts <- consist_df %>%
    count(Pathway, Direction, name = "n_pred") %>%
    filter(n_pred >= min_predictors)
  
  plot_df <- consist_df %>%
    semi_join(pwy_counts, by = c("Pathway", "Direction")) %>%
    group_by(Direction) %>%
    mutate(pwy_rank = dense_rank(Pathway)) %>%
    filter(pwy_rank <= max_pathways) %>%
    ungroup() %>%
    pivot_longer(cols = c(beta_1m, beta_6m),
                 names_to = "Timepoint",
                 values_to = "Beta") %>%
    mutate(
      Timepoint = recode(Timepoint, beta_1m = "1 month", beta_6m = "6 months"),
      Timepoint = factor(Timepoint, levels = c("1 month", "6 months")),
      Direction = factor(Direction, levels = c("Negative", "Positive"))
    )
  
  # Manual shortening of specific pathway labels
  pathway_labels <- c(
    "NAD salvage pathway V (PNC V cycle)"          = "NAD salvage pathway",
    "L-N&delta;-acetylornithine biosynthesis"      = "acetylornithine biosynthesis",
    "flavin biosynthesis I (bacteria and plants)"  = "flavin biosynthesis",
    "colanic acid building blocks biosynthesis"    = "colanic acid biosynthesis",
    "acetylene degradation"                        = "acetylene degradation"
  )
  
  plot_df <- plot_df %>%
    mutate(Pathway_short = dplyr::recode(Pathway_short, !!!pathway_labels))
  
  if (nrow(plot_df) == 0) {
    cat(sprintf("  [%s] No pathways meet min_predictors >= %d threshold\n",
                model_label, min_predictors))
    return(invisible(NULL))
  }
  
  # Order pathways: negative ones on top (sorted by beta), positive below
  # This keeps the two groups visually separated around the 0 line
  neg_order <- plot_df %>%
    filter(Direction == "Negative") %>%
    group_by(Pathway_short) %>%
    summarise(mean_beta = mean(Beta), .groups = "drop") %>%
    arrange(mean_beta) %>%          # most negative at top
    pull(Pathway_short)
  
  pos_order <- plot_df %>%
    filter(Direction == "Positive") %>%
    group_by(Pathway_short) %>%
    summarise(mean_beta = mean(Beta), .groups = "drop") %>%
    arrange(mean_beta) %>%          # least positive just above 0
    pull(Pathway_short)
  
  pwy_order <- unique(c(neg_order, pos_order))
  
  plot_df <- plot_df %>%
    mutate(Pathway_short = factor(Pathway_short, levels = pwy_order))
  
  p <- ggplot(plot_df,
              aes(x = Beta, y = Pathway_short,
                  color = predictor, shape = Timepoint,
                  group = interaction(Pathway_short, predictor))) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey30", linewidth = 0.8) +
    geom_line(aes(group = interaction(Pathway_short, predictor)),
              color = "grey50", linewidth = 1.2) +
    geom_point(size = 5.5, stroke = 1.2) +
    scale_shape_manual(values = c("1 month" = 16, "6 months" = 15)) +
    scale_color_brewer(palette = "Dark2") +
    labs(
      title    = NULL,
      subtitle = NULL,
      x        = "IQR-scaled \u03b2 (change in CLR-transformed pathway abundance\nper 1-IQR increase in log\u2082-transformed PFAS concentration)",
      y        = NULL,
      color    = "PFAS predictor",
      shape    = "Timepoint"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.y        = element_text(size = 13, color = "black"),
      axis.text.x        = element_text(size = 13, color = "black"),
      axis.title.x       = element_text(size = 13, face = "bold",
                                        color = "black", margin = margin(t = 10)),
      axis.line.x        = element_line(color = "black", linewidth = 0.5),
      axis.line.y        = element_line(color = "black", linewidth = 0.5),
      axis.ticks         = element_line(color = "black"),
      legend.position    = "right",
      legend.title       = element_text(size = 12, face = "bold"),
      legend.text        = element_text(size = 11),
      legend.key.size    = unit(0.5, "cm"),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white", color = NA),
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(size = 15, face = "bold"),
      plot.subtitle      = element_text(size = 12, color = "grey30"),
      plot.margin        = margin(10, 15, 10, 10)
    ) +
    scale_y_discrete(expand = expansion(add = 0.4))
  
  return(p)
}

# --- Calling code with dynamic height -----------------------------------------

p_dumbbell_cont <- make_dumbbell(consist_cont, "Continuous PFAS",
                                 min_predictors = 2, max_pathways = 25)
if (!is.null(p_dumbbell_cont)) {
  n_pwy <- consist_cont %>%
    count(Pathway, Direction) %>%
    filter(n >= 2) %>%
    nrow()
  plot_height <- max(4, n_pwy * 0.65)
  print(p_dumbbell_cont)
  ggsave(here::here("out_figures", "consistent_dumbbell_continuous.pdf"),
         plot = p_dumbbell_cont, width = 10, height = plot_height, device = cairo_pdf)
  cat("  Saved: consistent_dumbbell_continuous.pdf\n")
}

p_dumbbell_bin <- make_dumbbell(consist_bin, "Binary PFAS",
                                min_predictors = 2, max_pathways = 25)
if (!is.null(p_dumbbell_bin)) {
  n_pwy <- consist_bin %>%
    count(Pathway, Direction) %>%
    filter(n >= 2) %>%
    nrow()
  plot_height <- max(4, n_pwy * 0.65)
  print(p_dumbbell_bin)
  ggsave(here::here("out_figures", "consistent_dumbbell_binary.pdf"),
         plot = p_dumbbell_bin, width = 10, height = plot_height, device = cairo_pdf)
  cat("  Saved: consistent_dumbbell_binary.pdf\n")
}


# Save as excel in single file (all 4 scenarios) for supplemental
library(readr)
library(openxlsx)
library(here)

# Load the four pathway result files
cont_1m1m <- read_csv(here("out_files", "PATHWAY_CLR_continuous_1m_1m.csv"))
cont_1m6m <- read_csv(here("out_files", "PATHWAY_CLR_continuous_1m_6m.csv"))
bin_1m1m  <- read_csv(here("out_files", "PATHWAY_CLR_binary_1m_1m.csv"))
bin_1m6m  <- read_csv(here("out_files", "PATHWAY_CLR_binary_1m_6m.csv"))

# Create workbook
wb <- createWorkbook()

addWorksheet(wb, "PFAS1m_Pathway1m_Continuous")
writeData(wb, "PFAS1m_Pathway1m_Continuous", cont_1m1m)

addWorksheet(wb, "PFAS1m_Pathway6m_Continuous")
writeData(wb, "PFAS1m_Pathway6m_Continuous", cont_1m6m)

addWorksheet(wb, "PFAS1m_Pathway1m_Binary")
writeData(wb, "PFAS1m_Pathway1m_Binary", bin_1m1m)

addWorksheet(wb, "PFAS1m_Pathway6m_Binary")
writeData(wb, "PFAS1m_Pathway6m_Binary", bin_1m6m)

# Save workbook
saveWorkbook(
  wb,
  here("out_files", "Supplemental_Table_T2.xlsx"),
  overwrite = TRUE
)








