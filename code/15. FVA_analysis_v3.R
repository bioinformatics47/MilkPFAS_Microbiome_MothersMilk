# =============================================================================
# PFAS FVA ANALYSIS
#
# Includes:
#   - creation of the four 73-species FVA files
#   - SILVA/GTDB matching against the 73-species FVA list
#   - manual JEB00298 -> Bifidobacterium adolescentis correction
#   - matching audit and missing-species diagnostics
#   - uptake and secretion Wald tests
#   - Supplemental Tables T3 and T4
#   - Bifidobacterium capability verification
#
#  Code review: Daniel Fassler on August 3, 2026
#               Rechecked by Devendra Paudel on August 3, 2026
# =============================================================================

library(lme4)
library(car)
library(data.table)
library(openxlsx)
library(dplyr)
library(stringr)
library(here)

# =============================================================================
# Project paths
# =============================================================================

# setwd() removed for portability: run this script with the working
# directory set to the project root (the folder containing this "code" dir)

# All input files are stored directly in the PFAS working directory.
DIR_JOHANNES <- "insilico_data"
DIR_YUAN     <- "invitro_data"

TAB <- file.path(getwd(), "out_files")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(DIR_JOHANNES)) {
  stop("DIR_JOHANNES not found: ",
       normalizePath(DIR_JOHANNES, mustWork = FALSE))
}

if (!dir.exists(DIR_YUAN)) {
  stop("DIR_YUAN not found: ",
       normalizePath(DIR_YUAN, mustWork = FALSE))
}


# =============================================================================
# Part 1: Create _73 FVA files
# =============================================================================

species_to_remove_from_fva <- c(
  "Blautia coccoides",
  "Gordonibacter urolithinfaciens",
  "Lacticaseibacillus paracasei"
)

make_73_species_fva <- function(input_file) {
  
  dat <- read.csv(
    input_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (!"Species" %in% names(dat)) {
    stop(input_file, " does not contain a column named Species.")
  }
  
  n_before <- length(unique(dat$Species))
  
  dat_73 <- dat %>%
    filter(!(Species %in% species_to_remove_from_fva))
  
  n_after <- length(unique(dat_73$Species))
  
  # Insert _73 before .csv (works regardless of a _76 tag)
  output_file <- sub("\\.csv$", "_73.csv", input_file)
  
  if (output_file == input_file) {
    stop("Input file name does not end with .csv: ", input_file)
  }
  
  write.csv(
    dat_73,
    output_file,
    row.names = FALSE
  )
  
  cat("\nFile:", input_file, "\n")
  cat("Unique species before:", n_before, "\n")
  cat("Unique species after:", n_after, "\n")
  cat("Saved as:", output_file, "\n")
  
  if (n_before - n_after != 3) {
    warning(
      "Expected to remove 3 species (", n_before, " -> ", n_before - 3,
      "), but ended with ", n_after,
      ". Check species_to_remove_from_fva against this file's names."
    )
  }
  
  invisible(dat_73)
}

make_73_species_fva(file.path(DIR_JOHANNES, "FVA_Secretion_Capacity_HMO_BM_1month.csv"))
make_73_species_fva(file.path(DIR_JOHANNES, "FVA_Secretion_Capacity_HMO_Complex.csv"))
make_73_species_fva(file.path(DIR_JOHANNES, "FVA_Uptake_Capacity_HMO_BM_1month.csv"))
make_73_species_fva(file.path(DIR_JOHANNES, "FVA_Uptake_Capacity_HMO_Complex.csv"))


# =============================================================================
# Part 2: Load VMH and metadata
# =============================================================================

VMH_metabolites <- fread(
  file.path(DIR_JOHANNES, "VMH_Metabolites.tsv")
)

meta <- read.csv(
  file.path(DIR_YUAN, "formated_data.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Fix empty column names if present
bad_names <- is.na(names(meta)) | trimws(names(meta)) == ""

if (any(bad_names)) {
  names(meta)[bad_names] <- paste0("unnamed_column_", seq_len(sum(bad_names)))
}

names(meta) <- make.unique(names(meta))

meta <- meta[!(meta$StrainID == "JEB00285" & meta$median_value <= 1), ]
meta <- meta[meta$Drug %in% c("Veh", "PFOS-low", "PFOS-high"), ]


# =============================================================================
# Part 3: Helper functions
# =============================================================================

empty_to_na <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x[x %in% c("", "NA", "NaN", "nan", "NULL", "null")] <- NA_character_
  x
}

clean_species_name <- function(x) {
  x <- empty_to_na(x)
  
  # Example: Blautia_A obeum -> Blautia obeum
  x <- str_replace(
    x,
    "^([A-Za-z]+)_[A-Z]+\\s+",
    "\\1 "
  )
  
  # Example: species_A -> species
  x <- str_replace(
    x,
    "\\s+([a-z-]+)_[A-Z]+$",
    " \\1"
  )
  
  x <- str_squish(x)
  empty_to_na(x)
}

extract_gtdb_species <- function(x) {
  x <- empty_to_na(x)
  
  out <- rep(NA_character_, length(x))
  
  has_species_label <- !is.na(x) & str_detect(x, "s__")
  
  out[has_species_label] <- str_replace(
    x[has_species_label],
    "^.*s__",
    ""
  )
  
  # fallback if column already contains only a species name
  out[!has_species_label & !is.na(x)] <- x[!has_species_label & !is.na(x)]
  
  # remove possible trailing taxonomy
  out <- str_replace(out, ";.*$", "")
  
  clean_species_name(out)
}

clean_vmh_id <- function(x) {
  x <- as.character(x)
  
  # Remove exchange prefix
  x <- gsub("^EX_", "", x)
  
  # Remove extracellular compartment suffixes
  # Examples:
  # isocapr(e)     -> isocapr
  # 12ppd_S(e)     -> 12ppd_S
  # isocapr.e      -> isocapr
  # isocapr.e.     -> isocapr
  # isocapr[e]     -> isocapr
  x <- gsub("\\(e\\)$", "", x)
  x <- gsub("\\[e\\]$", "", x)
  x <- gsub("\\.e\\.?$", "", x)
  
  x
}


# =============================================================================
# Part 4: Read final 73-species FVA list
# =============================================================================

fva_reference <- read.csv(
  file.path(DIR_JOHANNES, "FVA_Secretion_Capacity_HMO_BM_1month_73.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"Species" %in% names(fva_reference)) {
  stop("FVA reference file does not contain a Species column.")
}

fva_species_73 <- fva_reference %>%
  mutate(
    Species = clean_species_name(Species)
  ) %>%
  filter(!is.na(Species)) %>%
  distinct(Species) %>%
  pull(Species)

cat(
  "\nUnique species in FVA reference file: ",
  length(fva_species_73),
  "\n",
  sep = ""
)

if (length(fva_species_73) != 73) {
  warning(
    "Expected 73 species in the FVA reference file, but found ",
    length(fva_species_73)
  )
}


# =============================================================================
# Part 5: Identify GTDB taxonomy column
# =============================================================================

gtdb_candidates <- c(
  "GTDB_08.rs214_taxonomy",
  "GTDB_08-rs214_taxonomy"
)

gtdb_col <- intersect(gtdb_candidates, names(meta))

if (length(gtdb_col) == 0) {
  stop(
    "Could not find GTDB taxonomy column. Expected one of: ",
    paste(gtdb_candidates, collapse = ", ")
  )
}

gtdb_col <- gtdb_col[1]

cat("\nUsing GTDB column: ", gtdb_col, "\n", sep = "")


# =============================================================================
# Part 6: Create final meta$Species merge column
# =============================================================================

meta <- meta %>%
  mutate(
    Species_SILVA_clean = clean_species_name(SILVA_Species),
    Species_GTDB_clean = extract_gtdb_species(.data[[gtdb_col]]),
    
    SILVA_in_FVA = Species_SILVA_clean %in% fva_species_73,
    GTDB_in_FVA = Species_GTDB_clean %in% fva_species_73,
    
    Species_merge_source = case_when(
      StrainID == "JEB00298" ~ "manual_JEB00298_GTDB",
      
      SILVA_in_FVA & GTDB_in_FVA &
        Species_SILVA_clean == Species_GTDB_clean ~ "SILVA_and_GTDB_agree",
      
      SILVA_in_FVA & !GTDB_in_FVA ~ "SILVA_only_in_FVA",
      
      GTDB_in_FVA & !SILVA_in_FVA ~ "GTDB_only_in_FVA",
      
      is.na(Species_GTDB_clean) & SILVA_in_FVA ~ "SILVA_used_GTDB_missing",
      
      is.na(Species_SILVA_clean) & GTDB_in_FVA ~ "GTDB_used_SILVA_missing",
      
      TRUE ~ "not_matched"
    ),
    
    # This is the final merge column used by the Wald code
    Species = case_when(
      StrainID == "JEB00298" ~ "Bifidobacterium adolescentis",
      
      SILVA_in_FVA & GTDB_in_FVA &
        Species_SILVA_clean == Species_GTDB_clean ~ Species_SILVA_clean,
      
      SILVA_in_FVA & !GTDB_in_FVA ~ Species_SILVA_clean,
      
      GTDB_in_FVA & !SILVA_in_FVA ~ Species_GTDB_clean,
      
      is.na(Species_GTDB_clean) & SILVA_in_FVA ~ Species_SILVA_clean,
      
      is.na(Species_SILVA_clean) & GTDB_in_FVA ~ Species_GTDB_clean,
      
      TRUE ~ NA_character_
    )
  )

# Keep only rows that can be merged to the final 73-species FVA files
meta <- meta %>%
  filter(!is.na(Species))

meta_adapted_file <- file.path(
  TAB,
  "formated_data_adapted.csv"
)

write.csv(
  meta,
  file = meta_adapted_file,
  row.names = FALSE,
  na = ""
)


# =============================================================================
# Part 7: Audit matching
# =============================================================================

species_matching_summary <- meta %>%
  count(
    Species_merge_source,
    name = "n_rows"
  ) %>%
  arrange(Species_merge_source)

species_check <- meta %>%
  distinct(
    Species,
    Species_merge_source
  ) %>%
  arrange(Species)

meta_species_not_in_fva <- setdiff(
  unique(meta$Species),
  fva_species_73
)

fva_species_not_in_meta <- setdiff(
  fva_species_73,
  unique(meta$Species)
)

cat("\nSpecies matching summary:\n")
print(species_matching_summary)

cat(
  "\nUnique final Species in meta: ",
  length(unique(meta$Species)),
  "\n",
  sep = ""
)

cat(
  "Species in meta but not in FVA: ",
  length(meta_species_not_in_fva),
  "\n",
  sep = ""
)

cat(
  "Species in FVA but not in meta: ",
  length(fva_species_not_in_meta),
  "\n",
  sep = ""
)


cat("\nSpecies present in metadata but absent from the FVA list:\n")
if (length(meta_species_not_in_fva) == 0) {
  cat("None\n")
} else {
  print(meta_species_not_in_fva)
}

cat("\nSpecies present in the FVA list but absent from metadata:\n")
if (length(fva_species_not_in_meta) == 0) {
  cat("None\n")
} else {
  print(fva_species_not_in_meta)
}

if (length(unique(meta$Species)) != 73) {
  warning(
    "Expected 73 unique species in meta$Species, but found ",
    length(unique(meta$Species)),
    ". Check species_check and fva_species_not_in_meta."
  )
}


if (length(fva_species_not_in_meta) > 0) {
  cat("\nMetadata rows potentially related to the missing FVA species:\n")

  possible_missing_matches <- meta %>%
    filter(
      Species_SILVA_clean %in% fva_species_not_in_meta |
      Species_GTDB_clean %in% fva_species_not_in_meta
    ) %>%
    select(
      StrainID,
      Drug,
      SILVA_Species,
      Species_SILVA_clean,
      all_of(gtdb_col),
      Species_GTDB_clean,
      Species_merge_source,
      Species
    ) %>%
    distinct()

  if (nrow(possible_missing_matches) == 0) {
    cat("No retained metadata rows matched the missing FVA species name.\n")
  } else {
    print(possible_missing_matches)
  }
}

write.xlsx(
  list(
    species_matching_summary = species_matching_summary,
    species_check = species_check,
    meta_species_not_in_fva = data.frame(Species = meta_species_not_in_fva),
    fva_species_not_in_meta = data.frame(Species = fva_species_not_in_meta)
  ),
  file = file.path(TAB, "species_merge_column_audit_73.xlsx"),
  overwrite = TRUE
)


# =============================================================================
# Part 8: Wald test
# =============================================================================

run_wald_test <- function(file,
                          flip_sign = FALSE,
                          output_file = NULL) {
  
  dat <- read.csv(
    file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (!"Species" %in% names(dat)) {
    stop("The FVA input file must contain a column named Species.")
  }
  
  dat$Species <- clean_species_name(dat$Species)
  
  binary_vars_all <- setdiff(names(dat), "Species")
  
  # Convert all FVA columns to numeric
  dat[, binary_vars_all] <- lapply(
    dat[, binary_vars_all, drop = FALSE],
    function(x) suppressWarnings(as.numeric(as.character(x)))
  )
  
  # Convert to binary capability
  if (flip_sign) {
    dat[, binary_vars_all] <- lapply(
      dat[, binary_vars_all, drop = FALSE],
      function(x) x * (-1)
    )
  }
  
  dat[, binary_vars_all] <- lapply(
    dat[, binary_vars_all, drop = FALSE],
    function(x) as.integer(x > 0)
  )
  
  merged <- merge(meta, dat, by = "Species")
  
  if (nrow(merged) == 0) {
    stop("No rows after merge. Check meta$Species and dat$Species.")
  }
  
  merged$Drug <- factor(merged$Drug)
  
  binary_vars <- setdiff(names(dat), "Species")
  
  # retain metabolites with prevalence 10-90%
  binary_vars <- binary_vars[
    sapply(dat[, binary_vars, drop = FALSE], function(x) {
      p1 <- mean(x == 1, na.rm = TRUE)
      p1 >= 0.1 & p1 <= 0.9
    })
  ]
  
  if (length(binary_vars) == 0) {
    stop("No metabolite variables passed the 10-90% prevalence filter.")
  }
  
  merged[, binary_vars] <- lapply(
    merged[, binary_vars, drop = FALSE],
    factor,
    levels = c(0, 1)
  )
  
  merged <- merged[complete.cases(merged[, binary_vars]), ]
  
  wald_results <- data.frame(
    Variable = character(),
    Chisq = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (var in binary_vars) {
    
    fit <- lmer(
      as.formula(
        paste0("log(median_value) ~ `", var, "` * Drug + (1|Species)")
      ),
      data = merged,
      REML = FALSE
    )
    
    coef_names <- names(fixef(fit))
    
    interaction_terms <- coef_names[
      grepl(":", coef_names, fixed = TRUE) &
        grepl("Drug", coef_names, fixed = TRUE)
    ]
    
    if (length(interaction_terms) > 0) {
      
      test_result <- linearHypothesis(
        fit,
        paste(interaction_terms, "= 0"),
        test = "Chisq"
      )
      
      wald_results <- rbind(
        wald_results,
        data.frame(
          Variable = var,
          Chisq = test_result$Chisq[2],
          p_value = test_result$`Pr(>Chisq)`[2],
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  # ---------------------------------------------------------------------------
  # Correct VMH ID cleaning
  # ---------------------------------------------------------------------------
  # Examples:
  # EX_isocapr(e)  -> isocapr
  # EX_12ppd_S(e)  -> 12ppd_S
  # EX_isocapr.e   -> isocapr
  # ---------------------------------------------------------------------------
  
  wald_results$VMH_ID <- clean_vmh_id(wald_results$Variable)
  
  wald_results$fullName <-
    VMH_metabolites$fullName[
      match(
        wald_results$VMH_ID,
        VMH_metabolites$abbreviation
      )
    ]
  
  wald_results <- wald_results[
    ,
    c("fullName", "VMH_ID", "Chisq", "p_value")
  ]
  
  if (!is.null(output_file)) {
    write.xlsx(
      wald_results,
      output_file,
      rowNames = FALSE,
      overwrite = TRUE
    )
  }
  
  return(wald_results)
}


# =============================================================================
# Part 9: Run Wald tests on the _73 files
# =============================================================================

# Run all four (both media x uptake/secretion). Keep results in memory;
# individual files are NOT written -- we consolidate into T3 and T4 below.

uptake_bm1     <- run_wald_test(
  file = file.path(DIR_JOHANNES, "FVA_Uptake_Capacity_HMO_BM_1month_73.csv"),
  flip_sign = TRUE
)
uptake_complex <- run_wald_test(
  file = file.path(DIR_JOHANNES, "FVA_Uptake_Capacity_HMO_Complex_73.csv"),
  flip_sign = TRUE
)
secretion_bm1     <- run_wald_test(
  file = file.path(DIR_JOHANNES, "FVA_Secretion_Capacity_HMO_BM_1month_73.csv"),
  flip_sign = FALSE
)
secretion_complex <- run_wald_test(
  file = file.path(DIR_JOHANNES, "FVA_Secretion_Capacity_HMO_Complex_73.csv"),
  flip_sign = FALSE
)


# =============================================================================
# Part 10: Build Supplemental Tables T3 (uptake) and T4 (secretion)
#          Each combines BOTH media (1-month human milk + complex medium)
# =============================================================================

combine_media <- function(bm1, complex) {
  out <- merge(
    bm1, complex,
    by = c("fullName", "VMH_ID"),
    all = TRUE,
    suffixes = c("_BM1month", "_Complex")
  )
  # rename p_value_* -> p_* for readability
  names(out)[names(out) == "p_value_BM1month"] <- "p_BM1month"
  names(out)[names(out) == "p_value_Complex"]  <- "p_Complex"
  out <- out[, c("fullName", "VMH_ID",
                 "Chisq_BM1month", "p_BM1month",
                 "Chisq_Complex",  "p_Complex")]
  # sort by human-milk p-value (primary medium), NA p's last
  out[order(out$p_BM1month, na.last = TRUE), ]
}

T3_uptake    <- combine_media(uptake_bm1,    uptake_complex)
T4_secretion <- combine_media(secretion_bm1, secretion_complex)

write.xlsx(
  T3_uptake,
  file.path(TAB, "SupplementalTable_T3_Uptake.xlsx"),
  rowNames = FALSE,
  overwrite = TRUE
)
write.xlsx(
  T4_secretion,
  file.path(TAB, "SupplementalTable_T4_Secretion.xlsx"),
  rowNames = FALSE,
  overwrite = TRUE
)

cat("\nWrote SupplementalTable_T3_Uptake.xlsx (", nrow(T3_uptake), " metabolites)\n", sep = "")
cat("Wrote SupplementalTable_T4_Secretion.xlsx (", nrow(T4_secretion), " metabolites)\n", sep = "")


# =============================================================================
# Part 11: Verify per-species Bifidobacterium capability claims (Results 2.9)
#   "Bifidobacterium species possessed ... absent D-alanine uptake, lactose
#    uptake capability, and absent cadaverine secretion"
#   Reads the RAW _73 FVA values (not binarized).
# =============================================================================

# Resolve an exchange-reaction column whether it ends in (e), [e], or nothing
get_ex <- function(df, id) {
  cands <- c(paste0("EX_", id, "(e)"),
             paste0("EX_", id, "[e]"),
             paste0("EX_", id))
  hit <- cands[cands %in% names(df)][1]
  if (is.na(hit)) return(rep(NA_real_, nrow(df)))
  df[[hit]]
}

bif_up  <- read.csv(
  file.path(DIR_JOHANNES, "FVA_Uptake_Capacity_HMO_BM_1month_73.csv"),
  check.names = FALSE
)
bif_sec <- read.csv(
  file.path(DIR_JOHANNES, "FVA_Secretion_Capacity_HMO_BM_1month_73.csv"),
  check.names = FALSE
)

bif_up  <- bif_up[grepl("Bifidobacterium",  bif_up$Species), ]
bif_sec <- bif_sec[grepl("Bifidobacterium", bif_sec$Species), ]

cat("\n=====================================================================\n")
cat("VERIFY: Bifidobacterium per-species capabilities \n")
cat("=====================================================================\n")

cat("\nD-alanine uptake     (claim: ABSENT  -> value ~ 0):\n")
print(data.frame(Species = bif_up$Species,
                 ala_D   = get_ex(bif_up, "ala_D")),
      row.names = FALSE)

cat("\nLactose uptake       (claim: PRESENT -> value < 0):\n")
print(data.frame(Species = bif_up$Species,
                 lcts    = get_ex(bif_up, "lcts")),
      row.names = FALSE)

cat("\nCadaverine secretion (claim: ABSENT  -> value ~ 0):\n")
print(data.frame(Species    = bif_sec$Species,
                 cadaverine = get_ex(bif_sec, "15dap")),
      row.names = FALSE)
