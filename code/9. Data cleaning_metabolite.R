# TITLE:   1. Data Cleaning.R
#
# PURPOSE: Read and clean breastmilk PFAS, Fecal metabolome and metadata
#
# DATE:    May 2025
#
# CODE REVIEW:
# Ellie Holzhausen (EAH) on April 27, 2026
# by Haonan Li (HL) on May 25, 2026
#
# ANNOTATION DEDUPLICATION APPROACH:
# Follows Li et al. (in preparation) supplementary methods exactly:
#
# Step 1 — Exact mz/RT deduplication on annotation output rows (BEFORE transposing):
#   Sub-case A: Same compound name appearing twice (redundant library entry)
#               → Keep first occurrence, remove duplicate row silently
#   Sub-case B: Different compound names sharing identical detected mz AND RT (co-eluters)
#               → Assign combined slash label to first row, remove subsequent rows
#               → Preserves transparency per MSI Level 1 annotation standards
#
#
# set up -----------------------------------------------------------------------

rm(list = ls())
options(scipen = 100)

library(tidyverse)
library(dplyr); library(tidyr); library(reshape)
library(purrr); library(stringr)
library(corrplot); library(lubridate)
library(MASS); library(pscl); library(naniar)
library(tableone); library(readxl); library(here)
library(tibble);  library(xMSanalyzer)
library(tidyr)


# Load meta trim data ----------------------------------------------------------
meta_trim <- read.csv(here::here("out_files/meta_trim.csv"))

# Load fecal metabolite data ---------------------------------------------------
Featuretable_c18 <- readRDS(here::here("input", "c18_rf_QRILC_badremove.rds"))
Featuretable_hilic <- readRDS(here::here("input", "hilic_rf_QRILC_badremove.rds"))

Confirmed_c18 <- read.csv(here::here("input", "c18_total_confirmed.csv"))
Confirmed_hilic <- read.csv(here::here("input", "hilic_total_confirmed.csv"))

# Shortnames_c18 <- read_excel(here::here("input", "c18_labs2.xlsx"))
# Shortnames_hilic <- read_excel(here::here("input", "hilic_labs2.xlsx"))

c18_preprocessed <- load(here::here("input", "c18_preprocessed_metabo.RData"))
hilic_preprocessed <- load(here::here("input", "hilic_preprocessed_metabo.RData"))


# Clean reference library files ------------------------------------------------
if ("X" %in% colnames(Confirmed_c18))   Confirmed_c18   <- dplyr::select(Confirmed_c18,   -"X")
if ("X" %in% colnames(Confirmed_hilic)) Confirmed_hilic <- dplyr::select(Confirmed_hilic, -"X")

Confirmed_c18$rt   <- as.numeric(Confirmed_c18$rt)
Confirmed_hilic$rt <- as.numeric(Confirmed_hilic$rt)

Confirmed_c18   <- Confirmed_c18[complete.cases(Confirmed_c18), ]
Confirmed_hilic <- Confirmed_hilic[complete.cases(Confirmed_hilic), ]

sum(is.na(Confirmed_c18))   # 0
sum(is.na(Confirmed_hilic)) # 0


# Clean feature files ----------------------------------------------------------
Featuretable_c18   <- data.frame(t(Featuretable_c18),   check.names = FALSE)
Featuretable_hilic <- data.frame(t(Featuretable_hilic), check.names = FALSE)

Featuretable_c18 <- cbind(
  Featuretable_c18[intersect(rownames(Featuretable_c18), rownames(c18.feat)), ],
  c18.feat[intersect(rownames(Featuretable_c18), rownames(c18.feat)), c("mz", "time")]
)

Featuretable_hilic <- cbind(
  Featuretable_hilic[intersect(rownames(Featuretable_hilic), rownames(hilic.feat)), ],
  hilic.feat[intersect(rownames(Featuretable_hilic), rownames(hilic.feat)), c("mz", "time")]
)



# ANNOTATION — C18-------------------------------------

Sig <- dplyr::select(Featuretable_c18, mz, time)
Ref <- dplyr::select(Confirmed_c18, c("mz", "rt"))
rownames(Ref) <- NULL

masteroverlap <- find.Overlapping.mzs(Sig, Ref, mz.thresh = 10, time.thresh = 50,
                                      alignment.tool = 'apLCMS')

Ref.sig   <- slice(Featuretable_c18, masteroverlap$index.A)
Ref.match <- slice(Confirmed_c18,    masteroverlap$index.B)
c18_ann   <- cbind(Ref.match, Ref.sig)

# Save C18 annotation objects before HILIC overwrites Ref.match / Ref.sig
c18_Ref_match <- Ref.match
c18_Ref_sig   <- Ref.sig


# =============================================================================
# ANNOTATION — HILIC
# =============================================================================

Sig <- dplyr::select(Featuretable_hilic, mz, time)
Ref <- dplyr::select(Confirmed_hilic, c("mz", "rt"))
rownames(Ref) <- NULL

masteroverlap <- find.Overlapping.mzs(Sig, Ref, mz.thresh = 10, time.thresh = 50,
                                      alignment.tool = 'apLCMS')

Ref.sig   <- slice(Featuretable_hilic, masteroverlap$index.A)
Ref.match <- slice(Confirmed_hilic,    masteroverlap$index.B)
hilic_ann <- cbind(Ref.match, Ref.sig)

# Save HILIC annotation objects
hilic_Ref_match <- Ref.match
hilic_Ref_sig   <- Ref.sig


# =============================================================================
# DIAGNOSTIC: Print all duplicate groups before any deduplication
# Follows Li et al. supplementary methods classification framework
# =============================================================================

c18_full <- data.frame(
  CNAME       = c18_Ref_match$CNAME,
  detected_mz = c18_Ref_sig$mz,
  detected_rt = c18_Ref_sig$time,
  ref_mz      = c18_Ref_match$mz,
  ref_rt      = c18_Ref_match$rt,
  stringsAsFactors = FALSE
)

dup_keys_c18   <- paste(c18_full$detected_mz, c18_full$detected_rt)
dup_mask_c18   <- dup_keys_c18 %in% dup_keys_c18[duplicated(dup_keys_c18)]
dup_c18        <- c18_full[dup_mask_c18, ]
dup_c18        <- dup_c18[order(dup_c18$detected_mz, dup_c18$detected_rt), ]
groups_c18_all <- unique(paste(dup_c18$detected_mz, dup_c18$detected_rt))

cat("Total duplicate groups in C18:", length(groups_c18_all), "\n")
cat("(These will be resolved by Step 1 deduplication below)\n\n")

for (i in seq_along(groups_c18_all)) {
  grp <- dup_c18[paste(dup_c18$detected_mz, dup_c18$detected_rt) == groups_c18_all[i], ]
  cat("--- C18 Group", i,
      "| detected mz:", unique(grp$detected_mz),
      "| detected RT:", unique(grp$detected_rt), "---\n")
  for (j in seq_len(nrow(grp))) {
    cat("  [", j, "]", grp$CNAME[j],
        "| ref mz:", round(grp$ref_mz[j], 5),
        "| ref RT:", round(grp$ref_rt[j], 1), "\n")
  }
  if (length(unique(grp$CNAME)) == 1) {
    cat("  >> SAME NAME — redundant library entry\n\n")
  } else {
    cat("  >> DIFFERENT NAMES — co-eluters or naming variants\n\n")
  }
}


hilic_full <- data.frame(
  CNAME       = hilic_Ref_match$CNAME,
  detected_mz = hilic_Ref_sig$mz,
  detected_rt = hilic_Ref_sig$time,
  ref_mz      = hilic_Ref_match$mz,
  ref_rt      = hilic_Ref_match$rt,
  stringsAsFactors = FALSE
)

dup_keys_hilic   <- paste(hilic_full$detected_mz, hilic_full$detected_rt)
dup_mask_hilic_d <- dup_keys_hilic %in% dup_keys_hilic[duplicated(dup_keys_hilic)]
dup_hilic        <- hilic_full[dup_mask_hilic_d, ]
dup_hilic        <- dup_hilic[order(dup_hilic$detected_mz, dup_hilic$detected_rt), ]
groups_hilic_all <- unique(paste(dup_hilic$detected_mz, dup_hilic$detected_rt))

cat("Total duplicate groups in HILIC:", length(groups_hilic_all), "\n\n")

for (i in seq_along(groups_hilic_all)) {
  grp <- dup_hilic[paste(dup_hilic$detected_mz, dup_hilic$detected_rt) == groups_hilic_all[i], ]
  cat("--- HILIC Group", i,
      "| detected mz:", unique(grp$detected_mz),
      "| detected RT:", unique(grp$detected_rt), "---\n")
  for (j in seq_len(nrow(grp))) {
    cat("  [", j, "]", grp$CNAME[j],
        "| ref mz:", round(grp$ref_mz[j], 5),
        "| ref RT:", round(grp$ref_rt[j], 1), "\n")
  }
  if (length(unique(grp$CNAME)) == 1) {
    cat("  >> SAME NAME — redundant library entry\n\n")
  } else {
    cat("  >> DIFFERENT NAMES — co-eluters or naming variants\n\n")
  }
}


# =============================================================================
# STEP 1 — EXACT mz/RT DEDUPLICATION — C18
# Removes duplicate rows where same detected mz AND RT appear more than once.
# After combined labels applied above, first occurrence of each detected peak
# carries the correct name (single compound name or combined slash label).
# =============================================================================

c18_feat_info <- data.frame(
  CNAME       = c18_Ref_match$CNAME,
  detected_mz = c18_Ref_sig$mz,
  detected_rt = c18_Ref_sig$time,
  stringsAsFactors = FALSE
)

cat("\n=== STEP 1: EXACT mz/RT DEDUPLICATION — C18 ===\n")
cat("C18 annotation rows before dedup:", nrow(c18_ann), "\n")

dup_mzrt_c18 <- duplicated(
  data.frame(mz = c18_Ref_sig$mz, time = c18_Ref_sig$time)
)

cat("C18 rows removed (duplicate detected peaks):", sum(dup_mzrt_c18), "\n")
cat("Removed (second+ occurrences of each detected peak):\n")
print(c18_Ref_match$CNAME[dup_mzrt_c18])

c18_ann <- c18_ann[!dup_mzrt_c18, ]
cat("C18 annotation rows after Step 1 dedup:", nrow(c18_ann), "\n\n")


# =============================================================================
# STEP 1 — EXACT mz/RT DEDUPLICATION — HILIC
# =============================================================================

hilic_feat_info <- data.frame(
  CNAME       = hilic_Ref_match$CNAME,
  detected_mz = hilic_Ref_sig$mz,
  detected_rt = hilic_Ref_sig$time,
  stringsAsFactors = FALSE
)

cat("\n=== STEP 1: EXACT mz/RT DEDUPLICATION — HILIC ===\n")
cat("HILIC annotation rows before dedup:", nrow(hilic_ann), "\n")

dup_mzrt_hilic <- duplicated(
  data.frame(mz = hilic_Ref_sig$mz, time = hilic_Ref_sig$time)
)

cat("HILIC rows removed (duplicate detected peaks):", sum(dup_mzrt_hilic), "\n")
cat("Removed (second+ occurrences of each detected peak):\n")
print(hilic_Ref_match$CNAME[dup_mzrt_hilic])

hilic_ann <- hilic_ann[!dup_mzrt_hilic, ]
cat("HILIC annotation rows after Step 1 dedup:", nrow(hilic_ann), "\n\n")


# =============================================================================
# TRANSPOSE C18
# =============================================================================

c18_ann <- as.data.frame(t(c18_ann))

# Exogenous compounds to exclude — C18 ----------------------------------------
# Classification per Li et al. framework:
# PFAS: same as exposure variables — circular association
# CHLOROBENZOATE: chlorinated pesticide degradation product — exogenous co-exposure
# ANILINE-2-SULFONIC ACID: industrial synthetic chemical, not in HMDB as endogenous

exclude_c18 <- c(
  "PFOS", "PFOA", "PFHXS", "PFNA", "PFDA",
  "CHLOROBENZOATE",
  "ANILINE-2-SULFONIC ACID"
)

# Save mz/rt lookup table for pathway enrichment
c18_mz_rt <- data.frame(
  CNAME = as.character(c18_ann["CNAME", ]),
  mz    = as.numeric(c18_ann["mz",    ]),
  rt    = as.numeric(c18_ann["time",  ])
) %>%
  filter(!is.na(mz), !is.na(rt), CNAME != "mz", CNAME != "time") %>%
  filter(!CNAME %in% exclude_c18) %>%
  dplyr::distinct(CNAME, .keep_all = TRUE)

cat("C18 mz/rt lookup — features saved:", nrow(c18_mz_rt), "\n")
write.csv(c18_mz_rt, here::here("out_files", "c18_mz_rt_lookup.csv"), row.names = FALSE)

# Clean up column names and drop metadata rows
colnames(c18_ann) <- c18_ann["CNAME", ]
c18_ann <- c18_ann[-c(1:5), ]

# Set up study ID, dyad_id, timepoint
c18_ann$studyID   <- rownames(c18_ann)
c18_ann$dyad_id   <- sub("X(\\d+)_.*", "\\1", c18_ann$studyID)
c18_ann$timepoint <- sub("X\\d+_(.*)",  "\\1", c18_ann$studyID)

c18_ann <- c18_ann %>%
  subset(timepoint == "01m" | timepoint == "06m") %>%
  mutate(timepoint = ifelse(timepoint == "01m", 1,
                            ifelse(timepoint == "06m", 6, NA))) %>%
  dplyr::select(-"studyID") %>%
  dplyr::select(c(dyad_id, timepoint, everything()))

fecalMetabolites_c18     <- colnames(c18_ann)[3:ncol(c18_ann)]
all_fecalMetabolites_c18 <- fecalMetabolites_c18

# Remove exogenous compounds
fecalMetabolites_c18 <- fecalMetabolites_c18[!fecalMetabolites_c18 %in% exclude_c18]

actually_excluded_c18 <- exclude_c18[exclude_c18 %in% all_fecalMetabolites_c18]
not_detected_c18      <- exclude_c18[!exclude_c18 %in% all_fecalMetabolites_c18]

cat("C18 compounds present and removed:", paste(actually_excluded_c18, collapse = ", "), "\n")
cat("C18 compounds in exclusion list but not in data:", paste(not_detected_c18, collapse = ", "), "\n")
cat("C18 metabolites remaining after exogenous removal:", length(fecalMetabolites_c18), "\n")

# Convert to numeric
c18_ann[, fecalMetabolites_c18] <- lapply(
  c18_ann[, fecalMetabolites_c18], function(x) as.numeric(x)
)

# Missingness check after annotation and dedup
na_check_c18 <- colSums(is.na(c18_ann[, fecalMetabolites_c18]))
cat("C18 metabolites with any NA:", sum(na_check_c18 > 0), "\n")
if (sum(na_check_c18 > 0) > 0) {
  cat("Metabolites with NAs:\n")
  print(na_check_c18[na_check_c18 > 0])
}


# =============================================================================
# STEP 2: INTENSITY-PROFILE DEDUPLICATION — C18
# Safety net: removes any remaining features with identical intensity vectors.
# After Step 1 applied combined labels and removed duplicate peaks,
# Step 2 catches edge cases where different peaks produce identical imputed values.
# Consistent with Holzhausen et al. published pipeline.
# =============================================================================

cat("\n=== STEP 2: INTENSITY-PROFILE DEDUPLICATION — C18 ===\n")
cat("C18 metabolites before Step 2:", length(fecalMetabolites_c18), "\n")

c18_met_only  <- c18_ann[, fecalMetabolites_c18]
t_c18         <- t(c18_met_only)
t_c18_dedup   <- t_c18[!duplicated(t_c18), ]
c18_met_dedup <- t(t_c18_dedup)

dropped_c18 <- setdiff(fecalMetabolites_c18, colnames(c18_met_dedup))
kept_c18    <- colnames(c18_met_dedup)

cat("C18 metabolites after Step 2:", length(kept_c18), "\n")
cat("C18 metabolites dropped by Step 2:", length(dropped_c18), "\n")
if (length(dropped_c18) > 0) {
  cat("Dropped:", paste(dropped_c18, collapse = ", "), "\n")
}
cat("\n")

non_met_cols         <- colnames(c18_ann)[!colnames(c18_ann) %in% fecalMetabolites_c18]
c18_ann              <- cbind(c18_ann[, non_met_cols], c18_met_dedup)
fecalMetabolites_c18 <- kept_c18

# Build C18 display name lookup from Step 1 diagnostic groups
# For each kept metabolite, find all co-eluting names that were removed in Step 1
# This matches Li et al. Loop J approach: combined labels for display, short names in data

c18_display_df <- data.frame(
  CNAME        = kept_c18,
  display_name = kept_c18,
  stringsAsFactors = FALSE
)

# Use the pre-dedup annotation info to find which names were removed per detected peak
# c18_full was built before Step 1 — contains all annotation rows including removed ones
for (kept in kept_c18) {
  # Find the detected mz/RT for this kept metabolite
  kept_row <- c18_full[c18_full$CNAME == kept & 
                         !duplicated(c18_full$CNAME), ]
  if (nrow(kept_row) == 0) next
  
  mz_val <- kept_row$detected_mz[1]
  rt_val <- kept_row$detected_rt[1]
  
  # Find all annotations at this exact detected mz/RT
  all_at_peak <- c18_full[
    abs(c18_full$detected_mz - mz_val) < 1e-8 &
      abs(c18_full$detected_rt - rt_val) < 1e-8, 
    "CNAME"
  ]
  
  # If more than one annotation at this peak — it is a co-eluter group
  unique_names <- unique(all_at_peak)
  if (length(unique_names) > 1) {
    idx <- which(c18_display_df$CNAME == kept)
    c18_display_df$display_name[idx] <- paste(unique_names, collapse = "; ")
  }
}

write.csv(c18_display_df,
          here::here("out_files", "c18_display_names.csv"),
          row.names = FALSE)
cat("C18 display name lookup saved:", nrow(c18_display_df), "metabolites\n")
cat("C18 co-eluter groups in display lookup:",
    sum(grepl(";", c18_display_df$display_name)), "\n")
cat("C18 co-eluter examples:\n")
print(c18_display_df[grepl(";", c18_display_df$display_name), ])

# =============================================================================
# TRANSPOSE HILIC
# =============================================================================

hilic_ann <- as.data.frame(t(hilic_ann))

# Exogenous compounds to exclude — HILIC --------------------------------------
# Phthalates/plasticizers: exogenous industrial compounds
# Pesticides: exogenous agricultural compounds
# Pharmaceuticals: external medication exposures
# Cocaine metabolites: illicit drug use
# Nicotine/cotinine: maternal smoking biomarkers
# Synthetic phenolics group (combined label): industrial/cosmetic compounds
# CAFFEINE: exogenous stimulant
# N,N-DIMETHYL-1,4-PHENYLENEDIAMINE: industrial synthetic amine, not in HMDB as endogenous

exclude_hilic <- c(
  # Phthalates / plasticizers
  "BIS(2-ETHYLHEXYL)PHTHALATE", "DIISOPROPYLPHTHALATE", "PHTHALIC ANHYDRIDE",
  # Pesticides
  "MALATHION", "METRIBUZIN", "PIRIMICARB",
  # Pharmaceuticals
  "CITALOPRAM", "RANITIDINE", "BUPROPION", "BUFURALOL",
  "MIDAZOLAM", "TOLBUTAMIDE", "BENSERAZIDE", "THEOPHYLLINE",
  # Cocaine and metabolites
  "ECGONINE", "HYDROXYBENZOYLECGONINE", "BENZOYLECGONINE", "METHYL ECGONINE",
  # Nicotine / cotinine
  "(S)-NICOTINE", "NICOTINE", "COTININE",
  # Additional exogenous compounds confirmed via HMDB
  "CAFFEINE",
  "N,N-DIMETHYL-1,4-PHENYLENEDIAMINE",
  "MELANIN"
)

# Save mz/rt lookup table for pathway enrichment
hilic_mz_rt <- data.frame(
  CNAME = as.character(hilic_ann["CNAME", ]),
  mz    = as.numeric(hilic_ann["mz",     ]),
  rt    = as.numeric(hilic_ann["time",   ])
) %>%
  filter(!is.na(mz), !is.na(rt), CNAME != "mz", CNAME != "time") %>%
  filter(!CNAME %in% exclude_hilic) %>%
  dplyr::distinct(CNAME, .keep_all = TRUE)

cat("HILIC mz/rt lookup — features saved:", nrow(hilic_mz_rt), "\n")
write.csv(hilic_mz_rt, here::here("out_files", "hilic_mz_rt_lookup.csv"), row.names = FALSE)

# Clean up column names and drop metadata rows
colnames(hilic_ann) <- hilic_ann["CNAME", ]
hilic_ann <- hilic_ann[-c(1:5), ]

# Set up study ID, dyad_id, timepoint
hilic_ann$studyID   <- rownames(hilic_ann)
hilic_ann$dyad_id   <- sub("X(\\d+)_.*", "\\1", hilic_ann$studyID)
hilic_ann$timepoint <- sub("X\\d+_(.*)", "\\1", hilic_ann$studyID)

hilic_ann <- hilic_ann %>%
  subset(timepoint == "01m" | timepoint == "06m") %>%
  mutate(timepoint = ifelse(timepoint == "01m", 1,
                            ifelse(timepoint == "06m", 6, NA))) %>%
  dplyr::select(-"studyID") %>%
  dplyr::select(c(dyad_id, timepoint, everything()))

fecalMetabolites_hilic     <- colnames(hilic_ann)[3:ncol(hilic_ann)]
all_fecalMetabolites_hilic <- fecalMetabolites_hilic

# Remove exogenous compounds
fecalMetabolites_hilic <- fecalMetabolites_hilic[!fecalMetabolites_hilic %in% exclude_hilic]

actually_excluded_hilic <- exclude_hilic[exclude_hilic %in% all_fecalMetabolites_hilic]
not_detected_hilic      <- exclude_hilic[!exclude_hilic %in% all_fecalMetabolites_hilic]

cat("HILIC compounds present and removed:", paste(actually_excluded_hilic, collapse = ", "), "\n")
cat("HILIC compounds in exclusion list but not in data:", paste(not_detected_hilic, collapse = ", "), "\n")
cat("HILIC metabolites remaining after exogenous removal:", length(fecalMetabolites_hilic), "\n")

# Convert to numeric
hilic_ann[, fecalMetabolites_hilic] <- lapply(
  hilic_ann[, fecalMetabolites_hilic], function(x) as.numeric(x)
)

# Missingness check
na_check_hilic <- colSums(is.na(hilic_ann[, fecalMetabolites_hilic]))
cat("HILIC metabolites with any NA:", sum(na_check_hilic > 0), "\n")
if (sum(na_check_hilic > 0) > 0) {
  cat("Metabolites with NAs:\n")
  print(na_check_hilic[na_check_hilic > 0])
}


# =============================================================================
# STEP 2: INTENSITY-PROFILE DEDUPLICATION — HILIC
# =============================================================================

cat("\n=== STEP 2: INTENSITY-PROFILE DEDUPLICATION — HILIC ===\n")
cat("HILIC metabolites before Step 2:", length(fecalMetabolites_hilic), "\n")

hilic_met_only  <- hilic_ann[, fecalMetabolites_hilic]
t_hilic         <- t(hilic_met_only)
t_hilic_dedup   <- t_hilic[!duplicated(t_hilic), ]
hilic_met_dedup <- t(t_hilic_dedup)

dropped_hilic <- setdiff(fecalMetabolites_hilic, colnames(hilic_met_dedup))
kept_hilic    <- colnames(hilic_met_dedup)

cat("HILIC metabolites after Step 2:", length(kept_hilic), "\n")
cat("HILIC metabolites dropped by Step 2:", length(dropped_hilic), "\n")
if (length(dropped_hilic) > 0) {
  cat("Dropped:", paste(dropped_hilic, collapse = ", "), "\n")
}
cat("\n")

non_met_cols_hilic  <- colnames(hilic_ann)[!colnames(hilic_ann) %in% fecalMetabolites_hilic]
hilic_ann           <- cbind(hilic_ann[, non_met_cols_hilic], hilic_met_dedup)
fecalMetabolites_hilic <- kept_hilic

# Build HILIC display name lookup programmatically from Step 2 dropped metabolites

hilic_display_df <- data.frame(
  CNAME        = kept_hilic,
  display_name = kept_hilic,
  stringsAsFactors = FALSE
)

for (kept in kept_hilic) {
  kept_row <- hilic_full[hilic_full$CNAME == kept &
                           !duplicated(hilic_full$CNAME), ]
  if (nrow(kept_row) == 0) next
  
  mz_val <- kept_row$detected_mz[1]
  rt_val <- kept_row$detected_rt[1]
  
  all_at_peak <- hilic_full[
    abs(hilic_full$detected_mz - mz_val) < 1e-8 &
      abs(hilic_full$detected_rt - rt_val) < 1e-8,
    "CNAME"
  ]
  
  unique_names <- unique(all_at_peak)
  if (length(unique_names) > 1) {
    idx <- which(hilic_display_df$CNAME == kept)
    hilic_display_df$display_name[idx] <- paste(unique_names, collapse = "; ")
  }
}

write.csv(hilic_display_df,
          here::here("out_files", "hilic_display_names.csv"),
          row.names = FALSE)
cat("HILIC display name lookup saved:", nrow(hilic_display_df), "metabolites\n")
cat("HILIC co-eluter groups in display lookup:",
    sum(grepl(";", hilic_display_df$display_name)), "\n")
cat("HILIC co-eluter examples:\n")
print(hilic_display_df[grepl(";", hilic_display_df$display_name), ])

# =============================================================================
# DIAGNOSTIC: Final metabolite lists after all deduplication and exclusion
# =============================================================================

cat("\n=== FINAL METABOLITE LISTS AFTER ALL PROCESSING ===\n")
cat("C18 final metabolite count:", length(fecalMetabolites_c18), "\n")
cat("HILIC final metabolite count:", length(fecalMetabolites_hilic), "\n")
cat("Combined total:", length(fecalMetabolites_c18) + length(fecalMetabolites_hilic), "\n\n")

cat("C18 final metabolites:\n")
print(sort(fecalMetabolites_c18))

cat("\nHILIC final metabolites:\n")
print(sort(fecalMetabolites_hilic))

# Check for overlapping names between platforms
cat("\nOverlapping metabolite names between platforms (before suffixing):\n")
print(intersect(fecalMetabolites_c18, fecalMetabolites_hilic))


# =============================================================================
# LOAD PFAS DATA
# =============================================================================

PFAS <- read.csv(here::here("out_files", "PFAS_1m6m_quantified.csv"))
if ("X" %in% colnames(PFAS)) {
  PFAS <- PFAS %>%
    dplyr::rename("timepoint" = Timepoint) %>%
    dplyr::select(-X)
} else {
  PFAS <- PFAS %>%
    dplyr::rename("timepoint" = Timepoint)
  message("Column 'X' does not exist. Skipping removal.")
}

PFAS <- PFAS %>%
  dplyr::select(dyad_id, timepoint, contains("_pgmL"))

mdl_values <- c(
  "N.EtFOSAA_pgmL" = 1.12, "N.MeFOSAA_pgmL" = 0.74, "PFBA_pgmL"   = 3.33,
  "PFDA_pgmL"      = 1.02, "PFDoA_pgmL"      = 0.74, "PFHpA_pgmL"  = 1.16,
  "PFHxA_pgmL"     = 3.93, "PFNA_pgmL"       = 0.66, "PFOA_pgmL"   = 2.28,
  "PFPeA_pgmL"     = 4.47, "PFTeDA_pgmL"     = 0.41, "PFUnA_pgmL"  = 1.16,
  "PFBS_pgmL"      = 7.00, "PFDoS_pgmL"      = 0.83, "PFHps_pgmL"  = 1.30,
  "PFHxS_pgmL"     = 0.62, "PFNS_pgmL"       = 1.44, "PFTrDA_pgmL" = 0.83,
  "PFOS_pgmL"      = 1.09, "PFPeAS_pgmL"     = 0.69
)

pfas_cols <- intersect(names(mdl_values), colnames(PFAS))

PFAS_long <- PFAS %>%
  dplyr::select(timepoint, all_of(pfas_cols)) %>%
  pivot_longer(cols = all_of(pfas_cols),
               names_to  = "PFAS_names",
               values_to = "value") %>%
  mutate(MDL = mdl_values[PFAS_names])

percent_below_mdl_by_tp <- PFAS_long %>%
  dplyr::group_by(PFAS_names, timepoint) %>%
  dplyr::summarize(
    n                 = sum(!is.na(value)),
    percent_below_MDL = mean(value < MDL, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  dplyr::mutate(percent_below_MDL = round(percent_below_MDL, 1))

eligible_pfas_1m <- percent_below_mdl_by_tp %>%
  filter(timepoint == 1, percent_below_MDL <= 25) %>%
  pull(PFAS_names)

eligible_pfas_1mDetect <- percent_below_mdl_by_tp %>%
  filter(timepoint == 1,
         percent_below_MDL > 25,
         percent_below_MDL < 80) %>%
  pull(PFAS_names)

PFAS_1m <- PFAS %>%
  filter(timepoint == 1) %>%
  dplyr::select(dyad_id, timepoint, all_of(as.character(eligible_pfas_1m)))

PFAS_1mDetect <- PFAS %>%
  filter(timepoint == 1) %>%
  dplyr::select(dyad_id, timepoint, all_of(as.character(eligible_pfas_1mDetect)))

unique(PFAS_1m$dyad_id)       %in% unique(meta_trim$dyad_id)
unique(PFAS_1mDetect$dyad_id) %in% unique(meta_trim$dyad_id)

PFAS_1m <- PFAS_1m %>%
  inner_join(meta_trim, by = c("dyad_id", "timepoint"))

PFAS_1mDetect <- PFAS_1mDetect %>%
  inner_join(meta_trim, by = c("dyad_id", "timepoint"))


# Create detect/non-detect columns ---------------------------------------------
create_detect_columns <- function(df, mdl_vals) {
  pfas_cols <- grep("_pgmL$", colnames(df), value = TRUE)
  for (pfas in pfas_cols) {
    mdl         <- mdl_vals[pfas]
    new_colname <- gsub("_pgmL$", "_Detect", pfas)
    df[[new_colname]] <- ifelse(df[[pfas]] >= mdl, "detect", "non-detect")
    cat(pfas, "->", new_colname, "(MDL =", mdl, ")\n")
  }
  return(df)
}

PFAS_1mDetect <- create_detect_columns(PFAS_1mDetect, mdl_values) %>%
  dplyr::select(-ends_with("_pgmL"))

grep("_Detect$", colnames(PFAS_1mDetect), value = TRUE)

create_n_detect <- function(df) {
  detect_cols    <- grep("_Detect$", colnames(df), value = TRUE)
  detect_numeric <- df[, detect_cols]
  detect_numeric[] <- lapply(detect_numeric, function(x) as.numeric(x == "detect"))
  df$n_detect <- rowSums(detect_numeric, na.rm = TRUE)
  cat("n_detect range:", min(df$n_detect, na.rm = TRUE), "-",
      max(df$n_detect, na.rm = TRUE), "\n")
  return(df)
}

PFAS_1mDetect <- create_n_detect(PFAS_1mDetect)


# Merge PFAS with fecal metabolomics -------------------------------------------

PFAS_1m       <- PFAS_1m       %>% mutate(dyad_id = as.character(dyad_id))
PFAS_1mDetect <- PFAS_1mDetect %>% mutate(dyad_id = as.character(dyad_id))
c18_ann       <- c18_ann       %>% mutate(dyad_id = as.character(dyad_id))
hilic_ann     <- hilic_ann     %>% mutate(dyad_id = as.character(dyad_id))

c18_1m   <- c18_ann   %>% filter(timepoint == 1)
c18_6m   <- c18_ann   %>% filter(timepoint == 6)
hilic_1m <- hilic_ann %>% filter(timepoint == 1)
hilic_6m <- hilic_ann %>% filter(timepoint == 6)

cat("C18  1m samples:", nrow(c18_1m),   "| 6m samples:", nrow(c18_6m),   "\n")
cat("HILIC 1m samples:", nrow(hilic_1m), "| 6m samples:", nrow(hilic_6m), "\n")

PFAS_1m_temp       <- PFAS_1m       %>% mutate(timepoint = 6)
PFAS_1mDetect_temp <- PFAS_1mDetect %>% mutate(timepoint = 6)

PFAS1m_c18_1m         <- PFAS_1m            %>% inner_join(c18_1m,   by = c("dyad_id", "timepoint"))
PFAS1m_c18_6m         <- PFAS_1m_temp       %>% inner_join(c18_6m,   by = c("dyad_id", "timepoint"))
PFAS1m_hilic_1m       <- PFAS_1m            %>% inner_join(hilic_1m, by = c("dyad_id", "timepoint"))
PFAS1m_hilic_6m       <- PFAS_1m_temp       %>% inner_join(hilic_6m, by = c("dyad_id", "timepoint"))
PFAS1mDetect_c18_1m   <- PFAS_1mDetect      %>% inner_join(c18_1m,   by = c("dyad_id", "timepoint"))
PFAS1mDetect_c18_6m   <- PFAS_1mDetect_temp %>% inner_join(c18_6m,   by = c("dyad_id", "timepoint"))
PFAS1mDetect_hilic_1m <- PFAS_1mDetect      %>% inner_join(hilic_1m, by = c("dyad_id", "timepoint"))
PFAS1mDetect_hilic_6m <- PFAS_1mDetect_temp %>% inner_join(hilic_6m, by = c("dyad_id", "timepoint"))


# Check missing values in covariates -------------------------------------------
covariates_to_check <- c("breastmilk_per_day", "SES_index_final",
                         "baby_birthweight_kg", "age_of_solid_foods",
                         "gestational_age_cat", "mode_of_delivery_cat")

for (nm in c("PFAS1m_c18_1m", "PFAS1m_c18_6m",
             "PFAS1m_hilic_1m", "PFAS1m_hilic_6m",
             "PFAS1mDetect_c18_1m", "PFAS1mDetect_c18_6m",
             "PFAS1mDetect_hilic_1m", "PFAS1mDetect_hilic_6m")) {
  df        <- get(nm)
  na_counts <- sapply(df[, covariates_to_check], function(x) sum(is.na(x)))
  total_na  <- sum(na_counts)
  cat(nm, "— total covariate NAs:", total_na, "\n")
  if (total_na > 0) print(na_counts[na_counts > 0])
}


# Covariate imputation ---------------------------------------------------------
get_mode <- function(x) {
  ux <- na.omit(x)
  ux[which.max(tabulate(match(ux, ux)))]
}

impute_covariates <- function(df) {
  if (sum(is.na(df$SES_index_final)) > 0)
    df$SES_index_final[is.na(df$SES_index_final)] <-
      median(df$SES_index_final, na.rm = TRUE)
  if (sum(is.na(df$baby_birthweight_kg)) > 0)
    df$baby_birthweight_kg[is.na(df$baby_birthweight_kg)] <-
      median(df$baby_birthweight_kg, na.rm = TRUE)
  if (sum(is.na(df$breastmilk_per_day)) > 0)
    df$breastmilk_per_day[is.na(df$breastmilk_per_day)] <-
      median(df$breastmilk_per_day, na.rm = TRUE)
  if (sum(is.na(df$age_of_solid_foods)) > 0)
    df$age_of_solid_foods[is.na(df$age_of_solid_foods)] <-
      median(df$age_of_solid_foods, na.rm = TRUE)
  if (sum(is.na(df$mode_of_delivery_cat)) > 0)
    df$mode_of_delivery_cat[is.na(df$mode_of_delivery_cat)] <-
      get_mode(df$mode_of_delivery_cat)
  if (sum(is.na(df$gestational_age_cat)) > 0)
    df$gestational_age_cat[is.na(df$gestational_age_cat)] <-
      get_mode(df$gestational_age_cat)
  return(df)
}

PFAS1m_c18_1m          <- impute_covariates(PFAS1m_c18_1m)
PFAS1m_c18_6m          <- impute_covariates(PFAS1m_c18_6m)
PFAS1m_hilic_1m        <- impute_covariates(PFAS1m_hilic_1m)
PFAS1m_hilic_6m        <- impute_covariates(PFAS1m_hilic_6m)
PFAS1mDetect_c18_1m    <- impute_covariates(PFAS1mDetect_c18_1m)
PFAS1mDetect_c18_6m    <- impute_covariates(PFAS1mDetect_c18_6m)
PFAS1mDetect_hilic_1m  <- impute_covariates(PFAS1mDetect_hilic_1m)
PFAS1mDetect_hilic_6m  <- impute_covariates(PFAS1mDetect_hilic_6m)


# PFAS processing --------------------------------------------------------------
process_pfas <- function(df, pfas_cols, mdl_vals) {
  for (pfas in pfas_cols) {
    mdl <- mdl_vals[pfas]
    df[[pfas]] <- ifelse(
      !is.na(df[[pfas]]) & df[[pfas]] < mdl,
      mdl / sqrt(2),
      df[[pfas]]
    )
  }
  df[pfas_cols] <- log2(df[pfas_cols])
  return(df)
}

pfas_continuous_cols <- grep("_pgmL$", colnames(PFAS1m_c18_1m), value = TRUE)
cat("Continuous PFAS columns:", pfas_continuous_cols, "\n")

PFAS1m_c18_1m   <- process_pfas(PFAS1m_c18_1m,   pfas_continuous_cols, mdl_values)
PFAS1m_c18_6m   <- process_pfas(PFAS1m_c18_6m,   pfas_continuous_cols, mdl_values)
PFAS1m_hilic_1m <- process_pfas(PFAS1m_hilic_1m, pfas_continuous_cols, mdl_values)
PFAS1m_hilic_6m <- process_pfas(PFAS1m_hilic_6m, pfas_continuous_cols, mdl_values)


# Metabolomics prevalence filter -----------------------------------------------
# Note: applied here for consistency; if data is fully imputed with no NAs,
# all metabolites will pass the 25% threshold and no metabolites will be dropped.
# The missingness check above will confirm this.

apply_prevalence_filter <- function(metabo_1m, metabo_6m, metab_cols,
                                    platform_name, threshold = 0.25) {
  combined   <- bind_rows(
    metabo_1m %>% dplyr::select(all_of(metab_cols)),
    metabo_6m %>% dplyr::select(all_of(metab_cols))
  )
  prevalence <- colMeans(!is.na(combined))
  keep_cols  <- names(prevalence[prevalence >= threshold])
  drop_cols  <- names(prevalence[prevalence <  threshold])
  cat(platform_name, "— metabolites before 25% filter:", length(metab_cols), "\n")
  cat(platform_name, "— metabolites after  25% filter:", length(keep_cols), "\n")
  cat(platform_name, "— metabolites dropped:", length(drop_cols), "\n")
  if (length(drop_cols) > 0) cat("Dropped:", paste(drop_cols, collapse = ", "), "\n")
  return(keep_cols)
}

keep_c18 <- apply_prevalence_filter(
  c18_1m, c18_6m, fecalMetabolites_c18, "C18"
)

keep_hilic <- apply_prevalence_filter(
  hilic_1m, hilic_6m, fecalMetabolites_hilic, "HILIC"
)

keep_cols_c18 <- function(df) {
  meta_cols <- colnames(df)[!colnames(df) %in% all_fecalMetabolites_c18]
  dplyr::select(df, all_of(c(meta_cols, keep_c18)))
}

keep_cols_hilic <- function(df) {
  meta_cols <- colnames(df)[!colnames(df) %in% all_fecalMetabolites_hilic]
  dplyr::select(df, all_of(c(meta_cols, keep_hilic)))
}

PFAS1m_c18_1m          <- keep_cols_c18(PFAS1m_c18_1m)
PFAS1m_c18_6m          <- keep_cols_c18(PFAS1m_c18_6m)
PFAS1mDetect_c18_1m    <- keep_cols_c18(PFAS1mDetect_c18_1m)
PFAS1mDetect_c18_6m    <- keep_cols_c18(PFAS1mDetect_c18_6m)
PFAS1m_hilic_1m        <- keep_cols_hilic(PFAS1m_hilic_1m)
PFAS1m_hilic_6m        <- keep_cols_hilic(PFAS1m_hilic_6m)
PFAS1mDetect_hilic_1m  <- keep_cols_hilic(PFAS1mDetect_hilic_1m)
PFAS1mDetect_hilic_6m  <- keep_cols_hilic(PFAS1mDetect_hilic_6m)


# NA check in PFAS after log2 transform ----------------------------------------
for (nm in c("PFAS1m_c18_1m", "PFAS1m_c18_6m",
             "PFAS1m_hilic_1m", "PFAS1m_hilic_6m")) {
  df           <- get(nm)
  pfas_present <- intersect(pfas_continuous_cols, colnames(df))
  na_counts    <- colSums(is.na(df[, pfas_present, drop = FALSE]))
  cat(nm, "— total PFAS NAs:", sum(na_counts), "\n")
  if (sum(na_counts) > 0) print(na_counts[na_counts > 0])
}

detect_cols_check <- grep("_Detect$", colnames(PFAS1mDetect_c18_1m), value = TRUE)
for (nm in c("PFAS1mDetect_c18_1m", "PFAS1mDetect_c18_6m",
             "PFAS1mDetect_hilic_1m", "PFAS1mDetect_hilic_6m")) {
  df     <- get(nm)
  na_det <- sum(sapply(df[detect_cols_check], function(x) sum(is.na(x))))
  cat(nm, "— detect column NAs:", na_det, "\n")
}


# Handle overlapping metabolite names between platforms ------------------------
overlap_metabs <- intersect(keep_c18, keep_hilic)
cat("Overlapping metabolite names between platforms:", length(overlap_metabs), "\n")
if (length(overlap_metabs) > 0) print(overlap_metabs)

colnames(PFAS1m_c18_1m)[colnames(PFAS1m_c18_1m)             %in% overlap_metabs] <- paste0(overlap_metabs, "_c18")
colnames(PFAS1m_c18_6m)[colnames(PFAS1m_c18_6m)             %in% overlap_metabs] <- paste0(overlap_metabs, "_c18")
colnames(PFAS1mDetect_c18_1m)[colnames(PFAS1mDetect_c18_1m) %in% overlap_metabs] <- paste0(overlap_metabs, "_c18")
colnames(PFAS1mDetect_c18_6m)[colnames(PFAS1mDetect_c18_6m) %in% overlap_metabs] <- paste0(overlap_metabs, "_c18")

colnames(PFAS1m_hilic_1m)[colnames(PFAS1m_hilic_1m)             %in% overlap_metabs] <- paste0(overlap_metabs, "_hilic")
colnames(PFAS1m_hilic_6m)[colnames(PFAS1m_hilic_6m)             %in% overlap_metabs] <- paste0(overlap_metabs, "_hilic")
colnames(PFAS1mDetect_hilic_1m)[colnames(PFAS1mDetect_hilic_1m) %in% overlap_metabs] <- paste0(overlap_metabs, "_hilic")
colnames(PFAS1mDetect_hilic_6m)[colnames(PFAS1mDetect_hilic_6m) %in% overlap_metabs] <- paste0(overlap_metabs, "_hilic")

keep_c18   <- ifelse(keep_c18   %in% overlap_metabs, paste0(keep_c18,   "_c18"),   keep_c18)
keep_hilic <- ifelse(keep_hilic %in% overlap_metabs, paste0(keep_hilic, "_hilic"), keep_hilic)

cat("Overlapping names after suffixing:", length(intersect(keep_c18, keep_hilic)), "\n")


# Save all 8 merged files ------------------------------------------------------
write.csv(PFAS1m_c18_1m,         here::here("out_files", "PFAS1m_c18_1m.csv"),         row.names = FALSE)
write.csv(PFAS1m_c18_6m,         here::here("out_files", "PFAS1m_c18_6m.csv"),         row.names = FALSE)
write.csv(PFAS1m_hilic_1m,       here::here("out_files", "PFAS1m_hilic_1m.csv"),       row.names = FALSE)
write.csv(PFAS1m_hilic_6m,       here::here("out_files", "PFAS1m_hilic_6m.csv"),       row.names = FALSE)
write.csv(PFAS1mDetect_c18_1m,   here::here("out_files", "PFAS1mDetect_c18_1m.csv"),   row.names = FALSE)
write.csv(PFAS1mDetect_c18_6m,   here::here("out_files", "PFAS1mDetect_c18_6m.csv"),   row.names = FALSE)
write.csv(PFAS1mDetect_hilic_1m, here::here("out_files", "PFAS1mDetect_hilic_1m.csv"), row.names = FALSE)
write.csv(PFAS1mDetect_hilic_6m, here::here("out_files", "PFAS1mDetect_hilic_6m.csv"), row.names = FALSE)

saveRDS(keep_c18,   here::here("out_files", "metabolite_cols_c18.rds"))
saveRDS(keep_hilic, here::here("out_files", "metabolite_cols_hilic.rds"))


# Save metabolite list for manual superclass annotation ------------------------
data.frame(
  CNAME    = c(keep_c18, keep_hilic),
  Platform = c(rep("C18",   length(keep_c18)),
               rep("HILIC", length(keep_hilic)))
) %>%
  distinct(CNAME, .keep_all = TRUE) %>%
  arrange(CNAME) %>%
  mutate(Super_Class = NA_character_,
         Class       = NA_character_) %>%
  write.csv(here::here("out_files", "metabolites_for_annotation.csv"),
            row.names = FALSE)

cat("\n=== SCRIPT COMPLETE ===\n")
cat("C18 final metabolites in MWAS:", length(keep_c18), "\n")
cat("HILIC final metabolites in MWAS:", length(keep_hilic), "\n")
cat("Combined total:", length(keep_c18) + length(keep_hilic), "\n")

# END --------------------------------------------------------------------------