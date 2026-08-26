# TITLE:   00_generate_synthetic_example_data.R
#
# PURPOSE: Generate SYNTHETIC example data that mimics the structure (column
#          names, data types, value ranges) of the real, IRB-restricted
#          Mother's Milk Study raw data files, so that reviewers can run the
#          analysis pipeline (scripts 0-16) end-to-end and see real output.
#
#          NO REAL PARTICIPANT DATA IS USED ANYWHERE IN THIS SCRIPT. All
#          subject IDs, PFAS concentrations, microbiome counts, and metadata
#          values below are randomly generated. Only non-human, public/
#          non-sensitive reference data (bacterial taxonomy dictionaries,
#          MetaCyc-style pathway identifiers) are reused for realism.
#
#          Real data underlying the published results are available as
#          described in the manuscript's Data Availability Statement
#          (NCBI accession for microbiome sequences; metabolomics data
#          deposition accession; other data available upon reasonable
#          request, subject to IRB data use agreement).
#
# USAGE:   Run once from the project root (the folder containing this
#          "code" directory) before running scripts 0-16:
#            Rscript code/00_generate_synthetic_example_data.R
#          This populates ./input/, and directly seeds a few ./out_files/
#          intermediates that stand in for the (not-reproduced) untargeted
#          metabolomics annotation pipeline (see script 9 note in README).

set.seed(2024)

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
  library(readr)
  library(vegan)
})

root <- here::here()
dir.create(file.path(root, "input", "missing_bmpercent"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "out_files"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "out_figures"), recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Synthetic cohort skeleton: dyads, timepoints, sample IDs
# =============================================================================

N_DYAD <- 40
dyad_ids <- 1:N_DYAD

# every dyad has 1-month and 6-month samples; first 5 dyads also have a
# 12-month sample, to exercise the "drop other timepoints" branch in script 1
manifest <- bind_rows(
  data.frame(dyad_id = dyad_ids, timepoint = 1L),
  data.frame(dyad_id = dyad_ids, timepoint = 6L),
  data.frame(dyad_id = dyad_ids[1:5], timepoint = 12L)
) %>%
  arrange(dyad_id, timepoint) %>%
  mutate(core_id = row_number(),
         mg_id   = paste0("MG", core_id),
         merge_id_dyad = sprintf("MM-%04d-%02d", dyad_id, timepoint))

n_control <- 10
control_ids <- paste0("L", 1:n_control)

cat("Synthetic cohort: ", N_DYAD, " dyads, ", nrow(manifest), " microbiome samples, ",
    n_control, " control samples\n", sep = "")

# =============================================================================
# 2. GORAN_MICROBIOME_MANIFEST 10-2022.xlsx
#    Two header rows (spanning header + true field names), matching the
#    real manifest layout that scripts 1/4/5/8 parse with colnames(IDs) <- IDs[1,]
# =============================================================================

wb <- createWorkbook()
addWorksheet(wb, "Sheet1")
writeData(wb, "Sheet1", t(c(NA, "Take from here", NA, NA, "Move to here", NA)),
          startRow = 1, colNames = FALSE)
writeData(wb, "Sheet1", t(c("CORE ID", "Old Together", "Box", "Position", "Box", "Position")),
          startRow = 2, colNames = FALSE)

old_together <- sprintf("%02d Month%d", manifest$timepoint, manifest$dyad_id)
manifest_sheet <- data.frame(
  core_id      = manifest$core_id,
  old_together = old_together,
  box1         = sample(1:10, nrow(manifest), replace = TRUE),
  position1    = paste0(sample(1:12, nrow(manifest), replace = TRUE), ".",
                         sample(LETTERS[1:8], nrow(manifest), replace = TRUE)),
  box2         = sample(1:10, nrow(manifest), replace = TRUE),
  position2    = paste0(sample(LETTERS[1:8], nrow(manifest), replace = TRUE),
                         sample(1:12, nrow(manifest), replace = TRUE))
)
writeData(wb, "Sheet1", manifest_sheet, startRow = 3, colNames = FALSE)
saveWorkbook(wb, file.path(root, "input", "GORAN_MICROBIOME_MANIFEST 10-2022.xlsx"), overwrite = TRUE)

# =============================================================================
# 3. Bracken count / relative-abundance tables at 6 taxonomic levels
#    Column headers = REAL taxonomy_id values sampled from the (public,
#    non-sensitive) taxonomy dictionaries, so the taxonomy join in script 1
#    behaves exactly as it does on real data.
# =============================================================================

read_tax_dict <- function(level) {
  f <- file.path(root, "input", paste0("taxonomyDictionary_brack_jan_", level, "_withLineage_bacteriaOnly.tsv"))
  d <- readr::read_tsv(f, show_col_types = FALSE)
  # exclude names with Newick-special characters (script 7 builds a Newick
  # tree from these names; a handful of real NCBI taxa, e.g. plant
  # phytoplasmas, contain apostrophes/parentheses that break tree parsing)
  bad_chars <- c("'", "\"", "(", ")", ",", ":", ";", "[", "]")
  has_bad_char <- Reduce(`|`, lapply(bad_chars, function(ch) grepl(ch, d$name, fixed = TRUE)))
  d[!has_bad_char, ]
}
read_tax_ids <- function(level) read_tax_dict(level)$taxonomy_id

# Genera that scripts 13/14 explicitly key on (PFAS-associated genus groups,
# and the in vitro PFOS growth-screening strain library). Real cohort
# species lists are gut-dominated, so a uniform random draw across all of
# bacterial taxonomy would rarely include these; we bias the synthetic
# species pool toward them so the genus-level matching logic in those
# scripts has real signal to work with.
gut_genera <- c("Bifidobacterium", "Streptococcus", "Escherichia", "Shigella", "Veillonella",
                "Bacteroides", "Phocaeicola", "Parabacteroides", "Blautia", "Dorea",
                "Anaerostipes", "Roseburia", "Coprococcus", "Lachnoclostridium", "Anaerobutyricum",
                "Clostridium", "Enterocloster", "Mediterraneibacter", "Ruminococcus", "Simiaoa",
                "Enterobacter", "Klebsiella")

pick_species_ids <- function(n_pick) {
  d <- read_tax_dict("species")
  gut_pattern <- paste0("^(", paste(gut_genera, collapse = "|"), ") ")
  is_gut <- grepl(gut_pattern, d$name)
  gut_pool   <- d$taxonomy_id[is_gut & !d$taxonomy_id %in% junk_ids$species]
  other_pool <- d$taxonomy_id[!is_gut & !d$taxonomy_id %in% junk_ids$species]
  n_gut <- min(length(gut_pool), round(n_pick * 0.6))
  c(sample(gut_pool, n_gut), sample(other_pool, n_pick - n_gut))
}

# taxon IDs that script 1 explicitly removes as non-bacterial contamination
junk_ids <- list(
  species = 9606,
  genus   = 9605,
  family  = 9604,
  order   = 9443,
  class   = c(183963, 2731619, 40674),
  phylum  = c(2731618, 28890, 7711)
)

# species needs a taxon pool large enough that, after the >=25% prevalence /
# >=50-count filter in script 1, a comfortably large number of species survive
n_taxa_per_level <- c(species = 260, genus = 60, family = 40, order = 25, class = 15, phylum = 8)

all_sample_ids <- c(control_ids, manifest$mg_id)
n_samples <- length(all_sample_ids)

# NOTE ON REAL FILE ORIENTATION: the real Bracken files are taxa-in-rows,
# samples-in-columns, with a (misleadingly named) "SampleID" first column
# that actually holds the taxonomy_id. Scripts 1/6 transpose this on load.
# We reproduce that exact on-disk layout here for fidelity.
make_bracken_tables <- function(level) {
  n_pick <- n_taxa_per_level[[level]]
  if (level == "species") {
    picked <- pick_species_ids(n_pick)
  } else {
    real_ids <- read_tax_ids(level)
    picked <- sample(real_ids[!real_ids %in% unlist(junk_ids)], n_pick)
  }
  taxa_ids <- c(picked, junk_ids[[level]])
  n_taxa <- length(taxa_ids)

  # power-law-ish relative weights so prevalence filtering has something to do
  weights <- sort(rexp(n_taxa, rate = 1), decreasing = TRUE)
  weights <- weights / sum(weights)

  total_reads <- round(runif(n_samples, 1.2e6, 3e6))
  # internal working matrix: samples (rows) x taxa (cols)
  counts <- t(sapply(seq_len(n_samples), function(i) {
    as.integer(rmultinom(1, total_reads[i], prob = weights))
  }))
  rownames(counts) <- all_sample_ids
  colnames(counts) <- as.character(taxa_ids)
  ra <- counts / rowSums(counts)

  # on-disk orientation: taxa (rows) x samples (cols)
  counts_df <- data.frame(SampleID = taxa_ids, t(counts), check.names = FALSE)
  colnames(counts_df)[-1] <- all_sample_ids
  ra_df <- data.frame(SampleID = taxa_ids, t(ra), check.names = FALSE)
  colnames(ra_df)[-1] <- all_sample_ids

  write.csv(counts_df, file.path(root, "input", paste0("counts_bracken_", level, ".csv")), row.names = FALSE)
  write.csv(ra_df, file.path(root, "input", paste0("relative_abundance_bracken_", level, ".csv")), row.names = FALSE)
  # return the samples x taxa matrices for internal use (e.g. alpha/beta diversity)
  invisible(list(counts = counts, ra = ra))
}

bracken_levels <- c("species", "genus", "family", "order", "class", "phylum")
bracken_out <- setNames(lapply(bracken_levels, make_bracken_tables), bracken_levels)

# CLR-transformed species abundance, keyed by merge_id_dyad, for injecting
# real taxa<->metabolite correlation into the synthetic metabolomics stand-in
# (section 9 below) -- otherwise every metabolite is pure noise and the
# taxa-metabolite / pathway-metabolite correlation scripts (11, 12, 14) have
# nothing to detect and skip their heatmaps.
species_clr <- {
  ra_mg <- bracken_out$species$ra[manifest$mg_id, , drop = FALSE]
  log_ra <- log(ra_mg + 1e-6)
  clr <- log_ra - rowMeans(log_ra)
  rownames(clr) <- manifest$merge_id_dyad
  clr
}
signal_taxa <- sample(colnames(species_clr), 8)

# =============================================================================
# 4. Alpha diversity (per sample) and beta diversity (Bray-Curtis distance)
#    Computed from the synthetic species relative-abundance matrix so the
#    values are internally consistent (not just random noise).
# =============================================================================

species_ra_mat <- bracken_out$species$ra

alpha_div <- data.frame(
  SampleID  = all_sample_ids,
  Shannon   = vegan::diversity(species_ra_mat, index = "shannon"),
  Richness  = rowSums(species_ra_mat > 0),
  Evenness  = vegan::diversity(species_ra_mat, index = "shannon") / log(rowSums(species_ra_mat > 0)),
  Simpson   = vegan::diversity(species_ra_mat, index = "simpson")
)
readr::write_tsv(alpha_div, file.path(root, "input",
  "alphaDiv_repeated_rarefied_100_data_readDepth_1000000_mothersMilk_replacementFALSE.tsv"))

beta_dist <- as.matrix(vegan::vegdist(species_ra_mat, method = "bray"))
beta_dist[upper.tri(beta_dist, diag = FALSE)] <- NA
beta_df <- data.frame(rn = rownames(beta_dist), beta_dist, check.names = FALSE)
colnames(beta_df)[-1] <- colnames(beta_dist)
readr::write_tsv(beta_df, file.path(root, "input",
  "betaDiv_bray_distance_repeated_rarefied_100_data_readDepth_1000000_mothersMilk_replacementFALSE.tsv"),
  col_names = TRUE, na = "")

# =============================================================================
# 5. HUMAnN pathway abundance tables (CPM); pathway names are real public
#    MetaCyc-style identifiers (not sensitive), abundances are synthetic.
# =============================================================================

pathway_names <- c(
  "UNMAPPED", "UNINTEGRATED",
  "PWY-6737: starch degradation V", "PWY-6737: starch degradation V|unclassified",
  "PWY-7111: pyruvate fermentation to isobutanol",
  "COA-PWY: coenzyme A biosynthesis I",
  "PWY-5484: glycolysis II", "PWY-6588: pyruvate fermentation to acetone",
  "ANAGLYCOLYSIS-PWY: glycolysis III", "PWY-6572: chondroitin sulfate degradation I",
  "P461-PWY: hexitol fermentation to lactate, formate, ethanol and acetate",
  "PWY-7237: myo-, chiro- and scillo-inositol degradation",
  "BIFIDOSHUNT-PWY: bifidobacterium shunt",
  "PWY-7456: mannan degradation",
  "UDPNAGSYN-PWY: UDP-N-acetylglucosamine biosynthesis I",
  "PWY-7013: L-1,2-propanediol degradation",
  "COLANSYN-PWY: colanic acid building blocks biosynthesis",
  "RHAMCAT-PWY: L-rhamnose degradation I"
)
pathway_short_names <- c(
  "UNMAPPED", "UNINTEGRATED",
  "Starch degradation V", "Starch degradation V (unclassified)",
  "Pyruvate ferm. to isobutanol",
  "CoA biosynthesis I",
  "Glycolysis II", "Pyruvate ferm. to acetone",
  "Glycolysis III", "Chondroitin sulfate degradation I",
  "Hexitol fermentation",
  "Inositol degradation",
  "Bifidobacterium shunt",
  "Mannan degradation",
  "UDP-GlcNAc biosynthesis I",
  "Propanediol degradation",
  "Colanic acid biosynthesis",
  "L-rhamnose degradation I"
)
write.csv(data.frame(Pathway = pathway_names, Short_Name = pathway_short_names),
          file.path(root, "out_files", "pathway_short_names.csv"), row.names = FALSE)

# Stand-in for PATHWAY_significant_with_species.csv (script 14). In the real
# pipeline this is a manually-curated cross-reference of which species
# contribute to which significant HUMAnN pathway (from stratified pathway
# output); no script produces it, so we build a small structurally-matching
# synthetic version here directly (real pathway/species names, random stats).
real_pwy <- pathway_names[-(1:2)]  # drop UNMAPPED/UNINTEGRATED
pwy_tax_dict <- read_tax_dict("species")
pwy_species_pool <- pwy_tax_dict[grepl(paste0("^(", paste(gut_genera, collapse = "|"), ") "), pwy_tax_dict$name), ]
categories <- c("Cell Wall", "Amino Acids", "Carbohydrates", "Vitamins & Cofactors",
                "Nucleotides", "Other", "Lipids", "SCFAs")
predictors <- c("N-detect", "Mixture", "PFBS", "PFHxS", "PFNA", "PFOA", "PFOS",
                 "N.MeFOSAA", "PFBA", "PFDA", "PFDoA", "PFHpA")
scenarios <- c("Binary 1m PFAS + 1m Pathway", "Binary 1m PFAS + 6m Pathway",
               "Continuous 1m PFAS + 1m Pathway", "Continuous 1m PFAS + 6m Pathway")

pwy_species_rows <- lapply(real_pwy, function(pwy) {
  n_sp <- sample(3:8, 1)
  sp <- pwy_species_pool[sample(nrow(pwy_species_pool), n_sp), ]
  beta <- rnorm(1, sd = 0.15)
  data.frame(
    predictor = sample(predictors, 1),
    pathway_clean = sub("^[A-Za-z0-9_-]+: ", "", pwy),
    betas = beta, ci_lo = beta - 0.1, ci_hi = beta + 0.1,
    p_val = runif(1, 0, 0.049), FDR = runif(1, 0, 0.09),
    scenario = sample(scenarios, 1),
    species_name = sp$name,
    genus = sub(" .*", "", sp$name),
    Direction = if (beta >= 0) "Positive" else "Negative",
    Category = sample(categories, 1)
  )
})
write.csv(bind_rows(pwy_species_rows),
          file.path(root, "out_files", "PATHWAY_significant_with_species.csv"), row.names = FALSE)

n_pwy <- length(pathway_names)
pwy_mat <- matrix(rexp(n_pwy * n_samples, rate = 1e-4), nrow = n_pwy)
pwy_mat[1, ] <- runif(n_samples, 5e4, 1.6e5)   # UNMAPPED
pwy_mat[2, ] <- runif(n_samples, 7e5, 8.7e5)   # UNINTEGRATED
colnames(pwy_mat) <- all_sample_ids
pwy_collapsed <- data.frame(PathwayName = pathway_names, pwy_mat, check.names = FALSE)
readr::write_tsv(pwy_collapsed, file.path(root, "input", "mothersmilk_humann_pathabundances_cpm_collapsed.tsv"))

# stratified version: reuse same rows, add "|unclassified" contributor already present
long_names <- paste0(pathway_names, "_stub_", seq_along(pathway_names) %% 4)
colnames(pwy_mat) <- paste0(all_sample_ids, "_stub.trimmed_clean_concat_Abundance")
pwy_full <- data.frame(`# Pathway` = pathway_names, pwy_mat, check.names = FALSE)
readr::write_tsv(pwy_full, file.path(root, "input", "mothersmilk_humann_pathabundances_cpm.tsv"))

# =============================================================================
# 6. Raw metadata CSV (REDCap-style). Real file has ~3,650 columns; only the
#    ~17 columns actually consumed by the analysis scripts are reproduced.
# =============================================================================

meta <- manifest %>%
  transmute(
    merge_id_dyad,
    dyad_id  = as.character(dyad_id),
    study_id = as.character(dyad_id),
    timepoint,
    mother_age = round(rnorm(n(), 31, 5)),
    prepreg_wt_kg = round(rnorm(n(), 70, 12), 1),
    prepreg_bmi_kgm2 = round(rnorm(n(), 26, 5), 1),
    baby_gender = sample(1:2, n(), replace = TRUE),
    baby_antibiotics = sample(1:4, n(), replace = TRUE, prob = c(0.75, 0.15, 0.05, 0.05)),
    mother_antibiotics = sample(1:4, n(), replace = TRUE, prob = c(0.75, 0.1, 0.05, 0.1)),
    mode_of_delivery = sample(1:2, n(), replace = TRUE, prob = c(0.75, 0.25)),
    breastmilk_per_day = sample(0:10, n(), replace = TRUE),
    age_in_days = timepoint * 30 + sample(-3:3, n(), replace = TRUE),
    age_solid_foods_clean = ifelse(timepoint >= 6, round(runif(n(), 4, 6), 1), NA),
    gestational_age_category = sample(c("<38", "38-40", "40-42"), n(), replace = TRUE, prob = c(0.15, 0.6, 0.25)),
    baby_birthweight_kg = round(rnorm(n(), 3.3, 0.5), 2),
    SES_index_final = round(rnorm(n(), 0, 1), 2)
  ) %>%
  # cleaned/derived covariates, mirroring script 1's transformations, so the
  # section-9 metabolomics stand-in (below) can reuse the SAME cleaned
  # column names/values that scripts 10-12/14 expect (matching out_files/
  # meta_trim.csv, which script 1 produces from the raw columns above)
  mutate(
    mode_of_delivery_cat = factor(ifelse(mode_of_delivery == 1, "Vaginal", "C-Section"),
                                   levels = c("Vaginal", "C-Section")),
    gestational_age_cat = factor(
      case_when(
        gestational_age_category == "<38"   ~ "Early",
        gestational_age_category == "38-40" ~ "Ontime",
        gestational_age_category == "40-42" ~ "Late"
      ),
      levels = c("Early", "Ontime", "Late")
    )
  )

write.csv(meta, file.path(root, "input", "mothersMilk_metadata_timepointsAsRows_updated051022_Temporary24mDiet.csv"),
          row.names = FALSE)

# =============================================================================
# 7. missing_bmpercent/bmpercent_missing.xlsx (24-hr diet recall lookup)
# =============================================================================

bm_missing <- data.frame(
  dyad_id = as.character(sample(dyad_ids, 8)),
  timepoint = 6,
  `1st recall BF/day` = sample(0:6, 8, replace = TRUE),
  `1st recall FF/day` = sample(0:3, 8, replace = TRUE),
  `2nd recall BF/day` = sample(0:6, 8, replace = TRUE),
  `2nd recall FF/day` = sample(0:3, 8, replace = TRUE),
  breastmilk_per_day = sample(1:6, 8, replace = TRUE),
  formula_per_day = sample(0:3, 8, replace = TRUE),
  check.names = FALSE
)
write.xlsx(bm_missing, file.path(root, "input", "missing_bmpercent", "bmpercent_missing.xlsx"))

# =============================================================================
# 8. Raw targeted-PFAS lab report (CLU0031_PFAS_Final_Report.xlsx)
#    Layout matches what script "0. Cleaning_targeted_PFAS.R" expects:
#    row2 (cols 6+) = compound names, row4 (cols 1-5) = meta field names,
#    row5+ = data.
# =============================================================================

pfas_compounds <- c("N-EtFOSAA", "N-MeFOSAA", "PFBA", "PFDA", "PFDoA", "PFHpA", "PFHxA",
                     "PFNA", "PFOA", "PFPeA", "PFTeDA", "PFUnA", "PFBS", "PFDoS", "PFHps",
                     "PFHxS", "PFNS", "PFTrDA", "PFOS", "PFPeAS")
mdl_values <- c(1.12, 0.74, 3.33, 1.02, 0.74, 1.16, 3.93, 0.66, 2.28, 4.47, 0.41, 1.16,
                7.00, 0.83, 1.30, 0.62, 1.44, 0.83, 1.09, 0.69)

sample_rows <- bind_rows(
  data.frame(dyad_id = dyad_ids, suffix = "BL"),
  data.frame(dyad_id = dyad_ids, suffix = "6M"),
  data.frame(dyad_id = dyad_ids[1:5], suffix = "12M")
)
sample_rows$Sample.ID <- paste0(sample_rows$dyad_id, "_", sample_rows$suffix)
n_pfas_rows <- nrow(sample_rows)

# Detection tier is assigned deterministically BY COMPOUND NAME, matching
# the exact compounds that downstream scripts (2-10) hardcode as literal
# column names (e.g. `pfas_cols <- c("PFBS_pgmL", "PFHxS_pgmL", "PFNA_pgmL",
# "PFOA_pgmL", "PFOS_pgmL")` in script 3, `pfas_vars_binary <- c("N.MeFOSAA_
# Detect", "PFBA_Detect", "PFDA_Detect", "PFDoA_Detect", "PFHpA_Detect")` in
# script 10). These names came from which real compounds happened to clear
# the <=25% / 25-80% below-MDL thresholds in the real cohort; we reproduce
# the same names/tiers here so the scripts run unmodified against synthetic
# data.
names(mdl_values) <- pfas_compounds
high_tier <- c("PFBS", "PFHxS", "PFNA", "PFOA", "PFOS")               # <=25% below MDL
mid_tier  <- c("N-MeFOSAA", "PFBA", "PFDA", "PFDoA", "PFHpA")          # 25-80% below MDL
tier <- ifelse(pfas_compounds %in% high_tier, "high",
        ifelse(pfas_compounds %in% mid_tier, "mid", "low"))
meanlog_mult <- c(high = 20, mid = 1, low = 1 / 20)

conc <- sapply(seq_along(mdl_values), function(i) {
  mdl <- mdl_values[i]
  round(rlnorm(n_pfas_rows, meanlog = log(mdl * meanlog_mult[[tier[i]]]), sdlog = 0.6), 3)
})
colnames(conc) <- pfas_compounds

meta_block <- data.frame(
  File.Name   = paste0("file_", seq_len(n_pfas_rows)),
  Sample.ID   = sample_rows$Sample.ID,
  Run_order   = seq_len(n_pfas_rows),
  Batch       = sample(1:3, n_pfas_rows, replace = TRUE),
  Sample_type = "Sample"
)

n_col_total <- 5 + length(pfas_compounds)
# NOTE: rows 1 and 3 need at least one non-blank cell, otherwise a fully-NA
# row is silently dropped when openxlsx writes the sheet (Excel does not
# store a fully empty leading/interior row). Script 0 never reads these
# cells, so the placeholder text is inert.
row1        <- c("Synthetic PFAS report (placeholder row)", rep(NA_character_, n_col_total - 1))
row2        <- c(rep(NA_character_, 5), pfas_compounds)
row3        <- c(rep(NA_character_, 5), "units: pg/mL", rep(NA_character_, n_col_total - 6))
row4        <- c(colnames(meta_block), rep(NA_character_, length(pfas_compounds)))
data_rows   <- as.matrix(cbind(meta_block, as.data.frame(conc)))
storage.mode(data_rows) <- "character"

sheet_mat <- rbind(row1, row2, row3, row4, data_rows)
dimnames(sheet_mat) <- NULL

wb2 <- createWorkbook()
addWorksheet(wb2, "Sheet1")
writeData(wb2, "Sheet1", sheet_mat, startRow = 1, colNames = FALSE)
saveWorkbook(wb2, file.path(root, "input", "CLU0031_PFAS_Final_Report.xlsx"), overwrite = TRUE)

# =============================================================================
# 9. Stand-in for the untargeted metabolomics pipeline (script 9 outputs)
#    Script 9 parses raw LC-MS annotation exports and RData preprocessing
#    objects that are highly specialized to this study's chemistry and are
#    NOT reproduced here (see README). Instead we synthesize its *outputs*
#    directly, with matching structure, so scripts 10-12 and 14 (which
#    consume these files) can still be run end-to-end.
# =============================================================================

make_metabolite_block <- function(platform, n_metab, meta_df) {
  metab_names <- paste0("METAB_", toupper(platform), "_", sprintf("%03d", seq_len(n_metab)))
  vals <- matrix(rnorm(nrow(meta_df) * n_metab, mean = 15, sd = 2), nrow = nrow(meta_df))
  colnames(vals) <- metab_names

  # Inject real correlation with several species CLR values into the first
  # few metabolite columns, so that IF a signal taxon and a signal
  # metabolite both independently clear the PFAS-significance pre-filter in
  # scripts 6/10, script 11's direct taxa<->metabolite correlation actually
  # finds real signal (pure-noise metabolites never survive its FDR
  # correction, so the heatmap would otherwise always be empty). Script 11
  # additionally requires a DENSE mutual signal (every qualifying taxon
  # correlated with several qualifying metabolites and vice versa), so each
  # signal metabolite here is a combination of MULTIPLE signal taxa (not
  # just one), forming a fully-connected taxa x metabolite signal clique.
  n_signal <- min(2 * length(signal_taxa), n_metab)
  clr_here <- species_clr[meta_df$merge_id_dyad, signal_taxa, drop = FALSE]
  combo <- rowMeans(clr_here)
  for (i in seq_len(n_signal)) {
    vals[, i] <- 15 + 4 * combo + rnorm(nrow(meta_df), sd = 1)
  }

  cbind(meta_df, as.data.frame(vals))
}

for (tp in c(1, 6)) {
  meta_tp <- meta %>% filter(timepoint == tp) %>%
    select(dyad_id, timepoint, merge_id_dyad, age_in_days, baby_gender, mother_age,
           baby_birthweight_kg, prepreg_bmi_kgm2, prepreg_wt_kg, baby_antibiotics,
           mother_antibiotics, mode_of_delivery_cat, gestational_age_cat,
           breastmilk_per_day, SES_index_final, age_solid_foods_clean)

  for (platform in c("c18", "hilic")) {
    n_metab <- if (platform == "c18") 50 else 40
    df <- make_metabolite_block(platform, n_metab, meta_tp)
    # PFAS columns actually used at 1-month exposure (per script 1's MDL screen)
    df$PFOA_pgmL <- round(rlnorm(nrow(df), log(3), 0.8), 3)
    df$PFOS_pgmL <- round(rlnorm(nrow(df), log(2), 0.8), 3)
    df$PFHxS_pgmL <- round(rlnorm(nrow(df), log(1), 0.8), 3)
    df$PFNA_pgmL <- round(rlnorm(nrow(df), log(1), 0.8), 3)
    df$PFBS_pgmL <- round(rlnorm(nrow(df), log(1), 0.8), 3)
    write.csv(df, file.path(root, "out_files", sprintf("PFAS1m_%s_%dm.csv", platform, tp)), row.names = FALSE)

    df_detect <- df
    # matches pfas_vars_binary hardcoded in script 10 (MWAS.R)
    for (p in c("N.MeFOSAA", "PFBA", "PFDA", "PFDoA", "PFHpA")) {
      df_detect[[paste0(p, "_Detect")]] <- sample(c("detect", "non-detect"), nrow(df_detect), replace = TRUE)
    }
    df_detect$n_detect <- rowSums(sapply(grep("_Detect$", names(df_detect)), function(i) df_detect[[i]] == "detect"))
    write.csv(df_detect, file.path(root, "out_files", sprintf("PFAS1mDetect_%s_%dm.csv", platform, tp)), row.names = FALSE)
  }
}

saveRDS(paste0("METAB_C18_", sprintf("%03d", 1:50)), file.path(root, "out_files", "metabolite_cols_c18.rds"))
saveRDS(paste0("METAB_HILIC_", sprintf("%03d", 1:40)), file.path(root, "out_files", "metabolite_cols_hilic.rds"))

write.csv(data.frame(CNAME = paste0("METAB_C18_", sprintf("%03d", 1:50)),
                      display_name = paste0("Synthetic metabolite C18-", 1:50)),
          file.path(root, "out_files", "c18_display_names.csv"), row.names = FALSE)

# Chemical superclass/class lookup (script 12). Category NAMES are reused
# from real, public ClassyFire-style chemical taxonomy (not sensitive);
# only the CNAME identifiers (our synthetic metabolite names) are assigned.
superclass_pairs <- data.frame(
  Super_Class = c("Lipids and lipid-like molecules", "Organic acids and derivatives",
                  "Organoheterocyclic compounds", "Nucleosides, nucleotides, and analogues",
                  "Phenylpropanoids and polyketides", "Benzenoids",
                  "Organic oxygen compounds", "Organic nitrogen compounds"),
  Class = c("Fatty Acyls", "Carboxylic acids and derivatives", "Isoquinolines and derivatives",
            "Purine nucleosides", "Acetophenones", "Phenols",
            "Organooxygen compounds", "Organonitrogen compounds")
)
cname_superclass <- bind_rows(
  data.frame(CNAME = paste0("METAB_C18_", sprintf("%03d", 1:50)), Platform = "C18"),
  data.frame(CNAME = paste0("METAB_HILIC_", sprintf("%03d", 1:40)), Platform = "HILIC")
) %>%
  bind_cols(superclass_pairs[sample(nrow(superclass_pairs), 90, replace = TRUE), ])
write.csv(cname_superclass, file.path(root, "out_files", "cname_superclass_lookup.csv"), row.names = FALSE)
write.csv(data.frame(CNAME = paste0("METAB_HILIC_", sprintf("%03d", 1:40)),
                      display_name = paste0("Synthetic metabolite HILIC-", 1:40)),
          file.path(root, "out_files", "hilic_display_names.csv"), row.names = FALSE)

# =============================================================================
# 10. Synthetic stand-in for invitro_data/formated_data.csv (scripts 13, 15, 16)
#    The real in vitro PFOS bacterial growth-screening data belongs to a
#    separate, not-yet-submitted study and cannot be shared here. Only this
#    one file is actually read by any script (confirmed by grepping "invitro_data"
#    across code/); the other real invitro_data files were unused and removed.
# =============================================================================

# Bias toward the fixed 73-species FVA reference list (insilico_data) so
# scripts 15/16's species-matching between the in vitro screening data and the
# metabolic models has realistic overlap (mirroring why the real study
# paired these two datasets in the first place). These are real, public
# species names (not sensitive) used only as labels for fabricated values.
fva_ref_species <- unique(read.csv(
  file.path(root, "insilico_data", "FVA_Secretion_Capacity_HMO_BM_1month_73.csv"),
  check.names = FALSE)$Species)
invitro_species_pool <- pwy_tax_dict[grepl(paste0("^(", paste(gut_genera, collapse = "|"), ") "),
                                         pwy_tax_dict$name), ]
# Script 13's figure specifically needs Bifidobacterium/Lachnospiraceae-genus
# overlap between the cohort species pool and the in vitro strain panel -- boost
# representation of exactly those genera here to make that overlap reliable.
bifido_lachno_genera <- c("Bifidobacterium", "Blautia", "Dorea", "Anaerostipes",
                          "Roseburia", "Coprococcus", "Lachnoclostridium",
                          "Anaerobutyricum", "Enterocloster", "Mediterraneibacter",
                          "Ruminococcus")
bifido_lachno_pool <- invitro_species_pool[grepl(
  paste0("^(", paste(bifido_lachno_genera, collapse = "|"), ") "), invitro_species_pool$name), ]
n_strains <- 175
n_extra <- n_strains - length(fva_ref_species)
n_anchor <- min(nrow(bifido_lachno_pool), round(n_extra * 0.4))
invitro_species_names <- c(
  fva_ref_species,
  sample(bifido_lachno_pool$name, n_anchor),
  sample(invitro_species_pool$name, n_extra - n_anchor, replace = TRUE)
)
invitro_strains <- data.frame(
  StrainID    = sprintf("JEB%05d", sample(10000:99999, n_strains)),
  species_name = sample(invitro_species_names),
  Strain_Info = paste0("DSM ", sample(1000:99999, n_strains))
)
invitro_strains$lineage <- pwy_tax_dict$taxonomic_lineage[match(invitro_strains$species_name, pwy_tax_dict$name)]
invitro_strains$gtdb <- paste0(invitro_strains$lineage, ";s__", invitro_strains$species_name)

drugs <- c("PFOS-high", "PFOS-low", "PFOSNO-high", "PFOSNO-low", "Veh")
invitro_rows <- lapply(seq_len(n_strains), function(i) {
  n_rep <- sample(4:8, 1)
  data.frame(
    StrainID   = invitro_strains$StrainID[i],
    Drug       = sample(drugs, n_rep, replace = TRUE),
    metric     = "AUC",
    date       = format(as.Date("2026-01-01") + sample(0:200, n_rep, replace = TRUE), "%Y%m%d"),
    median_value = round(runif(n_rep, 20, 60), 4),
    perc         = round(rnorm(n_rep, 100, 15), 4)
  ) %>%
    mutate(perc_inhibition = 100 - perc,
           SILVA_Species = invitro_strains$species_name[i],
           Strain_Info = invitro_strains$Strain_Info[i],
           `GTDB_08-rs214_taxonomy` = invitro_strains$gtdb[i])
})
invitro_formated <- bind_rows(invitro_rows)
write.csv(invitro_formated, file.path(root, "invitro_data", "formated_data.csv"), row.names = TRUE)

cat("\nDone. Synthetic example data written to:", file.path(root, "input"), "\n")
cat("Synthetic script-9 stand-in files written to:", file.path(root, "out_files"), "\n")
cat("Synthetic invitro_data stand-in written to:", file.path(root, "invitro_data"), "\n")
