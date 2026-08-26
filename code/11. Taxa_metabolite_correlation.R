# Title: TAXA–METABOLITE SPEARMAN CORRELATION HEATMAP
#
# PURPOSE: Compute spearman correlation between significant taxa (p<0.05) from pfas-taxa
# analysis with significant metabolites (p<0.05) from pfas-metabolite analysis
#
# DATE:    March 2026
#
# - Taxa selected from continuous PFAS-taxa CLR regression (p < 0.05, species)
# - Metabolites selected from MWAS (p < 0.05, C18 + HILIC combined)
# - Spearman correlation between per-sample CLR values and metabolite intensities
# - Bottom annotation: horizontal Spearman colour bar + PFAS direction legend
# - Top annotation: PFAS bars showing beta direction per taxon (empty PFAS removed)
# - Left annotation: PFAS bars showing beta direction per metabolite
# - Duplicate metabolites across platforms suffixed _C18 / _HILIC
# - Brackets [] removed from taxa names for correct alphabetical sorting
# - Separate cell/font size settings per scenario
#
#  Code Review
#  Ellie Holzhausen (EAH) on April 28, 2026
#  Haonan Li (HL) on May 25, 2026
#
#-------------------------------------------------------------------------------
# Set-up
rm(list = ls())

library(tidyverse)
library(dplyr)
library(here)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(gridExtra)

# Load metabolite column lists--------------------------------------------------
keep_c18   <- readRDS(here::here("out_files", "metabolite_cols_c18.rds"))
keep_hilic <- readRDS(here::here("out_files", "metabolite_cols_hilic.rds"))

display_c18   <- read.csv(here::here("out_files", "c18_display_names.csv"),
                          stringsAsFactors = FALSE)
display_hilic <- read.csv(here::here("out_files", "hilic_display_names.csv"),
                          stringsAsFactors = FALSE)
# PFAS variable names
pfas_vars_continuous <- c("PFBS", "PFHxS", "PFNA", "PFOA", "PFOS")
pfas_display <- c("PFBS"  = "PFBS",
                  "PFHxS" = "PFHxS",
                  "PFNA"  = "PFNA",
                  "PFOA"  = "PFOA",
                  "PFOS"  = "PFOS")


# Get significant taxa----------------------------------------------------------
get_sig_taxa_species <- function(taxa_results_file) {
  
  df <- read.csv(here::here("out_files", taxa_results_file)) %>%
    filter(taxa_level == "Species",
           predictor %in% pfas_vars_continuous,
           p_val < 0.05,
           !is.na(p_val),
           !is.na(betas)) %>%
    mutate(
      plot_beta = dplyr::coalesce(as.numeric(Beta_IQR), as.numeric(betas)),
      direction = ifelse(plot_beta > 0, "Positive", "Negative")
    )
  
  taxa_summary <- df %>%
    filter(!is.na(p_val)) %>%
    group_by(taxonomy_id, name) %>%
    summarise(best_pval = min(p_val, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(best_pval)) %>%
    mutate(
      clr_col = paste0(taxonomy_id, "_CLR"),
      # Remove brackets from display names e.g. [Clostridium] -> Clostridium
      name    = gsub("^\\[|\\]", "", name)
    )
  
  taxa_pfas_dir <- df %>%
    dplyr::select(taxonomy_id, name, predictor, direction) %>%
    distinct() %>%
    mutate(name = gsub("^\\[|\\]", "", name))
  
  list(summary = taxa_summary, pfas_dir = taxa_pfas_dir)
}


# Get significant metabolites--------------------------------------------------
# removing PFBS associated metabolites because there are no any PFBS associated taxa
get_sig_metabolites <- function(mwas_c18_file, mwas_hilic_file,
                                p_thresh = 0.05) {
  
  sig_c18 <- read.csv(here::here("out_files", mwas_c18_file)) %>%
    filter(P_value < p_thresh, 
           Metabolite %in% keep_c18,
           PFAS != "PFBS_pgmL") %>%
    pull(Metabolite) %>% unique()
  
  sig_hilic <- read.csv(here::here("out_files", mwas_hilic_file)) %>%
    filter(P_value < p_thresh, 
           Metabolite %in% keep_hilic,
           PFAS != "PFBS_pgmL") %>%
    pull(Metabolite) %>% unique()
  
  both       <- intersect(sig_c18, sig_hilic)
  c18_only   <- setdiff(sig_c18,   both)
  hilic_only <- setdiff(sig_hilic, both)
  
  met_labels <- c(
    c18_only,
    hilic_only,
    paste0(both, "_C18"),
    paste0(both, "_HILIC")
  )
  
  # Manual deduplication for known name mismatches across platforms
  name_aliases <- list(
    c("NICOTINAMIDE",  "NICOTINAMIDE(B3)"),
    c("ASPARAGINE",    "L-ASPARAGINE"),
    c("LAURIC ACID",   "FA12:0(LAURATE")
  )
  
  for (alias_pair in name_aliases) {
    in_c18   <- alias_pair[alias_pair %in% sig_c18]
    in_hilic <- alias_pair[alias_pair %in% sig_hilic]
    if (length(in_c18) > 0 && length(in_hilic) > 0) {
      c18_only   <- setdiff(c18_only,   in_c18)
      hilic_only <- setdiff(hilic_only, in_hilic)
      extra_labels <- c(paste0(in_c18,   "_C18"),
                        paste0(in_hilic, "_HILIC"))
      met_labels   <- c(met_labels, extra_labels)
      cat("  Alias dedup:", paste(alias_pair, collapse = " / "),
          "-> suffixed as C18/HILIC\n")
    }
  }
  
  cat("  Sig metabolites — C18:", length(sig_c18),
      "| HILIC:", length(sig_hilic),
      "| Combined:", length(met_labels), "\n")
  
  return(list(
    met_labels = met_labels,
    both       = both,
    c18_only   = c18_only,
    hilic_only = hilic_only
  ))
}


# Build per-sample matrix------------------------------------------------------
build_cor_matrix <- function(clr_file, metab_c18_file, metab_hilic_file,
                             sig_taxa, sig_mets, timepoint = "1m") {
  
  clr_df <- read.csv(here::here("out_files", clr_file), check.names = FALSE)
  
  clr_cols_avail <- intersect(sig_taxa$summary$clr_col, colnames(clr_df))
  clr_sub <- clr_df %>%
    dplyr::select(merge_id_dyad, all_of(clr_cols_avail))
  
  cat("  CLR taxa available:", length(clr_cols_avail),
      "of", nrow(sig_taxa$summary), "significant\n")
  
  metab_c18   <- read.csv(here::here("out_files", metab_c18_file),
                          check.names = FALSE)
  metab_hilic <- read.csv(here::here("out_files", metab_hilic_file),
                          check.names = FALSE)
  
  # Fix merge_id_dyad for 6m metabolite files
  if (timepoint == "6m") {
    # The 6m metabolite files store merge_id_dyad
    # ending in -01 (1m PFAS timepoint used during merging in data cleaning).
    # Recoding to -06 here is necessary to match the CLR file which uses
    # 6m microbiome timepoint IDs. Verified: no duplicates created after
    # recoding (0 duplicates across 102 samples), and 97/102 IDs successfully
    # match the 6m CLR file (5 non-matching = participants with metabolomics
    # but no microbiome data after QC filters).
    metab_c18 <- metab_c18 %>%
      mutate(merge_id_dyad = str_replace(merge_id_dyad, "-01$", "-06"))
    metab_hilic <- metab_hilic %>%
      mutate(merge_id_dyad = str_replace(merge_id_dyad, "-01$", "-06"))
    cat("  Recoded metabolite merge_id_dyad: -01 -> -06\n")
    cat("  Duplicates after recode:", sum(duplicated(metab_c18$merge_id_dyad)), "\n")
  }
  
  c18_cols   <- intersect(c(sig_mets$c18_only, sig_mets$both), keep_c18)
  hilic_cols <- intersect(c(sig_mets$hilic_only, sig_mets$both), keep_hilic)
  
  metab_c18_sub <- metab_c18 %>%
    dplyr::select(merge_id_dyad, any_of(c18_cols))
  
  metab_hilic_sub <- metab_hilic %>%
    dplyr::select(merge_id_dyad, any_of(hilic_cols))
  
  if (length(sig_mets$both) > 0) {
    metab_c18_sub <- metab_c18_sub %>%
      rename_with(~ paste0(., "_C18"),   any_of(sig_mets$both))
    metab_hilic_sub <- metab_hilic_sub %>%
      rename_with(~ paste0(., "_HILIC"), any_of(sig_mets$both))
  }
  
  merged <- clr_sub %>%
    inner_join(metab_c18_sub,   by = "merge_id_dyad") %>%
    inner_join(metab_hilic_sub, by = "merge_id_dyad")
  
  cat("  Merged samples:", nrow(merged), "\n")
  return(merged)
}

# Compute Spearman correlations------------------------------------------------
compute_spearman <- function(merged_df, sig_taxa, sig_mets) {
  
  clr_cols <- sig_taxa$summary$clr_col[
    sig_taxa$summary$clr_col %in% colnames(merged_df)]
  
  met_cols <- sig_mets$met_labels[sig_mets$met_labels %in% colnames(merged_df)]
  
  cat("  Metabolite cols matched in merged:", length(met_cols), "\n")
  cat("  CLR cols matched in merged:", length(clr_cols), "\n")
  cat("  Computing", length(clr_cols), "taxa x",
      length(met_cols), "metabolites =",
      length(clr_cols) * length(met_cols), "pairs\n")
  
  results_list <- list()
  idx <- 1L
  for (cc in clr_cols) {
    for (mc in met_cols) {
      x <- merged_df[[cc]]
      y <- merged_df[[mc]]
      ok <- sum(is.finite(x) & is.finite(y))
      if (ok < 3) {
        rho_val <- NA_real_; p_val <- NA_real_
      } else {
        ct      <- cor.test(x, y, method = "spearman", exact = FALSE)
        rho_val <- as.numeric(ct$estimate)
        p_val   <- as.numeric(ct$p.value)
      }
      results_list[[idx]] <- data.frame(
        clr_col   = cc,
        met_label = mc,
        rho       = rho_val,
        p_val     = p_val,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  
  results <- bind_rows(results_list)
  results$FDR <- p.adjust(results$p_val, method = "BH")
  results <- results %>%
    left_join(sig_taxa$summary %>% dplyr::select(clr_col, name),
              by = "clr_col")
  
  cat("  Significant pairs (FDR < 0.05):", sum(results$FDR < 0.05), "\n")
  cat("  Significant pairs (p < 0.05):",   sum(results$p_val < 0.05), "\n")
  
  return(results)
}


# Build PFAS annotation for taxa columns---------------------------------------
build_pfas_annotation <- function(sig_taxa, taxa_in_plot) {
  
  annot <- sig_taxa$pfas_dir %>%
    filter(name %in% taxa_in_plot) %>%
    mutate(predictor_label = pfas_display[predictor]) %>%
    dplyr::select(name, predictor_label, direction) %>%
    distinct() %>%
    pivot_wider(names_from  = predictor_label,
                values_from = direction,
                values_fill = NA_character_) %>%
    column_to_rownames("name")
  
  for (p in pfas_display) {
    if (!p %in% colnames(annot)) annot[[p]] <- NA_character_
  }
  
  annot <- annot[, as.character(pfas_display), drop = FALSE]
  return(annot)
}


# Build heatmap gtable (shared between view and save)--------------------------
# Returns pheatmap object — does NOT open a device
build_heatmap <- function(cor_results, sig_taxa,
                          mwas_c18_file, mwas_hilic_file,
                          scenario_label, tag,
                          min_sig_taxa  = 3,
                          min_sig_met   = 3,
                          p_threshold   = 0.05,
                          fdr_threshold = 0.10,
                          cellwidth     = 10,
                          cellheight    = 8.5,
                          fontsize_row  = 9,
                          fontsize_col  = 9,
                          bar_left      = 0.4,
                          bar_y         = 0.35) {
  # Build mwas_combined early — needed for direct PFAS filter and sorting
  mwas_c18_df   <- read.csv(here::here("out_files", mwas_c18_file))
  mwas_hilic_df <- read.csv(here::here("out_files", mwas_hilic_file))
  
  mwas_combined <- bind_rows(mwas_c18_df, mwas_hilic_df) %>%
    filter(P_value < p_threshold, !is.na(Estimate),
           PFAS != "PFBS_pgmL") %>%
    mutate(
      plot_est   = dplyr::coalesce(as.numeric(Beta_IQR), as.numeric(Estimate)),
      direction  = ifelse(plot_est > 0, "Positive", "Negative"),
      PFAS_clean = gsub("_pgmL$", "", PFAS)
    ) %>%
    filter(PFAS_clean %in% pfas_vars_continuous)
  
  # Iterative mutual filter — keeps looping until both conditions are
  # simultaneously satisfied:
  #   - every taxon has >= min_sig_taxa significant metabolite correlations
  #   - every metabolite has >= min_sig_met significant taxon correlations
  sig_pairs <- cor_results %>% filter(FDR < fdr_threshold)
  
  taxa_keep <- unique(sig_pairs$name)
  met_keep  <- unique(sig_pairs$met_label)
  
  repeat {
    taxa_prev <- taxa_keep
    met_prev  <- met_keep
    
    # Keep metabolites correlated with enough taxa
    met_keep <- sig_pairs %>%
      filter(name %in% taxa_keep) %>%
      count(met_label) %>%
      filter(n >= min_sig_met) %>%
      pull(met_label)
    
    # Keep taxa correlated with enough metabolites
    taxa_keep <- sig_pairs %>%
      filter(met_label %in% met_keep) %>%
      count(name) %>%
      filter(n >= min_sig_taxa) %>%
      pull(name)
    
    # Stop when nothing changes
    if (setequal(taxa_keep, taxa_prev) && setequal(met_keep, met_prev)) break
  }
  
  cat("  After iterative mutual filter (taxa >=", min_sig_taxa,
      "mets, mets >=", min_sig_met, "taxa):",
      length(taxa_keep), "taxa x", length(met_keep), "metabolites\n")
  
  # Final iterative filter combining direct PFAS requirement and thresholds
  # Guarantees every taxon has >= min_sig_taxa metabolites AND
  # every metabolite has >= min_sig_met taxa AND has a direct PFAS association
  repeat {
    taxa_prev <- taxa_keep
    met_prev  <- met_keep
    
    # Restrict to metabolites with direct PFAS association
    met_keep <- met_keep[met_keep %in% mwas_combined$Metabolite]
    
    # Re-apply metabolite threshold
    met_keep <- sig_pairs %>%
      filter(name %in% taxa_keep, met_label %in% met_keep) %>%
      count(met_label) %>%
      filter(n >= min_sig_met) %>%
      pull(met_label)
    
    # Keep only those with direct PFAS association
    met_keep <- met_keep[met_keep %in% mwas_combined$Metabolite]
    
    # Re-apply taxa threshold
    taxa_keep <- sig_pairs %>%
      filter(met_label %in% met_keep) %>%
      count(name) %>%
      filter(n >= min_sig_taxa) %>%
      pull(name)
    
    if (setequal(taxa_keep, taxa_prev) && setequal(met_keep, met_prev)) break
  }
  
  cat("  After final filter (direct PFAS + thresholds):",
      length(taxa_keep), "taxa x", length(met_keep), "metabolites\n")
  
  # Update mwas_combined to reflect final met_keep
  mwas_combined <- mwas_combined %>%
    filter(Metabolite %in% met_keep)
  
  if (length(taxa_keep) == 0 | length(met_keep) == 0) {
    cat("  Not enough significant pairs — skipping.\n")
    return(invisible(NULL))
  }
  
  cor_dedup <- cor_results %>%
    filter(name %in% taxa_keep, met_label %in% met_keep) %>%
    group_by(name, met_label) %>%
    slice_min(order_by = p_val, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # Sort taxa: positive group first then negative, alphabetically within each
  # Brackets [] removed from name for correct alphabetical sorting
  taxa_direction <- sig_taxa$pfas_dir %>%
    filter(name %in% taxa_keep) %>%
    count(name, direction) %>%
    group_by(name) %>%
    slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::select(name, dominant_dir = direction)
  
  taxa_order <- sig_taxa$summary %>%
    filter(name %in% taxa_keep) %>%
    left_join(taxa_direction, by = "name") %>%
    mutate(
      dominant_dir = replace_na(dominant_dir, "Positive"),
      # Remove brackets for sorting only — display name unchanged
      name_sort    = gsub("^\\[|\\]", "", name)
    ) %>%
    arrange(dominant_dir, name_sort) %>%
    pull(name)
  taxa_order <- c(taxa_order, setdiff(taxa_keep, taxa_order))
  
  # Sort metabolites: positive group first then negative, by max abs rho within
  met_direction <- mwas_combined %>%
    filter(Metabolite %in% met_keep) %>%
    count(Metabolite, direction) %>%
    group_by(Metabolite) %>%
    slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::rename(met_label    = Metabolite,
                  dominant_dir = direction)
  
  met_order <- cor_dedup %>%
    group_by(met_label) %>%
    summarise(max_rho = max(abs(rho), na.rm = TRUE), .groups = "drop") %>%
    left_join(met_direction, by = "met_label") %>%
    mutate(dominant_dir = replace_na(dominant_dir, "Positive")) %>%
    arrange(dominant_dir, desc(max_rho)) %>%
    pull(met_label)
  
  # Build rho matrix
  rho_wide <- cor_dedup %>%
    dplyr::select(name, met_label, rho) %>%
    pivot_wider(names_from = name, values_from = rho) %>%
    column_to_rownames("met_label") %>%
    as.matrix()
  
  taxa_order <- taxa_order[taxa_order %in% colnames(rho_wide)]
  met_order  <- met_order[met_order   %in% rownames(rho_wide)]
  rho_wide   <- rho_wide[met_order, taxa_order, drop = FALSE]
  
  # Significance matrix using FDR
  sig_wide <- cor_dedup %>%
    dplyr::select(name, met_label, FDR) %>%
    pivot_wider(names_from = name, values_from = FDR) %>%
    column_to_rownames("met_label") %>%
    as.matrix()
  sig_wide <- sig_wide[met_order, taxa_order, drop = FALSE]
  
  sig_labels <- matrix(
    ifelse(!is.na(sig_wide) & sig_wide < fdr_threshold, "*", ""),
    nrow = nrow(sig_wide), ncol = ncol(sig_wide),
    dimnames = dimnames(sig_wide)
  )
  
  # Taxa annotation — remove all-NA PFAS columns
  pfas_annot_col <- build_pfas_annotation(sig_taxa, taxa_order)
  pfas_annot_col <- pfas_annot_col[colnames(rho_wide), , drop = FALSE]
  pfas_annot_col <- pfas_annot_col[, colSums(!is.na(pfas_annot_col)) > 0,
                                   drop = FALSE]
  
  # Metabolite annotation
  annot_row <- mwas_combined %>%
    dplyr::select(Metabolite, PFAS_clean, direction) %>%
    distinct() %>%
    pivot_wider(names_from  = PFAS_clean,
                values_from = direction,
                values_fill = NA_character_) %>%
    column_to_rownames("Metabolite")
  
  for (p in pfas_vars_continuous) {
    if (!p %in% colnames(annot_row)) annot_row[[p]] <- NA_character_
  }
  annot_row <- annot_row[, pfas_vars_continuous, drop = FALSE]
  annot_row <- annot_row[rownames(rho_wide), , drop = FALSE]
  annot_row <- annot_row[, colSums(!is.na(annot_row)) > 0, drop = FALSE]
  
  # Manual display name overrides for plot labels
  # Applied after display lookup — covers co-eluters, long names, and duplicates
  manual_display_overrides <- c(
    # Co-eluters — combined short labels
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
    "GAMMA-LINOLENIC ACID"                      = "Linolenic Acid"
  )
  
  # Simple title-case cleaner for names not in manual overrides
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
  
  # Apply display names to row labels
  row_display_map <- setNames(rownames(rho_wide), rownames(rho_wide))
  
  for (nm in rownames(rho_wide)) {
    # Determine base name and suffix
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
    
    # Check manual override first
    if (base %in% names(manual_display_overrides)) {
      row_display_map[nm] <- paste0(manual_display_overrides[base], suffix)
    } else {
      # Fall back to display lookup then basic clean
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
  
  rownames(rho_wide)   <- row_display_map[rownames(rho_wide)]
  rownames(sig_wide)   <- row_display_map[rownames(sig_wide)]
  rownames(sig_labels) <- row_display_map[rownames(sig_labels)]
  rownames(annot_row)  <- row_display_map[rownames(annot_row)]
  
  cat("  Display names substituted:",
      sum(row_display_map != names(row_display_map)), "metabolites\n")
  cat("  Final row labels:\n")
  print(rownames(rho_wide))
  
  # Unified annotation colors
  dir_colors    <- c("Positive" = "purple3", "Negative" = "orange3")
  all_pfas_cols <- unique(c(colnames(pfas_annot_col), colnames(annot_row)))
  annot_colors  <- setNames(
    lapply(all_pfas_cols, function(p) dir_colors),
    all_pfas_cols
  )
  
  pal <- colorRampPalette(brewer.pal(11, "RdBu"))(100)
  
  # Build pheatmap silently — no device opened here
  p <- pheatmap(
    mat                  = rho_wide,
    color                = pal,
    breaks               = seq(-0.7, 0.7, length.out = 101),
    display_numbers      = sig_labels,
    fontsize_number      = 7,
    number_color         = "black",
    annotation_col       = pfas_annot_col,
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
    main                 = scenario_label,
    silent               = TRUE
  )
  
  # Make taxa (column) names italic
  gt <- p$gtable
  col_idx <- which(gt$layout$name == "col_names")
  if (length(col_idx) > 0) {
    gt$grobs[[col_idx]]$gp$fontface <- "italic"
  }
  p$gtable <- gt
  
  # Canvas sizing
  max_col_chars <- max(nchar(colnames(rho_wide)))
  max_row_chars <- max(nchar(rownames(rho_wide)))
  
  col_label_h <- max_col_chars * fontsize_col * 0.55 * sin(pi / 4) / 72 + 0.2
  
  # legend_h_in: inches reserved at bottom for horizontal legend strip
  legend_h_in <- 1.3
  
  fig_h <- max(5.0, nrow(rho_wide) * cellheight / 72 + col_label_h + 1.5 + legend_h_in)
  
  annot_row_w   <- ncol(annot_row) * cellwidth / 72 + 0.3
  row_label_w   <- max_row_chars * fontsize_row * 0.55 / 72 + 0.3
  fig_w         <- max(7.0, ncol(rho_wide) * cellwidth / 72 + annot_row_w + row_label_w + 0.55)
  
  return(list(p = p, fig_w = fig_w, fig_h = fig_h, legend_h_in = legend_h_in,
              bar_left = bar_left, bar_y = bar_y))
}


# Draw heatmap to a device (screen or file)------------------------------------
draw_heatmap <- function(heatmap_obj, legend_y = NULL) {
  
  if (is.null(heatmap_obj)) return(invisible(NULL))
  
  p           <- heatmap_obj$p
  legend_h_in <- heatmap_obj$legend_h_in
  
  grid.newpage()
  
  # Heatmap: full width, everything above the bottom legend strip
  pushViewport(viewport(
    x      = unit(0, "npc"),
    y      = unit(legend_h_in, "inches"),
    width  = unit(1, "npc"),
    height = unit(1, "npc") - unit(legend_h_in, "inches"),
    just   = c("left", "bottom")
  ))
  grid.draw(p$gtable)
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
  
  # ── Horizontal Spearman colour bar ────────────────────────────────────────
  bar_left  <- heatmap_obj$bar_left
  bar_right <- bar_left + 2.8
  bar_y     <- heatmap_obj$bar_y
  bar_h     <- 0.20
  each_w    <- (bar_right - bar_left) / n
  
  # Title above bar
  grid.text("Spearman Correlation",
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
  grid.text("-0.7", x = unit(bar_left, "inches"),
            y = unit(bar_y - bar_h / 2 - 0.08, "inches"),
            gp = gpar(fontsize = 9), just = c("center", "top"))
  grid.text("0",    x = unit((bar_left + bar_right) / 2, "inches"),
            y = unit(bar_y - bar_h / 2 - 0.08, "inches"),
            gp = gpar(fontsize = 9), just = c("center", "top"))
  grid.text("0.7",  x = unit(bar_right, "inches"),
            y = unit(bar_y - bar_h / 2 - 0.08, "inches"),
            gp = gpar(fontsize = 9), just = c("center", "top"))
  
  # ── Direction legend (right of colour bar) ────────────────────────────────
  dir_left <- bar_right + 0.8
  
  grid.text("PFAS-Taxa/Metabolite Direction",
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
            gp   = gpar(fontsize = 10),
            just = c("left", "center"))
  
  grid.rect(x      = unit(dir_left + 1.3, "inches"),
            y      = unit(bar_y + 0.02, "inches"),
            width  = bw, height = bh,
            just   = c("left", "center"),
            gp     = gpar(fill = "orange3", col = NA))
  grid.text("Negative",
            x    = unit(dir_left + 1.60, "inches"),
            y    = unit(bar_y + 0.02, "inches"),
            gp   = gpar(fontsize = 10),
            just = c("left", "center"))
  
  popViewport()
}


# Save heatmap to PDF----------------------------------------------------------
save_heatmap <- function(heatmap_obj, tag) {
  
  if (is.null(heatmap_obj)) {
    cat("  Nothing to save — heatmap object is NULL\n")
    return(invisible(NULL))
  }
  
  outfile <- here::here("out_figures",
                        paste0("heatmap_taxa_metab_", tag, ".pdf"))
  
  pdf(outfile,
      width  = heatmap_obj$fig_w,
      height = heatmap_obj$fig_h)
  
  draw_heatmap(heatmap_obj)
  
  dev.off()
  cat("  Saved:", outfile, "\n")
}


# Run 2 scenarios--------------------------------------------------------------
scenarios <- list(
  list(
    label            = NA,
    tag              = "1m",
    timepoint        = "1m",
    taxa_file        = "COMBINED_CLR_continuous_1m_1m.csv",
    clr_file         = "CLR_perSample_species_1m.csv",
    mwas_c18_file    = "MWAS_continuous_c18_1m.csv",
    mwas_hilic_file  = "MWAS_continuous_hilic_1m.csv",
    metab_c18_file   = "PFAS1m_c18_1m.csv",
    metab_hilic_file = "PFAS1m_hilic_1m.csv",
    min_sig_taxa     = 2,
    min_sig_met      = 2,
    p_threshold      = 0.05,
    fdr_threshold    = 0.05,
    cellwidth        = 12,
    cellheight       = 10,
    fontsize_row     = 10,
    fontsize_col     = 10,
    bar_left         = 1.5,
    bar_y            = 1.30
  ),
  list(
    label            = NA,
    tag              = "6m",
    timepoint        = "6m",
    taxa_file        = "COMBINED_CLR_continuous_1m_6m.csv",
    clr_file         = "CLR_perSample_species_6m.csv",
    mwas_c18_file    = "MWAS_continuous_c18_6m.csv",
    mwas_hilic_file  = "MWAS_continuous_hilic_6m.csv",
    metab_c18_file   = "PFAS1m_c18_6m.csv",
    metab_hilic_file = "PFAS1m_hilic_6m.csv",
    min_sig_taxa     = 5,
    min_sig_met      = 5,
    p_threshold      = 0.05,
    fdr_threshold    = 0.05,
    cellwidth        = 12,
    cellheight       = 9,
    fontsize_row     = 8,
    fontsize_col     = 8,
    bar_left         = 4.0,
    bar_y            = 1.30
  )
)


# STEP 1 — Build and view all plots interactively------------------------------
heatmap_objects <- list()

for (sc in scenarios) {
  
  cat("\n========================================\n")
  cat("Processing:", sc$label, "\n")
  cat("========================================\n")
  
  sig_taxa <- get_sig_taxa_species(sc$taxa_file)
  cat("  Significant taxa (p < 0.05):", nrow(sig_taxa$summary), "\n")
  
  sig_mets <- get_sig_metabolites(sc$mwas_c18_file, sc$mwas_hilic_file,
                                  p_thresh = sc$p_threshold)
  
  merged <- build_cor_matrix(
    sc$clr_file, sc$metab_c18_file, sc$metab_hilic_file,
    sig_taxa, sig_mets, timepoint = sc$timepoint
  )
  
  cor_results <- compute_spearman(merged, sig_taxa, sig_mets)
  
  write.csv(cor_results,
            here::here("out_files",
                       paste0("spearman_taxa_metab_", sc$tag, ".csv")),
            row.names = FALSE)
  
  hm <- build_heatmap(
    cor_results     = cor_results,
    sig_taxa        = sig_taxa,
    mwas_c18_file   = sc$mwas_c18_file,
    mwas_hilic_file = sc$mwas_hilic_file,
    scenario_label  = sc$label,
    tag             = sc$tag,
    min_sig_taxa    = sc$min_sig_taxa,
    min_sig_met     = sc$min_sig_met,
    p_threshold     = sc$p_threshold,
    fdr_threshold   = sc$fdr_threshold,
    cellwidth       = sc$cellwidth,
    cellheight      = sc$cellheight,
    fontsize_row    = sc$fontsize_row,
    fontsize_col    = sc$fontsize_col,
    bar_left        = sc$bar_left,
    bar_y           = sc$bar_y
  )
  
  heatmap_objects[[sc$tag]] <- hm
  
  cat("  Displaying:", sc$label, "\n")
  draw_heatmap(hm)
}


# STEP 2 — Save plots after reviewing------------------------------------------
for (sc in scenarios) {
  save_heatmap(heatmap_objects[[sc$tag]], sc$tag)
}
