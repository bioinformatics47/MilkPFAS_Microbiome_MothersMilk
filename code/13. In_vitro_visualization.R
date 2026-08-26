# header -----------------------------------------------------------------------
#
# TITLE:   In Vitro PFOS Growth Screening — Cohort Taxa Matching & Figures
#
# PURPOSE: Match cohort species (1m + 6m) to Yuan's in vitro PFOS bacterial
#          growth screening data (AUC-based, % inhibition vs vehicle).
#          Workflow:
#            (1) Extract unique species from cohort species file
#            (2) Annotate with genus/family via taxonomy dictionary
#            (3) Match to Yuan's strain library (Screening_bacteria_Yuan.xlsx)
#            (4) Merge with Yuan's AUC inhibition data (formated_data.csv)
#            (5) Filter to PFOS-low and PFOS-high only (drop PFOSNO)
#            (6) Generate manuscript figures
#
# INPUTS:
#   out_files/species.csv
#     - cohort species file (1m + 6m combined); filter for unique species
#   input/taxonomyDictionary_brack_jan_class_withLineage_bacteriaOnly.tsv
#     - taxonomy dictionary for genus/family annotation
#   invitro_data/Screening_bacteria_Yuan.xlsx
#     - Yuan's strain library with GTDB taxonomy (skip first row)
#   invitro_data/formated_data.csv
#     - Yuan's AUC-based growth inhibition data
#       cols: StrainID, Drug, median_value, perc, perc_inhibition,
#             SILVA_Species, Strain_Info, GTDB_08-rs214_taxonomy
#
#
# DATE: April 2026
# Ellie Holzhausen (EAH) on April 27, 2026
# Haonan Li (HL) on May 25, 2026
#
# set up -----------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(readxl)
library(ggplot2)
library(ggbeeswarm)
library(writexl)


# Set Paths---------------------------------------------------------------------
path_species      <- here::here("out_files/species.csv")
path_yuan_data    <- here::here("invitro_data", "formated_data.csv")
path_out <- here::here("out_figures")
path_out_files <- here::here("out_files")
dir.create(path_out, showWarnings = FALSE, recursive = TRUE)

# step 1: load cohort species and get unique taxon IDs ------------------------

species_raw <- read_csv(path_species, show_col_types = FALSE)

# species.csv is wide format: rows = samples, cols = taxon IDs (X562 etc.)
# extract unique numeric taxon IDs from column headers
taxon_ids <- names(species_raw) %>%
  str_subset("^X\\d+$") %>%
  str_remove("^X") %>%
  as.integer()

# step 2: get species names from taxonomy dictionary --------------------------
tax_dict     <- read_tsv(here::here("input", "taxonomyDictionary_brack_jan_species_withLineage_bacteriaOnly.tsv"))
# cols: name, taxonomy_id, taxonomic_lineage

cohort_species <- tax_dict %>%
  filter(taxonomy_id %in% taxon_ids) %>%
  dplyr::select(taxonomy_id, species_name = name)

message("  Matched: ", nrow(cohort_species), " / ", length(taxon_ids), " taxon IDs")

# step 3: read and clean Yuan data --------------------------------------------
# GTDB taxonomy is more reliable than SILVA_Species (17 SILVA errors found)
# GTDB species names have suffixes (_A, _B, _AP etc.) — strip before matching
# 120/175 strains missing GTDB — fall back to SILVA_Species for those

yuan_data <- read_csv(path_yuan_data, show_col_types = FALSE) %>%
  filter(Drug %in% c("Veh", "PFOS-low", "PFOS-high")) %>%
  mutate(
    Drug         = factor(Drug, levels = c("Veh", "PFOS-low", "PFOS-high")),
    is_floor     = median_value == 1.0,
    species_gtdb = str_extract(`GTDB_08-rs214_taxonomy`, "(?<=s__)[^;]+") %>%
      str_replace_all("_[A-Z0-9]+\\b", "") %>%
      str_trim(),
    species_match = if_else(!is.na(species_gtdb), species_gtdb, SILVA_Species)
  )

# step 3b: remove known outlier floor values for Blautia obeum (JEB00285)
# Two entries with median_value = 1.0 (assay floor) on different dates
# flagged by Yuan as failed runs — remove before analysis
yuan_data <- yuan_data %>%
  filter(!(StrainID == "JEB00285" & median_value == 1.0))

message("  Removed JEB00285 floor outliers. Remaining rows: ", nrow(yuan_data))

length(unique(yuan_data$StrainID)) #175 unique strains


# step 4: match cohort taxa to Yuan data --------------------------------------

yuan_matched <- yuan_data %>%
  filter(species_match %in% cohort_species$species_name)

message("  Strains matched: ",  n_distinct(yuan_matched$StrainID))
message("  Species matched: ",  n_distinct(yuan_matched$species_match))

# Save for supplementary table
write_xlsx(
  yuan_matched,
  path = file.path("out_files", "InVitro_strain_matched.xlsx")
)

# step 5: annotate with genus group and PFAS direction ------------------------

pfas_direction <- tibble(
  genus_group = c("Bifidobacterium", "Streptococcus", "Escherichia",
                  "Shigella", "Veillonella",
                  "Bacteroides", "Phocaeicola", "Parabacteroides"),
  pfas_assoc  = c("Negative", "Negative", "Negative",
                  "Negative", "Negative",
                  "Positive", "Positive", "Positive")
)

lachnospiraceae_genera <- c("Blautia", "Dorea", "Anaerostipes", "Roseburia",
                            "Coprococcus", "Lachnoclostridium", "Anaerobutyricum",
                            "Clostridium", "Enterocloster", "Mediterraneibacter",
                            "Ruminococcus", "Simiaoa")

yuan_annotated <- yuan_matched %>%
  mutate(
    # use species_match (GTDB-preferred) for genus extraction
    Genus = word(species_match, 1),
    genus_group = case_when(
      Genus %in% lachnospiraceae_genera ~ "Lachnospiraceae",
      TRUE                              ~ Genus
    ),
  ) %>%
  left_join(pfas_direction, by = "genus_group") %>%
  mutate(
    pfas_assoc = case_when(
      genus_group == "Lachnospiraceae" ~ "Positive",
      !is.na(pfas_assoc)              ~ pfas_assoc,
      TRUE                            ~ "Unknown"
    )
  )

message("  Rows in annotated dataset: ", nrow(yuan_annotated))
message("  Unique strains: ",            n_distinct(yuan_annotated$StrainID))

# step 6: summarise per strain and group --------------------------------------

message("\nStep 6: Summarizing per strain...")

strain_summary <- yuan_annotated %>%
  filter(Drug != "Veh") %>%
  group_by(StrainID, Strain_Info, Drug, species_match, genus_group, pfas_assoc) %>%
  dplyr::summarise(
    mean_inhibition = mean(perc - 100, na.rm = TRUE),
    n_reps          = dplyr::n(),
    any_floor       = any(is_floor),
    .groups = "drop"
  ) %>%
  filter(pfas_assoc != "Unknown")

# group-level summary
group_summary <- strain_summary %>%
  group_by(genus_group, pfas_assoc, Drug) %>%
  dplyr::summarise(
    mean_inh = mean(mean_inhibition, na.rm = TRUE),
    se_inh   = sd(mean_inhibition, na.rm = TRUE) / sqrt(dplyr::n()),
    n        = dplyr::n(),
    .groups  = "drop"
  )

message("  Groups in final dataset:")
print(dplyr::count(strain_summary, genus_group, pfas_assoc, Drug) %>% as.data.frame())

# summarise from strain level to species level (for per-species bar plots)
species_summary <- strain_summary %>%
  group_by(species_match, genus_group, pfas_assoc, Drug) %>%
  dplyr::summarise(
    mean_inhibition = mean(mean_inhibition, na.rm = TRUE),
    n_strains       = n(),
    any_floor       = any(any_floor),
    .groups = "drop"
  )

# # Figures --------------------------------------------------------------------
# # -- single faceted figure: all groups, strain level -------------------------
# 
# groups_of_interest <- c("Bifidobacterium", "Streptococcus", "Escherichia",
#                         "Bacteroides", "Lachnospiraceae")
# 
# # use strain_summary directly — strain level, not species mean
# plot_data_all <- strain_summary %>%
#   filter(genus_group %in% groups_of_interest) %>%
#   mutate(
#     strain_label = paste0(
#       str_replace(species_match, "^(\\w)\\w+ ", "\\1. "),
#       "\n(", Strain_Info, ")"
#     ),
#     genus_group  = factor(genus_group, levels = groups_of_interest),
#     Drug         = factor(Drug, levels = c("PFOS-low", "PFOS-high"))
#   ) %>%
#   arrange(genus_group, species_match) %>%
#   mutate(strain_label = factor(strain_label, levels = unique(strain_label)))
# 
# n_neg_strains <- plot_data_all %>%
#   filter(pfas_assoc == "Negative") %>%
#   pull(strain_label) %>%
#   unique() %>%
#   length()
# 
# fig_all <- ggplot(
#   plot_data_all,
#   aes(x = strain_label, y = mean_inhibition, fill = Drug)
# ) +
#   geom_hline(yintercept = 0,
#              linetype = "dashed", color = "grey40", linewidth = 0.4) +
#   geom_vline(xintercept = n_neg_strains + 0.5,
#              linetype = "solid", color = "grey60", linewidth = 0.5) +
#   geom_col(
#     position = position_dodge2(preserve = "single", padding = 0.1),
#     width = 0.65, alpha = 0.9
#   ) +
#   scale_fill_manual(
#     values = c("PFOS-low" = "#F59E0B", "PFOS-high" = "#EF4444"),
#     name   = "Condition"
#   ) +
#   scale_y_continuous(
#     breaks = c(-25, 0, 25, 50, 75, 100),
#     labels = function(x) paste0(x, "%"),
#     expand = expansion(mult = c(0.15, 0.1))
#   ) +
#   annotate("text", x = n_neg_strains / 2, y = 105,
#            label = "PFAS negative association",
#            size = 3.2, color = "grey30", fontface = "italic") +
#   annotate("text",
#            x = n_neg_strains + (n_distinct(plot_data_all$strain_label) - n_neg_strains) / 2,
#            y = 105,
#            label = "PFAS positive association",
#            size = 3.2, color = "grey30", fontface = "italic") +
#   labs(x = NULL, y = "% Growth inhibition relative to vehicle") +
#   theme_bw(base_size = 11) +
#   theme(
#     axis.text.x        = element_text(face = "bold.italic", size = 8,
#                                       angle = 40, hjust = 1),
#     axis.text.y        = element_text(size = 10),
#     axis.title.y       = element_text(size = 12, face = "bold"),
#     panel.grid.minor   = element_blank(),
#     panel.grid.major.x = element_blank(),
#     panel.grid.major.y = element_blank(),
#     legend.position    = "top"
#   )
# 
# ggsave(file.path(path_out, "fig_all_groups.pdf"),
#        fig_all, width = 16, height = 5)
# ggsave(file.path(path_out, "fig_all_groups.png"),
#        fig_all, width = 16, height = 5, dpi = 300)
# message("  Saved fig_all_groups")


# -- bifidobacterium and lachnospiraceae: strain level ------------------------
groups_subset <- c("Bifidobacterium", "Lachnospiraceae")

plot_data_subset <- strain_summary %>%
  filter(genus_group %in% groups_subset) %>%
  mutate(
    strain_label = paste0(
      str_replace(species_match, "^(\\w)\\w+ ", "\\1. "),
      "\n(", Strain_Info, ")"
    ),
    genus_group  = factor(genus_group, levels = groups_subset),
    Drug         = factor(
      ifelse(Drug == "PFOS-low", "Low PFOS", "High PFOS"),
      levels = c("Low PFOS", "High PFOS")
    )
  ) %>%
  arrange(genus_group, species_match) %>%
  mutate(strain_label = factor(strain_label, levels = unique(strain_label)))

# identify near-zero bars (tested but negligible response)
plot_data_nearzero <- plot_data_subset %>%
  filter(abs(mean_inhibition) < 0.5)
fig_subset <- ggplot(
  plot_data_subset,
  aes(x = strain_label, y = mean_inhibition, fill = Drug)
) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_col(
    position = position_dodge(width = 0.74, preserve = "single"),
    width = 0.68, color = "black", linewidth = 0.25
  ) +
  geom_text(                              # ← add this block here
    data = plot_data_nearzero,
    aes(label = "~0%", group = Drug),
    position = position_dodge(width = 0.74),
    vjust = -0.5,
    size  = 4.5,
    color = "black",
    fontface = "italic"
  ) +
  facet_grid(~ genus_group, scales = "free_x", space = "free_x") +
  facet_grid(~ genus_group, scales = "free_x", space = "free_x") +
  scale_fill_manual(
    values = c("Low PFOS" = "#2166AC", "High PFOS" = "#B2182B"),
    name   = NULL
  ) +
  scale_y_continuous(
    limits = c(-100, 25),
    breaks = seq(-100, 25, 25),
    labels = function(x) paste0(x, "%")
  ) +
  labs(x = NULL, y = "% Growth change relative to vehicle") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x        = element_text(face = "bold.italic", size = 7.3,
                                      angle = 40, hjust = 1, color = "black"),
    axis.text.y        = element_text(size = 10, face = "bold", color = "black"),
    axis.title.y       = element_text(size = 12, face = "bold"),
    axis.line          = element_line(color = "black", linewidth = 0.7),
    axis.ticks         = element_line(color = "black"),
    strip.text         = element_text(size = 12, face = "bold.italic"),
    strip.background   = element_blank(),
    panel.spacing      = unit(1.2, "lines"),
    legend.position    = c(0.5, 0.12),
    legend.direction   = "horizontal",
    legend.background  = element_rect(fill = "white", color = "grey80"),
    legend.text        = element_text(size = 9.5)
  )

ggsave(file.path(path_out, "fig_bifido_lach_strains.pdf"),
       fig_subset, width = 11, height = 5)
ggsave(file.path(path_out, "fig_bifido_lach_strains.png"),
       fig_subset, width = 11, height = 5, dpi = 300)
message("  Saved fig_bifido_lach")

# ── Matching summary statistics ───────────────────────────────────────────────
# Total strains in Yuan dataset (before cohort matching)
cat("Total strains in Yuan screening dataset:    ", n_distinct(yuan_data$StrainID), "\n")
cat("Total species in Yuan screening dataset:    ", n_distinct(yuan_data$species_match), "\n")

# After cohort matching (yuan_matched)
cat("\nStrains matched to cohort species:          ", n_distinct(yuan_matched$StrainID), "\n")
cat("Species matched to cohort species:          ", n_distinct(yuan_matched$species_match), "\n")

# After annotation and Unknown removal (strain_summary)
cat("\nStrains retained after filtering Unknown:   ", n_distinct(strain_summary$StrainID), "\n")
cat("Species retained after filtering Unknown:   ", n_distinct(strain_summary$species_match), "\n")

# Breakdown by genus group
cat("\nUnique strains per genus group (after filtering):\n")
strain_summary %>%
  group_by(genus_group, pfas_assoc) %>%
  summarise(n_strains  = n_distinct(StrainID),
            n_species  = n_distinct(species_match),
            .groups = "drop") %>%
  arrange(pfas_assoc, genus_group) %>%
  print(n = Inf)

# Figure 7 subset only (Bifido + Lachno)
cat("\nFigure 7 subset (Bifidobacterium + Lachnospiraceae only):\n")
strain_summary %>%
  filter(genus_group %in% c("Bifidobacterium", "Lachnospiraceae")) %>%
  group_by(genus_group) %>%
  summarise(n_strains = n_distinct(StrainID),
            n_species = n_distinct(species_match),
            .groups = "drop") %>%
  print()

# Data availability per group (how many strains had high vs low dose)
cat("\nData availability — strains with each dose condition:\n")
strain_summary %>%
  group_by(genus_group, Drug) %>%
  summarise(n_strains = n_distinct(StrainID), .groups = "drop") %>%
  pivot_wider(names_from = Drug, values_from = n_strains, values_fill = 0) %>%
  print()

cat("=============================================\n")

write.csv(cohort_species,
          here::here('out_files', 'cohort_species.csv'),
          row.names = FALSE)

