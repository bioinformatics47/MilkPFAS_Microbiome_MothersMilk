# TITLE:   1. Data Cleaning.R
#
# PURPOSE: Read and clean the meta data and microbiome data
#
# DATE:    January 12, 2024
# CODE REVIEW:
# Ellie Holzhausen (EAH) on April 21, 2026
# Haonan Li (HL) on May 25, 2026

# set up -----------------------------------------------------------------------

#clear the workspace
rm(list = ls())

#load libraries
library(plyr); library(tidyr); library(reshape)
library(purrr); library(stringr); library(lme4)
library(lmerTest); library(corrplot); library(lubridate)
library(pscl); library(dplyr); 
library(naniar); library(tibble); library(tableone); 
library(readxl); library(here); library(readr)
library(tidyverse)
library(ggplot2)

# read in the meta data and IDs
meta <- read.csv(here::here("input", "mothersMilk_metadata_timepointsAsRows_updated051022_Temporary24mDiet.csv"))
IDs <- read_excel(here::here("input", "GORAN_MICROBIOME_MANIFEST 10-2022.xlsx"))
# meta <- read.csv(here::here("input", "mothersMilk_metadata_timepointsAsRows_updated051022_Temporary24mDiet.csv"))
# IDs <- read_excel(here::here("input", "GORAN_MICROBIOME_MANIFEST 10-2022.xlsx"))

# Clean column names
colnames(IDs) <- IDs[1,]
#drop the first row
IDs <- IDs[-1,]

# Pull dyad_id and timepoint from `Old Together` column
IDs$dyad_id <- substr(IDs$`Old Together`, 9, 11)
IDs$timepoint <- substr(IDs$`Old Together`, 1, 2)

# We only want the 1 and 6 months data, so select irrelevant indices to drop
IDs_toremove <- IDs$`CORE ID`[which(IDs$timepoint %in% c("12", "18", "24", "36", "37"))]
IDs_toremove <- paste0("MG", IDs_toremove)
# there are 567 samples to remove from timepoints 12,18,24,36,37 timepoints

# read in the WGS count data
species <- read.csv(here::here("input", "counts_bracken_species.csv"))
genus <- read.csv(here::here("input", "counts_bracken_genus.csv"))
family <- read.csv(here::here("input", "counts_bracken_family.csv"))
order <- read.csv(here::here("input", "counts_bracken_order.csv"))
class <- read.csv(here::here("input", "counts_bracken_class.csv"))
phylum <- read.csv(here::here("input", "counts_bracken_phylum.csv"))


# Load Bracken relative abundance files (used for CLR analysis)
# Filtering decisions made on counts data above; relative abundances
# are subsetted to passing samples and taxa after all filtering steps
ra_species <- read.csv(here::here("input", "relative_abundance_bracken_species.csv"))
ra_genus   <- read.csv(here::here("input", "relative_abundance_bracken_genus.csv"))
ra_family  <- read.csv(here::here("input", "relative_abundance_bracken_family.csv"))
ra_order   <- read.csv(here::here("input", "relative_abundance_bracken_order.csv"))
ra_class   <- read.csv(here::here("input", "relative_abundance_bracken_class.csv"))
ra_phylum  <- read.csv(here::here("input", "relative_abundance_bracken_phylum.csv"))

# Transpose relative abundance files to match counts orientation
# rows = MG samples, cols = taxa (same structure as counts after transposing)
transpose_ra <- function(ra_df) {
  row.names(ra_df) <- paste0("X", ra_df$SampleID)
  ra_df <- dplyr::select(ra_df, -SampleID)
  data.frame(t(ra_df))
}

ra_species <- transpose_ra(ra_species)
ra_genus   <- transpose_ra(ra_genus)
ra_family  <- transpose_ra(ra_family)
ra_order   <- transpose_ra(ra_order)
ra_class   <- transpose_ra(ra_class)
ra_phylum  <- transpose_ra(ra_phylum)

# # Transpose counts data frames (same structure as RA files-so same function can be used)
# species <- transpose_ra(species)
# genus   <- transpose_ra(genus)
# family  <- transpose_ra(family)
# order   <- transpose_ra(order)
# class   <- transpose_ra(class)
# phylum  <- transpose_ra(phylum)

# need to transpose each of the data frames
row.names(species) <- paste0("X", species$SampleID)
species <- dplyr::select(species, -c("SampleID"))
species <- data.frame(t(species))

row.names(genus) <- paste0("X", genus$SampleID)
genus <- dplyr::select(genus, -c("SampleID"))
genus <- data.frame(t(genus))

row.names(family) <- paste0("X", family$SampleID)
family <- dplyr::select(family, -c("SampleID"))
family <- data.frame(t(family))

row.names(order) <- paste0("X", order$SampleID)
order <- dplyr::select(order, -c("SampleID"))
order <- data.frame(t(order))

row.names(class) <- paste0("X", class$SampleID)
class <- dplyr::select(class, -c("SampleID"))
class <- data.frame(t(class))

row.names(phylum) <- paste0("X", phylum$SampleID)
phylum <- dplyr::select(phylum, -c("SampleID"))
phylum <- data.frame(t(phylum))

# get the taxonomy data downloaded
species_tax <- read_tsv(here::here("input", "taxonomyDictionary_brack_jan_species_withLineage_bacteriaOnly.tsv"))
genus_tax <- read_tsv(here::here("input", "taxonomyDictionary_brack_jan_genus_withLineage_bacteriaOnly.tsv"))
family_tax <- read_tsv(here::here("input", "taxonomyDictionary_brack_jan_family_withLineage_bacteriaOnly.tsv"))
order_tax <- read_tsv(here::here("input", "taxonomyDictionary_brack_jan_order_withLineage_bacteriaOnly.tsv"))
class_tax <- read_tsv(here::here("input", "taxonomyDictionary_brack_jan_class_withLineage_bacteriaOnly.tsv"))
phylum_tax <- read_tsv(here::here("input", "taxonomyDictionary_brack_jan_phylum_withLineage_bacteriaOnly.tsv"))


# function to drop rare taxa ####
# Dropping rare taxa in metagenomics is common as it helps to reduce noise,
# improve computational efficiency and interpretation
## set threshold = 50 #### (meaning only taxa with total abundance, sum of counts across all samples, of 50 or more are only retained)
threshold = 50

drop_rare_taxa <- function(data, threshold){
  prop_lt_thresh <- colMeans(data < threshold) # Gives the proportion of data < threshold for each column
  cols_to_drop <- names(prop_lt_thresh[prop_lt_thresh >= 0.25]) #Find the species with >25% below-threshold data
  data[, !names(prop_lt_thresh) %in% cols_to_drop] #Drop species if > 25% of data has below threshold value)
}

# clean microbiome data --------------------------------------------------------

# drop the control samples (anything that is not "MG" is control)
species <- species[which(substr(row.names(species),1,2) == "MG"),]
genus <- genus[which(substr(row.names(genus), 1, 2) == "MG"),]
family <- family[which(substr(row.names(family), 1, 2) == "MG"),]
order <- order[which(substr(row.names(order), 1, 2) == "MG"),]
class <- class[which(substr(row.names(class), 1, 2) == "MG"),]
phylum <- phylum[which(substr(row.names(phylum), 1, 2) == "MG"),]


# drop the 12,18,24,36,37 month timepoints in the count data (keeping 1 and 6 months)
species <- species[-which(row.names(species) %in% IDs_toremove),]
genus <- genus[-which(row.names(genus) %in% IDs_toremove),]
family <- family[-which(row.names(family) %in% IDs_toremove),]
order <- order[-which(row.names(order) %in% IDs_toremove),]
class <- class[-which(row.names(class) %in% IDs_toremove),]
phylum <- phylum[-which(row.names(phylum) %in% IDs_toremove),]
# 358 observations after only keeping 1 and 6 months timepoint

### Avg reads / range per sample
aa <- rowSums(species)
bb <- rowSums(genus)
cc <- rowSums(family)
dd <- rowSums(order)
ee <- rowSums(class)
ff <- rowSums(phylum)
mean(c(aa,bb,cc,dd,ee,ff))
range(c(aa,bb,cc,dd,ee,ff))

# remove from all levels any observation with fewer than 1,000,000 total counts
# This removes noise and improve statistical robustness
species <- species[rowSums(species) >= 1000000, ]
genus   <- genus[rowSums(genus)   >= 1000000, ]
family  <- family[rowSums(family)  >= 1000000, ]
order   <- order[rowSums(order)   >= 1000000, ]
class   <- class[rowSums(class)   >= 1000000, ]
phylum  <- phylum[rowSums(phylum)  >= 1000000, ]
# dropped 1 observations with < 1 million reads, 357 observations remaining

# Checking again: avg reads / range per sample
aa <- rowSums(species)
bb <- rowSums(genus)
cc <- rowSums(family)
dd <- rowSums(order)
ee <- rowSums(class)
ff <- rowSums(phylum)
mean(c(aa,bb,cc,dd,ee,ff))
range(c(aa,bb,cc,dd,ee,ff))

## ==========================
## QC: Threshold exploration
## ==========================
# QC: species-only elbow & distributions
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(ggplot2); library(scales)
})

# --- helper to coerce to numeric matrix safely
as_num_mat <- function(x) { x <- as.data.frame(x); data.matrix(x) }

count_retained <- function(mat, threshold, prev_cut){
  mat <- as_num_mat(mat)
  prop_lt <- colMeans(mat < threshold, na.rm = TRUE)
  sum(prop_lt < prev_cut)
}

build_elbow_df <- function(mat, label,
                           thresholds = c(5,10,25,50,75,100,200,500,1000),
                           prev_grid  = c(0.10,0.25,0.50)){
  tidyr::expand_grid(threshold = thresholds, prev_cut = prev_grid) %>%
    mutate(n_taxa = purrr::map2_int(threshold, prev_cut, ~count_retained(mat, .x, .y)),
           level = label)
}

# ---- Elbow for species (combined 1m + 6m)
elbow_species <- build_elbow_df(species, "Species (1m + 6m)")

p_elbow <- ggplot(elbow_species, aes(x = threshold, y = n_taxa, 
                          color = factor(prev_cut), group = prev_cut)) +
  geom_line() +
  geom_point(size = 1.2) +
  geom_vline(xintercept = 50, linetype = "dashed") +   # chosen threshold
  scale_color_brewer(palette = "Set1",
                     name = "Prevalence cut\n(prop < threshold)") +
  scale_x_continuous(trans = "log10",
                     breaks = c(5, 10, 25, 50, 100, 200, 500, 1000)) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Count threshold (log10)",
    y = "Retained species",
    caption = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )

print(p_elbow)

# Save as PDF in out_figures/ folder
ggsave(
  filename = here::here("out_figures", "elbow_plot_species.pdf"),
  plot = p_elbow,
  width = 7,
  height = 5,
  units = "in"
)

## ==========================
## QC: Threshold exploration — ALL taxonomic levels
## ==========================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(ggplot2); library(scales)
})

as_num_mat <- function(x) { x <- as.data.frame(x); data.matrix(x) }

count_retained <- function(mat, threshold, prev_cut){
  mat <- as_num_mat(mat)
  prop_lt <- colMeans(mat < threshold, na.rm = TRUE)
  sum(prop_lt < prev_cut)
}

build_elbow_df <- function(mat, label,
                           thresholds = c(5,10,25,50,75,100,200,500,1000),
                           prev_grid  = c(0.10,0.25,0.50)){
  tidyr::expand_grid(threshold = thresholds, prev_cut = prev_grid) %>%
    mutate(n_taxa = purrr::map2_int(threshold, prev_cut, ~count_retained(mat, .x, .y)),
           level = label)
}

# Build elbow data for all taxonomic levels
elbow_all <- bind_rows(
  build_elbow_df(species, "Species"),
  build_elbow_df(genus,   "Genus"),
  build_elbow_df(family,  "Family"),
  build_elbow_df(order,   "Order"),
  build_elbow_df(class,   "Class"),
  build_elbow_df(phylum,  "Phylum")
)

# Set factor order so panels go from finest to coarsest resolution
elbow_all$level <- factor(elbow_all$level,
                          levels = c("Species","Genus","Family","Order","Class","Phylum"))

p_elbow_all <- ggplot(elbow_all, aes(x = threshold, y = n_taxa,
                                     color = factor(prev_cut), group = prev_cut)) +
  geom_line() +
  geom_point(size = 1.2) +
  geom_vline(xintercept = 50, linetype = "dashed") +
  facet_wrap(~ level, scales = "free_y", ncol = 3) +
  scale_color_brewer(palette = "Set1",
                     name = "Prevalence cut\n(prop < threshold)") +
  scale_x_continuous(trans = "log10",
                     breaks = c(5, 10, 25, 50, 100, 200, 500, 1000)) +
  labs(x = "Count threshold (log10)",
       y = "Retained taxa") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

print(p_elbow_all)

ggsave(
  filename = here::here("out_figures", "elbow_plot_all_levels.pdf"),
  plot = p_elbow_all,
  width = 12, height = 8, units = "in"
)

# remove rare taxa from the count data (helps reduce noise and focus on relevant taxa)
# NOTE: applying same filtering threshold at all taxa levels based on species level elbow plot is not technically correct. We have to look elbow plot at all levels and threshold.
# For this project, we will be using only species data and threshold is set based on elbow plot for species level
# species
species_norare <- drop_rare_taxa(species, threshold)
# genus
genus_norare <- drop_rare_taxa(genus, threshold)
#family
family_norare <- drop_rare_taxa(family, threshold)
#order
order_norare <- drop_rare_taxa(order, threshold)
#class
class_norare <- drop_rare_taxa(class, threshold)
#phylum
phylum_norare <- drop_rare_taxa(phylum, threshold)

# Make a table of the number of observations before and after filtering for each 
# tax. level (not included in manuscript, may add if reviewers ask)

stab <- data.frame(matrix(nrow = 6, ncol = 2))
colnames(stab) <- c("Before", "After")
row.names(stab) <- c("Species", "Genus", "Family", "Order", "Class", "Phylum")

# populate the table
stab["Species", "Before"] <- ncol(species)
stab["Species", "After"] <- ncol(species_norare)

stab["Genus", "Before"] <- ncol(genus)
stab["Genus", "After"] <- ncol(genus_norare)

stab["Family", "Before"] <- ncol(family)
stab["Family", "After"] <- ncol(family_norare)

stab["Order", "Before"] <- ncol(order)
stab["Order", "After"] <- ncol(order_norare)

stab["Class", "Before"] <- ncol(class)
stab["Class", "After"] <- ncol(class_norare)

stab["Phylum", "Before"] <- ncol(phylum)
stab["Phylum", "After"] <- ncol(phylum_norare)

# End table

# prep the ID data frame to merge
IDs$merge_id <- paste0("MG", IDs$`CORE ID`)

# make sure the order of the IDs data frame is the same as metagenome_t
IDs <- IDs[order(IDs$merge_id),]

# check if there are any IDs that aren't in the 1 and 6 month metagenome data
IDs$merge_id[!IDs$merge_id %in% row.names(species)] # IDs not found in metagenome data
length(IDs$merge_id[!IDs$merge_id %in% row.names(species)]) # 568 not found
# 567 belong to the IDs dropped from the removed timepoints and 1 removed due to low total counts

# make merge_id_dyad to add into metagenome_t
IDs$dyad_id2 <- str_pad(IDs$dyad_id, width = 4, pad = 0)
IDs$merge_id_dyad <- paste0("MM-", IDs$dyad_id2, "-", IDs$timepoint)

# Subset just the merge id dyad information
IDs_temp <- dplyr::select(IDs, c("merge_id_dyad", "merge_id"))

# Add merge id to microbiome data
species_norare$merge_id <- row.names(species_norare)
genus_norare$merge_id <- row.names(genus_norare)
family_norare$merge_id <- row.names(family_norare)
order_norare$merge_id <- row.names(order_norare)
class_norare$merge_id <- row.names(class_norare)
phylum_norare$merge_id <- row.names(phylum_norare)

# Add merge_id to RA files at same step as counts objects
# (rownames are MG sample IDs — same as counts after transposing)
ra_species$merge_id <- rownames(ra_species)
ra_genus$merge_id   <- rownames(ra_genus)
ra_family$merge_id  <- rownames(ra_family)
ra_order$merge_id   <- rownames(ra_order)
ra_class$merge_id   <- rownames(ra_class)
ra_phylum$merge_id  <- rownames(ra_phylum)

# # add the raw counts into the count data frame
# species_norare$counts <- rowSums(species)
# genus_norare$counts <- rowSums(genus)
# family_norare$counts <- rowSums(family)
# order_norare$counts <- rowSums(order)
# class_norare$counts <- rowSums(class)
# phylum_norare$counts <- rowSums(phylum)

# # Integrate the merge_id_dyad column
# species <- left_join(species_norare, IDs_temp, by = "merge_id")
# genus <- left_join(genus_norare, IDs_temp, by = "merge_id")
# family <- left_join(family_norare, IDs_temp, by = "merge_id")
# order <- left_join(order_norare, IDs_temp, by = "merge_id")
# class <- left_join(class_norare, IDs_temp, by = "merge_id")
# phylum <- left_join(phylum_norare, IDs_temp, by = "merge_id")
# 
# # check how many unique IDs we have for the 157 samples from 1 & 6 month microbiome
# unique_dyad_ids <- unique(substr(species$merge_id_dyad, 4, 7))
# length(unique_dyad_ids) #208
# 
# # Change rownames to indicate dyad IDs (6 months)
# row.names(species) <- species$merge_id_dyad
# row.names(genus) <- genus$merge_id_dyad
# row.names(family) <- family$merge_id_dyad
# row.names(order) <- order$merge_id_dyad
# row.names(class) <- class$merge_id_dyad
# row.names(phylum) <- phylum$merge_id_dyad
# 
# # Drop unnecessary columns (keep raw counts)
# species <- dplyr::select(species, -c("merge_id", "merge_id_dyad"))
# genus <- dplyr::select(genus, -c("merge_id", "merge_id_dyad"))
# family <- dplyr::select(family, -c("merge_id", "merge_id_dyad"))
# order <- dplyr::select(order, -c("merge_id", "merge_id_dyad"))
# class <- dplyr::select(class, -c("merge_id", "merge_id_dyad"))
# phylum <- dplyr::select(phylum, -c("merge_id", "merge_id_dyad"))
# 
# # set any observations which are 3 SDs above the mean to missing (outlier truncation)
# # It is done to identify and handle extreme outliers in multiple taxonomic data frames
# # don't check the last observation, since that is the total counts before filtering 
# for(i in 1:(ncol(species)-1)){
#   col_mean <- mean(species[,i], na.rm = T)
#   col_sd <- sd(species[,i], na.rm = T)
#   species[species[,i] > (col_mean + 3*col_sd), i] <- round(col_mean + 3*col_sd,0)
# }
# 
# for(i in 1:(ncol(genus)-1)){
#   col_mean <- mean(genus[,i], na.rm = T)
#   col_sd <- sd(genus[,i], na.rm = T)
#   genus[genus[,i] > (col_mean + 3*col_sd), i] <- round(col_mean + 3*col_sd,0)
# }
# 
# for(i in 1:(ncol(family)-1)){
#   col_mean <- mean(family[,i], na.rm = T)
#   col_sd <- sd(family[,i], na.rm = T)
#   family[family[,i] > (col_mean + 3*col_sd), i] <- round(col_mean + 3*col_sd,0)
# }
# 
# for(i in 1:(ncol(order)-1)){
#   col_mean <- mean(order[,i], na.rm = T)
#   col_sd <- sd(order[,i], na.rm = T)
#   order[order[,i] > (col_mean + 3*col_sd), i] <- round(col_mean + 3*col_sd,0)
# }
# 
# for(i in 1:(ncol(class)-1)){
#   col_mean <- mean(class[,i], na.rm = T)
#   col_sd <- sd(class[,i], na.rm = T)
#   class[class[,i] > (col_mean + 3*col_sd), i] <- round(col_mean + 3*col_sd,0)
# }
# 
# for(i in 1:(ncol(phylum)-1)){
#  col_mean <- mean(phylum[,i], na.rm = T)
#  col_sd <- sd(phylum[,i], na.rm = T)
#  phylum[phylum[,i] > (col_mean + 3*col_sd), i] <- round(col_mean + 3*col_sd,0)
# }

# Use Bracken relative abundance data (subsetted to passing samples/taxa)
# instead of counts data for downstream CLR analysis.
# merge_id_dyad is added via IDs_temp using the MG sample ID as key.

# add merge_id_dyad rownames to RA subset
add_merge_id_dyad <- function(ra_sub, IDs_temp) {
  ra_sub$merge_id <- rownames(ra_sub)
  df <- left_join(ra_sub, IDs_temp, by = "merge_id")
  rownames(df) <- df$merge_id_dyad
  df <- dplyr::select(df, -merge_id, -merge_id_dyad)
  return(df)
}

# Helper: subset RA file to passing samples and taxa from counts filtering
subset_bracken_ra <- function(ra_df, passing_mg_ids, passing_taxa_ids) {
  ra_samples <- intersect(passing_mg_ids, rownames(ra_df))
  cat("  Passing samples found in RA file:", length(ra_samples),
      "of", length(passing_mg_ids), "\n")
  ra_taxa <- intersect(passing_taxa_ids, colnames(ra_df))
  cat("  Passing taxa found in RA file:", length(ra_taxa),
      "of", length(passing_taxa_ids), "\n")
  ra_df[ra_samples, ra_taxa, drop = FALSE]
}

# Get passing sample IDs from counts-filtered _norare objects
# (these reflect read depth filter + rare taxa filter + non-bacterial removal)
passing_mg_ids_1m6m <- rownames(species_norare)
cat("Total passing MG samples (1m + 6m):", length(passing_mg_ids_1m6m), "\n")

# Subset each RA level — taxa IDs from each level's own _norare counts object
cat("\n--- Subsetting species RA ---\n")
ra_species_sub <- subset_bracken_ra(ra_species, passing_mg_ids_1m6m,
                                    colnames(species_norare))
cat("\n--- Subsetting genus RA ---\n")
ra_genus_sub   <- subset_bracken_ra(ra_genus,   passing_mg_ids_1m6m,
                                    colnames(genus_norare))
cat("\n--- Subsetting family RA ---\n")
ra_family_sub  <- subset_bracken_ra(ra_family,  passing_mg_ids_1m6m,
                                    colnames(family_norare))
cat("\n--- Subsetting order RA ---\n")
ra_order_sub   <- subset_bracken_ra(ra_order,   passing_mg_ids_1m6m,
                                    colnames(order_norare))
cat("\n--- Subsetting class RA ---\n")
ra_class_sub   <- subset_bracken_ra(ra_class,   passing_mg_ids_1m6m,
                                    colnames(class_norare))
cat("\n--- Subsetting phylum RA ---\n")
ra_phylum_sub  <- subset_bracken_ra(ra_phylum,  passing_mg_ids_1m6m,
                                    colnames(phylum_norare))

# Attach merge_id_dyad to RA subsets
# Overwrites species/genus etc. with RA data — rownames are now merge_id_dyad
species <- add_merge_id_dyad(ra_species_sub, IDs_temp)
genus   <- add_merge_id_dyad(ra_genus_sub,   IDs_temp)
family  <- add_merge_id_dyad(ra_family_sub,  IDs_temp)
order   <- add_merge_id_dyad(ra_order_sub,   IDs_temp)
class   <- add_merge_id_dyad(ra_class_sub,   IDs_temp)
phylum  <- add_merge_id_dyad(ra_phylum_sub,  IDs_temp)

# Verify
cat("Species samples after merge:", nrow(species), "\n")
cat("Species taxa after merge:",    ncol(species), "\n")
cat("Unique dyad IDs:", length(unique(substr(rownames(species), 4, 7))), "\n")

# NOTE: Winsorization (3SD truncation) is NOT applied to relative abundance data.
# The original winsorization was designed for raw counts. For Bracken relative
# abundances going into CLR, values are already bounded [0,1] and the 3SD
# threshold from count scale is not meaningful. CLR transformation handles
# compositional scale differences across samples.


# clean meta data ####
# Keep only IDs of 01 and 06 months
meta_breastfed <- meta[grep("-(01|06)$", meta$merge_id_dyad), ]
#View(meta_breastfed)

# Number of breastfeedings
# Recode breastmilk_per_day from categorical codes to actual feedings per day
# Code 0 = Never (0), Code 1 = less than daily (use 0.5 as approximation),
# Codes 2-9 = code - 1 feedings, Code 10 = 8+ (use 8 as conservative floor)

meta_breastfed <- meta_breastfed %>%
  mutate(breastmilk_per_day = case_when(
    breastmilk_per_day == 0  ~ 0,
    breastmilk_per_day == 1  ~ 0.5,   # less than once/day — adjust if preferred
    breastmilk_per_day == 10 ~ 8,     # 8 or more — conservative floor
    !is.na(breastmilk_per_day) ~ breastmilk_per_day - 1  # codes 2–9
  ))

sum(meta_breastfed$breastmilk_per_day == 0 & meta_breastfed$timepoint == 1, na.rm = TRUE) #7
sum(is.na(meta_breastfed$breastmilk_per_day) & meta_breastfed$timepoint == 1) #9
# At 1 month, 7 IDs have 0 and 9 IDs have NAs reported in breastmilk_per_day but we have PFAS data on them

# Check what IDs are these
# IDs that have 0's at 1 month
meta_breastfed$merge_id_dyad[
  meta_breastfed$timepoint == 1 &
    !is.na(meta_breastfed$breastmilk_per_day) &
    meta_breastfed$breastmilk_per_day == 0
]

# IDs that have NAs at 1 month
meta_breastfed$merge_id_dyad[
  is.na(meta_breastfed$breastmilk_per_day) & meta_breastfed$timepoint == 1
]

# So need to replace values on these IDs because these IDs lack data on 24hr diet recall as well
meta_breastfed$breastmilk_per_day[
  meta_breastfed$merge_id_dyad %in% c(
    # 7 with 0s
    "MM-0073-01","MM-0099-01","MM-0153-01","MM-0178-01","MM-0180-01","MM-0185-01","MM-0201-01",
    # 9 with NAs
    "MM-0019-01","MM-0020-01","MM-0026-01","MM-0027-01","MM-0033-01","MM-0203-01",
    "MM-0218-01","MM-0220-01","MM-0222-01"
  ) &
    meta_breastfed$timepoint == 1 &
    (is.na(meta_breastfed$breastmilk_per_day) | meta_breastfed$breastmilk_per_day == 0)
] <- median(
  meta_breastfed$breastmilk_per_day[
    meta_breastfed$timepoint == 1 &
      !is.na(meta_breastfed$breastmilk_per_day) &
      meta_breastfed$breastmilk_per_day > 0
  ],
  na.rm = TRUE
)


# Check again
sum((is.na(meta_breastfed$breastmilk_per_day) | meta_breastfed$breastmilk_per_day == 0) &
      meta_breastfed$timepoint == 1)

# Work on imputing for 6 months
sum(is.na(meta_breastfed$breastmilk_per_day) & meta_breastfed$timepoint == 6) #8
sum(meta_breastfed$breastmilk_per_day == 0 & meta_breastfed$timepoint == 6, na.rm = TRUE) #72

# Pull breastmilk_per_day from 24 hr diet recall for IDs that are zero for breastmilk/day but has PFAS value (ONLY 22 out of 72)
bm_missed <- read_excel(here::here("input", "missing_bmpercent/bmpercent_missing.xlsx")) %>%
  na.omit() %>% 
  dplyr::select(dyad_id, timepoint, breastmilk_per_day, formula_per_day) %>%
  mutate(
    dyad_id = as.character(dyad_id),  # keep dyad_id as character
    timepoint = as.integer(timepoint),
    merge_id_dyad = paste0("MM-", str_pad(dyad_id, width = 4, pad = "0"), "-", str_pad(timepoint, width = 2, pad = "0"))
  )

# Paste these values to meta_breastfed
meta_breastfed <- meta_breastfed %>%
  dplyr::left_join(bm_missed %>% dplyr::select(merge_id_dyad, breastmilk_per_day),
                   by = "merge_id_dyad",
                   suffix = c("_orig", "_impute")) %>%
  dplyr::mutate(
    breastmilk_per_day = ifelse(
      timepoint == 6 &
        breastmilk_per_day_orig == 0 &                # only overwrite 0’s
        !is.na(breastmilk_per_day_impute),            # only if we have recall values
      breastmilk_per_day_impute,                      # use recall value
      breastmilk_per_day_orig                         # otherwise keep original
    )
  ) %>%
  dplyr::select(-c(breastmilk_per_day_orig, breastmilk_per_day_impute))

# Check again
sum(is.na(meta_breastfed$breastmilk_per_day) & meta_breastfed$timepoint == 6) #8
sum(meta_breastfed$breastmilk_per_day == 0 & meta_breastfed$timepoint == 6, na.rm = TRUE) #50

# IDs that have 0's at 6 month
meta_breastfed$merge_id_dyad[
  meta_breastfed$timepoint == 6 &
    !is.na(meta_breastfed$breastmilk_per_day) &
    meta_breastfed$breastmilk_per_day == 0
]

# IDs that have NAs at 6 month
meta_breastfed$merge_id_dyad[
  is.na(meta_breastfed$breastmilk_per_day) & meta_breastfed$timepoint == 6
]

# Impute median values for NAs at 6 months (These IDs do not have diet recall records. Likely left study)
meta_breastfed$breastmilk_per_day[
  meta_breastfed$merge_id_dyad %in% c(
    "MM-0043-06","MM-0046-06","MM-0053-06","MM-0106-06",
    "MM-0118-06","MM-0144-06","MM-0152-06","MM-0183-06"
  ) &
    meta_breastfed$timepoint == 6 &
    is.na(meta_breastfed$breastmilk_per_day)
] <- median(
  meta_breastfed$breastmilk_per_day[
    meta_breastfed$timepoint == 6 &
      !is.na(meta_breastfed$breastmilk_per_day) &
      meta_breastfed$breastmilk_per_day > 0
  ],
  na.rm = TRUE
)

# check again
sum(is.na(meta_breastfed$breastmilk_per_day) & meta_breastfed$timepoint == 6)
sum(meta_breastfed$breastmilk_per_day == 0 & meta_breastfed$timepoint == 6, na.rm = TRUE) #50 (For these 50, these could be true zeros as infant start solid foods by 6 months)

summary(meta_breastfed$breastmilk_per_day)

# Make BMI categories
meta_breastfed <- meta_breastfed %>%
  mutate(
    prepreg_bmi_cat = case_when(
      prepreg_bmi_kgm2 < 25              ~ "Normal",
      prepreg_bmi_kgm2 >= 25 & prepreg_bmi_kgm2 < 30 ~ "Overweight",
      prepreg_bmi_kgm2 >= 30              ~ "Obesity",
      TRUE                                ~ NA_character_
    ),
    prepreg_bmi_cat = factor(prepreg_bmi_cat, levels = c("Normal", "Overweight", "Obesity"))
  )

## Make gestational age categories
meta_breastfed <- meta_breastfed %>%
  mutate(
    gestational_age_cat = case_when(
      gestational_age_category == "<38" ~ "Early",
      gestational_age_category %in% c("38-40", "40") ~ "Ontime",
      gestational_age_category %in% c("40-42", ">40") ~ "Late",
      TRUE ~ NA_character_
    ),
    gestational_age_cat = factor(
      gestational_age_cat,
      levels = c("Early", "Ontime", "Late")
    )
  )

# Populate age of solid foods where it is NA
# Convert the age_solid_foods_clean to numeric 
meta_breastfed$age_solid_foods_clean <- as.numeric(meta_breastfed$age_solid_foods_clean) # Warning here about NAs is fine. NAs are due to non-numeric values like 2-May and solid food names
# for 2 IDs. 2-May is an error in the file because the baby's birthdate is 3rd may.

# Create a variable indicating at what month solids foods were introduced
meta_breastfed <- meta_breastfed %>%
  dplyr::group_by(dyad_id) %>% # Group by participant
  dplyr::mutate(age_of_solid_foods = ifelse(is.na(age_solid_foods_clean), # Populate NAs with available age of solid foods value
                                     dplyr::first(na.omit(age_solid_foods_clean)),
                                     age_solid_foods_clean))

# mode_of_delivery_cat
meta_breastfed$mode_of_delivery_cat <-factor(ifelse(meta_breastfed$mode_of_delivery == 1,
                                                    "Vaginal", ifelse(meta_breastfed$mode_of_delivery == 2, "C-Section", NA)),
                                             levels = c("Vaginal","C-Section"))
summary(meta_breastfed$mode_of_delivery_cat)
# Vaginal C-Section      NA's 
#     307       102        12 

# Antibiotics
table(meta_breastfed$baby_antibiotics)
#   1   2   3   4 
# 369  34   4   2  

meta_breastfed$baby_antibiotics[meta_breastfed$baby_antibiotics == 1] <- "No"
meta_breastfed$baby_antibiotics[meta_breastfed$baby_antibiotics %in% c(2,3)] <- "Yes"
meta_breastfed$baby_antibiotics[meta_breastfed$baby_antibiotics == 4] <- NA
table(meta_breastfed$baby_antibiotics)
#  No Yes
# 369  38

table(meta_breastfed$mother_antibiotics)
#   1   2   3   4 
# 358  31   8  12 

meta_breastfed$mother_antibiotics[meta_breastfed$mother_antibiotics == 1] <- "No"
meta_breastfed$mother_antibiotics[meta_breastfed$mother_antibiotics %in% c(2,3)] <- "Yes"
meta_breastfed$mother_antibiotics[meta_breastfed$mother_antibiotics == 4] <- NA
table(meta_breastfed$mother_antibiotics)
#  No Yes
# 358  39

# Make Baby gender to categorical
meta_breastfed <- meta_breastfed %>%
  mutate(baby_gender = recode(baby_gender, `1` = "Female", `2` = "Male"))

table(meta_breastfed$baby_gender)
# Female   Male 
# 223    197 

meta_trim <- dplyr::select(meta_breastfed, c("merge_id_dyad", "timepoint", "dyad_id",
                                             "age_in_days", "baby_gender", "mother_age",
                                             "mode_of_delivery_cat", "baby_birthweight_kg", "prepreg_bmi_cat", "prepreg_wt_kg",
                                             "baby_antibiotics","mother_antibiotics",
                                             "gestational_age_cat", "breastmilk_per_day",
                                             "SES_index_final","age_of_solid_foods")) 

# Save meta_trim file
write.csv(meta_trim, here::here("out_files", "meta_trim.csv"))

# check that all the taxa have labels in the tax files
print(colnames(species)[-which(colnames(species) %in% paste0("X", species_tax$taxonomy_id))]) # 9606 is homo sapiens, remove
print(colnames(genus)[-which(colnames(genus) %in% paste0("X", genus_tax$taxonomy_id))]) # 9605 is genus homo, remove
print(colnames(family)[-which(colnames(family) %in% paste0("X", family_tax$taxonomy_id))]) # 9604 is Hominidae, great apes, remove
print(colnames(order)[-which(colnames(order) %in% paste0("X", order_tax$taxonomy_id))]) # 9443 is primates, remove
print(colnames(class)[-which(colnames(class) %in% paste0("X", class_tax$taxonomy_id))]) # 183963 belong to Archaea, 2731619 belong to virus, 40674 is homosapiens, remove all
print(colnames(phylum)[-which(colnames(phylum) %in% paste0("X", phylum_tax$taxonomy_id))]) # 2731618 belong to virus, 28890 belong to Archaea, 7711 belong to chordata, remove all

# remove irrelevant taxa
species <- species[, !colnames(species) %in% c("X9606")]
genus   <- genus[, !colnames(genus) %in% c("X9605")]
family  <- family[, !colnames(family) %in% c("X9604")]
order   <- order[, !colnames(order) %in% c("X9443")]
class   <- class[, !colnames(class) %in% c("X183963",  "X2731619", "X40674")]
phylum  <- phylum[, !colnames(phylum) %in% c("X2731618", "X28890", "X7711")]

# save the microbiome outputs
write.csv(species, here::here("out_files", "species.csv"))
write.csv(genus, here::here("out_files", "genus.csv"))
write.csv(family, here::here("out_files", "family.csv"))
write.csv(order, here::here("out_files", "order.csv"))
write.csv(class, here::here("out_files", "class.csv"))
write.csv(phylum, here::here("out_files", "phylum.csv"))

# Load quantified PFAS data ----------------------------------------------------
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


# Filter only needed variables
PFAS <- PFAS %>%
  dplyr::select(dyad_id, timepoint, contains("_pgmL"))

# These MDL values are coming from quantification report provided by Doug walker lab (Document is located here: JHSPH-Shares/EHE/Restricted/ECLIPSE/Lab Projects/PFAS R01/Reprocessed PFAS Data June 2025/CLU0031-Study_Report_LC_Summary_Report.docx)
mdl_values <- c(
  "N.EtFOSAA_pgmL" = 1.12, "N.MeFOSAA_pgmL" = 0.74, "PFBA_pgmL" = 3.33,
  "PFDA_pgmL" = 1.02, "PFDoA_pgmL" = 0.74, "PFHpA_pgmL"  = 1.16,
  "PFHxA_pgmL" = 3.93, "PFNA_pgmL" = 0.66, "PFOA_pgmL"  = 2.28,
  "PFPeA_pgmL"  = 4.47, "PFTeDA_pgmL" = 0.41, "PFUnA_pgmL" = 1.16,
  "PFBS_pgmL" = 7.00, "PFDoS_pgmL" = 0.83, "PFHps_pgmL" = 1.30,
  "PFHxS_pgmL" = 0.62, "PFNS_pgmL"  = 1.44, "PFTrDA_pgmL" = 0.83,
  "PFOS_pgmL" = 1.09, "PFPeAS_pgmL"    = 0.69
)

pfas_cols <- intersect(names(mdl_values), colnames(PFAS))

PFAS_long <- PFAS %>%
  dplyr::select(timepoint, all_of(pfas_cols)) %>%
  pivot_longer(
    cols = all_of(pfas_cols),
    names_to = "PFAS_names",
    values_to = "value"
  ) %>%
  mutate(MDL = mdl_values[PFAS_names])

# percent below MDL by timepoint
percent_below_mdl_by_tp <- PFAS_long %>%
  dplyr::group_by(PFAS_names, timepoint) %>%
  dplyr::summarize(
    n = sum(!is.na(value)),
    percent_below_MDL = mean(value < MDL, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  dplyr::mutate(percent_above_MDL = 100 - percent_below_MDL)


# Make a plot showing percent of samples below MDL------------------------------
# Convert to factor and ensure labels match the values in your scale_fill_manual
percent_below_mdl_by_tp$timepoint <- factor(percent_below_mdl_by_tp$timepoint, 
                                            levels = c(1, 6), 
                                            labels = c("1 Month", "6 Month"))

# Keep only 1 month
percent_below_mdl_by_tp_1m <- percent_below_mdl_by_tp %>%
  dplyr::filter(timepoint == "1 Month")

# Order by average percent below MDL across timepoints (descending)
pfas_order <- percent_below_mdl_by_tp_1m %>%
  dplyr::group_by(PFAS_names) %>%
  dplyr::summarize(mean_below = mean(percent_below_MDL), .groups = "drop") %>%
  dplyr::arrange((mean_below)) %>%
  dplyr::pull(PFAS_names)

percent_below_mdl_by_tp_1m <- percent_below_mdl_by_tp_1m %>%
  mutate(PFAS_names = factor(PFAS_names, levels = pfas_order))

ggplot(percent_below_mdl_by_tp_1m, aes(x = PFAS_names, y = percent_above_MDL, fill = timepoint)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), width = 0.5) +
  geom_hline(yintercept = 75, linetype = "dotted", color = "black", linewidth = 0.8) +
  geom_hline(yintercept = 25, linetype = "dotted", color = "black", linewidth = 0.8) +
  scale_fill_manual(values = c("1 Month" = "#2166AC")) +
  scale_x_discrete(labels = function(x) {
    x <- gsub("_pgmL", "", x)
    x <- gsub("N\\.MeFOSAA", "N-MeFOSAA", x)
    x <- gsub("N\\.EtFOSAA", "N-EtFOSAA", x)
    x <- gsub("PFHps", "PFHpS", x)
    x
  }) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(
    x = NULL,
    y = "% Samples Above MDL",
    fill = "Timepoint"
  ) +
  theme_classic() +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.8),  # keep only x & y axis
    
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    
    legend.position = "none",
    
    panel.grid = element_blank(),
    panel.border = element_blank(),
    
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(here::here("out_figures", "PFAS_percent_above_MDL.pdf"), 
       width = 10, height = 6, dpi = 300)


# PFAS that have ≤25% of values below MDL
eligible_pfas_1m <- percent_below_mdl_by_tp %>%
  filter(timepoint == "1 Month", percent_below_MDL <= 25) %>%
  pull(PFAS_names)
eligible_pfas_6m <- percent_below_mdl_by_tp %>%
  filter(timepoint == "6 Month", percent_below_MDL <= 25) %>%
  pull(PFAS_names)

# Select PFAS that have >25 and <80% of values below MDL
eligible_pfas_1mDetect <- percent_below_mdl_by_tp %>%
  filter(timepoint == "1 Month",
         percent_below_MDL > 25,
         percent_below_MDL < 80) %>%
  pull(PFAS_names)

eligible_pfas_6mDetect <- percent_below_mdl_by_tp %>%
  filter(timepoint == "6 Month",
         percent_below_MDL > 25,
         percent_below_MDL < 80) %>%
  pull(PFAS_names)


# Keep only PFAS that has ≤25% of values below MDL
PFAS_1m <- PFAS %>%
  filter(timepoint == 1) %>%
  select(dyad_id, timepoint, all_of(as.character(eligible_pfas_1m)))

PFAS_6m <- PFAS %>%
  filter(timepoint == 6) %>%
  select(dyad_id, timepoint, all_of(as.character(eligible_pfas_6m)))

# Now also keep PFAS that has 25-80% of values below MDL
PFAS_1mDetect <- PFAS %>%
  filter(timepoint == 1) %>%
  select(dyad_id, timepoint, all_of(as.character(eligible_pfas_1mDetect)))

PFAS_6mDetect <- PFAS %>%
  filter(timepoint == 6) %>%
  select(dyad_id, timepoint, all_of(as.character(eligible_pfas_6mDetect)))

# Make sure all participants w/ PFAS data are found in meta data
unique(PFAS_1m$dyad_id) %in% unique(meta_trim$dyad_id) # All true
unique(PFAS_6m$dyad_id) %in% unique(meta_trim$dyad_id)
unique(PFAS_1mDetect$dyad_id) %in% unique(meta_trim$dyad_id)
unique(PFAS_6mDetect$dyad_id) %in% unique(meta_trim$dyad_id)

# Merge PFAS and meta_trim 
PFAS_1m <- PFAS_1m %>%
  inner_join(meta_trim, by = c("dyad_id", "timepoint"))
PFAS_6m <- PFAS_6m %>%
  inner_join(meta_trim, by = c("dyad_id", "timepoint"))
PFAS_1mDetect <- PFAS_1mDetect %>%
  inner_join(meta_trim, by = c("dyad_id", "timepoint"))
PFAS_6mDetect <- PFAS_6mDetect %>%
  inner_join(meta_trim, by = c("dyad_id", "timepoint"))


# Create detect/non-detect for select PFAS files
# Function to create detect/non-detect columns
create_detect_columns <- function(df, mdl_vals) {
  
  # Identify PFAS columns (those ending in _pgmL)
  pfas_cols <- grep("_pgmL$", colnames(df), value = TRUE)
  
  for (pfas in pfas_cols) {
    
    # Get MDL value for this PFAS
    mdl <- mdl_vals[pfas]
    
    # Create new column name: remove _pgmL and add _Detect
    # e.g., "PFHpA_pgmL" -> "PFHpA_Detect"
    new_colname <- gsub("_pgmL$", "_Detect", pfas)
    
    # Create detect/non-detect column
    df[[new_colname]] <- ifelse(df[[pfas]] >= mdl, "detect", "non-detect")
    
    cat(pfas, "-> ", new_colname, "(MDL =", mdl, ")\n")
  }
  
  return(df)
}

# Apply to both datasets
PFAS_1mDetect <- create_detect_columns(PFAS_1mDetect, mdl_values)
PFAS_6mDetect <- create_detect_columns(PFAS_6mDetect, mdl_values)

# Verify the new columns
grep("_Detect$", colnames(PFAS_1mDetect), value = TRUE)
grep("_Detect$", colnames(PFAS_6mDetect), value = TRUE)

# Create n_detect (count of detected PFAS)
create_n_detect <- function(df) {
  
  detect_cols <- grep("_Detect$", colnames(df), value = TRUE)
  
  # Convert detect/non-detect to numeric (1 = detect, 0 = non-detect)
  detect_numeric <- df[, detect_cols]
  detect_numeric[] <- lapply(detect_numeric, function(x) as.numeric(x == "detect"))
  
  # Number of PFAS detected per participant
  df$n_detect <- rowSums(detect_numeric, na.rm = TRUE)
  
  # Print summary
  cat("n_detect range:", min(df$n_detect, na.rm = TRUE), "-", max(df$n_detect, na.rm = TRUE), "\n")
  
  return(df)
}

PFAS_1mDetect <- create_n_detect(PFAS_1mDetect)
PFAS_6mDetect <- create_n_detect(PFAS_6mDetect)

# Prepare data for cross-sectional analysis-------------------------------------
# Bring rowname on species file to column
species <- species %>%
  rownames_to_column(var = "merge_id_dyad")

# Select species of 1 month and 6 months
species_1m <- species %>% filter(substr(merge_id_dyad, nchar(merge_id_dyad)-1, nchar(merge_id_dyad)) == "01")
species_6m <- species %>% filter(substr(merge_id_dyad, nchar(merge_id_dyad)-1, nchar(merge_id_dyad)) == "06")

# Now merge 1 month PFAS with 1 month microbiome
PFAS1m_micro1m <- inner_join(
  PFAS_1m,
  species_1m,
  by = "merge_id_dyad"
)

PFAS1mDetect_micro1m <- inner_join(
  PFAS_1mDetect,
  species_1m,
  by = "merge_id_dyad"
)

# To merge 1m PFAS with 6m microbiome, timepoint need to be changed on merge_id_dyad
PFAS_1m_temp <- PFAS_1m %>%
  mutate(
    merge_id_dyad = sub("-01$", "-06", merge_id_dyad)
  )

PFAS_1mDetect_temp <- PFAS_1mDetect %>%
  mutate(
    merge_id_dyad = sub("-01$", "-06", merge_id_dyad)
  )

# Merge 1 month PFAS with 6 month microbiome
PFAS1m_micro6m <- inner_join(
  PFAS_1m_temp,  # timepoint changed to -06
  species_6m,
  by = "merge_id_dyad"
)

PFAS1mDetect_micro6m <- inner_join(
  PFAS_1mDetect_temp,  # timepoint changed to -06
  species_6m,
  by = "merge_id_dyad"
)

# Save PFAS_meta file here to create 1m and 6m PFAS boxplot, correlation plot, and get table for median and % contribution
write.csv(PFAS1m_micro1m, here::here("out_files", "PFAS_1m_meta.csv"))
write.csv(PFAS1m_micro6m, here::here("out_files", "PFAS_6m_meta.csv"))
write.csv(PFAS1mDetect_micro1m, here::here("out_files", "PFAS_1mDetect_meta.csv"))
write.csv(PFAS1mDetect_micro6m, here::here("out_files", "PFAS_6mDetect_meta.csv"))


# Replace 0s with pseudo-zero in detect files for actual PFAS data (for sensitivity analysis- suggested by Lida/Tanya)
detect_pfas_cols_1m <- grep("_pgmL$", colnames(PFAS1mDetect_micro1m), value = TRUE)
detect_pfas_cols_6m <- grep("_pgmL$", colnames(PFAS1mDetect_micro6m), value = TRUE)

PFAS1mDetect_micro1m[detect_pfas_cols_1m] <- lapply(
  PFAS1mDetect_micro1m[detect_pfas_cols_1m],
  function(x) ifelse(x == 0, 0.001, x)
)
PFAS1mDetect_micro6m[detect_pfas_cols_6m] <- lapply(
  PFAS1mDetect_micro6m[detect_pfas_cols_6m],
  function(x) ifelse(x == 0, 0.001, x)
)

# Check how many NAs in covariates of interest so that we can impute them
vars_to_check <- c(
  "breastmilk_per_day",
  "mode_of_delivery_cat",
  "baby_birthweight_kg",
  "gestational_age_cat",
  "SES_index_final",
  "age_of_solid_foods"
)

pfas_list <- list(
  PFAS1m_micro1m = PFAS1m_micro1m,
  PFAS1m_micro6m = PFAS1m_micro6m,
  PFAS1mDetect_micro1m = PFAS1mDetect_micro1m,
  PFAS1mDetect_micro6m = PFAS1mDetect_micro6m
)

na_counts <- lapply(pfas_list, function(df) {
  sapply(df[, vars_to_check], function(x) sum(is.na(x)))
})

na_counts


# Impute missing values on covariates of interest --------------------
# Mode function
get_mode <- function(x) {
  ux <- na.omit(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Function to impute missing values in a dataframe
impute_covariates <- function(df) {
  
  # Impute categorical variables with mode
  if (sum(is.na(df$gestational_age_cat)) > 0) {
    df$gestational_age_cat[is.na(df$gestational_age_cat)] <- 
      get_mode(df$gestational_age_cat)
  }
  
  # Impute continuous variables with median
  if (sum(is.na(df$SES_index_final)) > 0) {
    df$SES_index_final[is.na(df$SES_index_final)] <- 
      median(df$SES_index_final, na.rm = TRUE)
  }
  
  if (sum(is.na(df$age_of_solid_foods)) > 0) {
    df$age_of_solid_foods[is.na(df$age_of_solid_foods)] <- 
      median(df$age_of_solid_foods, na.rm = TRUE)
  }
  
  return(df)
}

# Apply imputation to each dataset
PFAS1m_micro1m <- impute_covariates(PFAS1m_micro1m)
PFAS1m_micro6m <- impute_covariates(PFAS1m_micro6m)
PFAS1mDetect_micro1m <- impute_covariates(PFAS1mDetect_micro1m)
PFAS1mDetect_micro6m <- impute_covariates(PFAS1mDetect_micro6m)

# Verify imputation worked
vars_to_check <- c(
  "breastmilk_per_day",
  "mode_of_delivery_cat",
  "baby_birthweight_kg",
  "gestational_age_cat",
  "SES_index_final",
  "age_of_solid_foods"
)

pfas_list <- list(
  PFAS1m_micro1m = PFAS1m_micro1m,
  PFAS1m_micro6m = PFAS1m_micro6m,
  PFAS1mDetect_micro1m = PFAS1mDetect_micro1m,
  PFAS1mDetect_micro6m = PFAS1mDetect_micro6m
)

na_counts_after <- lapply(pfas_list, function(df) {
  sapply(df[, vars_to_check], function(x) sum(is.na(x)))
})

na_counts_after


# Impute values below MDL with MDL/sqrt(2) for 1 month pfas-1m micro data
for(pfas in eligible_pfas_1m) {
  mdl <- mdl_values[pfas]
  PFAS1m_micro1m[[pfas]] <- ifelse(
    PFAS1m_micro1m[[pfas]] < mdl,
    mdl / sqrt(2),
    PFAS1m_micro1m[[pfas]]
  )
}

# Impute values below MDL with MDL/sqrt(2) for 1 month pfas-6m micro data
for(pfas in eligible_pfas_1m) {
  mdl <- mdl_values[pfas]
  PFAS1m_micro6m[[pfas]] <- ifelse(
    PFAS1m_micro6m[[pfas]] < mdl,
    mdl / sqrt(2),
    PFAS1m_micro6m[[pfas]]
  )
}

# Drop participants with missing PFAS values (1-2 per compound, <2% missingness)
# Complete case analysis is appropriate at this level of missingness
PFAS1m_micro1m <- PFAS1m_micro1m %>%
  filter(if_all(all_of(eligible_pfas_1m), ~ !is.na(.)))

PFAS1m_micro6m <- PFAS1m_micro6m %>%
  filter(if_all(all_of(eligible_pfas_1m), ~ !is.na(.)))

PFAS1mDetect_micro1m <- PFAS1mDetect_micro1m %>%
  filter(if_all(all_of(detect_pfas_cols_1m), ~ !is.na(.)))

PFAS1mDetect_micro6m <- PFAS1mDetect_micro6m %>%
  filter(if_all(all_of(detect_pfas_cols_6m), ~ !is.na(.)))

# Verify
cat("After complete case exclusion:\n")
cat("PFAS1m_micro1m:", nrow(PFAS1m_micro1m), "\n")
cat("PFAS1m_micro6m:", nrow(PFAS1m_micro6m), "\n")
cat("PFAS1mDetect_micro1m:", nrow(PFAS1mDetect_micro1m), "\n")
cat("PFAS1mDetect_micro6m:", nrow(PFAS1mDetect_micro6m), "\n")

# Define PFAS variables
pfas_vars1 <- colnames(PFAS1m_micro1m)[3:7]
pfas_vars2 <- colnames(PFAS1m_micro6m)[3:7]

# Define meta_vars as all other columns except those in pfas_vars
meta_vars1 <- setdiff(colnames(PFAS1m_micro1m), pfas_vars1)
meta_vars2 <- setdiff(colnames(PFAS1m_micro6m), pfas_vars2)

# Log2 transform PFAS data
PFAS1m_micro1m[pfas_vars1] <- log2(PFAS1m_micro1m[pfas_vars1])
PFAS1m_micro6m[pfas_vars2] <- log2(PFAS1m_micro6m[pfas_vars2])
# For pfas in binary file too (it is a part of sensitivity analysis)
PFAS1mDetect_micro1m[detect_pfas_cols_1m] <- log2(PFAS1mDetect_micro1m[detect_pfas_cols_1m])
PFAS1mDetect_micro6m[detect_pfas_cols_6m] <- log2(PFAS1mDetect_micro6m[detect_pfas_cols_6m])

# Select only needed variables
microbeStart <- 22
microbeEnd <- 161

PFAS1m_micro1m <- PFAS1m_micro1m %>%
  select(
    merge_id_dyad,
    mode_of_delivery_cat,
    baby_birthweight_kg,
    gestational_age_cat,
    breastmilk_per_day,
    SES_index_final,
    age_of_solid_foods,
    contains("pgmL"),
    microbeStart:microbeEnd) #Warning here is okay

PFAS1m_micro6m <- PFAS1m_micro6m %>%
  select(merge_id_dyad,
         mode_of_delivery_cat,
         baby_birthweight_kg,
         gestational_age_cat,
         breastmilk_per_day,
         SES_index_final,
         age_of_solid_foods,
         contains("pgmL"),
         microbeStart:microbeEnd)

# Similar for detects/non-detects
microbeStart <- 28
microbeEnd <- 167

PFAS1mDetect_micro1m <- PFAS1mDetect_micro1m %>%
  select(
    merge_id_dyad,
    mode_of_delivery_cat,
    baby_birthweight_kg,
    gestational_age_cat,
    breastmilk_per_day,
    SES_index_final,
    age_of_solid_foods,
    contains("_pgmL"),
    contains("detect"),
    microbeStart:microbeEnd)

PFAS1mDetect_micro6m <- PFAS1mDetect_micro6m %>%
  select(merge_id_dyad,
         mode_of_delivery_cat,
         baby_birthweight_kg,
         gestational_age_cat,
         breastmilk_per_day,
         SES_index_final,
         age_of_solid_foods,
         contains("_pgmL"),
         contains("detect"),
         microbeStart:microbeEnd)

# Check if any NAs
sum(is.na(PFAS1m_micro1m))
sum(is.na(PFAS1m_micro6m))
sum(is.na(PFAS1mDetect_micro1m))
sum(is.na(PFAS1mDetect_micro6m))

# # Remove rows if any NAs
# PFAS1m_micro1m <- PFAS1m_micro1m %>%
#   na.omit()
# PFAS1m_micro6m <- PFAS1m_micro6m %>%
#   na.omit()
# PFAS1mDetect_micro1m <- PFAS1mDetect_micro1m %>%
#   na.omit()
# PFAS1mDetect_micro6m <- PFAS1mDetect_micro6m %>%
#   na.omit()

# Save files for downstream analysis
write.csv(PFAS1m_micro1m, here::here("out_files", "PFAS1m_micro1m_species.csv"))
write.csv(PFAS1m_micro6m, here::here("out_files", "PFAS1m_micro6m_species.csv"))
write.csv(PFAS1mDetect_micro1m, here::here("out_files", "PFAS1mDetect_micro1m_species.csv"))
write.csv(PFAS1mDetect_micro6m, here::here("out_files", "PFAS1mDetect_micro6m_species.csv"))

# Prepare files for other taxonomic levels (genus through phylum) --------------
# Extract the already-processed PFAS + covariate columns from species files
# (log2 transformed PFAS, imputed covariates) — just swap out the microbe and count columns

# Extract PFAS + covariate columns only (no taxa, no counts column)
pfas_meta_1m        <- PFAS1m_micro1m       %>% dplyr::select(-matches("^X\\d+"))
pfas_meta_6m        <- PFAS1m_micro6m       %>% dplyr::select(-matches("^X\\d+"))
pfas_meta_1mDetect  <- PFAS1mDetect_micro1m %>% dplyr::select(-matches("^X\\d+"))
pfas_meta_6mDetect  <- PFAS1mDetect_micro6m %>% dplyr::select(-matches("^X\\d+"))

other_levels <- list(
  genus  = genus,
  family = family,
  order  = order,
  class  = class,
  phylum = phylum
)

for (lvl_name in names(other_levels)) {

  taxa_df <- other_levels[[lvl_name]] %>%
    rownames_to_column(var = "merge_id_dyad")

  taxa_1m <- taxa_df %>% filter(grepl("-01$", merge_id_dyad))
  taxa_6m <- taxa_df %>% filter(grepl("-06$", merge_id_dyad))

  # Adjust merge_id_dyad on 6m PFAS side to match 6m microbiome rows
  pfas_meta_6m_temp       <- pfas_meta_6m       %>% mutate(merge_id_dyad = sub("-01$", "-06", merge_id_dyad))
  pfas_meta_6mDetect_temp <- pfas_meta_6mDetect %>% mutate(merge_id_dyad = sub("-01$", "-06", merge_id_dyad))

  out_1m       <- inner_join(pfas_meta_1m,        taxa_1m, by = "merge_id_dyad")
  out_6m       <- inner_join(pfas_meta_6m_temp,   taxa_6m, by = "merge_id_dyad")
  out_1mDetect <- inner_join(pfas_meta_1mDetect,  taxa_1m, by = "merge_id_dyad")
  out_6mDetect <- inner_join(pfas_meta_6mDetect_temp, taxa_6m, by = "merge_id_dyad")

  assign(paste0("PFAS1m_micro1m_",       lvl_name), out_1m)
  assign(paste0("PFAS1m_micro6m_",       lvl_name), out_6m)
  assign(paste0("PFAS1mDetect_micro1m_", lvl_name), out_1mDetect)
  assign(paste0("PFAS1mDetect_micro6m_", lvl_name), out_6mDetect)

  write.csv(out_1m,       here::here("out_files", paste0("PFAS1m_micro1m_",       lvl_name, ".csv")), row.names = FALSE)
  write.csv(out_6m,       here::here("out_files", paste0("PFAS1m_micro6m_",       lvl_name, ".csv")), row.names = FALSE)
  write.csv(out_1mDetect, here::here("out_files", paste0("PFAS1mDetect_micro1m_", lvl_name, ".csv")), row.names = FALSE)
  write.csv(out_6mDetect, here::here("out_files", paste0("PFAS1mDetect_micro6m_", lvl_name, ".csv")), row.names = FALSE)
}

# END