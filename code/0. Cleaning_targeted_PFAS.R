# TITLE:   1. Data Cleaning.R
#
# PURPOSE: Read and clean the meta data and microbiome data
#
# DATE:    June 2025
# CODE REVIEW:

# set up -----------------------------------------------------------------------

#clear the workspace
rm(list = ls())

#load libraries
library(readxl)
library(tidyverse)
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)

# Read targeted PFAS with no headers
file_path <- here::here("input", "CLU0031_PFAS_Final_Report.xlsx")
PFAS <- read_excel(file_path, sheet = "Sheet1", col_names = FALSE)

# Extract metadata colnames from Row 4 (Excel row 5)
meta_names <- as.character(PFAS[4, 1:5])

# Extract compound names from Row 2 (Excel row 3), cols 6+
compound_names <- as.character(PFAS[2, 6:ncol(PFAS)])
compound_names <- paste0(compound_names, "_pgmL")

# Final column names
final_colnames <- c(meta_names, compound_names)

# Drop top 4 rows and set colnames
PFAS <- PFAS[-c(1:4), ]
colnames(PFAS) <- final_colnames
rownames(PFAS) <- NULL

# Convert compound measurement columns to numeric
compound_cols <- (length(meta_names) + 1):ncol(PFAS)
PFAS[compound_cols] <- lapply(PFAS[compound_cols], as.numeric)

# Add a flag column to record omitted compounds
PFAS$Omitted_Compounds <- ""

# Define outlier internal standards---------------------------------------------
# These are coming from Quantification report from Doug's Lab. These sample IDs PFAS has their internal standard peak outside 3SD
compounds_to_omit <- list(
  `92_6M`     = c("PFBS", "PFPeAS"),
  `MA-072_BL` = c("PFBA"),
  `MA-036_BL` = c("PFUnA", "PFTeDA"),
  `205_BL`    = c("PFTeDA"),
  `14_BL`     = c("PFBA"),
  `166_BL`    = c("PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA", "PFUnA", "PFDoA", "PFTeDA", "PFBA", "PFBS", "PFPeAS"), #last two were missed. so I added
  `212_6M`    = c("PFBA"),
  `106_BL`    = c("PFPeA", "PFOA", "PFHxA", "PFHpA", "PFBA"),
  `79_BL`     = c("PFBA"),
  `MA-108_BL` = c("PFBA"),
  `MA-051_BL` = c("PFBA"),
  `24_BL`     = c("PFBA"),
  `MA-081_BL` = c("PFOS", "PFNS", "PFTrDA", "PFBS", "PFHxS", "PFHps", "PFPeAS", "PFDoS")
)

# Set NAs and update flags according to Doug's documentation--------------------
for (sample_id in names(compounds_to_omit)) {
  compounds <- compounds_to_omit[[sample_id]]
  row_index <- which(PFAS$Sample.ID == sample_id)  # corrected from df to PFAS
  
  flagged_compounds <- c()  # track what was set to NA
  
  for (compound in compounds) {
    col_name <- paste0(compound, "_pgmL")
    if (col_name %in% colnames(PFAS)) {
      PFAS[row_index, col_name] <- NA
      flagged_compounds <- c(flagged_compounds, compound)
    } else {
      warning(sprintf("Compound column '%s' not found for sample '%s'", col_name, sample_id))
    }
  }
  
  # Update the flag column with a comma-separated list of omitted compounds
  if (length(flagged_compounds) > 0) {
    PFAS$Omitted_Compounds[row_index] <- paste(flagged_compounds, collapse = ", ")
  }
}

# Add a column identifying the study source
PFAS$Study <- ifelse(grepl("^MA", PFAS$Sample.ID), "MAMITAS", "Mothers Milk")

# Convert timepoint labels to numeric values
PFAS$Timepoint <- sub(".*_", "", PFAS$Sample.ID)  # extract
PFAS$Timepoint <- recode(PFAS$Timepoint,
                       "BL" = 1,
                       "6M" = 6,
                       "12M" = 12,
                       "24M" = 24)

# View updated dataframe
head(PFAS)

# Create dyad_id to merge later
PFAS$dyad_id <- sub("_.*", "", PFAS$Sample.ID)

# # Save MM+MA (Elsendero milk pfas)
# write.csv(PFAS, here::here("output", "PFAS_1m6m_elsendero_quantified.csv"))

# Now keep only Mother's milk data
PFAS <- PFAS %>%
  dplyr::filter(Study == "Mothers Milk")

# Save the cleaned file ready to merge
write.csv(PFAS, here::here("out_files", "PFAS_1m6m_quantified.csv"))



