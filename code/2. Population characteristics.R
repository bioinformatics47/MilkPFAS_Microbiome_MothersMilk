# TITLE:   1. Population characteristics.R
#
# PURPOSE: Extract population data
#
# DATE:    January 12, 2024
# CODE REVIEW:
# Reviewed by Ellie Holzhausen (EAH) on April 23, 2026
# by Haonan Li (HL) on May 25, 2026
#
# set up -----------------------------------------------------------------------
rm(list = ls())

library(tidyverse)
library(dplyr)
library(tableone)
library(kableExtra)
library(flextable)
library(officer)
library(here)
library(knitr)
library(ggplot2)
library(scales)
library(ggpubr)
library(openxlsx)


# Get analytic sample IDs from the analysis-ready file
analytic_ids <- read.csv(here::here("out_files", "PFAS1m_micro1m_species.csv")) %>%
  pull(merge_id_dyad)

# Load full meta and restrict to analytic sample
pop <- read.csv(here::here("out_files", "PFAS_1m_meta.csv")) %>%
  select(-X) %>%
  filter(merge_id_dyad %in% analytic_ids)

tableDF <- pop %>% 
  select(-merge_id_dyad)

# Recode and factor categorical variables
tableDF <- tableDF %>%
  mutate(
    prepreg_bmi_cat = factor(prepreg_bmi_cat,
                             levels = c("Normal", "Overweight", "Obesity")),
    mode_of_delivery_cat = factor(mode_of_delivery_cat,
                                  levels = c("Vaginal", "C-Section")),
    gestational_age_cat = factor(gestational_age_cat,
                                 levels = c("Early", "Ontime", "Late"))
  )

table_vars <- c("mother_age",
                "SES_index_final",
                "mother_antibiotics",
                "prepreg_bmi_cat",
                "mode_of_delivery_cat",
                "breastmilk_per_day",
                "age_in_days",
                "baby_gender",
                "baby_birthweight_kg",
                "baby_antibiotics",
                "gestational_age_cat",
                "age_of_solid_foods")

cont_vars <- c("age_in_days", "breastmilk_per_day", "mother_age", 
               "prbaby_birthweight_kg", "SES_index_final", 
               "age_of_solid_foods")

cat_vars <- c("baby_gender", "mode_of_delivery_cat", "baby_antibiotics",
              "mother_antibiotics", "gestational_age_cat", "prepreg_bmi_cat")

mom_table <- CreateTableOne(vars = table_vars,
                            data = tableDF,         
                            factorVars = cat_vars)

mom_print <- print(mom_table,
                   quote = FALSE,
                   noSpaces = TRUE,
                   printToggle = FALSE,
                   showAllLevels = TRUE)

mom_df <- as.data.frame(mom_print)
if (ncol(mom_df) > 2 && colnames(mom_df)[ncol(mom_df)] == "") {
  mom_df <- mom_df[, -ncol(mom_df)]
}

colnames(mom_df) <- c("Level", "1-Month")

mom_df$Characteristics <- c(
  "N",
  "Maternal Age (Years) (SD)",
  "Socio-economic status (SES) Index (SD)",
  "Maternal Antibiotics (N, %)",          # appears on first level row (No)
  "",                                      # Yes level
  "Pre-pregnancy BMI Category (N, %)",    # appears on Normal level row
  "",                                      # Overweight level
  "",                                      # Obesity level
  "Mode of Delivery (N, %)",              # appears on Vaginal level row
  "",                                      # C-Section level
  "Breastmilk per Day (SD)",
  "Infant Age (Days) (SD)",
  "Infant Sex (N, %)",                    # appears on Female level row
  "",                                      # Male level
  "Infant Birthweight (kg) (SD)",
  "Infant Antibiotics (N, %)",            # appears on No level row
  "",                                      # Yes level
  "Gestational Age (N, %)",              # appears on Early level row
  "",                                      # Ontime level
  "",                                      # Late level
  "Age of Solid Foods Started (SD)"
)

# Helper function to format min - max range
fmt_range <- function(x, digits = 2) {
  paste0(round(min(x, na.rm = TRUE), digits), " - ", round(max(x, na.rm = TRUE), digits))
}

# Build range vector: one entry per row of mom_df
# Rows for continuous variables get min - max; all others get ""
# Row order must match the Characteristics vector above
range_vec <- c(
  "",                                           # N
  fmt_range(tableDF$mother_age),                # Maternal Age
  fmt_range(tableDF$SES_index_final),           # SES Index
  "",                                           # Maternal Antibiotics No
  "",                                           # Yes
  "",                                           # BMI Normal
  "",                                           # Overweight
  "",                                           # Obesity
  "",                                           # Mode of Delivery Vaginal
  "",                                           # C-Section
  fmt_range(tableDF$breastmilk_per_day),        # Breastmilk per Day
  fmt_range(tableDF$age_in_days, digits = 0),   # Infant Age
  "",                                           # Infant Sex Female
  "",                                           # Male
  fmt_range(tableDF$baby_birthweight_kg),       # Infant Birthweight
  "",                                           # Infant Antibiotics No
  "",                                           # Yes
  "",                                           # Gestational Age Early
  "",                                           # Ontime
  "",                                           # Late
  fmt_range(tableDF$age_of_solid_foods)         # Age of Solid Foods
)

# Add range
mom_df$Range <- range_vec

# Reorder columns: Characteristics first
mom_df <- mom_df %>%
  relocate(Characteristics, .before = Level)

mom_df <- mom_df[, colSums(is.na(mom_df) | mom_df == "") != nrow(mom_df)]

# Quick kable preview
mom_df %>%
  kable("pipe", align = "c") %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE)

# Create flextable
ft <- mom_df %>%
  flextable() %>%
  theme_vanilla() %>%
  fontsize(size = 10, part = "all") %>%
  bold(part = "header") %>%
  
  align(align = "center", part = "header") %>%
  align(align = "left", j = 1:2, part = "all") %>%
  align(align = "center", j = 3:4, part = "all") %>%  # 1-Month and Range centered
  
  bg(j = 1:ncol(mom_df), 
     bg = ifelse(seq_len(nrow(mom_df)) %% 2 == 0, "#f5f5f5", "white"), 
     part = "body") %>%
  
  border_remove() %>%
  hline_top(border = fp_border(color = "black", width = 2), part = "header") %>%
  hline_bottom(border = fp_border(color = "black", width = 2), part = "header") %>%
  hline_bottom(border = fp_border(color = "black", width = 2), part = "body") %>%
  
  autofit() %>%
  width(j = 1, width = 3.5) %>%   # Characteristics
  width(j = 2, width = 1.5) %>%   # Level
  width(j = 3, width = 1.8) %>%   # 1-Month
  width(j = 4, width = 1.8) %>%   # Range
  fit_to_width(max_width = 7.5)

# Export to Word doc
doc <- read_docx() %>%
  body_add_par("") %>%
  body_add_flextable(value = ft) %>%
  body_add_par("")

print(doc, target = here::here("out_figures", "descriptive_table_1month_portrait.docx"))


# Calculate Median, IQR, % contribution of PFAS --------------------------------

# Analytic IDs for continuous PFAS
analytic_ids_cont <- read.csv(here::here("out_files", "PFAS1m_micro1m_species.csv")) %>%
  pull(merge_id_dyad)
analytic_dyads_cont <- unique(sub("MM-(\\d+)-\\d+", "\\1", analytic_ids_cont)) %>% as.numeric()

# Analytic IDs for binary PFAS
analytic_ids_detect <- read.csv(here::here("out_files", "PFAS1mDetect_micro1m_species.csv")) %>%
  pull(merge_id_dyad)
analytic_dyads_detect <- unique(sub("MM-(\\d+)-\\d+", "\\1", analytic_ids_detect)) %>% as.numeric()

# Load and filter each to its own analytic sample
PFAS_meta1 <- read.csv(here::here("out_files", "PFAS_1m_meta.csv")) %>% 
  dplyr::select(-X) %>%
  filter(dyad_id %in% analytic_dyads_cont) %>%
  select(dyad_id, timepoint, contains("_pgmL"))

PFAS_meta2 <- read.csv(here::here("out_files", "PFAS_1mDetect_meta.csv")) %>% 
  dplyr::select(-X) %>%
  filter(dyad_id %in% analytic_dyads_detect) %>%
  select(dyad_id, timepoint, contains("_pgmL"))

PFAS_meta <- PFAS_meta1 %>%
  full_join(PFAS_meta2, by = c("dyad_id", "timepoint"))

# DEFINE MDL VALUES
mdl_values <- c(
  "N.EtFOSAA_pgmL" = 1.12, "N.MeFOSAA_pgmL" = 0.74, "PFBA_pgmL" = 3.33,
  "PFDA_pgmL" = 1.02, "PFDoA_pgmL" = 0.74, "PFHpA_pgmL"  = 1.16,
  "PFHxA_pgmL" = 3.93, "PFNA_pgmL" = 0.66, "PFOA_pgmL"  = 2.28,
  "PFPeA_pgmL"  = 4.47, "PFTeDA_pgmL" = 0.41, "PFUnA_pgmL" = 1.16,
  "PFBS_pgmL" = 7.00, "PFDoS_pgmL" = 0.83, "PFHps_pgmL" = 1.30,
  "PFHxS_pgmL" = 0.62, "PFNS_pgmL"  = 1.44, "PFTrDA_pgmL" = 0.83,
  "PFOS_pgmL" = 1.09, "PFPeAS_pgmL" = 0.69
)

mdl_df <- tibble(PFAS = names(mdl_values), MDL = mdl_values)

# 1-MONTH PFAS SUMMARY ---------------------------------------------------------

pfas_1m_long <- PFAS_meta %>%
  dplyr::filter(timepoint == 1) %>%
  tidyr::pivot_longer(
    cols = contains("_pgmL"),
    names_to = "PFAS",
    values_to = "concentration"
  ) %>%
  dplyr::left_join(mdl_df, by = "PFAS") %>%
  dplyr::mutate(
    # Detected-only: sub-MDL → NA. This is the universe for all summary stats.
    concentration_detected = if_else(concentration > MDL, concentration, NA_real_)
  )

# Summary stats: detected values only, flag low-detection compounds
pfas_1m_summary <- pfas_1m_long %>%
  dplyr::group_by(PFAS, MDL) %>%
  dplyr::summarise(
    n_total      = sum(!is.na(concentration)),
    n_detected   = sum(!is.na(concentration_detected)),
    pct_detected = round(n_detected / n_total * 100, 1),
    # Stats computed on detected values only
    Median = median(concentration_detected, na.rm = TRUE),
    Q1     = quantile(concentration_detected, 0.25, na.rm = TRUE),
    Q3     = quantile(concentration_detected, 0.75, na.rm = TRUE),
    Max    = max(concentration_detected, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    # Suppress median/IQR when >50% of values are below MDL because
    # the sample median itself would fall below MDL and is not estimable
    # from detected values alone
    Median = if_else(pct_detected < 50, NA_real_, Median),
    Q1     = if_else(pct_detected < 50, NA_real_, Q1),
    Q3     = if_else(pct_detected < 50, NA_real_, Q3),
    # Format IQR as a single string for the table
    IQR_range = if_else(
      is.na(Q1), "< MDL",
      paste0(round(Q1, 2), " - ", round(Q3, 2))
    ),
    Median_display = if_else(is.na(Median), "< MDL", as.character(round(Median, 2)))
  )

# % Contribution uses MDL/sqrt(2) for non-detects to enable rowSums.
# For compounds with high detection (>75%), median contribution reflects
# observed data. For low-detection compounds, contributions are dominated
# by the MDL/sqrt(2) constant and should be interpreted cautiously.
pfas_1m_wide_imputed <- PFAS_meta %>%
  dplyr::filter(timepoint == 1) %>%
  dplyr::select(dyad_id, timepoint, all_of(names(mdl_values)[names(mdl_values) %in% names(.)]))

# Apply MDL/sqrt(2) imputation per column
for (col in intersect(names(mdl_values), names(pfas_1m_wide_imputed))) {
  mdl_val <- mdl_values[[col]]
  pfas_1m_wide_imputed[[col]] <- if_else(
    pfas_1m_wide_imputed[[col]] <= mdl_val,
    mdl_val / sqrt(2),          # impute sub-MDL with MDL/√2
    pfas_1m_wide_imputed[[col]] # keep detected values as-is
  )
}

pfas_1m_contribution <- pfas_1m_wide_imputed %>%
  dplyr::mutate(
    total_PFAS = rowSums(dplyr::select(., contains("_pgmL")), na.rm = TRUE)
  ) %>%
  tidyr::pivot_longer(
    cols = contains("_pgmL"),
    names_to = "PFAS",
    values_to = "concentration_imputed"
  ) %>%
  dplyr::mutate(
    pct_contribution = (concentration_imputed / total_PFAS) * 100
  ) %>%
  dplyr::group_by(PFAS) %>%
  dplyr::summarise(
    percent_contribution = round(median(pct_contribution, na.rm = TRUE), 2),
    .groups = "drop"
  )

# Combine all summaries
pfas_1m_final <- pfas_1m_summary %>%
  dplyr::left_join(pfas_1m_contribution, by = "PFAS") %>%
  dplyr::select(
    PFAS,
    MDL,
    Median_display,
    IQR_range,
    Max,
    percent_contribution,
    n_detected,
    n_total,
    pct_detected
  ) %>%
  dplyr::rename(
    `Median (pg/mL)`       = Median_display,
    `IQR (pg/mL)`          = IQR_range,
    `Max (pg/mL)`          = Max,
    `% Contribution`       = percent_contribution,
    `N above MDL`          = n_detected,
    `N`                    = n_total,
    `% Above MDL`          = pct_detected
  )

write.xlsx(pfas_1m_final,
           here::here("out_files", "pfas_1m_summary.xlsx"), rowNames = FALSE)


# 1m PFAS Boxplots--------------------------------------------------------------

# Prepare data — keep original column names for MDL matching
pfas_cols_in_data <- PFAS_meta %>%
  select(contains("_pgmL")) %>%
  names()

# Now keep only the MDL entries that match real columns
existing_mdl <- mdl_values[names(mdl_values) %in% pfas_cols_in_data]

# Start from the wide PFAS data at 1 month only
PFAS_1m_wide <- PFAS_meta %>%
  filter(timepoint == 1) %>%
  select(dyad_id, timepoint, all_of(names(existing_mdl))) %>%
  mutate(across(everything(), as.numeric))

# Apply MDL filter: ≤ MDL → NA  (only keep values truly above MDL)
for (col in names(existing_mdl)) {
  PFAS_1m_wide[[col]] <- ifelse(
    PFAS_1m_wide[[col]] > existing_mdl[[col]],
    PFAS_1m_wide[[col]],
    NA_real_
  )
}

pfas_display_names <- c(
  "PFBS_pgmL"      = "PFBS",
  "PFHxS_pgmL"     = "PFHxS",
  "PFNA_pgmL"      = "PFNA",
  "PFOA_pgmL"      = "PFOA",
  "PFOS_pgmL"      = "PFOS",
  "N.MeFOSAA_pgmL" = "N-MeFOSAA",
  "PFBA_pgmL"      = "PFBA",
  "PFDA_pgmL"      = "PFDA",
  "PFDoA_pgmL"     = "PFDoA",
  "PFHpA_pgmL"     = "PFHpA"
  # add others if you later include more PFAS with enough detections
)

continuously_measured <- c("PFBS", "PFHxS", "PFNA", "PFOA", "PFOS")
binary_classified    <- c("N-MeFOSAA", "PFBA", "PFDA", "PFDoA", "PFHpA")

desired_order <- c(continuously_measured, binary_classified)


# Pivot → filter above MDL → log2 → factor levels
df_long_1m <- PFAS_1m_wide %>%
  # pivot all PFAS columns that exist
  pivot_longer(
    cols = any_of(names(mdl_values)),
    names_to = "PFAS_col",
    values_to = "concentration_pg_mL",
    values_drop_na = FALSE
  ) %>%
  # only keep detected values
  filter(!is.na(concentration_pg_mL)) %>%
  # create nice display name
  mutate(
    PFAS = recode(PFAS_col, !!!pfas_display_names),
    log2_conc = log2(concentration_pg_mL),
    PFAS_type = if_else(
      PFAS %in% continuously_measured,
      "Continuously measured PFAS",
      "Binary classified PFAS"
    ),
    PFAS_type = factor(
      PFAS_type,
      levels = c("Continuously measured PFAS", "Binary classified PFAS")
    ),
    PFAS = factor(PFAS, levels = desired_order)
  )

# Plot
pfas_boxplot_1m <- ggplot(df_long_1m, aes(x = PFAS, y = log2_conc, fill = PFAS_type)) +
  geom_boxplot(
    outlier.size = 1.8,
    outlier.shape = 21,
    outlier.fill = "black",
    outlier.colour = "black"
  ) +
  scale_fill_manual(
    values = c("Binary classified PFAS"  = "#d7191c",
               "Continuously measured PFAS" = "#2c7bb6"
               )
  ) +
  scale_y_continuous(breaks = seq(-5, 10, by = 2)) +
  labs(
    x    = NULL,
    y    = expression("log"[2]*" concentration (pg/mL)"),
    fill = "PFAS group",
    title = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.y       = element_text(size = 20, face = "bold", color = "black"),
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
    axis.text.y        = element_text(size = 14),
    legend.position    = "top",
    legend.box.spacing = unit(-8, "pt"),
    legend.title       = element_blank(),
    legend.text        = element_text(size = 14),
    panel.grid.major.y = element_blank()
  )
print(pfas_boxplot_1m)

# Save
ggsave(
  here::here("out_figures", "PFAS_1m_boxplot_above_MDL_only.pdf"),
  pfas_boxplot_1m,
  width  = 11,
  height = 7,
  dpi    = 400,
  bg     = "white"
)

#END

