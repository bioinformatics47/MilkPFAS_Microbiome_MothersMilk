#-------------------------------------------------------------------------------
# Script 14: Mechanistic Integration — PFAS → Taxa → Gene Pathways → Metabolites
#
# PURPOSE:
#   Establish mechanistic links across three biological layers:
#
#   AXIS 1 (HMO disruption — the core story):
#     PFAS → ↓ Bifidobacterium → ↓ HMO-degradation pathways (sucrose deg., Bifido shunt,
#            UDP-GlcNAc biosyn.) → ↓ fecal HMO metabolites (fucose, rhamnose, GlcNAc)
#
#   AXIS 2 (Lachnospiraceae overgrowth):
#     PFAS → ↑ Lachnospiraceae (Dorea, Blautia, Roseburia, Coprococcus) →
#            does not restore HMO metabolism → negative correlation with HMO metabolites
#
#   AXIS 3 (Enterobacteriaceae expansion):
#     PFAS → ↑ Enterobacter/Klebsiella → ↑ L-rhamnose degradation I,
#            colanic acid biosynthesis, LPS → no benefit to HMO metabolite profile
#
#   KEY NOTE (Script 8 Figure 4 arc ring):
#     Escherichia appears as dominant contributor to POSITIVELY-associated pathways
#     because E. coli is the MetaCyc primary reference organism (richest gene catalog).
#     However, Escherichia is NEGATIVELY associated with PFAS in all taxa regression
#     models. It is NOT included as a positive-PFAS taxon in this script.
#     The true PFAS-positive community is Lachnospiraceae + Enterobacter/Klebsiella
#     + Bacteroides/Parabacteroides + Faecalibacterium prausnitzii.
#
# OUTPUTS:
#   out_figures/14_combined_fig15_fig16_v9.pdf/png  — Combined taxa-metabolite +
#                                                      pathway-metabolite bubble plots
#   out_figures/14_network_flow_unified.pdf/png      — Unified network flow diagram
#
# DATA REQUIREMENTS:
#   out_files/COMBINED_CLR_continuous_1m_6m.csv
#   out_files/PATHWAY_significant_with_species.csv
#   out_files/spearman_taxa_metab_6m.csv
#   out_files/spearman_taxa_metab_cont_1m_1m.csv
#   out_files/spearman_taxa_metab_1m.csv
#   out_files/pathway_metabolite_spearman_6m.csv
#   out_files/pathway_metabolite_spearman_1m.csv
#
#   Reviewed by Ellie Holzhausen (EAH) on April 23, 2026
#   Haonan Li (HL) on May 25, 2026
#
#-------------------------------------------------------------------------------

rm(list = ls())
# setwd() removed for portability: run this script with the working
# directory set to the project root (the folder containing this "code" dir)

library(here)
library(tidyverse)
library(ggplot2)
library(ggtext)
library(cowplot)
library(egg)
library(grid)

# ── Helper: save figures ───────────────────────────────────────────────────────
save_fig <- function(fname, w, h, expr) {
  pdf(here("out_figures", paste0(fname, ".pdf")), width = w, height = h)
  expr
  dev.off()
  png(here("out_figures", paste0(fname, ".png")), width = w * 150,
      height = h * 150, res = 150)
  expr
  dev.off()
  message("Saved: ", fname)
}

# ── Colour palette ─────────────────────────────────────────────────────────────
LABEL_NEG  <- "maroon4"
LABEL_POS  <- "#3A7EC2"
BUBBLE_NEG <- "orange2"
BUBBLE_POS <- "blue"
TAXA_NEG_COL <- "orange2"
TAXA_POS_COL <- "blue"
COL_BIF    <- "#D9622B"
COL_LAC    <- "#3A7EC2"

# ── Load data ──────────────────────────────────────────────────────────────────
pwy_species <- read.csv(here("out_files", "PATHWAY_significant_with_species.csv"))

pwy_met_1m_path <- here("out_files", "pathway_metabolite_spearman_1m.csv")
pwy_met_1m <- if (file.exists(pwy_met_1m_path)) {
  read.csv(pwy_met_1m_path)
} else {
  message("NOTE: pathway_metabolite_spearman_1m.csv not found — run Script 12 first.")
  NULL
}

# ── Define key taxa ────────────────────────────────────────────────────────────
bifido_sp <- c("Bifidobacterium bifidum", "Bifidobacterium dentium",
               "Bifidobacterium breve", "Bifidobacterium animalis",
               "Bifidobacterium pseudocatenulatum", "Bifidobacterium scardovii",
               "Bifidobacterium longum", "Bifidobacterium pseudolongum",
               "Bifidobacterium catenulatum", "Bifidobacterium actinocoloniiforme",
               "Bifidobacterium saguini", "Bifidobacterium thermophilum")

lachno_sp <- c("Dorea longicatena", "Dorea formicigenerans",
               "Blautia obeum", "Blautia wexlerae", "Blautia sp. SC05B48",
               "Coprococcus catus", "Coprococcus comes",
               "Roseburia hominis", "Roseburia intestinalis",
               "Ruminococcus torques", "Anaerobutyricum hallii",
               "Enterocloster bolteae")

entero_sp <- c("Enterobacter asburiae", "Enterobacter roggenkampii",
               "Enterobacter kobei", "Enterobacter ludwigii",
               "Enterobacter hormaechei", "Enterobacter cloacae")

# ==============================================================================
# 15.  PANEL A: TAXA GROUP × METABOLITE BUBBLE  (1m vs 6m)
#      Fill = mean Spearman rho (within taxa group); size = n significant species
# ==============================================================================

assign_grp15 <- function(nm) {
  g <- sub(" .*", "", nm)
  case_when(
    g == "Bifidobacterium"                       ~ "Bifidobacterium",
    g %in% c("Blautia","Dorea","Roseburia",
             "Coprococcus","Anaerobutyricum",
             "Enterocloster","Ruminococcus",
             "Eubacterium")                      ~ "Lachnospiraceae",
    g == "Enterobacter"                          ~ "Enterobacter",
    TRUE                                         ~ NA_character_
  )
}

story_mets15 <- tribble(
  ~met_id,                                    ~display,              ~pfas_dir,
  "6-DEOXY-L-GALACTOSE (FUCOSE)",             "Fucose/L-Rhamnose",   "DOWN",
  "N-ACETYL-D-GLUCOSAMINE",                   "GlcNAc",              "DOWN",
  "N-ACETYLNEURAMINATE",                      "Neu5Ac",              "DOWN",
  "CHOLINE",                                  "Choline",             "DOWN",
  "CHOLESTEROL",                              "Cholesterol",         "DOWN",
  "N-ACETYL-L-ASPARTIC ACID",                 "N-Acetylaspartate",   "DOWN",
  "4-PYRIDOXATE",                             "4-Pyridoxate",        "UP",
  "TRYPTAMINE",                               "Tryptamine",          "UP",
  "CADAVERINE",                               "Cadaverine",          "UP",
  "FA5:0(VALERATE, ISOVALERATE, OTHERS)",     "Valerate",            "UP",
)

met_order15 <- c(
  "Cholesterol", "Choline", "N-Acetylaspartate",
  "Fucose/L-Rhamnose", "GlcNAc", "Neu5Ac",
  "4-Pyridoxate", "Tryptamine", "Cadaverine", "Valerate"
)

grp_pal15 <- c("Bifidobacterium" = LABEL_NEG,
               "Lachnospiraceae" = LABEL_POS,
               "Enterobacter"    = LABEL_POS)

met_col_vec <- ifelse(
  story_mets15$pfas_dir[match(met_order15, story_mets15$display)] == "DOWN",
  LABEL_NEG, LABEL_POS
)

make_bubble15 <- function(spearman_df, tp) {
  spearman_df %>%
    filter(met_label %in% story_mets15$met_id) %>%
    mutate(group = assign_grp15(name)) %>%
    filter(!is.na(group)) %>%
    left_join(story_mets15 %>% select(met_id, display, pfas_dir),
              by = c("met_label" = "met_id")) %>%
    mutate(sig = p_val < 0.05) %>%
    group_by(group, display, pfas_dir) %>%
    summarise(mean_rho = mean(rho, na.rm = TRUE),
              n_sig    = sum(sig,  na.rm = TRUE),
              .groups  = "drop") %>%
    mutate(
      timepoint    = tp,
      display      = factor(display, levels = met_order15),
      group        = factor(group,   levels = names(grp_pal15)),
      pfas_dir_lab = factor(ifelse(pfas_dir == "DOWN", "PFAS negative", "PFAS positive"),
                            levels = c("PFAS negative", "PFAS positive"))
    )
}

bub15 <- bind_rows(
  make_bubble15(read.csv(here("out_files","spearman_taxa_metab_1m.csv")), "1-month"),
  make_bubble15(read.csv(here("out_files","spearman_taxa_metab_6m.csv")),         "6-month")
) %>% mutate(timepoint = factor(timepoint, levels = c("1-month","6-month")))

# Panel A plot (no legend)
p_A_v5 <- ggplot(bub15 %>% filter(n_sig > 0),
                 aes(x = group, y = display, fill = mean_rho, size = n_sig)) +
  geom_point(shape = 21, colour = "black", stroke = 0.4, alpha = 0.9) +
  geom_point(aes(colour = pfas_dir_lab), shape = NA) +
  scale_fill_gradient2(
    low = BUBBLE_NEG, mid = "white", high = BUBBLE_POS,
    midpoint = 0, limits = c(-0.8, 0.8),
    name = "Mean Spearman rho \n(within taxa)"
  ) +
  scale_colour_manual(
    values = c("PFAS negative" = LABEL_NEG, "PFAS positive" = LABEL_POS),
    guide  = "none"
  ) +
  scale_size_continuous(
    name   = "No. sig. species",
    range  = c(1, 9),
    breaks = c(1, 5, 10)
  ) +
  scale_y_discrete(
    limits = met_order15,
    expand = expansion(add = c(0.4, 0.4))
  ) +
  scale_x_discrete(expand = expansion(add = c(0.3, 0.3))) +
  facet_wrap(~ timepoint, nrow = 1) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(size = 12, angle = 45, hjust = 1,
                                    face = "bold.italic",
                                    colour = unname(grp_pal15[levels(bub15$group)])),
    axis.text.y      = element_text(size = 12, colour = met_col_vec, face = "bold"),
    strip.background = element_rect(fill = "black"),
    strip.text       = element_text(colour = "white", face = "bold", size = 12),
    panel.grid.major = element_line(colour = "grey88"),
    panel.grid.minor = element_blank(),
    legend.position  = "none"
  )

# ==============================================================================
# 16.  PANEL B: PATHWAY × METABOLITE BUBBLE  (1m vs 6m)
#      Fill = Spearman rho (individual pathway-metabolite pairs)
# ==============================================================================

story_pwy_map <- list(
  `HMO &\nGlycan` = c(
    "Bifidobacterium shunt",
    "L-rhamnose degradation I",
    "UDP-N-acetyl-D-glucosamine biosynthesis I",
    "sucrose degradation III",
    "sucrose degradation IV",
    "stachyose degradation"
  ),
  `Choline &\nMembrane` = c(
    "CDP-diacylglycerol biosynthesis I",
    "CDP-diacylglycerol biosynthesis II",
    "phosphatidylglycerol biosynthesis I"
  ),
  `Cholesterol &\nLipid` = c(
    "mevalonate pathway I",
    "palmitate biosynthesis",
    "fatty acid &beta;-oxidation VI"
  ),
  `Polyamine &\nNitrogen` = c(
    "L-lysine biosynthesis III",
    "urea cycle"
  ),
  `B-Vitamin\nCatabolism` = c(
    "NAD salvage pathway V",
    "folate transformations III"
  ),
  `SCFA &\nFermentation` = c(
    "acetyl-CoA fermentation to butanoate II",
    "superpathway of glycerol degradation to 1,3-propanediol",
    "isopropanol biosynthesis"
  ),
  `LPS &\nBiofilm` = c(
    "colanic acid building blocks biosynthesis"
  )
)

shorten_pwy16 <- function(x) {
  x %>%
    gsub("^Bifidobacterium shunt$",                    "Bifidobacterium shunt", .) %>%
    gsub("^L-rhamnose degradation I$",                 "L-Rhamnose deg. I", .) %>%
    gsub("UDP-N-acetyl-D-glucosamine biosynthesis I",  "UDP-GlcNAc biosyn. I", .) %>%
    gsub("^sucrose degradation III$",                  "Sucrose deg. III", .) %>%
    gsub("^sucrose degradation IV$",                   "Sucrose deg. IV", .) %>%
    gsub("^stachyose degradation$",                    "Stachyose deg.", .) %>%
    gsub("CDP-diacylglycerol biosynthesis I$",         "CDP-DAG biosyn. I", .) %>%
    gsub("CDP-diacylglycerol biosynthesis II",         "CDP-DAG biosyn. II", .) %>%
    gsub("phosphatidylglycerol biosynthesis I",        "PG biosyn. I", .) %>%
    gsub("mevalonate pathway I",                       "Mevalonate pathway I", .) %>%
    gsub("^palmitate biosynthesis$",                   "Palmitate biosyn.", .) %>%
    gsub("fatty acid &beta;-oxidation VI",             "FA b-oxidation VI", .) %>%
    gsub("L-lysine biosynthesis III",                  "L-Lys biosyn. III", .) %>%
    gsub("^urea cycle$",                               "Urea cycle", .) %>%
    gsub("NAD salvage pathway V",                      "NAD salvage V", .) %>%
    gsub("folate transformations III.*",               "Folate transform. III", .) %>%
    gsub("acetyl-CoA fermentation to butanoate II",    "Ac-CoA to butanoate II", .) %>%
    gsub("superpathway of glycerol degradation to 1,3-propanediol", "Glycerol to 1,3-PDO", .) %>%
    gsub("isopropanol biosynthesis.*",                 "Isopropanol biosyn.", .) %>%
    gsub("colanic acid building blocks biosynthesis",  "Colanic acid biosyn.", .)
}

all_story_pwy <- unique(unlist(story_pwy_map, use.names = FALSE))

pwy_axis_df <- data.frame(
  pathway_clean = all_story_pwy,
  met_axis      = rep(names(story_pwy_map),
                      sapply(story_pwy_map, length))[seq_along(all_story_pwy)],
  stringsAsFactors = FALSE
)
pwy_axis_df <- pwy_axis_df[!duplicated(pwy_axis_df$pathway_clean), ]

strip_pwy_id <- function(x) {
  x %>%
    sub("^[A-Z0-9][A-Z0-9.-]+-PWY[0-9A-Z.-]* ", "", .) %>%
    sub("^P[0-9]+-PWY ", "", .) %>%
    sub("^PWY[A-Z0-9]*-[0-9]+ ", "", .) %>%
    sub(" \\(.*\\)$", "", .)
}

met_axis_order16 <- c("LPS &\nBiofilm", "SCFA &\nFermentation",
                      "B-Vitamin\nCatabolism", "Polyamine &\nNitrogen",
                      "Cholesterol &\nLipid", "Choline &\nMembrane",
                      "HMO &\nGlycan")

# ── Biological category colours ────────────────────────────────────────────────
theme_colors <- c(
  "HMO metabolism"       = "#D9622B",
  "Membrane & Choline"   = "#3A7EC2",
  "Cholesterol & Lipid"  = "#B8860B",
  "Polyamine & Nitrogen" = "#8B4789",
  "B-Vitamin"            = "#2E8B57",
  "SCFA & Fermentation"  = "#B22222",
  "LPS & Biofilm"        = "#2F4F4F"
)

pwy_theme_map <- tribble(
  ~pwy_short,               ~theme,
  "Bifidobacterium shunt",  "HMO metabolism",
  "Stachyose deg.",         "HMO metabolism",
  "Sucrose deg. III",       "HMO metabolism",
  "Sucrose deg. IV",        "HMO metabolism",
  "UDP-GlcNAc biosyn. I",   "HMO metabolism",
  "L-Rhamnose deg. I",      "HMO metabolism",
  "CDP-DAG biosyn. I",      "Membrane & Choline",
  "CDP-DAG biosyn. II",     "Membrane & Choline",
  "PG biosyn. I",           "Membrane & Choline",
  "Mevalonate pathway I",   "Cholesterol & Lipid",
  "Palmitate biosyn.",      "Cholesterol & Lipid",
  "FA b-oxidation VI",      "Cholesterol & Lipid",
  "L-Lys biosyn. III",      "Polyamine & Nitrogen",
  "Urea cycle",             "Polyamine & Nitrogen",
  "NAD salvage V",          "B-Vitamin",
  "Folate transform. III",  "B-Vitamin",
  "Ac-CoA to butanoate II", "SCFA & Fermentation",
  "Glycerol to 1,3-PDO",    "SCFA & Fermentation",
  "Isopropanol biosyn.",    "SCFA & Fermentation",
  "Colanic acid biosyn.",   "LPS & Biofilm"
)

group_taxa16 <- function(genus) {
  case_when(
    genus == "Bifidobacterium"                      ~ "Bifidobacterium",
    genus %in% c("Blautia","Dorea","Roseburia",
                 "Coprococcus","Anaerobutyricum",
                 "Enterocloster","Ruminococcus",
                 "Eubacterium")                     ~ "Lachnospiraceae",
    genus == "Enterobacter"                         ~ "Enterobacter",
    TRUE                                            ~ "Other"
  )
}

# ── Build pathway-metabolite data ──────────────────────────────────────────────
build_pm16 <- function(fname, tp) {
  read.csv(here("out_files", fname)) %>%
    mutate(pathway_clean = strip_pwy_id(Pathway)) %>%
    filter(met_base %in% story_mets15$met_id,
           pathway_clean %in% all_story_pwy,
           p_val < 0.05) %>%
    left_join(story_mets15 %>% select(met_id, display, pfas_dir),
              by = c("met_base" = "met_id")) %>%
    left_join(pwy_axis_df, by = "pathway_clean") %>%
    mutate(pwy_short  = shorten_pwy16(pathway_clean),
           met_axis   = factor(met_axis, levels = met_axis_order16),
           display    = factor(display,  levels = met_order15),
           timepoint  = tp) %>%
    filter(!is.na(display)) %>%
    mutate(display = droplevels(display))
}

pm16 <- bind_rows(
  build_pm16("pathway_metabolite_spearman_1m.csv", "1-month"),
  build_pm16("pathway_metabolite_spearman_6m.csv", "6-month")
)

cat("Pathways NOT matched from story_pwy_map:\n")
print(all_story_pwy[!all_story_pwy %in% unique(pm16$pathway_clean)])

pm16 <- pm16 %>%
  group_by(timepoint, display) %>%
  mutate(n_pwy_corr = n()) %>%
  ungroup() %>%
  filter(n_pwy_corr >= 1) %>%
  select(-n_pwy_corr)

# ── Build taxa contribution data ───────────────────────────────────────────────
sp16_all <- pwy_species %>%
  filter(grepl("Continuous 1m PFAS + 1m Pathway", scenario, fixed = TRUE) |
           grepl("Continuous 1m PFAS + 6m Pathway", scenario, fixed = TRUE),
         pathway_clean %in% all_story_pwy, p_val < 0.05) %>%
  mutate(taxa_group = group_taxa16(genus),
         pwy_short  = shorten_pwy16(pathway_clean)) %>%
  left_join(pwy_axis_df, by = "pathway_clean") %>%
  mutate(met_axis = factor(met_axis, levels = met_axis_order16)) %>%
  group_by(pwy_short, Direction, taxa_group, met_axis) %>%
  summarise(n_species = n_distinct(species_name), .groups = "drop")

sp16_story <- sp16_all %>%
  filter(taxa_group != "Other") %>%
  mutate(taxa_group = factor(taxa_group,
                             levels = c("Bifidobacterium",
                                        "Lachnospiraceae",
                                        "Enterobacter")))

pm16 <- pm16 %>%
  filter(!is.na(pwy_short), !is.na(display)) %>%
  mutate(
    display   = factor(as.character(display),
                       levels = met_order15[met_order15 %in% as.character(display)]),
    timepoint = factor(timepoint,
                       levels = intersect(c("1-month","6-month"), unique(timepoint)))
  )

sp16_story <- sp16_story %>%
  filter(pwy_short %in% unique(pm16$pwy_short))

pathways_with_taxa <- unique(as.character(sp16_story$pwy_short))
pathways_with_taxa <- pathways_with_taxa[!is.na(pathways_with_taxa)]

if (length(pathways_with_taxa) == 0) {
  message(
    "\nNo pathways survived the story-pathway + taxa-overlap filter on this run.\n",
    "This final figure requires the specific PFAS<->Bifidobacterium/Lachnospiraceae/\n",
    "Enterobacter mechanistic-story pathways (hardcoded above from the manuscript's\n",
    "real findings) to be statistically significant -- essentially impossible to hit\n",
    "by chance on random synthetic data. This is expected on synthetic data, not an\n",
    "error; see code_submission_1/README.txt. Skipping the rest of script 14."
  )
  quit(save = "no", status = 0)
}

pm16       <- pm16 %>% filter(as.character(pwy_short) %in% pathways_with_taxa)
sp16_story <- sp16_story %>% filter(as.character(pwy_short) %in% pathways_with_taxa)

# ── PFAS direction per pathway (for y-axis label colors) ──────────────────────
pwy_dom_dir16 <- tribble(
  ~pwy_short,               ~pfas_pwy_dir,
  "Bifidobacterium shunt",  "PFAS negative",
  "Stachyose deg.",         "PFAS negative",
  "Sucrose deg. III",       "PFAS negative",
  "Sucrose deg. IV",        "PFAS negative",
  "UDP-GlcNAc biosyn. I",   "PFAS negative",
  "L-Rhamnose deg. I",      "PFAS negative",
  "CDP-DAG biosyn. I",      "PFAS positive",
  "CDP-DAG biosyn. II",     "PFAS positive",
  "PG biosyn. I",           "PFAS positive",
  "Palmitate biosyn.",      "PFAS negative",
  "FA b-oxidation VI",      "PFAS positive",
  "L-Lys biosyn. III",      "PFAS positive",
  "Urea cycle",             "PFAS positive",
  "NAD salvage V",          "PFAS positive",
  "Folate transform. III",  "PFAS positive",
  "Glycerol to 1,3-PDO",    "PFAS positive",
  "Colanic acid biosyn.",   "PFAS positive"
)

pwy_ord16 <- pm16 %>%
  group_by(pwy_short, met_axis) %>%
  summarise(mean_rho = mean(rho), .groups = "drop") %>%
  left_join(pwy_theme_map, by = "pwy_short") %>%
  mutate(met_axis = factor(met_axis, levels = met_axis_order16)) %>%
  arrange(met_axis, mean_rho) %>%
  pull(pwy_short) %>% unique() %>% rev()

pm16$pwy_short       <- factor(as.character(pm16$pwy_short), levels = pwy_ord16)
sp16_story$pwy_short <- factor(as.character(sp16_story$pwy_short), levels = pwy_ord16)

pwy_dir_for_col <- pwy_dom_dir16$pfas_pwy_dir[
  match(as.character(pwy_ord16), as.character(pwy_dom_dir16$pwy_short))
]
pwy_col_vec16 <- ifelse(
  is.na(pwy_dir_for_col), "grey40",
  ifelse(pwy_dir_for_col == "PFAS negative", "maroon4", "#3A7EC2")
)

# ── Metabolite HTML labels (colored by PFAS direction) ────────────────────────
met_html_labels <- setNames(
  paste0('<span style="color:', met_col_vec, '">', met_order15, '</span>'),
  met_order15
)

# ── Pathway HTML labels (colored by PFAS direction) ───────────────────────────
pwy_html_labels_v7 <- setNames(
  paste0('<span style="color:', pwy_col_vec16, '">',
         as.character(pwy_ord16), '</span>'),
  as.character(pwy_ord16)
)

# ── Biological category sidebar ────────────────────────────────────────────────
sidebar_df_v7 <- data.frame(
  pwy_short = factor(pwy_ord16, levels = pwy_ord16),
  cat       = pwy_theme_map$theme[match(pwy_ord16, pwy_theme_map$pwy_short)],
  x         = 1,
  stringsAsFactors = FALSE
) %>%
  mutate(cat = tidyr::replace_na(cat, "Other"))

p_sidebar_v7 <- ggplot(sidebar_df_v7,
                       aes(x = x, y = pwy_short, fill = cat)) +
  geom_tile(width = 0.9, height = 0.88, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(theme_colors, "Other" = "grey80"),
                    guide  = "none") +
  scale_y_discrete(limits = pwy_ord16,
                   expand = expansion(add = c(0.4, 0.4))) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(0, 1, 0, 1))

p_sidebar_v8_padded <- p_sidebar_v7 +
  theme(plot.margin = margin(t = 23, r = 1, b = 0, l = 8))

# ── Shared bubble theme ────────────────────────────────────────────────────────
bubble_theme_v7 <- theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_markdown(size = 11, angle = 45, hjust = 1,
                                        face = "bold"),
    axis.text.y      = element_markdown(size = 11, face = "bold"),
    strip.background = element_rect(fill = "black"),
    strip.text       = element_text(colour = "white", face = "bold", size = 12),
    panel.grid.major = element_line(colour = "grey88"),
    panel.grid.minor = element_blank(),
    legend.position  = "none"
  )

# ── Panel B: 1-month (only Neu5Ac, 4-Pyridoxate, Tryptamine) ─────────────────
mets_1m_show     <- c("Neu5Ac", "4-Pyridoxate", "Tryptamine")
met_html_labels_1m <- met_html_labels[mets_1m_show]

pm16_1m_v7 <- pm16 %>% filter(timepoint == "1-month")
pm16_1m_v7$timepoint <- factor("1-month")

p_B_1m_v7 <- ggplot(pm16_1m_v7,
                    aes(x = display, y = pwy_short,
                        fill = rho, size = abs(rho))) +
  geom_point(shape = 21, colour = "grey30", stroke = 0.4, alpha = 0.9) +
  scale_fill_gradient2(
    low      = BUBBLE_NEG, mid = "white", high = BUBBLE_POS,
    midpoint = 0, limits = c(-0.8, 0.8),
    name     = "Mean Spearman rho\n(within taxa)"
  ) +
  scale_size_continuous(range = c(1.5, 7), guide = "none") +
  scale_x_discrete(
    limits = mets_1m_show,
    labels = met_html_labels_1m,
    expand = expansion(add = c(0.5, 0.5))
  ) +
  scale_y_discrete(
    limits = pwy_ord16,
    labels = pwy_html_labels_v7,
    expand = expansion(add = c(0.4, 0.4))
  ) +
  facet_wrap(~ timepoint, nrow = 1) +
  labs(x = NULL, y = NULL) +
  bubble_theme_v7

# ── Panel B: 6-month (all metabolites, no y-axis labels) ─────────────────────
pm16_6m_v7 <- pm16 %>% filter(timepoint == "6-month")
n_met_6m16 <- n_distinct(pm16_6m_v7$display)

met_html_labels_6m <- met_html_labels[
  met_order15[met_order15 %in% as.character(pm16_6m_v7$display)]
]

p_B_6m_v7 <- ggplot(pm16_6m_v7,
                    aes(x = display, y = pwy_short,
                        fill = rho, size = abs(rho))) +
  geom_point(shape = 21, colour = "grey30", stroke = 0.4, alpha = 0.9) +
  scale_fill_gradient2(
    low      = TAXA_NEG_COL, mid = "white", high = TAXA_POS_COL,
    midpoint = 0, limits = c(-0.8, 0.8),
    name     = "Spearman \u03c1\n(pathway\u2013metabolite)"
  ) +
  scale_size_continuous(range = c(1.5, 7), guide = "none") +
  scale_x_discrete(
    limits = met_order15[met_order15 %in% as.character(pm16_6m_v7$display)],
    labels = met_html_labels_6m,
    expand = expansion(add = c(0.5, 0.5))
  ) +
  scale_y_discrete(
    limits = pwy_ord16,
    expand = expansion(add = c(0.4, 0.4))
  ) +
  facet_wrap(~ timepoint, nrow = 1) +
  labs(x = NULL, y = NULL) +
  bubble_theme_v7 +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

# Keep for legend extraction
p_B_6m_v7_for_legend <- p_B_6m_v7 + theme(legend.position = "right")

# ── Figure dimensions ──────────────────────────────────────────────────────────
w_panelA   <- 2 * 1.2 + 2.5
sidebar_w_v7 <- 0.25
w_1m_v7    <- length(mets_1m_show) * 0.85 + 0.5
w_6m_v7    <- n_met_6m16 * 0.85 + 0.5
h_plots    <- max(length(met_order15) * 0.40 + 2.5,
                  length(pwy_ord16) * 0.22 + 2.8)

# ── Panel A with left margin fix ──────────────────────────────────────────────
p_A_v9 <- p_A_v5 +
  theme(plot.margin = margin(t = 0, r = 20, b = 0, l = 15))

# ── egg::ggarrange for aligned assembly ───────────────────────────────────────
fig_panels_v9 <- egg::ggarrange(
  p_A_v9,
  p_sidebar_v8_padded,
  p_B_1m_v7,
  p_B_6m_v7,
  nrow   = 1,
  widths = c(w_panelA, sidebar_w_v7, w_1m_v7, w_6m_v7)
)

# ── Legends ────────────────────────────────────────────────────────────────────
leg_title_sz <- 13
leg_text_sz  <- 12
bar_w        <- unit(3.5, "cm")
bar_h        <- unit(0.6, "cm")

leg_base_theme <- theme(
  legend.position  = "bottom",
  legend.direction = "horizontal",
  legend.title     = element_text(size = leg_title_sz, face = "bold"),
  legend.text      = element_text(size = leg_text_sz)
)

leg_meanrho_v9 <- get_legend(
  p_A_v5 +
    scale_fill_gradient2(
      low      = BUBBLE_NEG, mid = "white", high = BUBBLE_POS,
      midpoint = 0, limits = c(-0.8, 0.8),
      name     = "Mean Spearman rho\n(within taxa)"
    ) +
    guides(
      fill   = guide_colorbar(direction      = "horizontal",
                              title.position = "top",
                              barwidth       = bar_w,
                              barheight      = bar_h),
      size   = "none",
      colour = "none"
    ) +
    leg_base_theme +
    theme(legend.position = "bottom")
)

leg_nsigspp_v9 <- get_legend(
  p_A_v5 +
    guides(
      fill   = "none",
      size   = guide_legend(title          = "No. sig. species",
                            title.position = "top"),
      colour = "none"
    ) +
    leg_base_theme +
  theme(legend.position = "bottom")
)

leg_dir_v9 <- get_legend(
  ggplot(
    data.frame(
      dir   = factor(c("PFAS negative", "PFAS positive"),
                     levels = c("PFAS negative", "PFAS positive")),
      dummy = 1
    ),
    aes(x = dummy, y = dummy, colour = dir)
  ) +
    geom_point(size = 5, shape = 15) +
    scale_colour_manual(
      name   = "PFAS\u2013Taxa/Pathway/Metabolite Direction",
      values = c("PFAS negative" = LABEL_NEG,
                 "PFAS positive" = LABEL_POS),
      guide  = guide_legend(
        nrow = 2, title.position = "top",
        override.aes = list(size   = 5, shape = 22,
                            fill   = c(LABEL_NEG, LABEL_POS),
                            colour = NA)
      )
    ) +
    theme_void() +
    leg_base_theme +
    theme(legend.key.spacing.y = unit(-2, "pt"))
)

leg_rho_v9 <- get_legend(
  p_B_6m_v7_for_legend +
    scale_fill_gradient2(
      low      = TAXA_NEG_COL, mid = "white", high = TAXA_POS_COL,
      midpoint = 0, limits = c(-0.8, 0.8),
      name = "Spearman rho\n(pathway-metabolite)"
    ) +
    guides(
      fill = guide_colorbar(direction      = "horizontal",
                            title.position = "top",
                            barwidth       = bar_w,
                            barheight      = bar_h)
    ) +
    leg_base_theme
)

leg_bio_v9 <- get_legend(
  ggplot(
    data.frame(
      cat   = factor(names(theme_colors), levels = names(theme_colors)),
      dummy = 1
    ),
    aes(x = dummy, y = dummy, fill = cat)
  ) +
    geom_tile() +
    scale_fill_manual(
      values = theme_colors,
      name   = "Biological Categories",
      guide  = guide_legend(ncol = 3, title.position = "top")
    ) +
    theme_void() +
    leg_base_theme +
    theme(legend.key.size = unit(0.55, "cm"))
)

legends_v9 <- plot_grid(
  leg_meanrho_v9,
  leg_nsigspp_v9,
  leg_dir_v9,
  leg_rho_v9,
  leg_bio_v9,
  nrow       = 1,
  rel_widths = c(0.9, 0.75, 1.1, 0.9, 1.8)
)

# ── Final combined figure (Fig 15 + 16) ───────────────────────────────────────
fig_v9 <- plot_grid(
  fig_panels_v9,
  legends_v9,
  ncol        = 1,
  rel_heights = c(h_plots, 1.5)
)

save_fig("14_combined_fig15_fig16_v9",
         w = w_panelA + sidebar_w_v7 + w_1m_v7 + w_6m_v7,
         h = h_plots + 2.2,
         { print(fig_v9) })

message("Saved: 14_combined_fig15_fig16_v9.pdf/png")

# ==============================================================================
# 17.  UNIFIED NETWORK FLOW: 1m (top) + 6m (bottom)
# ==============================================================================

X1 <- 1.5;  X2 <- 5.5;  X3 <- 10.5;  X4 <- 15.5
YB <- 7.5;  YL <- 3.0;  YC_1m <- 13.0

pfas_nd <- data.frame(
  label = c("PFOS", "PFOA", "PFNA", "PFHxS"),
  x     = X1,
  y     = c(13.0, 10.0, 7.0, 3.5),
  stringsAsFactors = FALSE
)

taxa_nd_1m <- data.frame(
  label = "Bifidobacterium\n (B. bifidum; depleted at 1m)",
  x = X2, y = YC_1m, col = COL_BIF, stringsAsFactors = FALSE
)

path_nd_1m <- data.frame(
  label = "Bifidobacterium shunt\nUDP-GlcNAc biosyn.\nSuccrose deg. III/IV",
  x = X3, y = YC_1m, col = COL_BIF, stringsAsFactors = FALSE
)

met_nd_1m <- data.frame(
  label = "Neu5Ac (Sialic acid)\n(HMO Metabolism)",
  x = X4, y = YC_1m, col = COL_BIF, stringsAsFactors = FALSE
)

taxa_nd <- data.frame(
  label = c(
    "Bifidobacterium spp.\n(B. longum, B. breve; depleted at 6m)",
    "Lachnospiraceae (Dorea, Blautia, Roseburia) &\nEnterobacter spp.\n(enriched at 6m)"
  ),
  x   = rep(X2, 2),
  y   = c(YB, YL),
  col = c(COL_BIF, COL_LAC),
  stringsAsFactors = FALSE
)

path_nd <- data.frame(
  label = c(
    "(Bifidobacterium shunt\nL-Rhamnose deg.\nGlcNAc/Neu5Ac deg.)",
    "Folate transform. III\nNAD salvage V\n FA B-oxidation VIn\nPalmitate Biosyn."
  ),
  x   = rep(X3, 2),
  y   = c(YB, YL),
  col = c(COL_BIF, COL_LAC),
  stringsAsFactors = FALSE
)

met_nd <- data.frame(
  label = c(
    "Fucose (HMO Metabolism)\nL-Rhamnose (HMO Metabolism)\nGlcNAc (HMO Metabolism\nNeu5Ac (HMO Metabolism)\nCholine (Acetylcholine Precursor)\nCholesterol (Contributes to Myelination)",
    "Tryptamine (Tryptophan shunting)\nCadaverine (Polyamine, Neurotoxic at high levels)\n4-Pyridoxate (Vit. B6 status marker)\nValerate (Gut barrier)"
  ),
  x   = rep(X4, 2),
  y   = c(YB, YL),
  col = c(COL_BIF, COL_LAC),
  stringsAsFactors = FALSE
)

edges_horiz <- data.frame(
  x    = c(X2, X3,   X2, X3),
  y    = c(YB, YB,   YL, YL),
  xend = c(X3, X4,   X3, X4),
  yend = c(YB, YB,   YL, YL),
  col  = c(COL_BIF, COL_BIF, COL_LAC, COL_LAC),
  lwd  = c(1.7, 1.7, 1.5, 1.5),
  stringsAsFactors = FALSE
)

edges_1m_horiz <- data.frame(
  x    = c(X2, X3),
  y    = rep(YC_1m, 2),
  xend = c(X3, X4),
  yend = rep(YC_1m, 2),
  col  = rep(COL_BIF, 2),
  lwd  = rep(1.7, 2),
  stringsAsFactors = FALSE
)

edges_pfas_1m <- data.frame(
  x    = rep(X1, 4),
  y    = pfas_nd$y,
  xend = rep(X2, 4),
  yend = rep(YC_1m, 4),
  col  = rep(COL_BIF, 4),
  lwd  = c(2.0, 1.8, 1.5, 1.5),
  stringsAsFactors = FALSE
)

edges_pfas_6m <- data.frame(
  x    = rep(X1, 6),
  y    = pfas_nd$y[c(1, 2, 3, 4, 1, 2)],
  xend = rep(X2, 6),
  yend = c(rep(YB, 4), rep(YL, 2)),
  col  = c(rep(COL_BIF, 4), rep(COL_LAC, 2)),
  lwd  = c(2.0, 1.7, 1.5, 1.4, 1.8, 1.6),
  stringsAsFactors = FALSE
)

transition_arrow <- data.frame(
  x = X4, y = YC_1m - 1.0,
  xend = X4, yend = YB + 1.2
)

p_network <- ggplot() +
  annotate("rect", xmin = 0.5, xmax = 17.5, ymin = 11.8, ymax = 14.5,
           fill = COL_BIF, alpha = 0.05) +
  annotate("rect", xmin = 0.5, xmax = 17.5, ymin = 5.5, ymax = 9.8,
           fill = COL_BIF, alpha = 0.05) +
  annotate("rect", xmin = 0.5, xmax = 17.5, ymin = 0.8, ymax = 4.8,
           fill = COL_LAC, alpha = 0.04) +
  annotate("segment",
           x = 0.6, xend = 17.5, y = 11.2, yend = 11.2,
           colour = "grey55", linewidth = 0.5, linetype = "dashed") +
  annotate("text", x = 9.0, y = 11.55,
           label = "developmental transition (1m \u2192 6m)",
           colour = "grey45", fontface = "italic", size = 3.8) +
  geom_segment(data = transition_arrow,
               aes(x = x, y = y, xend = xend, yend = yend),
               colour = "grey50", linewidth = 0.9, linetype = "dashed",
               arrow = arrow(length = unit(0.20, "cm"), type = "closed")) +
  annotate("text", x = X4 + 0.6, y = (YC_1m + YB) / 2,
           label = "Sialylated \u2192\nFucosylated\nHMO shift",
           colour = "grey40", fontface = "italic", size = 3.0, hjust = 0) +
  geom_curve(data = edges_pfas_1m,
             aes(x = x, y = y, xend = xend, yend = yend,
                 colour = col, linewidth = lwd),
             curvature = -0.12,
             arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
             alpha = 0.65, show.legend = FALSE) +
  geom_segment(data = edges_1m_horiz,
               aes(x = x, y = y, xend = xend, yend = yend,
                   colour = col, linewidth = lwd),
               arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
               alpha = 0.65, show.legend = FALSE) +
  geom_curve(data = edges_pfas_6m,
             aes(x = x, y = y, xend = xend, yend = yend,
                 colour = col, linewidth = lwd),
             curvature = 0.12,
             arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
             alpha = 0.65, show.legend = FALSE) +
  geom_segment(data = edges_horiz,
               aes(x = x, y = y, xend = xend, yend = yend,
                   colour = col, linewidth = lwd),
               arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
               alpha = 0.65, show.legend = FALSE) +
  scale_linewidth_identity() +
  scale_colour_identity() +
  geom_label(data = pfas_nd,
             aes(x = x, y = y, label = label),
             fill = "grey22", colour = "white", fontface = "bold",
             size = 4.5, label.padding = unit(0.30, "lines"),
             label.r = unit(0.18, "lines"), show.legend = FALSE) +
  geom_label(data = taxa_nd_1m,
             aes(x = x, y = y, label = label, fill = col),
             colour = "white", fontface = "bold", size = 4.0,
             label.padding = unit(0.40, "lines"),
             label.r = unit(0.12, "lines"), show.legend = FALSE) +
  geom_label(data = path_nd_1m,
             aes(x = x, y = y, label = label, fill = col),
             colour = "white", size = 3.5,
             label.padding = unit(0.45, "lines"),
             label.r = unit(0.10, "lines"), show.legend = FALSE) +
  geom_label(data = met_nd_1m,
             aes(x = x, y = y, label = label, fill = col),
             colour = "white", size = 3.8,
             label.padding = unit(0.45, "lines"),
             label.r = unit(0.10, "lines"), show.legend = FALSE) +
  geom_label(data = taxa_nd,
             aes(x = x, y = y, label = label, fill = col),
             colour = "white", fontface = "bold", size = 4.0,
             label.padding = unit(0.40, "lines"),
             label.r = unit(0.12, "lines"), show.legend = FALSE) +
  geom_label(data = path_nd,
             aes(x = x, y = y, label = label, fill = col),
             colour = "white", size = 3.5,
             label.padding = unit(0.48, "lines"),
             label.r = unit(0.10, "lines"), show.legend = FALSE) +
  geom_label(data = met_nd,
             aes(x = x, y = y, label = label, fill = col),
             colour = "white", size = 3.8,
             label.padding = unit(0.48, "lines"),
             label.r = unit(0.10, "lines"), show.legend = FALSE) +
  scale_fill_identity() +
  annotate("rect",
           xmin = c(0.6, 3.8, 8.5, 13.5),
           xmax = c(2.4, 7.2, 12.5, 17.4),
           ymin = 15.0, ymax = 15.8,
           fill = "grey22", alpha = 0.90) +
  annotate("text",
           x = c(X1, X2, X3, X4), y = 15.4,
           label = c("PFAS\n(1m breast milk)",
                     "Gut Microbiome\nTaxa",
                     "Gene Functional\nPathways",
                     "Fecal\nMetabolites"),
           fontface = "bold", size = 4.5, colour = "white", hjust = 0.5) +
  annotate("text", x = 17.8, y = YC_1m,
           label = "1\nmonth", angle = -90,
           fontface = "bold", colour = "grey40", size = 4.0) +
  annotate("text", x = 17.8, y = (YB + YL) / 2,
           label = "6\nmonths", angle = -90,
           fontface = "bold", colour = "grey40", size = 4.0) +
  annotate("segment",
           x = 0.6, xend = 1.5, y = 0.15, yend = 0.15,
           colour = COL_BIF, linewidth = 1.6) +
  annotate("text", x = 1.7, y = 0.15,
           label = "Negatively PFAS-associated (Bifidobacterium)",
           colour = COL_BIF, hjust = 0, size = 3.8) +
  annotate("segment",
           x = 8.5, xend = 9.4, y = 0.15, yend = 0.15,
           colour = COL_LAC, linewidth = 1.6) +
  annotate("text", x = 9.6, y = 0.15,
           label = "Positively PFAS-associated (Lachnospiraceae & Enterobacter)",
           colour = COL_LAC, hjust = 0, size = 3.8) +
  coord_cartesian(xlim = c(0.4, 18.2), ylim = c(0.0, 16.2)) +
  theme_void(base_size = 13) +
  theme(plot.margin = margin(5, 15, 10, 5))

save_fig("14_network_flow_unified", w = 18, h = 13,
         { print(p_network) })

message("===== Script 14 complete =====")
message("  14_combined_fig15_fig16_v9.pdf/png")
message("  14_network_flow_unified.pdf/png")