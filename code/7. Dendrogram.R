# header -----------------------------------------------------------------------
#
# TITLE:   7. Circular Phylogenetic Dendrogram with PFAS Results
#
# PURPOSE: 6 circular cladograms. then run each
#          scenario block in Section 8 individually to preview before saving.
#
# Code Review
# Ellie Holzhausen (EAH) on April 27, 2026
# Hoanan Li (HL) on
#
# set up -----------------------------------------------------------------------
rm(list = ls())

library(tidyverse)
library(dplyr)
library(stringr)
library(here)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(RColorBrewer)
library(scales)

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
mutate <- dplyr::mutate


# Load ALL 6 combined result files----------------------------------------------
pfas_order_continuous <- c("PFBS", "PFHxS", "PFNA", "PFOA", "PFOS")
pfas_order_binary     <- c("N.MeFOSAA", "PFBA", "PFDA", "PFDoA", "PFHpA")

results_files <- list(
  cont_1m1m      = "out_files/COMBINED_CLR_continuous_1m_1m.csv",
  cont_1m6m      = "out_files/COMBINED_CLR_continuous_1m_6m.csv",
  bin_1m1m       = "out_files/COMBINED_CLR_binary_1m_1m.csv",
  bin_1m6m       = "out_files/COMBINED_CLR_binary_1m_6m.csv",
  semiquant_1m1m = "out_files/COMBINED_CLR_semiquant_1m_1m.csv",
  semiquant_1m6m = "out_files/COMBINED_CLR_semiquant_1m_6m.csv",
  sens_cont      = "out_files/COMBINED_CLR_sensitivity_continuous_1m_6m.csv",
  sens_bin       = "out_files/COMBINED_CLR_sensitivity_binary_1m_6m.csv",
  sens_cont_bf   = "out_files/COMBINED_CLR_sensitivity_cont_1m_6m_bf.csv",
  sens_bin_bf    = "out_files/COMBINED_CLR_sensitivity_bin_1m_6m_bf.csv"
)

all_results <- lapply(results_files, function(f) read.csv(here::here(f)))

# Extract each result as a named object in the global environment
list2env(all_results, envir = .GlobalEnv)

# Species metadata from ALL scenarios for complete tree skeleton (~140 species)
species_meta_raw <- bind_rows(
  lapply(all_results, function(df)
    df %>% filter(taxa_level == "Species") %>%
      distinct(taxonomy_id, taxa_name, taxa_full)
  )
) %>%
  distinct(taxonomy_id, .keep_all = TRUE) %>%
  filter(!is.na(taxa_full), taxa_full != "")


# Diagnostic: inspect taxa_full format
cat(head(species_meta_raw$taxa_full, 3), sep = "\n")
cat(head(species_meta_raw$taxa_name, 3), sep = "\n")


# Convert taxa_full to double-dash Newick format--------------------
# FORMAT A = semicolon-separated

FORMAT <- "A"

convert_to_double_dash <- function(taxa_full_vec, format = "A") {
  unlist(lapply(taxa_full_vec, function(tf) {
    if (is.na(tf) || tf == "") return(NA_character_)
    if (format == "A") {
      parts <- str_split(tf, ";")[[1]] %>% trimws()
      # Split any part with embedded _x__ rank boundary
      # e.g. "g__Genus_s__Species name" -> two parts
      parts <- unlist(lapply(parts, function(p) {
        hits <- gregexpr("_[dpcofgs]__", p)[[1]]
        if (hits[1] == -1) return(p)
        cuts <- c(1, hits + 1)
        ends <- c(hits - 1, nchar(p))
        mapply(function(s, e) substr(p, s, e), cuts, ends)
      }))
      paste(parts, collapse = "--")
    } else {
      gsub("_(d|p|c|o|f|g|s)__", "--\\1__", tf)
    }
  }))
}

species_meta <- species_meta_raw %>%
  mutate(
    taxa_full_dd = convert_to_double_dash(taxa_full, format = FORMAT),
    display_name = taxa_name %>% gsub("_", " ", .) %>% gsub("\\[|\\]", "", .)
  ) %>%
  filter(!is.na(taxa_full_dd)) %>%
  # Fix brackets in taxa_full_dd (e.g. s__[Clostridium]_innocuum)
  mutate(taxa_full_dd = gsub("\\[|\\]", "", taxa_full_dd))

# Check how it looks
cat(head(species_meta$taxa_full_dd, 3), sep = "\n")

dd_count <- str_count(species_meta$taxa_full_dd, "--")
print(table(dd_count))
if (any(dd_count != 6)) {
  cat("WARNING — species with wrong rank count:\n")
  print(species_meta$taxa_full_dd[dd_count != 6])
}

# Check if removing brackets creates any duplicate taxa names
before_removal <- species_meta_raw %>% 
  filter(grepl("\\[|\\]", taxa_name)) %>%
  pull(taxa_name)

after_removal <- gsub("\\[|\\]", "", before_removal)

cat("Taxa with brackets:", length(before_removal), "\n")
cat("Duplicates created by removal:", 
    sum(duplicated(after_removal)), "\n")
print(data.frame(original = before_removal, 
                 cleaned  = after_removal))

# Build taxonomy rank name vectors and phylo_level------------------------------
species_names_full <- species_meta$taxa_full_dd

extract_rank_names <- function(full_names, stop_before) {
  unique(sapply(full_names, function(x) {
    parts <- str_split(x, "--")[[1]]
    idx   <- grep(paste0("^", stop_before, "__"), parts)
    if (length(idx) == 0) return(NA_character_)
    paste(parts[1:(idx - 1)], collapse = "--")
  })) %>% na.omit() %>% as.character()
}

genus_names  <- extract_rank_names(species_names_full, "s")
family_names <- extract_rank_names(species_names_full, "g")
order_names  <- extract_rank_names(species_names_full, "f")
class_names  <- extract_rank_names(species_names_full, "o")
phylum_names <- extract_rank_names(species_names_full, "c")

names_all   <- c(species_names_full, phylum_names, class_names,
                 order_names, family_names, genus_names)
phylo_level <- sapply(str_split(names_all, "--"), length)

cat("\nTaxonomy level counts:\n")
print(table(phylo_level))


# Build Newick tree bottom-up---------------------------------------
l1 <- list()
for (n in 7:3) {
  l2       <- list()
  parents  <- names_all[phylo_level == (n - 1)]
  children <- names_all[phylo_level == n]
  for (m in parents) {
    matched <- children[grep(paste0("^", m, "--"), children, fixed = FALSE)]
    if (length(matched) == 0) next
    tips <- if (n == 7) matched else { t <- l1[matched]; t[!sapply(t, is.null)] }
    if (length(tips) > 0)
      l2[[m]] <- paste0("(", paste(tips, collapse = ","), ")", m)
  }
  l1 <- l2
  cat("Level", n, ":", length(l1), "nodes\n")
}

# Strip node labels from Newick string — long lineage strings break ape parser
l3       <- paste0("(", paste(l1, collapse = ","), ")Bacteria;")
l3_clean <- gsub("\\)([^,;()]+)", ")", l3)
l3_clean <- paste0(substr(l3_clean, 1, nchar(l3_clean) - 1), "Bacteria;")

tree <- ape::read.tree(text = l3_clean)
tree$node.label <- NULL

cat("\nTree: Tips =", ape::Ntip(tree), "| Nodes =", ape::Nnode(tree), "\n")
if (ape::Ntip(tree) < 10) stop("Tree parsing failed — check l3_clean")


# Map tip labels to species display names---------------------------------------
tip_lookup <- species_meta %>%
  select(taxa_full_dd, display_name, taxonomy_id) %>%
  distinct()

tree_display <- tree
tree_display$tip.label <- ifelse(
  tree_display$tip.label %in% tip_lookup$taxa_full_dd,
  tip_lookup$display_name[match(tree_display$tip.label, tip_lookup$taxa_full_dd)],
  tree_display$tip.label
)
tree_display$node.label <- NULL

cat("First 6 tip labels:\n")
print(head(tree_display$tip.label, 6))


# Preview with tip labels (topology check only) # Warning here is okay
p_preview <- ggtree(tree_display, layout = "fan", open.angle = 20,
                    branch.length = "none", size = 0.40, color = "black") +
  geom_tiplab(size = 5, fontface = "bold.italic", offset = 0.3) +
  theme(plot.background  = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

print(p_preview)

# p_preview <- ggtree(tree_display, layout = "fan", open.angle = 20,
#                     branch.length = "none", size = 0.40, color = "black") +
#   geom_tiplab(size = 7, fontface = "bold.italic", offset = 0.2) +  # bigger text, closer to tree
#   theme(plot.background  = element_rect(fill = "white", color = NA),
#         panel.background = element_rect(fill = "white", color = NA)) +
#   xlim(0, 8)  # compress horizontal space — reduces tree width
# 
# print(p_preview)
# ggsave(here::here("out_figures", "dendrogram_preview_tree_only.pdf"),
#        plot = p_preview, width = 20, height = 22,  # reduced from 25x28
#        device = cairo_pdf, units = "in")


# Extract tree layout data (tip y positions + node x/y)-------------------------
# Warning about deprecated aesthetic mapping is expected ggtree behavior
# and does not affect the output — tree_layout_full is created correctly
tree_layout_full <- ggtree(tree_display, layout = "fan",
                           branch.length = "none")$data

x_max <- max(tree_layout_full$x[tree_layout_full$isTip], na.rm = TRUE)

# Tip y positions for symbol overlay
tip_y_pos <- tree_layout_full %>%
  filter(isTip) %>%
  select(display_name = label, y_tip = y)

# All node positions (both tips and internal) for mixture dots
node_xy <- tree_layout_full %>%
  select(node, x, y) %>%
  filter(!is.na(node))

cat("x_max:", round(x_max, 3), "\n")
cat("Total nodes in layout:", nrow(node_xy), "\n")


# Phylum MRCA nodes and p_base--------------------------------------------------
# NO text/labels on p_base — final plots are clean of all text
node_num <- rep(NA, length(phylum_names))
tip_num  <- rep(NA, length(phylum_names))

for (i in seq_along(phylum_names)) {
  ph_spp <- species_meta %>%
    filter(str_detect(taxa_full_dd, fixed(phylum_names[i]))) %>%
    pull(display_name)
  tips_m <- tree_display$tip.label[tree_display$tip.label %in% ph_spp]
  
  if (length(tips_m) == 0) {
    node_num[i] <- NA; tip_num[i] <- NA
  } else if (length(tips_m) == 1) {
    node_num[i] <- which(tree_display$tip.label == tips_m); tip_num[i] <- 1
  } else {
    nd <- try(ggtree::MRCA(tree_display, tips_m), silent = TRUE)
    if (inherits(nd, "try-error") || is.null(nd)) {
      node_num[i] <- NA; tip_num[i] <- NA
    } else {
      node_num[i] <- nd; tip_num[i] <- length(tips_m)
    }
  }
}

dd3 <- data.frame(phylum_names, node_num, tip_num) %>%
  mutate(node_num = as.numeric(node_num),
         Phylum   = str_extract(phylum_names, "p__[^-]+")) %>%
  filter(!is.na(node_num)) %>%
  arrange(node_num)

cat("Phyla with MRCA nodes:", nrow(dd3), "\n")
print(dd3 %>% select(Phylum, tip_num))

n_phyla     <- nrow(dd3)


phylum_cols <- setNames(
  RColorBrewer::brewer.pal(max(3, min(n_phyla, 12)), "Set3")[seq_len(n_phyla)],
  dd3$Phylum
)

# Base tree: NO text of any kind (Phylum shading commented out for now)---------
p_base <- ggtree(tree_display, layout = "fan", open.angle = 20,
                 branch.length = "none", size = 0.7, color = "grey30")
# Warning here is expected ggtree behavior (p_base is created correctly)

# for (i in seq_len(nrow(dd3))) {
#   p_base <- p_base +
#     geom_hilight(node  = dd3$node_num[i],
#                  fill  = phylum_cols[dd3$Phylum[i]],
#                  alpha = 0.15, extend = 0.5)
#   show.legend = FALSE
# }

# Phylum color legend via invisible points (no text on tree)
# phylum_legend_df <- dd3 %>%
#   mutate(label_clean = str_replace(Phylum, "p__", ""))
# 
# p_base <- p_base +
#   geom_point(data        = phylum_legend_df,
#              mapping     = aes(x = -Inf, y = -Inf, fill = label_clean),
#              shape       = 22, size = 8,
#              inherit.aes = FALSE, alpha = 0) +
#   scale_fill_manual(
#     values = setNames(phylum_cols, str_replace(names(phylum_cols), "p__", "")),
#     name   = "Phylum",
#     guide  = guide_legend(order = 0,
#                           override.aes = list(alpha = 0.6,
#                                               size   = 8,
#                                               shape  = 22))
#   )
  

# Helper functions--------------------------------------------------------------

# Get MRCA node — extracts actual rank+name from higher taxa taxa_full
# Higher taxa format: "d__X;p__X;...;g___g__Actinomyces" (last seg = g___g__Name)
get_mrca_node <- function(taxa_full_raw, taxa_level, species_meta, tree_display) {
  
  # Extract actual rank identifier from last semicolon segment
  # e.g. "g___g__Actinomyces"   -> "g__Actinomyces"
  #      "g___f__Micrococcaceae" -> "f__Micrococcaceae"
  #      "g___o__Nostocales"     -> "o__Nostocales"
  #      "g___c__Alphaproteobacteria" -> "c__Alphaproteobacteria"
  last_seg <- trimws(tail(str_split(taxa_full_raw, ";")[[1]], 1))
  actual_id <- str_match(last_seg, "[a-z]___([a-z]__[^;]+)")[1, 2]
  
  if (is.na(actual_id)) {
    cat("   Could not extract actual_id from last_seg:", last_seg, "\n")
    return(NA_integer_)
  }
  actual_id <- trimws(actual_id)
  
  # Find all species in species_meta whose lineage contains this rank+name
  member_species <- species_meta %>%
    filter(str_detect(taxa_full_dd, fixed(actual_id))) %>%
    pull(display_name)
  
  tips_in_tree <- tree_display$tip.label[tree_display$tip.label %in% member_species]
  
  if (length(tips_in_tree) == 0) return(NA_integer_)
  if (length(tips_in_tree) == 1)
    return(which(tree_display$tip.label == tips_in_tree))
  
  nd <- try(ggtree::MRCA(tree_display, tips_in_tree), silent = TRUE)
  if (inherits(nd, "try-error") || is.null(nd)) return(NA_integer_)
  return(nd)
}

# Mixture / n_detect dots — ALL taxa levels
prep_dot_nodes <- function(results_df, predictor_name,
                           tip_lookup, species_meta,
                           tree_display, node_xy, FORMAT) {
  
  empty <- data.frame(node = integer(), x = numeric(), y = numeric(),
                      Direction = character(), sig_type = character(),
                      taxa_level = character(), stringsAsFactors = FALSE)
  
  dat <- results_df %>%
    filter(predictor == predictor_name, !is.na(p_val), p_val < 0.05)
  
  if (nrow(dat) == 0) {
    cat("   No p<0.05 results for:", predictor_name, "\n")
    return(empty)
  }
  
  cat("   p<0.05 by taxa_level for", predictor_name, ":\n")
  print(table(dat$taxa_level))
  
  dat <- dat %>%
    mutate(
      betas     = as.numeric(betas),
      Beta_IQR  = if ("Beta_IQR" %in% names(dat)) as.numeric(Beta_IQR) else NA_real_,
      plot_beta = dplyr::coalesce(Beta_IQR, betas),
      Direction = case_when(plot_beta > 0 ~ "Positive",
                            plot_beta < 0 ~ "Negative",
                            TRUE          ~ NA_character_),
      sig_type  = ifelse(!is.na(FDR) & FDR < 0.05, "FDR", "nominal")
    )
  
  rows_out <- lapply(seq_len(nrow(dat)), function(i) {
    row <- dat[i, ]
    
    if (row$taxa_level == "Species") {
      # Match species via display_name -> tip label
      dn  <- tip_lookup$display_name[tip_lookup$taxonomy_id == row$taxonomy_id]
      if (length(dn) == 0 || is.na(dn[1])) return(NULL)
      nid <- which(tree_display$tip.label == dn[1])
    } else {
      # Higher taxa: extract rank+name from raw taxa_full last segment
      nid <- get_mrca_node(row$taxa_full, row$taxa_level, species_meta, tree_display)
    }
    
    if (length(nid) == 0 || is.na(nid)) return(NULL)
    xy <- node_xy %>% filter(node == nid)
    if (nrow(xy) == 0) return(NULL)
    
    data.frame(node       = nid,
               x          = xy$x[1],
               y          = xy$y[1],
               Direction  = row$Direction,
               sig_type   = row$sig_type,
               taxa_level = row$taxa_level,
               stringsAsFactors = FALSE)
  })
  
  result <- bind_rows(rows_out)
  if (nrow(result) == 0) return(empty)
  
  cat("   Dot nodes successfully placed by level:\n")
  print(table(result$taxa_level))
  result %>% distinct()
}

# Wide beta matrix for gheatmap (species level, ALL betas)
make_pfas_wide <- function(results_df, pfas_vars, tip_lookup) {
  results_df %>%
    filter(taxa_level == "Species",
           predictor %in% pfas_vars,
           taxonomy_id %in% tip_lookup$taxonomy_id) %>%
    left_join(tip_lookup %>% select(taxonomy_id, display_name),
              by = "taxonomy_id") %>%
    mutate(
      predictor = factor(predictor, levels = pfas_vars),
      betas     = as.numeric(betas),
      Beta_IQR  = if ("Beta_IQR" %in% names(.)) as.numeric(Beta_IQR) else NA_real_,
      plot_beta = dplyr::coalesce(Beta_IQR, betas)
    ) %>%
    select(display_name, predictor, plot_beta) %>%
    pivot_wider(names_from = predictor, values_from = plot_beta) %>%
    as.data.frame() %>%
    { rownames(.) <- .$display_name; . } %>%
    select(-display_name)
}

# Symbol matrix — * for ALL p<0.05 (includes FDR<0.05 automatically)
make_pfas_annot_wide <- function(results_df, pfas_vars, tip_lookup) {
  results_df %>%
    filter(taxa_level == "Species",
           predictor %in% pfas_vars,
           taxonomy_id %in% tip_lookup$taxonomy_id) %>%
    left_join(tip_lookup %>% select(taxonomy_id, display_name),
              by = "taxonomy_id") %>%
    mutate(
      predictor = factor(predictor, levels = pfas_vars),
      symbol    = ifelse(!is.na(p_val) & p_val < 0.05, "*", NA_character_)
    ) %>%
    filter(!is.na(symbol)) %>%
    select(display_name, predictor, symbol) %>%
    pivot_wider(names_from = predictor, values_from = symbol) %>%
    as.data.frame() %>%
    { rownames(.) <- .$display_name; . } %>%
    select(-display_name)
}

# Master plot builder-----------------------------------------------------------
make_dendro_plot <- function(results_df,
                             pfas_vars,
                             is_binary,
                             color_scheme    = NULL,
                             p_base,
                             tip_lookup,
                             species_meta,
                             tree_display,
                             node_xy,
                             tip_y_pos,
                             x_max,
                             FORMAT,
                             col_width       = 0.06,
                             col_gap         = 0.07,
                             starting_offset = 0.3) {
  
  mix_colors <- c("Positive" = "#3333FF", "Negative" = "red")
  
  # Detect correct mixture predictor — "Mixture" for continuous, "N-detect" for binary
  mix_predictor <- ifelse("Mixture" %in% results_df$predictor, "Mixture", "N-detect")
  cat("  Mixture predictor:", mix_predictor, "\n")
  
  mix_nodes <- prep_dot_nodes(results_df, mix_predictor,
                              tip_lookup, species_meta,
                              tree_display, node_xy, FORMAT)
  cat("  Mixture nodes p<0.05:", nrow(mix_nodes),
      "| FDR:", sum(mix_nodes$sig_type == "FDR", na.rm = TRUE), "\n")
  if (nrow(mix_nodes) > 0) print(table(mix_nodes$taxa_level))
  
  if (nrow(mix_nodes) > 0) {
    p <- p_base %<+% mix_nodes +
      geom_point(aes(x = x, y = y, color = Direction),
                 shape = 16, size = 7.0, alpha = 0.95, na.rm = TRUE) +
      scale_colour_manual(values       = mix_colors,
                          na.translate = FALSE,
                          name         = "PFAS Mixture",
                          guide        = guide_legend(order = 1,
                                                      override.aes = list(size = 6, shape = 16)))
  } else {
    p <- p_base +
      scale_colour_manual(values       = mix_colors,
                          na.translate = FALSE,
                          name         = "PFAS Mixture",
                          guide        = guide_legend(order = 1,
                                                      override.aes = list(size = 6, shape = 16)))
  }
  
  # Individual PFAS heatmap
  pfas_wide       <- make_pfas_wide(results_df, pfas_vars, tip_lookup)
  pfas_annot_wide <- make_pfas_annot_wide(results_df, pfas_vars, tip_lookup)
  
  scale_limit <- ceiling(max(abs(as.vector(as.matrix(pfas_wide))),
                             na.rm = TRUE) * 10) / 10
  cat("  Heatmap scale: ±", scale_limit, "\n")
  
  col_offset <- starting_offset
  
  # Force-clear any discrete fill before gheatmap
  p <- p + new_scale_fill()
  
  # Determine heatmap colors based on color_scheme / is_binary
  heat_low  <- case_when(
    !is.null(color_scheme) && color_scheme == "semiquant" ~ "#8B008B",
    is_binary                                             ~ "#D55E00",
    TRUE                                                  ~ "orange"
  )
  heat_high <- case_when(
    !is.null(color_scheme) && color_scheme == "semiquant" ~ "#0072B2",
    is_binary                                             ~ "#009E73",
    TRUE                                                  ~ "#3333FF"
  )
  
  for (pfas in pfas_vars) {
    col_df <- pfas_wide %>% select(all_of(pfas))
    p <- suppressMessages(
      gheatmap(p, col_df,
               offset            = col_offset,
               width             = col_width,
               colnames_offset_y = -3,
               colnames_angle    = 45,
               hjust             = 0,
               font.size         = 2.5,
               color             = "grey8") +
        scale_fill_gradient2(
          low      = heat_low,
          mid      = "white",
          high     = heat_high,
          midpoint = 0, na.value = "grey92",
          limits   = c(-scale_limit, scale_limit),
          name     = ifelse(is_binary, 
                            "Individual PFAS Estimate",
                            "Individual PFAS Estimate (IQR scaled)"),
          guide    = guide_colorbar(order = 4,
                                    barheight = unit(1.2, "in"),
                                    barwidth  = unit(0.25, "in"))
        )
    )
    col_offset <- col_offset + col_width * x_max + col_gap
  }
  
  # Extract actual cell x-centers from gheatmap tile layers
  pfas_x_centers <- list()
  for (lyr in p$layers) {
    ld <- tryCatch(lyr$data, error = function(e) NULL)
    if (!is.null(ld) && is.data.frame(ld) &&
        "variable" %in% names(ld) && "x" %in% names(ld) && nrow(ld) > 0) {
      var_vals <- unique(as.character(ld$variable))
      for (v in var_vals) {
        if (v %in% pfas_vars && !v %in% names(pfas_x_centers)) {
          pfas_x_centers[[v]] <- mean(ld$x[ld$variable == v], na.rm = TRUE)
        }
      }
    }
  }
  cat("  Heatmap x-centers extracted:", length(pfas_x_centers), "columns\n")
  
  # Symbol overlay: * = p<0.05, centered in cells
  if (nrow(pfas_annot_wide) > 0 && length(pfas_x_centers) > 0) {
    sym_long <- pfas_annot_wide %>%
      tibble::rownames_to_column("display_name") %>%
      pivot_longer(-display_name, names_to = "predictor", values_to = "symbol") %>%
      filter(!is.na(symbol)) %>%
      mutate(x_center = unlist(pfas_x_centers)[predictor]) %>%
      left_join(tip_y_pos, by = "display_name") %>%
      filter(!is.na(x_center), !is.na(y_tip))
    
    if (nrow(sym_long) > 0) {
      p <- p + geom_text(data        = sym_long,
                         mapping     = aes(x = x_center, y = y_tip, label = symbol),
                         size        = 7.5,
                         fontface    = "bold",
                         color       = "grey10",
                         vjust       = 0.6,
                         hjust       = 0.5,
                         inherit.aes = FALSE)
    }
  }
  
  p + theme(
    legend.position  = "right",
    legend.key.size  = unit(0.9, "cm"),
    legend.text      = element_text(size = 14),
    legend.title     = element_text(size = 16, face = "bold"),
    legend.box       = "vertical",
    legend.spacing.y = unit(0.3, "cm"),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
}

n_tips   <- ape::Ntip(tree_display)
fig_size <- max(14, n_tips * 0.10)


# Preview the plots and save----------------------------------------------------

# PLOT 1: Continuous PFAS | 1m breastmilk → 1m microbiome
p1 <- make_dendro_plot(
  results_df = all_results[["cont_1m1m"]], pfas_vars = pfas_order_continuous,
  is_binary = FALSE, p_base = p_base, tip_lookup = tip_lookup,
  species_meta = species_meta, tree_display = tree_display,
  node_xy = node_xy, tip_y_pos = tip_y_pos, x_max = x_max, FORMAT = FORMAT
)
print(p1) # we will see warning messages here. we can safely ignore them
ggsave(here::here("out_figures","dendrogram_continuous_1m_1m.pdf"), plot=p1,
       width=fig_size+5, height=fig_size*1.5, device=cairo_pdf, units="in")
ggsave(here::here("out_figures", "dendrogram_continuous_1m_1m.png"), plot = p1,
       width = fig_size + 5, height = fig_size * 1.5, dpi = 600, units = "in")


# PLOT 2: Continuous PFAS | 1m breastmilk → 6m microbiome
p2 <- make_dendro_plot(
  results_df = all_results[["cont_1m6m"]], pfas_vars = pfas_order_continuous,
  is_binary = FALSE, p_base = p_base, tip_lookup = tip_lookup,
  species_meta = species_meta, tree_display = tree_display,
  node_xy = node_xy, tip_y_pos = tip_y_pos, x_max = x_max, FORMAT = FORMAT
)
print(p2)
ggsave(here::here("out_figures","dendrogram_continuous_1m_6m.pdf"), plot=p2,
       width=fig_size+5, height=fig_size, device=cairo_pdf, units="in")
ggsave(here::here("out_figures", "dendrogram_continuous_1m_6m.png"), plot = p2,
       width = fig_size + 5, height = fig_size, dpi = 600, units = "in")


# PLOT 3: Binary PFAS | 1m breastmilk → 1m microbiome
p3 <- make_dendro_plot(
  results_df = all_results[["bin_1m1m"]], pfas_vars = pfas_order_binary,
  is_binary = TRUE, p_base = p_base, tip_lookup = tip_lookup,
  species_meta = species_meta, tree_display = tree_display,
  node_xy = node_xy, tip_y_pos = tip_y_pos, x_max = x_max, FORMAT = FORMAT
)
print(p3)
ggsave(here::here("out_figures","dendrogram_binary_1m_1m.pdf"), plot=p3,
       width=fig_size+5, height=fig_size, device=cairo_pdf, units="in")
ggsave(here::here("out_figures", "dendrogram_binary_1m_1m.png"), plot = p3,
       width = fig_size + 5, height = fig_size, dpi = 600, units = "in")

# PLOT 4: Binary PFAS | 1m breastmilk → 6m microbiome
p4 <- make_dendro_plot(
  results_df = all_results[["bin_1m6m"]], pfas_vars = pfas_order_binary,
  is_binary = TRUE, p_base = p_base, tip_lookup = tip_lookup,
  species_meta = species_meta, tree_display = tree_display,
  node_xy = node_xy, tip_y_pos = tip_y_pos, x_max = x_max, FORMAT = FORMAT
)
print(p4)
ggsave(here::here("out_figures","dendrogram_binary_1m_6m.pdf"), plot=p4,
       width=fig_size+5, height=fig_size, device=cairo_pdf, units="in")
ggsave(here::here("out_figures", "dendrogram_binary_1m_6m.png"), plot = p4,
       width = fig_size + 5, height = fig_size, dpi = 600, units = "in")

# PLOT 5: Semi-quantitative PFAS | 1m breastmilk → 1m microbiome
p5 <- make_dendro_plot(
  results_df   = all_results[["semiquant_1m1m"]], pfas_vars = pfas_order_binary,
  is_binary    = FALSE, color_scheme = "semiquant", p_base = p_base,
  tip_lookup   = tip_lookup, species_meta = species_meta,
  tree_display = tree_display, node_xy = node_xy, tip_y_pos = tip_y_pos,
  x_max = x_max, FORMAT = FORMAT
)
print(p5)
ggsave(here::here("out_figures", "dendrogram_semiquant_1m_1m.pdf"), plot = p5,
       width = fig_size + 5, height = fig_size, device = cairo_pdf, units = "in")
ggsave(here::here("out_figures", "dendrogram_semiquant_1m_1m.png"), plot = p5,
       width = fig_size + 5, height = fig_size, dpi = 600, units = "in")

# PLOT 6: Semi-quantitative PFAS | 1m breastmilk → 6m microbiome
p6 <- make_dendro_plot(
  results_df   = all_results[["semiquant_1m6m"]], pfas_vars = pfas_order_binary,
  is_binary    = FALSE, color_scheme = "semiquant", p_base = p_base,
  tip_lookup   = tip_lookup, species_meta = species_meta,
  tree_display = tree_display, node_xy = node_xy, tip_y_pos = tip_y_pos,
  x_max = x_max, FORMAT = FORMAT
)
print(p6)
ggsave(here::here("out_figures", "dendrogram_semiquant_1m_6m.pdf"), plot = p6,
       width = fig_size + 5, height = fig_size, device = cairo_pdf, units = "in")
ggsave(here::here("out_figures", "dendrogram_semiquant_1m_6m.png"), plot = p6,
       width = fig_size + 5, height = fig_size, dpi = 600, units = "in")

# PLOT 7: Sensitivity | Continuous PFAS | 1m → 6m
p7 <- make_dendro_plot(
  results_df = all_results[["sens_cont"]], pfas_vars = pfas_order_continuous,
  is_binary = FALSE, p_base = p_base, tip_lookup = tip_lookup,
  species_meta = species_meta, tree_display = tree_display,
  node_xy = node_xy, tip_y_pos = tip_y_pos, x_max = x_max, FORMAT = FORMAT
)
print(p7)
ggsave(here::here("out_figures","dendrogram_sensitivity_continuous_1m_6m.pdf"), plot=p7,
       width=fig_size+5, height=fig_size, device=cairo_pdf, units="in")
ggsave(here::here("out_figures", "dendrogram_sensitivity_continuous_1m_6m.png"), plot = p7,
       width = fig_size + 5, height = fig_size, dpi = 600, units = "in")

# PLOT 8: Sensitivity | Binary PFAS | 1m → 6m
p8 <- make_dendro_plot(
  results_df = all_results[["sens_bin"]], pfas_vars = pfas_order_binary,
  is_binary = TRUE, p_base = p_base, tip_lookup = tip_lookup,
  species_meta = species_meta, tree_display = tree_display,
  node_xy = node_xy, tip_y_pos = tip_y_pos, x_max = x_max, FORMAT = FORMAT
)
print(p8)
ggsave(here::here("out_figures","dendrogram_sensitivity_binary_1m_6m.pdf"), plot=p8,
       width=fig_size+5, height=fig_size, device=cairo_pdf, units="in")
ggsave(here::here("out_figures", "dendrogram_sensitivity_binary_1m_6m.png"), plot = p8,
       width = fig_size + 5, height = fig_size, dpi = 600, units = "in")

# We dont need to run below Codes as I will have beta correlation plot for below scenarios and We decided to not include individual dendrograms plots for now
# # PLOT 9: Sensitivity BF | Continuous PFAS | 1m → 6m (BF ≥1/day at 6m)
# p9 <- make_dendro_plot(
#   results_df = all_results[["sens_cont_bf"]], pfas_vars = pfas_order_continuous,
#   is_binary = FALSE, p_base = p_base, tip_lookup = tip_lookup,
#   species_meta = species_meta, tree_display = tree_display,
#   node_xy = node_xy, tip_y_pos = tip_y_pos, x_max = x_max, FORMAT = FORMAT
# )
# print(p9)
# ggsave(here::here("out_figures", "dendrogram_sensitivity_cont_1m_6m_bf.pdf"), plot = p9,
#        width = fig_size + 5, height = fig_size, device = cairo_pdf, units = "in")
# ggsave(here::here("out_figures", "dendrogram_sensitivity_cont_1m_6m_bf.png"), plot = p9,
#        width = fig_size + 5, height = fig_size, dpi = 600, units = "in")
# 
# # PLOT 10: Sensitivity BF | Binary PFAS | 1m → 6m (BF ≥1/day at 6m)
# p10 <- make_dendro_plot(
#   results_df = all_results[["sens_bin_bf"]], pfas_vars = pfas_order_binary,
#   is_binary = TRUE, p_base = p_base, tip_lookup = tip_lookup,
#   species_meta = species_meta, tree_display = tree_display,
#   node_xy = node_xy, tip_y_pos = tip_y_pos, x_max = x_max, FORMAT = FORMAT
# )
# print(p10)
# ggsave(here::here("out_figures", "dendrogram_sensitivity_bin_1m_6m_bf.pdf"), plot = p10,
#        width = fig_size + 5, height = fig_size, device = cairo_pdf, units = "in")
# ggsave(here::here("out_figures", "dendrogram_sensitivity_bin_1m_6m_bf.png"), plot = p10,
#        width = fig_size + 5, height = fig_size, dpi = 600, units = "in")
# 
# # Supplementary: Individual PFAS cladograms-------------------------------------
# # One plot per PFAS, dots at species tips + MRCA nodes for higher taxa
# make_indiv_pfas_plot <- function(results_df,
#                                  pfas_name,
#                                  p_base,
#                                  tip_lookup,
#                                  species_meta,
#                                  tree_display,
#                                  node_xy,
#                                  FORMAT) {
#   
#   pos_col <- "#3333FF"  # red = positive
#   neg_col <- "red"  # green = negative
#   
#   # Get significant results for this PFAS at ALL taxa levels
#   dat <- results_df %>%
#     filter(predictor == pfas_name, !is.na(p_val), p_val < 0.05) %>%
#     mutate(
#       betas     = as.numeric(betas),
#       Beta_IQR  = if ("Beta_IQR" %in% names(.)) as.numeric(Beta_IQR) else NA_real_,
#       plot_beta = dplyr::coalesce(Beta_IQR, betas),
#       Direction = case_when(plot_beta > 0 ~ "Positive",
#                             plot_beta < 0 ~ "Negative",
#                             TRUE          ~ NA_character_),
#       sig_type  = ifelse(!is.na(FDR) & FDR < 0.05, "FDR", "nominal")
#     )
#   
#   cat("  ", pfas_name, "— p<0.05 rows:", nrow(dat), "\n")
#   if (nrow(dat) > 0) print(table(dat$taxa_level))
#   
#   empty <- data.frame(node = integer(), x = numeric(), y = numeric(),
#                       Direction = character(), sig_type = character(),
#                       taxa_level = character(), stringsAsFactors = FALSE)
#   
#   if (nrow(dat) == 0) {
#     cat("  No significant results for", pfas_name, "— returning empty plot\n")
#     nodes_df <- empty
#   } else {
#     rows_out <- lapply(seq_len(nrow(dat)), function(i) {
#       row <- dat[i, ]
#       if (row$taxa_level == "Species") {
#         dn  <- tip_lookup$display_name[tip_lookup$taxonomy_id == row$taxonomy_id]
#         if (length(dn) == 0 || is.na(dn[1])) return(NULL)
#         nid <- which(tree_display$tip.label == dn[1])
#       } else {
#         nid <- get_mrca_node(row$taxa_full, row$taxa_level, species_meta, tree_display)
#       }
#       if (length(nid) == 0 || is.na(nid)) return(NULL)
#       xy <- node_xy %>% filter(node == nid)
#       if (nrow(xy) == 0) return(NULL)
#       data.frame(node       = nid,
#                  x          = xy$x[1],
#                  y          = xy$y[1],
#                  Direction  = row$Direction,
#                  sig_type   = row$sig_type,
#                  taxa_level = row$taxa_level,
#                  stringsAsFactors = FALSE)
#     })
#     nodes_df <- bind_rows(rows_out) %>% distinct()
#     cat("  Dots placed:", nrow(nodes_df), "\n")
#   }
#   
#   p <- p_base %<+% nodes_df +
#     geom_point(aes(x = x, y = y, color = Direction),
#                shape = 16, size = 4.5, alpha = 0.95, na.rm = TRUE) +
#     scale_colour_manual(values       = c("Positive" = pos_col,
#                                          "Negative" = neg_col),
#                         na.translate = FALSE,
#                         name         = pfas_name,
#                         guide        = guide_legend(order = 1,
#                                                     override.aes = list(size = 5,
#                                                                         shape = 16))) +
#     labs(title = paste0(pfas_name, " — significant associations (p < 0.05)")) +
#     theme(
#       plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
#       legend.position  = "right",
#       legend.key.size  = unit(0.6, "cm"),
#       legend.text      = element_text(size = 11),
#       legend.title     = element_text(size = 13, face = "bold"),
#       plot.background  = element_rect(fill = "white", color = NA),
#       panel.background = element_rect(fill = "white", color = NA)
#     )
#   
#   return(p)
# }
# 
# # Continuous PFAS supplementary plots (1m → 6m, main prospective scenario)
# supp_cont_plots <- lapply(pfas_order_continuous, function(pfas) {
#   cat("\n--- Continuous:", pfas, "---\n")
#   make_indiv_pfas_plot(
#     results_df   = all_results[["cont_1m6m"]],
#     pfas_name    = pfas,
#     p_base       = p_base,
#     tip_lookup   = tip_lookup,
#     species_meta = species_meta,
#     tree_display = tree_display,
#     node_xy      = node_xy,
#     FORMAT       = FORMAT
#   )
# })
# names(supp_cont_plots) <- pfas_order_continuous
# 
# # Binary PFAS supplementary plots (1m → 6m)
# supp_bin_plots <- lapply(pfas_order_binary, function(pfas) {
#   cat("\n--- Binary:", pfas, "---\n")
#   make_indiv_pfas_plot(
#     results_df   = all_results[["bin_1m6m"]],
#     pfas_name    = pfas,
#     p_base       = p_base,
#     tip_lookup   = tip_lookup,
#     species_meta = species_meta,
#     tree_display = tree_display,
#     node_xy      = node_xy,
#     FORMAT       = FORMAT
#   )
# })
# names(supp_bin_plots) <- pfas_order_binary
# 
# # Preview any individual plot
# suppressWarnings(print(supp_cont_plots[["PFOS"]]))
# suppressWarnings(print(supp_bin_plots[["PFBA"]]))
# 
# # Save all supplementary plots
# # Continuous
# for (pfas in pfas_order_continuous) {
#   ggsave(here::here("out_figures", paste0("supp_dendro_continuous_", pfas, ".pdf")),
#          plot   = supp_cont_plots[[pfas]],
#          width  = fig_size + 3, height = fig_size,
#          device = cairo_pdf, units = "in")
# }
# 
# # Binary
# for (pfas in pfas_order_binary) {
#   ggsave(here::here("out_figures", paste0("supp_dendro_binary_", pfas, ".pdf")),
#          plot   = supp_bin_plots[[pfas]],
#          width  = fig_size + 3, height = fig_size,
#          device = cairo_pdf, units = "in")
# }

#END
