#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(quantmod)
library(pracma)
library(hdf5r)
library(Matrix)
library(zoo)
library(scales)
library(Seurat)

# --- Argument Parsing ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript CIRI_unified_complete.R <directory> <variable_guides (1 or 2)> <h5matrix> [thresh_a] [thresh_i]")
}

dir <- args[1]
variable_guides <- as.integer(args[2])
h5matrix <- args[3]

# Optional Thresholds
input_thresh_a <- if (length(args) >= 4) as.numeric(args[4]) else -1
input_thresh_i <- if (length(args) >= 5) as.numeric(args[5]) else -1

if (!variable_guides %in% c(1, 2)) {
  stop("Error: variable_guides must be 1 or 2.")
}

setwd(dir)
print(paste0("Working directory: ", dir))

# --- Configuration ---
#a_genes = c("MYOD_1", "MYOD_2")
#i_genes = c("NANOG","OCT4","SOX2")

# --- Load Annotations ---
ann = read.csv(file.path(dir, "guides.csv"), header = F)
names(ann) = c("feature", "type", "fixed")
ann$seq = NULL

a_genes = ann[ann$type == "a" & ann$fixed == "f", "feature"]
i_genes = ann[ann$type == "i" & ann$fixed == "f", "feature"]

# --- Load 10x Data ---
CIRI_file <- H5File$new(h5matrix, mode = "r")
feat <- CIRI_file[["matrix/features"]]
names <- feat[["name"]][]   
types <- feat[["feature_type"]][]
crispr_idx <- which(types == "CRISPR Guide Capture")

data    <- CIRI_file[["matrix/data"]][]
indices <- CIRI_file[["matrix/indices"]][]
indptr  <- CIRI_file[["matrix/indptr"]][]
shape   <- CIRI_file[["matrix/shape"]][]
barcodes <- CIRI_file[["matrix/barcodes"]][]

mat <- new("dgCMatrix", x = as.numeric(data), i = indices, p = indptr, Dim = shape)
rownames(mat) <- names
colnames(mat) <- barcodes

# Subset CRISPR
crispr_mat <- mat[crispr_idx, , drop = FALSE]
CIRI_df <- as.data.frame(t(as.matrix(crispr_mat)))

CIRI_long <- CIRI_df %>%
  tibble::rownames_to_column("cell_barcode") %>%
  pivot_longer(cols = -cell_barcode, names_to = "feature", values_to = "umi")

CIRI_long = left_join(CIRI_long, ann, by = "feature")

# ==============================================================================
# COMMON PLOTTING (QC)
# ==============================================================================

# Total UMIs per guide
plot_df <- CIRI_long %>%
  group_by(feature, type) %>%
  summarise(total_umi = sum(umi), .groups = "drop") %>%
  arrange(type, feature) %>%
  mutate(feature = factor(feature, levels = unique(feature)))

p = ggplot(plot_df, aes(x = feature, y = total_umi, fill = type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Total UMIs per guide", x = "Guide", y = "Total UMIs") +
  scale_y_continuous(labels = label_number(accuracy = 1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(dir, "total_umixguide.pdf"), p)

p = ggplot(CIRI_long, aes(x = umi)) + geom_histogram(binwidth = 1)
ggsave(file.path(dir, "cellranger_UMIxguide_hist.pdf"), p)

umis_total <- CIRI_long %>%
  group_by(cell_barcode) %>%
  summarise(total_umis = sum(umi), .groups = "drop")

CIRI_long <- CIRI_long %>%
  left_join(umis_total, by = "cell_barcode") %>%
  mutate(percentage = (umi / total_umis) * 100)
write.csv(CIRI_long, file.path(dir, "CIRI_long.csv"))

tmp = distinct(CIRI_long[, c("cell_barcode", "total_umis")])
p = ggplot(tmp, aes(x = total_umis)) + geom_histogram(binwidth = 1)
ggsave(file.path(dir, "cellranger_UMIxcell_hist.pdf"), p)

dir.create(file.path(dir, "single_guide_plots"), showWarnings = FALSE)
for (guide in unique(CIRI_long$feature)) {
  df = filter(CIRI_long, feature == guide)
  p = ggplot(df, aes(x = umi)) + geom_histogram(binwidth = 1)
  ggsave(file.path(dir, "single_guide_plots", paste0("single_guide_hist_", guide, ".pdf")), p)
  p = ggplot(df, aes(x = percentage)) + geom_histogram(binwidth = 1)
  ggsave(file.path(dir, "single_guide_plots", paste0("single_guide_hist_percentage_", guide, ".pdf")), p)
}

p = ggplot(CIRI_long, aes(x = umi)) + geom_histogram(binwidth = 1) + facet_wrap(vars(feature))
ggsave(file.path(dir, "single_guide_hist_facet.pdf"), p, width = 10, height = 10)

p = ggplot(CIRI_long, aes(x = percentage)) + geom_histogram(binwidth = 1) + facet_wrap(vars(feature))
ggsave(file.path(dir, "single_guide_hist_percentage_facet.pdf"), p, width = 10, height = 10)

p = ggplot(CIRI_long, aes(x = percentage, y = umi)) +
  geom_count(alpha = 0.1, colour = "blue", fill = "blue") +
  scale_y_log10() + facet_wrap(vars(feature))
ggsave(file.path(dir, "2D_facet.pdf"), p, width = 10, height = 10)

# --- Fixed vs Variable Logic ---
crispra <- CIRI_long %>% filter(type == "a") %>% mutate(fixed = feature %in% a_genes) %>%
  group_by(cell_barcode) %>% summarise(total_umi = sum(umi), fixed_umi = sum(umi[fixed == TRUE]), percentage = (fixed_umi / total_umi) * 100, feature = "CRISPRa", .groups = "drop") %>% filter(total_umi > 0)

crispri <- CIRI_long %>% filter(type == "i") %>% mutate(fixed = feature %in% i_genes) %>%
  group_by(cell_barcode) %>% summarise(total_umi = sum(umi), fixed_umi = sum(umi[fixed == TRUE]), percentage = (fixed_umi / total_umi) * 100, feature = "CRISPRi", .groups = "drop") %>% filter(total_umi > 0)

CRISPRai <- bind_rows(crispra, crispri)
CRISPRai_plot <- CRISPRai %>% mutate(variable_umi = total_umi - fixed_umi)

p <- ggplot(CRISPRai_plot, aes(x = variable_umi, y = fixed_umi, color = feature)) +
  geom_point(alpha = 0.1, size = 0.5) + theme_minimal() +
  labs(title = "Fixed vs Variable UMI per Cell", x = "Variable UMI", y = "Fixed UMI", color = "Feature") +
  scale_x_continuous(labels = label_number(accuracy = 1)) + scale_y_continuous(labels = label_number(accuracy = 1)) + theme(legend.position = "top")
ggsave(file.path(dir, "CRISPRai_fixed_vs_variable_scatter.pdf"), p, width = 12, height = 10)

p = ggplot(crispra, aes(x = percentage)) + geom_histogram(binwidth = 1)
ggsave(file.path(dir, "CRISPRa_hist_percentage_fixed_vs_variable.pdf"), p)

p = ggplot(crispri, aes(x = percentage)) + geom_histogram(binwidth = 1)
ggsave(file.path(dir, "CRISPRi_hist_percentage_fixed_vs_variable.pdf"), p)

p = ggplot(CRISPRai, aes(x = percentage, y = total_umi)) +
  geom_count(alpha = 0.1, colour = "blue", fill = "blue") + scale_y_log10() + facet_wrap(vars(feature))
ggsave(file.path(dir, "2D_CRISPR_facet_fixed_vs_variable.pdf"), p, width = 20, height = 10)


# --- Detailed Fixed Guide Analysis ---
crispra_df <- CIRI_long %>% filter(feature %in% a_genes) %>%
  group_by(cell_barcode) %>% summarise(feature = "CRISPRa", umi = sum(umi), percentage = sum(percentage), .groups = "drop") %>% filter(umi > 0)

crispri_df <- CIRI_long %>% filter(feature %in% i_genes) %>%
  group_by(cell_barcode) %>% summarise(feature = "CRISPRi", umi = sum(umi), percentage = sum(percentage), .groups = "drop") %>% filter(umi > 0)

CRISPR <- bind_rows(crispra_df, crispri_df)

p = ggplot(crispra_df, aes(x = umi)) + geom_histogram(binwidth = 1)
ggsave(file.path(dir, "CRISPRa_hist.pdf"), p)
p = p + coord_cartesian(xlim = c(0,800) ,ylim = c(0, 80)) 
ggsave(file.path(dir, "CRISPRa_hist_80_800.pdf"), p)
p = ggplot(crispra_df, aes(x = percentage)) + geom_histogram(binwidth = 1)
ggsave(file.path(dir, "CRISPRa_hist_percentage.pdf"), p)

p = ggplot(crispri_df, aes(x = umi)) + geom_histogram(binwidth = 1)
ggsave(file.path(dir, "CRISPRi_hist.pdf"), p)
p = p + coord_cartesian(xlim = c(0,1000) ,ylim = c(0, 15)) 
ggsave(file.path(dir, "CRISPRi_hist_15_1000.pdf"), p)
p = ggplot(crispri_df, aes(x = percentage)) + geom_histogram(binwidth = 1)
ggsave(file.path(dir, "CRISPRi_hist_percentage.pdf"), p)

p = ggplot(CRISPR, aes(x = percentage, y = umi)) +
  geom_count(alpha = 0.1, colour = "blue", fill = "blue") + scale_y_log10() + facet_wrap(vars(feature))
ggsave(file.path(dir, "2D_CRISPR_facet.pdf"), p, width = 20, height = 10)

# ==============================================================================
# THRESHOLD CALCULATION & PLOTTING (DYNAMIC N)
# ==============================================================================

calculate_threshold_and_plot <- function(df, type_label, manual_thresh) {
  
  auto_thresh <- 32
  y_fit <- NULL 
  x_fit <- NULL
  
  try({
    # --- Dynamic N Calculation ---
    max_val <- max(df$umi, na.rm = TRUE)
    target_n <- (max_val * 1.5) * 2 
    calc_n <- 2^ceiling(log2(target_n))
    final_n <- max(1024, min(calc_n, 16384))
    
    cat(paste0("Dynamic Density N: ", final_n, " (Max UMI: ", max_val, ")\n"))
    
    # --- DENSITY ESTIMATION ---
    d <- density(df$umi, adjust = 0.1, n = final_n)
    
    x_fit <- d$x
    y_fit <- d$y
    d1 <- diff(y_fit) 
    main_peak_idx <- which.max(y_fit)
    
    if (main_peak_idx < length(x_fit)) {
      search_indices <- (main_peak_idx):(length(x_fit) - 2)
      valley_found <- FALSE
      for (i in search_indices) {
        if (d1[i] < 0 && d1[i+1] >= 0) {
          auto_thresh <- x_fit[i]
          valley_found <- TRUE
          break 
        }
      }
      
      if (valley_found) {
        cat(paste0("Density Fit (", type_label, "): Valley found at ", round(auto_thresh, 2), "\n"))
      } else {
        cat(paste0("Density Fit (", type_label, "): No valley found. Using default.\n"))
      }
    }
  }, silent = TRUE)
  
  final_thresh <- if(manual_thresh != -1) manual_thresh else auto_thresh
  print(paste("Final threshold for", type_label, ":", final_thresh))
  
  # 4. Standard Plotting
  p = ggplot(df, aes(x = umi)) + 
    geom_histogram(binwidth = 1, fill = "gray70") + 
    geom_vline(xintercept = final_thresh, col = "red", linetype = "dashed") +
    theme_minimal() +
    labs(title = paste("Threshold:", type_label), subtitle = paste("Cutoff:", round(final_thresh, 2)))
  ggsave(file.path(dir, paste0("CRISPR", type_label, "_hist_threshold.pdf")), p)
  
  # 5. Curve Verification Plot
  if (!is.null(y_fit)) {
    scale_factor <- length(df$umi) * 1
    curve_df <- data.frame(x = x_fit, y = y_fit * scale_factor)
    max_x_plot <- max(quantile(df$umi, 0.99), final_thresh * 2)
    
    p_curve = ggplot() +
      geom_histogram(data = df, aes(x = umi), binwidth = 1, fill = "gray80", color = NA) +
      geom_line(data = curve_df, aes(x = x, y = y), color = "blue", alpha = 0.8) +
      geom_vline(xintercept = final_thresh, col = "red", linetype = "dashed") +
      coord_cartesian(xlim = c(0, max_x_plot)) +
      theme_minimal() +
      labs(title = paste("Density Fit Check:", type_label), 
           subtitle = paste0("Blue = Density (Scaled), n=", final_n))
    
    ggsave(file.path(dir, paste0("CRISPR", type_label, "_hist_with_curve.pdf")), p_curve)
  }
  
  return(final_thresh)
}

# Apply to CRISPRa
thresh_a_val <- calculate_threshold_and_plot(crispra_df, "a", input_thresh_a)
crispra_df <- mutate(crispra_df, TorF = umi > thresh_a_val)

p = ggplot(crispra_df, aes(x = percentage, y = umi, colour = TorF)) +
  geom_count(alpha = 0.2) + scale_y_log10() + geom_hline(yintercept=thresh_a_val, colour = "red")
ggsave(file.path(dir, "2D_crispra_threshold.pdf"), p, width = 20, height = 10)

# Apply to CRISPRi
thresh_i_val <- calculate_threshold_and_plot(crispri_df, "i", input_thresh_i)
crispri_df <- mutate(crispri_df, TorF = umi > thresh_i_val)

p = ggplot(crispri_df, aes(x = percentage, y = umi, colour = TorF)) +
  geom_count(alpha = 0.2) + scale_y_log10() + geom_hline(yintercept=thresh_i_val, colour = "red")
ggsave(file.path(dir, "2D_crispri_threshold.pdf"), p, width = 20, height = 10)

# --- Common Summary Merge ---
common = merge(crispra_df, crispri_df, all = T)
common_summary <- common %>%
  filter(feature %in% c("CRISPRa", "CRISPRi")) %>%
  group_by(cell_barcode, feature) %>%
  summarise(TorF = any(TorF), .groups = "drop") %>%
  pivot_wider(names_from = feature, values_from = TorF, values_fill = FALSE) %>%
  mutate(CIRI = CRISPRa & CRISPRi)

common = merge(common, common_summary)
p = ggplot(common, aes(x = percentage, y = umi, colour = CIRI)) +
  geom_count(alpha = 0.2) + scale_y_log10() + facet_wrap(vars(feature))
ggsave(file.path(dir, "2D_CIRI_threshold.pdf"), p, width = 20, height = 10)
write.csv(common, file.path(dir, "all_cells_CIRI.csv"))

# --- Prepare Data for Branching ---
all = common
all_wide <- all %>%
  pivot_wider(id_cols = c(cell_barcode, CRISPRa, CRISPRi, CIRI), names_from = feature, values_from = c(umi, percentage), names_sep = "_") %>%
  rename(umi_a = umi_CRISPRa, umi_i = umi_CRISPRi, percentage_a = percentage_CRISPRa, percentage_i = percentage_CRISPRi) %>%
  filter(CRISPRa | CRISPRi)

OG_CIRI = filter(CIRI_long, cell_barcode %in% all_wide$cell_barcode)
CIRI_sec = filter(OG_CIRI, !feature %in% a_genes & !feature %in% i_genes)

# ==============================================================================
# BRANCHING LOGIC
# ==============================================================================

# This variable will hold the intermediate table before final count aggregation
assigned_wide <- NULL

if (variable_guides == 1) {
  print("--- Mode 1: Single Guide Logic ---")
  
  df_ratio <- CIRI_sec %>%
    group_by(cell_barcode, type) %>%
    summarise(sorted_umis = list(sort(umi, decreasing = TRUE)), .groups = "drop") %>%
    mutate(
      top1 = sapply(sorted_umis, function(x) if(length(x) >= 1) x[1] else NA),
      top2 = sapply(sorted_umis, function(x) if(length(x) >= 2) x[2] else NA),
      ratio = ifelse(!is.na(top2) & top2 > 0, top1 / top2, Inf)
    ) %>%
    select(cell_barcode, type, top1, top2, ratio)
  
  # Plots specific to Mode 1
  p = ggplot(df_ratio, aes(x = ratio, fill = type)) +
    geom_histogram(binwidth = 1, position = "identity") +
    facet_wrap(~type, scales = "free_y") + theme_minimal() + coord_cartesian(xlim = c(0, 500)) +
    labs(title = "Distribution of UMI ratio (top1 / top2)", x = "Ratio (Top1 / Top2)", y = "Count")
  ggsave(file.path(dir, "ratio_secondary_guides.pdf"), p, width = 20, height = 10)
  
  p = ggplot(df_ratio, aes(x = ratio, fill = type)) +
    geom_histogram(binwidth = 0.1, position = "identity") +
    facet_wrap(~type, scales = "free_y") + theme_minimal() + coord_cartesian(xlim = c(0, 20)) +
    labs(title = "Distribution of UMI ratio (top1 / top2)", x = "Ratio (Top1 / Top2)", y = "Count")
  ggsave(file.path(dir, "ratio_secondary_guides_zoom.pdf"), p, width = 20, height = 10)
  
  p = ggplot(df_ratio, aes(x = top1, fill = type)) +
    geom_histogram(binwidth = 1, position = "identity") +
    facet_wrap(~type, scales = "free_y") + theme_minimal() +
    labs(title = "Distribution of UMIs of top guide", x = "UMIs", y = "Count")
  ggsave(file.path(dir, "top_guide_hist.pdf"), p, width = 20, height = 10)
  
  p = ggplot(df_ratio, aes(x = top2, fill = type)) +
    geom_histogram(binwidth = 1, position = "identity") +
    facet_wrap(~type, scales = "free_y") + theme_minimal() +
    labs(title = "Distribution of UMIs of second top guide", x = "UMIs", y = "Count")
  ggsave(file.path(dir, "second_top_guide_hist.pdf"), p, width = 20, height = 10)
  
  tmp_plot <- df_ratio %>% pivot_longer(cols = c(top1, top2), names_to = "rank", values_to = "UMIs")
  p <- ggplot(tmp_plot, aes(x = UMIs, fill = rank)) +
    geom_histogram(binwidth = 1, alpha = 0.3, position = "identity") +
    facet_wrap(~type, scales = "free_y") + theme_minimal() +
    coord_cartesian(xlim = c(0, 2000), ylim = c(0, 1000))+
    labs(title = "Distribution of UMIs of top guides", x = "UMIs", y = "Count")
  ggsave(file.path(dir, "top_guides_overlay_hist.pdf"), p, width = 20, height = 10)
  
  assigned = merge(df_ratio, distinct(common[, c("cell_barcode", "CIRI", "CRISPRa", "CRISPRi")]), by = c("cell_barcode"))
  
  assigned <- assigned %>%
    mutate(pass_filter = top1 >= 10 & ratio >= 5) %>%
    group_by(cell_barcode) %>%
    filter(
      (CIRI[1] == TRUE & all(pass_filter[type %in% c("a","i")])) |
        (CRISPRa[1] == TRUE & CRISPRi[1] == FALSE & any(pass_filter[type == "a"])) |
        (CRISPRi[1] == TRUE & CRISPRa[1] == FALSE & any(pass_filter[type == "i"]))
    ) %>%
    ungroup()
  
  assigned_filtered = filter(assigned, pass_filter == T)
  
  CIRI_sec$top1 = CIRI_sec$umi
  assigned_filtered = left_join(assigned_filtered, CIRI_sec[, c("cell_barcode","feature","top1", "type")], by = c("cell_barcode", "top1", "type"))
  
  counts = as.data.frame(table(distinct(assigned_filtered[c("cell_barcode", "feature")])$feature)) %>% arrange(desc(Freq))
  write.csv(counts, file = file.path(dir, "features_count.csv"))
  
  counts_CIRI = as.data.frame(table(distinct(filter(assigned_filtered, CIRI == T)[c("cell_barcode", "feature")])$feature)) %>% arrange(desc(Freq))
  write.csv(counts_CIRI, file = file.path(dir, "features_count_CIRI.csv"))
  
  # Prepare for Final Join
  assigned_wide <- assigned_filtered %>%
    mutate(feature_col = paste0("feature_", type), umi_col = paste0("top1_", type)) %>%
    pivot_wider(id_cols = cell_barcode, names_from = type, values_from = c(feature, top1)) %>%
    rename(feature_a = feature_a, feature_i = feature_i) 
  
} else {
  print("--- Mode 2: Dual Guide Logic ---")
  
  CIRI_sec <- CIRI_sec %>% mutate(gene = str_extract(feature, "^[^_]+"))
  
  top1_df <- CIRI_sec %>%
    group_by(cell_barcode, type) %>%
    arrange(desc(umi), feature) %>% slice(1) %>% ungroup() %>%
    select(cell_barcode, type, top1_umi = umi, top1_feature = feature, top1_gene = gene)
  
  ranked_wide <- CIRI_sec %>%
    left_join(top1_df, by = c("cell_barcode", "type")) %>%
    group_by(cell_barcode, type) %>%
    arrange(desc(umi), desc(gene == top1_gene), feature) %>%
    mutate(rank = row_number()) %>%
    filter(rank <= 3) %>%
    select(cell_barcode, type, rank, umi, feature) %>%
    pivot_wider(names_from = rank, values_from = c(umi, feature), names_prefix = "top") %>%
    mutate(ratio = ifelse(!is.na(umi_top2) & umi_top2 > 0, (umi_top1 + umi_top2) / umi_top3, Inf)) %>%
    ungroup()
  
  tmp <- ranked_wide %>% 
    select(cell_barcode, type, umi_top1, umi_top2, umi_top3) %>%
    pivot_longer(cols = c(umi_top1, umi_top2, umi_top3), names_to = "rank", values_to = "UMIs")
  
  p <- ggplot(tmp, aes(x = UMIs, fill = rank)) +
    geom_histogram(binwidth = 1, alpha = 0.3, position = "identity") +
    facet_wrap(~type, scales = "free_y") + theme_minimal() +
    coord_cartesian(xlim = c(0, 2000), ylim = c(0, 1000)) +
    labs(title = "Distribution of UMIs of top 3 guides", x = "UMIs", y = "Count")
  ggsave(file.path(dir, "top3_guides_overlay_hist.pdf"), p, width = 20, height = 10)
  
  assigned <- merge(ranked_wide, distinct(common[, c("cell_barcode", "CIRI", "CRISPRa", "CRISPRi")]), by = "cell_barcode") %>%
    mutate(pass_filter = (umi_top1 + umi_top2 >= 4 & ratio >= 10)) %>%
    group_by(cell_barcode) %>%
    filter(
      (CIRI[1] == TRUE & all(pass_filter[type %in% c("a","i")])) |
        (CRISPRa[1] == TRUE & CRISPRi[1] == FALSE & any(pass_filter[type == "a"])) |
        (CRISPRi[1] == TRUE & CRISPRa[1] == FALSE & any(pass_filter[type == "i"]))
    ) %>%
    ungroup() %>%
    filter(pass_filter == TRUE)
  
  assigned <- assigned %>%
    filter(!(type == "a" & CRISPRa == FALSE), !(type == "i" & CRISPRi == FALSE)) %>%
    mutate(
      gene1 = str_extract(feature_top1, "^[^_]+"), gene2 = str_extract(feature_top2, "^[^_]+"),
      num1 = str_extract(feature_top1, "(?<=_)\\d+"), num2 = str_extract(feature_top2, "(?<=_)\\d+"),
      letter1 = str_extract(feature_top1, "[A-Za-z]$"), letter2 = str_extract(feature_top2, "[A-Za-z]$")
    ) %>%
    filter(
      (is.na(feature_top1) & !is.na(feature_top2)) |
        (!is.na(feature_top1) & is.na(feature_top2)) |
        (gene1 == gene2 & num1 == num2 & ((letter1 == "A" & letter2 == "B") | (letter1 == "B" & letter2 == "A")))
    ) %>%
    select(-gene1, -gene2, -num1, -num2, -letter1, -letter2)
  
  collapse_feats <- function(f1, f2) {
    feats <- sort(na.omit(c(f1, f2)))
    if(length(feats) == 0) return(NA_character_)
    paste(feats, collapse = ";")
  }
  
  assigned_filtered <- assigned %>% mutate(feature_set = mapply(collapse_feats, feature_top1, feature_top2))
  
  counts <- assigned_filtered %>% select(cell_barcode, feature_set) %>% filter(!is.na(feature_set)) %>% distinct(cell_barcode, feature_set) %>% count(feature_set, name = "Freq") %>% arrange(desc(Freq))
  write.csv(counts, file = file.path(dir, "features_count.csv"), row.names = FALSE)
  
  counts_CIRI <- assigned_filtered %>% filter(CIRI == TRUE) %>% select(cell_barcode, feature_set) %>% filter(!is.na(feature_set)) %>% distinct(cell_barcode, feature_set) %>% count(feature_set, name = "Freq") %>% arrange(desc(Freq))
  write.csv(counts_CIRI, file = file.path(dir, "features_count_CIRI.csv"), row.names = FALSE)
  
  # Prepare for Final Join
  assigned_summary <- assigned_filtered %>%
    group_by(cell_barcode, type) %>%
    summarise(
      feature = if(length(na.omit(feature_set))>0) na.omit(feature_set)[1] else NA_character_,
      top1 = if(all(is.na(umi_top1))) NA_integer_ else max(umi_top1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(id_cols = cell_barcode, names_from = type, values_from = c(feature, top1)) %>%
    rename(feature_a = feature_a, feature_i = feature_i)
  
  assigned_wide <- assigned_summary
}

# ==============================================================================
# FINAL JOIN AND COUNTS (COMMON TO BOTH MODES)
# ==============================================================================

# Join with all_wide to get CIRI/CRISPRa/CRISPRi status
assigned_wide <- left_join(assigned_wide, all_wide, by = "cell_barcode")

# Apply Final Naming/Type Logic
assigned_wide <- assigned_wide %>%
  mutate(
    feature_a = case_when(CRISPRa == TRUE ~ feature_a, TRUE ~ NA_character_),
    feature_i = case_when(CRISPRi == TRUE ~ feature_i, TRUE ~ NA_character_),
    new_names = paste0(cell_barcode, "-", feature_a, "-", feature_i),
    final_type = case_when(
      CIRI == TRUE ~ "CIRI",
      CIRI == FALSE & CRISPRa == TRUE ~ "CRISPRa",
      CIRI == FALSE & CRISPRi == TRUE ~ "CRISPRi",
      TRUE ~ NA_character_
    )
  )

# Calculate Counts (Now safe because CIRI column exists)
assigned_counts <- assigned_wide %>%
  group_by(feature_a, feature_i) %>% summarise(n_cells = n(), .groups = "drop") %>% arrange(desc(n_cells))
write.csv(assigned_counts, file = file.path(dir, "combinations_counts.csv"), row.names = FALSE)

assigned_counts_CIRI <- assigned_wide %>% filter(CIRI == TRUE) %>%
  group_by(feature_a, feature_i) %>% summarise(n_cells = n(), .groups = "drop") %>% arrange(desc(n_cells))
write.csv(assigned_counts_CIRI, file = file.path(dir, "combinations_counts_CIRI.csv"), row.names = FALSE)

write.csv(assigned_wide, file.path(dir, "annotation_data.csv"), row.names = FALSE)
print("Analysis Complete.")

