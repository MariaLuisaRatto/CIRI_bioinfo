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
if (length(args) < 2) {
  stop("Usage: Rscript CIRI_unified.R <directory> <variable_guides (1 or 2)>")
}

dir <- args[1]
variable_guides <- as.integer(args[2])

if (!variable_guides %in% c(1, 2)) {
  stop("Error: variable_guides must be 1 or 2.")
}

setwd(dir) # Ensure we are working in the target directory for saving outputs

# --- Configuration (Hardcoded as per original scripts) ---
a_genes = c("MYOD_1", "MYOD_2")
i_genes = c("NANOG","OCT4","SOX2")

# --- Load Annotations ---
ann = read.csv("./guides.csv", header = F)
names(ann) = c("feature", "seq", "type")
ann$seq = NULL

tmp_a = data.frame(a_genes, rep("a", length(a_genes)))
tmp_i = data.frame(i_genes, rep("i", length(i_genes)))
colnames(tmp_a) = c("feature", "type")
names(tmp_i) = c("feature", "type")
ann = rbind(ann, tmp_a, tmp_i)

# --- Load 10x Data ---
CIRI_file <- H5File$new("./filtered_feature_bc_matrix.h5", mode = "r")
feat <- CIRI_file[["matrix/features"]]
names <- feat[["name"]][]   
types <- feat[["feature_type"]][]
crispr_idx <- which(types == "CRISPR Guide Capture")

data    <- CIRI_file[["matrix/data"]][]
indices <- CIRI_file[["matrix/indices"]][]
indptr  <- CIRI_file[["matrix/indptr"]][]
shape   <- CIRI_file[["matrix/shape"]][]
barcodes <- CIRI_file[["matrix/barcodes"]][]

# Reconstruct Sparse Matrix
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

# --- Basic QC Plotting (Common) ---
# (Simplified for brevity, ensuring key outputs are generated)
umis_total <- CIRI_long %>%
  group_by(cell_barcode) %>%
  summarise(total_umis = sum(umi), .groups = "drop")

CIRI_long <- CIRI_long %>%
  left_join(umis_total, by = "cell_barcode") %>%
  mutate(percentage = (umi / total_umis) * 100)

write.csv(CIRI_long, "./CIRI_long.csv")

# --- Threshold Detection (Common Logic) ---
get_threshold <- function(df_subset, type_label) {
  hist_data <- hist(df_subset$umi, breaks = 2000, plot = FALSE)
  counts <- hist_data$counts
  mids <- hist_data$mids
  smooth_counts <- zoo::rollmean(counts, k = 2, fill = NA)
  peaks <- findpeaks(smooth_counts, minpeakheight = 10, minpeakdistance = 10)
  
  threshold <- 22.5 # Default fallback
  
  if (!is.null(peaks)) {
    main_peak_idx <- peaks[which.max(peaks[, 1]), 2]
    right_peaks <- peaks[peaks[,2] > main_peak_idx, , drop = FALSE]
    
    if (nrow(right_peaks) > 0) {
      next_peak_idx <- right_peaks[which.min(right_peaks[,2]), 2]
      valley_range <- (main_peak_idx:next_peak_idx)
      valley_idx <- valley_range[which.min(smooth_counts[valley_range])]
      threshold <- mids[valley_idx]
    }
  }
  return(threshold)
}

# Process CRISPRa
crispra_df <- CIRI_long %>% filter(feature %in% a_genes) %>%
  group_by(cell_barcode) %>% summarise(feature="CRISPRa", umi=sum(umi), percentage=sum(percentage))
thresh_a <- get_threshold(filter(crispra_df, umi > 0), "a")
crispra_df <- mutate(crispra_df, TorF = umi > thresh_a)

# Process CRISPRi
crispri_df <- CIRI_long %>% filter(feature %in% i_genes) %>%
  group_by(cell_barcode) %>% summarise(feature="CRISPRi", umi=sum(umi), percentage=sum(percentage))
thresh_i <- get_threshold(filter(crispri_df, umi > 0), "i")
crispri_df <- mutate(crispri_df, TorF = umi > thresh_i)

# Merge Logic (Common)
common = merge(crispra_df, crispri_df, all = T)
common_summary <- common %>%
  filter(feature %in% c("CRISPRa", "CRISPRi")) %>%
  group_by(cell_barcode, feature) %>%
  summarise(TorF = any(TorF), .groups = "drop") %>%
  pivot_wider(names_from = feature, values_from = TorF, values_fill = FALSE) %>%
  mutate(CIRI = CRISPRa & CRISPRi)

# Prepare Base Data for Assignment
all_wide <- common_summary %>%
  left_join(filter(crispra_df, feature=="CRISPRa"), by="cell_barcode") %>%
  rename(umi_a = umi, perc_a = percentage) %>%
  left_join(filter(crispri_df, feature=="CRISPRi"), by="cell_barcode") %>%
  rename(umi_i = umi, perc_i = percentage) %>%
  filter(CRISPRa | CRISPRi)

OG_CIRI = filter(CIRI_long, cell_barcode %in% all_wide$cell_barcode)
CIRI_sec = filter(OG_CIRI, !feature %in% a_genes & !feature %in% i_genes)

# ==========================================
# BRANCHING LOGIC
# ==========================================

if (variable_guides == 1) {
  print("--- Running Single Guide Assignment Logic ---")
  
  df_ratio <- CIRI_sec %>%
    group_by(cell_barcode, type) %>%
    summarise(sorted_umis = list(sort(umi, decreasing = TRUE)), .groups = "drop") %>%
    mutate(
      top1 = sapply(sorted_umis, function(x) if(length(x) >= 1) x[1] else NA),
      top2 = sapply(sorted_umis, function(x) if(length(x) >= 2) x[2] else NA),
      ratio = ifelse(!is.na(top2) & top2 > 0, top1 / top2, Inf)
    )
  
  assigned <- merge(df_ratio, distinct(common_summary[, c("cell_barcode", "CIRI", "CRISPRa", "CRISPRi")]), by = "cell_barcode")
  
  # Single Guide Filters: Top1 >= 10, Ratio >= 5
  assigned <- assigned %>%
    mutate(pass_filter = top1 >= 10 & ratio >= 5) %>%
    group_by(cell_barcode) %>%
    filter(
      (CIRI[1] == TRUE & all(pass_filter[type %in% c("a","i")])) |
        (CRISPRa[1] == TRUE & CRISPRi[1] == FALSE & any(pass_filter[type == "a"])) |
        (CRISPRi[1] == TRUE & CRISPRa[1] == FALSE & any(pass_filter[type == "i"]))
    ) %>%
    ungroup() %>%
    filter(pass_filter == TRUE)
  
  # Join Feature Names
  assigned <- left_join(assigned, CIRI_sec[, c("cell_barcode","feature","umi", "type")], 
                        by = c("cell_barcode", "type", "top1"="umi"))
  
  assigned_wide <- assigned %>%
    pivot_wider(id_cols = cell_barcode, names_from = type, values_from = c(feature, top1)) %>%
    rename(feature_a = feature_a, feature_i = feature_i)
  
} else {
  print("--- Running Dual (2) Guide Assignment Logic ---")
  
  # 1. Deterministic Ranking
  CIRI_sec <- CIRI_sec %>% mutate(gene = str_extract(feature, "^[^_]+"))
  
  # Compute top1
  top1_df <- CIRI_sec %>%
    group_by(cell_barcode, type) %>%
    arrange(desc(umi), feature) %>% slice(1) %>% ungroup() %>%
    select(cell_barcode, type, top1_umi = umi, top1_feature = feature, top1_gene = gene)
  
  # Compute top2/3 with preference for same gene
  ranked_wide <- CIRI_sec %>%
    left_join(top1_df, by = c("cell_barcode", "type")) %>%
    group_by(cell_barcode, type) %>%
    arrange(desc(umi), desc(gene == top1_gene), feature) %>%
    mutate(rank = row_number()) %>%
    filter(rank <= 3) %>%
    pivot_wider(names_from = rank, values_from = c(umi, feature), names_prefix = "top") %>%
    mutate(ratio = ifelse(!is.na(umi_top2) & umi_top2 > 0, (umi_top1 + umi_top2) / umi_top3, Inf)) %>%
    ungroup()
  
  # 2. Dual Guide Filters: (Top1+Top2) >= 4, Ratio >= 10
  assigned <- merge(ranked_wide, distinct(common_summary[, c("cell_barcode", "CIRI", "CRISPRa", "CRISPRi")]), by = "cell_barcode") %>%
    mutate(pass_filter = (umi_top1 + umi_top2 >= 4 & ratio >= 10)) %>%
    group_by(cell_barcode) %>%
    filter(
      (CIRI[1] == TRUE & all(pass_filter[type %in% c("a","i")])) |
        (CRISPRa[1] == TRUE & CRISPRi[1] == FALSE & any(pass_filter[type == "a"])) |
        (CRISPRi[1] == TRUE & CRISPRa[1] == FALSE & any(pass_filter[type == "i"]))
    ) %>%
    ungroup() %>%
    filter(pass_filter == TRUE)
  
  # 3. Consistency Check (A/B logic)
  assigned <- assigned %>%
    mutate(
      gene1 = str_extract(feature_top1, "^[^_]+"), gene2 = str_extract(feature_top2, "^[^_]+"),
      num1 = str_extract(feature_top1, "(?<=_)\\d+"), num2 = str_extract(feature_top2, "(?<=_)\\d+"),
      letter1 = str_extract(feature_top1, "[A-Za-z]$"), letter2 = str_extract(feature_top2, "[A-Za-z]$")
    ) %>%
    filter(
      (is.na(feature_top1) & !is.na(feature_top2)) |
        (!is.na(feature_top1) & is.na(feature_top2)) |
        (gene1 == gene2 & num1 == num2 & ((letter1 == "A" & letter2 == "B") | (letter1 == "B" & letter2 == "A")))
    )
  
  # 4. Collapse Features
  collapse_feats <- function(f1, f2) {
    feats <- sort(na.omit(c(f1, f2)))
    if(length(feats) == 0) return(NA_character_)
    paste(feats, collapse = ";")
  }
  
  assigned_filtered <- assigned %>% mutate(feature_set = mapply(collapse_feats, feature_top1, feature_top2))
  
  assigned_wide <- assigned_filtered %>%
    group_by(cell_barcode, type) %>%
    summarise(feature = if(length(na.omit(feature_set))>0) na.omit(feature_set)[1] else NA_character_, .groups="drop") %>%
    pivot_wider(id_cols = cell_barcode, names_from = type, values_from = feature) %>%
    rename(feature_a = feature_a, feature_i = feature_i)
}

# --- Final Merge and Save ---
final_output <- left_join(assigned_wide, all_wide, by = "cell_barcode")

# Add "type" classification column
final_output <- final_output %>%
  mutate(
    final_type = case_when(
      CIRI == TRUE ~ "CIRI",
      CIRI == FALSE & CRISPRa == TRUE ~ "CRISPRa",
      CIRI == FALSE & CRISPRi == TRUE ~ "CRISPRi",
      TRUE ~ NA_character_
    )
  )

write.csv(final_output, "./annotation_data.csv", row.names = FALSE)
print("Analysis Complete.")