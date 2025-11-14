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

a_genes = c("MYOD_1", "MYOD_2")
i_genes = c("NANOG","OCT4","SOX2")

ann = read.csv("./guides.csv", header = F)
names(ann) = c("feature", "seq", "type")
#seq is not correct in this file but it is not actually used here
ann$seq = NULL

tmp_a = data.frame(a_genes, rep("a", length(a_genes)))
tmp_i = data.frame(i_genes, rep("i", length(i_genes)))

colnames(tmp_a) = c("feature", "type")
names(tmp_i) = c("feature", "type")

ann = rbind(ann, tmp_a, tmp_i)

#LOAD
CIRI <- H5File$new("./filtered_feature_bc_matrix.h5", mode = "r")
CIRI[["matrix"]]$ls()

feat <- CIRI[["matrix/features"]]
names <- feat[["name"]][]   
types <- feat[["feature_type"]][]

crispr_idx <- which(types == "CRISPR Guide Capture")
guide_names <- names[crispr_idx]

# 2. Load sparse matrix components
data   <- CIRI[["matrix/data"]][]
indices <- CIRI[["matrix/indices"]][]
indptr <- CIRI[["matrix/indptr"]][]
shape <- CIRI[["matrix/shape"]][]  # [features, barcodes]

data <- as.numeric(data)

# Reconstruct full sparse matrix
mat <- new("dgCMatrix",
           x = data,
           i = indices,
           p = indptr,
           Dim = shape)

# Add row and column names
barcodes <- CIRI[["matrix/barcodes"]][]
rownames(mat) <- names
colnames(mat) <- barcodes

# Subset only CRISPR guides
crispr_mat <- mat[crispr_idx, , drop = FALSE]
df_crispr <- as.data.frame(as.matrix(crispr_mat))
CIRI = as.data.frame(t(df_crispr))

CIRI_long <- CIRI %>%
  tibble::rownames_to_column("cell_barcode") %>%
  pivot_longer(
    cols = -cell_barcode,
    names_to = "feature",
    values_to = "umi"
  )

#create final complete df 
CIRI_long = left_join(CIRI_long, ann)

#plotting 
plot_df <- CIRI_long %>%
  group_by(feature, type) %>%
  summarise(total_umi = sum(umi), .groups = "drop")

plot_df <- plot_df %>%
  arrange(type, feature) %>%
  mutate(feature = factor(feature, levels = unique(feature)))

p = ggplot(plot_df, aes(x = feature, y = total_umi, fill = type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "Total UMIs per guide",
    x = "Guide",
    y = "Total UMIs"
  ) +
  scale_y_continuous(labels = label_number(accuracy = 1)) + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
ggsave("./total_umixguide.pdf", p)

#hist of UMI x guide 
p = ggplot(CIRI_long, aes(x = umi)) +
  geom_histogram(binwidth = 1)
ggsave("./cellranger_UMIxguide_hist.pdf", p)

# Total UMIs per cell 
umis_total <- CIRI_long %>%
  group_by(cell_barcode) %>%
  summarise(total_umis = sum(umi), .groups = "drop")

# Add percentage to CIRI_long
CIRI_long <- CIRI_long %>%
  left_join(umis_total, by = "cell_barcode") %>%
  mutate(percentage = ( umi / total_umis) * 100 ) 

write.csv(CIRI_long, "./CIRI_long.csv")

#hist of UMI x cell
tmp = distinct(CIRI_long[, c("cell_barcode", "total_umis")])
p = ggplot(tmp, aes(x = total_umis)) +
  geom_histogram(binwidth = 1)
ggsave("./cellranger_UMIxcell_hist.pdf", p)

#UMIxguide and percentage for each guide
for (guide in unique(CIRI_long$feature)) {
  df = filter(CIRI_long, feature == guide)
  p = ggplot(df, aes(x = umi)) + 
    geom_histogram(binwidth = 1)
  ggsave(paste0("./single_guide_hist_", guide, ".pdf"), p)
  
  p = ggplot(df, aes(x = percentage))+ 
    geom_histogram(binwidth = 1)
  ggsave(paste0("./single_guide_hist_percentage_", guide, ".pdf"), p)
}

p = ggplot(CIRI_long, aes(x = umi)) + 
  geom_histogram(binwidth = 1) + 
  facet_wrap(vars(feature))
ggsave(paste0("./single_guide_hist_facet.pdf"), p, width = 10, height = 10)

p = ggplot(df, aes(x = percentage))+ 
  geom_histogram(binwidth = 1)
ggsave(paste0("./single_guide_hist_percentage_", guide, ".pdf"), p)

p = ggplot(CIRI_long, aes(x = percentage))+ 
  geom_histogram(binwidth = 1)+ 
  facet_wrap(vars(feature))
ggsave(paste0("./single_guide_hist_percentage_facet.pdf"), p, width = 10, height = 10)


p = ggplot(CIRI_long, aes(x = percentage, y = umi))+
  geom_count(alpha = 0.1, colour = "blue", fill = "blue") +
  scale_y_log10() +
  facet_wrap(vars(feature))
ggsave(paste0("./2D_facet.pdf"), p, width = 10, height = 10)


unique(CIRI_long$feature)

# consider CRISPRa and CRISPRi separately -> fixed guides over total 
crispra <- CIRI_long %>%
  filter(type == "a") %>%
  mutate(fixed = feature %in% a_genes) %>%
  group_by(cell_barcode) %>%
  summarise(
    total_umi = sum(umi),  # All UMIs in the cell
    fixed_umi = sum(umi[fixed == TRUE]),  # Only fixed guides
    percentage = (fixed_umi / total_umi) * 100,  # Fixed UMI % of total UMI
    feature = "CRISPRa",
    .groups = "drop"
  )

crispra = filter(crispra, total_umi > 0)

crispri <- CIRI_long %>%
  filter(type == "i") %>%
  mutate(fixed = feature %in% i_genes) %>%
  group_by(cell_barcode) %>%
  summarise(
    total_umi = sum(umi),  # All UMIs in the cell
    fixed_umi = sum(umi[fixed == TRUE]),  # Only fixed guides
    percentage = (fixed_umi / total_umi) * 100,  # Fixed UMI % of total UMI
    feature = "CRISPRi",
    .groups = "drop"
  )

crispri = filter(crispri, total_umi > 0)

CRISPRai <- bind_rows(crispra, crispri)

# Create copy of CRISPRai and add variable_umi
CRISPRai_plot <- CRISPRai %>%
  mutate(variable_umi = total_umi - fixed_umi)

# Scatter plot: variable_umi vs fixed_umi, colored by feature (CRISPRa / CRISPRi)
p <- ggplot(CRISPRai_plot, aes(x = variable_umi, y = fixed_umi, color = feature)) +
  geom_point(alpha = 0.1, size = 0.5) +
  theme_minimal() +
  labs(
    title = "Fixed vs Variable UMI per Cell",
    x = "Variable UMI",
    y = "Fixed UMI",
    color = "Feature"
  ) +
  scale_x_continuous(labels = label_number(accuracy = 1)) +
  scale_y_continuous(labels = label_number(accuracy = 1)) +
  theme(
    legend.position = "top"
  )

ggsave("./CRISPRai_fixed_vs_variable_scatter.pdf", p, width = 12, height = 10)

p = ggplot(crispra, aes(x = percentage)) + 
  geom_histogram(binwidth = 1)
ggsave(paste0("./CRISPRa_hist_percentage_fixed_vs_variable.pdf"), p)

p = ggplot(crispri, aes(x = percentage)) + 
  geom_histogram(binwidth = 1)
ggsave(paste0("./CRISPRi_hist_percentage_fixed_vs_variable.pdf"), p)


p = ggplot(CRISPRai, aes(x = percentage, y = total_umi))+
  geom_count(alpha = 0.1, colour = "blue", fill = "blue") +
  scale_y_log10() +
  facet_wrap(vars(feature))
ggsave(paste0("./2D_CRISPR_facet_fixed_vs_variable.pdf"), p, width = 20, height = 10)


# consider only CRISPRa and CRISPRi fixed guides 
crispra_df <- CIRI_long %>%
  filter(feature %in% a_genes) %>%
  group_by(cell_barcode) %>%
  summarise(feature = "CRISPRa", umi = sum(umi), percentage = sum(percentage), .groups = "drop")

crispra_df = filter(crispra_df, umi > 0)

crispri_df <- CIRI_long %>%
  filter(feature %in% i_genes) %>%
  group_by(cell_barcode) %>%
  summarise(feature = "CRISPRi", umi = sum(umi), percentage = sum(percentage), .groups = "drop")

crispri_df = filter(crispri_df, umi > 0)


CRISPR <- bind_rows(crispra_df, crispri_df)

p = ggplot(crispra_df, aes(x = umi)) + 
  geom_histogram(binwidth = 1)
ggsave(paste0("./CRISPRa_hist.pdf"), p)

p = p + coord_cartesian(xlim = c(0,800) ,ylim = c(0, 80)) 
ggsave(paste0("./CRISPRa_hist_80_800.pdf"), p)

p = ggplot(crispra_df, aes(x = percentage)) + 
  geom_histogram(binwidth = 1)
ggsave(paste0("./CRISPRa_hist_percentage.pdf"), p)


p = ggplot(crispri_df, aes(x = umi)) + 
  geom_histogram(binwidth = 1)
ggsave(paste0("./CRISPRi_hist.pdf"), p)

p = p + coord_cartesian(xlim = c(0,1000) ,ylim = c(0, 15)) 
ggsave(paste0("./CRISPRi_hist_15_1000.pdf"), p)

p = ggplot(crispri_df, aes(x = percentage)) + 
  geom_histogram(binwidth = 1)
ggsave(paste0("./CRISPRi_hist_percentage.pdf"), p)


p = ggplot(CRISPR, aes(x = percentage, y = umi))+
  geom_count(alpha = 0.1, colour = "blue", fill = "blue") +
  scale_y_log10() +
  facet_wrap(vars(feature))
ggsave(paste0("./2D_CRISPR_facet.pdf"), p, width = 20, height = 10)


# COMMON NO FILTER 
length(unique(CIRI$cell_barcode))
length(unique(CIRI_long$cell_barcode))
length(unique(crispra_df$cell_barcode))
length(unique(crispri_df$cell_barcode))

common_a = filter(crispra_df, cell_barcode %in% crispri_df$cell_barcode)
length(unique(common_a$cell_barcode))
common_i = filter(crispri_df, cell_barcode %in% crispra_df$cell_barcode)
length(unique(common_i$cell_barcode))

p = ggplot(common_a, aes(x = umi)) + 
  geom_histogram(binwidth = 1)
ggsave(paste0("./CRISPRa_hist_common.pdf"), p)

p = p + coord_cartesian(xlim = c(0,800) ,ylim = c(0, 80)) 
ggsave(paste0("./CRISPRa_hist_80_800_common.pdf"), p)

p = ggplot(common_a, aes(x = percentage)) + 
  geom_histogram(binwidth = 1)
ggsave(paste0("./CRISPRa_hist_percentage_common.pdf"), p)


p = ggplot(common_i, aes(x = umi)) + 
  geom_histogram(binwidth = 1)
ggsave(paste0("./CRISPRi_hist_common.pdf"), p)

p = p + coord_cartesian(xlim = c(0,1000) ,ylim = c(0, 15)) 
ggsave(paste0("./CRISPRi_hist_15_1000_common.pdf"), p)

p = ggplot(common_i, aes(x = percentage)) + 
  geom_histogram(binwidth = 1)
ggsave(paste0("./CRISPRi_hist_percentage_common.pdf"), p)


p = ggplot(CRISPR, aes(x = percentage, y = umi))+
  geom_count(alpha = 0.1, colour = "blue", fill = "blue") +
  scale_y_log10() +
  facet_wrap(vars(feature))
ggsave(paste0("./2D_CRISPR_facet.pdf"), p, width = 20, height = 10)

#FIND thresholds a
# histogram
hist_data <- hist(crispra_df$umi, breaks = 2000, plot = FALSE)
counts <- hist_data$counts
mids <- hist_data$mids

# smooth counts
smooth_counts <- zoo::rollmean(counts, k = 2, fill = NA)

# find all peaks
peaks <- findpeaks(smooth_counts, minpeakheight = 10, minpeakdistance = 10)

if (!is.null(peaks)) {
  # main peak = highest peak
  main_peak_idx <- peaks[which.max(peaks[, 1]), 2]
  main_peak_x <- mids[main_peak_idx]
  
  # candidate secondary peaks to the right
  right_peaks <- peaks[peaks[,2] > main_peak_idx, , drop = FALSE]
  
  if (nrow(right_peaks) > 0) {
    # take the closest right-side peak
    next_peak_idx <- right_peaks[which.min(right_peaks[,2]), 2]
    
    # find valley between peaks
    valley_range <- (main_peak_idx:next_peak_idx)
    valley_idx <- valley_range[which.min(smooth_counts[valley_range])]
    main_valley_x <- mids[valley_idx]
    
    cat("Main peak at:", main_peak_x, "\n")
    cat("Next peak at:", mids[next_peak_idx], "\n")
    cat("Valley between peaks (threshold):", main_valley_x, "\n")
  } else {
    cat("No secondary peak found to the right of main peak\n")
  }
} else {
  cat("No peaks found\n")
}
#main_valley_x = 250

p = ggplot(crispra_df, aes(x = umi)) + 
  geom_histogram(binwidth = 1)+
  geom_vline(xintercept = main_valley_x, col = "red")#+ 
  #coord_cartesian(xlim = c(0,1000) ,ylim = c(0, 15)) 
ggsave(paste0("./CRISPRa_hist_threshold.pdf"), p)

p = ggplot(crispra_df, aes(x = umi)) + 
  geom_histogram(bins = 2000)+
  geom_vline(xintercept = main_valley_x, col = "red")#+ 
#coord_cartesian(xlim = c(0,1000) ,ylim = c(0, 15)) 
ggsave(paste0("./CRISPRa_hist_threshold_2000bins.pdf"), p)

crispra_df = mutate(crispra_df, TorF = case_when(umi > main_valley_x ~ TRUE, .default = F))


p = ggplot(crispra_df, aes(x = percentage, y = umi, colour = TorF))+
  geom_count(alpha = 0.2) +
  scale_y_log10()+ 
  geom_hline(yintercept=main_valley_x, colour = "red")
ggsave(paste0("./2D_crispra_threshold.pdf"), p, width = 20, height = 10)

sum(crispra_df$TorF)

##############
#FIND thresholds i
# histogram
hist_data <- hist(crispri_df$umi, breaks = 2000, plot = FALSE)
counts <- hist_data$counts
mids <- hist_data$mids

# smooth counts
smooth_counts <- zoo::rollmean(counts, k = 2, fill = NA)

# find all peaks
peaks <- findpeaks(counts, minpeakheight = 10, minpeakdistance = 10)

if (!is.null(peaks)) {
  # main peak = highest peak
  main_peak_idx <- peaks[which.max(peaks[, 1]), 2]
  main_peak_x <- mids[main_peak_idx]
  
  # candidate secondary peaks to the right
  right_peaks <- peaks[peaks[,2] > main_peak_idx, , drop = FALSE]
  
  if (nrow(right_peaks) > 0) {
    # take the closest right-side peak
    next_peak_idx <- right_peaks[which.min(right_peaks[,2]), 2]
    
    # find valley between peaks
    valley_range <- (main_peak_idx:next_peak_idx)
    valley_idx <- valley_range[which.min(smooth_counts[valley_range])]
    main_valley_x <- mids[valley_idx]
    
    cat("Main peak at:", main_peak_x, "\n")
    cat("Next peak at:", mids[next_peak_idx], "\n")
    cat("Valley between peaks (threshold):", main_valley_x, "\n")
  } else {
    cat("No secondary peak found to the right of main peak\n")
  }
} else {
  cat("No peaks found\n")
}
#main_valley_x = 300

p = ggplot(crispri_df, aes(x = umi)) + 
  geom_histogram(binwidth = 1)+
  geom_vline(xintercept = main_valley_x, col = "red")#+ 
#coord_cartesian(xlim = c(0,1000) ,ylim = c(0, 15)) 
ggsave(paste0("./CRISPRi_hist_threshold.pdf"), p)

p = ggplot(crispri_df, aes(x = umi)) + 
  geom_histogram(bins = 2000)+
  geom_vline(xintercept = main_valley_x, col = "red")#+ 
#coord_cartesian(xlim = c(0,1000) ,ylim = c(0, 15)) 
ggsave(paste0("./CRISPRi_hist_threshold_2000bins.pdf"), p)


crispri_df = mutate(crispri_df, TorF = case_when(umi > main_valley_x ~ TRUE, .default = F))


p = ggplot(crispri_df, aes(x = percentage, y = umi, colour = TorF))+
  geom_count(alpha = 0.2) +
  scale_y_log10()+ 
  geom_hline(yintercept=main_valley_x, colour = "red")
ggsave(paste0("./2D_crispri_threshold.pdf"), p, width = 20, height = 10)

sum(crispri_df$TorF)

a_list = unique(c(filter(crispra_df, TorF == T)$cell_barcode))
i_list = unique(c(filter(crispri_df, TorF == T)$cell_barcode))

common = intersect(a_list, i_list)

common = merge(crispra_df, crispri_df, all = T)

common_summary <- common %>%
  filter(feature %in% c("CRISPRa", "CRISPRi")) %>%
  group_by(cell_barcode, feature) %>%
  summarise(TorF = any(TorF), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = feature,
    values_from = TorF,
    values_fill = FALSE
  ) %>%
  mutate(
    CIRI = case_when(
      CRISPRa & CRISPRi ~ TRUE,
      TRUE ~ FALSE
    )
  )

sum(common_summary$CIRI)
common = merge(common, common_summary)

p = ggplot(common, aes(x = percentage, y = umi, colour = CIRI))+
  geom_count(alpha = 0.2) +
  scale_y_log10()+ 
  facet_wrap(vars(feature))
ggsave(paste0("./2D_CIRI_threshold.pdf"), p, width = 20, height = 10)

write.csv(common, "./all_cells_CIRI.csv")


#### check variable I guide distribution 
#Consider also only a and only i!!!
all = common
all_wide <- all %>%
  pivot_wider(
    id_cols = c(cell_barcode, CRISPRa, CRISPRi, CIRI),
    names_from = feature,
    values_from = c(umi, percentage),
    names_sep = "_"
  ) %>%
  rename(
    umi_a = umi_CRISPRa,
    umi_i = umi_CRISPRi,
    percentage_a = percentage_CRISPRa,
    percentage_i = percentage_CRISPRi
  )

all_wide <- all_wide %>%
  filter(CRISPRa | CRISPRi)

OG_CIRI = filter(CIRI_long, cell_barcode %in% all_wide$cell_barcode)
CIRI_sec = filter(OG_CIRI, ! feature %in% a_genes & ! feature %in% i_genes)

df_ratio <- CIRI_sec %>%
  group_by(cell_barcode, type) %>%
  summarise(
    sorted_umis = list(sort(umi, decreasing = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    top1 = sapply(sorted_umis, function(x) if(length(x) >= 1) x[1] else NA),
    top2 = sapply(sorted_umis, function(x) if(length(x) >= 2) x[2] else NA),
    ratio = ifelse(!is.na(top2) & top2 > 0, top1 / top2, Inf)
  ) %>%
  select(cell_barcode, type, top1, top2, ratio)

p = ggplot(df_ratio, aes(x = ratio, fill = type)) +
  geom_histogram(binwidth = 1, position = "identity") +
  facet_wrap(~type, scales = "free_y") +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 500)) +
  labs(
    title = "Distribution of UMI ratio (top1 / top2)",
    x = "Ratio (Top1 / Top2)",
    y = "Count"
  )
ggsave(paste0("./ratio_secondary_guides.pdf"), p, width = 20, height = 10)

p = ggplot(df_ratio, aes(x = ratio, fill = type)) +
  geom_histogram(binwidth = 0.1, position = "identity") +
  facet_wrap(~type, scales = "free_y") +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 20)) +
  labs(
    title = "Distribution of UMI ratio (top1 / top2)",
    x = "Ratio (Top1 / Top2)",
    y = "Count"
  )
ggsave(paste0("./ratio_secondary_guides_zoom.pdf"), p, width = 20, height = 10)

p = ggplot(df_ratio, aes(x = top1, fill = type)) +
  geom_histogram(binwidth = 1, position = "identity") +
  facet_wrap(~type, scales = "free_y") +
  theme_minimal() +
  #coord_cartesian(xlim = c(0, 500)) +
  labs(
    title = "Distribution of UMIs of top guide",
    x = "UMIs",
    y = "Count"
  )
ggsave(paste0("./top_guide_hist.pdf"), p, width = 20, height = 10)

p = ggplot(df_ratio, aes(x = top2, fill = type)) +
  geom_histogram(binwidth = 1, position = "identity") +
  facet_wrap(~type, scales = "free_y") +
  theme_minimal() +
  #coord_cartesian(xlim = c(0, 500)) +
  labs(
    title = "Distribution of UMIs of second top guide",
    x = "UMIs",
    y = "Count"
  )
ggsave(paste0("./second_top_guide_hist.pdf"), p, width = 20, height = 10)

# Reshape data: put top1 and top2 into one column
tmp <- df_ratio %>%
  pivot_longer(cols = c(top1, top2), 
               names_to = "rank", 
               values_to = "UMIs")

# Now plot with alpha for overlapping histograms
p <- ggplot(tmp, aes(x = UMIs, fill = rank)) +
  geom_histogram(binwidth = 1, alpha = 0.3, position = "identity") +
  facet_wrap(~type, scales = "free_y") +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 2000), ylim = c(0, 1000))+
  labs(
    title = "Distribution of UMIs of top guides",
    x = "UMIs",
    y = "Count"
  )

ggsave("./top_guides_overlay_hist.pdf", p, width = 20, height = 10)

### Filter for ratio >= 5 and UMI >= 10 
assigned = merge(df_ratio, distinct(common[, c("cell_barcode", "CIRI", "CRISPRa", "CRISPRi")]), by = c("cell_barcode"))

assigned <- assigned %>%
  # First, apply per-row filters
  mutate(pass_filter = top1 >= 10 & ratio >= 5) %>%
  
  # Then, group by cell_barcode to apply the conditional logic
  group_by(cell_barcode) %>%
  filter(
    # Case 1: CIRI TRUE -> both type a and i must pass
    (CIRI[1] == TRUE & all(pass_filter[type %in% c("a","i")])) |
      
      # Case 2: CRISPRa TRUE & CRISPRi FALSE -> only type a must pass
      (CRISPRa[1] == TRUE & CRISPRi[1] == FALSE & any(pass_filter[type == "a"])) |
      
      # Case 3: CRISPRi TRUE & CRISPRa FALSE -> only type i must pass
      (CRISPRi[1] == TRUE & CRISPRa[1] == FALSE & any(pass_filter[type == "i"]))
  ) %>%
  ungroup()

assigned_filtered = filter(assigned, pass_filter == T)

#Consider cell_barcodes are repeated in long df!!
sum(assigned_filtered$CIRI)
sum(assigned_filtered$CIRI == FALSE & assigned_filtered$CRISPRa == TRUE)
sum(assigned_filtered$CIRI == FALSE & assigned_filtered$CRISPRi == TRUE)

CIRI_sec$top1 = CIRI_sec$umi
assigned_filtered = left_join(assigned_filtered, CIRI_sec[, c("cell_barcode","feature","top1", "type")], by = c("cell_barcode", "top1", "type"))

counts = as.data.frame(table(distinct(assigned_filtered[c("cell_barcode", "feature")])$feature)) %>%
  arrange(desc(Freq))
write.csv(counts, file = "./features_count.csv")

counts_CIRI = as.data.frame(table(distinct(filter(assigned_filtered, CIRI == T)[c("cell_barcode", "feature")])$feature)) %>%
  arrange(desc(Freq))
write.csv(counts_CIRI, file = "./features_count_CIRI.csv")

assigned_wide <- assigned_filtered %>%
  mutate(feature_col = paste0("feature_", type),
         umi_col = paste0("top1_", type)) %>%
  pivot_wider(
    id_cols = cell_barcode,
    names_from = type,
    values_from = c(feature, top1)
  )

assigned_counts <- assigned_wide %>%
  group_by(feature_a, feature_i) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  arrange(desc(n_cells))
write.csv(assigned_counts, file = "./combinations_counts.csv")

assigned_wide = left_join(assigned_wide, all_wide, by = "cell_barcode")

assigned_wide <- assigned_wide %>%
  mutate(
    feature_a = case_when(CRISPRa == TRUE ~ feature_a, TRUE ~ NA_character_),
    feature_i = case_when(CRISPRi == TRUE ~ feature_i, TRUE ~ NA_character_)
  )

#new names and update matrix 
assigned_wide = mutate(assigned_wide, new_names = paste0(cell_barcode, "-", feature_a, "-", feature_i))

assigned_wide <- assigned_wide %>%
  mutate(
    type = case_when(
      CIRI == TRUE ~ "CIRI",
      CIRI == FALSE & CRISPRa == TRUE ~ "CRISPRa",
      CIRI == FALSE & CRISPRi == TRUE ~ "CRISPRi",
      TRUE ~ NA_character_
    )
  )

assigned_counts_CIRI <- filter(assigned_wide, CIRI == T) %>%
  group_by(feature_a, feature_i) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  arrange(desc(n_cells))
write.csv(assigned_counts_CIRI, file = "./combinations_counts_CIRI.csv")

sum(assigned_wide$CIRI)
sum(assigned_wide$CIRI == FALSE & assigned_wide$CRISPRa == TRUE)
sum(assigned_wide$CIRI == FALSE & assigned_wide$CRISPRi == TRUE)

write.csv(assigned_wide, "annotation_data.csv")

library(rCASC)
h5tocsv(
  group = "docker",
  file = "./filtered_feature_bc_matrix.h5",
  type = "h5",
  version = "5"
)

mitoRiboUmi(
  group = "docker",
  scratch.folder = "/Users/marialuisaratto/scripts/scratch/",
  file = "/Users/marialuisaratto/scripts/CIRI/filtered_feature_bc_matrix.csv",
  separator = ",",
  gtf.name = "Homo_sapiens.GRCh38.101.gtf",
  bio.type = "protein_coding",
  umiXgene = 3
)

scannobyGtf(
  group = "docker",
  file = "/Users/marialuisaratto/scripts/CIRI/RiboMito/filtered_feature_bc_matrix.csv",
  gtf.name= "Homo_sapiens.GRCh38.101.gtf",
  biotype = "protein_coding",
  mt =F,
  ribo.proteins =F,
  umiXgene = 3,
  riboStart.percentage = 0,
  riboEnd.percentage = 100,
  mitoStart.percentage = 1,
  mitoEnd.percentage = 100,
  thresholdGenes = 250
)

scannobyGtf(
  group = "docker",
  file = "/30tb/3tb/data/ratto/AB11_screening_CRISPR/EB_analysis/scratch/filtered_feature_bc_matrix.csv",
  gtf.name= "Homo_sapiens.GRCh38.101.gtf",
  biotype = "protein_coding",
  mt =F,
  ribo.proteins =F,
  umiXgene = 3,
  riboStart.percentage = 0,
  riboEnd.percentage = 100,
  mitoStart.percentage = 0,
  mitoEnd.percentage = 20,
  thresholdGenes = 250
)

exp = read.csv("./filtered_annotated_filtered_feature_bc_matrix.csv", row.names = 1)

exp = exp[, colSums(exp != 0) > 0]

old_names = c(colnames(exp))
old_names = sub("\\.", "-", old_names)
colnames(exp) = old_names
exp = exp[, colnames(exp) %in% assigned_wide$cell_barcode]

# Create a named vector: old -> new
name_map <- setNames(assigned_wide$new_names, assigned_wide$cell_barcode)

# Replace column names in exp
colnames(exp) <- name_map[colnames(exp)]

write.csv(exp, "annotated_matrix.csv")







#secondary analysis load 
suppressMessages(library(devtools))
suppressMessages(library(monocle3))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(gtools))
suppressMessages(library(dplyr))
set.seed(1234597698)

dir = "."
file = "annotated_matrix.csv"
res = 3e-4
#exp$X = rownames(exp)

#LOAD ANNOTATED GENE EXP
#exp = read.csv("/30tb/3tb/data/ratto/testing/annotated_silencing_matrix_complete_all_samples.csv", header = T)
exp = read.csv(paste0(dir,"/", file), header = T)
exp = exp[, colSums(exp != 0) > 0]
exp <- separate(
  exp,
  col="X",
  into=c("a","b"),
  sep = ":",
  remove = TRUE,
  convert = FALSE,
  extra = "warn",
  fill = "warn"
)
exp$a <- NULL
names <- make.unique(exp$b, sep=".")
rownames(exp) <- names
exp$b <- NULL

#CREATE ANNOTATION DATAFRAME
data = data.frame()
data = as.data.frame(colnames(exp))
colnames(data) = c("nomi")
rownames(data)=data$nomi
data = separate(data, nomi, into = c("cellID","sample", "guide_a", "guide_i"), sep = "\\.", remove = F, convert = T)
data = mutate(data, comb = paste(guide_a, guide_i, sep = "_"))
data <- data %>%
  mutate(gene_a = sapply(strsplit(guide_a, "_"), `[`, 1))
data <- data %>%
  mutate(gene_i = sapply(strsplit(guide_i, "_"), `[`, 1))
data = mutate(data, gene_comb = paste(gene_a, gene_i, sep = "_"))

data <- data %>%
  mutate(
    type = case_when(
      !is.na(guide_a) & !is.na(guide_i) ~ "CIRI",
      !is.na(guide_a) &  is.na(guide_i) ~ "CRISPRa",
      is.na(guide_a) & !is.na(guide_i) ~ "CRISPRi",
      TRUE ~ NA_character_
    )
  )

#update names
data = mutate(data, nomi = paste(nomi, gene_a, gene_i, gene_comb, type, sep = "."))
rownames(data) = data$nomi

#ORDER IN SAME WAY MATRIX AND ANNOTATION 
data_complete = exp[, order(colnames(exp))]
data = data[order(row.names(data)), ]
colnames(data_complete) = row.names(data)

#save data
write.csv(data_complete, paste(dir, "/espression_data.csv",sep=""))
write.csv(data, paste(dir, "/cell_metadata.csv",sep=""))

#CREATE MONOCLE DATASET
data_matrix = as.matrix(data_complete)
starting_cds <- new_cell_data_set(expression_data = data_matrix,
                                  cell_metadata = data,
                                  gene_metadata = NULL)
# Save dataframe as an RData file
save(starting_cds, file = paste0(dir, "/starting_cds.RData"))

#NORMALIZE for size and log
norm = normalized_counts(
  starting_cds#,norm_method = "log",pseudocount = 1
)

cds <- new_cell_data_set(expression_data = norm,
                         cell_metadata = data,
                         gene_metadata = NULL)

cds <- preprocess_cds(cds, num_dim = 20)
p = plot_pc_variance_explained(cds)
ggsave(p, filename = paste0(dir,"/PCA.pdf"),
       width = 5, height = 5)
#cds <- align_cds(cds, alignment_group = "replicate", useNames = TRUE)
cds <- reduce_dimension(cds, reduction_method="UMAP", umap.fast_sgd = FALSE,cores=1,n_sgd_threads=1)

UMAP = as.data.frame(reducedDims(cds)$UMAP)
names(UMAP) = c("x", "y")
UMAP$nomi = rownames(UMAP)
UMAP = separate(UMAP, nomi,into = c("cellID","sample", "guide_a", "guide_i", "gene_a", "gene_i", "gene_comb", "type"), sep = "\\.", remove = F, convert = T)
UMAP$sample = as.factor(UMAP$sample)
p = ggplot(UMAP, aes(x, y, color = sample, alpha = 0.1)) +
  geom_point(size = 0.4)
ggsave(p, filename = paste0(dir,"/UMAP.pdf"),
       width = 6, height = 6)

#Plot by type
p = ggplot(UMAP, aes(x, y, color = type, alpha = 0.1)) +
  geom_point(size = 0.4)
ggsave(p, filename = paste0(dir,"/UMAP_type.pdf"),
       width = 6, height = 6)

p = ggplot(UMAP, aes(x, y, color = type, alpha = 0.1)) +
  geom_point(size = 0.4) + 
  facet_wrap(vars(type))
ggsave(p, filename = paste0(dir,"/UMAP_type_facet.pdf"),
       width = 18, height = 6)

#Plot by guide 
p = ggplot(UMAP, aes(x, y, color = guide_a, alpha = 0.1)) +
  geom_point(size = 0.4)+ 
  facet_wrap(vars(type))
ggsave(p, filename = paste0(dir,"/UMAP_guide_a.pdf"),
       width = 18, height = 6)

p = ggplot(UMAP, aes(x, y, color = guide_i, alpha = 0.1)) +
  geom_point(size = 0.4)+ 
  facet_wrap(vars(type))
ggsave(p, filename = paste0(dir,"/UMAP_guide_i.pdf"),
       width = 18, height = 6)

#Plot by gene
p = ggplot(UMAP, aes(x, y, color = gene_a, alpha = 0.1)) +
  geom_point(size = 0.4)+ 
  facet_wrap(vars(type))
ggsave(p, filename = paste0(dir,"/UMAP_gene_a.pdf"),
       width = 18, height = 6)

p = ggplot(UMAP, aes(x, y, color = gene_i, alpha = 0.1)) +
  geom_point(size = 0.4)+ 
  facet_wrap(vars(type))
ggsave(p, filename = paste0(dir,"/UMAP_gene_i.pdf"),
       width = 18, height = 6)

#Plot CIRI by gene comb
tmp = filter(UMAP, type == "CIRI")
p = ggplot(tmp, aes(x, y, color = gene_i, alpha = 0.1)) +
  geom_point(size = 0.4)+ 
  facet_wrap(vars(gene_a))
ggsave(p, filename = paste0(dir,"/UMAP_facet_a_color_i.pdf"),
       width = 12, height = 9)

p = ggplot(tmp, aes(x, y, color = gene_a, alpha = 0.1)) +
  geom_point(size = 0.4)+ 
  facet_wrap(vars(gene_i))
ggsave(p, filename = paste0(dir,"/UMAP_facet_i_color_a.pdf"),
       width = 12, height = 9)

#PLOT GENE EXPRESSION
genes = "genes_of_interest.txt"
glist = c(unique(t(read.delim(paste0(dir, "/", genes), sep = ",", header = F))))
#put gene names in $gene_short_name to solve bug
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name

p = plot_cells(cds,
               genes = glist,
               label_cell_groups = T,
               show_trajectory_graph = FALSE)+
  labs(x = "UMAP 1", y = "UMAP 2", title = "")#+
#scale_color_viridis(option="E", discrete=F)

ggsave(p, filename = paste0(dir,"/UMAP_gene_expression.pdf"),
       width = 10, height = 10)

p = plot_cells(cds,
               genes = c(unique(data$gene_a), unique(data$gene_i)),
               label_cell_groups = T,
               show_trajectory_graph = FALSE)+
  labs(x = "UMAP 1", y = "UMAP 2", title = "")#+
#scale_color_viridis(option="E", discrete=F)

ggsave(p, filename = paste0(dir,"/UMAP_gene_expression_CRISPR.pdf"),
       width = 10, height = 10)

#CLUSTERING
cds <- cluster_cells(cds, resolution=res, random_seed = 42)
p = plot_cells(cds)
ggsave(p, filename = paste0(dir,"/UMAP_clustering.pdf"),
       width = 5, height = 5)
colData(cds)$clusters = clusters(cds)


#trajectories
print("Learning trajectories...")
cds <- learn_graph(cds)
p = plot_cells(cds,
               color_cells_by = "cluster",
               label_groups_by_cluster=FALSE,
               label_leaves=FALSE,
               label_branch_points=FALSE)
ggsave(p, filename = paste0(dir,"/UMAP_clustering_trajectories.pdf"),
       width = 5, height = 5)

# Save dataas an RData file
save(cds, file = paste0(dir, "/processed_cds.RData"))

#MARKER GENES
print("Finding clusters marker genes...")
marker_test_res <- top_markers(cds, group_cells_by="cluster",
                               reference_cells=1000, cores=8)

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
#save data
write.csv(marker_test_res, paste(dir, "/gene_markers.csv",sep=""))
write.csv(top_specific_markers, paste(dir, "/top1_gene_markers.csv",sep=""))

p = plot_cells(cds,
               genes = unique(top_specific_markers$gene_id),
               label_cell_groups = T,
               show_trajectory_graph = FALSE)+
  labs(x = "UMAP 1", y = "UMAP 2", title = "")#+
#scale_color_viridis(option="E", discrete=F)

ggsave(p, filename = paste0(dir,"/UMAP_gene_expression_top_markers.pdf"),
       width = 10, height = 10)

print("Done!")



#### enrichment / depletion 
min_cells_cluster = 40

#calculate cells x cluster and filter if too few
colData(cds)$clusters = clusters(cds)
df = as.data.frame(colData(cds))
cellsxcluster = distinct(df[,c("clusters", "nomi")]) %>% 
  group_by(clusters) %>%
  mutate(cellsxcluster = n()) 

cellsxcluster = distinct(cellsxcluster[, c("clusters", "cellsxcluster")])
print(cellsxcluster)
print(paste0("Min cells per clusters: ", min_cells_cluster))
print(paste0("Removing clusters: ", filter(cellsxcluster, cellsxcluster <= min_cells_cluster)[1]))
cellsxcluster = filter(cellsxcluster, cellsxcluster > min_cells_cluster)
df = filter(df, clusters %in% cellsxcluster$clusters)


# Calculate the number of cells per cluster and condition a 
cluster_condition_counts_a <- df %>%
  group_by(sample, clusters, gene_a) %>%
  summarise(count = n(), .groups = "drop")

cluster_condition_percentages_a <- cluster_condition_counts_a %>%
  group_by(sample, clusters) %>%
  mutate(total = sum(count),
         percentage = (count / total) * 100) %>%
  ungroup()

# Plot
p <- ggplot(cluster_condition_percentages_a, aes(x = factor(sample), y = percentage, fill = gene_a)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
  theme_minimal() +  
  facet_grid(~clusters)

ggsave(p, filename = paste0(dir, "/gene_a_in_clusters.pdf"), width = 20, height = 10)


# Calculate the number of cells per cluster and condition i
cluster_condition_counts_i <- df %>%
  group_by(sample, clusters, gene_i) %>%
  summarise(count = n(), .groups = "drop")

cluster_condition_percentages_i <- cluster_condition_counts_i %>%
  group_by(sample, clusters) %>%
  mutate(total = sum(count),
         percentage = (count / total) * 100) %>%
  ungroup()

# Plot
p <- ggplot(cluster_condition_percentages_i, aes(x = factor(sample), y = percentage, fill = gene_i)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
  theme_minimal() +  
  facet_grid(~clusters)

ggsave(p, filename = paste0(dir, "/gene_i_in_clusters.pdf"), width = 20, height = 10)

# Other way around 
# Counts per sample, gene_a, and cluster
gene_a_cluster_counts <- df %>%
  group_by(sample, gene_a, clusters) %>%
  summarise(count = n(), .groups = "drop")

# Percentages of clusters within each sample + gene_a
gene_a_cluster_percentages <- gene_a_cluster_counts %>%
  group_by(sample, gene_a) %>%
  mutate(total = sum(count),
         percentage = (count / total) * 100) %>%
  ungroup()

p <- ggplot(gene_a_cluster_percentages, 
            aes(x = factor(sample), y = percentage, fill = clusters)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
  theme_minimal() +
  facet_grid(~gene_a) +
  labs(x = "Sample", y = "Percentage of cells", fill = "Cluster")

ggsave(p, filename = paste0(dir, "/clusters_within_gene_a.pdf"), 
       width = 20, height = 10)

# Counts per sample, gene_i, and cluster
gene_i_cluster_counts <- df %>%
  group_by(sample, gene_i, clusters) %>%
  summarise(count = n(), .groups = "drop")

# Percentages of clusters within each sample + gene_i
gene_i_cluster_percentages <- gene_i_cluster_counts %>%
  group_by(sample, gene_i) %>%
  mutate(total = sum(count),
         percentage = (count / total) * 100) %>%
  ungroup()

p <- ggplot(gene_i_cluster_percentages, 
            aes(x = factor(sample), y = percentage, fill = clusters)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
  theme_minimal() +
  facet_grid(~gene_i) +
  labs(x = "Sample", y = "Percentage of cells", fill = "Cluster")

ggsave(p, filename = paste0(dir, "/clusters_within_gene_i.pdf"), 
       width = 20, height = 10)

#REGRESSION
print("regression...")
# Fit models for perturbations GENE
gene_fits <- fit_models(cds, model_formula_str = paste0("~gene_a"))
print("genes fitted for genes a...")
# Extract coefficient tables and adjust with BH
fit_coefs <- coefficient_table(gene_fits)

regression = fit_coefs %>% filter(term != "(Intercept)")
regression = regression %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)
write.csv(regression, file= paste0(dir,"/regression_filtered_genes_per_gene_a.csv"), row.names=FALSE)

gene_fits <- fit_models(cds, model_formula_str = paste0("~gene_i"))
print("genes fitted for genes i...")
# Extract coefficient tables and adjust with BH
fit_coefs <- coefficient_table(gene_fits)

regression = fit_coefs %>% filter(term != "(Intercept)")
regression = regression %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)
write.csv(regression, file= paste0(dir,"/regression_filtered_genes_per_gene_i.csv"), row.names=FALSE)

gene_fits <- fit_models(cds, model_formula_str = paste0("~gene_a + gene_i"))
print("genes fitted for genes a + i...")
# Extract coefficient tables and adjust with BH
fit_coefs <- coefficient_table(gene_fits)

regression = fit_coefs %>% filter(term != "(Intercept)")
regression = regression %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)
write.csv(regression, file= paste0(dir,"/regression_filtered_genes_per_gene_a_and_i.csv"), row.names=FALSE)

