suppressMessages(library(devtools))
suppressMessages(library(monocle3))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(gtools))
suppressMessages(library(dplyr))
suppressMessages(library(ggrepel))
suppressMessages(library(viridis))
suppressMessages(library(colorspace))
set.seed(1234597698)

dir = "."
ccs = c("3")
group_of_interest = paste("Group", paste(ccs, collapse = "_"), sep = "_")

load(paste0(dir, "/processed_cds.Rdata"))


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


# Function to plot cluster percentages per gene with facets by sample
plot_cluster_percentages <- function(df, gene_col, output_file) {
  p <- ggplot(df, aes_string(x = gene_col, y = "percentage", fill = "clusters")) +
    geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
    # Add total cell counts above the bars
    geom_text(data = df,
              aes_string(x = gene_col, y = "101", label = "total"),
              inherit.aes = FALSE,
              size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_grid(~sample) +
    labs(x = gene_col, y = "Percentage of cells", fill = "Cluster",
         title = paste("Clusters within", gene_col))
  
  ggsave(p, filename = output_file, width = 20, height = 10)
}

df = mutate(df, guide_a = paste0(guide_a1, ";", guide_a2))
df = mutate(df, guide_i = paste0(guide_i1, ";", guide_i2))

# ---- Plot gene_a ----
plot_cluster_percentages(
  df = gene_a_cluster_percentages,
  gene_col = "gene_a",
  output_file = paste0(dir, "/clusters_within_gene_a.pdf")
)

# ---- Plot gene_i ----
plot_cluster_percentages(
  df = gene_i_cluster_percentages,
  gene_col = "gene_i",
  output_file = paste0(dir, "/clusters_within_gene_i.pdf")
)

# ---- GUIDE LEVEL----
# guide_a
guide_a_cluster_counts <- df %>%
  group_by(sample, guide_a, clusters) %>%
  summarise(count = n(), .groups = "drop")

guide_a_cluster_percentages <- guide_a_cluster_counts %>%
  group_by(sample, guide_a) %>%
  mutate(total = sum(count),
         percentage = (count / total) * 100) %>%
  ungroup()

# guide_i
guide_i_cluster_counts <- df %>%
  group_by(sample, guide_i, clusters) %>%
  summarise(count = n(), .groups = "drop")

guide_i_cluster_percentages <- guide_i_cluster_counts %>%
  group_by(sample, guide_i) %>%
  mutate(total = sum(count),
         percentage = (count / total) * 100) %>%
  ungroup()

# ---- Plot ----
# guide_a
plot_cluster_percentages(
  df = guide_a_cluster_percentages,
  gene_col = "guide_a",
  output_file = paste0(dir, "/clusters_within_guide_a.pdf")
)

# guide_i
plot_cluster_percentages(
  df = guide_i_cluster_percentages,
  gene_col = "guide_i",
  output_file = paste0(dir, "/clusters_within_guide_i.pdf")
)





####################
#### only consider muscle vs rest 
colData(cds)$clusters = clusters(cds)
df = as.data.frame(colData(cds))

cellsxcluster = distinct(df[,c("clusters", "nomi")]) %>% 
  group_by(clusters) %>%
  mutate(cellsxcluster = n()) %>%
  distinct(clusters, cellsxcluster)

print(cellsxcluster)

# ---- Collapse clusters into 2 groups ----
df <- df %>%
  mutate(cluster_group = ifelse(clusters %in% ccs, group_of_interest, "Group_other"))

# ---- Counts per sample, gene_comb, and cluster_group ----
gene_comb_cluster_counts <- df %>%
  group_by(sample, gene_comb, cluster_group) %>%
  summarise(count = n(), .groups = "drop")

gene_comb_cluster_percentages <- gene_comb_cluster_counts %>%
  group_by(sample, gene_comb) %>%
  mutate(total = sum(count),
         percentage = (count / total) * 100) %>%
  ungroup()

# ---- Order gene_comb within each sample by Group_4_10 percentage ----
order_df <- gene_comb_cluster_percentages %>%
  filter(cluster_group == group_of_interest) %>%
  select(sample, gene_comb, percentage)

gene_comb_cluster_percentages <- gene_comb_cluster_percentages %>%
  left_join(order_df, by = c("sample", "gene_comb"), suffix = c("", "_grp4")) %>%
  mutate(percentage_grp4 = ifelse(is.na(percentage_grp4), 0, percentage_grp4),
         gene_comb = fct_reorder2(gene_comb, sample, percentage_grp4))   # reorder per sample

# ---- Identify comb with at least 10 cells per sample ----
totals <- gene_comb_cluster_percentages %>%
  group_by(sample, gene_comb) %>%
  summarise(total_cells = sum(count), .groups = "drop")

totals$gene_comb <- factor(totals$gene_comb, levels = levels(gene_comb_cluster_percentages$gene_comb))

# ---- Plot with total cell counts on top ----
p <- ggplot(gene_comb_cluster_percentages, 
            aes(x = factor(sample), y = percentage, fill = cluster_group)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
  geom_text(data = totals, 
            aes(x = factor(sample), y = 100, label = total_cells), 
            inherit.aes = FALSE, size = 3, angle = 90) +
  theme_minimal() +
  facet_grid(~gene_comb, scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(x = "Sample", y = "Percentage of cells", fill = "Cluster group")

ggsave(p, filename = paste0(dir, "/cluster_groups_within_gene_comb_totals_ordered.pdf"), 
       width = 30, height = 7)

p <- ggplot(gene_comb_cluster_percentages, 
            aes(x = gene_comb, y = percentage, fill = cluster_group)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
  geom_text(aes(x = gene_comb, y = 105, label = total), 
            inherit.aes = FALSE, size = 3) +
  theme_minimal() +
  facet_grid(rows = vars(sample), scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_text(angle = 90, vjust = 0.5),
        axis.text.x = element_text(angle = 90)) +
  labs(x = "Gene combination", y = "Percentage of cells", fill = "Cluster group")

ggsave(p, filename = paste0(dir, "/cluster_groups_within_gene_comb_totals_facet_sample_ordered.pdf"), 
       width = 30, height = 7)

# ---- Filter data for plotting ----
gene_to_keep <- totals %>%
  group_by(sample, gene_comb) %>%
  summarise(comb_total = sum(total_cells), .groups = "drop") %>%
  filter(comb_total >= 10)

plot_data <- gene_comb_cluster_percentages %>%
  filter(gene_comb %in% gene_to_keep$gene_comb)

totals_filtered <- totals %>%
  filter(gene_comb %in% gene_to_keep$gene_comb)

# ---- Plot with totals on top ----
p <- ggplot(plot_data, 
            aes(x = factor(sample), y = percentage, fill = cluster_group)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
  # Add total cell counts above the bar
  geom_text(data = totals_filtered, 
            aes(x = factor(sample), y = 100, label = total_cells), 
            inherit.aes = FALSE, size = 3, angle = 90) +
  theme_minimal() +
  facet_grid(~gene_comb, scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(x = "Sample", y = "Percentage of cells", fill = "Cluster group")

ggsave(p, filename = paste0(dir, "/cluster_groups_within_gene_comb_totals_filtered_10.pdf"), 
       width = 30, height = 7)

# ---- Loop over samples ----
for (s in unique(gene_to_keep$sample)) {
  print(s)
  # Keep only data for this sample + valid gene_comb
  gene_keep_sample <- gene_to_keep %>% filter(sample == s) %>% pull(gene_comb)
  
  plot_data <- gene_comb_cluster_percentages %>%
    filter(sample == s, gene_comb %in% gene_keep_sample)
  
  totals_filtered <- totals %>%
    filter(sample == s, gene_comb %in% gene_keep_sample)
  
  # ---- Order gene_comb by percentage in Group_4_10 ----
  # Get Group_4_10 percentage per gene_comb (0 if missing)
  order_df <- plot_data %>%
    filter(cluster_group == group_of_interest) %>%
    select(gene_comb, percentage)
  
  # Make a complete ordering vector
  all_gene_combs <- unique(plot_data$gene_comb)
  order_vec <- all_gene_combs %>%
    setNames(all_gene_combs) %>%
    sapply(function(x) {
      if (x %in% order_df$gene_comb) {
        order_df$percentage[order_df$gene_comb == x]
      } else {
        0  # assign 0 if missing in Group_4_10
      }
    })
  
  # Order descending by percentage in Group_4_10
  levels_ordered <- names(sort(order_vec, decreasing = TRUE))
  
  # Assign factor levels safely
  plot_data$gene_comb <- factor(plot_data$gene_comb, levels = levels_ordered)
  totals_filtered$gene_comb <- factor(totals_filtered$gene_comb, levels = levels_ordered)
  
  # ---- Plot ----
  p <- ggplot(plot_data, 
              aes(x = gene_comb, y = percentage, fill = cluster_group)) +
    geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
    geom_text(data = totals_filtered, 
              aes(x = gene_comb, y = 100, label = total_cells), 
              inherit.aes = FALSE, size = 3, angle = 90) +
    theme_minimal() +
    labs(x = "Gene combination", 
         y = "Percentage of cells", 
         fill = "Cluster group",
         title = paste("Sample:", s)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  ggsave(p, filename = paste0(dir, "/cluster_groups_gene_comb_", s, "_ordered.pdf"), 
         width = 12, height = 4)
}


#### at guide level ###########
# ---- Counts per sample, gene_comb, and cluster_group ----
comb_cluster_counts <- df %>%
  group_by(sample, comb, cluster_group) %>%
  summarise(count = n(), .groups = "drop")

# ---- Percentages of cluster_groups within each sample + gene_comb ----
comb_cluster_percentages <- comb_cluster_counts %>%
  group_by(sample, comb) %>%
  mutate(total = sum(count),
         percentage = (count / total) * 100) %>%
  ungroup()
# ---- Identify comb with at least 10 cells per sample ----
totals <- comb_cluster_percentages %>%
  group_by(sample, comb) %>%
  summarise(total_cells = sum(count), .groups = "drop")

comb_to_keep <- totals %>%
  group_by(sample, comb) %>%
  summarise(comb_total = sum(total_cells), .groups = "drop") %>%
  filter(comb_total >= 10)

# ---- Loop over samples ----
for (s in unique(comb_to_keep$sample)) {
  
  # Keep only data for this sample + valid comb
  comb_keep_sample <- comb_to_keep %>% filter(sample == s) %>% pull(comb)
  
  plot_data <- comb_cluster_percentages %>%
    filter(sample == s, comb %in% comb_keep_sample)
  
  totals_filtered <- totals %>%
    filter(sample == s, comb %in% comb_keep_sample)
  
  # ---- Order comb by percentage in Group_4_10 ----
  order_df <- plot_data %>%
    filter(cluster_group == group_of_interest) %>%
    arrange(desc(percentage))
  
  all_combs <- unique(plot_data$comb)
  missing_combs <- setdiff(all_combs, order_df$comb)
  
  if (length(missing_combs) > 0) {
    order_df <- bind_rows(order_df, 
                          data.frame(comb = missing_combs, percentage = 0))
  }
  
  plot_data$comb <- factor(plot_data$comb, 
                           levels = order_df %>% arrange(desc(percentage)) %>% pull(comb))
  totals_filtered$comb <- factor(totals_filtered$comb, levels = levels(plot_data$comb))
  
  # ---- Plot ----
  p <- ggplot(plot_data, 
              aes(x = comb, y = percentage, fill = cluster_group)) +
    geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
    geom_text(data = totals_filtered, 
              aes(x = comb, y = 100, label = total_cells), 
              inherit.aes = FALSE, size = 3, angle = 90) +
    theme_minimal() +
    labs(x = "Combination", 
         y = "Percentage of cells", 
         fill = "Cluster group",
         title = paste("Sample:", s)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  ggsave(p, filename = paste0(dir, "/cluster_groups_comb_", s, "_ordered.pdf"), 
         width = 12, height = 6)
}



##HEATMAP
samples <- unique(totals$sample)

for (sample_id in samples) {
  message("Processing sample: ", sample_id)
  
  # ---- Step 1: Filter totals for this sample ----
  totals_filtered <- totals %>%
    filter(sample == sample_id, total_cells >= 10)
  
  # Get combinations passing filter in this sample
  comb_pass <- totals_filtered %>% pull(comb)
  
  if (length(comb_pass) == 0) {
    warning("No combinations with >= 10 cells in sample ", sample_id, ". Skipping.")
    next
  }
  
  # ---- Step 2: Compute percentages for this sample ----
  percentages_sample <- comb_cluster_percentages %>%
    filter(sample == sample_id,
           comb %in% comb_pass,
           cluster_group == group_of_interest) %>%
    select(sample, comb, percentage)
  
  # Skip if empty
  if (nrow(percentages_sample) == 0) {
    warning("No data for ", group_of_interest, " in sample ", sample_id)
    next
  }
  
  # ---- Step 3: Prepare for heatmap ----
  percentages_wide <- percentages_sample %>%
    separate(comb, into = c("guide_a", "guide_i"),
             sep = "-", remove = FALSE, fill = "right") %>%
    left_join(totals_filtered, by = "comb")
  
  # ---- Step 4: Heatmap ----
  p <- ggplot(percentages_wide, aes(x = guide_a, y = guide_i, fill = percentage)) +
    geom_tile(color = "white") +
    scale_fill_viridis(option = "cividis") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 8)
    ) +
    labs(
      x = "Guide A",
      y = "Guide I",
      fill = paste0("% cells in ", group_of_interest),
      title = paste0("Sample: ", sample_id)
    )
  
  # ---- Step 5: Save output ----
  out_path <- file.path(dir, paste0("heatmap_guideA_guideI_", group_of_interest, "_", sample_id, ".pdf"))
  ggsave(p, filename = out_path, width = 10, height = 8)
  message("Saved heatmap for ", sample_id, " → ", out_path)
}

###SCATER PLOT
# ---- Step 1: Identify combs passing filter in both samples ----
# Filter totals for cells > 10
totals_filtered_min <- totals %>% filter(total_cells >= 10)

# Count number of samples passing filter per comb
comb_pass <- totals_filtered_min %>%
  group_by(comb) %>%
  summarise(n_samples = n_distinct(sample), .groups = "drop") %>%
  filter(n_samples == 2) %>%  # keep only comb present in both samples
  pull(comb)

# ---- Step 2: Compute percentages for Group_4_10 per comb per sample ----
percentages_group <- comb_cluster_percentages %>%
  filter(comb %in% comb_pass, cluster_group == group_of_interest) %>%
  select(sample, comb, percentage)

# ---- Step 3: Pivot to wide format (sample 1 vs sample 2) ----
percentages_wide <- percentages_group %>%
  pivot_wider(names_from = sample, values_from = percentage, values_fill = 0)

# Rename columns for clarity (assuming your samples are named "S1" and "S2")
colnames(percentages_wide) <- c("comb", "sample1", "sample2")
percentages_wide <- percentages_wide %>%
  separate(comb, into = c("guide_a", "guide_i"),
           sep = "-", remove = FALSE, fill = "right")

percentages_wide = left_join(percentages_wide, totals_filtered_min)

# ---- Step 4: Scatter plot ----
p <- ggplot() +
  geom_point(data = filter(percentages_wide, sample == 1), mapping = aes(x = sample1, y = sample2, size = total_cells), alpha = 0.2, colour = "red") +
  geom_point(data = filter(percentages_wide, sample == 2), mapping = aes(x = sample1, y = sample2, size = total_cells), alpha = 0.2, colour = "green") +
  geom_text_repel(data = filter(percentages_wide, sample == 1), mapping = aes(x = sample1, y = sample2, label = comb), size = 1.5, max.overlaps = Inf, inherit.aes = F) +
  theme_minimal() +
  labs(x = "Percentage in Sample 1", 
       y = "Percentage in Sample 2", 
       title = paste0(group_of_interest, " percentages per combination"))

ggsave(p, filename = paste0(dir, "/scatter_", group_of_interest, "_samples.pdf"),
       width = 6, height = 6)


##Plot subclusters 
## Function to subset clusters and generate UMAP plots
analyze_clusters <- function(cds, clusters_to_keep, group_name, dir, res = 1e-2, g = "SOX2") {
  
  # ---- Subset CDS ----
  cds_sub <- cds[, colData(cds)$clusters %in% clusters_to_keep]
  
  # ---- Save original clusters ----
  colData(cds_sub)$clusters_main <- colData(cds_sub)$clusters
  
  # ---- Re-run preprocessing & subclustering ----
  cds_sub <- preprocess_cds(cds_sub, num_dim = 50)
  cds_sub <- reduce_dimension(cds_sub, reduction_method = "UMAP")
  cds_sub <- cluster_cells(cds_sub, resolution = res, random_seed = 42)
  
  # Save subclusters separately
  colData(cds_sub)$clusters_sub <- clusters(cds_sub)
  
  # ---- Plot by main clusters ----
  p_main <- plot_cells(
    cds_sub,
    color_cells_by = "clusters_main",
    show_trajectory_graph = FALSE,
    label_cell_groups = TRUE,
    label_leaves = FALSE,
    label_branch_points = FALSE
  )
  ggsave(p_main, filename = paste0(dir, "/UMAP_", group_name, "_mainclusters.pdf"),
         width = 5, height = 4)
  
  # ---- Plot by subclusters ----
  p_sub <- plot_cells(
    cds_sub,
    color_cells_by = "clusters_sub",
    show_trajectory_graph = FALSE,
    label_cell_groups = TRUE,
    label_leaves = FALSE,
    label_branch_points = FALSE
  )
  ggsave(p_sub, filename = paste0(dir, "/UMAP_", group_name, "_subclusters.pdf"),
         width = 5, height = 4)
  
  
  # ---- Extract UMAP coordinates ----
  UMAP_sub <- as.data.frame(reducedDims(cds_sub)$UMAP)
  names(UMAP_sub) <- c("x", "y")
  UMAP_sub$nomi <- rownames(UMAP_sub)
  UMAP_sub <- separate(UMAP_sub, nomi,
                       into = c("cellID","sample", "guide_a1", "guide_a2", "guide_i1", "guide_i2", "gene_a", "gene_i", "gene_comb", "type"),
                       sep = "\\.", remove = FALSE, convert = TRUE)
  UMAP_sub$sample <- as.factor(UMAP_sub$sample)
  
  # ---- Base UMAP plots ----
  p <- ggplot(UMAP_sub, aes(x, y, color = sample, alpha = 0.1)) +
    geom_point(size = 0.4) +
    scale_color_manual(values = viridis::turbo(length(unique(UMAP_sub$sample)))) + theme_minimal()
  ggsave(p, filename = paste0(dir, "/UMAP_", group_name, ".pdf"),
         width = 5, height = 4)
  
  # ---- By type ----
  p <- ggplot(UMAP_sub, aes(x, y, color = type, alpha = 0.1)) +
    geom_point(size = 0.4)+
    scale_color_manual(values = viridis::turbo(length(unique(UMAP_sub$type)))) + theme_minimal()
  ggsave(p, filename = paste0(dir, "/UMAP_type_", group_name, ".pdf"),
         width = 5, height = 4)
  
  p <- ggplot(UMAP_sub, aes(x, y, color = type, alpha = 0.1)) +
    geom_point(size = 0.4) + facet_wrap(vars(type))+
    scale_color_manual(values = viridis::turbo(length(unique(UMAP_sub$type)))) + theme_minimal()
  ggsave(p, filename = paste0(dir, "/UMAP_type_facet_", group_name, ".pdf"),
         width = 18, height = 6)
  
  UMAP_sub = mutate(UMAP_sub, guide_a = paste(guide_a1, guide_a1, sep = ";"))
  UMAP_sub = mutate(UMAP_sub, guide_i = paste(guide_i1, guide_i1, sep = ";"))
  
  # ---- By guides ----
  p <- ggplot(UMAP_sub, aes(x, y, color = guide_a, alpha = 0.1)) +
    geom_point(size = 0.4) + facet_wrap(vars(type))+
    scale_color_manual(values = viridis::turbo(length(unique(UMAP_sub$guide_a)))) + theme_minimal()
  ggsave(p, filename = paste0(dir, "/UMAP_guide_a_", group_name, ".pdf"),
         width = 18, height = 6)
  
  p <- ggplot(UMAP_sub, aes(x, y, color = guide_i, alpha = 0.1)) +
    geom_point(size = 0.4) + facet_wrap(vars(type))+
    scale_color_manual(values = viridis::turbo(length(unique(UMAP_sub$guide_i)))) + theme_minimal()
  ggsave(p, filename = paste0(dir, "/UMAP_guide_i_", group_name, ".pdf"),
         width = 9, height = 3)
  
  # ---- By genes ----
  p <- ggplot(UMAP_sub, aes(x, y, color = gene_a, alpha = 0.1)) +
    geom_point(size = 0.4) + facet_wrap(vars(type))+
    scale_color_manual(values = viridis::turbo(length(unique(UMAP_sub$gene_a)))) + theme_minimal()
  ggsave(p, filename = paste0(dir, "/UMAP_gene_a_", group_name, ".pdf"),
         width = 9, height = 3)
  
  p <- ggplot(UMAP_sub, aes(x, y, color = gene_a, alpha = 0.1)) +
    geom_point(size = 0.4) + facet_wrap(vars(type, gene_a))+
    scale_color_manual(values = viridis::turbo(length(unique(UMAP_sub$gene_a)))) + theme_minimal()
  ggsave(p, filename = paste0(dir, "/UMAP_gene_a_facet_", group_name, ".pdf"),
         width = 9, height = 9)
  
  p <- ggplot(UMAP_sub, aes(x, y, color = gene_i, alpha = 0.1)) +
    geom_point(size = 0.4) + facet_wrap(vars(type))+
    scale_color_manual(values = viridis::turbo(length(unique(UMAP_sub$gene_i)))) + theme_minimal()
  ggsave(p, filename = paste0(dir, "/UMAP_gene_i_", group_name, ".pdf"),
         width = 9, height = 3)
  
  p <- ggplot(UMAP_sub, aes(x, y, color = gene_i, alpha = 0.1)) +
    geom_point(size = 0.4) + facet_wrap(vars(type, gene_i))+
    scale_color_manual(values = viridis::turbo(length(unique(UMAP_sub$gene_i)))) + theme_minimal()
  ggsave(p, filename = paste0(dir, "/UMAP_gene_i_facet_", group_name, ".pdf"),
         width = 9, height = 9)
  
  # ---- By gene combinations ----
  p <- ggplot(UMAP_sub, aes(x, y, color = gene_i, alpha = 0.1)) +
    geom_point(size = 0.4) + facet_wrap(vars(gene_a))+
    scale_color_manual(values = viridis::turbo(length(unique(UMAP_sub$gene_i)))) + theme_minimal()
  ggsave(p, filename = paste0(dir, "/UMAP_facet_a_color_i_", group_name, ".pdf"),
         width = 6, height = 5)
  
  p <- ggplot(UMAP_sub, aes(x, y, color = gene_a, alpha = 0.1)) +
    geom_point(size = 0.4) + facet_wrap(vars(gene_i))+
    scale_color_manual(values = viridis::turbo(length(unique(UMAP_sub$gene_a)))) + theme_minimal()
  ggsave(p, filename = paste0(dir, "/UMAP_facet_i_color_a_", group_name, ".pdf"),
         width = 6, height = 5)
  
  # ---- Marker genes for subclusters ----
  message("Finding marker genes for subclusters...")
  marker_test_res <- top_markers(
    cds_sub,
    group_cells_by = "clusters_sub",
    reference_cells = 1000,
    cores = 8
  )
  
  top_specific_markers <- marker_test_res %>%
    filter(fraction_expressing >= 0.10) %>%
    group_by(cell_group) %>%
    top_n(1, pseudo_R2)
  
  top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
  
  # Save data
  write.csv(marker_test_res,
            paste0(dir, "/gene_markers_", group_name, ".csv"),
            row.names = FALSE)
  write.csv(top_specific_markers,
            paste0(dir, "/top1_gene_markers_", group_name, ".csv"),
            row.names = FALSE)
  
  # Plot expression of top markers
  p <- plot_cells(
    cds_sub,
    genes = unique(top_specific_marker_ids),
    label_cell_groups = TRUE,
    show_trajectory_graph = FALSE
  ) + labs(x = "UMAP 1", y = "UMAP 2", title = paste0("Top markers - ", group_name))
  
  ggsave(p, filename = paste0(dir, "/UMAP_top_markers_", group_name, ".pdf"),
         width = 9, height = 9)
  
  #PLOT GENE EXPRESSION
  genes = "genes_of_interest.txt"
  glist = c(unique(t(read.delim(paste0(dir, "/", genes), sep = ",", header = F))))
  #put gene names in $gene_short_name to solve bug
  rowData(cds_sub)$gene_name <- rownames(cds_sub)
  rowData(cds_sub)$gene_short_name <- rowData(cds_sub)$gene_name
  
  p = plot_cells(cds_sub,
                 genes = glist,
                 label_cell_groups = T,
                 show_trajectory_graph = FALSE)+
    labs(x = "UMAP 1", y = "UMAP 2", title = "")#+
  #scale_color_viridis(option="E", discrete=F)
  
  ggsave(p, filename = paste0(dir,"/UMAP_gene_expression_", group_name, ".pdf"),
         width = 10, height = 10)
  
  p = plot_cells(cds_sub,
                 genes = c(unique(UMAP_sub$gene_a), unique(UMAP_sub$gene_i)),
                 label_cell_groups = T,
                 show_trajectory_graph = FALSE)+
    labs(x = "UMAP 1", y = "UMAP 2", title = "")
  
  ggsave(p, filename = paste0(dir,"/UMAP_gene_expression_CRISPR_", group_name, ".pdf"),
         width = 10, height = 10)
  
  #pseudotime
  cds_sub <- learn_graph(cds_sub)
  
  # Helper function to identify the root principal node with the lowest TTN expression
  get_root_lowest_TTN <- function(cds, gene = g) {
    # Extract expression matrix
    exprs_mat <- exprs(cds)
    
    # Check that the gene exists
    if (!(gene %in% rownames(exprs_mat))) {
      stop(paste("Gene", gene, "not found in expression matrix."))
    }
    
    # Get expression of TTN
    ttn_expr <- exprs_mat[gene, ]
    
    # Identify the cell(s) with the lowest TTN expression (CHNAGED TO MAX FOR SOX2)
    min_cells <- names(ttn_expr[ttn_expr == max(ttn_expr, na.rm = TRUE)])
    
    # Find the principal graph projection (which vertex each cell maps to)
    closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    
    # Subset only the vertices for the cells with lowest TTN expression
    vertices <- closest_vertex[min_cells, , drop = FALSE]
    
    # Find the most common vertex among these cells
    root_pr_node <- names(which.max(table(vertices)))
    root_pr_node = paste0("Y_", root_pr_node)
    
    return(root_pr_node)
  }
  
  # Then use it when ordering cells:
  cds_sub <- order_cells(cds_sub, root_pr_nodes = get_root_lowest_TTN(cds_sub))
  #cds_sub <- order_cells(cds_sub)
  pt = as.data.frame(pseudotime(cds_sub))
  names(pt) = c("pseudotime")
  write.csv(pt, paste0(dir, "/pseudotime_", group_name, ".csv"))
  p_sub = plot_cells(cds_sub,
                     color_cells_by = "pseudotime",
                     label_cell_groups=FALSE,
                     label_leaves=FALSE,
                     label_branch_points=FALSE,
                     graph_label_size=1.5)
  ggsave(p_sub, filename = paste0(dir, "/UMAP_", group_name, "_pseudotime.pdf"),
         width = 5, height = 4)
  
  # ---- ENRICHMENT / DEPLETION (subclusters) ----
  print("Running enrichment/depletion (subclusters)...")
  df_sub <- as.data.frame(colData(cds_sub))
  cellsxsubcluster <- df_sub %>%
    distinct(clusters_sub, nomi) %>%
    group_by(clusters_sub) %>%
    mutate(cellsxsubcluster = n()) %>%
    distinct(clusters_sub, cellsxsubcluster)
  
  # Gene A
  gene_a_subcluster_counts <- df_sub %>%
    group_by(sample, gene_a, clusters_sub) %>%
    summarise(count = n(), .groups = "drop")
  gene_a_subcluster_percentages <- gene_a_subcluster_counts %>%
    group_by(sample, gene_a) %>%
    mutate(total = sum(count),
           percentage = (count / total) * 100) %>%
    ungroup()
  
  p7 <- ggplot(gene_a_subcluster_percentages, 
               aes(x = factor(sample), y = percentage, fill = clusters_sub)) +
    geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
    # Add total cell counts above the bar
    geom_text(aes(x = factor(sample), y = 100, label = total), 
              inherit.aes = FALSE, size = 3, angle = 90) +
    theme_minimal() +
    theme(strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), paper = "white") +
    facet_grid(~gene_a) +
    labs(x = "Sample", y = "Percentage of cells", fill = "Subcluster")+
    scale_fill_manual(values = lighten(viridis::turbo(length(unique(gene_a_subcluster_percentages$clusters_sub))), amount = 0.4))
  ggsave(p7, filename = paste0(dir, "/subclusters_within_gene_a_", group_name, ".pdf"), 
         width = 10, height = 5)
  
  df_sub = mutate(df_sub, guide_a = paste(guide_a1, guide_a1, sep = ";"))
  df_sub = mutate(df_sub, guide_i = paste(guide_i1, guide_i1, sep = ";"))
  
  # ---- By guide_a ----
  guide_a_subcluster_counts <- df_sub %>%
    group_by(sample, guide_a, clusters_sub) %>%
    summarise(count = n(), .groups = "drop")
  
  guide_a_subcluster_percentages <- guide_a_subcluster_counts %>%
    group_by(sample, guide_a) %>%
    mutate(total = sum(count),
           percentage = (count / total) * 100) %>%
    ungroup()
  
  p_guide_a <- ggplot(guide_a_subcluster_percentages, 
                      aes(x = factor(sample), y = percentage, fill = clusters_sub)) +
    geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
    # Add total cell counts above the bar
    geom_text(aes(x = factor(sample), y = 100, label = total), 
              inherit.aes = FALSE, size = 3, angle = 90) +
    theme_minimal() +
    theme(strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), paper = "white") +
    facet_grid(~guide_a) +
    labs(x = "Sample", y = "Percentage of cells", fill = "Subcluster")+
    scale_fill_manual(values = lighten(viridis::turbo(length(unique(guide_a_subcluster_percentages$clusters_sub))), amount = 0.4))
  
  ggsave(p_guide_a, filename = paste0(dir, "/subclusters_within_guide_a_", group_name, ".pdf"), 
         width = 10, height = 5)
  
  # Gene I
  gene_i_subcluster_counts <- df_sub %>%
    group_by(sample, gene_i, clusters_sub) %>%
    summarise(count = n(), .groups = "drop")
  gene_i_subcluster_percentages <- gene_i_subcluster_counts %>%
    group_by(sample, gene_i) %>%
    mutate(total = sum(count),
           percentage = (count / total) * 100) %>%
    ungroup()
  
  p8 <- ggplot(gene_i_subcluster_percentages, 
               aes(x = factor(sample), y = percentage, fill = clusters_sub)) +
    geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
    # Add total cell counts above the bar
    geom_text(aes(x = factor(sample), y = 100, label = total), 
              inherit.aes = FALSE, size = 3, angle = 90) +
    theme_minimal() +
    facet_grid(~gene_i) +
    theme(strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), paper = "white") +
    labs(x = "Sample", y = "Percentage of cells", fill = "Subcluster")+
    scale_fill_manual(values = lighten(viridis::turbo(length(unique(gene_i_subcluster_percentages$clusters_sub))), amount = 0.4))
  ggsave(p8, filename = paste0(dir, "/subclusters_within_gene_i_", group_name, ".pdf"), 
         width = 10, height = 5)
  
  # ---- By guide_i ----
  guide_i_subcluster_counts <- df_sub %>%
    group_by(sample, guide_i, clusters_sub) %>%
    summarise(count = n(), .groups = "drop")
  
  guide_i_subcluster_percentages <- guide_i_subcluster_counts %>%
    group_by(sample, guide_i) %>%
    mutate(total = sum(count),
           percentage = (count / total) * 100) %>%
    ungroup()
  
  p_guide_i <- ggplot(guide_i_subcluster_percentages, 
                      aes(x = factor(sample), y = percentage, fill = clusters_sub)) +
    geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
    # Add total cell counts above the bar
    geom_text(aes(x = factor(sample), y = 100, label = total), 
              inherit.aes = FALSE, size = 3, angle = 90) +
    theme_minimal() +
    theme(strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), paper = "white") +
    facet_grid(~guide_i) +
    labs(x = "Sample", y = "Percentage of cells", fill = "Subcluster")+
    scale_fill_manual(values = lighten(viridis::turbo(length(unique(guide_i_subcluster_percentages$clusters_sub))), amount = 0.4))
  
  ggsave(p_guide_i, filename = paste0(dir, "/subclusters_within_guide_i_", group_name, ".pdf"), 
         width = 10, height = 5)
  
  #PLOT GENE EXPRESSION
  #put gene names in $gene_short_name to solve bug
  rowData(cds_sub)$gene_name <- rownames(cds_sub)
  rowData(cds_sub)$gene_short_name <- rowData(cds_sub)$gene_name
  
  p = plot_cells(cds_sub,
                 genes = glist,
                 label_cell_groups = T,
                 show_trajectory_graph = FALSE)+
    labs(x = "UMAP 1", y = "UMAP 2", title = "")#+
  #scale_color_viridis(option="E", discrete=F)
  
  ggsave(p, filename = paste0(dir,"/UMAP_gene_expression_", group_name, ".pdf"),
         width = 10, height = 10)
  
  # ---- Return UMAP data for downstream use ----
  return(UMAP_sub)
}

# Usage
UMAP_muscle <- analyze_clusters(cds, clusters_to_keep = ccs, group_name = group_of_interest, dir = dir, res = 0.001e-2)

