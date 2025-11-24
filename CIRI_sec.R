#secondary analysis
suppressMessages(library(devtools))
suppressMessages(library(monocle3))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(gtools))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(viridis))
suppressMessages(library(ggrepel))
suppressMessages(library(colorspace))
set.seed(1234597698)

dir = "."
ccs = c("5")
group_of_interest = paste("Group", paste(ccs, collapse = "_"), sep = "_")
g = "SOX2"

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
  scale_fill_manual(values = viridis::turbo(length(unique(cluster_condition_percentages_a$gene_a)))) +
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
  scale_fill_manual(values = viridis::turbo(length(unique(cluster_condition_percentages_i$gene_i)))) +
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
    scale_fill_manual(values = viridis::turbo(length(unique(df$clusters)))) +
    labs(x = gene_col, y = "Percentage of cells", fill = "Cluster",
         title = paste("Clusters within", gene_col))
  
  ggsave(p, filename = output_file, width = 20, height = 10)
}

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
  scale_fill_manual(values = lighten(viridis::turbo(2), amount = 0.4)) +
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
  scale_fill_manual(values = lighten(viridis::turbo(2), amount = 0.4)) +
  labs(x = "Gene combination", y = "Percentage of cells", fill = "Cluster group")

ggsave(p, filename = paste0(dir, "/cluster_groups_within_gene_comb_totals_facet_sample_ordered.pdf"), 
       width = 30, height = 7)


# ---- Identify gene_comb with at least 20 cells total ----
gene_to_keep <- totals %>%
  group_by(gene_comb) %>%
  summarise(gene_total = sum(total_cells), .groups = "drop") %>%
  filter(gene_total >= 20) %>%
  pull(gene_comb)

# ---- Filter data for plotting ----
plot_data <- gene_comb_cluster_percentages %>%
  filter(gene_comb %in% gene_to_keep)

totals_filtered <- totals %>%
  filter(gene_comb %in% gene_to_keep)

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
  scale_fill_manual(values = lighten(viridis::turbo(2), amount = 0.4)) +
  labs(x = "Sample", y = "Percentage of cells", fill = "Cluster group")

ggsave(p, filename = paste0(dir, "/cluster_groups_within_gene_comb_totals_filtered_20.pdf"), 
       width = 30, height = 7)

# ---- Identify gene_comb with at least 10 cells per sample ----
gene_to_keep <- totals %>%
  group_by(sample, gene_comb) %>%
  summarise(gene_total = sum(total_cells), .groups = "drop") %>%
  filter(gene_total >= 10)

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
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    scale_fill_manual(values = lighten(viridis::turbo(2), amount = 0.4))
  
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
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
    scale_fill_manual(values = lighten(viridis::turbo(2), amount = 0.4)) 
  
  ggsave(p, filename = paste0(dir, "/cluster_groups_comb_", s, "_ordered.pdf"), 
         width = 12, height = 6)
}


# ---- Step 1: Identify combs passing filter in both samples ----
# Filter totals for cells > 10
totals_filtered_min <- totals %>% filter(total_cells > 10)

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
  pivot_wider(names_from = sample, values_from = percentage)

# Rename columns for clarity (assuming your samples are named "S1" and "S2")
colnames(percentages_wide) <- c("comb", "sample1", "sample2")
percentages_wide <- percentages_wide %>%
  separate(comb, into = c("guide_a", "guide_i"),
           sep = "-", remove = FALSE, fill = "right")

# ---- Step 4: Scatter plot ----
percentages_wide = left_join(percentages_wide, totals_filtered_min)

p <- ggplot() +
  geom_point(data = filter(percentages_wide, sample == 1), mapping = aes(x = sample1, y = sample2, size = total_cells), alpha = 0.4, colour = "#f0f921") +
  geom_point(data = filter(percentages_wide, sample == 2), mapping = aes(x = sample1, y = sample2, size = total_cells), alpha = 0.4, colour = "#0d0887") +
  geom_text_repel(data = filter(percentages_wide, sample == 1), mapping = aes(x = sample1, y = sample2, label = comb), size = 1.5, max.overlaps = Inf, inherit.aes = F) +
  theme_minimal() +
  labs(x = "Percentage in Sample 1", 
       y = "Percentage in Sample 2", 
       title = paste0(group_of_interest, " percentages per combination"))

ggsave(p, filename = paste0(dir, "/scatter_", group_of_interest, "_samples.pdf"),
       width = 6, height = 6)

##HEATMAP
p <- ggplot(percentages_wide, aes(x = guide_a, y = guide_i, fill = sample1)) +
  geom_tile(color = "white") +
  scale_fill_viridis(option="cividis")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Guide A", y = "Guide I", fill = paste0("% cells in ", group_of_interest))

ggsave(p, filename = paste0(dir, "/heatmap_guideA_guideI_", group_of_interest, "_sample1.pdf"),
       width = 10, height = 8)

p <- ggplot(percentages_wide, aes(x = guide_a, y = guide_i, fill = sample2)) +
  geom_tile(color = "white") +
  scale_fill_viridis(option="cividis")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Guide A", y = "Guide I", fill = paste0("% cells in ", group_of_interest))

ggsave(p, filename = paste0(dir, "/heatmap_guideA_guideI_",group_of_interest, "_sample2.pdf"),
       width = 10, height = 8)
