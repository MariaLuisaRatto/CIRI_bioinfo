# load 
suppressMessages(library(devtools))
suppressMessages(library(monocle3))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(gtools))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(viridis))
suppressMessages(library(ggrepel))
suppressMessages(library(ggsignif))
suppressMessages(library(ggtext))
set.seed(1234597698)

dir = "/Users/marialuisaratto/scripts/CIRI12_starting_dataset/"

load(paste0(dir, "processed_cds.RData"))
norm_mat = normalized_counts(cds)

#CREATE ANNOTATION DATAFRAME
data = data.frame()
data = as.data.frame(colnames(norm_mat))
colnames(data) = c("nomi")
rownames(data)=data$nomi
# Check number of pieces in the first column name
data = data.frame(nomi = colnames(norm_mat))
rownames(data) = data$nomi
num_pieces = length(strsplit(colnames(norm_mat)[1], "\\.")[[1]])

if(num_pieces == 8){
  data = data.frame(nomi = colnames(norm_mat))
  rownames(data) = data$nomi
  
  data = separate(data, nomi, into = c("cellID","sample","guide_a","guide_i", "gene_a", "gene_i", "comb", "type"), sep = "\\.", remove = FALSE, convert = TRUE)
} else if(num_pieces == 10){
  data = separate(
    data,
    nomi,
    into = c("cellID", "sample", "guide_a1", "guide_a2", "guide_i1", "guide_i2", "gene_a", "gene_i", "comb", "type"),
    sep = "\\.",
    remove = FALSE,
    convert = TRUE,
    fill = "right"
  )
  
  data = mutate(data, guide_a = paste(guide_a1, guide_a2, sep = ";"))
  data = mutate(data, guide_i = paste(guide_i1, guide_i2, sep = ";"))
}

#ORDER IN SAME WAY MATRIX AND ANNOTATION 
data_complete = norm_mat[, order(colnames(norm_mat))]
data = data[order(row.names(data)), ]
colnames(data_complete) = row.names(data)

## plot gene expression of a and i genes
# create lookup to classify genes
gene_class = tibble(
  gene = c(data$gene_a, data$gene_i),
  class = c(rep("a", length(data$gene_a)),
            rep("i", length(data$gene_i)))
) %>%
  filter(!is.na(gene)) %>%
  distinct()

build_violin_df = function(gene){
  if(!(gene %in% rownames(norm_mat))){
    warning(paste("Gene not found:", gene))
    return(NULL)
  }
  
  idx = grep(paste0("^", gene, "$"), rownames(norm_mat))
  if(length(idx) == 0){
    warning(paste("Gene not found:", gene))
    return(NULL)
  }
  
  expr = norm_mat[idx, ]
  
  tmp = tibble(
    cell = names(expr),
    expr = as.numeric(expr),
    perturbed = ifelse(
      data$gene_a %in% gene | data$gene_i %in% gene,
      "perturbed",
      "other"
    ),
    gene = gene
  )
  
  # add a or i tag
  tmp = left_join(tmp, gene_class, by = "gene")
  
  tmp
}

all_targets = unique(c(data$gene_a, data$gene_i))
all_targets = all_targets[!is.na(all_targets)]

violin_df = map_df(all_targets, build_violin_df)
violin_df = violin_df %>%
  mutate(perturbed = factor(perturbed, levels = c("perturbed", "other")))

# colored facet labels
lab_fun = function(g){
  cls = gene_class$class[match(g, gene_class$gene)]
  col = ifelse(cls == "a", "darkgreen", "darkred")
  paste0("<span style='color:", col, "'>", g, "</span>")
}

# Compute t-test per gene
stat_table = violin_df %>%
  group_by(gene) %>%
  summarise(
    test = list(wilcox.test(expr ~ perturbed)),
    p_value = test[[1]]$p.value,
    statistic = test[[1]]$statistic,
    direction = ifelse(
      median(expr[perturbed=="perturbed"]) > median(expr[perturbed=="other"]),
      "perturbed higher",
      "other higher"
    ),
    .groups = "drop"
  )

# Create a label: always p-value, add direction only if p < 0.05
stat_table = stat_table %>%
  mutate(
    label = ifelse(
      p_value < 0.05,
      paste0("p=", signif(p_value, 2), "\n", direction),
      paste0("p=", signif(p_value, 2))
    )
  )

# Add annotations to the plot
p = ggplot(violin_df, aes(x = perturbed, y = expr, fill = perturbed)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.12, outlier.size = 0.3) +
  facet_wrap(~ gene, scales = "free_y", labeller = labeller(gene = lab_fun)) +
  theme_bw() +
  theme(strip.text = element_markdown()) +
  labs(x = "cell group", y = "normalized expression") +
  geom_text(
    data = stat_table,
    aes(
      x = 1.5, 
      y = max(violin_df$expr), 
      label = label
    ),
    inherit.aes = FALSE,
    size = 2
  )

ggsave(p, filename = paste0(dir, "/violin_CRISPR_targets.pdf"),
       width = 10, height = 8)

#### guide level
### CRIPRa
# Subset cells in the cds object
crispra_cells = colData(cds)$type == "CRISPRa"
cds_crispra = cds[, crispra_cells]

# Now get normalized counts only for CRISPRa cells
norm_mat = normalized_counts(cds_crispra)

# Also subset metadata if you use it later
data_a = filter(data, type == "CRISPRa")

# Build the violin dataframe with guide info
build_violin_df_guide = function(gene){
  if(!(gene %in% rownames(norm_mat))){
    warning(paste("Gene not found:", gene))
    return(NULL)
  }
  
  idx = grep(paste0("^", gene, "$"), rownames(norm_mat))
  if(length(idx) == 0){
    warning(paste("Gene not found:", gene))
    return(NULL)
  }
  
  expr = norm_mat[idx, ]
  
  tmp = tibble(
    cell = names(expr),
    expr = as.numeric(expr),
    perturbed = ifelse(
      data_a$gene_a %in% gene,
      "perturbed",
      "other"
    ),
    gene = gene
  )
  
  # add a or i tag
  tmp = left_join(tmp, gene_class, by = "gene")
  
  tmp
}

violin_df_guides = map_df(unique(data_a$gene_a), build_violin_df_guide)
#violin_df_guides = separate(violin_df_guides, cell, into = c("cellID","sample", "guide_a", "guide_i"), sep = "\\.", remove = F, convert = T)
# Check number of pieces in the first column name
num_pieces = length(strsplit(violin_df_guides$cell[1], "\\.")[[1]])

if(num_pieces == 8){
  violin_df_guides = separate(violin_df_guides, cell, into = c("cellID","sample","guide_a","guide_i", "gene_a", "gene_i", "comb", "type"), sep = "\\.", remove = FALSE, convert = TRUE)
} else if(num_pieces == 10){
  violin_df_guides = separate(
    violin_df_guides, cell,
    into = c("cellID", "sample", "guide_a1", "guide_a2", "guide_i1", "guide_i2", "gene_a", "gene_i", "comb", "type"),
    sep = "\\.",
    remove = FALSE,
    convert = TRUE,
    fill = "right"
  )
  
  violin_df_guides = mutate(violin_df_guides, guide_a = paste(guide_a1, guide_a2, sep = ";"))
  violin_df_guides = mutate(violin_df_guides, guide_i = paste(guide_i1, guide_i2, sep = ";"))
}

violin_df_guides = violin_df_guides[, c("cell", "gene", "class", "guide_a", "expr", "perturbed")]
violin_df_guides = violin_df_guides %>%
  mutate(perturbed = factor(perturbed, levels = c("perturbed", "other")))


violin_df_guides_filtered = violin_df_guides %>%
  group_by(guide_a, perturbed) %>%
  filter(n() >= 20) %>%   # keep only groups with at least 10 observations per perturbed level
  ungroup() %>%
  group_by(guide_a) %>%
  filter(n_distinct(perturbed) == 2) %>%  # ensure both perturbed and control exist
  ungroup()

# ##power calc
# # 1. Compute per-gene minimal n
# genes <- unique(violin_df$gene)
# 
# # Sweep n to find minimal per-group sample size that reaches target power
# find_n_empirical <- function(ctrl_values, fold = 1.3, alpha = 0.05, target_power = 0.8,
#                              n_sim = 2000, n_grid = seq(10, 500, by = 5), seed = 1234,
#                              return_table = FALSE) {
#   # Returns the minimal n and achieved power. If return_table=TRUE returns full table.
#   powers <- sapply(n_grid, function(n) {
#     power_wilcox_empirical(ctrl_values = ctrl_values, fold = fold,
#                            n_per_group = n, n_sim = n_sim, alpha = alpha, seed = seed)
#   })
#   df <- data.frame(n = n_grid, power = powers)
#   idx <- which(powers >= target_power)
#   if(length(idx) == 0) {
#     result <- list(min_n = NA, achieved_power = max(powers), table = df)
#   } else {
#     min_n <- n_grid[min(idx)]
#     result <- list(min_n = min_n, achieved_power = powers[min(idx)], table = df)
#   }
#   if(return_table) return(result)
#   result[c("min_n","achieved_power")]
# }
# 
# # Simulation-based power (empirical ctrl distribution)
# power_wilcox_empirical <- function(ctrl_values, fold = 1.3, n_per_group = 50,
#                                    n_sim = 2000, alpha = 0.05, seed = 1234) {
#   # ctrl_values: numeric vector of expression from "other" cells (empirical control)
#   # fold: multiplicative change applied to control to create "perturbed" (1.3 = +30%)
#   # n_per_group: sample size per group to test
#   # n_sim: number of simulation replicates
#   # alpha: significance threshold for each test
#   set.seed(seed)
#   res <- replicate(n_sim, {
#     x <- sample(ctrl_values, size = n_per_group, replace = TRUE)
#     y <- sample(ctrl_values * fold, size = n_per_group, replace = TRUE)
#     p <- wilcox.test(x, y, exact = FALSE)$p.value
#     as.integer(p < alpha)
#   })
#   mean(res)
# }
# 
# pergene <- lapply(genes, function(g) {
#   ctrl <- violin_df %>% filter(gene == g, perturbed == "other") %>% pull(expr)
#   
#   # skip genes with insufficient control distribution
#   if(length(ctrl) < 20) {
#     return(tibble(gene = g, min_n = NA))
#   }
#   
#   # skip genes with >90% zeros (cannot detect 30% increase reliably)
#   if(mean(ctrl == 0) > 0.9) {
#     return(tibble(gene = g, min_n = NA))
#   }
#   
#   # compute required n for +30% shift
#   res <- find_n_empirical(ctrl_values = ctrl,
#                           fold = 1.3,
#                           alpha = 0.05,
#                           target_power = 0.8,
#                           n_sim = 1500,
#                           n_grid = seq(10, 400, 5),
#                           seed = 42)
#   
#   print(paste("Gene:", g))
#   print(res$power_curve)
#   
#   tibble(gene = g, min_n = res$min_n)
# })
# 
# pergene_df <- bind_rows(pergene)
# valid <- pergene_df %>% filter(!is.na(min_n))
# 
# global_cutoff <- quantile(valid$min_n, 0.8)   # 80% of genes powered
# global_cutoff
# 
# ######

# Now compute t-test per guide_a safely
stat_table_guides = violin_df_guides_filtered %>%
  group_by(guide_a) %>%
  summarise(
    test = list(wilcox.test(expr ~ perturbed)),
    p_value = test[[1]]$p.value,
    statistic = test[[1]]$statistic,
    mean_perturbed = mean(expr[perturbed == "perturbed"]),
    mean_other     = mean(expr[perturbed == "other"]),
    direction = ifelse(
      mean_perturbed > mean_other,
      "perturbed higher",
      "other higher"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("p=", signif(p_value, 2), "\n", direction)
  )

# Compute max expression per guide for annotation
stat_table_guides = stat_table_guides %>%
  left_join(
    violin_df_guides_filtered %>% 
      group_by(guide_a) %>% 
      summarise(max_expr = max(expr), .groups = "drop"),
    by = "guide_a"
  ) %>%
  mutate(y_pos = max_expr - 0.1 * max_expr)

# Compute medians for each guide_a
median_table_guides_a = violin_df_guides_filtered %>%
  group_by(guide_a, perturbed) %>%
  summarise(median_expr = median(expr), .groups = "drop") %>%
  pivot_wider(
    names_from = perturbed,
    values_from = median_expr,
    names_prefix = "median_"
  )

stat_export_guides_a <- stat_table_guides %>%
  select(guide_a, p_value, mean_perturbed, mean_other) %>%
  left_join(median_table_guides_a, by = "guide_a")

write.table(
  stat_export_guides_a,
  file = paste0(dir, "CRISPR_violin_stats_guides_a.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Plot
p_guides = ggplot(violin_df_guides_filtered, aes(x = perturbed, y = expr, fill = perturbed)) +
  geom_violin() +
  geom_boxplot(width = 0.12, outlier.size = 0.3) +
  facet_wrap(~ guide_a, scales = "free_y") +
  theme_bw() +
  theme(strip.text = element_markdown()) +
  labs(x = "cell group", y = "normalized expression") +
  geom_text(
    data = stat_table_guides,
    aes(
      x = 1.5,
      y = y_pos,
      label = label
    ),
    inherit.aes = FALSE,
    size = 2
  )

ggsave(p_guides, filename = paste0(dir, "/violin_CRISPR_targets_guide_a_colored.pdf"),
       width = 10, height = 8)

### CRIPRi
# Subset cells in the cds object
crispri_cells = colData(cds)$type == "CRISPRi"
cds_crispri = cds[, crispri_cells]

# Now get normalized counts only for CRISPRa cells
norm_mat = normalized_counts(cds_crispri)

# Also subset metadata if you use it later
data_i = filter(data, type == "CRISPRi")

# Build the violin dataframe with guide info
build_violin_df_guide = function(gene){
  if(!(gene %in% rownames(norm_mat))){
    warning(paste("Gene not found:", gene))
    return(NULL)
  }
  
  idx = grep(paste0("^", gene, "$"), rownames(norm_mat))
  if(length(idx) == 0){
    warning(paste("Gene not found:", gene))
    return(NULL)
  }
  
  expr = norm_mat[idx, ]
  
  tmp = tibble(
    cell = names(expr),
    expr = as.numeric(expr),
    perturbed = ifelse(
      data_i$gene_i %in% gene,
      "perturbed",
      "other"
    ),
    gene = gene
  )
  
  # add a or i tag
  tmp = left_join(tmp, gene_class, by = "gene")
  
  tmp
}

violin_df_guides = map_df(unique(data_i$gene_i), build_violin_df_guide)
#violin_df_guides = separate(violin_df_guides, cell, into = c("cellID","sample", "guide_a", "guide_i"), sep = "\\.", remove = F, convert = T)
# Check number of pieces in the first column name
num_pieces = length(strsplit(violin_df_guides$cell[1], "\\.")[[1]])

if(num_pieces == 8){
  violin_df_guides = separate(violin_df_guides, cell, into = c("cellID","sample","guide_a","guide_i", "gene_a", "gene_i", "comb", "type"), sep = "\\.", remove = FALSE, convert = TRUE)
} else if(num_pieces == 10){
  violin_df_guides = separate(
    violin_df_guides, cell,
    into = c("cellID", "sample", "guide_a1", "guide_a2", "guide_i1", "guide_i2", "gene_a", "gene_i", "comb", "type"),
    sep = "\\.",
    remove = FALSE,
    convert = TRUE,
    fill = "right"
  )
  
  violin_df_guides = mutate(violin_df_guides, guide_a = paste(guide_a1, guide_a2, sep = ";"))
  violin_df_guides = mutate(violin_df_guides, guide_i = paste(guide_i1, guide_i2, sep = ";"))
}

violin_df_guides = violin_df_guides[, c("cell", "gene", "class", "guide_i", "expr", "perturbed")]
violin_df_guides = violin_df_guides %>%
  mutate(perturbed = factor(perturbed, levels = c("perturbed", "other")))

violin_df_guides_filtered = violin_df_guides %>%
  group_by(guide_i, perturbed) %>%
  filter(n() >= 20) %>%   # keep only groups with at least 10 observations per perturbed level
  ungroup() %>%
  group_by(guide_i) %>%
  filter(n_distinct(perturbed) == 2) %>%  # ensure both perturbed and control exist
  ungroup()

# Now compute t-test per guide_a safely
stat_table_guides = violin_df_guides_filtered %>%
  group_by(guide_i) %>%
  summarise(
    test = list(wilcox.test(expr ~ perturbed)),
    p_value = test[[1]]$p.value,
    statistic = test[[1]]$statistic,
    mean_perturbed = mean(expr[perturbed == "perturbed"]),
    mean_other     = mean(expr[perturbed == "other"]),
    direction = ifelse(
      mean_perturbed > mean_other,
      "perturbed higher",
      "other higher"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("p=", signif(p_value, 2), "\n", direction)
  )

# Compute max expression per guide for annotation
stat_table_guides = stat_table_guides %>%
  left_join(
    violin_df_guides_filtered %>% 
      group_by(guide_i) %>% 
      summarise(max_expr = max(expr), .groups = "drop"),
    by = "guide_i"
  ) %>%
  mutate(y_pos = max_expr - 0.1 * max_expr)

median_table_guides_i = violin_df_guides_filtered %>%
  group_by(guide_i, perturbed) %>%
  summarise(median_expr = median(expr), .groups = "drop") %>%
  pivot_wider(
    names_from = perturbed,
    values_from = median_expr,
    names_prefix = "median_"
  )

stat_export_guides_i <- stat_table_guides %>%
  select(guide_i, p_value, mean_perturbed, mean_other) %>%
  left_join(median_table_guides_i, by = "guide_i")

write.table(
  stat_export_guides_i,
  file = paste0(dir, "CRISPR_violin_stats_guides_i.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Plot
p_guides = ggplot(violin_df_guides_filtered, aes(x = perturbed, y = expr, fill = perturbed)) +
  geom_violin() +
  geom_boxplot(width = 0.12, outlier.size = 0.3) +
  facet_wrap(~ guide_i, scales = "free_y") +
  theme_bw() +
  theme(strip.text = element_markdown()) +
  labs(x = "cell group", y = "normalized expression") +
  geom_text(
    data = stat_table_guides,
    aes(
      x = 1.5,
      y = y_pos,
      label = label
    ),
    inherit.aes = FALSE,
    size = 2
  )

ggsave(p_guides, filename = paste0(dir, "/violin_CRISPR_targets_guide_i_colored.pdf"),
       width = 10, height = 8)

print("Done!")
