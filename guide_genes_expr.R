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
norm_mat <- normalized_counts(cds)

#CREATE ANNOTATION DATAFRAME
data = data.frame()
data = as.data.frame(colnames(norm_mat))
colnames(data) = c("nomi")
rownames(data)=data$nomi
# Check number of pieces in the first column name
data <- data.frame(nomi = colnames(norm_mat))
rownames(data) <- data$nomi
num_pieces <- length(strsplit(colnames(norm_mat)[1], "\\.")[[1]])

if(num_pieces == 8){
  data <- data.frame(nomi = colnames(norm_mat))
  rownames(data) <- data$nomi
  
  data <- separate(data, nomi, into = c("cellID","sample","guide_a","guide_i", "gene_a", "gene_i", "comb", "type"), sep = "\\.", remove = FALSE, convert = TRUE)
} else if(num_pieces == 10){
  data <- separate(
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
gene_class <- tibble(
  gene = c(data$gene_a, data$gene_i),
  class = c(rep("a", length(data$gene_a)),
            rep("i", length(data$gene_i)))
) %>%
  filter(!is.na(gene)) %>%
  distinct()

build_violin_df <- function(gene){
  if(!(gene %in% rownames(norm_mat))){
    warning(paste("Gene not found:", gene))
    return(NULL)
  }
  
  idx <- grep(paste0("^", gene, "$"), rownames(norm_mat))
  if(length(idx) == 0){
    warning(paste("Gene not found:", gene))
    return(NULL)
  }
  
  expr <- norm_mat[idx, ]
  
  tmp <- tibble(
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
  tmp <- left_join(tmp, gene_class, by = "gene")
  
  tmp
}

all_targets <- unique(c(data$gene_a, data$gene_i))
all_targets <- all_targets[!is.na(all_targets)]

violin_df <- map_df(all_targets, build_violin_df)
violin_df <- violin_df %>%
  mutate(perturbed = factor(perturbed, levels = c("perturbed", "other")))

# colored facet labels
lab_fun <- function(g){
  cls <- gene_class$class[match(g, gene_class$gene)]
  col <- ifelse(cls == "a", "darkgreen", "darkred")
  paste0("<span style='color:", col, "'>", g, "</span>")
}

# Compute t-test per gene
stat_table <- violin_df %>%
  group_by(gene) %>%
  summarise(
    test = list(t.test(expr ~ perturbed)),
    p_value = test[[1]]$p.value,
    statistic = test[[1]]$statistic,
    direction = ifelse(statistic > 0, "perturbed higher", "other higher"),
    .groups = "drop"
  )

# Create a label: always p-value, add direction only if p < 0.05
stat_table <- stat_table %>%
  mutate(
    label = ifelse(
      p_value < 0.05,
      paste0("p=", signif(p_value, 2), "\n", direction),
      paste0("p=", signif(p_value, 2))
    )
  )

# Add annotations to the plot
p <- ggplot(violin_df, aes(x = perturbed, y = expr, fill = perturbed)) +
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
crispra_cells <- colData(cds)$type == "CRISPRa"
cds_crispra <- cds[, crispra_cells]

# Now get normalized counts only for CRISPRa cells
norm_mat <- normalized_counts(cds_crispra)

# Also subset metadata if you use it later
data_a = filter(data, type == "CRISPRa")

# Build the violin dataframe with guide info
build_violin_df_guide <- function(gene){
  if(!(gene %in% rownames(norm_mat))){
    warning(paste("Gene not found:", gene))
    return(NULL)
  }
  
  idx <- grep(paste0("^", gene, "$"), rownames(norm_mat))
  if(length(idx) == 0){
    warning(paste("Gene not found:", gene))
    return(NULL)
  }
  
  expr <- norm_mat[idx, ]
  
  tmp <- tibble(
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
  tmp <- left_join(tmp, gene_class, by = "gene")
  
  tmp
}

violin_df_guides <- map_df(unique(data$gene_i), build_violin_df_guide)
#violin_df_guides = separate(violin_df_guides, cell, into = c("cellID","sample", "guide_a", "guide_i"), sep = "\\.", remove = F, convert = T)
# Check number of pieces in the first column name
num_pieces <- length(strsplit(violin_df_guides$cell[1], "\\.")[[1]])

if(num_pieces == 8){
  violin_df_guides <- separate(violin_df_guides, cell, into = c("cellID","sample","guide_a","guide_i", "gene_a", "gene_i", "comb", "type"), sep = "\\.", remove = FALSE, convert = TRUE)
} else if(num_pieces == 10){
  violin_df_guides <- separate(
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
violin_df_guides <- violin_df_guides %>%
  mutate(perturbed = factor(perturbed, levels = c("perturbed", "other")))


violin_df_guides_filtered <- violin_df_guides %>%
  group_by(guide_a, perturbed) %>%
  filter(n() >= 20) %>%   # keep only groups with at least 10 observations per perturbed level
  ungroup() %>%
  group_by(guide_a) %>%
  filter(n_distinct(perturbed) == 2) %>%  # ensure both perturbed and control exist
  ungroup()

# Now compute t-test per guide_a safely
stat_table_guides <- violin_df_guides_filtered %>%
  group_by(guide_a) %>%
  summarise(
    test = list(t.test(expr ~ perturbed)),
    p_value = test[[1]]$p.value,
    statistic = test[[1]]$statistic,
    direction = ifelse(statistic > 0, "perturbed higher", "other higher"),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("p=", signif(p_value, 2), "\n", direction)
  )

# Compute max expression per guide for annotation
stat_table_guides <- stat_table_guides %>%
  left_join(
    violin_df_guides_filtered %>% 
      group_by(guide_a) %>% 
      summarise(max_expr = max(expr), .groups = "drop"),
    by = "guide_a"
  ) %>%
  mutate(y_pos = max_expr - 0.1 * max_expr)

# Plot
p_guides <- ggplot(violin_df_guides_filtered, aes(x = perturbed, y = expr, fill = perturbed)) +
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

### CRIPRa
# Subset cells in the cds object
crispri_cells <- colData(cds)$type == "CRISPRi"
cds_crispri <- cds[, crispri_cells]

# Now get normalized counts only for CRISPRa cells
norm_mat <- normalized_counts(cds_crispri)

# Also subset metadata if you use it later
data_i = filter(data, type == "CRISPRi")

# Build the violin dataframe with guide info
build_violin_df_guide <- function(gene){
  if(!(gene %in% rownames(norm_mat))){
    warning(paste("Gene not found:", gene))
    return(NULL)
  }
  
  idx <- grep(paste0("^", gene, "$"), rownames(norm_mat))
  if(length(idx) == 0){
    warning(paste("Gene not found:", gene))
    return(NULL)
  }
  
  expr <- norm_mat[idx, ]
  
  tmp <- tibble(
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
  tmp <- left_join(tmp, gene_class, by = "gene")
  
  tmp
}

violin_df_guides <- map_df(unique(data$gene_i), build_violin_df_guide)
#violin_df_guides = separate(violin_df_guides, cell, into = c("cellID","sample", "guide_a", "guide_i"), sep = "\\.", remove = F, convert = T)
# Check number of pieces in the first column name
num_pieces <- length(strsplit(violin_df_guides$cell[1], "\\.")[[1]])

if(num_pieces == 8){
  violin_df_guides <- separate(violin_df_guides, cell, into = c("cellID","sample","guide_a","guide_i", "gene_a", "gene_i", "comb", "type"), sep = "\\.", remove = FALSE, convert = TRUE)
} else if(num_pieces == 10){
  violin_df_guides <- separate(
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
violin_df_guides <- violin_df_guides %>%
  mutate(perturbed = factor(perturbed, levels = c("perturbed", "other")))

violin_df_guides_filtered <- violin_df_guides %>%
  group_by(guide_i, perturbed) %>%
  filter(n() >= 20) %>%   # keep only groups with at least 10 observations per perturbed level
  ungroup() %>%
  group_by(guide_i) %>%
  filter(n_distinct(perturbed) == 2) %>%  # ensure both perturbed and control exist
  ungroup()

# Now compute t-test per guide_a safely
stat_table_guides <- violin_df_guides_filtered %>%
  group_by(guide_i) %>%
  summarise(
    test = list(t.test(expr ~ perturbed)),
    p_value = test[[1]]$p.value,
    statistic = test[[1]]$statistic,
    direction = ifelse(statistic > 0, "perturbed higher", "other higher"),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("p=", signif(p_value, 2), "\n", direction)
  )

# Compute max expression per guide for annotation
stat_table_guides <- stat_table_guides %>%
  left_join(
    violin_df_guides_filtered %>% 
      group_by(guide_i) %>% 
      summarise(max_expr = max(expr), .groups = "drop"),
    by = "guide_i"
  ) %>%
  mutate(y_pos = max_expr - 0.1 * max_expr)

# Plot
p_guides <- ggplot(violin_df_guides_filtered, aes(x = perturbed, y = expr, fill = perturbed)) +
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
