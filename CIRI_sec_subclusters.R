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

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Error: No directory provided. Please supply the input directory path.")
}

dir = args[1]
ccs = c(strsplit(args[2], split = "-"))
#dir = "/Users/marialuisaratto/scripts/CIRI12_starting_dataset/"
#ccs = c("4")
group_of_interest = paste("Group", paste(ccs, collapse = "_"), sep = "_")
print(group_of_interest)
g = "SOX2"
g = args[3]

#file = args[4]
#load(paste0(dir, file))
load(paste0(dir, "/processed_cds.RData"))
clusters_to_keep = ccs
group_name = args[4]
#group_name = "muscle"
#res = 1e-3
res = as.numeric(args[5])
#res = 1e-4
#example command 
#Rscript ../CIRI_sec_subclusters.R /analysis/data/ 3 SOX2 muscle 1e-4
  
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
) +
  labs(caption = paste("res =", res))
ggsave(p_sub, filename = paste0(dir, "/UMAP_", group_name, "_subclusters.pdf"),
       width = 5, height = 4)

print("Subclustering done")


# ---- Extract UMAP coordinates ----
UMAP_sub <- as.data.frame(reducedDims(cds_sub)$UMAP)
names(UMAP_sub) <- c("x", "y")
UMAP_sub$nomi <- rownames(UMAP_sub)
num_pieces <- length(strsplit(UMAP_sub$nomi, "\\.")[[1]])
#AAACCAACAGGATTAA_30_AACCGGAG_SMAD2_TCTATT_1_SMAD2.2_78__0___singleSample

#if(num_pieces < 10){
  UMAP_sub <- separate(UMAP_sub, nomi, into = c("cellID","sample","guide_a","guide_i","comb", "gene_a", "gene_i", "gene_comb", "type"), sep = "\\.", remove = FALSE, convert = TRUE)
# } else if(num_pieces == 10){
#   # Split into max 6 parts; if fewer, fill with NA instead of shifting
#   UMAP_sub <- separate(
#     UMAP_sub,
#     nomi,
#     into = c("cellID", "sample", "guide_a1", "guide_a2", "guide_i1", "guide_i2","gene_a", "gene_i", "gene_comb", "type"),
#     sep = "\\.",
#     remove = FALSE,
#     convert = TRUE,
#     fill = "right"   # <-- important: pads missing with NA instead of shifting
#   )
#   
#   UMAP_sub = mutate(UMAP_sub, comb = paste0(guide_a1,";", guide_a2,"-", guide_i1, ";", guide_i2))
#   UMAP_sub <- UMAP_sub %>%
#     mutate(gene_a = sapply(strsplit(guide_a1, "_"), `[`, 1))
#   UMAP_sub <- UMAP_sub %>%
#     mutate(gene_i = sapply(strsplit(guide_i1, "_"), `[`, 1))
#   UMAP_sub = mutate(UMAP_sub, gene_comb = paste0(gene_a, "-", gene_i))
#   
#   UMAP_sub = mutate(UMAP_sub, guide_a = paste(guide_a1, guide_a2, sep = ";"))
#   UMAP_sub = mutate(UMAP_sub, guide_i = paste(guide_i1, guide_i2, sep = ";"))
# }


UMAP_sub$sample <- as.factor(UMAP_sub$sample)
print(head(UMAP_sub))

print("UMAPs...")
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
print("Plotting gene expression...")
genes = "genes_of_interest.txt"
glist <- readLines(paste0(dir, "/", genes))
glist <- trimws(glist)
glist <- glist[glist != ""]
glist <- unique(glist)
#put gene names in $gene_short_name to solve bug
old_names <- rownames(cds_sub)
new_names <- make.unique(old_names)

rownames(cds_sub) <- new_names
rowData(cds_sub)$gene_short_name <- new_names
#rowData(cds_sub)$gene_name <- rownames(cds_sub)
#rowData(cds_sub)$gene_short_name <- rowData(cds_sub)$gene_name

valid_genes <- intersect(
  glist,
  rowData(cds_sub)$gene_short_name
)

valid_genes  # sanity check
any(duplicated(colData(cds_sub)$clusters_sub))

colData(cds_sub)$clusters_sub <- factor(
  as.character(colData(cds_sub)$clusters_sub)
)

p = plot_cells(cds_sub,
               genes = valid_genes,
               label_cell_groups = T,
               show_trajectory_graph = FALSE)+
  labs(x = "UMAP 1", y = "UMAP 2", title = "")#+
#scale_color_viridis(option="E", discrete=F)

ggsave(p, filename = paste0(dir,"/UMAP_gene_expression_", group_name, ".pdf"),
       width = 10, height = 10)

# p = plot_cells(cds_sub,
#                genes = c(unique(UMAP_sub$gene_a), unique(UMAP_sub$gene_i)),
#                label_cell_groups = T,
#                show_trajectory_graph = FALSE)+
#   labs(x = "UMAP 1", y = "UMAP 2", title = "")
# 
# ggsave(p, filename = paste0(dir,"/UMAP_gene_expression_CRISPR_", group_name, ".pdf"),
#        width = 10, height = 10)

## ---- PSEUDOTIME ----
print("Running pseudotime...")
cds_sub <- learn_graph(cds_sub)

# Helper function to identify the root principal node with the lowest TTN expression
get_root_lowest_TTN <- function(cds, gene) {
    
    # sanity check
    if (!gene %in% rownames(cds)) {
      stop(paste("Gene", gene, "not found in CDS"))
    }
    
    gene_expr <- as.numeric(exprs(cds)[gene, ])
    names(gene_expr) <- colnames(cds)  # ensure names exist
    
    if (all(is.na(gene_expr))) stop("Gene expression is all NA")
    
    max_val <- max(gene_expr, na.rm = TRUE)
    min_cells <- names(gene_expr[gene_expr == max_val])
    
    if (length(min_cells) == 0) stop("No cells found with max expression")
    
    # Find the principal graph projection (which vertex each cell maps to)
    closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    vertices <- closest_vertex[min_cells, , drop = FALSE]
    
    # Find the most common vertex among these cells
    root_pr_node <- names(which.max(table(vertices)))
    root_pr_node = paste0("Y_", root_pr_node)
    print(root_pr_node)
    return(root_pr_node)
}

cds_sub <- order_cells(
  cds_sub,
  root_pr_nodes = get_root_lowest_TTN(cds_sub, gene = g)
)
# Save dataas an RData file
save(cds_sub, file = paste0(dir, "/processed_cds_", group_name, ".RData"))

#pseudotime for each sample 
samples <- unique(colData(cds_sub)$sample)
for (s in samples) {
  
  message(paste("Processing sample:", s))
  
  ## ---- subset CDS ----
  cds_sample <- cds_sub[, colData(cds_sub)$sample %in% s]
  group_name_s <- paste0(group_name, "_", s)
  
  save(cds_sample, file = paste0(dir, "/cds_sample_", s, "_ordered.RData"))
  
  pt <- as.data.frame(pseudotime(cds_sample))
  colnames(pt) <- "pseudotime"
  write.csv(
    pt,
    file = paste0(dir, "/pseudotime_", group_name_s, ".csv")
  )
  
  p_sub <- plot_cells(
    cds_sample,
    color_cells_by = "pseudotime",
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    graph_label_size = 1.5
  ) +
    ggtitle(paste("Pseudotime - sample", s))
  
  ggsave(
    p_sub,
    filename = paste0(dir, "/UMAP_", group_name_s, "_pseudotime.pdf"),
    width = 5,
    height = 4
  )
  
  ## ---- ENRICHMENT / DEPLETION ----
  df_sub <- as.data.frame(colData(cds_sample))
  
  ## ---- Gene A ----
  gene_a_counts <- df_sub %>%
    group_by(gene_a, clusters_sub) %>%
    summarise(count = n(), .groups = "drop")
  
  gene_a_percentages <- gene_a_counts %>%
    group_by(gene_a) %>%
    mutate(
      total = sum(count),
      percentage = count / total * 100
    ) %>%
    ungroup()
  
  p_gene_a <- ggplot(
    gene_a_percentages,
    aes(x = gene_a, y = percentage, fill = clusters_sub)
  ) +
    geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
    geom_text(
      aes(y = 100, label = total),
      #inherit.aes = FALSE,
      size = 3,
      angle = 90
    ) +
    #facet_grid(~gene_a) +
    labs(
      title = paste("Subclusters within gene_a - sample", s),
      x = "gene_a",
      y = "Percentage of cells",
      fill = "Subcluster"
    ) +
    theme_minimal() +
    theme(
      strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    scale_fill_manual(
      values = lighten(
        viridis::turbo(length(unique(gene_a_percentages$clusters_sub))),
        amount = 0.4
      )
    )
  
  ggsave(
    p_gene_a,
    filename = paste0(dir, "/subclusters_gene_a_", group_name_s, ".pdf"),
    width = 10,
    height = 5
  )
  
  ## ---- guide_a ----
  guide_a_counts <- df_sub %>%
    group_by(guide_a, clusters_sub) %>%
    summarise(count = n(), .groups = "drop")
  
  guide_a_percentages <- guide_a_counts %>%
    group_by(guide_a) %>%
    mutate(
      total = sum(count),
      percentage = count / total * 100
    ) %>%
    ungroup()
  
  p_guide_a <- ggplot(
    guide_a_percentages,
    aes(x = guide_a, y = percentage, fill = clusters_sub)
  ) +
    geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
    geom_text(
      aes(y = 100, label = total),
      #inherit.aes = FALSE,
      size = 3,
      angle = 90
    ) +
    #facet_grid(~guide_a) +
    labs(
      title = paste("Subclusters within guide_a - sample", s),
      x = "guide_a",
      y = "Percentage of cells",
      fill = "Subcluster"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    scale_fill_manual(
      values = lighten(
        viridis::turbo(length(unique(guide_a_percentages$clusters_sub))),
        amount = 0.4
      )
    )
  
  ggsave(
    p_guide_a,
    filename = paste0(dir, "/subclusters_guide_a_", group_name_s, ".pdf"),
    width = 10,
    height = 5
  )
  
  ## ---- Gene I ----
  gene_i_counts <- df_sub %>%
    group_by(gene_i, clusters_sub) %>%
    summarise(count = n(), .groups = "drop")
  
  gene_i_percentages <- gene_i_counts %>%
    group_by(gene_i) %>%
    mutate(
      total = sum(count),
      percentage = count / total * 100
    ) %>%
    ungroup()
  
  p_gene_i <- ggplot(
    gene_i_percentages,
    aes(x = gene_i, y = percentage, fill = clusters_sub)
  ) +
    geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
    geom_text(
      aes(y = 100, label = total),
      #inherit.aes = FALSE,
      size = 3,
      angle = 90
    ) +
    #facet_grid(~gene_i) +
    labs(
      title = paste("Subclusters within gene_i - sample", s),
      x = "gene_i",
      y = "Percentage of cells",
      fill = "Subcluster"
    ) +
    theme_minimal() +
    theme(
      strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    scale_fill_manual(
      values = lighten(
        viridis::turbo(length(unique(gene_i_percentages$clusters_sub))),
        amount = 0.4
      )
    )
  
  ggsave(
    p_gene_i,
    filename = paste0(dir, "/subclusters_gene_i_", group_name_s, ".pdf"),
    width = 10,
    height = 5
  )
  
  ## ---- guide_i ----
  guide_i_counts <- df_sub %>%
    group_by(guide_i, clusters_sub) %>%
    summarise(count = n(), .groups = "drop")
  
  guide_i_percentages <- guide_i_counts %>%
    group_by(guide_i) %>%
    mutate(
      total = sum(count),
      percentage = count / total * 100
    ) %>%
    ungroup()
  
  p_guide_i <- ggplot(
    guide_i_percentages,
    aes(x = guide_i, y = percentage, fill = clusters_sub)
  ) +
    geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
    geom_text(
      aes(y = 100, label = total),
      #inherit.aes = FALSE,
      size = 3,
      angle = 90
    ) +
    #facet_grid(~guide_i) +
    labs(
      title = paste("Subclusters within guide_i - sample", s),
      x = "guide_i",
      y = "Percentage of cells",
      fill = "Subcluster"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    scale_fill_manual(
      values = lighten(
        viridis::turbo(length(unique(guide_i_percentages$clusters_sub))),
        amount = 0.4
      )
    )
  
  ggsave(
    p_guide_i,
    filename = paste0(dir, "/subclusters_guide_i_", group_name_s, ".pdf"),
    width = 10,
    height = 5
  )
}


# 
# # Then use it when ordering cells:
# cds_sub <- order_cells(cds_sub, root_pr_nodes = get_root_lowest_TTN(cds_sub))
# pt = as.data.frame(pseudotime(cds_sub))
# names(pt) = c("pseudotime")
# write.csv(pt, paste0(dir, "/pseudotime_", group_name, ".csv"))
# p_sub = plot_cells(cds_sub,
#                    color_cells_by = "pseudotime",
#                    label_cell_groups=FALSE,
#                    label_leaves=FALSE,
#                    label_branch_points=FALSE,
#                    graph_label_size=1.5)
# ggsave(p_sub, filename = paste0(dir, "/UMAP_", group_name, "_pseudotime.pdf"),
#        width = 5, height = 4)
# 
# # ---- ENRICHMENT / DEPLETION (subclusters) ----
# print("Running enrichment/depletion (subclusters)...")
# df_sub <- as.data.frame(colData(cds_sub))
# cellsxsubcluster <- df_sub %>%
#   distinct(clusters_sub, nomi) %>%
#   group_by(clusters_sub) %>%
#   mutate(cellsxsubcluster = n()) %>%
#   distinct(clusters_sub, cellsxsubcluster)
# 
# # Gene A
# gene_a_subcluster_counts <- df_sub %>%
#   group_by(sample, gene_a, clusters_sub) %>%
#   summarise(count = n(), .groups = "drop")
# gene_a_subcluster_percentages <- gene_a_subcluster_counts %>%
#   group_by(sample, gene_a) %>%
#   mutate(total = sum(count),
#          percentage = (count / total) * 100) %>%
#   ungroup()
# 
# p7 <- ggplot(gene_a_subcluster_percentages, 
#              aes(x = factor(sample), y = percentage, fill = clusters_sub)) +
#   geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
#   # Add total cell counts above the bar
#   geom_text(aes(x = factor(sample), y = 100, label = total), 
#             inherit.aes = FALSE, size = 3, angle = 90) +
#   theme_minimal() +
#   theme(strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), paper = "white") +
#   facet_grid(~gene_a) +
#   labs(x = "Sample", y = "Percentage of cells", fill = "Subcluster")+
#   scale_fill_manual(values = lighten(viridis::turbo(length(unique(gene_a_subcluster_percentages$clusters_sub))), amount = 0.4))
# ggsave(p7, filename = paste0(dir, "/subclusters_within_gene_a_", group_name, ".pdf"), 
#        width = 10, height = 5)
# 
# # ---- By guide_a ----
# guide_a_subcluster_counts <- df_sub %>%
#   group_by(sample, guide_a, clusters_sub) %>%
#   summarise(count = n(), .groups = "drop")
# 
# guide_a_subcluster_percentages <- guide_a_subcluster_counts %>%
#   group_by(sample, guide_a) %>%
#   mutate(total = sum(count),
#          percentage = (count / total) * 100) %>%
#   ungroup()
# 
# p_guide_a <- ggplot(guide_a_subcluster_percentages, 
#                     aes(x = factor(sample), y = percentage, fill = clusters_sub)) +
#   geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
#   # Add total cell counts above the bar
#   geom_text(aes(x = factor(sample), y = 100, label = total), 
#             inherit.aes = FALSE, size = 3, angle = 90) +
#   theme_minimal() +
#   theme(strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), paper = "white") +
#   facet_grid(~guide_a) +
#   labs(x = "Sample", y = "Percentage of cells", fill = "Subcluster")+
#   scale_fill_manual(values = lighten(viridis::turbo(length(unique(guide_a_subcluster_percentages$clusters_sub))), amount = 0.4))
# 
# ggsave(p_guide_a, filename = paste0(dir, "/subclusters_within_guide_a_", group_name, ".pdf"), 
#        width = 10, height = 5)
# 
# # Gene I
# gene_i_subcluster_counts <- df_sub %>%
#   group_by(sample, gene_i, clusters_sub) %>%
#   summarise(count = n(), .groups = "drop")
# gene_i_subcluster_percentages <- gene_i_subcluster_counts %>%
#   group_by(sample, gene_i) %>%
#   mutate(total = sum(count),
#          percentage = (count / total) * 100) %>%
#   ungroup()
# 
# p8 <- ggplot(gene_i_subcluster_percentages, 
#              aes(x = factor(sample), y = percentage, fill = clusters_sub)) +
#   geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
#   # Add total cell counts above the bar
#   geom_text(aes(x = factor(sample), y = 100, label = total), 
#             inherit.aes = FALSE, size = 3, angle = 90) +
#   theme_minimal() +
#   facet_grid(~gene_i) +
#   theme(strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), paper = "white") +
#   labs(x = "Sample", y = "Percentage of cells", fill = "Subcluster")+
#   scale_fill_manual(values = lighten(viridis::turbo(length(unique(gene_i_subcluster_percentages$clusters_sub))), amount = 0.4))
# ggsave(p8, filename = paste0(dir, "/subclusters_within_gene_i_", group_name, ".pdf"), 
#        width = 10, height = 5)
# 
# # ---- By guide_i ----
# guide_i_subcluster_counts <- df_sub %>%
#   group_by(sample, guide_i, clusters_sub) %>%
#   summarise(count = n(), .groups = "drop")
# 
# guide_i_subcluster_percentages <- guide_i_subcluster_counts %>%
#   group_by(sample, guide_i) %>%
#   mutate(total = sum(count),
#          percentage = (count / total) * 100) %>%
#   ungroup()
# 
# p_guide_i <- ggplot(guide_i_subcluster_percentages, 
#                     aes(x = factor(sample), y = percentage, fill = clusters_sub)) +
#   geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
#   # Add total cell counts above the bar
#   geom_text(aes(x = factor(sample), y = 100, label = total), 
#             inherit.aes = FALSE, size = 3, angle = 90) +
#   theme_minimal() +
#   theme(strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), paper = "white") +
#   facet_grid(~guide_i) +
#   labs(x = "Sample", y = "Percentage of cells", fill = "Subcluster")+
#   scale_fill_manual(values = lighten(viridis::turbo(length(unique(guide_i_subcluster_percentages$clusters_sub))), amount = 0.4))
# 
# ggsave(p_guide_i, filename = paste0(dir, "/subclusters_within_guide_i_", group_name, ".pdf"), 
#        width = 10, height = 5)


# Usage
#UMAP_muscle <- analyze_clusters(cds, clusters_to_keep = ccs, group_name = group_of_interest, dir = dir, res = 0.05e-2)

#UMAP_noRA <- analyze_clusters(cds, clusters_to_keep = c("5","2"), group_name = "noRA", dir = dir, res =  0.05e-2)

#UMAP_RA <- analyze_clusters(cds, clusters_to_keep = c("1","3"), group_name = "RA", dir = dir, res =  0.05e-2)

#UMAP_misc <- analyze_clusters(cds, clusters_to_keep = c("6","7","8","9"), group_name = "misc", dir = dir, res =  0.05e-2)
