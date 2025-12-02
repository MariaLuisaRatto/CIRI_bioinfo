#secondary analysis load 
suppressMessages(library(devtools))
suppressMessages(library(monocle3))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(gtools))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(viridis))
suppressMessages(library(ggrepel))
set.seed(1234597698)


dir = "."
file = "annotated_matrix.csv"
#res = 0.5e-4
#for 1 guide
res = 0.1e-4
#for 2 guides
#exp$X = rownames(exp)

#LOAD ANNOTATED GENE EXP
#exp = read.csv("/30tb/3tb/data/ratto/testing/annotated_silencing_matrix_complete_all_samples.csv", header = T)
exp = read.csv(paste0(dir,"/", file), header = T)
rownames(exp) = exp$X
exp$X = NULL
#exp$X = rownames(exp)
# exp = exp[, colSums(exp != 0) > 0]
# exp <- separate(
#   exp,
#   col="X",
#   into=c("a","b"),
#   sep = ":",
#   remove = TRUE,
#   convert = FALSE,
#   extra = "warn",
#   fill = "warn"
# )
# exp$a <- NULL
# names <- make.unique(exp$b, sep=".")
# rownames(exp) <- names
# exp$b <- NULL


#CREATE ANNOTATION DATAFRAME
data = data.frame()
data = as.data.frame(colnames(exp))
colnames(data) = c("nomi")
rownames(data)=data$nomi

num_pieces <- length(strsplit(data$nomi, "\\.")[[1]])
if(num_pieces <= 4){
  data = separate(data, nomi, into = c("cellID","sample", "guide_a", "guide_i"), sep = "\\.", remove = F, convert = T)
  data = mutate(data, comb = paste(guide_a, guide_i, sep = "-"))
  data <- data %>%
    mutate(gene_a = sapply(strsplit(guide_a, "_"), `[`, 1))
  data <- data %>%
    mutate(gene_i = sapply(strsplit(guide_i, "_"), `[`, 1))
  data = mutate(data, gene_comb = paste(gene_a, gene_i, sep = "-"))
} else if(num_pieces  > 4){
  # Split into max 6 parts; if fewer, fill with NA instead of shifting
  data <- separate(
    data,
    nomi,
    into = c("cellID", "sample", "guide_a1", "guide_a2", "guide_i1", "guide_i2"),
    sep = "\\.",
    remove = FALSE,
    convert = TRUE,
    fill = "right"   # <-- important: pads missing with NA instead of shifting
  )
  
  data = mutate(data, comb = paste0(guide_a1,";", guide_a2,"-", guide_i1, ";", guide_i2))
  data <- data %>%
    mutate(gene_a = sapply(strsplit(guide_a1, "_"), `[`, 1))
  data <- data %>%
    mutate(gene_i = sapply(strsplit(guide_i1, "_"), `[`, 1))
  data = mutate(data, gene_comb = paste0(gene_a, "-", gene_i))
  
  data = mutate(data, guide_a = paste(guide_a1, guide_a2, sep = ";"))
  data = mutate(data, guide_i = paste(guide_i1, guide_i2, sep = ";"))
}

data <- data %>%
  mutate(
    type = case_when(
      !is.na(gene_a) & !is.na(gene_i) ~ "CIRI",
      !is.na(gene_a) &  is.na(gene_i) ~ "CRISPRa",
      is.na(gene_a) & !is.na(gene_i) ~ "CRISPRi",
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

p = ggplot(UMAP, aes(x, y, color = gene_a, alpha = 0.1)) +
  geom_point(size = 0.4)+ 
  facet_wrap(vars(type, gene_a))
ggsave(p, filename = paste0(dir,"/UMAP_gene_a_facet.pdf"),
       width = 18, height = 18)

p = ggplot(UMAP, aes(x, y, color = gene_i, alpha = 0.1)) +
  geom_point(size = 0.4)+ 
  facet_wrap(vars(type))
ggsave(p, filename = paste0(dir,"/UMAP_gene_i.pdf"),
       width = 18, height = 6)

p = ggplot(UMAP, aes(x, y, color = gene_i, alpha = 0.1)) +
  geom_point(size = 0.4)+ 
  facet_wrap(vars(type, gene_i))
ggsave(p, filename = paste0(dir,"/UMAP_gene_i_facet.pdf"),
       width = 18, height = 18)

#Plot by gene comb
p = ggplot(UMAP, aes(x, y, color = gene_i, alpha = 0.1)) +
  geom_point(size = 0.4)+ 
  facet_wrap(vars(gene_a))
ggsave(p, filename = paste0(dir,"/UMAP_facet_a_color_i.pdf"),
       width = 12, height = 9)

p = ggplot(UMAP, aes(x, y, color = gene_a, alpha = 0.1)) +
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

genes = "cc_genes.txt"
glist = c(unique(t(read.delim(paste0(dir, "/", genes), sep = ",", header = F))))

p = plot_cells(cds,
               genes = glist,
               label_cell_groups = T,
               show_trajectory_graph = FALSE)+
  labs(x = "UMAP 1", y = "UMAP 2", title = "")#+
#scale_color_viridis(option="E", discrete=F)

ggsave(p, filename = paste0(dir,"/UMAP_gene_expression_cc.pdf"),
       width = 10, height = 10)

genes = "retinoic_acid_genes.txt"
glist = c(unique(t(read.delim(paste0(dir, "/", genes), sep = ",", header = F))))

p = plot_cells(cds,
               genes = glist,
               label_cell_groups = T,
               show_trajectory_graph = FALSE)+
  labs(x = "UMAP 1", y = "UMAP 2", title = "")#+
#scale_color_viridis(option="E", discrete=F)

ggsave(p, filename = paste0(dir,"/UMAP_gene_expression_RA.pdf"),
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

#############

