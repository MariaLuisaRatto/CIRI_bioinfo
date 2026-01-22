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
library(biomaRt)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Error: No directory provided. Please supply the input directory path.")
}


dir <- args[1]

message(paste("Loading data from:", dir))

#Load 10X data
file = args[2]
data <- Read10X_h5(file.path(dir, file), unique.features = TRUE)
data_seurat <- CreateSeuratObject(counts = data$`Gene Expression`, assay = "RNA")

# Calculate percentages ribo mito
data_seurat[["percent.mt"]] <- PercentageFeatureSet(data_seurat, pattern = "^MT-")
data_seurat[["percent.rb"]] <- PercentageFeatureSet(data_seurat, pattern = "^RPS|^RPL")

# Extract metadata
qc_data <- data_seurat@meta.data
# Fix color scale range (use global limits based on the dataset)
color_limits <- range(qc_data$nCount_RNA, na.rm = TRUE)

# Scatter plot
p = ggplot(qc_data, aes(x = percent.mt, y = percent.rb, colour = nCount_RNA)) +
  geom_point(alpha = 0.6, size = 1) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Mitochondrial gene percentage",
    y = "Ribosomal gene percentage",
    title = "Mitochondrial vs Ribosomal content per cell",
    colour = "nCount_RNA"
  ) +
  scale_color_viridis_c(
    option = "turbo",
    trans = "log10",
    limits = color_limits,
    oob = scales::squish
  )
ggsave(p, filename = paste0(dir, "RiboMito.pdf"), width = 10, height = 10)

#filter ngenes x cell 
# Filter cells with < 250 detected genes
data_seurat <- subset(data_seurat, subset = nFeature_RNA >= 250)

# Calculate mitochondrial and ribosomal gene percentages
data_seurat[["percent.mt"]] <- PercentageFeatureSet(data_seurat, pattern = "^MT-")
data_seurat[["percent.rb"]] <- PercentageFeatureSet(data_seurat, pattern = "^RPS|^RPL")

# Extract metadata
qc_data <- data_seurat@meta.data

# Scatter plot
p <- ggplot(qc_data, aes(x = percent.mt, y = percent.rb, colour = nCount_RNA)) +
  geom_point(alpha = 0.6, size = 1) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Mitochondrial gene percentage",
    y = "Ribosomal gene percentage",
    title = "Mitochondrial vs Ribosomal content per cell",
    colour = "nCount_RNA"
  ) +
  scale_color_viridis_c(
    option = "turbo",
    trans = "log10",
    limits = color_limits,
    oob = scales::squish
  )

ggsave(p, filename = paste0(dir, "RiboMito_filtered.pdf"), width = 10, height = 10)

# Get protein-coding genes
mart <- useMart(
  biomart = "ensembl", 
  dataset = "hsapiens_gene_ensembl", 
  host = "https://sep2025.archive.ensembl.org" 
)
all_coding_genes <- getBM(
  mart = mart,
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "biotype",
  values = list(biotype = "protein_coding")
)

# Fallback if symbol is missing
all_coding_genes$mix <- ifelse(all_coding_genes$hgnc_symbol == "", all_coding_genes$ensembl_gene_id, all_coding_genes$hgnc_symbol)

# Filter to protein-coding genes only
data_seurat <- subset(data_seurat, features = intersect(rownames(data_seurat), all_coding_genes$mix))

# Extract and export matrix
mat <- GetAssayData(data_seurat, assay = "RNA", layer = "counts")

# --- Filters ---
#  Remove mitochondrial and ribosomal genes
keep_genes <- !grepl("^MT-|^RPS|^RPL", rownames(mat), ignore.case = TRUE)
mat <- mat[keep_genes, ]

# ⃣ Remove genes with <3 total UMIs across all cells
gene_totals <- Matrix::rowSums(mat)
mat <- mat[gene_totals >= 3, ]

# Optional sanity check
print(dim(mat))

gc()                      # garbage collect
rm(list = setdiff(ls(), c("mat", "dir")))  # remove everything except needed
gc()                      # collect again

mat <- as.data.frame(mat)
print(dim(mat))
mat = mat[, colSums(mat != 0) > 0]
print(dim(mat))
mat = as.data.frame(mat)
fwrite(mat, file = paste0(dir, "mingenes_proteincoding_filtered_feature_bc_matrix.csv"), row.names = TRUE)


# #library(rCASC)
# h5tocsv(
#   group = "docker",
#   file = "/30tb/3tb/data/ratto/AB012/Results_aggr/outs/count/filtered_feature_bc_matrix.h5",
#   type = "h5",
#   version = "5"
# )
# 
# mitoRiboUmi(
#   group = "docker",
#   scratch.folder = "/2tb/torino/ratto/scratch/",
#   file = "/2tb/torino/ratto/RiboMito/filtered_feature_bc_matrix.csv",
#   separator = ",",
#   gtf.name = "Homo_sapiens.GRCh38.101.gtf",
#   bio.type = "protein_coding",
#   umiXgene = 3
# )
# 
# scannobyGtf(
#   group = "docker",
#   file = "/30tb/3tb/data/ratto/AB012/Results_aggr/outs/count/filtered_feature_bc_matrix.csv",
#   gtf.name= "Homo_sapiens.GRCh38.110.gtf",
#   biotype = "protein_coding",
#   mt =F,
#   ribo.proteins =F,
#   umiXgene = 3,
#   riboStart.percentage = 0,
#   riboEnd.percentage = 100,
#   mitoStart.percentage = 1,
#   mitoEnd.percentage = 100,
#   thresholdGenes = 250
# )
# library(rCASC)
# scannobyGtf(
#   group = "docker",
#   file = "/2tb/torino/ratto/RiboMito/filtered_feature_bc_matrix.csv",
#   gtf.name= "Homo_sapiens.GRCh38.101.gtf",
#   biotype = "protein_coding",
#   mt =F,
#   ribo.proteins =F,
#   umiXgene = 3,
#   riboStart.percentage = 0,
#   riboEnd.percentage = 100,
#   mitoStart.percentage = 0,
#   mitoEnd.percentage = 20,
#   thresholdGenes = 250
# )

## see /30tb/3tb/data/ratto/rCASC

#exp = read.csv("./annotated_filtered_feature_bc_matrix.csv", row.names = 1)

old_names = c(colnames(mat))
old_names = sub("\\.", "-", old_names)
colnames(mat) = old_names
assigned_wide = read.csv(paste0(dir, "annotation_data.csv"))
mat = mat[, colnames(mat) %in% assigned_wide$cell_barcode]

# Create a named vector: old -> new
assigned_wide = mutate(assigned_wide, new_names = paste0(cell_barcode, "-", feature_a, "-", feature_i))
name_map <- setNames(assigned_wide$new_names, assigned_wide$cell_barcode)

# Replace column names in exp
colnames(mat) <- name_map[colnames(mat)]

#write.csv(mat, "annotated_matrix.csv")
fwrite(mat, file = paste0(dir, "annotated_matrix.csv"), row.names = TRUE)

