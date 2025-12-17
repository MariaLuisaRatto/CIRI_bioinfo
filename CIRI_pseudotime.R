#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#ARGS
#dir = "."
dir = as.character(args[1])
print(paste("Folder:", dir))
#file = "cds_sample_1_ordered.RData"
file = args[2]
print(paste("File:", file))
#pseudotime = "pseudotime_muscle_1.csv"
#pseudotime = "signature_values_csv/sarcomere_core.csv"
pseudotime = args[3]
print(paste("Pseudotime file:", pseudotime))
#analysis = "pseudotime_sample_1"
analysis = args[4]
all = args[5]
#all = T
#all = as.logical(args[4])
print(paste("All:", all))

print("Starting...")


suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(gtools))
suppressMessages(library(dplyr))
suppressMessages(library(monocle3))
suppressMessages(library(viridis))
suppressMessages(library(colorspace))

if (file.exists(paste0(dir,"/", file)) == F) {
  stop(paste("Error: cds file does not exist at path:", paste0(dir,"/", file)))
}

if (file.exists(paste0(dir, "/", pseudotime)) == F) {
  stop(paste("Error: Pseudotime ile does not exist at path:", paste0(dir,"/", pseudotime)))
}



# Create a temporary environment to load the data into
temp_env <- new.env()
load(paste0(dir, "/", file), envir = temp_env)

# Grab the first object in that environment and assign it to 'cds'
cds <- temp_env[[ls(temp_env)[1]]]
pt = read.csv(paste0(dir, "/", pseudotime), row.names = 1)
analysis = names(pt)
names(pt) = c("pseudotime")

#CUMOLATIVE FREQ
pt$nomi = rownames(pt)
num_pieces <- length(strsplit(pt$nomi, "\\.")[[1]])
#AAACCAACAGGATTAA_30_AACCGGAG_SMAD2_TCTATT_1_SMAD2.2_78__0___singleSample

if(num_pieces == 8){
  pt <- separate(pt, nomi, into = c("cellID","sample","guide_a","guide_i", "gene_a", "gene_i", "gene_comb", "type"), sep = "\\.", remove = FALSE, convert = TRUE)
} else if(num_pieces == 10){
  pt <- separate(
    pt,
    nomi,
    into = c("cellID", "sample", "guide_a1", "guide_a2", "guide_i1", "guide_i2", "gene_a", "gene_i", "gene_comb", "type"),
    sep = "\\.",
    remove = FALSE,
    convert = TRUE,
    fill = "right"
  )
  
  pt = mutate(pt, guide_a = paste(guide_a1, guide_a2, sep = ";"))
  pt = mutate(pt, guide_i = paste(guide_i1, guide_i2, sep = ";"))
}

pt$sample = as.factor(pt$sample)


# gene_a
if (length(unique(pt$gene_a)) < 21) {
  ggplot(pt, aes(pseudotime, colour = gene_a)) + 
    stat_ecdf(geom = "step") +
    theme_classic(base_size = 14) +
    ylab("Cumulative frequency") +
    xlab(analysis) +
    ggtitle(paste("Cumulative frequency on", analysis))+
    scale_color_manual(values = viridis::turbo(length(unique(pt$gene_a))))+
    facet_grid(vars(type))
  ggsave(filename = paste0(dir, "/cumulative_frequency_", analysis, "_gene_a.pdf"),
         width = 10, height = 14)
} else {
  print("Too many genes_a to plot cumulative_frequency_pseudotime_gene_a.pdf")
}

# gene_i
if (length(unique(pt$gene_i)) < 21) {
  ggplot(pt, aes(pseudotime, colour = gene_i)) + 
    stat_ecdf(geom = "step") +
    theme_classic(base_size = 14) +
    ylab("Cumulative frequency") +
    xlab(analysis) +
    ggtitle(paste("Cumulative frequency on", analysis))+
    scale_color_manual(values = viridis::turbo(length(unique(pt$gene_i))))+
    facet_grid(vars(type))
  ggsave(filename = paste0(dir, "/cumulative_frequency_", analysis, "_gene_i.pdf"),
         width = 10, height = 14)
} else {
  print("Too many genes_i to plot cumulative_frequency_pseudotime_gene_i.pdf")
}

# guide_a
if (length(unique(pt$guide_a)) < 21) {
  ggplot(pt, aes(pseudotime, colour = guide_a)) + 
    stat_ecdf(geom = "step") +
    theme_classic(base_size = 14) +
    ylab("Cumulative frequency") +
    xlab(analysis) +
    ggtitle(paste("Cumulative frequency on", analysis)) +
    scale_color_manual(values = viridis::turbo(length(unique(pt$guide_a))))+
    facet_grid(vars(type))
  ggsave(filename = paste0(dir, "/cumulative_frequency_", analysis, "_guide_a.pdf"),
         width = 7, height = 10)
} else {
  print("Too many guide_a to plot cumulative_frequency_pseudotime_guide_a.pdf")
}

# guide_i
if (length(unique(pt$guide_i)) < 21) {
  ggplot(pt, aes(pseudotime, colour = guide_i)) + 
    stat_ecdf(geom = "step") +
    theme_classic(base_size = 14) +
    ylab("Cumulative frequency") +
    xlab(analysis) +
    ggtitle(paste("Cumulative frequency on", analysis)) +
    scale_color_manual(values = viridis::turbo(length(unique(pt$guide_i))))+
    facet_grid(vars(type))
  ggsave(filename = paste0(dir, "/cumulative_frequency_", analysis, "_guide_i.pdf"),
         width = 7, height = 10)
} else {
  print("Too many guide_i to plot cumulative_frequency_pseudotime_guide_i.pdf")
}


# --- KS TEST based on control samples ---

control_genes_a = c("NTCa", "NA")
control_genes_i = c("NTCi", "NA")
control_guides_a = control_genes_a
control_guides_i = control_genes_i
pt$gene_a[is.na(pt$gene_a)] <- "NA"
pt$gene_i[is.na(pt$gene_i)] <- "NA"
pt$guide_a[is.na(pt$guide_a)] <- "NA"
pt$guide_i[is.na(pt$guide_i)] <- "NA"


run_ks_test <- function(data, type, target, control_genes, dir, all, t_label) {
  cat("\n--- Running KS test for:", target, "in type:", type, "(", t_label, ") ---\n")
  
  ks_stat <- data.frame(stringsAsFactors = FALSE)
  
  # Ensure target column is atomic
  data <- data %>% mutate(!!sym(target) := as.character(!!sym(target)))
  
  # Filter by type
  tmp2 <- data %>% filter(type == !!t_label)
  
  # Skip if not enough unique values
  n_unique <- length(unique(tmp2[[target]]))
  if (n_unique <= 1) {
    message("Skipping ", target, " in ", t_label, ": only ", n_unique, " unique value.")
    return(NULL)
  }
  
  for (cval in unique(tmp2[[target]])) {
    tryCatch({
      tmp_TET <- tmp2 %>% filter(!!sym(target) == cval)
      if (nrow(tmp_TET) < 2 || all(is.na(tmp_TET$pseudotime))) next
      if (length(unique(tmp_TET$pseudotime)) <= 1) next
      
      # Compare vs each control gene (if all == TRUE)
      if (isTRUE(all)) {
        for (CG in control_genes) {
          tryCatch({
            if (CG == "NA") {
              # Take "NA" controls from all types
              tmp_CTR <- data %>% filter(!!sym(target) == "NA")
            } else {
              # Take control from same type
              tmp_CTR <- tmp2 %>% filter(!!sym(target) == CG)
            }
            
            if (nrow(tmp_CTR) < 2 || all(is.na(tmp_CTR$pseudotime))) next
            if (length(unique(tmp_CTR$pseudotime)) <= 1) next
            
            ks_g <- ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative = "greater")
            ks_l <- ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative = "less")
            
            ks_stat <- rbind(
              ks_stat,
              data.frame(
                gene = cval,
                pval = ks_g$p.value,
                statistic = -ks_g$statistic[["D^+"]],
                analysis = paste0("vs", CG),
                type = "greater"
              ),
              data.frame(
                gene = cval,
                pval = ks_l$p.value,
                statistic = ks_l$statistic[["D^-"]],
                analysis = paste0("vs", CG),
                type = "less"
              )
            )
          }, error = function(e) message("Skipping ", cval, " vs ", CG, ": ", e$message))
        }
      }
      
      # Compare vs all controls combined
      if (length(control_genes) >= 1) {
        if ("NA" %in% control_genes) {
          tmp_CTR <- data %>% filter(!!sym(target) == "NA")
          tmp_CTR2 <- tmp2 %>% filter(!!sym(target) %in% setdiff(control_genes, "NA"))
          tmp_CTR <- rbind(tmp_CTR, tmp_CTR2)
        } else {
          tmp_CTR <- tmp2 %>% filter(!!sym(target) %in% control_genes)
        }
        
        if (nrow(tmp_CTR) >= 2 && !all(is.na(tmp_CTR$pseudotime))) {
          ks_g <- ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative = "greater")
          ks_l <- ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative = "less")
          
          ks_stat <- rbind(
            ks_stat,
            data.frame(
              gene = cval,
              pval = ks_g$p.value,
              statistic = -ks_g$statistic[["D^+"]],
              analysis = paste0("vsControl_", target),
              type = "greater"
            ),
            data.frame(
              gene = cval,
              pval = ks_l$p.value,
              statistic = ks_l$statistic[["D^-"]],
              analysis = paste0("vsControl_", target),
              type = "less"
            )
          )
        }
      }
      
    }, error = function(e) message("ERROR for ", cval, ": ", e$message))
  }
  
  if (nrow(ks_stat) == 0) return(NULL)
  
  # Adjust p-values
  ks_stat$padj <- p.adjust(ks_stat$pval, method = "BH")
  ks_stat$sig <- ks_stat$padj <= 0.05
  
  ks_stat_fil <- ks_stat %>%
    group_by(gene, analysis) %>%
    filter(padj == min(padj)) %>%
    ungroup() %>%
    distinct()
  
  # Output files
  out_csv <- file.path(dir, paste0(t_label, "_", target, "_ks_statistics_all", all, analysis, ".csv"))
  out_pdf <- file.path(dir, paste0(t_label, "_volcano_ks_stats_", target, "_all", all, analysis, ".pdf"))
  
  write.csv(ks_stat, out_csv, row.names = FALSE)
  message("Wrote: ", out_csv)
  
  # Volcano plot
  plot_title <- paste("KS test statistics for", target, "- type:", t_label, "_", analysis)
  
  p <- ggplot(ks_stat_fil, aes(x = statistic, y = -log10(padj))) +
    geom_point(aes(colour = analysis), size = 2, alpha = 0.3) +
    ggrepel::geom_text_repel(aes(label = gene), size = 4, max.overlaps = 20) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "red") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    ggtitle(plot_title) +
    theme_minimal() +
    labs(x = "KS statistic", y = "-Log10 adj p-value") +
    scale_color_viridis_d()
  
  ggsave(filename = out_pdf, plot = p, width = 14, height = 8)
  message("Saved: ", out_pdf)
  
  return(ks_stat)
}
# ---- Main loop ----

targets <- c("gene_a", "gene_i", "guide_a", "guide_i")

for (t in unique(pt$type)) {
  cat("\n==============================\n")
  cat("Processing type:", t, "\n")
  cat("==============================\n")
  
  for (target in targets) {
    control_set <- if (grepl("_a$", target)) control_genes_a else control_genes_i
    run_ks_test(data = pt, type = t, target = target,
                control_genes = control_set, dir = dir,
                all = all, t_label = t)
  }
}

#CIRI guide comb volcano plot
pt = mutate(pt, guide_comb = paste(guide_a, guide_i, sep = "-"))
pt_CIRI = filter(pt, type == "CIRI")
# Remove guide_comb groups with fewer than 10 cells inside CIRI
pt_CIRI = pt_CIRI %>%
  group_by(guide_comb, sample) %>%
  filter(n() >= 8) %>%
  ungroup()

message("Remaining guide_comb in CIRI after filtering: ",
        length(unique(pt_CIRI$guide_comb)))

###############################################
# --- CUMULATIVE FREQUENCY PER SAMPLE (guide_comb) ---
###############################################

samples = unique(pt$sample)

for (s in samples) {
  tmp = filter(pt_CIRI, sample == s)
  
  p = ggplot(tmp, aes(pseudotime, colour = guide_i)) +
      stat_ecdf(geom = "step") +
      theme_classic(base_size = 14) +
      ylab("Cumulative frequency") +
      xlab("monocle pseudotime") +
      ggtitle(paste("Cumulative frequency on", analysis, "sample", s)) +
      facet_grid(vars(guide_a))+
      scale_color_manual(values = viridis::turbo(length(unique(tmp$guide_i)))) 
    
    ggsave(
      filename = paste0(dir, "/cumulative_frequency_", analysis, "_guide_comb_facet_sample_", s, ".pdf"),
      plot = p,
      width = 8, height = 15
    )
  
    print(length(unique(tmp$guide_comb)))
    #if (length(unique(tmp$guide_comb)) < 26) {
      p = ggplot(tmp, aes(pseudotime, colour = guide_i, linetype = guide_a)) +
        stat_ecdf(geom = "step") +
        theme_classic(base_size = 14) +
        ylab("Cumulative frequency") +
        xlab("monocle pseudotime") +
        ggtitle(paste("Cumulative frequency on", analysis, "sample", s)) +
        scale_color_manual(values = viridis::turbo(length(unique(tmp$guide_i))))
      
      ggsave(
        filename = paste0(dir, "/cumulative_frequency_", analysis, "_guide_comb_sample_", s, ".pdf"),
        plot = p,
        width = 14, height = 10
      )
   # } else {
    #  message("Too many guide_comb for sample ", s)
    #}
}

###############################################
# --- KS-TEST and VOLCANO for guide_comb (CIRI only) ---
###############################################

control_comb = c("NTCa-NTCi")

ks_df = data.frame()

for (gc in unique(pt_CIRI$guide_comb)) {
  tmp_case = filter(pt_CIRI, guide_comb == gc)
  
  if (nrow(tmp_case) < 2 || length(unique(tmp_case$pseudotime)) <= 1)
    next
  
  tmp_ctrl = filter(pt_CIRI, guide_comb %in% control_comb)
  
  if (nrow(tmp_ctrl) < 2 || length(unique(tmp_ctrl$pseudotime)) <= 1)
    next
  
  ks_g = ks.test(tmp_case$pseudotime, tmp_ctrl$pseudotime, alternative = "greater")
  ks_l = ks.test(tmp_case$pseudotime, tmp_ctrl$pseudotime, alternative = "less")
  
  ks_df = rbind(
    ks_df,
    data.frame(
      guide_comb = gc,
      pval = ks_g$p.value,
      stat = -ks_g$statistic[["D^+"]],
      direction = "greater"
    ),
    data.frame(
      guide_comb = gc,
      pval = ks_l$p.value,
      stat = ks_l$statistic[["D^-"]],
      direction = "less"
    )
  )
}

if (nrow(ks_df) > 0) {
  ks_df$padj = p.adjust(ks_df$pval, method = "BH")
  ks_df$neglog = -log10(ks_df$padj)
  ks_df$sig = ks_df$padj <= 0.05
  
  write.csv(
    ks_df,
    paste0(dir, "/CIRI_guide_comb_ks_statistics.csv"),
    row.names = FALSE
  )
  
  library(ggrepel)
  
  p_volc = ggplot(ks_df, aes(x = stat, y = neglog)) +
    geom_point(aes(colour = direction), alpha = 0.5, size = 2) +
    geom_text_repel(aes(label = guide_comb), size = 3, max.overlaps = 20) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "red") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    theme_minimal(base_size = 14) +
    labs(
      title = "KS volcano plot for guide_comb (CIRI)",
      x = "KS statistic",
      y = "-log10 adj p-value"
    ) +
    scale_color_viridis_d()
  
  ggsave(
    filename = paste0(dir, "/CIRI_volcano_guide_comb.pdf"),
    plot = p_volc,
    width = 12, height = 7
  )
}


# --- KS TEST based on control genes_a and genes_i ---
ks_stat = data.frame()
tmp2 = pt

# loop for gene_a
for (c in unique(pt$gene_a)) {
  tryCatch({
    print(c)
    tmp_TET = filter(tmp2, gene_a == c)
    
    if (all == T) {
      for (CG in control_genes_a) {
        tryCatch({
          print(CG)
          tmp_CTR = filter(tmp2, gene_a == CG)
          
          ks_g = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="greater")
          ks_stat = rbind(ks_stat, data.frame(gene = c, pval = ks_g$p.value, statistic = -ks_g$statistic[["D^+"]], analysis = paste0("vs", CG), type = "greater"))
          
          ks_l = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="less")
          ks_stat = rbind(ks_stat, data.frame(gene = c, pval = ks_l$p.value, statistic = ks_l$statistic[["D^-"]], analysis = paste0("vs", CG), type = "less"))
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }
    }
    
    # control together
    tmp_CTR = filter(tmp2, gene_a %in% control_genes_a)
    ks_g = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="greater")
    ks_l = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="less")
    
    ks_stat = rbind(ks_stat, 
                    data.frame(gene = c, pval = ks_g$p.value, statistic = -ks_g$statistic[["D^+"]], analysis="vsControlGenes_a", type="greater"),
                    data.frame(gene = c, pval = ks_l$p.value, statistic = ks_l$statistic[["D^-"]], analysis="vsControlGenes_a", type="less"))
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

ks_stat$padj <- p.adjust(ks_stat$pval, method = "BH")
ks_stat$sig <- ks_stat$padj <= 0.05

ks_stat_fil = ks_stat %>%
  group_by(gene, analysis) %>%
  filter(padj == min(padj)) %>%
  ungroup() %>% 
  distinct()

write.csv(ks_stat, file= paste0(dir, "/genes_a_ks_statistics_vscontrol_genes_all",all, "_", analysis, ".csv"), row.names=FALSE)

p = ggplot(data=ks_stat_fil, aes(x=statistic, y=log10(padj)*-1)) +
  geom_point(aes(colour = analysis), size=2, alpha = 0.3) +
  ggrepel::geom_text_repel(
    data=ks_stat_fil,
    aes(label = gene),
    size = 4,
    max.overlaps=10
  ) +
  geom_hline(yintercept = log10(0.05)*-1, linetype="dotted", color="red") +
  geom_vline(xintercept = 0, linetype="dotted") +
  ggtitle("KS test statistics  genes_a") +
  theme_minimal() +
  labs(x = "KS statistic summary", y = "-Log10 adj p-value") +
  scale_color_manual(values = viridis::turbo(length(unique(ks_stat_fil$analysis))))

ggsave(p, filename = paste0(dir, "/volcano_ks_stats_vs_control_genes_a_all",all, "_", analysis, ".pdf"),
       width = 14, height = 8)



# loop for gene_i
ks_stat = data.frame()
tmp2 = pt
for (c in unique(pt$gene_i)) {
  tryCatch({
    print(c)
    tmp_TET = filter(tmp2, gene_i == c)
    
    if (all == T) {
      for (CG in control_genes_i) {
        tryCatch({
          print(CG)
          tmp_CTR = filter(tmp2, gene_i == CG)
          
          ks_g = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="greater")
          ks_stat = rbind(ks_stat, data.frame(gene = c, pval = ks_g$p.value, statistic = -ks_g$statistic[["D^+"]], analysis = paste0("vs", CG), type = "greater"))
          
          ks_l = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="less")
          ks_stat = rbind(ks_stat, data.frame(gene = c, pval = ks_l$p.value, statistic = ks_l$statistic[["D^-"]], analysis = paste0("vs", CG), type = "less"))
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }
    }
    
    # control together
    tmp_CTR = filter(tmp2, gene_i %in% control_genes_i)
    ks_g = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="greater")
    ks_l = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="less")
    
    ks_stat = rbind(ks_stat, 
                    data.frame(gene = c, pval = ks_g$p.value, statistic = -ks_g$statistic[["D^+"]], analysis="vsControlgenes_i", type="greater"),
                    data.frame(gene = c, pval = ks_l$p.value, statistic = ks_l$statistic[["D^-"]], analysis="vsControlgenes_i", type="less"))
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# loop for gene_i
for (c in unique(pt$gene_i)) {
  tryCatch({
    print(c)
    tmp_TET = filter(tmp2, gene_i == c)
    
    # control together
    tmp_CTR = filter(tmp2, gene_i %in% control_genes_i)
    ks_g = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="greater")
    ks_l = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="less")
    
    ks_stat = rbind(ks_stat, 
                    data.frame(gene = c, pval = ks_g$p.value, statistic = -ks_g$statistic[["D^+"]], analysis="vsControlgenes_i", type="greater"),
                    data.frame(gene = c, pval = ks_l$p.value, statistic = ks_l$statistic[["D^-"]], analysis="vsControlgenes_i", type="less"))
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

ks_stat$padj <- p.adjust(ks_stat$pval, method = "BH")
ks_stat$sig <- ks_stat$padj <= 0.05

ks_stat_fil = ks_stat %>%
  group_by(gene, analysis) %>%
  filter(padj == min(padj)) %>%
  ungroup() %>% 
  distinct()

write.csv(ks_stat, file= paste0(dir, "/genes_i_ks_statistics_vscontrol_genes_all",all, "_", analysis, ".csv"), row.names=FALSE)

p = ggplot(data=ks_stat_fil, aes(x=statistic, y=log10(padj)*-1)) +
  geom_point(aes(colour = analysis), size=2, alpha = 0.3) +
  ggrepel::geom_text_repel(
    data=ks_stat_fil,
    aes(label = gene),
    size = 4,
    max.overlaps=10
  ) +
  geom_hline(yintercept = log10(0.05)*-1, linetype="dotted", color="red") +
  geom_vline(xintercept = 0, linetype="dotted") +
  ggtitle("KS test statistics  genes_i") +
  theme_minimal() +
  labs(x = "KS statistic summary", y = "-Log10 adj p-value") +
  scale_color_manual(values = viridis::turbo(length(unique(ks_stat_fil$analysis))))

ggsave(p, filename = paste0(dir, "/volcano_ks_stats_vs_control_genes_i_all",all, "_", analysis, ".pdf"),
       width = 14, height = 8)



# # --- KS TEST based on guide_a ---
# ks_stat = data.frame()
# tmp2 = pt
# for (c in unique(pt$guide_a)) {
#   tryCatch({
#     print(c)
#     tmp_TET = filter(tmp2, guide_a == c)
#     
#     if (all == T) {
#       for (CG in control_genes_a) {
#         tryCatch({
#           print(CG)
#           tmp_CTR = filter(tmp2, guide_a == CG)
#           
#           ks_g = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="greater")
#           ks_stat = rbind(ks_stat, data.frame(gene = c, pval = ks_g$p.value, statistic = -ks_g$statistic[["D^+"]], analysis = paste0("vs", CG), type = "greater"))
#           
#           ks_l = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="less")
#           ks_stat = rbind(ks_stat, data.frame(gene = c, pval = ks_l$p.value, statistic = ks_l$statistic[["D^-"]], analysis = paste0("vs", CG), type = "less"))
#         }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
#       }
#     }
#     
#     # control together
#     tmp_CTR = filter(tmp2, guide_a %in% control_genes_a)
#     ks_g = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="greater")
#     ks_l = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="less")
#     
#     ks_stat = rbind(ks_stat, 
#                     data.frame(gene = c, pval = ks_g$p.value, statistic = -ks_g$statistic[["D^+"]], analysis="vsControlgenes_a", type="greater"),
#                     data.frame(gene = c, pval = ks_l$p.value, statistic = ks_l$statistic[["D^-"]], analysis="vsControlgenes_a", type="less"))
#     
#   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
# }
# 
# # loop for guide_a
# for (c in unique(pt$guide_a)) {
#   tryCatch({
#     print(c)
#     tmp_TET = filter(tmp2, guide_a == c)
#     
#     # control together
#     tmp_CTR = filter(tmp2, guide_a %in% control_genes_a)
#     ks_g = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="greater")
#     ks_l = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="less")
#     
#     ks_stat = rbind(ks_stat, 
#                     data.frame(gene = c, pval = ks_g$p.value, statistic = -ks_g$statistic[["D^+"]], analysis="vsControlgenes_a", type="greater"),
#                     data.frame(gene = c, pval = ks_l$p.value, statistic = ks_l$statistic[["D^-"]], analysis="vsControlgenes_a", type="less"))
#     
#   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
# }
# 
# ks_stat$padj <- p.adjust(ks_stat$pval, method = "BH")
# ks_stat$sig <- ks_stat$padj <= 0.05
# 
# ks_stat_fil = ks_stat %>%
#   group_by(gene, analysis) %>%
#   filter(padj == min(padj)) %>%
#   ungroup() %>% 
#   distinct()
# 
# write.csv(ks_stat, file= paste0(dir, "/genes_i_ks_statistics_vscontrol_guide_all",all, "_", analysis, ".csv"), row.names=FALSE)
# 
# p = ggplot(data=ks_stat_fil, aes(x=statistic, y=log10(padj)*-1)) +
#   geom_point(aes(colour = analysis), size=2, alpha = 0.3) +
#   ggrepel::geom_text_repel(
#     data=ks_stat_fil,
#     aes(label = gene),
#     size = 4,
#     max.overlaps=10
#   ) +
#   geom_hline(yintercept = log10(0.05)*-1, linetype="dotted", color="red") +
#   geom_vline(xintercept = 0, linetype="dotted") +
#   ggtitle("KS test statistics  guide_a") +
#   theme_minimal() +
#   labs(x = "KS statistic summary", y = "-Log10 adj p-value") +
#   scale_color_manual(values = viridis::turbo(length(unique(ks_stat_fil$analysis))))
# 
# ggsave(p, filename = paste0(dir, "/volcano_ks_stats_vs_control_guide_a_all",all, "_", analysis, ".pdf"),
#        width = 14, height = 8)
# 
# 
# # --- KS TEST based on guide_i ---
# ks_stat = data.frame()
# tmp2 = pt
# for (c in unique(pt$guide_i)) {
#   tryCatch({
#     print(c)
#     tmp_TET = filter(tmp2, guide_i == c)
#     
#     if (all == T) {
#       for (CG in control_genes_i) {
#         tryCatch({
#           print(CG)
#           tmp_CTR = filter(tmp2, guide_i == CG)
#           
#           ks_g = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="greater")
#           ks_stat = rbind(ks_stat, data.frame(gene = c, pval = ks_g$p.value, statistic = -ks_g$statistic[["D^+"]], analysis = paste0("vs", CG), type = "greater"))
#           
#           ks_l = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="less")
#           ks_stat = rbind(ks_stat, data.frame(gene = c, pval = ks_l$p.value, statistic = ks_l$statistic[["D^-"]], analysis = paste0("vs", CG), type = "less"))
#         }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
#       }
#     }
#     
#     # control together
#     tmp_CTR = filter(tmp2, guide_i %in% control_genes_i)
#     ks_g = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="greater")
#     ks_l = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="less")
#     
#     ks_stat = rbind(ks_stat, 
#                     data.frame(gene = c, pval = ks_g$p.value, statistic = -ks_g$statistic[["D^+"]], analysis="vsControlgenes_i", type="greater"),
#                     data.frame(gene = c, pval = ks_l$p.value, statistic = ks_l$statistic[["D^-"]], analysis="vsControlgenes_i", type="less"))
#     
#   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
# }
# 
# # loop for guide_i
# for (c in unique(pt$guide_i)) {
#   tryCatch({
#     print(c)
#     tmp_TET = filter(tmp2, guide_i == c)
#     
#     # control together
#     tmp_CTR = filter(tmp2, guide_i %in% control_genes_i)
#     ks_g = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="greater")
#     ks_l = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="less")
#     
#     ks_stat = rbind(ks_stat, 
#                     data.frame(gene = c, pval = ks_g$p.value, statistic = -ks_g$statistic[["D^+"]], analysis="vsControlgenes_i", type="greater"),
#                     data.frame(gene = c, pval = ks_l$p.value, statistic = ks_l$statistic[["D^-"]], analysis="vsControlgenes_i", type="less"))
#     
#   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
# }
# 
# ks_stat$padj <- p.adjust(ks_stat$pval, method = "BH")
# ks_stat$sig <- ks_stat$padj <= 0.05
# 
# ks_stat_fil = ks_stat %>%
#   group_by(gene, analysis) %>%
#   filter(padj == min(padj)) %>%
#   ungroup() %>% 
#   distinct()
# 
# write.csv(ks_stat, file= paste0(dir, "/genes_i_ks_statistics_vscontrol_guide_all",all, "_", analysis, ".csv"), row.names=FALSE)
# 
# ggplot(data=ks_stat_fil, aes(x=statistic, y=log10(padj)*-1)) +
#   geom_point(aes(colour = analysis), size=2, alpha = 0.3) +
#   ggrepel::geom_text_repel(
#     data=ks_stat_fil,
#     aes(label = gene),
#     size = 4,
#     max.overlaps=10
#   ) +
#   geom_hline(yintercept = log10(0.05)*-1, linetype="dotted", color="red") +
#   geom_vline(xintercept = 0, linetype="dotted") +
#   ggtitle("KS test statistics  guide_i") +
#   theme_minimal() +
#   labs(x = "KS statistic summary", y = "-Log10 adj p-value") +
#   scale_color_manual(values = viridis::turbo(length(unique(ks_stat_fil$analysis))))
# 
# ggsave(filename = paste0(dir, "/volcano_ks_stats_vs_control_guide_i_all",all, "_", analysis, ".pdf"),
#        width = 14, height = 8)



#gene and pseudotime correlation
cds <- cds[
  ,
  colData(cds) %>%
    row.names
]
pr_graph_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=1)
pr_graph_test_res = pr_graph_test_res[order(pr_graph_test_res$morans_I, decreasing = T), ]
write.csv(pr_graph_test_res, file= paste0(dir,"regression_pseudotime.csv"))
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

p = plot_cells(cds, genes=pr_deg_ids[1:30],
               show_trajectory_graph=FALSE,
               label_cell_groups=FALSE,
               label_leaves=FALSE)
ggsave(p, filename = paste0(dir, "/correlated_pseudotime_gene_exp", analysis, ".pdf"),
       width = 14, height = 14) 


###############################################
# --- KS-TEST: CIRI guide_comb vs CRISPRa NTCa (NTCa-NA) ---
###############################################

# Define the specific control target
control_target_comb = "NTCa_1A;NTCa_1B-NA;NA"

# 1. Extract control cells from the GLOBAL pt object 
# (We use 'pt' instead of 'pt_CIRI' to ensure we capture the NTCa-NA cells 
# even if they are labeled with a different 'type' like 'CRISPRa')
tmp_ctrl_ntca = filter(pt, guide_comb == control_target_comb)

message("\n--- Starting analysis: CIRI vs ", control_target_comb, " ---")
message("Number of control cells found (", control_target_comb, "): ", nrow(tmp_ctrl_ntca))

# Check if control exists
if (nrow(tmp_ctrl_ntca) < 5) {
  warning("WARNING: Not enough control cells found for ", control_target_comb, ". Skipping this comparison.")
} else {
  
  ks_df_vs_ntca = data.frame()
  
  # 2. Iterate through the filtered CIRI guide combinations
  # (Using pt_CIRI which is already filtered for n >= 8 and type == CIRI)
  for (gc in unique(pt_CIRI$guide_comb)) {
    
    # Skip if the current group IS the control (sanity check)
    if (gc == control_target_comb) next
    
    tmp_case = filter(pt_CIRI, guide_comb == gc)
    
    # Quality check: ensure enough cells and variance
    if (nrow(tmp_case) < 2 || length(unique(tmp_case$pseudotime)) <= 1) next
    
    # Run KS Test (Greater)
    ks_g = ks.test(tmp_case$pseudotime, tmp_ctrl_ntca$pseudotime, alternative = "greater")
    # Run KS Test (Less)
    ks_l = ks.test(tmp_case$pseudotime, tmp_ctrl_ntca$pseudotime, alternative = "less")
    
    # Store results
    ks_df_vs_ntca = rbind(
      ks_df_vs_ntca,
      data.frame(
        guide_comb = gc,
        pval = ks_g$p.value,
        stat = -ks_g$statistic[["D^+"]],
        direction = "greater"
      ),
      data.frame(
        guide_comb = gc,
        pval = ks_l$p.value,
        stat = ks_l$statistic[["D^-"]],
        direction = "less"
      )
    )
  }
  
  # 3. Save and Plot
  if (nrow(ks_df_vs_ntca) > 0) {
    # Adjust p-values
    ks_df_vs_ntca$padj = p.adjust(ks_df_vs_ntca$pval, method = "BH")
    ks_df_vs_ntca$neglog = -log10(ks_df_vs_ntca$padj)
    ks_df_vs_ntca$sig = ks_df_vs_ntca$padj <= 0.05
    
    # Output CSV
    out_csv_ntca = paste0(dir, "/CIRI_guide_comb_vs_", control_target_comb, "_ks_statistics.csv")
    write.csv(ks_df_vs_ntca, out_csv_ntca, row.names = FALSE)
    message("Saved CSV: ", out_csv_ntca)
    
    ks_df_vs_ntca = separate(ks_df_vs_ntca, guide_comb, into = c("guide_a", "guide_i"), sep = "-", remove = F)
    ks_df_vs_ntca = separate(ks_df_vs_ntca, guide_a, into = c("guide_a1", "guide_a2"), sep = ";", remove = F)
    ks_df_vs_ntca = separate(ks_df_vs_ntca, guide_a1, into = c("gene_a"), sep = "_", remove = F)
    
    # Volcano Plot
    p_volc_ntca = ggplot(ks_df_vs_ntca, aes(x = stat, y = neglog)) +
      geom_point(aes(colour = gene_a), alpha = 0.5, size = 2) +
      ggrepel::geom_text_repel(aes(label = guide_comb), size = 3, max.overlaps = 20) +
      geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "red") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      theme_minimal(base_size = 14) +
      labs(
        title = paste("KS Volcano: CIRI guide_comb vs", control_target_comb),
        x = "KS statistic",
        y = "-log10 adj p-value"
      ) +
      scale_color_viridis_d()
    
    out_pdf_ntca = paste0(dir, "/CIRI_volcano_guide_comb_vs_", control_target_comb, ".pdf")
    ggsave(filename = out_pdf_ntca, plot = p_volc_ntca, width = 12, height = 7)
    message("Saved Plot: ", out_pdf_ntca)
  } else {
    message("No significant results or no comparisons run for NTCa analysis.")
  }
}



