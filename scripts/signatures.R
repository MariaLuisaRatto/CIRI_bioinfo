library(monocle3)
library(dplyr)
library(ggplot2)


dir = "."
ccs = c("4")
group_name = "muscle"
group_of_interest = paste("Group", paste(ccs, collapse = "_"), sep = "_")
load(paste0(dir, "/processed_cds_", group_name, ".RData"))

# ============================================================
# HUMAN SKELETAL MUSCLE & iPSC-MYOGENIC SIGNATURES (v2025)
# Integrated from Chemello et al. 2025, De Micheli et al. 2020,
# Murgia et al. 2021, Hicks et al. 2025, Front Cell Dev Biol 2025,
# and classic literature (Schiaffino & Reggiani 2011; Bodine 2001).
# ============================================================

gs <- list(
  
  # ===============================
  # CELL CYCLE / RIBOSOME TONE
  # ===============================
  cellcycle_s = c("MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2",
                  "MCM6","CDCA7","DTL","PRIM1","UHRF1","HELLS","RFC2","RPA2","NASP",
                  "RAD51","CHEK1","ORC6","MCM3","MCM7","MLF1IP","BRIP1","E2F8"),
  # Ref: Seurat S.Score; Whitfield et al., Cell 2002; Trapnell et al., Nat Biotech 2014
  
  cellcycle_g2m = c("HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80",
                    "CKS2","NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","SMC4",
                    "CCNB2","CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B",
                    "GTSE1","KIF20B","HJURP","CDCA3","CDC20","TTK","CDC45","CDC6",
                    "EXO1","TIPIN","DSCC1","BLM","CASP8AP2","USP1","CLSPN","POLA1",
                    "CHAF1B","BRCA1","EXOSC2","USP16","POLD3","MSH2","ATAD2","RAD51AP1","RRM2"),
  # Ref: Seurat G2M.Score; Whitfield et al., Cell 2002; Trapnell 2014
  
  ribosome_tone = c("RPLP0","RPL3","RPL5","RPL7","RPL11","RPL13A","RPL23A",
                    "RPS3","RPS6","RPS8","RPS12","RPS18","RPS27","NCL","FBL",
                    "NPM1","EIF4EBP1","MYC","DDIT4"),
  # Ref: Murgia et al. 2021; Chemello et al. 2025 (fiber proteome core)
  
  # ===============================
  # FIBER-TYPE / METABOLIC PROGRAMS
  # ===============================
  fiber_type_I_slow = c("MYH7","TNNI1","TNNT1","MYL2","ATP2A2","CKMT2","MB",
                        "COX7A1","COX6A2","PPARGC1A","SDHA","CS","SLC25A4",
                        "MYBPC1","MYL3","PDK4","ESRRA","NDUFA10"),
  # Ref: Schiaffino & Reggiani 2011; Chemello et al. 2025; Murgia et al. 2021
  
  fiber_type_IIa_fast = c("MYH2","TNNI2","TNNT3","MYL1","ATP2A1","ACTN3","TPM2",
                          "PVALB","LDHA","PKM","PFKM","ENO1","PGAM2","ALDOA",
                          "CKM","DMD","CASQ1"),
  # Ref: Murgia et al. 2021; Bloemberg & Quadrilatero 2012; Chemello 2025
  
  fiber_type_IIx_fast = c("MYH1","ACTN3","MYLK2","TNNT3","TNNI2","ATP2A1","TPM2",
                          "CKM","LDHA","ENO1","PFKM","MYH4"),
  # Ref: Murgia 2021; Uhlén et al., HPA 2021
  
  developmental_embryonic_perinatal = c("MYH3","MYH8","MYL4","TNNI1","TNNT1",
                                        "MYOG","MYOD1","MEF2C","DES","MYMK",
                                        "MYMX","CHRNG","ACTA1","CKM"),
  # Ref: Chal et al. 2015; Hicks et al. 2025; Sampath et al. 2018
  
  sarcomere_core = c("TTN","NEB","ACTN2","TCAP","MYBPC1","MYBPC2","MYL1","MYL2",
                     "TNNI1","TNNI2","TNNT1","TNNT3","TPM1","TPM2","MYOM1","MYOM2","CSRP3"),
  # Ref: Lundberg et al. 2019; GO:0030016 "myofibril"
  
  calcium_excitation_contraction = c("RYR1","CACNA1S","ATP2A1","ATP2A2","CASQ1",
                                     "CASQ2","SLN","PLN","TRDN","STAC3","JPH1",
                                     "CAMK2A","CAMK2B","CAV3"),
  # Ref: Protasi et al. 2021; Jorgensen et al. 2023 (hiPSC EC coupling)
  
  oxidative_mitochondrial = c("SDHA","SDHB","UQCRC1","NDUFA9","NDUFS1","COX4I1",
                              "COX6A2","ATP5F1A","ATP5MC1","CKMT2","SLC25A4",
                              "PPARGC1A","TFAM","CS","ACADM","NDUFB8","CYCS","PDHA1"),
  # Ref: Lin et al. Nature 2002; Chemello et al. 2025; Murgia 2021
  
  glycolytic_fast = c("HK2","PFKM","ALDOA","TPI1","GAPDH","PGK1","PGAM2","ENO1",
                      "PKM","LDHA","PDK1","SLC2A4","PFKFB3","ENO3","MYL1"),
  # Ref: Schiaffino & Reggiani 2011; Chemello 2025
  
  # ===============================
  # STEM / MYOGENIC LINEAGES
  # ===============================
  satellite_cell_core = c("PAX7","VCAM1","ITGA7","CXCR4","MYF5","CDH15","SDC4",
                          "EGFR","NCAM1","SPRY1","MUSTN1","NOS1"),
  # Ref: De Micheli et al. 2020; Yin et al. 2013; GSE143437 atlas
  
  myoblast_proliferation_cycle = c("MYOD1","MYF5","MKI67","TOP2A","PCNA","CCNB1",
                                   "CDK1","E2F1","NUSAP1","HMGB2","UBE2C","AURKB","RRM2"),
  # Ref: Hicks et al. 2025; Trapnell et al. 2014 pseudotime
  
  myotube_differentiation_fusion = c("MYOG","MEF2C","MYMK","MYMX","MYH3","MYH8",
                                     "CKM","DES","DNM2","ACTA1"),
  # Ref: Chal et al. 2015; Sampath et al. 2018; Hicks 2025
  
  myogenic_regulators_axis = c("PAX3","PAX7","MYF5","MYOD1","MYOG","MEF2C","MYF6"),
  # Ref: Weintraub 1991; Rudnicki et al. 1993
  
  # ===============================
  # NMJ / ECM / NICHE
  # ===============================
  neuromuscular_junction = c("AGRN","LRP4","MUSK","DOK7","RAPSN","COLQ","CHRNA1",
                             "CHRNB1","CHRND","CHRNE","CHRNG","ACHE","UTRN"),
  # Ref: Burden et al. 2018; Cantor et al. 2022 (hiPSC NMJ coculture)
  
  muscle_ecm_basement = c("LAMA2","LAMC1","COL6A1","COL6A2","COL6A3","COL1A1",
                          "COL1A2","COL3A1","FN1","SPP1","LUM","DCN","DAG1",
                          "ITGA7","THY1","COL4A1","COL4A2","COL5A1"),
  # Ref: De Micheli et al. 2020; Saini et al. 2021
  
  fibro_adipogenic_progenitors = c("PDGFRA","COL1A1","COL3A1","DCN","LUM","CXCL14",
                                   "THY1","GPC6","FBN1","DPT"),
  # Ref: De Micheli 2020; Front Cell Dev Biol 2025
  
  atrophy_stress = c("FBXO32","TRIM63","DDIT4","GADD45A","ATF4","FOXO3","PSMA1",
                     "PSMB1","SESN1","SESN3","EIF4EBP1","PPP1R15A"),
  # Ref: Bodine et al. 2001; Gomes et al. 2001; Milan et al. 2021
  
  # ===============================
  # PLURIPOTENCY → MYOGENIC (hiPSC/hESC)
  # ===============================
  pluripotency_core = c("POU5F1","NANOG","SOX2","PRDM14","DPPA2","DPPA4","ZFP42",
                        "DNMT3A","DNMT3B"),
  # Ref: Theunissen et al. 2016; Nichols & Smith 2012
  
  naive_pluripotency = c("KLF17","KLF4","TFCP2L1","DNMT3L","DPPA3","TBX3","ESRRB","SUSD2"),
  # Ref: Collier et al. 2022; Pastor et al. 2018
  
  mesoderm_specification = c("MESP1","PDGFRA","KDR","ISL1","HAND1","HAND2","GATA4",
                             "TBX5","NKX2-5","MEF2C","FGF8","EOMES","TBXT"),
  # Ref: Chan et al. 2013; Loh et al. 2016; Hicks 2025
  
  early_myogenic_commitment = c("PAX3","PAX7","MYF5","MYOD1","MET","ITGA7","FGFR4","LBX1","MYOG"),
  # Ref: Chal et al. 2015; Hicks et al. 2025
  # --- Pluripotency axes ---
  core_pluripotency = c("POU5F1","NANOG","SOX2","PRDM14","DPPA2","DPPA4","DNMT3A","DNMT3B","ZFP42"),
  naive_pluripotency = c("KLF17","KLF4","TFCP2L1","DNMT3L","DPPA3","TBX3","ESRRB","SUSD2"),
  formative_pluripotency = c("OTX2","DPPA2","DPPA4","POU3F1","ZIC2","ZIC3","SOX3"),  # OCT6 = POU3F1
  primed_pluripotency = c("OTX2","FGF5","ZIC1","ZIC3","ZIC5","SOX11","DNMT3A","DNMT3B","CD24","EPCAM"),
  # --- Epithelial vs EMT/stress tone ---
  epithelial = c("CDH1","EPCAM","CLDN6","CLDN3","TJP1","ITGA6","ITGB1"),
  emt_stress = c("VIM","CD44","ITGA5","ITGAV","SNAI1","SNAI2","ZEB1","ZEB2","TWIST1","MALAT1","FN1"),
  # --- Pathways / signaling readouts ---
  activin_nodal = c("NODAL","LEFTY1","LEFTY2","CER1","SMAD7","TDGF1"),
  bmp_axis = c("BMP4","ID1","ID2","ID3","MSX1","DLX5"),
  fgf_erk = c("FGF2","ETV4","ETV5","SPRY2","DUSP6"),
  wnt_beta_catenin = c("WNT3","WNT3A","AXIN2","LEF1","TCF7","DKK1"),
  notch_axis = c("HES1","HES5"),
  hippo_yap = c("CTGF","CYR61","AMOTL2"),
  # --- Lineage priming / early fate modules ---
  neuroectoderm = c("SOX1","PAX6","LHX2","HES1","HES5","OTX1","NES"),
  neural_border_crest_placodal = c("TFAP2A","TFAP2C","PAX3","PAX7","SOX9","SOX10","DLX5","DLX6","MSX1","GATA2","GATA3","SIX1","EYA1","DACH1","DACH2","KRT8","KRT18","KRT19"),
  mesendoderm = c("TBXT","EOMES","MIXL1","GSC","NODAL","FGF8"),
  cardiac_mesoderm_early = c("MESP1","PDGFRA","KDR","ISL1","HAND1","HAND2","GATA4","NKX2-5","TBX5","MEF2C"),
  endoderm = c("SOX17","FOXA2","GATA6","HHEX","CER1"),
  extraembryonic_amnion_te = c("TFAP2C","GATA2","GATA3","KRT7","PEG10","TACSTD2","CGA","CGB","ISL1","KRT8","KRT18","KRT19"),
  ribosome_tone = c("RPLP0","RPL3","RPL5","RPL7","RPL11","RPL13A","RPL23A","RPS3","RPS6","RPS8","RPS12","RPS18","RPS27","NCL","FBL","NPM1","EIF4EBP1","MYC","DDIT4"),
  # --- Stress responses ---
  oxidative_stress = c("SOD2","PRDX1","PRDX2","TXN"),
  er_upr = c("HSPA5","ATF4","DDIT3","XBP1","HSP90B1"),
  dna_damage_p53 = c("TP53","CDKN1A","GADD45A"),
  apoptosis_axis = c("BAX","BCL2","BCL2L11","PMAIP1"),
  # --- Adhesion/surface (handy for sort validation) ---
  adhesion_surface = c("PODXL","B3GALT5","EPCAM","SUSD2","CD24","ITGA6","ITGB1")
)


# ============================
# FILTER SIGNATURES
# ============================
available = rownames(cds_sub)

gs_filtered = lapply(gs, function(v) intersect(v, available))

empty = names(gs_filtered)[sapply(gs_filtered, length) == 0]
if(length(empty) > 0) {
  message("Empty signatures: ", paste(empty, collapse = ", "))
}

# ============================
# CONVERT TO gene_group_df
# ============================

gene_group_df = do.call(
  rbind,
  lapply(names(gs_filtered), function(sig) {
    data.frame(
      gene_id = gs_filtered[[sig]],
      gene_group = sig,
      stringsAsFactors = FALSE
    )
  })
)

# ============================
# AGGREGATE
# ============================
# Compute aggregated signatures (signature x cell)
agg = aggregate_gene_expression(
  cds = cds_sub,
  gene_group_df = gene_group_df
)

# Transpose to cell x signature
agg_t = t(agg)
agg_t = agg_t[colnames(cds_sub), , drop = FALSE]

sig_df = as.data.frame(agg_t)

# Add signatures to colData
for(sig in colnames(sig_df)){
  colData(cds_sub)[[sig]] <- sig_df[[sig]]
}

# Create output folder if needed
outdir = "signature_umaps"
if(!dir.exists(outdir)) dir.create(outdir)

# # Save one PNG per signature
# for(sig in colnames(sig_df)){
#   p = plot_cells(
#     cds_sub,
#     color_cells_by = sig,
#     show_trajectory_graph = FALSE
#   ) + ggtitle(sig)
#   
#   ggsave(
#     filename = file.path(outdir, paste0(sig, ".pdf")),
#     plot = p,
#     width = 6,
#     height = 5
#   )
# }


#with ggplot 
umap = reducedDims(cds_sub)$UMAP
umap_df = data.frame(
  cell = rownames(umap),
  UMAP_1 = umap[,1],
  UMAP_2 = umap[,2]
)

agg = aggregate_gene_expression(cds_sub, gene_group_df)
agg_t = t(agg)
agg_t = agg_t[colnames(cds_sub), , drop = FALSE]
sig_df = as.data.frame(agg_t)
sig_df$cell = rownames(sig_df)
plot_df = merge(umap_df, sig_df, by = "cell")

for(sig in colnames(sig_df)[colnames(sig_df) != "cell"]){
  
  p = ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = .data[[sig]])) +
    geom_point(size = 0.5) +
    theme_classic() +
    labs(color = sig, title = sig)+
    ggplot2::scale_color_viridis_c(option = "E")
  
  ggsave(
    filename = file.path(outdir, paste0(sig, ".pdf")),
    plot = p,
    width = 6,
    height = 5
  )
}

csv_outdir = "signature_values_csv"
if(!dir.exists(csv_outdir)) dir.create(csv_outdir)

for(sig in colnames(sig_df)[colnames(sig_df) != "cell"]){
  
  df_sig = data.frame(value = plot_df[[sig]])
  rownames(df_sig) = plot_df$cell
  colnames(df_sig) = sig
  
  write.csv(
    df_sig,
    file = file.path(csv_outdir, paste0(sig, ".csv")),
    row.names = TRUE
  )
}




