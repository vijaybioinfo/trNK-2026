
#################################################################
#################################################################
##### Figures Conpendium for KoMa TrNK cells in Lung cancer #####
#################################################################
#################################################################


resources = c(
  "/home/ciro/scripts/handy_functions/devel/file_reading.R", # readfile
  "/home/ciro/scripts/handy_functions/devel/utilities.R",
  # filters_columns is.file.finished show_commas (cluster_reports)
  "/home/ciro/scripts/clustering/R/plotting.R", # cluster_reports
  "/home/ciro/scripts/clustering/R/utilities.R", # get_top_n
  "/home/ciro/scripts/handy_functions/devel/filters.R", # sample_even
  "/home/ciro/scripts/handy_functions/devel/plots.R", # plot_pct getlegend mytheme
  "/home/ciro/scripts/handy_functions/R/stats_summary_table.R",
  "/home/fcastaneda/bin/clustering/R/plotting.R" #plots of clustering pipeline
)
for(i in resources){ source(i) }
source("/home/ciro/scripts/figease/figease.R")
#source("/home/ciro/scripts/figease/figease.R")
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)
library(Seurat)
library(parallel)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(ggridges)
library(matrixStats)
library(ggforce)
library(cowplot)
library(patchwork)
library(tidyverse)
library(phenoptr)
library(scales)


source("/home/donaldo/bioadhoc/spatial/KoMa/tumor_merge/scripts/KoMa_NK_paper_subclustering/figures/global_v3.R")

setwd("/mnt/cephfs/vd-vijay/donaldo/spatial/KoMa/tumor_merge/results/figures/KoMa_NKpaper_figures_v3")


{ cat(redb("### Metadata and matrix for co-clustering Single-cell/Spatial ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

  km_rds_spa <- readRDS("/home/donaldo/bioadhoc/spatial/KoMa/tumor_merge/results/clutering/KoMa_NK_paper/KoMa_trNKpaper_run2_slide1_samples_selected_REDO/.object_init_mean0.01_pct40_pc25.rds")
  
  km_rds_spa_sub <- subset(km_rds_spa, subset = Xenium_snn_res.0.5 == '2')
  km_rds_spa_sub@meta.data$orig.data_class <- "Xenium"
  km_rds_spa_sub@meta.data$orig.celltype <- paste0("X_", km_rds_spa_sub@meta.data$Xenium_snn_res.0.5)
  km_rds_spa_sub@meta.data$orig.sub_celltype <- paste0("X_", km_rds_spa_sub@meta.data$Xenium_snn_res.0.5)
  table(km_rds_spa_sub@meta.data$orig.celltype, km_rds_spa_sub@meta.data$Xenium_snn_res.0.5)

  identical(rownames(km_rds_spa_sub@meta.data), colnames(km_rds_spa_sub))

  #### ------------ Single-Cell object (only immune clusters)
  #rds selected TNK_3025
  source("/home/donaldo/bioadhoc/KoMa/scripts/figures/global.R")
  setwd("/home/donaldo/bioadhoc/spatial/KoMa/tumor_merge/results/figures/KoMa_NKpaper_figures_v3")
  aaaa <- "/mnt/cephfs/vd-vijay/donaldo/spatial/KoMa/tumor_merge/results/figures/KoMa_NKpaper_figures_v3"
  if(getwd() != aaaa) {stop(cat("working dir is wrong!", "\n"))} else {cat('working directory is good!', "\n")}

  TNK_3025$orig.data_class <- "single_cell"
  TNK_3025$origlib <- TNK_3025$orig.donors.sort.batch.tag
  md_keep_sc <- TNK_3025@meta.data[, c(4, 7, 18, 19, 31:35, 45:46)]
  TNK_KM_sc <- TNK_3025
  TNK_KM_sc@meta.data <- md_keep_sc
  identical(rownames(TNK_KM_sc@meta.data), colnames(TNK_KM_sc))

  ##### ---------------------------- Integration ---------------------------- #####

  # To accomplish data transfer across cells of either two assays (scRNA-seq and imaging-based ST), we first need to find cells that are similar between the two datasets. 
  # -> To this end, we use the common genes of both assays to integrate both datasets into a joint embedding space.
  # -> to find the nearest scRNA-seq cell neighbors of each Xenium cell in the joint embedding space, assuming that these cells will have similar transcriptional profiles.


  # Integration of sc and Xenium object
  # genes present in both datasets
  gs_TNK <- intersect(rownames(km_rds_spa_sub), rownames(TNK_KM_sc))
  # 3329 genes

  obj_TNK_scxn <- lapply(list(km_rds_spa_sub, TNK_KM_sc), function(obj_TNK_scxn) {
      # obj_TNK_scxn <- km_rds_spa_sub
      obj_TNK_scxn <- subset(obj_TNK_scxn, features = gs_TNK) #keep intersect markers only
      if ("Xenium" %in% unique(obj_TNK_scxn$orig.data_class)) {
      y_SH <- GetAssayData(obj_TNK_scxn, assay = "Xenium", slot = "counts")
      } else {
      y_SH <- GetAssayData(obj_TNK_scxn, slot = "counts") } # sparse count matrix
      ob <- CreateSeuratObject(
        counts = y_SH,
        meta.data = obj_TNK_scxn@meta.data) #new object with only intersect genes matrix
      return(ob)
  })

  obj_TNK_scxn <- merge(obj_TNK_scxn[[1]], obj_TNK_scxn[[2]])
  #joining counts layers
  obj_TNK_scxn <- JoinLayers(obj_TNK_scxn)
  table(obj_TNK_scxn$orig.celltype, obj_TNK_scxn$orig.sub_celltype, useNA='always')
  table(obj_TNK_scxn$orig.data_class)

  aver_counts <- GetAssayData(obj_TNK_scxn, assay = "RNA", slot = "counts")
  aver_counts2 <- as.matrix(aver_counts)
  obj_TNK_scxn$percent.mt <- 10
  table(obj_TNK_scxn$orig.data_class)

  obj_TNK_scxn@meta.data <- obj_TNK_scxn@meta.data[,c(1:14,27:37)]

  # saveRDS(obj_TNK_scxn@meta.data, "/home/donaldo/bioadhoc/spatial/KoMa/tumor_merge/sc_spatial_aggr/paper/info/seurat_objects/concat_seurat_object_metadata_matching_genes_populations_July31_2026_NK.rds")
  # saveRDS(aver_counts2, "/home/donaldo/bioadhoc/spatial/KoMa/tumor_merge/sc_spatial_aggr/paper/info/seurat_objects/concat_counts_mtx_matching_genes_populations_July31_2026_Tcells_NK.rds")

}

{ cat(redb("### SC/XN co-clustering annotations for T/NKpopulations and transfering to general clustering XN ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  
  # KM_tnk_v3 - res 1.6
  # KM_tnk_v3_c10 - 
  KM_tnk_v3 <- readRDS("/mnt/BioScratch/donaldo/KoMA_donors_spe_SC_XN_match_markers_populations_Tcells_NK_v3_run2_object_lock_mean0.01_pct20_pc20.rds"  )
  KM_tnk_v3$RNA_nn <- NULL

  counts_TNK <- GetAssayData(KM_tnk_v3, assay = "RNA", slot = "counts")
  counts_TNK2 <- as.matrix(counts_TNK)

  KM_tnk_v3@meta.data$orig.sub_celltype <- ifelse(KM_tnk_v3@meta.data$orig.sub_celltype == 'CD4_TFH_CTL', 'CD4_TFH', ifelse(KM_tnk_v3@meta.data$orig.sub_celltype == 'CD8_CTL_GZMKp', 'CD8_GZMKp', KM_tnk_v3@meta.data$orig.sub_celltype))

  # saveRDS(counts_TNK2, "/home/donaldo/bioadhoc/spatial/KoMa/tumor_merge/results/clutering/KoMa_NK_paper/KoMA_donors_spe_SC_XN_match_markers_populations_Tcells_NK_v3_run2/seurat_mean0.01_pct20_pc20_res1.6/concat_gene_exp_matrix_20pct20pc_16.rds")
  # saveRDS(KM_tnk_v3@meta.data, "/home/donaldo/bioadhoc/spatial/KoMa/tumor_merge/results/clutering/KoMa_NK_paper/KoMA_donors_spe_SC_XN_match_markers_populations_Tcells_NK_v3_run2/seurat_mean0.01_pct20_pc20_res1.6/metadata_20pct20pc_16.rds")

  # 0,17,18 -> CD8_GZMKp 
  # 1,7 CD8_TRM_TCR_eng
  # 2 -> CD4_TCM, CD8_TCM, CD4_TRM
  # 19 -> CD4_TRM
  # 3 -> CD8_IFNR
  # 4,13,15 -> CD8_TRM
  # 4 -> CD4_TFH
  # 6,8,10 -> CD4_TFR_CD4_THIFNR_CD4_Tregs_CD8_IFNR
  # 9 -> CD8_TRM_TCR_eng_CD4_TFH
  # 11 -> CD8_TEM
  # 12 -> NK CD16n (to do tNK tr and circ sub-clustering)
  # 14 -> NK CD16p
  # 16 -> CD8_KIRp

  col_annot_res1.6  <- c ("0" = 'CD8_GZMKp', #
                         "1" = 'CD8_TRM_TCR_eng', #
                         "2" = 'CD4_TCM_CD8_TCM_CD4_TRM', #
                         "3" = 'CD8_IFNR', #
                         "4" = 'CD8_TRM', #
                         "5" = 'CD4_TFH', #
                         "6" = 'CD4_TFR_CD4_THIFNR_CD4_Tregs_CD8_IFNR', #
                         "7" = 'CD8_TRM_TCR_eng', #
                         "8" = 'CD4_TFR_CD4_THIFNR_CD4_Tregs_CD8_IFNR', #
                         "9" = 'CD8_TRM_TCR_eng_CD4_TFH', #
                         "10" = 'CD4_TFR_CD4_THIFNR_CD4_Tregs_CD8_IFNR', #
                         "11" = 'CD8_TEM', #
                         "12" = 'NK_CD16n', #
                         "13" = 'CD8_TRM', #
                         "14" = 'NK_CD16n', #
                         "15" = 'CD8_TRM', #
                         "16" = 'CD8_KIRp', #
                         "17" = 'CD8_GZMKp', #
                         "18" = 'CD8_GZMKp', #
                         "19" = 'CD4_TRM' #
                         )

  KM_tnk_v3@meta.data$orig.sub_celltype_coclust <- col_annot_res1.6[as.character(KM_tnk_v3$RNA_snn_res.1.6)]
  table(KM_tnk_v3@meta.data$orig.sub_celltype_coclust, KM_tnk_v3$RNA_snn_res.1.6, useNA='always')

  ##### correction in annotations for clusters

  # c12 (NK_CD16n)
  # KM_tnk_v3_c12 - res 0.3
  KM_tnk_v3_c12$orig.celltype_clean <- ifelse(KM_tnk_v3_c12$RNA_snn_res.0.3 == '0', 'NK_CD16n_residency_markers', 'NK_CD16n_circulation_markers')
  table(KM_tnk_v3_c12$orig.celltype_clean, KM_tnk_v3_c12$RNA_snn_res.0.3, KM_tnk_v3_c12$orig.data_class)
  table(KM_tnk_v3_c12$orig.celltype_clean, KM_tnk_v3_c12$orig.sub_celltype)

  # c2 (CD4_TCM, CD8_TCM, CD4_TRM)
  # KM_tnk_v3_c2 - res 0.5
  col_annot_c2  <- c ("0" = 'CD4_TCM',
                         "1" = 'CD8_TCM',
                         "2" = 'CD4_TRM',
                         "3" = 'CD4_TCM')

  KM_tnk_v3_c2@meta.data$orig.celltype_clean <- col_annot_c2[as.character(KM_tnk_v3_c2$RNA_snn_res.0.5)]
  table(KM_tnk_v3_c2$orig.celltype_clean, KM_tnk_v3_c2$RNA_snn_res.0.5, KM_tnk_v3_c2$orig.data_class)

  # c6 & c8 & c10 (CD4_TFR, CD4_THIFNR, CD4_Tregs, CD8_IFNR)
  # KM_tnk_v3_c6_8_10 - res 0.4
  col_annot_c6_8_10  <- c ("0" = 'CD4_Tregs',
                         "1" = 'CD4_TFR',
                         "2" = 'CD4_TFR',
                         "3" = 'CD8_IFNR',
                         "4" = 'CD4_THIFNR',
                         "5" = 'CD4_THIFNR'
                         )

  KM_tnk_v3_c6_8_10@meta.data$orig.celltype_clean <- col_annot_c6_8_10[as.character(KM_tnk_v3_c6_8_10$RNA_snn_res.0.4)]
  table(KM_tnk_v3_c6_8_10$orig.celltype_clean, KM_tnk_v3_c6_8_10$RNA_snn_res.0.4, KM_tnk_v3_c6_8_10$orig.data_class)

  # c9 clean-up (CD8_TRM_TCR_eng, CD4_TFH)
  # KM_tnk_v3_c9 - res 0.2
  KM_tnk_v3_c9$orig.celltype_clean <- ifelse(KM_tnk_v3_c9$RNA_snn_res.0.2 == '0', 'CD8_TRM_TCR_eng', 'CD4_TFH')
  table(KM_tnk_v3_c9$orig.celltype_clean, KM_tnk_v3_c9$RNA_snn_res.0.2, KM_tnk_v3_c9$orig.data_class)

  #transder clean annotations to seurat object with rest of T/NK populations
  seu_cleanups <- list(KM_tnk_v3_c12, KM_tnk_v3_c2, KM_tnk_v3_c6_8_10, KM_tnk_v3_c9)

  annots_to_update <- lapply(seq_along(seu_cleanups), function(x) {
    # x <- 1
    rds_cu <- seu_cleanups[[x]]
    rds_cu$cellname <- rownames(rds_cu@meta.data)
    rds_cu2 <- rds_cu@meta.data[rds_cu$orig.data_class == 'Xenium', c('cellname', 'orig.celltype_clean')]
    return(rds_cu2)
  })

  clean_labels_join <- Reduce(rbind, annots_to_update)
  colnames(clean_labels_join)[2] <- "orig.sub_celltype_coclust_clean"

  #joining corrected annotations fro Tc/NK populations to object
  KM_tnk_v3@meta.data$cellname <- rownames(KM_tnk_v3@meta.data)
  KM_tnk_v3@meta.data <- left_join(KM_tnk_v3@meta.data, clean_labels_join, by = 'cellname')
  rownames(KM_tnk_v3@meta.data) <- KM_tnk_v3@meta.data$cellname

  KM_tnk_v3@meta.data$orig.sub_celltype_coclust_clean2 <- case_when(
    is.na(KM_tnk_v3@meta.data$orig.sub_celltype_coclust_clean) & KM_tnk_v3@meta.data$orig.data_class == 'Xenium' ~ KM_tnk_v3@meta.data$orig.sub_celltype_coclust, 
    is.na(KM_tnk_v3@meta.data$orig.sub_celltype_coclust_clean) & KM_tnk_v3@meta.data$orig.data_class == 'single_cell' ~ KM_tnk_v3@meta.data$orig.sub_celltype,
    TRUE ~ KM_tnk_v3@meta.data$orig.sub_celltype_coclust_clean
    )

    KM_tnk_v3@meta.data = joindf(KM_tnk_v3@meta.data,
    as.data.frame(KM_tnk_v3@reductions$umap@cell.embeddings))

  # saveRDS(KM_tnk_v3, "/home/donaldo/bioadhoc/spatial/KoMa/tumor_merge/info/KoMa_NKpaper/seurat_objects/KoMa_single_cell_only_coclustering_single_TFH_Jul31_2026.rds")


  ###### General clustering annotation ######
  #labeling general populations
  cl_gen <- c("0"='Tumor_cells', #
              "1"='CAF', #
              "2"='Tcells_NK', #
              "3"='Myeloid', #
              "4"='Plasma_cells', #
              "5"='Tumor_cells', #
              "6"='Tumor_cells', #
              "7"='Endothelial_cells', #
              "8"='Plasma_cells', #
              "9"='Tumor_cells', #
              "10"='Tumor_cells', #
              "11"='Endothelial_cells', #
              "12"='Tumor_cells', #
              "13"='Tumor_cells', #
              "14"='Bcells', #
              "15"='proliferating_cells', #
              "16"='Tumor_cells' #
              )

  km_rds_spa@meta.data$orig.celltype <- cl_gen[as.character(km_rds_spa@meta.data$Xenium_snn_res.0.5)]
  table(km_rds_spa@meta.data$Xenium_snn_res.0.5, km_rds_spa@meta.data$orig.celltype, useNA='always')

  km_rds_spa@meta.data <- left_join(km_rds_spa@meta.data, KM_tnk_v3@meta.data[, c('orig.sub_celltype_coclust_clean2', 'cellname')], by = 'cellname')
  rownames(km_rds_spa@meta.data) <- km_rds_spa@meta.data$cellname
  km_rds_spa@meta.data$orig.sub_celltype_coclust_clean2 <- case_when(
    is.na(km_rds_spa@meta.data$orig.sub_celltype_coclust_clean2) & km_rds_spa@meta.data$orig.celltype == 'Tcells_NK' ~ 'tc_bq',
    is.na(km_rds_spa@meta.data$orig.sub_celltype_coclust_clean2) & km_rds_spa@meta.data$orig.celltype != 'Tcells_NK' ~ km_rds_spa@meta.data$orig.celltype, 
    TRUE ~ km_rds_spa@meta.data$orig.sub_celltype_coclust_clean2
  )

  table(km_rds_spa@meta.data$orig.sub_celltype_coclust_clean2, km_rds_spa@meta.data$orig.sub_celltype_coclust_clean2, useNA='alway')
  table(km_rds_spa@meta.data$orig.sub_celltype_coclust_clean2, km_rds_spa@meta.data$orig.celltype, useNA='alway')

  km_rds_spa_clean <- subset(km_rds_spa, subset = orig.sub_celltype_coclust_clean2 != 'tc_bq')

  saveRDS(km_rds_spa, "/home/donaldo/bioadhoc/spatial/KoMa/tumor_merge/info/KoMa_NKpaper/seurat_objects/KoMa_coclustering_general_single_TFH_Jul31_wth_bq_2026.rds")
  saveRDS(km_rds_spa_clean, "/home/donaldo/bioadhoc/spatial/KoMa/tumor_merge/info/KoMa_NKpaper/seurat_objects/KoMa_coclustering_general_single_TFH_Jul31_clean_2026.rds")

}

{ cat(redb("### Secondary global variables ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

  colors_gen <- c('Tumor_cells' = '#5D75AE',
              'Tcells_NK' = '#E75883',
              'Bcells' = '#FFEB72',
              'proliferating_cells' = '#5D95FB',
              'CAF' = '#24aed2',
              'Plasma_cells' = '#92cea6',
              'Myeloid' = '#fbb6ca',
              'Endothelial_cells' = '#d2b6ff',
              'tc_bq' = 'gray90')

  #glasbey I use it to generate color_blind friendly paletters
  # library(reticulate)
  # glasbey <- import("glasbey")
  
  # pal <- glasbey$create_palette( palette_size = 40L, colorblind_safe = TRUE)
  # pal <- pal[c(1,2,4,7,8,11,12,15,16,19,21,25,27,29:34,36:40)]
  # pal <- c(pal, "gray90")
  # pops_tnk_co <- c("CD4_TFH", "CD4_Tregs", "CD4_THIFNR", "CD4_TFR", "CD4_TRM", "CD4_TCM", "CD4_CTL_GZMKp", "CD4_CTL_GZMKn", "CD8_TRM_TCR_eng", "CD8_TCM", "CD8_GZMKp", "CD8_TRM", "CD8_IFNR", "CD8_TEM", "CD8_MAIT", "CD8_KIRp", "NK_CD16n_residency_markers", "NK_CD16p", "NK_CD16n_circulation_markers", "NK_CD16p_mature", "NK_ILC3", "NK_IFN", "NK_CD16p_active", "tc_bq")
  # names(pal) <- pops_tnk_co
  # colors_coclust_pops <- pal

  # colors_coclust_pops - > color blind friendly
  colors_coclust_pops <- c('CD4_TFH' = "#0075ff",
                        'CD4_Tregs' = "#69ffba",
                        'CD4_THIFNR' = "#927155",
                        'CD4_TFR' = "#08d2ff",
                        'CD4_TRM' = "#5daa8e",
                        'CD4_TCM' = "#ffb2a2",
                        'CD4_CTL_GZMKp' = "#00650c",
                        'CD4_CTL_GZMKn' = "#0096a6",
                        'CD8_TRM_TCR_eng' = "#7996ff",
                        'CD8_TCM' = "#00dfc2",
                        'CD8_GZMKp' = "#ca2800",
                        'CD8_TRM' = "#000096",
                        'CD8_IFNR' = "#006986",
                        'CD8_TEM' = "#d792aa",
                        'CD8_MAIT' = "#9a001c",
                        'CD8_KIRp' = "#96698a",
                        'NK_CD16n_residency_markers' = "#ff316d",
                        'NK_CD16p' = "#69c66d",
                        'NK_CD16n_circulation_markers' = "#ffe779",
                        'NK_CD16p_mature' = "#aa9ec6",
                        'NK_ILC3' = "#008671",
                        'NK_IFN' = "#9600ca",
                        'NK_CD16p_active' = "#a2f3ff",
                        'tc_bq' = "gray90")


  col_data_class <- c('Xenium' = '#ce1c45', 'single_cell' = '#d2b6ff')

  col_pop_gen <- c('CD8' = '#fbb6ca', 'CD4' = '#24aed2', 'NK' = '#ffaa49')


    #Visualization
    cells<-data.frame(
      cell_id =km_rds_spa_clean$cellname,
      group = km_rds_spa_clean@meta.data[,"orig.sub_celltype_coclust_clean2"]
    )

    dir.create("xenium_explorer_tables/")

    write.csv(cells, "xenium_explorer_tables/cells_classifications_to_xenium_explorer.csv")

}

{ cat(redb("### dim.reduction UMAP resolution chosen sc and xn object ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  
  dir.create("dim.reduction")
  dir.create("dim.reduction/UMAPs_paper_NK_KM")

  fconfigs = list(
    list(result_id = "dim.reduction/UMAPs_paper_NK_KM/UMAP_reductions_KoMa_sc_xn_tnk_coclustering", 
      edata = GetAssayData(KM_tnk_v3, assay = "RNA", slot = "data"), metadata = "KM_tnk_v3@meta.data",
      axis_x = "UMAP_1", axis_y = "UMAP_2",
      vars2col = c("orig.sub_celltype_coclust_clean2"),
      colors = colors_coclust_pops,
      size = c(6,6),
      redu = "UMAP_1",
      facet = "orig.data_class"
    ),
    list(result_id = "dim.reduction/UMAPs_paper_NK_KM/UMAP_reductions_KoMa_sc_xn_tnk_coclustering_XN_populations_v2", 
      edata = GetAssayData(KM_tnk_v3, assay = "RNA", slot = "data"), metadata = "KM_tnk_v3@meta.data[!KM_tnk_v3$orig.sub_celltype_coclust_clean2 %in% c('CD4_CTL_GZMKn','CD4_CTL_GZMKp', 'CD8_MAIT', 'NK_CD16p_active', 'NK_CD16p_mature', 'NK_IFN', 'NK_ILC3'),]",
      axis_x = "UMAP_1", axis_y = "UMAP_2",
      vars2col = c("orig.sub_celltype_coclust_clean2"),
      colors = colors_coclust_pops,
      size = c(6,6),
      redu = "UMAP_1",
      facet = "orig.data_class"
    ),
    list(result_id = "dim.reduction/UMAPs_paper_NK_KM/UMAP_reductions_KoMa_general_clustering", 
      edata = GetAssayData(km_rds_spa_clean, assay = "Xenium", slot = "data"), metadata = "km_rds_spa_clean@meta.data[km_rds_spa_clean$orig.sub_celltype_coclust_clean2 != 'tc_bq', ]",
      axis_x = "umap_1", axis_y = "umap_2",
      vars2col = c("orig.celltype"),
      colors = colors_gen,
      size = c(6,6),
      redu = "umap_1"
    )
  )
  
  pp_main_dr = fig_plot_base(
    fconfigs[2],
    theme_extra = function(x){
      p <- x + geom_point(size = 0.3, na.rm = TRUE) + labs(x = "Dim 1", y = "Dim 2") + theme_classic() +
      guides(col = guide_legend(ncol = 1, override.aes = list(size = 6), label.position = "left"))
    })



}

{ cat(redb("### KoMa: Dotplot of Expression ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

  source('/mnt/BioAdHoc/Groups/vd-vijay/kmatsumura/R24_NK/paper_development/job_scripts/functions_collection.0.3.R')

  ##### function from KoMa dotplot
  dot.plot_xenium <- function(
      seurat.obj, features,
      assay = NULL,                   # <- 追加：Xenium / RNA / SCT など
      layer = NULL,                   # <- 追加：v5の "counts"/"data"/"scale.data" 相当
      slot = 'data',                  # 互換用（v4はslot、v5はlayer優先）
      do.norm = FALSE, ensembl = FALSE,
      groups.tag, groups.order=NULL, groups.of.int=NULL, groups.label=NULL,
      filter.tags=NULL, groups.to.filter=NULL, keep=TRUE, na.rm=TRUE, feature.thold=NULL,
      this.color.scale=NULL, col.min=NULL, col.max=NULL,
      scale.by='radius', dot.scale=6, size.min=NA, size.max=NA,
      file.name=NULL
  ){
    suppressPackageStartupMessages({
      library(data.table); library(tidyr); library(stringr); library(ggplot2)
    })
    
    # ---- color scale default
    if(is.null(this.color.scale)){
      this.color.scale <- c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
    } else {
      if(length(this.color.scale)==1){
        if(this.color.scale %in% c('hot.and.cold','cold.and.hot','cold.n.hot')) this.color.scale <- NA
      }
    }
    
    # ---- size scale
    scale.func <- switch(
      EXPR=scale.by,
      size=scale_size,
      radius=scale_radius,
      stop("'scale.by' must be either 'size' or 'radius'")
    )
    
    # ---- assay default: Xenium優先
    if (is.null(assay)) {
      assay <- if ("Xenium" %in% Assays(seurat.obj)) "Xenium" else DefaultAssay(seurat.obj)
    }
    DefaultAssay(seurat.obj) <- assay
    
    # ---- features check
    features <- unique(features)
    feats_use <- if (ensembl) names(features) else features
    if(!all(feats_use %in% rownames(seurat.obj))){
      missing <- feats_use[!feats_use %in% rownames(seurat.obj)]
      stop("Features not found in object: ", paste(missing, collapse=", "))
    }
    
    # ---- groups.tag existence
    if(!(groups.tag %in% colnames(seurat.obj@meta.data))){
      stop("groups.tag not defined in meta.data: ", groups.tag)
    }
    
    # ---- filter.tags checks
    if(!is.null(filter.tags)){
      ok <- filter.tags %in% colnames(seurat.obj@meta.data)
      if(!all(ok)) stop("Some filter.tags not in meta.data: ", paste(filter.tags[!ok], collapse=", "))
      
      if(length(filter.tags)>1){
        if(!is.list(groups.to.filter) || length(filter.tags)!=length(groups.to.filter)){
          stop("When multiple tags are provided, groups.to.filter must be a list with same length.")
        }
      } else {
        if(!is.list(groups.to.filter)){
          if(!is.vector(groups.to.filter)) stop('groups.to.filter must be a list or vector.')
          groups.to.filter <- list(groups.to.filter)
        }
      }
      names(groups.to.filter) <- paste0('filter.tag.', seq_along(filter.tags))
    }
    
    # ---- helper: safe matrix fetch (v4/v5)
    .get_expr_dt <- function(obj, feats, assay, slot, layer){
      # v5優先：layer 指定があれば layer を使う
      # layerがNULLなら、slotで取りに行く（両対応）
      if (!is.null(layer)) {
        m <- GetAssayData(obj, assay = assay, layer = layer)
        m <- m[feats, , drop=FALSE]
        return(as.data.table(t(as.matrix(m))))
      } else {
        # slot path
        # FetchDataはslotを受け取るが、Assay5で不安定なことがあるので GetAssayData を優先
        m <- tryCatch(
          GetAssayData(obj, assay = assay, slot = slot),
          error = function(e) NULL
        )
        if (is.null(m)) {
          # 最後の手段
          df <- FetchData(obj, vars = feats, slot = slot)
          return(as.data.table(df))
        } else {
          m <- m[feats, , drop=FALSE]
          return(as.data.table(t(as.matrix(m))))
        }
      }
    }
    
    # ---- fetch expression (cells x features)
    feature.vals <- .get_expr_dt(
      obj = seurat.obj,
      feats = feats_use,
      assay = assay,
      slot = slot,
      layer = layer
    )
    setDT(feature.vals)
    colnames(feature.vals) <- feats_use
    
    # ---- tags
    if(!is.null(filter.tags)){
      tag.dt <- data.table(seurat.obj@meta.data[, c(groups.tag, filter.tags), drop=FALSE])
      colnames(tag.dt) <- c("groups.tag", paste0("filter.tag.", seq_along(filter.tags)))
    } else {
      tag.dt <- data.table(groups.tag = seurat.obj@meta.data[, groups.tag])
    }
    
    # ---- join
    tmp.data <- cbind(tag.dt, feature.vals)
    
    # ---- NA drop
    if(na.rm) tmp.data <- tmp.data[!is.na(groups.tag)]
    
    # ---- groups.of.int
    if(!is.null(groups.of.int)){
      tmp.data[!groups.tag %in% groups.of.int, groups.tag := 'Rest']
      if(!is.null(groups.label)){
        tmp.data[groups.tag %in% groups.of.int, groups.tag := as.character(groups.label)]
        new.lvls <- c(as.character(groups.label), 'Rest')
        tmp.data$groups.tag <- factor(as.character(tmp.data$groups.tag), levels=new.lvls)
      }
    }
    
    # ---- filtering
    if(!is.null(filter.tags)){
      if(na.rm){
        for(idx in seq_along(filter.tags)){
          tmp.col <- paste0('filter.tag.', idx)
          tmp.data <- tmp.data[!is.na(get(tmp.col))]
        }
      }
      if(!is.null(groups.to.filter)){
        for(idx in seq_along(filter.tags)){
          tmp.col <- paste0('filter.tag.', idx)
          if(keep) tmp.data <- tmp.data[get(tmp.col) %in% groups.to.filter[[tmp.col]]]
          else     tmp.data <- tmp.data[!get(tmp.col) %in% groups.to.filter[[tmp.col]]]
        }
      }
      # drop filter columns
      keep.cols <- colnames(tmp.data)[!str_detect(colnames(tmp.data), '^filter\\.tag\\.\\d+$')]
      tmp.data <- tmp.data[, ..keep.cols]
    }
    
    # ---- long format
    tmp.data <- as.data.table(tidyr::gather(tmp.data, key='feature', value='exp.val', -groups.tag))
    
    # ---- threshold
    if(!is.null(feature.thold)){
      if(!is.numeric(feature.thold)) stop("feature.thold must be numeric")
      tmp.data <- tmp.data[exp.val > feature.thold]
    }
    
    # ---- summary stats
    tmp.data <- tmp.data[, .(
      exp.mean = mean(exp.val, na.rm=TRUE),
      burst.freq = .SD[exp.val > 0, .N] * 100 / .N
    ), by=.(groups.tag, feature)]
    
    # ---- z-score row-wise
    if(do.norm){
      wide <- tidyr::spread(tmp.data[, .(groups.tag, feature, exp.mean)], key=groups.tag, value=exp.mean)
      wide <- as.data.frame(wide)
      rownames(wide) <- wide$feature
      wide$feature <- NULL
      wide2 <- t(scale(t(as.matrix(wide))))
      wide2 <- as.data.frame(wide2)
      wide2$feature <- rownames(wide2)
      wide2 <- as.data.table(wide2)
      wide2 <- as.data.table(tidyr::gather(wide2, key='groups.tag', value='exp.mean', -feature))
      tmp.data <- merge(tmp.data[, .(groups.tag, feature, burst.freq)], wide2, by=c("feature","groups.tag"))
    }
    
    # ---- clamp color
    if(!is.null(col.min) && !is.null(col.max) && col.min > col.max) stop("col.min must be <= col.max")
    if(!is.null(col.min)) tmp.data[exp.mean < col.min, exp.mean := col.min]
    if(!is.null(col.max)) tmp.data[exp.mean > col.max, exp.mean := col.max]
    
    # ---- hot/cold scale (only meaningful after z-score)
    if(do.norm && all(is.na(this.color.scale))){
      col.breaks <- seq(min(tmp.data$exp.mean), max(tmp.data$exp.mean), length.out=100)
      mid.point <- which.min(abs(col.breaks - 0))
      this.color.scale.1 <- colorRampPalette(c('blue','mediumblue','white'))(mid.point)
      this.color.scale.2 <- colorRampPalette(c('white','yellow','gold'))(100-(mid.point+1))
      this.color.scale <- c(this.color.scale.1, this.color.scale.2)
    }
    
    # ---- clamp size
    if(!is.na(size.min) && !is.na(size.max) && size.min > size.max) stop("size.min must be <= size.max")
    if(!is.na(size.min)) tmp.data[burst.freq < size.min, burst.freq := size.min]
    if(!is.na(size.max)) tmp.data[burst.freq > size.max, burst.freq := size.max]
    
    # ---- factor order (mixedsortの代わりに base の sort)
    if(is.null(groups.order)){
      new.lvls <- sort(unique(as.character(tmp.data$groups.tag)))
      tmp.data$groups.tag <- factor(as.character(tmp.data$groups.tag), levels=new.lvls)
    } else {
      if(!(all(groups.order %in% unique(as.character(tmp.data$groups.tag))) &&
          all(unique(as.character(tmp.data$groups.tag)) %in% groups.order))){
        stop("groups.order must match observed groups exactly.")
      }
      tmp.data$groups.tag <- factor(as.character(tmp.data$groups.tag), levels=as.character(groups.order))
    }
    
    # ---- feature order
    new.lvls <- feats_use
    tmp.data$feature <- factor(tmp.data$feature, levels=rev(new.lvls))
    
    # ---- plot
    tmp.ggplot <- ggplot(tmp.data, aes(x=groups.tag, y=feature, size=burst.freq, fill=exp.mean)) +
      geom_point(stroke=0.7, shape=21, color='black') +
      scale_fill_gradientn(colors=this.color.scale) +
      scale.func(range=c(0, dot.scale), limits=c(size.min, size.max)) +
      labs(x='Groups', y='Gene', fill='Average expression', size='Percent expressed') +
      theme_classic()
    
    if(is.null(file.name)){
      return(tmp.ggplot)
    } else {
      grDevices::pdf(file=file.name)
      print(tmp.ggplot)
      grDevices::dev.off()
      return(invisible(NULL))
    }
  }

  km_rds_spa_clean@meta.data$orig.group_coclust <- case_when(km_rds_spa_clean@meta.data$orig.sub_celltype_coclust_clean2 %in% c("CD4_TCM", "CD4_TFR", "CD4_THIFNR", "CD4_TRM", "CD4_Tregs") ~ 'other_CD4', km_rds_spa_clean@meta.data$orig.sub_celltype_coclust_clean2 %in% c("CD8_GZMKp", "CD8_IFNR", "CD8_TCM", "CD8_TEM", "CD8_KIRp", "CD8_TRM", "CD8_TRM_TCR_eng") ~ 'CD8_pops', km_rds_spa_clean@meta.data$orig.sub_celltype_coclust_clean2 %in% c("NK_CD16n_circulation_markers", "NK_CD16n_residency_markers", "NK_CD16p") ~ 'NK_pops', km_rds_spa_clean@meta.data$orig.sub_celltype_coclust_clean2 %in% c("Myeloid", "Plasma_cells", "Bcells") ~ "other_immune", TRUE ~ km_rds_spa_clean@meta.data$orig.sub_celltype_coclust_clean2)


  table(km_rds_spa_clean@meta.data$orig.sub_celltype_coclust_clean2, km_rds_spa_clean@meta.data$orig.group_coclust, useNA='always')
  

  grp_order <- c("CD4_TFH", "other_CD4", "CD8_pops", "NK_pops", "other_immune", "proliferating_cells", "Endothelial_cells", "CAF", "Tumor_cells")
  ####IL21 expression
  ##Dot plot
  tmp.ggplot <- dot.plot_xenium(
    seurat.obj = km_rds_spa_clean,
    assay = "Xenium",           # ここ重要
    layer = "data",             # v5なら "counts" / "data" を明示が安全
    features = c('IL21'),
    groups.tag = "orig.group_coclust",
    do.norm = TRUE,
    groups.order = grp_order,
    this.color.scale='hot.and.cold'
  )
  tmp.ggplot
  # Output.
  #tmp.width <- (length(these.pops[[t.subset]]) * 0.5) + 0.3
  
  tmp.height <- 5
  # tmp.width <- 5
  #tmp.height <- (length(these.markers[[t.subset]]) * 0.4) + 0.3
  tmp.width <- (length(unique(km_rds_spa_clean@meta.data$orig.group_coclust)) * 0.4) + 0.3

  tmp.ggplot <- tmp.ggplot +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )
  tmp.ggplot 
  tmp.lab <- paste0('legend_IL21', '_Markers_DotPlot_Opt-A')
  publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=paste0(getwd(), "/IL21_expression"), file.name=tmp.lab, type='pdf', blank.comp=blank.complement.2, do.legend=TRUE,
    width=tmp.width, height=tmp.height, c.width=tmp.width+2, c.height=tmp.height, legend.height=5, legend.width=0.5
  )
  tmp.height <- (1 * 0.4) + 0.3
  tmp.lab <- paste0('IL21', '_Markers_DotPlot_Opt-A')
  publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=paste0(getwd(), "/IL21_expression"), file.name=tmp.lab, type='pdf', blank.comp=blank.complement.2, do.legend=TRUE,
    width=tmp.width, height=tmp.height, c.width=tmp.width+2, c.height=tmp.height, legend.height=5, legend.width=0.5
  )

  #dot plot all populations

  grep_order2 <- c("NK_CD16p", "NK_CD16n_circulation_markers", "NK_CD16n_residency_markers", "CD8_KIRp", "CD8_TRM", "CD8_TRM_TCR_eng", "CD8_GZMKp", "CD8_TCM", "CD8_TEM", "CD8_IFNR", "CD4_THIFNR", "CD4_TFH", "CD4_TCM", "CD4_TFR", "CD4_TRM", "CD4_Tregs", "Myeloid", "Plasma_cells", "Bcells", "proliferating_cells", "Endothelial_cells", "CAF", "Tumor_cells")

  km_rds_spa_clean_sub <- subset(km_rds_spa_clean, subset = orig.celltype == 'Tcells_NK')

  grep_order3 <- c("NK_CD16p", "NK_CD16n_circulation_markers", "NK_CD16n_residency_markers", "CD8_KIRp", "CD8_TRM_TCR_eng", "CD8_TRM", "CD8_GZMKp", "CD8_TCM", "CD8_TEM", "CD8_IFNR", "CD4_THIFNR", "CD4_Tregs", "CD4_TFR", "CD4_TFH", "CD4_TCM", "CD4_TRM")

  tmp.ggplot_all <- dot.plot_xenium(
    seurat.obj = km_rds_spa_clean_sub,
    assay = "Xenium",           # ここ重要
    layer = "data",             # v5なら "counts" / "data" を明示が安全
    features = c("FGFBP2", "FCGR3A", "TRDC", "KLRC1", "KLRD1", "ITGAE", "ZNF683", "ITGA1", "CD8A", "GZMK", "TCF7", "S1PR5", "IFIT1", "MX1", "CD4", "FOXP3", "CTLA4", "RTKN2", "IKZF2", "PTPN13", "CXCL13", "TCF7", "S1PR1", "LEF1", "MAL", "CD40LG"),
    groups.tag = "orig.sub_celltype_coclust_clean2",
    do.norm = TRUE,
    groups.order = grep_order3,
    this.color.scale='hot.and.cold'
  )
  # tmp.ggplot_all
  # Output.
  #tmp.width <- (length(these.pops[[t.subset]]) * 0.5) + 0.3
  tmp.lab <- "TNK_subset"
  tmp.ggplot_all <- tmp.ggplot_all +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    ) + scale_size(range = c(2, 8), breaks = c(10, 30, 70))

  tmp.ggplot_all <- tmp.ggplot_all +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  scale_size(
    range = c(0, 12),
    breaks = c(10, 30, 80),
    limits = c(0, 80)
  )

  pdf(paste0(getwd(), "/IL21_expression/", tmp.lab, ".pdf"), height = 10.5, width = 8.2)
  print(tmp.ggplot_all)
  dev.off()

  km_rds_spa_clean@meta.data$orig.group_coclust2 <- case_when(km_rds_spa_clean@meta.data$orig.sub_celltype_coclust_clean2 %in% c("CD4_TCM", "CD4_TFR", "CD4_THIFNR", "CD4_TRM", "CD4_Tregs", "CD4_TFH") ~ 'CD4_pops', km_rds_spa_clean@meta.data$orig.sub_celltype_coclust_clean2 %in% c("CD8_GZMKp", "CD8_IFNR", "CD8_TCM", "CD8_TEM", "CD8_KIRp", "CD8_TRM", "CD8_TRM_TCR_eng") ~ 'CD8_pops', km_rds_spa_clean@meta.data$orig.sub_celltype_coclust_clean2 %in% c("NK_CD16n_circulation_markers", "NK_CD16n_residency_markers", "NK_CD16p") ~ 'NK_pops', TRUE ~ km_rds_spa_clean@meta.data$orig.sub_celltype_coclust_clean2)

  grp_order2 <- c("NK_pops", "CD8_pops", "CD4_pops", "Myeloid", "Plasma_cells", "Bcells", "proliferating_cells", "Endothelial_cells", "CAF", "Tumor_cells")


  tmp.ggplot_all_gen <- dot.plot_xenium(
    seurat.obj = km_rds_spa_clean,
    assay = "Xenium",           # ここ重要
    layer = "data",             # v5なら "counts" / "data" を明示が安全
    features = c(c("TRDC", "KLRD1", "CD3E", "CD8A", "CD8B", "CD4"), rev(c("EPCAM", "COL5A1", "FLT4", "CD34", "PECAM1", "TYMS", "MKI67", "MS4A1", "CD79A", "MZB1", "FCGR2A", "CD14"))),
    groups.tag = "orig.group_coclust2",
    do.norm = TRUE,
    groups.order = grp_order2,
    this.color.scale='hot.and.cold'
  )
  # tmp.ggplot_all_gen
  # Output.
  #tmp.width <- (length(these.pops[[t.subset]]) * 0.5) + 0.3
  tmp.lab <- "General_subset"
  tmp.ggplot_all_gen <- tmp.ggplot_all_gen +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    ) + scale_size(range = c(2, 8), breaks = c(10, 30, 70))

  tmp.ggplot_all_gen <- tmp.ggplot_all_gen +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  scale_size(
    range = c(0, 12),
    breaks = c(10, 30, 80),
    limits = c(0, 80)
  )

  pdf(paste0(getwd(), "/IL21_expression/", tmp.lab, ".pdf"), height = 7.5, width = 6)
  print(tmp.ggplot_all_gen)
  dev.off()

}

{cat(redb("### Phenoptr distance matrix ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  
  #donor information per TMA
  donor_TMA_id <- c('B12' = "1568",
                    'B13' = "1568",
                    'D10' = '1620',
                    'D11' = '1620',
                    'D12' = '1621',
                    'D13' = '1621',
                    'E8' = '1636')

  km_rds_spa_clean@meta.data$orig.donor <- donor_TMA_id[km_rds_spa_clean@meta.data$orig.TMA_core]

  table(km_rds_spa_clean$orig.donor, km_rds_spa_clean$orig.TMA_core)

  ## Calculate distances from T cells to all other celltypes
  source("/home/fcastaneda/bin/clustering/R/stats_summary_table.R")

  ## Calculating the distance matrix
    slide_xe <- km_rds_spa_clean
    slide_xe$donor_TMA <- paste0(slide_xe$orig.TMA_core, "_", slide_xe$orig.donor)
    # cells coordinates (based on fov)
    coords <- Seurat::GetTissueCoordinates(slide_xe, image = "fov")
    colnames(coords)[3] <- "cellname"
    #sanity check for matrix
    all(coords$cellname %in% slide_xe@meta.data$cellname)
    # coords with metadata info for phenoptr
    coords2 <- left_join(coords, slide_xe@meta.data[, c("cellname", "orig.sub_celltype_coclust_clean2", "donor_TMA")], by = 'cellname')
    #renaming columns for function
    colnames(coords2) <- c("Cell X Position", "Cell Y Position", "Cell ID", "Phenotype", "sample")
    coords2 <- coords2 %>% select("Cell ID", "Cell X Position", "Cell Y Position", "Phenotype", "sample") 
    rownames(coords2) <- coords2$'Cell ID'

    head(coords2)
    # Distance matrix generated per donor
     list_dist_mtx <- lapply(unique(coords2$sample), function(donor) {
        # donor <- unique(coords2$sample)[1]
        donor_distances <- coords2 %>% filter(sample == donor) 
        dist_mtx <- distance_matrix(donor_distances)
        return(dist_mtx)
    })

    names(list_dist_mtx) <- unique(coords2$sample)
    #celltypes in rds
    celltypes_comp_orig <- unique(slide_xe$orig.sub_celltype_coclust_clean2)
    from_celltypes_comp_orig <-  unique(slide_xe$orig.sub_celltype_coclust_clean2)
    # celltypes_comp_orig <- "CD4_TFH"
    dir.create("distance_matrix")

    #matrix per celltype
    all_list <- lapply(celltypes_comp_orig, function(from_celltype){ 
        
        # from_celltype <-  "Tumor_cells"
        #selecting based (from) celltype for matrix
        cat("based (from) celltype: ", from_celltype, "\n")
          
        matrix_donor <- lapply(names(list_dist_mtx), function(donor){   
          # donor <- names(list_dist_mtx)[1]

          cat("sample: ", donor, "\n")
          #sample specific matrix
          sample_matrix <- list_dist_mtx[[donor]]
          # all celltypes but based 'from' celltype
          celltypes_comp <- as.character(from_celltypes_comp_orig)
          celltypes_comp <- grep(from_celltype, celltypes_comp, value=TRUE, invert=TRUE)

          distance_to_closest_cell <- lapply(celltypes_comp, function(to_celltype){ 
            # to_celltype <- "CD8"

            # performing matrix to any other celltype but based (from) celltype
            # cat(to_celltype, "\n")

            # from cell type
            from_celltype_cells <- slide_xe@meta.data %>% filter(orig.sub_celltype_coclust_clean2 == from_celltype  & donor_TMA ==  donor) %>% pull(cellname)
            # to cell type
            to_celltype_cells <- slide_xe@meta.data %>% filter(orig.sub_celltype_coclust_clean2 == to_celltype & donor_TMA ==  donor) %>% pull(cellname)

            if(length(from_celltype_cells) > 1 & length(to_celltype_cells) > 1) {
              #constructing matrix per to_celltype based on from celltype
              sample_matrix2 <- sample_matrix[to_celltype_cells, from_celltype_cells]
              to_ret <- data.frame(
                distance = rowMins(sample_matrix2),
                to = to_celltype, 
                cell= names(rowMins(sample_matrix2)),
                donor_sample = donor,
                donor_id = sub("_[^_]*$", "", donor))

            } else {
              to_ret <- data.frame(
                distance = 0,
                to = to_celltype, 
                cell= "no_cells",
                donor_sample = donor,
                donor_id = sub("_[^_]*$", "", donor))
            }
            return(to_ret) 
          })

          distance_from_tcells <- Reduce(rbind, distance_to_closest_cell)
          return(distance_from_tcells)

        })

        donor_matrix <-  Reduce(rbind,matrix_donor)
        return(donor_matrix)

    })

    names(all_list) <- celltypes_comp_orig

  lapply(seq_along(names(all_list)), function(x) {
    # x <- 1
    cp_km <- names(all_list)[[x]]
    m_df <- all_list[[cp_km]]

    if(cp_km %in% m_df$to) {stop(cat("something went wrong :c", "\n"))} else {cat("all cool for ", cp_km, "!", "\n")}
    
    write.csv(m_df, paste0(getwd(), "/distance_matrix/distance_matrix_", cp_km, ".csv"))
    
    return(NULL)

  })

}

{cat(redb("### Distance to closest immune celltype for each structural celltype per group and as a whole ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

  # For TFH and trNK

    colors_dist <- c('CD4_TFH' = "#0075ff",
                        'CD4_Tregs' = "#69ffba",
                        'CD4_THIFNR' = "#927155",
                        'CD4_TFR' = "#08d2ff",
                        'CD4_TRM' = "#5daa8e",
                        'CD4_TCM' = "#ffb2a2",
                        'CD4_CTL_GZMKp' = "#00650c",
                        'CD4_CTL_GZMKn' = "#0096a6",
                        'CD8_TRM_TCR_eng' = "#7996ff",
                        'CD8_TCM' = "#00dfc2",
                        'CD8_GZMKp' = "#ca2800",
                        'CD8_TRM' = "#000096",
                        'CD8_IFNR' = "#006986",
                        'CD8_TEM' = "#d792aa",
                        'CD8_MAIT' = "#9a001c",
                        'CD8_KIRp' = "#96698a",
                        'NK_CD16n_residency_markers' = "#ff316d",
                        'NK_CD16p' = "#69c66d",
                        'NK_CD16n_circulation_markers' = "#ffe779",
                        'NK_CD16p_mature' = "#aa9ec6",
                        'NK_ILC3' = "#008671",
                        'NK_IFN' = "#9600ca",
                        'NK_CD16p_active' = "#a2f3ff",
                        'Tumor_cells' = '#5D75AE',
                        'Bcells' = '#FFEB72',
                        'proliferating_cells' = '#5D95FB',
                        'CAF' = '#24aed2',
                        'Plasma_cells' = '#92cea6',
                        'Myeloid' = '#fbb6ca',
                        'Endothelial_cells' = '#d2b6ff')

  ############################################################  
  ########### Distance and quartile analysis  ################
  ############################################################  

  ## To calculate the distance to populations from a given celltype 

  distance_files <- list.files("/home/donaldo/bioadhoc/spatial/KoMa/tumor_merge/results/figures/KoMa_NKpaper_figures_v3/distance_matrix/", pattern="distance_matrix_", recursive=TRUE, full.names=TRUE)

  celltypes_comp_orig2 <- c("CD4_TFH", "NK_CD16n_residency_markers")

  dir.create("distance_matrix/distance_per_pop")

  lapply(celltypes_comp_orig2, function(per_celltype){ 
      # per_celltype <- celltypes_comp_orig2[2]

      cat(per_celltype, "\n")

      #selecting specific matrix
      files_celltype <- grep(per_celltype, distance_files, value=TRUE)
      distance_to_celltypes <- read.csv(files_celltype, row.names = 1)
      #removing 'never' for disease status and no_cells found
      distance_to_celltypes <- distance_to_celltypes %>% filter(cell != "no_cells")
      #documented this one for distance to all cells.
      # distance_to_celltypes <- distance_to_celltypes %>% filter(to %in% c("CD8_TRM", "CD8_TEM", "CD8_TEFF", "CD8_GZMK", "CD8_TCM", "CD8_circKIR", "CD8_KIRp", "CD8_IFNR", "CD8_MAIT", "CD4_TRM", "CD4_Treg", "CD4_TCM", "CD4_TFH", "CD4_THIFNR", "CD4_conv_CTL"))

      # distance_to_celltypes <-Reduce(rbind, lapply(files_celltype, function(distances_files){
      #   distances_slides <- read.csv(distances_files)
      #   distances_slides$slide <- str_extract(distances_files, "slide1|slide2")
      #   return(distances_slides)
      #   }))

      if(per_celltype == "NK_CD16n_residency_markers") {
        order_ct_plot <- c('CD4_TFH', 'CD4_TCM', 'CD4_TFR', 'CD4_THIFNR', 'CD4_TRM', 'CD4_Tregs', 'CD8_GZMKp', 'CD8_IFNR', 'CD8_TCM', 'CD8_TEM', 'CD8_KIRp', 'CD8_TRM', 'CD8_TRM_TCR_eng', 'NK_CD16n_circulation_markers', 'NK_CD16p', 'Myeloid', 'CAF', 'Endothelial_cells', 'Plasma_cells', 'Bcells', 'Tumor_cells', 'proliferating_cells')
      } else {
        order_ct_plot <- c('NK_CD16n_residency_markers', 'CD4_TCM', 'CD4_TFR', 'CD4_THIFNR', 'CD4_TRM', 'CD4_Tregs', 'CD8_GZMKp', 'CD8_IFNR', 'CD8_TCM', 'CD8_TEM', 'CD8_KIRp', 'CD8_TRM', 'CD8_TRM_TCR_eng', 'NK_CD16n_circulation_markers', 'NK_CD16p', 'Myeloid', 'CAF', 'Endothelial_cells', 'Plasma_cells', 'Bcells', 'Tumor_cells', 'proliferating_cells')
      }

      # #plots per celltype without considering disease group
      # pdf(file = paste0(getwd(), "/distance_matrix/distance_per_pop/distance_per_celltype_fropm_", per_celltype, ".pdf"), width = 10)
      # distance_plot2 <- ggplot(distance_to_celltypes, aes(x = factor(to, levels = order_ct_plot), y = distance, fill=to)) + 
      #       geom_jitter(position = position_jitter(0.2), alpha=0.1, size=0.1) + 
      #       geom_violin(width = 1) +
      #       scale_fill_manual(values=colors_dist) +
      #       geom_boxplot(width = 0.1, outlier.shape = NA) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(paste0("Distance from ", per_celltype)) + scale_y_continuous(trans = log2_trans()) #facet_wrap_paginate(~ orig.disease2_wo_batch, ncol = 1, nrow = 1, scales = "free")
      # print(distance_plot2)
      # dev.off()

      # #keeping only cell populations with at least 10 cells and 2 samples
      cell_count_ct_tma <- table(km_rds_spa_clean$orig.sub_celltype_coclust_clean2, km_rds_spa_clean$orig.TMA_core) %>% data.frame()

      # cell_count_sub <- cell_count_ct_tma %>% group_by(Var1) %>% filter(sum(Freq >= 10) >= 2) %>% ungroup() %>% distinct(Var1) %>% pull(Var1) %>% as.character()

      # pdf(file = paste0(getwd(), "/distance_matrix/distance_per_pop/distance_per_celltype_fropm_", per_celltype, "_10_count_2_samples_thershold.pdf"), width = 10)
      # distance_plot_thr <- ggplot(distance_to_celltypes %>% filter(to %in% cell_count_sub), aes(x = factor(to, levels = order_ct_plot), y = distance, fill=to)) + 
      #       geom_jitter(position = position_jitter(0.2), alpha=0.1, size=0.1) + 
      #       geom_violin(width = 1) +
      #       scale_fill_manual(values=colors_dist) +
      #       geom_boxplot(width = 0.1, outlier.shape = NA) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(paste0("Distance from ", per_celltype)) + scale_y_continuous(trans = log2_trans()) #facet_wrap_paginate(~ orig.disease2_wo_batch, ncol = 1, nrow = 1, scales = "free")
      # print(distance_plot_thr)
      # dev.off()

      cell_count_sub2 <- cell_count_ct_tma %>% group_by(Var1) %>% filter(sum(Freq >= 10) >= 3) %>% ungroup() %>% distinct(Var1) %>% pull(Var1) %>% as.character()

      # pdf(file = paste0(getwd(), "/distance_matrix/distance_per_pop/distance_per_celltype_fropm_", per_celltype, "_10_count_3_samples_thershold.pdf"), width = 10)
      # distance_plot_thr <- ggplot(distance_to_celltypes %>% filter(to %in% cell_count_sub2), aes(x = factor(to, levels = order_ct_plot), y = distance, fill=to)) + 
      #       geom_jitter(position = position_jitter(0.2), alpha=0.1, size=0.1) + 
      #       geom_violin(width = 1) +
      #       scale_fill_manual(values=colors_dist) +
      #       geom_boxplot(width = 0.1, outlier.shape = NA) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(paste0("Distance from ", per_celltype)) + scale_y_continuous(trans = log2_trans()) #facet_wrap_paginate(~ orig.disease2_wo_batch, ncol = 1, nrow = 1, scales = "free")
      # print(distance_plot_thr)
      # dev.off()

      # pdf(file = paste0(getwd(), "/distance_matrix/distance_per_pop/distance_per_celltype_fropm_", per_celltype, "_10_count_3_samples_thershold_mod.pdf"), width = 12)
      # distance_plot_thr <- ggplot(distance_to_celltypes %>% filter(to %in% cell_count_sub2), aes(x = factor(to, levels = order_ct_plot), y = distance, fill=to)) + 
      #       geom_jitter(position = position_jitter(0.2), alpha=0.1, size=0.1) + 
      #       geom_violin(width = 1, alpha = 0.8) +
      #       scale_fill_manual(values=colors_dist) +
      #       # scale_color_manual(values=colors_dist) +
      #       geom_boxplot(width = 0.15, outlier.shape = NA) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(paste0("Distance from ", per_celltype)) + scale_y_continuous(trans = log2_trans(), breaks = c(1, 10, 50, 100, 200, 500), labels = c(1, 10, 50, 100, 200, 500)) + coord_cartesian(ylim = c(NA, 500)) #+ geom_hline(yintercept = 100, linetype = "dotted", color = "red") + #facet_wrap_paginate(~ orig.disease2_wo_batch, ncol = 1, nrow = 1, scales = "free")
      # print(distance_plot_thr)
      # dev.off()

      pdf(file = paste0(getwd(), "/distance_matrix/distance_per_pop/distance_per_celltype_fropm_", per_celltype, "_10_count_3_samples_thershold_mod_TFHonly.pdf"), width = 3, height = 5)
      distance_plot_thr <- ggplot(distance_to_celltypes %>% filter(to == 'CD4_TFH'), aes(x = factor(to, levels = order_ct_plot), y = distance, fill=to)) + 
            # geom_jitter(position = position_jitter(0.2), alpha=0.5, size=0.2) + 
            geom_violin(width = 1, alpha = 0.6) +
            scale_fill_manual(values=colors_dist) +
            # scale_color_manual(values=colors_dist) +
            geom_boxplot(width = 0.15, outlier.shape = NA) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(paste0("Distance from ", per_celltype)) + scale_y_continuous(trans = log2_trans(), breaks = c(1, 10, 50, 100, 200, 500), labels = c(1, 10, 50, 100, 200, 500)) + coord_cartesian(ylim = c(NA, 500)) #+ geom_hline(yintercept = 100, linetype = "dotted", color = "red") + #facet_wrap_paginate(~ orig.disease2_wo_batch, ncol = 1, nrow = 1, scales = "free")
      print(distance_plot_thr)
      dev.off()

    return(NULL)

  })

}
