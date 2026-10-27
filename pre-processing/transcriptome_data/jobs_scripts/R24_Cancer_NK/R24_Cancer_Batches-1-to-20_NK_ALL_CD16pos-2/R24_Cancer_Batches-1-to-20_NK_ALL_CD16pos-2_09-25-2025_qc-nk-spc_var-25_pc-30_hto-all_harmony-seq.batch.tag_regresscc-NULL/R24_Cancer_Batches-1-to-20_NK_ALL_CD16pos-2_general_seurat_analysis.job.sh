#!/bin/sh
#SBATCH --job-name=R24_Cancer_Batches-1-to-20_NK_ALL_CD16pos-2_SA
#SBATCH --output=/mnt/BioAdHoc/Groups//vd-vijay/vfajardo/R24/jobs_scripts/seurat_analysis//R24_Cancer_NK/R24_Cancer_Batches-1-to-20_NK_ALL_CD16pos-2/R24_Cancer_Batches-1-to-20_NK_ALL_CD16pos-2_09-25-2025_qc-nk-spc_var-25_pc-30_hto-all_harmony-seq.batch.tag_regresscc-NULL/R24_Cancer_Batches-1-to-20_NK_ALL_CD16pos-2_general_seurat_analysis.out.txt
#SBATCH --error=/mnt/BioAdHoc/Groups//vd-vijay/vfajardo/R24/jobs_scripts/seurat_analysis//R24_Cancer_NK/R24_Cancer_Batches-1-to-20_NK_ALL_CD16pos-2/R24_Cancer_Batches-1-to-20_NK_ALL_CD16pos-2_09-25-2025_qc-nk-spc_var-25_pc-30_hto-all_harmony-seq.batch.tag_regresscc-NULL/R24_Cancer_Batches-1-to-20_NK_ALL_CD16pos-2_general_seurat_analysis.err.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=250g
#SBATCH --time=160:00:00

echo -e "\n\n######### ------------- Job to Run Seurat Analysis ------------ #########\n\n"
echo -e "### --------------------------- Arguments --------------------------- ###"
echo -e "### --------------------- Job general arguments --------------------- ###"
echo -e "Project name: R24_Cancer_Batches-1-to-20_NK_ALL_CD16pos-2"
echo -e "Job output directory: /mnt/BioAdHoc/Groups//vd-vijay/vfajardo/R24/jobs_scripts/seurat_analysis//R24_Cancer_NK/R24_Cancer_Batches-1-to-20_NK_ALL_CD16pos-2/R24_Cancer_Batches-1-to-20_NK_ALL_CD16pos-2_09-25-2025_qc-nk-spc_var-25_pc-30_hto-all_harmony-seq.batch.tag_regresscc-NULL"
echo -e "### ------------------- Seurat analysis arguments ------------------- ###"
echo -e "Seurat script version (absolute path): /path/to/home/scripts/seurat_analysis/general_seurat_analysis.2.4.1.R"
echo -e "Project name: R24_Cancer_Batches-1-to-20_NK_ALL_CD16pos-2"
echo -e "Absolute path to count matrix file: "
echo -e "Feature ID: ensembl"
echo -e "Minimum and maximum thresholds for counts number per cell: 1000 and 6000"
echo -e "Minimum and maximum thresholds for genes number per cell: 600 and 2800"
echo -e "Maximum mitochondrial genes percentage per cell: 10"
echo -e "Normalization method: LogNormalize"
echo -e "Variable features selection method: vst"
echo -e "Number of most variable genes to be considered: 25"
echo -e "Resolution argument for cluster analysis: 'c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)'"
echo -e "PCs to be taken into account for dimensionality reduction: 'c(30)'\n\n"
echo -e "Mean cutoff: 0.01"
echo -e "Should likely low quality cells be filtered out? TRUE"
echo -e "Absolute path to markers file (or NULL value): /path/to/home/scripts/seurat_analysis/general_data/NK-cell_markers_ensembl.1.0.RData"
echo -e "### --------------------- Annotations arguments --------------------- ###"
echo -e "Should annotations be added? TRUE (If FALSE, this chunk's variables' values don't really matter)."
echo -e "Absolute path to annotations table: /path/to/user/sequencing_data/05-14-2022/aggr/data/R24_Cancer_Batches-1-to-20_NK_ALL_aggr_table_annotated.1.0.csv"
echo -e "Lane ID: 'species.tag;cell.type.tag;donors.sort.batch.tag;donors.no.tag;tissue.tag'"
echo -e "10X Chromium Batch ID: chrom.batch.tag"
echo -e "Sequencing batch ID: seq.batch.tag"
echo -e "HTO ID: hashtag.tag"
echo -e "Demuxlet ID: "
echo -e "Overlay ID: donor.tag"
echo -e "### ---------------------- Filtering arguments ---------------------- ###"
echo -e "Should there be any kind of preprocess filetring? TRUE (If FALSE, this chunk's variables' values don't really matter)."
echo -e "Tags-related criteria: --TagsCriteria /mnt/BioAdHoc/Groups//vd-vijay/vfajardo/R24/paper_developments/R24_Cancer/data_for_subsetting/custom_subsets//R24_Cancer_NK//R24_Cancer_Batches-1-to-20_NK_ALL_CD16pos-1/R24_Cancer_Batches-1-to-20_NK_ALL_CD16pos-1_09-24-2025_qc-nk-spc_var-30_pc-30_hto-all_harmony-seq.batch.tag_regresscc-NULL/CD16pos-2/TagsSubsetCriteria.csv"
echo -e "Features-related criteria: "

echo -e "Here starts analysis with R package seurat...\n\n"

module load R/3.6.1

Rscript /path/to/home/scripts/seurat_analysis/general_seurat_analysis.2.4.1.R \
--ReportsPath /mnt/BioAdHoc/Groups//vd-vijay/vfajardo/R24/seurat_analysis//R24_Cancer_NK/R24_Cancer_Batches-1-to-20_NK_ALL_CD16pos-2/R24_Cancer_Batches-1-to-20_NK_ALL_CD16pos-2_09-25-2025_qc-nk-spc_var-25_pc-30_hto-all_harmony-seq.batch.tag_regresscc-NULL --PrjName R24_Cancer_Batches-1-to-20_NK_ALL_CD16pos-2 --RAM 250 \
--DataFile /path/to/user/R24/seurat_analysis/R24_Cancer_NK/R24_Cancer_Batches-1-to-20_NK_ALL_CD16pos-1/R24_Cancer_Batches-1-to-20_NK_ALL_CD16pos-1_09-24-2025_qc-nk-spc_var-30_pc-30_hto-all_harmony-seq.batch.tag_regresscc-NULL/seurat_objects/SeuratObjectForPrjR24_Cancer_Batches-1-to-20_NK_ALL_CD16pos-1_WithArgs_NoPCs_30.RDS --Raw10x /path/to/user/sequencing_data/05-14-2022/aggr/R24_Cancer_Batches-1-to-20_NKAndB_ALL/outs/filtered_feature_bc_matrix --InputType seurat --FeatureID ensembl \
--FilterOut TRUE --minCounts 1000 --maxCounts 6000 --minFeatures 600 --maxFeatures 2800 --maxMP 10 --IsFivePrime FALSE \
--GenNormMethod LogTransform --NormMethod LogNormalize \
--FVFsMethod vst --FeatsForDSA 25 --MeanCutoff 0.01 --PropCutoff 0.001 \
--PCs 'c(30)' --ForHarmony 'c("seq.batch.tag")' --VarsToRegress 'c("nCount_RNA", "percent.mt")'  --DimReduction 'c("umap")' --Resolution 'c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)' \
--MarkersFile /path/to/home/scripts/seurat_analysis/general_data/NK-cell_markers_ensembl.1.0.RData --DEA FALSE \
--DoAnnotate TRUE --AggrTable /path/to/user/sequencing_data/05-14-2022/aggr/data/R24_Cancer_Batches-1-to-20_NK_ALL_aggr_table_annotated.1.0.csv --DonorsMetaData /mnt/BioAdHoc/Groups//vd-vijay/vfajardo/R24/paper_developments/R24_Cancer/donors_metadata/GeneralMetaDataProcessedTableForR24CancerBatches-1-to-20.1.0.csv  --MergeBy donor.tag --LaneID 'species.tag;cell.type.tag;donors.sort.batch.tag;donors.no.tag;tissue.tag' --MultAnnPerLane TRUE --ChromBatchID chrom.batch.tag --SeqBatchID seq.batch.tag --HTOID hashtag.tag --OverlayID donor.tag  --NewMeta /path/to/user/R24/paper_developments/R24_Cancer/paper_items/labels_for_transfer/NK_v2.csv \
--DoSubset TRUE --TagsCriteria /mnt/BioAdHoc/Groups//vd-vijay/vfajardo/R24/paper_developments/R24_Cancer/data_for_subsetting/custom_subsets//R24_Cancer_NK//R24_Cancer_Batches-1-to-20_NK_ALL_CD16pos-1/R24_Cancer_Batches-1-to-20_NK_ALL_CD16pos-1_09-24-2025_qc-nk-spc_var-30_pc-30_hto-all_harmony-seq.batch.tag_regresscc-NULL/CD16pos-2/TagsSubsetCriteria.csv 

echo -e "Job completed!\nCheck for errors if any."
