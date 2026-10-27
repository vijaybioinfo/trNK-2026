#!/bin/sh
#SBATCH --job-name=R24_Cancer_Batches-1-to-20_NK_ALL_SA
#SBATCH --output=/path/to/user/R24/jobs_scripts/seurat_analysis//R24_Cancer_NK/R24_Cancer_Batches-1-to-20_NK_ALL/R24_Cancer_Batches-1-to-20_NK_ALL_09-24-2025_qc-nk-spc_var-30_pc-30_hto-all_harmony-seq.batch.tag_regresscc-NULL/R24_Cancer_Batches-1-to-20_NK_ALL_general_seurat_analysis.out.txt
#SBATCH --error=/path/to/user/R24/jobs_scripts/seurat_analysis//R24_Cancer_NK/R24_Cancer_Batches-1-to-20_NK_ALL/R24_Cancer_Batches-1-to-20_NK_ALL_09-24-2025_qc-nk-spc_var-30_pc-30_hto-all_harmony-seq.batch.tag_regresscc-NULL/R24_Cancer_Batches-1-to-20_NK_ALL_general_seurat_analysis.err.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=300g
#SBATCH --time=120:00:00

echo -e "\n\n######### ------------- Job to Run Seurat Analysis ------------ #########\n\n"
echo -e "### --------------------------- Arguments --------------------------- ###"
echo -e "### --------------------- Job general arguments --------------------- ###"
echo -e "Project name: R24_Cancer_Batches-1-to-20_NK_ALL"
echo -e "Job output directory: /path/to/user/R24/jobs_scripts/seurat_analysis//R24_Cancer_NK/R24_Cancer_Batches-1-to-20_NK_ALL/R24_Cancer_Batches-1-to-20_NK_ALL_09-24-2025_qc-nk-spc_var-30_pc-30_hto-all_harmony-seq.batch.tag_regresscc-NULL"
echo -e "### ------------------- Seurat analysis arguments ------------------- ###"
echo -e "Seurat script version (absolute path): /path/to/home/scripts/seurat_analysis/general_seurat_analysis.2.3.9.R"
echo -e "Project name: R24_Cancer_Batches-1-to-20_NK_ALL"
echo -e "Absolute path to count matrix file: "
echo -e "Feature ID: ensembl"
echo -e "Minimum and maximum thresholds for counts number per cell: 1000 and 6000"
echo -e "Minimum and maximum thresholds for genes number per cell: 600 and 2800"
echo -e "Maximum mitochondrial genes percentage per cell: 10"
echo -e "Normalization method: LogNormalize"
echo -e "Variable features selection method: vst"
echo -e "Number of most variable genes to be considered: 30"
echo -e "Resolution argument for cluster analysis: 'c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)'"
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
echo -e "Tags-related criteria: --TagsCriteria /path/to/user/R24/paper_developments/R24_Cancer/data_for_subsetting/subset_hto-all/TagsSubsetCriteria.csv"
echo -e "Features-related criteria: "

echo -e "Here starts analysis with R package seurat...\n\n"

module load R/3.6.1

Rscript /path/to/home/scripts/seurat_analysis/general_seurat_analysis.2.3.9.R --ReportsPath /path/to/user/R24/seurat_analysis//R24_Cancer_NK/R24_Cancer_Batches-1-to-20_NK_ALL/R24_Cancer_Batches-1-to-20_NK_ALL_09-24-2025_qc-nk-spc_var-30_pc-30_hto-all_harmony-seq.batch.tag_regresscc-NULL --PrjName R24_Cancer_Batches-1-to-20_NK_ALL --DataFile /path/to/user/sequencing_data/05-14-2022/reanalyze/R24_Cancer_Batches-1-to-20_NK_ALL/outs/filtered_feature_bc_matrix --Raw10x /path/to/user/sequencing_data/05-14-2022/reanalyze/R24_Cancer_Batches-1-to-20_NK_ALL/outs/filtered_feature_bc_matrix --InputType matrix   --FeatureID ensembl --minCounts 1000 --maxCounts 6000 --minFeatures 600 --maxFeatures 2800 --maxMP 10 --GenNormMethod LogTransform --NormMethod LogNormalize --FVFsMethod vst --FeatsForDSA 30 --PCs 'c(30)' --ForHarmony 'c("seq.batch.tag")' --Resolution 'c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)' --IsFivePrime FALSE --MeanCutoff 0.01 --PropCutoff 0.001 --VarsToRegress 'c("nCount_RNA", "percent.mt")'  --FilterOut TRUE --MarkersFile /path/to/home/scripts/seurat_analysis/general_data/NK-cell_markers_ensembl.1.0.RData --DEA FALSE --DimReduction 'c("umap")' --RAM 300 --DoAnnotate TRUE --AggrTable /path/to/user/sequencing_data/05-14-2022/aggr/data/R24_Cancer_Batches-1-to-20_NK_ALL_aggr_table_annotated.1.0.csv --DonorsMetaData /path/to/user/R24/paper_developments/R24_Cancer/donors_metadata/GeneralMetaDataProcessedTableForR24CancerBatches-1-to-20.1.0.csv  --MergeBy donor.tag --LaneID 'species.tag;cell.type.tag;donors.sort.batch.tag;donors.no.tag;tissue.tag' --MultAnnPerLane TRUE --ChromBatchID chrom.batch.tag --SeqBatchID seq.batch.tag --HTOID hashtag.tag --OverlayID donor.tag  --DoSubset TRUE --TagsCriteria /path/to/user/R24/paper_developments/R24_Cancer/data_for_subsetting/subset_hto-all/TagsSubsetCriteria.csv 

echo -e "Job completed!\nCheck for errors if any."
