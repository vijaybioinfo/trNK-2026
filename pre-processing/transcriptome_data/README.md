Transcriptome data from single-cell RNA-seq.
===========

This process was applied to each aggregated dataset containing both NK and B cells, either separately or together. The process was first applied to the aggregated datasets containing both cell types but separately according to sample type--human lungs or tumors--to allow for lineage deconvolution on the basis of their transcriptomic profiles. Then, the good-quality cell barcodes assigned with the NK cell identity from each type of sample were merged into a separate aggregated dataset using cellranger reanalyze—-cell barcode subset and job script provided in the folder `./reanalyze`.

This section comprises the following major steps:

* Preprocessing of reads by mapping them to the reference genome and collapsing them into UMI count matrices (10x's cellranger; see further details in enclosing folder).
* Donor deconvolution based on hashtag oligonucleotide (HTO) data (see further details in enclosing folder).
* Quality control (QC).
* Data normalization and batch effect removal.
* Clustering and dimensionality reduction.
* Population definition based on differential expression analysis and manual assignment of B or NK cell subsets.

The first step is described in the methods section of our manuscript in detail. The second step is described in `./../hto_data` in this repo. Here we provide the code for the rest of the steps.

---
# Scripts and job scripts
We have performed the standard scRNA-seq analysis for each aggregated dataset by applying the [Seurat toolkit](https://satijalab.org/seurat/) (shoutout to the Seurat developers and maintainers!).<br>
Specifically, several rounds of cleaning were applied to complete effective transcriptome-based immune lineage deconvolution-to distinguish between NK and B cells. Again, this was done first separately for humans lungs and tumors, and later on, the aggregation of both sample types for NK cells only was subjected to further clustering analysis. Supplemental details on the rounds of analysis are provided in [an excel file](https://github.com/vijaybioinfo-org/NK_2021/preprocessing/supplements/CleanUpStepTracker.xlsx).<br>
For this section to be completed, we have developed a workflow that is described in detail in [another repo](https://github.com/VicenteFR/Seurat-based_scRNA-seq_Analysis_v2.3).<br>
The specific job scripts we wrote and ran are provided in separate folders under [jobs_scripts](https://github.com/vijaybioinfo-org/NK_2021/preprocessing/jobs_scripts).

### **Important disclaimers**
- The results will not be fully reproducible if the versions of the dependencies that were used during the analysis are not the same that you have installed to your machine.
  - Indeed, we don't expect the UMAP plots nor clustering results to be fully reproducible unless the exact same versions are used. Nonetheless, pretty similar results can be obtained that, importantly, reflect the biology that were captured during our analysis and we present in our manuscript.
  - For you to obtain similar results, you might need to tweak the contamination clusters that are removed during the process.
- The data we report on in our manuscript are huge! Therefore, we required tons of resources to be spent to process these data. Specifically, we ran these analyses in a cluster that allows for up to 300 GB in RAM (with as many processors as there might be available since it helps speed up the process). We haven't attempted to run this analysis in a local machine (and you shouldn't either).

---
# System requirements and installation
* R (version $\geq$ v3.6).
* Seurat v3.2.2
* MAST v1.12.
* Versions of other dependencies are specified in the out files under [jobs_scripts](https://github.com/vijaybioinfo-org/NK_2021/preprocessing/jobs_scripts).
* Please follow the instructions in [this repo](https://github.com/VicenteFR/Seurat-based_scRNA-seq_Analysis_v2.3) to install and set up the workflow. Do install the recommended dependencies for full reproducibility.

---
# Data requirements
To obtain the raw data for the results presented in our study, you must download it from our GEO submission (**to be provided**). The input for this section is, among other files, the following:
* The aggregated gene by cell UMI count matrix for each aggregated dataset either as the cellranger aggr/reanalyze output or the seurat object from a previous round of clustering, as indicated in [this excel file](https://github.com/vijaybioinfo-org/NK_2021/preprocessing/supplements/CleanUpStepTracker.xlsx).
* Cell-wise donor assignments—output from the HTO data preprocessing steps, as described in `./../hto_data` in this repo.