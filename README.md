# IL-21 sustains cytotoxic activity in tissue-resident NK cells within tumors
------------

Description
------------

Repo for our manuscript entitled as above. This repository self-contains code as well as links to other repos with further code used to process and analyze the data produced in our study. We also provide the code written to generate the associated figures and tables.

### Summary
Adoptive natural killer (NK) cell-based therapies are less effective against solid tumors than hematological malignancies, largely because transferred NK cells have a limited capacity to infiltrate and persist within the tumor tissue microenvironment (TME). Potentially more persistent NK subtypes such as tissue-resident NK cells (trNK cells) and infection-induced adaptive NK cells have been identified. However, the spectrum of subtypes, and the mechanisms underlying their cytotoxic function and development within human tumors remain poorly understood, limiting the design of more effective NK cell-based therapies. Here, we profiled naturally occurring trNK cells from tumor and tumor-adjacent uninvolved lung tissue from >100 patients with early-stage lung cancer to identify mechanisms that support the generation of cytotoxic trNK cells within the TME. We identified distinct TCF7-expressing progenitor-like and adaptive trNK cell states, with adaptive trNK cells showing enhanced cytotoxicity and clonal expansion. Tumor-infiltrating cytotoxic trNK cells were enriched for gene signatures linked to IL-21 signaling, and spatial analyses revealed their close association with IL21-expressing T follicular helper cells. TGF-beta-mediated suppression of trNK cell cytotoxicity was reversed by IL-21 treatment ex vivo while preserving tissue-resident features, thereby enhancing the anti-tumor activity of trNK cells in vitro and in vivo. These findings establish IL-21 as an important driver of cytotoxic activity in trNK cells and provide a framework for developing NK cell therapies that combine durable tissue residency with enhanced anti-tumor cytotoxicity.

Major requirements
------------

For scRNA-seq data the project has been processed following modules/programs:  
• [R](https://cran.r-project.org/) (v4.1.3)  
• [Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/latest) (v3.1.0 or v7.1.0, per indicated)  
• [Scrublet](https://github.com/swolock/scrublet/blob/master/README.md) (v0.2.3)  
• [Seurat](https://satijalab.org/seurat/) (v4.3.0)  
• [vdjtools](https://vdjtools-doc.readthedocs.io/en/master/) (v1.2.1)  
• [Cellchat](https://github.com/jinworks/CellChat) (v2.2.0)

Regarding Spatial transcriptomics analysis, data has been processed following modules/programs:  
• [R](https://cran.r-project.org/) (v4.3.3)  
• [Xenium Ranger](https://www.10xgenomics.com/support/software/xenium-ranger/latest) (v4.0.1.1)  
• [Seurat](https://satijalab.org/seurat/) (v5.1.0)  
• [phenoptr](https://akoyabio.github.io/phenoptr/) (v0.3.2)

Raw data
------------

The raw and processed data files from our study can be downloaded from the following sources.
* **scRNA-seq data**: dbGaP, accession number: *TBD*
* **Bulk RNA-seq data**: GEO, accession number: *TBD*
* **Spatial transcriptomics data**: GEO, accession number: *TBD*

Raw data processing
------------

For more specific information about the data generation and processing, please check either the relevant folders in this repo or, if needed, the Methods section of our manuscript.  

> Relevant scripts and jobs scripts for both 5'v2 scRNA-seq and Spatial Transcriptomics data are all located in: `./pre-processing`

### Single-cell RNA-seq
* Preprocesing was performed with 10x's cellranger (v3.1.0) as generally described (see Methods section of publication).
* Quality control (QC), unbiased clustering, dimensionality reduction and cluster annotation were performed with the workflow detailed in [this repo](https://github.com/VicenteFR/Seurat-based_scRNA-seq_Analysis_v2.3)
> Relevant scripts are located in: `./pre-processing/transcriptome_data`

### Hashtag oligonucleotide (HTO) by sequencing data.
* Standard processing scripts in this repo.
> Relevant scripts are located in: `./pre-processing/hto_data`

### Spatial Transcriptomics
* To generate the clustering of Spatial Transcriptomics data just pull [our in-house pipeline](https://github.com/vijaybioinfo/Xenium-clustering/) using Seurat.  

### mtSNV pipeline
* To generate the results for the identification of informative mtSNVs, use [our in-house pipeline](https://github.com/vijaybioinfo/mtDNA_variant_detection).

### Bulk RNA-seq pipeline
* To generate the results for the Bulk RNA-seq, use [our in-house pipeline](https://github.com/vijaybioinfo/bulk_rnaseq_pipeline_2026.git).  


Figures & tables
------------
> Relevant scripts are located in: `./figures` and `./tables`

DGEA - You can follow [our in-house pipeline](https://github.com/vijaybioinfo/dgea)  
CellChat - You can use our in-house cellchat implementation.  


Usage & Citation
--------------
Run the following command to clone this repo:
```bash
git clone https://github.com/vijaybioinfo/trNK-2026.git
```
Please cite the following manuscript if you are using this repository:
> *Cite to be provided upon publication.*</br>
*Link to be provided upon publication*

Maintainers
-----------

Current maintainers:
* Kosuke Matsumura (kmatsumura@lji.org)
* Donaldo Sosa-Garcia (donaldo@lji.org)
* Ignacio Ramírez Bernabe (iramirez@lji.org)
* Vicente Fajardo Rosas (vfajardo@lji.org)

Vijayanand Lab.  
La Jolla Institute for Immunology La Jolla, CA 92037, USA

Contact
-----------
Please email any of the current maintainers.