###########    -------   Get aggr cell subsets    --------    ###########

# Version.
# Full version: 0.1
# --> Version: 0
#   * No stable version reached yet.
# --> Subversion: 1
#   * First subversion.


### -------------------------- Description -------------------------- ###
# Script to obtain the aggr cell barcodes for the DICE tissue datasets.


### --------------------------- Libraries --------------------------- ###
library(Seurat)
library(data.table)
library(stringr)
# library(tidyr)


### --------------------------- Functions --------------------------- ###
# coll.version <- 0.5
gen.data.path <- '/path/to/user/R24/paper_developments/R24_Cancer/paper_items'
# coll.file <- paste0(gen.data.path, '/jobs_scripts/functions_collection.', coll.version, '.R')
# file.exists(coll.file)
# source(coll.file)


### ----------------------- General Arguments ----------------------- ###

# ---> General definitions.
# @ Seed.
set.seed(seed=1)
# @ Cell types.
cell.type.sets <- list(
    'B'='B', 
    'CD4'='CD4',
    'CD8'='CD8',
    'NK'='NK',
    'MYELOID'='Myeloid',
    'NKAndB'=c('NK', 'B')
)
cell.type.aggr.ids <- c(
    'B'='NKAndB', 
    'CD4'='CD4',
    'CD8'='CD8',
    'NK'='NKAndB',
    'MYELOID'='Myeloid'
)
tissue.types <- c('hlty'='Normal-lung', 'tumor'='Tumor-lung', 'all'='ALL')
# @ Cluster labels for each dataset.
obj.extended.names <- c(
  hlty.cd4='Healthy-CD4',
  hlty.cd8='Healthy-CD8',
  hlty.nk='Healthy-NK',
  hlty.b='Healthy-B',
  hlty.myeloid='Healthy-Myeloid',
# obj.extended.names <- c(
  tumor.cd4='Tumor-CD4',
  tumor.cd8='Tumor-CD8',
  tumor.nk='Tumor-NK',
  tumor.b='Tumor-B',
  tumor.myeloid='Tumor-Myeloid'
)
clust.labs <- c(
  hlty.cd4='RNA_snn_res.0.3',
  hlty.cd8='RNA_snn_res.0.2',
  hlty.nk='RNA_snn_res.0.1',
  hlty.b='RNA_snn_res.0.1',
  hlty.myeloid='RNA_snn_res.0.2',
# clust.labs <- c(
  tumor.cd4='RNA_snn_res.0.4',
  tumor.cd8='RNA_snn_res.0.2',
  tumor.nk='RNA_snn_res.0.3',
  tumor.b='RNA_snn_res.0.1',
  tumor.myeloid='RNA_snn_res.0.2'
)
# ---> Path definitions.
# Individual aggr table file paths
aggr.table.path <- '/path/to/user/sequencing_data/05-14-2022/aggr/data'
aggr.data.path <- '/mnt/hpcscratch/vfajardo/sequencing_data/05-14-2022/aggr'
# Reports path
reports.path <- paste0(gen.data.path, '/custom_aggrs')
if(!dir.exists(reports.path)) dir.create(reports.path)
# ---> File definitions.
# Seurat objects.
obj.files <- paste0(gen.data.path, '/SeuratObj_', obj.extended.names, '.RDS')
names(obj.files) <- names(obj.extended.names)
# Raw data directory.
# tenx.raw.data <- '/mnt/bioadhoc-temp/Groups/vd-vijay/vfajardo/sequencing_data/05-14-2022/aggr/R24_Cancer_Batches-1-to-20_CD4_Normal-lung/outs/filtered_feature_bc_matrix'
# ---> Check directories and files.
if(!all(dir.exists(gen.data.path))) stop(paste0('Following paths must be already defined:\n', gen.data.path, '\n', module.sigs.path))
reports.path <- gen.data.path
essential.files <- c(obj.files)
essential.files <- essential.files[!unlist(lapply(X=essential.files, FUN=file.exists))]
if(length(essential.files) > 0) stop(paste0('Next essential files are not appropriately defined -make sure you\'ve got their names right-:\n', paste0(essential.files, collapse='\n'), '\n'))


### ------------------------- Data Loading -------------------------- ###

# ---> Seurat objects' metadata.
meta.data <- lapply(X=obj.files, FUN=function(seurat.obj.file){
    seurat.obj <- readRDS(file=seurat.obj.file)
    tmp.data <- as.data.table(seurat.obj@meta.data)
    tmp.data[, barcode:=Cells(seurat.obj)]
    return(tmp.data)
})
meta.data <- rbindlist(l=meta.data, use.names=TRUE, fill=TRUE, idcol='data.set')
meta.data[, 
    tissue.type:=tissue.types[
        str_extract(string=data.set, pattern='hlty|tumor')
    ]
]
meta.data[, 
    cell.type:=toupper(
        str_replace(string=data.set, pattern='hlty\\.|tumor\\.', replacement='')
    )
]

# ---> Raw data directory.
# feature.info <- load.feature.info.file(tenx.raw.data=tenx.raw.data)

# ---> Define aggr tables.
# Define cell and tissue types.
aggr.pffx.lab <- 'R24_Cancer_Batches-1-to-20_'
aggr.tables <- lapply(X=names(cell.type.aggr.ids), FUN=function(cell.type){
    aggr.files <- paste0(aggr.table.path, '/', aggr.pffx.lab, cell.type.aggr.ids[cell.type], '_', tissue.types, '_aggr_table_annotated.1.0.csv')
    names(aggr.files) <- tissue.types
    tmp.data <- lapply(X=aggr.files, FUN=fread)
    return(tmp.data)
})
names(aggr.tables) <- names(cell.type.aggr.ids)


### ------------------------- Main program -------------------------- ###

# ---> Harmonize cell barcode between general and tissue-dependent aggregations (process per cell type).
meta.data[, gen.barcode:=barcode]
meta.data[, tissue.suffix:=str_replace(
    string=str_extract(string=barcode, pattern='-\\d+$'),
    pattern='-',
    replacement=''
)]
tmp.data <- lapply(X=names(cell.type.aggr.ids), FUN=function(cell.type){
    # Retrieve barcode suffixes.
    tmp.vals <- setdiff(x=tissue.types, y='ALL')
    tmp.data.1 <- lapply(X=tmp.vals, FUN=function(tmp.tissue){
        tmp.data.1 <- aggr.tables[[cell.type]][['ALL']][, .(
            sample.id=library_id,
            all.suffix=as.character(1:.N)
        )]
        tmp.data.2 <- aggr.tables[[cell.type]][[tmp.tissue]][, .(
            sample.id=library_id,
            tissue.suffix=as.character(1:.N)
        )]
        tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by='sample.id')
    })
    names(tmp.data.1) <- tmp.vals
    tmp.data.1 <- rbindlist(l=tmp.data.1, use.names=TRUE, idcol='tissue.type')
    setorderv(x=tmp.data.1, col='all.suffix')
    tmp.data.1[, cell.type:=cell.type]
    # Replace w/ updated suffixes.
    tmp.val <- cell.type
    tmp.data.2 <- meta.data[cell.type==tmp.val]
    tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by=c('tissue.suffix', 'tissue.type', 'cell.type'), all=TRUE, sort=FALSE)
    return(tmp.data)
})
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE)
tmp.check <- tmp.data[, .N]==meta.data[, .N]
if(!tmp.check) stop('Unexpected error. Same number of cells before and after barcode suffix replacement expected.\n')
tmp.data[, gen.barcode:=str_replace(string=gen.barcode, pattern=tissue.suffix, replacement=all.suffix)]
meta.data <- tmp.data

# ---> Output final cell of barcodes per cell type.
for(aggr in names(cell.type.sets)){
    data.sets <- cell.type.sets[[aggr]]
    set.pttn <- lapply(X=names(tissue.types), FUN=function(tissue.type) paste0(tissue.type, '.', tolower(data.sets)))
    set.pttn <- unlist(set.pttn)
    tmp.data <- meta.data[data.set %in% set.pttn]
    tmp.data <- tmp.data[, .(Barcode=gen.barcode)]
    tmp.file.name <- paste0(reports.path, '/', aggr, '_ALL.csv')
    fwrite(file=tmp.file.name, x=tmp.data)
}