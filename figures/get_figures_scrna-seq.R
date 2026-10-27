############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    -----   Preliminary figure stack    -----    ############
############    ------------ DICE NK project  -----------    ############
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############

## Author Kosuke Matsumura
# ---> About the script.
# Version: 1
# Subversion: 12_3
# Best module: R/4.2.2


############Description-----------------------------------------    ############
### -------------------------- Description -------------------------- ###
############Description-----------------------------------------    ############

# Script to get all bioinformatics-based panel figures for 
# any of the main and supplementary figures for the paper of the DICE trNK project.
# Specifically for the analysis of scRNA-seq data.
rm(list = ls())
gc()


############Libraries-----------------------------------------    ############
###--------------------------- Libraries --------------------------- ###
############libraries-----------------------------------------    ############

############Libraries-----------------------------------------    ############
###--------------------------- Libraries --------------------------- ###
############libraries-----------------------------------------    ############
library(Seurat)
library(data.table)
library(tidyr)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)
library(fgsea)
library(pheatmap)
library(GSVA)
# library(FSA)
# library(RColorBrewer)
# library(forcats)


############Functions-----------------------------------------    ############
### --------------------------- Functions --------------------------- ###
############Functions-----------------------------------------    ############
#Load handy fuctions
coll.version <- 0.3
gen.data.path <- '/home/kmatsumura/BioAdhoc/R24_NK/paper_development/job_scripts'
coll.file <- paste0(gen.data.path, '/functions_collection.', coll.version, '.R')
if(file.exists(coll.file)) source(coll.file) else stop('Faulty function collection file.\n')
total.workers <- parallel::detectCores(logical = FALSE)
#total.workers <- try(expr=as.integer(system(command='echo $PBS_NUM_PPN', intern=TRUE)))
if(is.na(total.workers)) total.workers <- try(expr=as.integer(system(command='echo $SLURM_JOB_CPUS_PER_NODE', intern=TRUE)))
if(!is.numeric(total.workers)) stop('Problem while defining number of processors/CPUs available per node.\n')


##Define FeaturePlot function
# Name: Get UMAP gene expression ggplot.
# Description:
# Provided a seurat object, the name of a file to save it to and any further pdf function arguments, this function will output the publication quality pdf of such a plot in both versions: 1) 'Complete' version: version with all labels; and 2) 'Blank' version: version to be used in the actual paper figures file.
# Arguments ------------------------------>
#
# Function:
get.umap.feature.gg <- function(
    obj.name,
    feature,
    cell.subset = NULL,
    assay = "RNA",
    slot = "data",
    dot.size = 1,
    hmap.center = 0,
    hmap.col.scale =   c("lightgrey",
                         "blue"),
    ordered = TRUE
){
  #convert the features to emsamble id
  tmp.feats <- feature.info[feature.info$name %in% feature,]
  tmp.feats <- tmp.feats$id
  
  obj <- srt.objs.list[[obj.name]]
  if (is.null(cell.subset)) {
    cell.subset <- Cells(obj)
  }
  # ---- fetch UMAP + expression ----
  meta.data <- as.data.table(cbind(
    obj@reductions$umap@cell.embeddings[cell.subset, ],
    expr = GetAssayData(
      obj,
      assay = assay,
      slot = slot
    )[tmp.feats, cell.subset] |> as.numeric()
  ))
  # ---- FeaturePlot-like ordering (important) ----
  if (ordered == TRUE){
    meta.data <- meta.data[order(meta.data[["expr"]])]
  }
  # ---- Color scale ----
  col.breaks <- seq(from=min(range(meta.data$expr, na.rm=TRUE)), to=max(range(meta.data$expr, na.rm=TRUE)), length.out=100)
  mid.point <- hmap.center
  mid.point <- which.min(abs(col.breaks - mid.point))
  #hmap.col.scale.1 <- colorRampPalette(c('#8FD4C8', '#9FDAD0', '#ABE0D6', '#BFE8E0', '#CEEBE5', '#DEF2EE', '#EAF8F5', '#FFFFFF'))(mid.point)
  #hmap.col.scale.2 <- colorRampPalette(c('#FFFFFF', '#F6EDF4', '#ECDEEC', '#E2C7E2', '#DBBADC', '#CFA5CF', '#C58FC6', '#BC7FBE'))(100-(mid.point+1))
  
  #hmap.col.scale.1 <- colorRampPalette(c('#8FD4C8', '#9FDAD0', '#ABE0D6', '#BFE8E0', '#CEEBE5', '#DEF2EE', '#EAF8F5', '#FFFFFF'))(mid.point)
  #hmap.col.scale.2 <- colorRampPalette(c('#FFFFFF', '#F6EDF4', '#ECDEEC', '#E2C7E2', '#DBBADC', '#CFA5CF', '#C58FC6', '#BC7FBE'))(100-(mid.point+1))
  
  # c("lightgrey", "#4134f7", "blue")
  
  #hmap.col.scale.1 <- colorRampPalette(
  #  c('#FFFFBF',
  #    '#FEE08B',
  #    '#FEC44F',
  #    '#FDAE61')
  #)(mid.point)
  
  # High 
  #hmap.col.scale.2 <- colorRampPalette(
  #  c('#FC8D59',
  #    '#D7301F',
  #    '#B22222',
  #    '#8B0000')
  #)(100 - (mid.point + 1))
  
  
  
  # ---- ggplot ----
  tmp.ggplot <- ggplot(data=meta.data, aes_string(x='UMAP_1', y='UMAP_2', color = "expr")) +
    geom_point(alpha=1, size=dot.size) +
    scale_color_gradientn(colors=hmap.col.scale) +
    scale_y_continuous(expand=expansion(add=0)) + scale_x_continuous(expand=expansion(add=0)) +
    labs(x='UMAP 1', y='UMAP 2', col='', fill=feature)
  return(tmp.ggplot)
}


##Define FGSEA function 
# Name: Do FGSEA
# Description:
# Given a seurat object and a tag listing a value of interest (defined in the object's metadata and provided to this function as 'tag.of.int'), this function will apply FGSEA comparing the signature enrichment between a set of values defined for the tag of interest (for each one independently) and the rest of the cells. Ranking metrics are calculated independently (see function above, 'get.rank.metrics') and same values are valid for this function.
# NOTE: If you wish to use this function for the development of other papers, pay attention to how files are being output here (see argument 'files.preffix').
# Arguments ------------------------------>
# metrics.src - A seurat object with features for which to calculate metrics. In its metadata the groups of interest are defined under the name provided through argument 'tag.of.int'.
# tag.of.int - Tag name (as defined in the seurat object's metadata) listing value of interest (to be taken as reference group).
# vals.to.depict - Character, lists the sets of values defined for the tag of interest for which an output file (GSEA plot) should be output.
# metric - Metric for ranking features to be calculated. Accepted metrics are: signal.to.noise, t.test, ratio.of.classes, diff.of.classes, log.ratio.of.classes (see function above, 'get.rank.metrics')
# files.preffix - Preffix to add to every plot to be output by this function.
# Function:

do.fgsea <- function(
    metrics.src, modules.feats, metric='signal.to.noise', cell.subset=NULL,
    tag.of.int, vals.to.depict=NULL, 
    output.path, files.preffix, full.reports=FALSE,
    do.heatmap=FALSE, signature.order=NULL, ref.order=NULL,
    col.metadata=NULL, ann.colors=ann.colors, just.get.res=FALSE
){
  # I added follwing three lines
  modules.feats <- lapply(X = modules.feats, FUN = function(tmp.feats){
    tmp.feats <- feature.info[feature.info$name %in% tmp.feats,]
    tmp.feats <- tmp.feats$id
    return(tmp.feats)
  })
  # ---> Predefinitions.
  files.preffix <- paste0(files.preffix, '_FGSEA-')
  # ---> Subset seurat object if necessary.
  if(!is.null(cell.subset)) metrics.src <- subset(x=metrics.src, cells=cell.subset)
  # ---> Metrics calculation for expresion values.
  # Define tag values.
  tag.values <- gtools::mixedsort(unique(metrics.src@meta.data[, tag.of.int]))
  tag.values <- tag.values[!is.na(tag.values)]
  if(is.factor(tag.values)) tag.values <- as.character(tag.values)
  # Confirm values to be depicted that were provided are indeed defined within the seurat object.
  if(is.null(vals.to.depict) | do.heatmap){
    vals.to.depict <- tag.values
  }else{
    to.check <- all(vals.to.depict %in% tag.values)
    if(!to.check){
      to.check <- vals.to.depict %in% tag.values
      to.check <- vals.to.depict[!to.check]
      tmp.err <- paste0('Next values were requested to be depicted, but they have not been previously defined for tag \'', tag.of.int, '\':\n',
                        paste0(to.check, collapse='\n'), '\nPlease make sure to provide valid values for the tag of interest.\n')
      stop(tmp.err)
    }
  }
  # Obtain metrics for each relevant tag value taken as a reference.
  val.metrics <- lapply(X=vals.to.depict, FUN=function(val.of.int) get.rank.metrics(seurat.obj=metrics.src, tag.of.int=tag.of.int, val.of.int=val.of.int, metric=metric))
  names(val.metrics) <- vals.to.depict
  # Remove values for which metrics calculation was not successful.
  to.check <- is.na(val.metrics)
  to.check <- names(to.check)[to.check]
  if(length(to.check)>0){
    tmp.warning <- paste('Metrics calculation was not successful for next values because their group size was too small:\n\t', paste0(to.check, collapse='\t\n'), '\nThey will be removed from further processes.\n')
    warning(tmp.warning)
    val.metrics <- val.metrics[!is.na(val.metrics)]
    to.check <- length(val.metrics)<1
    if(to.check){ tmp.error <- paste0('Seems that metrics calculation failed for all of the possible value from the tag of interest (', tag.of.int, '), so the program cannot continue any further.\n'); stop(tmp.error)}
  }
  # Save calculated metrics.
  if(full.reports){
    to.output <- rbindlist(l=lapply(X=val.metrics, FUN=function(x) data.table(gene=names(x), metric=x)), use.names=TRUE, idcol='tag.value')
    to.output <- spread(data=to.output, key=tag.value, value=metric)
    tmp.file.name <- paste0(output.path, '/', files.preffix, 'RankingMetricsPerTagValueAcrossGenes_Tag-', tag.of.int, '.csv')
    fwrite(file=tmp.file.name, x=to.output, na='NaN', quote=FALSE)
  }
  # ---> FGSEA
  cols.to.keep <- c('pathway', 'pval', 'padj', 'ES', 'NES', 'nMoreExtreme', 'size')
  fgsea.results <- lapply(X=vals.to.depict, FUN=function(tag.value){
    # FGSEA
    set.of.metrics <- val.metrics[[tag.value]]
    set.of.metrics <- set.of.metrics[!is.nan(set.of.metrics)]
    fgsea.results <- fgsea(pathways=modules.feats, stats=set.of.metrics, minSize=5, maxSize=500, nperm=10000)
    cat("FGSEA result size: ", nrow(fgsea.results), "\n")
    #fgsea.results <- fgsea(pathways=modules.feats, stats=set.of.metrics, minSize=10, maxSize=500, nproc=total.workers, nperm=10000)
    # Save plain results.
    fgsea.results <- fgsea.results[, ..cols.to.keep]
    if(full.reports){
      tmp.file.name <- paste0(output.path, '/', files.preffix, 'FGSEAResultsForTagValue_Tag-', tag.of.int, '_Val-', tag.value, '.csv')
      fwrite(file=tmp.file.name, x=fgsea.results, na=NA, quote=FALSE)
    }
    # Return general results.
    return(fgsea.results)
  })
  names(fgsea.results) <- vals.to.depict
  # ---> Just get results
  if(just.get.res) return(fgsea.results)
  # ---> Individual enrichment plots.
  if(do.heatmap){
    # Tidy data.
    tmp.data <- rbindlist(l=fgsea.results, use.names=TRUE, idcol='ref.tag')
    tmp.data <- tmp.data[, .(ref.tag, signature=pathway, p.adj=padj, NES)]
    tmp.data[p.adj>0.05, NES:=NA]
    tmp.data[, p.adj:=NULL]
    tmp.data <- spread(data=tmp.data, key=ref.tag, value=NES)
    tmp.data <- as.data.frame(tmp.data, stringsAsFactors=FALSE)
    row.names(tmp.data) <- tmp.data$signature; tmp.data$signature <- NULL
    # Set color scale and breaks for heatmap.
    col.breaks <- seq(from=min(range(tmp.data, na.rm=TRUE)), to=max(range(tmp.data, na.rm=TRUE)), length.out=100)
    mid.point <- which.min(abs(col.breaks - 0))
    hmap.col.scale.1 <- colorRampPalette(c('#8FD4C8', '#9FDAD0', '#ABE0D6', '#BFE8E0', '#CEEBE5', '#DEF2EE', '#EAF8F5', '#FFFFFF'))(mid.point)
    hmap.col.scale.2 <- colorRampPalette(c('#FFFFFF', '#F6EDF4', '#ECDEEC', '#E2C7E2', '#DBBADC', '#CFA5CF', '#C58FC6', '#BC7FBE'))(100-(mid.point+1))
    hmap.col.scale <- c(hmap.col.scale.1, hmap.col.scale.2)
    # Set order of data.
    if(!is.null(signature.order)){
      if(is.null(names(signature.order))) names(signature.order) <- signature.order
      tmp.data <- tmp.data[signature.order, ]
      signature.order <- setNames(object=names(signature.order), nm=signature.order)
      row.names(tmp.data) <- signature.order[row.names(tmp.data)]
    }
    if(!is.null(ref.order)){
      if(is.null(names(ref.order))) names(ref.order) <- ref.order
      tmp.data <- tmp.data[, ref.order]
      ref.order <- setNames(object=names(ref.order), nm=ref.order)
      colnames(tmp.data) <- ref.order[colnames(tmp.data)]
    }
    # @ Heatmap.
    # Plot 'complete' version
    tmp.width <- (ncol(tmp.data)*1) + 0.6
    tmp.height <- (nrow(tmp.data)*1) + 0.8
    tmp.file.name <- paste0(output.path, '/', files.preffix, 'SignaturesNESHeatmap.C.pdf')
    pheatmap(
      mat=tmp.data,
      color=hmap.col.scale, breaks=col.breaks, scale='none',
      na_col='#b3b3b3', border_color='black',
      cluster_rows=FALSE, cluster_cols=FALSE,
      annotation_col=col.metadata, annotation_names_col=FALSE, annotation_colors=ann.colors,
      show_rownames=TRUE, show_colnames=TRUE, legend=TRUE, annotation_legend=TRUE,
      filename=tmp.file.name, height=tmp.height, width=tmp.width+1.5
    )
    # Plot 'blank' version
    tmp.file.name <- paste0(output.path, '/', files.preffix, 'SignaturesNESHeatmap.B.pdf')
    pheatmap(
      mat=tmp.data,
      color=hmap.col.scale, breaks=col.breaks, scale='none',
      na_col='#b3b3b3', border_color='black',
      cluster_rows=FALSE, cluster_cols=FALSE,
      annotation_col=col.metadata, annotation_names_col=FALSE, annotation_colors=ann.colors,
      show_rownames=FALSE, show_colnames=FALSE, legend=FALSE, annotation_legend=FALSE,
      filename=tmp.file.name, height=tmp.height, width=tmp.width
    )
    # @ Color legend only.
    tmp.data <- as.data.frame(tmp.data)
    tmp.data$signature <- row.names(tmp.data)
    tmp.data <- gather(data=tmp.data, key='pop', value='nes', -`signature`)
    tmp.data <- tmp.data[!is.na(tmp.data$nes), ]
    # W/ labels
    tmp.ggplot <- ggplot(data=tmp.data, aes(x=nes, y=nes, col=nes)) + 
      geom_point() +
      scale_color_gradientn(colors=hmap.col.scale, breaks=scales::pretty_breaks(n=5)) +
      theme(
        legend.ticks=element_line(color='black', linewidth=0.6),
        legend.ticks.length=unit(0.22, "cm"),
        legend.frame=element_rect(color='black', linewidth=0.6)
      )
    tmp.ggplot <- get_legend(p=tmp.ggplot)
    tmp.file.name <- paste0(output.path, '/', files.preffix, 'SignaturesNESHeatmap.L1.pdf')
    pdf(file=tmp.file.name, height=2, width=2)
    print(as_ggplot(tmp.ggplot))
    dev.off()
    # W/out labels
    tmp.ggplot <- ggplot(data=as.data.frame(tmp.data), aes(x=nes, y=nes, col=nes)) + 
      geom_point() +
      scale_color_gradientn(colors=hmap.col.scale, breaks=scales::pretty_breaks(n=5), name=NULL, labels=NULL) +
      theme(
        legend.ticks=element_line(color='black', linewidth=0.6),
        legend.ticks.length=unit(0.22, "cm"),
        legend.frame=element_rect(color='black', linewidth=0.6)
      )
    tmp.ggplot <- get_legend(p=tmp.ggplot)
    tmp.file.name <- paste0(output.path, '/', files.preffix, 'SignaturesNESHeatmap.L2.pdf')
    pdf(file=tmp.file.name, height=1.25, width=0.3)
    print(as_ggplot(tmp.ggplot))
    dev.off()
  }else{
    for(tag.value in vals.to.depict){
      for(tmp.module in names(modules.feats)){
        # Check pathway was actually processed.
        tmp.check <- tmp.module %in% fgsea.results[[tag.value]][, pathway]
        if(!tmp.check) next
        # Retrieve stats, module signatures and further labels.
        tmp.lab <- toupper(str_replace_all(string=tmp.module, pattern='_', replacement=' '))
        tmp.pathway <- modules.feats[[tmp.module]]
        tmp.stats <- val.metrics[[tag.value]]
        tmp.stats <- tmp.stats[!is.nan(tmp.stats)]
        fgsea.stats <- fgsea.results[[tag.value]][pathway==tmp.module, round(x=c(ES, NES, padj), digits=6)]
        tmp.caption <- paste0('ES: ', fgsea.stats[1], '; NES: ', fgsea.stats[2], '; adj. P value: ', fgsea.stats[3])
        # Get plot.
        tmp.ggplot <- plot.fgsea.enrichment(pathway=tmp.pathway, stats=tmp.stats) + labs(title=tmp.lab, x='Rank', y='Enrichment score', caption=tmp.caption)
        # Output enrichment plot.
        tmp.file.name <- paste0(files.preffix, 'EnrichmentPlotForSignature-', tmp.module, '_Tag-', tag.of.int, '_Val-', tag.value)
        publish.plot(
          tmp.ggplot=tmp.ggplot, output.path=output.path, file.name=tmp.file.name, type='pdf', blank.comp=blank.complement.3, do.legend=FALSE, height=7, width=7
        )
      }
    }
  }
}

# ---------------------------------------->
# Name: Module scoring.
# Dependencies:
#   Seurat and stringr, ggplot2, stringr, data.table, R_handy_functions
# Functions dependencies:
#   translate.ids, among others.
# Description:
# You can find an extended explanation on the AddModuleScore function from Seurat in the script associated to this module function.
# This program takes a seurat object and a directory where a set of files listing gene signatures should be and calculates module score for each signature.
# Arguments ------------------------------>
# seurat.obj - seurat object.
# module.feats.dir - String for the absolute path to the directory storing one or more gene-signature lists. Each file should have a csv format (one column which must be named 'feature') with rows listing the set of features (usually gene names) that define the module and should be part of the seurat object; if they're not, they will be filtered out. Each file name must follow the next pattern '<Signature name -spaced substituted by underscores->_signature.csv
# signal.thold <- Signal threshold to pick positive and negative cells for the signal (i.e., to define module tag). Value ignored if flag FPR is activated (i.e., if it's a value different than NULL).
# is.ensembl - Logical indicating whether input seurat object has ENSEMBL IDs defined as gene IDs instead of common genen names.
# reports.path - Absolute path to where filtered gene signatures should be stored.
# Value ----------------------------->
# Seurat object with module scores calculated.
# Function:

get.module.scores <- function(seurat.obj, module.feats.dir, signal.thold, is.ensembl=FALSE, feature.info.path, reports.path=NULL){
  cat('\n\n')
  cat('############    -----------   Module Scoring    -----------    ############\n')
  cat('############    -Based on module-defining genes expression-    ############\n')
  
  ### --------------------------- Arguments --------------------------- ###
  this.color.scale <- c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000') # provided by Ciro.
  
  cat('\n\n')
  ### --------------------------- Load data --------------------------- ###
  cat('### --------------------------- Load data --------------------------- ###\n')
  # Module-defining genes file.
  feature.info <- read.csv(file = feature.info.path,  stringsAsFactors=FALSE)
  module.feats.files <- list.files(path=module.feats.dir, recursive=FALSE, include.dirs=FALSE, pattern='signature.csv$', full.names=TRUE)
  modules.feats <- lapply(X=module.feats.files, FUN=function(tmp.file){
    tmp.feats <- read.csv(file=tmp.file, stringsAsFactors=FALSE)
    if(!'feature' %in% colnames(tmp.feats)) stop(paste0('File ', tmp.file, ' not appropriately defined. Column listing the gene features should be named \'feature\'.\n\n'))
    tmp.feats <- tmp.feats$feature
    return(tmp.feats)
  })
  # I added follwing three lines
  modules.feats <- lapply(X = modules.feats, FUN = function(tmp.feats){
    tmp.feats <- feature.info[feature.info$name %in% tmp.feats,]
    tmp.feats <- tmp.feats$id
    return(tmp.feats)
  })
  
  names(modules.feats) <- list.files(path=module.feats.dir, recursive=FALSE, include.dirs=FALSE, pattern='signature.csv$', full.names=FALSE)
  names(modules.feats) <- str_replace(string=names(modules.feats), pattern='_signature.csv$', replacement='')
  names(modules.feats) <- str_replace_all(string=names(modules.feats), pattern='_', replacement='.')
  modules.names <- names(modules.feats)
  cat('Module-defining features loaded...\n')
  
  ### ----------------------- Module definition ----------------------- ###
  cat('### ----------------------- Module definition ----------------------- ###\n')
  data.path <- if(!is.null(reports.path)) paste0(reports.path, '/data') else NULL
  if(!is.null(data.path)) if(!dir.exists(data.path)) dir.create(data.path)
  
  # ---> Module-defining genes.
  # Make sure that all module-defining genes are part of the seurat object.
  cat('Filtering non-defined features...\n')
  modules.feats <- lapply(X=names(modules.feats), FUN=function(tmp.module){
    module.feats <- modules.feats[[tmp.module]]
    # Translate IDs if necessary.
    if(is.ensembl) module.feats <- translate.ids(ids=module.feats, ensembl=FALSE) else names(module.feats) <- module.feats
    # Define IDs as a data frame.
    module.feats <- module.feats[!is.na(names(module.feats))]
    module.feats <- data.frame(features=names(module.feats), original.features=module.feats, present=(names(module.feats) %in% rownames(seurat.obj)), stringsAsFactors=FALSE)
    module.feats <- module.feats[!is.na(module.feats$features), ]
    # Make sure there are any features present to continue.
    print(module.feats)
    tmp.check <- all(!module.feats$present)
    if(tmp.check){
      tmp.warning <- paste0('There are not enough features from module ', tmp.module, '. It will be skipped.\n')
      warning(tmp.warning)
      return(NA)
    }
    # Then, output final list of module-defining genes.
    if(!is.null(data.path)){
      tmp.file.name <- paste0(data.path, '/', str_replace_all(string=toupper(tmp.module), pattern='\\.', replacement='_'), '_ModuleDefiningGenesInExpressionData.csv')
      write.csv(x=module.feats, file=tmp.file.name, quote=FALSE, row.names=FALSE)
    }
    # Check at least 90% of the original siganture features are present in the seurat object.
    if(sum(module.feats$present)/nrow(module.feats) < 0.85) warning(paste0('\nFor module/signature ', tmp.module, ', we couldn\'t recover for the seurat object more than 85% of genes in the original list. Something may be wrong either with the list or the seurat object. Please be careful about that when interpreting the results.\n'))
    # Output final features.
    return(module.feats$features[module.feats$present])
  })
  modules.names <- modules.names[!is.na(modules.feats)]
  modules.feats <- modules.feats[!is.na(modules.feats)]
  names(modules.feats) <- modules.names
  cat('List of module-defining features depicting presence in the gene expression data output.\n')
  
  cat('\n\n')
  ### -------------------- Module score inference --------------------- ###
  cat('### -------------------- Module score inference --------------------- ###\n')
  # ---> Module score.
  # Add module score across cells.
  seurat.obj <- AddModuleScore(object=seurat.obj, features=modules.feats, name=modules.names)
  # Fix module names given by seurat by removing the suffix number.
  modules.names <- paste0(modules.names, '.score')
  names(modules.names) <- paste0(str_replace(string=modules.names, pattern='\\.score', replacement=''), 1:length(modules.names))
  for(tmp.module in names(modules.names)){
    colnames(seurat.obj@meta.data) <- str_replace(string=colnames(seurat.obj@meta.data), pattern=tmp.module, replacement=modules.names[tmp.module])
  }
  names(modules.names) <- str_replace(string=names(modules.names), pattern='\\d+$', replacement='')
  names(modules.names) <- str_replace_all(string=names(modules.names), pattern='\\.', replacement='_')
  
  return(seurat.obj)
}

# ---------------------------------------->
#Define the function again
# Name: Dot plot.
# Of note: This function's idea was based on Seurat's output dot plot, but it was developed independently up to the previous version (see past function's file version). The size scale function implementation was copied from Seurat's function.
# Updates:
# @     Compared to the function in version 1.4, this one:
# 1. If required, applies a row-wise normalization by z-scoring.
# 2. Can set boundaries for both scales, color and size, through arguments col.min/col.max and size.min/size.max.
# 3. If required, sets a hot and cold color scale taking into account the minimum and maximum color scale values by centering at zero. This won't work if there's only negative or only positive values. (This needs further revising).
# 4. Sets a scale by taking into consideration either the size or the radius of the points in the plot. Original seurat's function does apply a radius scaling and it seems to be the best option for this kind of plots.
# Description:
# Given a seurat object, this function will:
# 1. Fetch the appropriate data columns from it.
# 2. If required (groups.of.int not NULL), keep the original identity of the groups indicated only, setting the others' identity to 'rest'
# 3. If indicated (filter.tag set to a -valid- value other than NULL), filter out or forward (keep set to F or T, respectively) a set of groups (groups.to.filter) of a tag (filter.tag).
# 4. Calculate expression mean and burst frequency for each gene across groups.
# 5. Normalize the mean expression values across the groups of interest (i.e., per gene or in a row-wise manner).
# 6. Set color and size scales by setting boundaries and plugging the appropriate arguments.
# 7. Generate a dot plot (based on seurat's style) as a ggplot, which will be either returned (when file.name set to NULL) or output (when a valid file.name provided).
# Arguments ------------------------->
#     @ About data.
#         Data specificities.
# seurat.obj - seurat object. The tags indicated should already be defined. No default.
# features - Features for which to output expression values' maetrics. All should be defined as part of a valid data slot (choose slot accordingly if any prefered -'data' as default-). No default.
# slot - Valid data slot to fetch data from. Valid values are the usual ones: counts, data or scale.data. Default: data.
# do.norm - Logical, indicates whether a row-wise normalization should be applied through z-scoring. Default: FALSE.
# ensembl - Logical, defines whether input features' vector (defined through argument 'fetaures') lists ENSEMBL IDs (then expecting the gene names to be defined as the vector's names). Default: FALSE.
#         Data filtering and grouping.
# groups.tag - Tag defined in the seurat object metadata whose groups will be the basis for the x axis of the product. No default.
# groups.of.int - Character listing the set of tag groups whose identity should be kept as originally defined, otherwise they'll be set to 'rest'. This step is skipped is it's set to NULL. Default: NULL.
# groups.order - Order to show the groups in the x axis. When an order is provided, all groups defined in the object (post filtering, if any) are expected to be in the order vector and the other way around. If none provided (i.e., NULL value), the function will attempt to sort them based on their names (in an alphanumeric order). Default: NULL
# filter.tags - Chracter vector, which should list the tag or tags in the seurat object metadata whose groups will be used to filter out cells. When filtering should be applied based on multiple tags, their groups to be considered to keep or exclude cells must be passed as chracter vectors embedded within a list; then, this character vector must have the same length as that of the list provided through the argument 'groups.to.filter'. No cell filtering is applied when it's set to NULL. Default: NULL.
# groups.to.filter - List of vectors for the set of tag groups that'll be kept ot removed according to the argument 'keep' value. Must have the same length as the vector provided in 'filter.tags'. This step is skipped is it's set to NULL. Default: NULL.
# keep - Logical indicating if groups defined in 'groups.to.filter' should be filter out or forward when set to FALSE and TRUE, respectively. Default: TRUE
# na.rm - Logical, indicates whether NA values should be discarded for the tag listing groups. Default: TRUE.
# feature.thold - Numeric, must be a number put as a threshold to filter forward only certain cells.
#     @ About plot and output.
#         Color scale.
# this.color.scale - Character, indicates color scale to use to represent average expression. If set to NULL, a color scale that's broadly used in Vijay's lab will be used instead. If a length-1 vector that equals 'hot.and.cold', 'hot.n.cold' or 'cold.and.hot' is input, a hot and cold color scale is infered by taking into account the range of average expression values seen in the data (post normalization if required).
# col.min - Numeric, lower boundary for color scale. Must be smaller than col.max (if provided) to input an appropriate color range. If NULL, no boundary is set. Default: NULL.
# col.max - Numeric, upper boundary for color scale. Must be larger than col.min (if provided) to input an appropriate color range. If NULL, no boundary is set. Default: NULL.
#         Color scale.
# scale.by - Character, indicates whether a radius scaling or a standard size scaling should be applied. Accepted values are 'radius' and 'size', respectively. Default: 'radius'.
# dot.scale - Character, second value in the argument 'range' of 'scale_radius' or 'scale_size'. Default: 6.
# size.min - Numeric, lower boundary for size scale. Must be smaller than size.max (if provided) to input an appropriate color range. If NULL, no boundary is set. Default: NULL.
# size.max - Numeric, upper boundary for size scale. Must be larger than col.min (if provided) to input an appropriate color range. If NULL, no boundary is set. Default: NULL.
#         Output.
# file.name - Though NULL indicates that the ggplot should be returned instead, this value indicates the name of the file to output plot to. Default: NULL.
# Value ----------------------------->
# Either a ggplot object or NULL value and the plot is output as a file (when file.name is different than NULL)
# Function:

dot.plot <- function(
    # @   About data.
  # Data specificities.
  seurat.obj, features, slot='data', do.norm=FALSE, ensembl=FALSE,
  # Data filtering and grouping.
  groups.tag, groups.order=NULL, groups.of.int=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=TRUE, na.rm=TRUE, feature.thold=NULL,
  # @   About plot and output.
  # Color scale
  this.color.scale=NULL, col.min=NULL, col.max=NULL,
  # Size scale.
  scale.by='radius', dot.scale=6, size.min=NA, size.max=NA,
  # File.
  file.name=NULL,
  #Labels
  named.labels=NULL
){
  # ---> Arguments and checks.
  # Set color scale.
  if(is.null(this.color.scale)){
    this.color.scale <- c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
  }else{
    if(length(this.color.scale)==1){
      if(this.color.scale=='hot.and.cold' | this.color.scale=='cold.and.hot' | this.color.scale=='cold.n.hot') this.color.scale <- NA
    }
  }
  # Set size scale.
  scale.func <- switch(
    EXPR=scale.by,
    size=scale_size,
    radius=scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  # Features must be unambigously defined.
  features <- unique(features)
  # Define CPMs-related arguments.
  if(is.null(slot)){
    cpm <- FALSE
  }else{
    if(slot == 'cpm'){
      slot <- 'counts'
      cpm <- TRUE
    }else{
      cpm <- FALSE
    }
  }
  # Filter tags and groups to filter dimensions must agree with each other.
  if(!is.null(filter.tags)){
    if(length(filter.tags)>1){
      if(!is.list(groups.to.filter) | length(filter.tags)!=length(groups.to.filter)) stop(paste0('When multiple tags are provided for cell filtering, the argument groups.to.filter must be a list with the same dimension length.\n'))
    }else{
      if(!is.list(groups.to.filter)){
        if(!is.vector(groups.to.filter)) stop('Argument "groups.to.filter" must be either a list (with the same length as that for the argument "filter.tags") or a vector.\n') else groups.to.filter <- list(groups.to.filter)
      }
    }
    names(groups.to.filter) <- paste0('filter.tag.', 1:length(filter.tags))
  }
  # ---> Fetch data.
  # Look for feature.
  # First, look for it in the metadata if slot is set to NULL.
  if(is.null(slot)){
    if(!all(features %in% colnames(seurat.obj@meta.data))){
      features <- features[!features %in% colnames(seurat.obj@meta.data)]
      stop(paste0('Following features could not be found in the seurat object\'s metadata:\n', paste0(features, collapse='\n'), '\n'))
    }
    feature.vals <- data.table(seurat.obj@meta.data[, features])
  }else{
    # If set to a value different than NULL, make sure it's a valid slot and try to find the feature there.
    if(!(slot %in% c('counts', 'data', 'scale.data'))) stop('Not valid slot. Should be anyone among "counts", "data", "scale.data".\n')
    if(!all(features %in% rownames(seurat.obj))){
      features <- features[!features %in% rownames(seurat.obj)]
      stop(paste0('Following features could not be found in the seurat object\'s feature set:\n', paste0(features, collapse='\n'), '\n'))
    }
    # Proceed if all good, though process branches based on CPM possible requirement.
    if(cpm){
      feature.vals <- RelativeCounts(data=GetAssayData(object=seurat.obj, slot=slot, assay='RNA'), scale=1e6, verbose=FALSE)
      feature.vals <- log2(feature.vals[features, ] + 1)
      feature.vals <- as.data.table(t(feature.vals))
    }else{
      feature.vals <- as.data.table(FetchData(object=seurat.obj, vars=features, slot=slot))
    }
  }
  # Fetch the data for the rest of the tags after checking everything exists.
  if(!(groups.tag %in% colnames(seurat.obj@meta.data))) stop('Tag for groups of interest not defined in seurat object\'s metadata.\n')
  if(!is.null(filter.tags)){
    tmp.assess <- filter.tags %in% colnames(seurat.obj@meta.data)
    if(!all(tmp.assess)){
      tmp.assess <- filter.tags[!tmp.assess]
      stop(paste0('Tags for filtering set to a value different than NULL, but some of them not defined in seurat object\'s metadata, as follows.\n', paste0(tmp.assess, collapse='\n'), '\n'))
    }
    tmp.data <- data.table(seurat.obj@meta.data[, filter.tags])
    colnames(tmp.data) <- paste0('filter.tag.', 1:ncol(tmp.data))
    tmp.data <- cbind(data.table(groups.tag=seurat.obj@meta.data[, groups.tag]), tmp.data)
  }else{
    tmp.data <- data.table(groups.tag=seurat.obj@meta.data[, groups.tag])
  }
  # Join altogether.
  tmp.data <- cbind(tmp.data, feature.vals)
  # ---> Groups tag.
  # Filter out non-desired cells according to groups.tag
  if(na.rm) tmp.data <- tmp.data[!is.na(groups.tag)]
  # Set rest group for the groups tag if necessary.
  if(!is.null(groups.of.int)){
    tmp.data[!groups.tag %in% groups.of.int, groups.tag:='Rest']
    # Similarly, merge groups of interest if their label is provided.
    if(!is.null(groups.label)){
      tmp.data[groups.tag %in% groups.of.int, groups.tag:=as.character(groups.label)]
      # Set factors leaving 'Rest' label right at the end of the levels' order.
      new.lvls <- c(as.character(groups.label), 'Rest')
      tmp.data$groups.tag <- factor(x=as.character(tmp.data$groups.tag), levels=new.lvls)
    }
  }
  # ---> Filter out non-desired cells.
  if(!is.null(filter.tags)){
    # Filter all cells with not available data regarding all column-to-filter information.
    if(na.rm){
      for(idx in 1:length(filter.tags)){
        tmp.col <- paste0('filter.tag.', idx)
        tmp.data <- tmp.data[!is.na(get(tmp.col))]
      }
    }
    # Then, when requested, keep or remove input values.
    if(!is.null(groups.to.filter)){
      for(idx in 1:length(filter.tags)){
        tmp.col <- paste0('filter.tag.', idx)
        if(keep) tmp.data <- tmp.data[get(tmp.col) %in% groups.to.filter[[tmp.col]]] else tmp.data <- tmp.data[!get(tmp.col) %in% groups.to.filter[[tmp.col]]]
      }
    }
  }
  # ---> Tidy data up!
  if(!is.null(filter.tags)){
    tmp.cols <- colnames(tmp.data)[!str_detect(string=colnames(tmp.data), pattern='^filter.tag.\\d+$')]
    tmp.data <- tmp.data[, ..tmp.cols]
  }
  tmp.data <- as.data.table(gather(data=tmp.data, key='feature', value='exp.val', -groups.tag))
  # ---> Apply feature threshold.
  if(!is.null(feature.thold)){
    if(!is.numeric(feature.thold)) stop('Not valid feature threshold provided. Must be a numeric value.\n')
    tmp.data <- tmp.data[exp.val > feature.thold]
  }
  # ---> Values for color and size scales.
  tmp.data <- tmp.data[, .(exp.mean=mean(exp.val, na.rm=TRUE), burst.freq=.SD[exp.val > 0, .N]*100/.N), by=.(groups.tag, feature)]
  # ---> Row-wise normalization.
  # Apply only if required.
  if(do.norm){
    tmp.data.2 <- as.data.frame(spread(data=tmp.data[, .(groups.tag, feature, exp.mean)], key=groups.tag, value=exp.mean))
    row.names(tmp.data.2) <- tmp.data.2$feature; tmp.data.2$feature <- NULL
    tmp.data.2 <- as.data.frame(t(scale(x=t(tmp.data.2))))
    tmp.data.2$feature <- row.names(tmp.data.2); tmp.data.2 <- as.data.table(tmp.data.2)
    tmp.data.2 <- gather(data=tmp.data.2, key='groups.tag', value='exp.mean', -feature)
    tmp.data <- merge(x=tmp.data[, .(groups.tag, feature, burst.freq)], y=tmp.data.2, by=c('feature', 'groups.tag'))
  }
  # ---> Set boundaries for scales and get final color definitions.
  # Color scale.
  if(!is.null(col.min) & !is.null(col.max)) if(col.min > col.max) stop('Boundaries for color scale were provided, but the minimum value to be set is larger than the maximum one. Provide a valid range if any is required.\n')
  if(!is.null(col.min)) tmp.data[exp.mean < col.min, exp.mean:=col.min]
  if(!is.null(col.max)) tmp.data[exp.mean > col.max, exp.mean:=col.max]
  # Warm and cold color scale is only available when row-wise normalization is applied
  if(do.norm & all(is.na(this.color.scale))){
    col.breaks <- seq(from=min(range(tmp.data$exp.mean)), to=max(range(tmp.data$exp.mean)), length.out=100)
    mid.point <- which.min(abs(col.breaks - 0))
    this.color.scale.1 <- colorRampPalette(c('blue', 'mediumblue', 'white'))(mid.point)
    this.color.scale.2 <- colorRampPalette(c('white', 'yellow', 'gold'))(100-(mid.point+1))
    this.color.scale <- c(this.color.scale.1, this.color.scale.2)
  }
  # Size scale.
  if(!is.na(size.min) & !is.na(size.max)) if(size.min > col.max) stop('Boundaries for size scale were provided, but the minimum value to be set is larger than the maximum one. Provide a valid range if any is required.\n')
  if(!is.na(size.min)) tmp.data[burst.freq<size.min, burst.freq:=size.min]
  if(!is.na(size.max)) tmp.data[burst.freq>size.max, burst.freq:=size.max]
  # ---> Set factors.
  # For groups.
  if(is.null(groups.order)){
    new.lvls <- mixedsort(x=unique(as.character(tmp.data$groups.tag)))
    tmp.data$groups.tag <- factor(x=as.character(tmp.data$groups.tag), levels=new.lvls)
  }else{
    if(!(all(groups.order %in% unique(as.character(tmp.data$groups.tag))) & all(unique(as.character(tmp.data$groups.tag)) %in% groups.order))) stop('When an order is provided, all groups defined in the object (post filtering, if any) are expected to be in the order vector and the other way around.\n')
    tmp.data$groups.tag <- factor(x=as.character(tmp.data$groups.tag), levels=as.character(groups.order))
  }
  # For features (same order as the one provided)
  new.lvls <- if(ensembl) names(features) else features
  tmp.data$feature <- factor(x=tmp.data$feature, levels=rev(new.lvls))
  # ---> Dot plot.
  tmp.ggplot <- ggplot(data=tmp.data, aes(x=groups.tag, y=feature, size=burst.freq, fill=exp.mean)) +
    geom_point(stroke=0.7, shape=21, col='black') +
    # geom_point() +
    scale_fill_gradientn(colors=this.color.scale) +
    scale.func(range=c(0, dot.scale), limits=c(size.min, size.max)) +
    {if (!is.null(named.labels)) scale_y_discrete(labels = named.labels)} +
    labs(x='Groups', y='Gene', fill='Average expression', size='Percent expressed') +
    theme(panel.background=element_blank(), panel.border=element_blank(), axis.line=element_line())
  # ---> Return.
  # If reports path is set to null, return ggplot.
  if(is.null(file.name)){
    return(tmp.ggplot)
  }else{
    pdf(file=file.name)
    print(tmp.ggplot)
    dev.off()
    return(NA)
  }
}


# ---------------------------------------->
#Define the function for GSVA
# Name: get.exp.data, process_cluster_condition
# 
####################GSVA function
# ── 1) Get exp data ─────────────────────────────────────
get.exp.data <- function(seurat.obj,
                         tag.of.interest,
                         group,
                         norm.data = "NOT-cpms") {
  
  exp.data <- if (norm.data == "cpms") {
    RelativeCounts(seurat.obj@assays$RNA@counts,
                   scale   = 1e6,
                   verbose = TRUE)
  } else {
    GetAssayData(seurat.obj, slot = "data", assay = "RNA")
  }
  
  grp.cells <- Cells(seurat.obj)[
    seurat.obj@meta.data[[tag.of.interest]] == group &
      !is.na(seurat.obj@meta.data[[tag.of.interest]])
  ]
  
  if (length(grp.cells) < 10) return(NULL)
  exp.data[, grp.cells, drop = FALSE]
}

# ── 2) Make mean table of donor and tissue ────────────────
process_cluster_condition <- function(seurat.obj,
                                      cluster.id,
                                      tissue.condition,
                                      cluster.col,
                                      donor.tag   = "full.donor.id.tag",
                                      norm.data   = "NOT-cpms") {
  
  cell.use <- Cells(seurat.obj)[
    seurat.obj@meta.data[[cluster.col]] == cluster.id &
      seurat.obj@meta.data$tissue.tag   == tissue.condition
  ]
  if (length(cell.use) == 0) return(NULL)
  
  sub.obj <- subset(seurat.obj, cells = cell.use)
  
  donors <- gtools::mixedsort(unique(na.omit(
    as.character(sub.obj@meta.data[[donor.tag]])
  )))
  if (length(donors) == 0) return(NULL)
  
  gene.stats <- lapply(donors, function(d) {
    exp <- get.exp.data(sub.obj, donor.tag, d, norm.data)
    if (is.null(exp)) return(NULL)
    data.table::data.table(
      donor = d,
      t(Matrix::rowMeans(exp))       # transpose
    )
  })
  gene.stats <- data.table::rbindlist(gene.stats, fill = TRUE)
  
  gene.stats   # wide format（row:donor、col:gene）
}

# ---------------------------------------->
#Define the function for UMAP
# Name: get.exp.data, process_cluster_condition
# 
####################UMAP function
#Make mean table of donor and tissue ────────────────
process_condition <- function(seurat.obj,
                              tissue.condition,
                              donor.tag   = "full.donor.id.tag",
                              norm.data   = "NOT-cpms") {
  
  cell.use <- Cells(seurat.obj)[
    seurat.obj@meta.data$tissue.tag   == tissue.condition
  ]
  if (length(cell.use) == 0) return(NULL)
  
  sub.obj <- subset(seurat.obj, cells = cell.use)
  
  donors <- gtools::mixedsort(unique(na.omit(
    as.character(sub.obj@meta.data[[donor.tag]])
  )))
  if (length(donors) == 0) return(NULL)
  
  gene.stats <- lapply(donors, function(d) {
    exp <- get.exp.data(sub.obj, donor.tag, d, norm.data)
    if (is.null(exp)) return(NULL)
    data.table::data.table(
      donor = d,
      t(Matrix::rowMeans(exp))  #Transpose
    )
  })
  gene.stats <- data.table::rbindlist(gene.stats, fill = TRUE)
  
  gene.stats   # wide format（row:donor、col:gene）
}



############General Arguments-----------------------------------------    ############
### ----------------------- General Arguments ----------------------- ###
############General Arguments-----------------------------------------    ############
# ---> General definitions.
# @ Paper.
main.prj <- 'DICE'
this.prj <- 'NK'
this.figure <- 'figure_stack_v12_4'
# @ Seed.
set.seed(seed=1)
# @ Dataset labels.
obj.extended.names <- c(
  CD16neg = 'SeuratObjectForPrjR24_Cancer_Batches-1-to-20_NK_ALL_CD16neg-2_WithArgs_NoPCs_30',
  CD16pos = 'SeuratObjectForPrjR24_Cancer_Batches-1-to-20_NK_ALL_CD16pos-2_WithArgs_NoPCs_30',
  trNK = 'SeuratObjectForPrjR24_Cancer_Batches-1-to-20_NK_ALL_CD16neg-2-TR-1_WithArgs_NoPCs_30',
  aggre = 'SeuratObjectForPrjR24_Cancer_Batches-1-to-20_NK_ALL_WithArgs_NoPCs_30'
)
gen.cell.types <- c(
  CD16neg ='NK',
  CD16pos ='NK' ,
  trNK  ='NK' ,
  aggre = 'NK'
)
# @ Cluster labels for each dataset.
clust.labs <- c(
  CD16neg='RNA_snn_res.0.6',
  CD16pos='RNA_snn_res.0.6',
  trNK = 'RNA_snn_res.0.4',
  aggre = 'RNA_snn_res.0.7'
)
# @ For tag analysis.
# min.no.cells.per.group <- 50 # @ Internal parameter for tag specific analysis (part of figure ?).
# colors.table <- NULL
# ---> Cluster identities
clusts.defs <- list(
  `CD16neg`=c(
    `0`='trNK',
    `3`='CD56_bright',
    `4`='CD160',
    `6`='ILC3'    
  ),
  `CD16pos`=c(
    `0`='CX3CR1',
    `2`='EGR1',
    `3`='AREG',
    `4`='Adaptive',
    `5`='TNFRSF9'
  ),
  `trNK`=c(
    `0`='Adaptive_trNK',
    `1`='ENTPD1',
    `2`='Cytotoxic',
    `3`='CREM',
    `4`='TCF7'
  ),
  `aggre`=c(
    `0`='CD16pos',
    `3`='CD16neg'
  )
)
# @ Color per cluster per dataset.
# See definitions at: https://docs.google.com/spreadsheets/d/1_lNI4M7P-LBuwfQWEgag4pHaycMZMNcmzGv_k_e7x5U/edit#gid=0
clusts.cols <- list(
  `CD16neg`=c(
    `0`='#DF594C', # trNK
    `3`='#F5C04A', # CD56_bright
    `4`='#006A8E', # CD160
    `6`='#4CAF50' # ILC3
  ),
  `CD16pos`=c(
    `0`='#8B5FE0', # CX3CR1 hi
    `2`='#DC5C9F', # EGR1
    `3`='#4E9FD1', # AREG
    `4`='#A985D8', # Adaptive
    `5`='#63B76C' # TNFRSF9
  ),
  `trNK`=c(
    `0`="#7BC8A4", 
    `1`="#A78BD3", 
    `2`="#F08A72", 
    `3`="#F2C14E", 
    `4`="#6BAED6" 
  ),
  `aggre`=c(
    `0`='#5D72E0', 
    `3`='#D4824B'
  )
)
# @ Other terms' colors.
signatures.col.scale <- c('#ffffff', '#ffffe0', '#ffffad', '#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
cont.col.scale.1 <- c('#FFFFFF', '#414487', '#2A788E', '#22A884', '#FDE725', '#FFEA00')
#Number of the cells for subseting
subset.n <- 50000
# ---> Path definitions.
# General data and reports paths.
prj.gen.path <- paste0('/home/kmatsumura/BioAdhoc/R24_NK/paper_item')
gen.data.path <- paste0('/home/kmatsumura/BioAdhoc/R24_NK/paper_item/seuratobject')
# final.figs.path <- paste0('/home/kmatsumura/BioAdhoc/R24_NK/paper_development/fig')
final.figs.path <- paste0('/path/to/user/R24/paper_developments/DICE-Cancer-NK/final_figures')
this.fig.path <- paste0(final.figs.path, '/', this.figure)
if(!dir.exists(this.fig.path)) dir.create(this.fig.path)
reports.path <- paste0(this.fig.path, '/', this.figure, '_panels')
if(!dir.exists(reports.path)) dir.create(reports.path)
# ---> File definitions.
# DICE cohort seurat objects.
objs.files.list <- paste0('/home/kmatsumura/BioAdhoc/R24_NK/paper_item/seuratobject/',obj.extended.names,'.RDS')
names(objs.files.list) <- names(obj.extended.names)
# Gene sets' files (for module scoring & FGSE analyses).
module.feats.dir <- paste0(prj.gen.path, '/signatures')
# Feature Info file (for conversion gene name to emsamble ID).
feature.info.path <- paste0(prj.gen.path, '/FeaturesInfo.csv')
# CD4 expression
CD4.exp.path <-paste0(prj.gen.path, '/CD4_tumor_exp.RDS')
# Patient metadata
metadata.path <-paste0(prj.gen.path, '/metadata.csv')
# Count table path
CD16neg.data.path <- paste0(prj.gen.path, '/CD16neg_tumor_raw.csv')
CD16pos.data.path <- paste0(prj.gen.path, '/CD16pos_tumor_raw.csv')
CD45.data.path <- paste0(prj.gen.path, '/DICE CD45 Full Summary.csv')
# TIL data 
#TIL.info.path <- paste0(prj.gen.path, '/TILSummary.csv')
#CD8_CD3_cell_density
CD8.data.path <- paste0(prj.gen.path, '/DICETumour_CD3CD8meanSE.csv')
#MHCpos
MHC.data.path <- paste0(prj.gen.path, '/trNKMHCExpressionPerPatient.csv')
#facs data
facs.data.path <- paste0(prj.gen.path, '/FACSCellCountInfo.csv')


############Data Loading    -----------------------------------------    ############
### ------------------------- Data Loading -------------------------- ###
############Data Loading    -----------------------------------------    ############
# Main seurat objects.
srt.objs.list <- lapply(X=objs.files.list, FUN=readRDS)
names(srt.objs.list) <- names(objs.files.list)

#Module files
module.feats.files <- list.files(path=paste0(module.feats.dir), recursive=FALSE, include.dirs=FALSE, pattern='signature.csv$', full.names=TRUE)
modules.feats <- lapply(X=module.feats.files, FUN=function(tmp.file){
  tmp.feats <- read.csv(file=tmp.file, stringsAsFactors=FALSE)
  if(!'feature' %in% colnames(tmp.feats)) stop(paste0('File ', tmp.file, ' not appropriately defined. Column listing the gene features should be named \'feature\'.\n\n'))
  tmp.feats <- tmp.feats$feature
  # Replace underscores with points in the gene names just as Seurat does when one creates a seurat object.
  tmp.check <- any(str_detect(string=tmp.feats, pattern='_'))
  if(tmp.check){
    tmp.warning <- paste0('Feature names cannot have underscores (\'_\'), replacing with dashes (\'-\'). This for file: ', basename(tmp.file))
    warning(tmp.warning)
    tmp.feats <- str_replace_all(string=tmp.feats, pattern='_', replacement='-')
  }
  # Return.
  return(tmp.feats)
})
names(modules.feats) <- list.files(path=paste0(module.feats.dir), recursive=FALSE, include.dirs=FALSE, pattern='signature.csv$', full.names=FALSE)
names(modules.feats) <- str_replace(string=names(modules.feats), pattern='_signature.csv$', replacement='')
names(modules.feats) <- str_replace_all(string=names(modules.feats), pattern='_', replacement='.')
modules.names <- names(modules.feats)

# Expression data of CD4
CD4.exp <- readRDS(CD4.exp.path)

# load feature info
feature.info <- read.csv(file = feature.info.path,  stringsAsFactors=FALSE)

# Patient metadata
metadata <- read.csv(metadata.path)
colnames(metadata) <- c('Sample', 'donor.id.tag', 'age', 'Gender', 'Staging')

# Count table
# Count table path
CD16neg.data <- read.csv(CD16neg.data.path)
CD16neg.data[is.na(CD16neg.data)] <- 0
setDT(CD16neg.data)
CD16neg.data[, trNK := tr_adaptive + unnamed + CD16neg_TNFRSF9 + tr_CD39]

CD45.data<- read.csv(CD45.data.path)
facs.data <- read.csv(facs.data.path)
#CD8_CD3_cell_density
CD8.data <- read.csv(CD8.data.path)
#MHCpos
MHC.data<- read.csv(MHC.data.path)

############Data preprocessing-----------------------------------------    ############
### ---------------------- Data preprocessing ----------------------- ###
############Data preprocessing-----------------------------------------    ############
# ---> Seurat objects, basic processing
# Combination instructions.
comb.insts <- list(
  `CD16neg`=list(
    `0`=c('0', '1', '2', '10'),
    `3`=c('3', '5'),
    `6`=c('6', '7', '9')
  ),
  `CD16pos`=list(
    `0`=c('0', '1')
  ),
  `trNK`=list(
    `2`=c('2', '5')
  ),
  `aggre`=list(
    # NOTE VFR: Important change. We didn't consider clusters 12, 13 and 14, since they represent contaminating lineages like B cells (12) or T cells (13). 14 is simply too small.
    # `0`=c('0', '1','2', '4', '9', '10', '11', '12', '14'),
    `0`=c('0', '1', '2', '4', '9', '10'),
    # `3`=c('3', '5','6', '7', '8', '13')
    `3`=c('3', '5', '6', '7', '8')
  )
)

# Process per seurat object
for(tmp.lab in names(srt.objs.list)){
  # Combine clusters (if necessary)
  if(tmp.lab %in% names(comb.insts)){
    for(tmp.comb in comb.insts[[tmp.lab]]){
      tmp.new.val <- tmp.comb[1]
      tmp.data <- as.character(srt.objs.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]])
      tmp.data[tmp.data %in% tmp.comb] <- tmp.new.val
      srt.objs.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]] <- tmp.data
    }
  }
  # Keep only relevant clusters
  clusts.to.keep <- names(clusts.cols[[tmp.lab]])
  tmp.data <- Cells(srt.objs.list[[tmp.lab]])[
    as.character(srt.objs.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]]) %in% clusts.to.keep
  ]
  srt.objs.list[[tmp.lab]] <- subset(x=srt.objs.list[[tmp.lab]], cells=tmp.data)
  # Set populations as factors.
  srt.objs.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]] <- factor(
    x=as.character(srt.objs.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]]),
    levels=names(clusts.defs[[tmp.lab]])
  )
}

#Define the colom containing sub clusterinformation
t.subset <- 'aggre'
# Extract relevant cell barcodes for the parent aggr NK Srt Obj
cell.barcodes.CD16pos_CX3CR1 <- rownames(srt.objs.list[['CD16pos']]@meta.data[srt.objs.list[['CD16pos']]@meta.data[clust.labs[['CD16pos']]] == names(clusts.defs[['CD16pos']][clusts.defs[['CD16pos']] == 'CX3CR1']),])
cell.barcodes.CD16pos_EGR1 <- rownames(srt.objs.list[['CD16pos']]@meta.data[srt.objs.list[['CD16pos']]@meta.data[clust.labs[['CD16pos']]] == names(clusts.defs[['CD16pos']][clusts.defs[['CD16pos']] == 'EGR1']),])
cell.barcodes.CD16pos_AREG <- rownames(srt.objs.list[['CD16pos']]@meta.data[srt.objs.list[['CD16pos']]@meta.data[clust.labs[['CD16pos']]] == names(clusts.defs[['CD16pos']][clusts.defs[['CD16pos']] == 'AREG']),])
cell.barcodes.CD16pos_Adaptive <- rownames(srt.objs.list[['CD16pos']]@meta.data[srt.objs.list[['CD16pos']]@meta.data[clust.labs[['CD16pos']]] == names(clusts.defs[['CD16pos']][clusts.defs[['CD16pos']] == 'Adaptive']),])
cell.barcodes.CD16pos_TNFRSF9 <- rownames(srt.objs.list[['CD16pos']]@meta.data[srt.objs.list[['CD16pos']]@meta.data[clust.labs[['CD16pos']]] == names(clusts.defs[['CD16pos']][clusts.defs[['CD16pos']] == 'TNFRSF9']),])

cell.barcodes.CD16neg_CD160 <- rownames(srt.objs.list[['CD16neg']]@meta.data[srt.objs.list[['CD16neg']]@meta.data[clust.labs[['CD16neg']]] == names(clusts.defs[['CD16neg']][clusts.defs[['CD16neg']] == 'CD160']),])
cell.barcodes.CD16neg_trNK <- rownames(srt.objs.list[['CD16neg']]@meta.data[srt.objs.list[['CD16neg']]@meta.data[clust.labs[['CD16neg']]] == names(clusts.defs[['CD16neg']][clusts.defs[['CD16neg']] == 'trNK']),])
cell.barcodes.CD16neg_ILC3 <- rownames(srt.objs.list[['CD16neg']]@meta.data[srt.objs.list[['CD16neg']]@meta.data[clust.labs[['CD16neg']]] == names(clusts.defs[['CD16neg']][clusts.defs[['CD16neg']] == 'ILC3']),])
cell.barcodes.CD16neg_CD56_bright <- rownames(srt.objs.list[['CD16neg']]@meta.data[srt.objs.list[['CD16neg']]@meta.data[clust.labs[['CD16neg']]] == names(clusts.defs[['CD16neg']][clusts.defs[['CD16neg']] == 'CD56_bright']),])


srt.objs.list[[t.subset]]@meta.data$barcode <- rownames(srt.objs.list[[t.subset]]@meta.data)

# Create a column defining the desired populations in the metadata
srt.objs.list[[t.subset]]@meta.data$label_transfer_tag <- 'Rest' 
#aggr.Parent_NK.srt.objs.list[[t.subset]]_mod@meta.data$pop.info<-ifelse(aggr.Parent_NK.srt.objs.list[[t.subset]]_mod@meta.data$barcode %in% cell.barcodes.CD16neg_TR.Clust5, "CD16neg_TR.Clust5", "Rest")

srt.objs.list[[t.subset]]@meta.data$label_transfer_tag <-ifelse(srt.objs.list[[t.subset]]@meta.data$barcode %in% cell.barcodes.CD16pos_CX3CR1, "0", srt.objs.list[[t.subset]]@meta.data$label_transfer_tag) 
srt.objs.list[[t.subset]]@meta.data$label_transfer_tag <-ifelse(srt.objs.list[[t.subset]]@meta.data$barcode %in% cell.barcodes.CD16pos_EGR1, "0", srt.objs.list[[t.subset]]@meta.data$label_transfer_tag) 
srt.objs.list[[t.subset]]@meta.data$label_transfer_tag <-ifelse(srt.objs.list[[t.subset]]@meta.data$barcode %in% cell.barcodes.CD16pos_AREG, "0", srt.objs.list[[t.subset]]@meta.data$label_transfer_tag) 
srt.objs.list[[t.subset]]@meta.data$label_transfer_tag <-ifelse(srt.objs.list[[t.subset]]@meta.data$barcode %in% cell.barcodes.CD16pos_Adaptive, "0", srt.objs.list[[t.subset]]@meta.data$label_transfer_tag) 
srt.objs.list[[t.subset]]@meta.data$label_transfer_tag <-ifelse(srt.objs.list[[t.subset]]@meta.data$barcode %in% cell.barcodes.CD16pos_TNFRSF9, "0", srt.objs.list[[t.subset]]@meta.data$label_transfer_tag) 

srt.objs.list[[t.subset]]@meta.data$label_transfer_tag <-ifelse(srt.objs.list[[t.subset]]@meta.data$barcode %in% cell.barcodes.CD16neg_CD160, "1", srt.objs.list[[t.subset]]@meta.data$label_transfer_tag) 
srt.objs.list[[t.subset]]@meta.data$label_transfer_tag <-ifelse(srt.objs.list[[t.subset]]@meta.data$barcode %in% cell.barcodes.CD16neg_trNK, "2", srt.objs.list[[t.subset]]@meta.data$label_transfer_tag) 
srt.objs.list[[t.subset]]@meta.data$label_transfer_tag <-ifelse(srt.objs.list[[t.subset]]@meta.data$barcode %in% cell.barcodes.CD16neg_ILC3, "3", srt.objs.list[[t.subset]]@meta.data$label_transfer_tag) 
srt.objs.list[[t.subset]]@meta.data$label_transfer_tag <-ifelse(srt.objs.list[[t.subset]]@meta.data$barcode %in% cell.barcodes.CD16neg_CD56_bright, "4", srt.objs.list[[t.subset]]@meta.data$label_transfer_tag) 

# Update definitions for "aggre" object.
clust.labs['aggre'] <- 'label_transfer_tag'
clusts.defs[['aggre']] <- c(
  `0`='CD16pos',
  #`1`='EGR1',
  #`2`='AREG',
  #`3`='Adaptive',
  #`4`='TNFRSF9',
  `1`='CD160',
  `2`='trNK',
  `3`='ILC3',
  `4`='CD56_bright'
)
# Sanity checks.
tmp.data <- as.data.table(srt.objs.list[['aggre']]@meta.data)
table(tmp.data[, .(label_transfer_tag, RNA_snn_res.0.7)])
# Cells left unlabeled ("Rest") in the original aggr shouldn't be present the in the final NK cell subset aggrs.
tmp.vals <- c('0', '3')
for(tmp.val in tmp.vals){
  print(tmp.data[
    RNA_snn_res.0.7==tmp.val & label_transfer_tag=='Rest',
    .(
      pos.freq.rel=sum(barcode %in% Cells(srt.objs.list[['CD16pos']])) / .N,
      neg.freq.rel=sum(barcode %in% Cells(srt.objs.list[['CD16neg']])) / .N
    )
  ])
}
# Estimate of cells eliminated during processing for each subset.
#   CD16pos NK, total: 4,102
#       HSP cluster (cluster 5) in round 1 of CD16pos subclustering: 4102
#   CD16neg NK, total: 3,970
#       HSP cluster (cluster 4) in round 1 of CD16neg subclustering: 3087
#       Small undefinable cluster (cluster 10) in round 1 of CD16neg subclustering: 379
#       B cell cluster (cluster 8) in round 2 of CD16neg subclustering: 504
# Retrieve number of cells in the original NK cell aggregation that were lost during the subprocessing steps, and compare against estimated values.


# ---> Obtain cell subset of each dataset that's equal across cell types.
subset.list <- lapply(X=srt.objs.list, FUN=function(seurat.obj){
  if(subset.n>=length(Cells(seurat.obj))){
    return(Cells(seurat.obj))
  }else{
    sample(x=Cells(seurat.obj), size=subset.n)
  }
})

# ---> Donor metadata to be consistently displayed.
# Categorize variables.
# Patient metadata
metadata <- read.csv(metadata.path)
colnames(metadata) <- c('Sample', 'donor.id.tag', 'age', 'Gender', 'Staging')

metadata <- metadata %>% left_join(CD16neg.data, by = 'donor.id.tag')
metadata <- metadata %>% left_join(CD45.data, by = 'donor.id.tag')
facs.data <- facs.data %>% 
  filter(tissue.type == 'CT')
metadata <- metadata %>% left_join(facs.data, by = 'donor.id.tag')
metadata <- metadata %>% left_join(CD8.data, by = 'donor.id.tag')
metadata <- metadata %>% left_join(MHC.data, by = 'donor.id.tag')

#metadata <- metadata %>% left_join(TIL.info, by = 'donor.id.tag')

ordered_donors <- metadata[order(metadata$trNK), ]$donor.id.tag
metadata.trNK <- metadata[, c('donor.id.tag', 'trNK')]
metadata.CD45 <- metadata[, c('donor.id.tag', 'mean_CD45permm')]
metadata.CD8 <- metadata[, c('donor.id.tag', 'CD3_CD8_meancells_per_mm2')]
metadata.MHC <- metadata[, c('donor.id.tag', 'MHC_pos_Panck_pos_pct')]

#metadata <- metadata[, c('donor.id.tag', 'age', 'Gender', 'Staging')]
donor.meta.hmap.trNK <- copy(metadata.trNK)
donor.meta.hmap.CD45 <- copy(metadata.CD45)
donor.meta.hmap.CD8 <- copy(metadata.CD8)
donor.meta.hmap.MHC <-  copy(metadata.MHC)

donor.meta.hmap <- copy(metadata)
setDT(donor.meta.hmap)
# @ Set donor order
#tmp.vals <- factor(x=as.character(donor.meta.hmap$donor.id.tag), levels=rev(custom.donor.ord))
#set(x=donor.meta.hmap, j='donor.id.tag', value=tmp.vals)

# @ Gender
donor.meta.lvls <- list(gender=c('F', 'M'))

# @ Age
# 4 categories. Chosen according to quartiles, then adjusted for ease in communication (according to 5-y bins). <66, 66-70, 71-75, >75
# donor.meta.hmap[, quantile(x=age, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)]
# donor.meta.hmap[, quantile(x=age, probs=c(0.05, 0.1, 0.15, 0.20, 0.25), na.rm=TRUE)]
# donor.meta.hmap[, quantile(x=age, probs=c(0.75, 0.8, 0.9, 0.95, 0.99), na.rm=TRUE)]
# NOTE VFR: The groups should be adjusted according to distribution of the current cohort.
donor.meta.hmap[,
                age:=ifelse(
                  test=age<66,
                  yes='<66',
                  no=ifelse(
                    test=age<71,
                    yes='66-70',
                    no=ifelse(
                      test=age<76,
                      yes='71-75',
                      no='>75'
                    )
                  )
                )
]

donor.meta.lvls[['age']] <-  c('<66', '66-70', '71-75', '>75')

donor.meta.hmap[, Gender := as.character(Gender)]
donor.meta.hmap[, Staging := as.character(Staging)]
#donor.meta.hmap[, til.score := as.character(til.score)]
# @ Tidy format.
donor.meta.hmap <- as.data.table(gather(data=donor.meta.hmap, 
                                        key='donor.var', 
                                        value='var.val', 
                                        -`donor.id.tag`))
donor.meta.hmap[!is.na(var.val), fill.val := paste(donor.var, var.val, sep='.')]
donor.meta.hmap$donor.id.tag <- factor(donor.meta.hmap$donor.id.tag, levels = ordered_donors)

donor.meta.hmap.trNK <- as.data.table(gather(data=donor.meta.hmap.trNK, 
                                             key='donor.var', 
                                             value='var.val', 
                                             -`donor.id.tag`))
#donor.meta.hmap[!is.na(var.val), fill.val := paste(donor.var, var.val, sep='.')]
donor.meta.hmap.trNK$donor.id.tag <- factor(donor.meta.hmap.trNK$donor.id.tag, levels = ordered_donors)

donor.meta.hmap.CD45 <- as.data.table(gather(data=donor.meta.hmap.CD45, 
                                             key='donor.var', 
                                             value='var.val', 
                                             -`donor.id.tag`))
donor.meta.hmap[!is.na(var.val), fill.val := paste(donor.var, var.val, sep='.')]
donor.meta.hmap.CD45$donor.id.tag <- factor(donor.meta.hmap.trNK$donor.id.tag, levels = ordered_donors)

donor.meta.hmap.CD8 <- as.data.table(gather(data=donor.meta.hmap.CD8, 
                                            key='donor.var', 
                                            value='var.val', 
                                            -`donor.id.tag`))
#donor.meta.hmap[!is.na(var.val), fill.val := paste(donor.var, var.val, sep='.')]
donor.meta.hmap.CD8$donor.id.tag <- factor(donor.meta.hmap.CD8$donor.id.tag, levels = ordered_donors)

donor.meta.hmap.MHC <- as.data.table(gather(data=donor.meta.hmap.MHC, 
                                            key='donor.var', 
                                            value='var.val', 
                                            -`donor.id.tag`))
#donor.meta.hmap[!is.na(var.val), fill.val := paste(donor.var, var.val, sep='.')]
donor.meta.hmap.MHC$donor.id.tag <- factor(donor.meta.hmap.MHC$donor.id.tag, levels = ordered_donors)

# @ Define colors.
donor.meta.cols <- c(
  'Gender.F'='#BE82BF',
  'Gender.M'='#8FD5C9',
  'age.<66'='#93F6FF',
  'age.66-70'='#91D2F3',
  'age.71-75'='#99C1FF',
  'age.>75'='#919EF3',
  'Staging.1a2' = '#F7FBFF',
  'Staging.1a3' = '#C6DBEF',
  'Staging.1b'  = '#6BAED6',
  'Staging.2a'  = '#3182BD',
  'Staging.2b'  = '#08519C',
  'Staging.3a'  = '#756BB1',
  'Staging.3b'  = '#54278F',
  'Staging.4a'  = '#3F007D'
  #'til.score.TIL HIGH' = '#B22222',
  #'til.score.TIL LOW' = '#F08080',
  #'til.score.TIL MOD' = 'white'
)

# ----> Output directory.
fig.phen.comps.path <- paste0(reports.path, '/figure_on_phen_comps')
if(!dir.exists(fig.phen.comps.path)) dir.create(fig.phen.comps.path)
# ---> Comments: Main dataset for this section, CD4 T cell set.
main.obj.name <- NULL



############Figure1-----------------------------------------    ############
###--------------------------- Figure1 --------------------------- ###
############Figure1-----------------------------------------    ############
############Main Figures-----------------------------------------    ############
###--------------------------- Main Figures --------------------------- ###
############Main Figures-----------------------------------------    ############
############Donor meta data generation-----------------------------------------    ############
tmp.reports.path <- paste0(fig.phen.comps.path, '/tile')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
##Save the source data
metadata_save <- metadata[,c('donor.id.tag', 'trNK','Sample', 'donor.id.tag', 'age', 'Gender', 'Staging',
   'mean_CD45permm', 'CD3_CD8_meancells_per_mm2', 'MHC_pos_Panck_pos_pct')]
metadata
metadata_save %>% 
  write.csv(paste0(tmp.reports.path, '/source_data_tile.csv'))

# @ Version w/ full set of selected variables.donor.id.tag

tmp.ggplot <- ggplot(data=donor.meta.hmap, aes(x=donor.id.tag, y= factor(donor.var, levels = c("til.score","Staging","Gender", "age")), fill=fill.val)) +
  geom_tile(color='black', linewidth=0.1, linetype=1, alpha=0.8) +
  coord_fixed(ratio = 1.6)+
  scale_fill_manual(values=donor.meta.cols) +
  labs(x='', y='', fill='Value') +
  theme(axis.text.x=element_text(angle=45))
tmp.lab <- paste0('ptmatadata')
publish.plot(
  tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
  blank.comp=blank.complement.1, do.legend=FALSE, do.rotate=TRUE, height=5, width=10
)


# ---> Individual heatmaps for continuous variables.
tmp.vars <- list(
  `trNK_FreqRel`=donor.meta.hmap.trNK,
  `CD45`=donor.meta.hmap.CD45,
  `CD8`=donor.meta.hmap.CD8,
  `MHC`=donor.meta.hmap.MHC
)
for(tmp.var in names(tmp.vars)){
  tmp.ggplot <- ggplot(data=tmp.vars[[tmp.var]], aes(x=donor.id.tag, y=donor.var, fill=var.val)) +
    geom_tile(color='black', linewidth=0.1, linetype=1, alpha=0.8) +
    coord_fixed(ratio = 2)+
    scale_fill_gradientn(colors = cont.col.scale.1) +
    labs(x='', y='', fill='Value') +
    theme(axis.text.x=element_text(angle=45))
  tmp.lab <- paste0(tmp.var)
  publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.1, do.legend=TRUE, do.rotate=TRUE, height=5, width=10
  )
}
# Sanity checks
all(as.character(donor.meta.hmap$donor.id.tag) == as.character(donor.meta.hmap.trNK$donor.id.tag))
all(as.character(donor.meta.hmap$donor.id.tag) == as.character(donor.meta.hmap.CD45$donor.id.tag))
all(as.character(donor.meta.hmap$donor.id.tag) == as.character(donor.meta.hmap.CD8$donor.id.tag))
all(as.character(donor.meta.hmap$donor.id.tag) == as.character(donor.meta.hmap.MHC$donor.id.tag))



color_lists <- list(
  Staging_1 = c('Staging.1a2' = '#F7FBFF', 'Staging.1a3' = '#C6DBEF', 
                'Staging.1b'  = '#6BAED6', 'Staging.2a'  = '#3182BD'),
  Staging = c('Staging.2b'  = '#08519C', 'Staging.3a'  = '#756BB1',
              'Staging.3b'  = '#54278F', 'Staging.4a'  = '#3F007D'),
  Age     = c('age.<66'='#93F6FF', 'age.66-70'='#91D2F3', 'age.71-75'='#99C1FF', 'age.>75'='#919EF3'),
  Gender  = c('Gender.F'='#BE82BF', 'Gender.M'='#8FD5C9')
  #til_score = c(  'til.score.TIL HIGH' = '#B22222','til.score.TIL LOW' = '#F08080','til.score.TIL MOD' = 'white')
)

# Function for legend
make_custom_legend <- function(colors, title) {
  df <- data.frame(
    label = factor(names(colors), levels = rev(names(colors))),
    x = 1
  )
  
  ggplot(df, aes(x = x, y = label, fill = label)) +
    geom_tile(color = "black", size = 1) +
    #geom_text(aes(label = label), size = 4) + # 文字入れが必要な場合
    scale_fill_manual(values = colors) +
    labs(title = title) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    coord_fixed(ratio = 0.4) # 
}

# Generate plot
legend_plots <- lapply(names(color_lists), function(nm) {
  make_custom_legend(color_lists[[nm]], nm)
})

unit_h <- 0.5 #
for(i in seq_along(legend_plots)) {
  title_nm <- names(color_lists)[i]
  row_count <- length(color_lists[[i]])
  
  # 
  dynamic_height <- (row_count * unit_h) + 0.5 
  
  ggsave(
    filename = paste0(tmp.reports.path,"/Legend_", title_nm, ".tiff"),
    plot = legend_plots[[i]],
    width = 2,           # 
    height = dynamic_height, # 
    dpi = 300
  )
}


############UMAP plots depicting clusters.-----------------------------------------    ############
###--------------------------- UMAP plots depicting clusters.--------------------------- ###
############UMAP plots depicting clusters.-----------------------------------------    ############
tmp.reports.path <- paste0(fig.phen.comps.path, '/UMAP')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

for(t.subset in names(obj.extended.names)){
  cat(t.subset, '\n')
  tmp.vals <- c(
    'CD16neg',
    'CD16pos'
  )
  dot.size <- if(t.subset %in% tmp.vals) 0.05 else 0.8
  # Whole dataset.
  tmp.ggplot <- get.umap.gg(obj.name=t.subset, attempted.format='tiff', cell.subset=NULL, dot.size=dot.size)
  tmp.lab <- paste0(obj.extended.names[t.subset], '_PopsOnUMAP_Whole')
  publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='tiff',
    blank.comp=blank.complement.1, do.legend=FALSE,
    height=10, width=10
  )
  # Subset
  tmp.ggplot <- get.umap.gg(obj.name=t.subset, attempted.format='tiff', cell.subset=subset.list[[t.subset]], dot.size=dot.size)
  tmp.lab <- paste0(obj.extended.names[t.subset], '_PopsOnUMAP_Subset')
  publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='tiff',
    blank.comp=blank.complement.1, do.legend=FALSE,
    height=10, width=10
  )
}

## Without density map
for(t.subset in names(obj.extended.names)){
  ## Fetch data
  cat(t.subset, '\n')
  tmp.vals <- c(
    'CD16neg',
    'CD16pos'
  )
  dot.size <- if(t.subset %in% tmp.vals) 0.05 else 0.8
  tmp.tag <- 'tissue.tag'
  
  #Get meta data
  tmp.data.1 <- as.data.frame(srt.objs.list[[t.subset]]@meta.data)
  tmp.data.1 <- tmp.data.1[!is.na(tmp.tag)]
  tmp.data.1$barcode <- row.names(tmp.data.1)
  
  #Get UMAP data
  tmp.data.2 <- as.data.frame(srt.objs.list[[t.subset]]@reductions$umap@cell.embeddings)
  tmp.data.2$barcode <- row.names(tmp.data.2)
  
  #Define the tag (in this case, I used tissue.tag)
  tmp.types <- as.character(unique(tmp.data.1$tissue.tag))
  print(tmp.types)
  
  #merge tabel
  tmp.data <- merge(
    x=tmp.data.1, y=tmp.data.2,
    by='barcode'
  )
  alpha.scale <- c(
    `no`=0.1,
    `yes`=1
  )
  size.scale <- c(
    `no`=0.1,
    `yes`=1.5
  )
  
  xlim_global <- range(tmp.data$UMAP_1, na.rm=TRUE)
  ylim_global <- range(tmp.data$UMAP_2, na.rm=TRUE)
  
  
  for (tmp.type in tmp.types){
    tmp.ggplot <-ggplot(data = tmp.data[tmp.data[tmp.tag] == tmp.type,], aes_string(x='UMAP_1', y = 'UMAP_2', col=clust.labs[[t.subset]] #color=col.val,
    ) )+geom_point(alpha=1, size=dot.size)+scale_color_manual(values=clusts.cols[[t.subset]])+
      labs(x='UMAP 1', y='UMAP 2', col='', fill='Cluster')+
      #geom_density_2d(
      #  data=tmp.data[tmp.data[tmp.tag] == tmp.type,],
      #  linewidth=1.5, color='black'
      #)
      scale_x_continuous(limits=xlim_global, expand=expansion(add=0)) +
      scale_y_continuous(limits=ylim_global, expand=expansion(add=0)) +
      coord_fixed() 
    tmp.lab <- paste0(obj.extended.names[t.subset], '_wo_density', tmp.type)
    publish.plot(
      tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='tiff',
      blank.comp=blank.complement.1, do.legend=FALSE,
      height=10, width=10
    )
  }
}


############Featureplot-----------------------------------------    ############
### ----------------------- Featureplot ----------------------- ###
############Featureplot-----------------------------------------    ############
tmp.reports.path <- paste0(fig.phen.comps.path, '/Featureplot')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
t.subset <- 'aggre'

#Define genesetes to be ploted
genes <- c('FCGR3A')
tn.tags <- c('tumor.lung', 'normal.lung')
for (tn.tag in tn.tags){
  cell.subset <- rownames(srt.objs.list[[t.subset]]@meta.data[srt.objs.list[[t.subset]]@meta.data$tissue.tag == tn.tag,])
  for (tmp.gene in genes){
    tmp.ggplot <- get.umap.feature.gg(obj.name=t.subset,feature = tmp.gene, cell.subset=cell.subset, hmap.center =2, ordered = FALSE)
    tmp.lab <- paste0(obj.extended.names[t.subset],tn.tag ,tmp.gene,'_PopsOnUMAP_Whole')
    tmp.ggplot
    publish.plot(
      tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='tiff',
      blank.comp=blank.complement.1, do.legend=TRUE,
      height=10, width=10
    )
  } 
}


tmp.reports.path <- paste0(fig.phen.comps.path, '/Featureplot')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
t.subset <- 'CD16neg'
#Define genesetes to be ploted
genes <- c('ITGAE', 'ITGA1', 'ZNF683')
tn.tags <- c('tumor.lung', 'normal.lung')

for (tn.tag in tn.tags){
  cell.subset <- rownames(srt.objs.list[[t.subset]]@meta.data[srt.objs.list[[t.subset]]@meta.data$tissue.tag == tn.tag,])
  for (tmp.gene in genes){
    if (tmp.gene == 'ZNF683'){
      tmp.ggplot <- get.umap.feature.gg(obj.name=t.subset,feature = tmp.gene, cell.subset=cell.subset, hmap.center =1, ordered = FALSE)
    }else{
      tmp.ggplot <- get.umap.feature.gg(obj.name=t.subset,feature = tmp.gene, cell.subset=cell.subset, hmap.center =1, ordered = TRUE)
    }
    tmp.lab <- paste0(obj.extended.names[t.subset],tn.tag ,tmp.gene,'Whole')
    publish.plot(
      tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='tiff',
      blank.comp=blank.complement.1, do.legend=TRUE,
      height=10, width=10
    )
  } 
}


############FGSEA-----------------------------------------    ############
### ----------------------- FGSEA ----------------------- ###
###########FGSEA-----------------------------------------    #############
tmp.reports.path <- paste0(fig.phen.comps.path, '/fgsea')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
#Do GSEA trNK signature Nature Communications volume 10, Article number: 3841 (2019) 
####FGSEA aubtract setting
# @ Seeting definition
# Reference group
# Cluster definition pattern
# Signature
# Cluster definition pattern for sbtraction
gsea.settings <- list(
  `CD16neg`=list(
    `trNK`=c('trNK', '^trNK$', 'trNK')
  )
)

for(data.set in names(gsea.settings)){
  cat(paste0(obj.extended.names[data.set], '\n'))
  this.reports.path <- paste0(tmp.reports.path, '/', obj.extended.names[data.set])
  if(!dir.exists(this.reports.path)) dir.create(this.reports.path)
  for(tmp.eval in names(gsea.settings[[data.set]])){
    tmp.ref <- gsea.settings[[data.set]][[tmp.eval]][1]
    tmp.pttn <- gsea.settings[[data.set]][[tmp.eval]][2]
    tmp.mod <- gsea.settings[[data.set]][[tmp.eval]][3]
    srt.objs.list[[data.set]]@meta.data[, 'tmp.tag'] <- ifelse(
      test=str_detect(
        string=clusts.defs[[data.set]][
          as.character(srt.objs.list[[data.set]]@meta.data[, clust.labs[data.set]])
        ],
        pattern=tmp.pttn
      ),
      yes=tmp.ref,
      no='Rest'
    )
    do.fgsea(
      metrics.src=srt.objs.list[[data.set]],
      modules.feats=modules.feats[tmp.mod], 
      metric='signal.to.noise', cell.subset=subset.list[[data.set]],
      tag.of.int='tmp.tag',
      output.path=this.reports.path,
      vals.to.depict=tmp.ref,
      files.preffix=paste0(obj.extended.names[data.set])
    )
    srt.objs.list[[data.set]]@meta.data[, 'tmp.tag'] <- NULL
  }
}


###########Dot plot-----------------------------------------    ############
### ----------------------- Dot plot ----------------------- ###
###########Dot plot-----------------------------------------    ############

tmp.reports.path <- paste0(fig.phen.comps.path, '/dotplot')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

# Define genes of interest.
these.markers <- list(
  `CD16pos`=c(
    # ---> Adaptive
    'KLRC2', 'CD3E',
    # ---> TNFRSF9
    'TNFRSF9', 'CCL4L2',
    # ---> CX3CR1
    'CX3CR1', 'PTGDS',
    # ---> EGR1
    'EGR1',
    # ---> AREG
    'AREG', 'DUSP2'
  ),
  `CD16neg`=c(
    # ---> CD160
    #'CD160', 'CCL3',
    # ---> ILC3
    #'KIT', 'AHR', 'IL7R',
    # ---> Circ
    #'SELL', 'S1PR5', 'KLF2',
    # ---> Tissue resident
    'ITGAE', 'ITGA1', 'ZNF683'
  ),
  `trNK`=c(
    # ---> GZMK
    'GZMK', 'CCL3', 'CCL4',
    # ---> Adaptive tr
    'KLRC2', 'KIR2DL3','KIR2DL1', 'KIR3DL1', 
    # ---> ieILC1
    'S100A13', 'NCR2', 'ENTPD1', 'TMIGD2',
    # ---> CREM
    'LMNA', 'IFNGR1', 'CREM', 'SATB1',
    # ---> XCL1
    'IL7R', 'XCL1', 'TCF7', 'LTB'
  ),
  `aggre`= rev(c(
    'ITGAE', 'ITGA1', 'ZNF683', 'CD160',  'SELL','KIT'
    # ---> GZMK
    #'GZMK', 'CCL3', 'CCL4',
    # ---> Adaptive tr
    #'KLRC2', 'KIR2DL3','KIR2DL1', 'KIR3DL1', 
    # ---> ieILC1
    #'S100A13', 'NCR2', 'ENTPD1', 'TMIGD2',
    # ---> CREM
    #'LMNA', 'IFNGR1', 'CREM', 'SATB1',
    # ---> XCL1
    #'IL7R', 'XCL1', 'TCF7', 'LTB'
  ))
)

# Define clusters' order.
these.pops <- list(
  `CD16pos`=c('4', '5', '0', '2', '3'),
  `CD16neg`=c('4', '6', '3', '0'),
  `trNK`=c('2', '0', '1', '3', '4'),
  `aggre`=c('2', '1', '4', '3', '0')
)

t.subset <- c("aggre")

# Process per lineage.

#if(t.subset == "aggre"){
#  next
#}
# Plot
tmp.markers <- lapply(X = these.markers, FUN = function(tmp.feats){
  tmp.ids <- feature.info[feature.info$name %in% tmp.feats,]
  tmp.ids <- tmp.ids[match(tmp.feats, tmp.ids$name), ]
  tmp.ids <- tmp.ids$id
  tmp.ids <- rev(tmp.ids)
  return(tmp.ids)
})
tmp.gene.names <- lapply(X = these.markers, FUN = function(tmp.feats){
  tmp.names <- feature.info[feature.info$name %in% tmp.feats,]
  tmp.names <- tmp.names[match(tmp.feats, tmp.names$name), ]
  tmp.names <- tmp.names$name
  tmp.names <- rev(tmp.names)
  return(tmp.names)
})

named.labels <-  setNames(tmp.gene.names[[t.subset]], tmp.markers[[t.subset]])

tmp.srt.obj <- subset(srt.objs.list[[t.subset]], subset = label_transfer_tag != 'Rest')


tmp.ggplot <- dot.plot(
  seurat.obj=tmp.srt.obj, features=tmp.markers[[t.subset]], slot='data', ensembl=FALSE,
  # do.norm=FALSE,
  do.norm=TRUE,
  groups.tag=clust.labs[t.subset], groups.order=these.pops[[t.subset]], groups.of.int=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=FALSE, na.rm=TRUE, feature.thold=NULL,
  # this.color.scale=signatures.col.scale,
  this.color.scale='hot.and.cold',
  col.min=NULL, col.max=NULL,
  scale.by='radius', dot.scale=12,
  size.min=0, size.max=NA,
  file.name=NULL,
  named.labels = named.labels
)

tmp.ggplot <- tmp.ggplot + theme(legend.position='bottom')
tmp.ggplot <- tmp.ggplot #+ coord_flip()
tmp.ggplot
# Output.
#tmp.width <- (length(these.pops[[t.subset]]) * 0.5) + 0.3
tmp.height <- (length(these.markers[[t.subset]]) * 0.4) + 0.3
# tmp.width <- 5
#tmp.height <- (length(these.markers[[t.subset]]) * 0.4) + 0.3

tmp.width <- (length(these.pops[[t.subset]]) * 0.4) + 0.3
tmp.lab <- paste0(obj.extended.names[t.subset], '_Markers_DotPlot_Opt-A')
publish.plot(
  tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf', blank.comp=blank.complement.2, do.legend=TRUE,
  width=tmp.width, height=tmp.height, c.width=tmp.width+2, c.height=tmp.height, legend.height=5, legend.width=0.5
)


tmp.reports.path <- paste0(fig.phen.comps.path, '/dotplot')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

# Define genes of interest.
these.markers <- list(
  `CD16pos`=c(
    # ---> Adaptive
    'KLRC2', 'CD3E',
    # ---> TNFRSF9
    'TNFRSF9', 'CCL4L2',
    # ---> CX3CR1
    'CX3CR1', 'PTGDS',
    # ---> EGR1
    'EGR1',
    # ---> AREG
    'AREG', 'DUSP2'
  ),
  `CD16neg`=c(
    # ---> CD160
    #'CD160', 'CCL3',
    # ---> ILC3
    #'KIT', 'AHR', 'IL7R',
    # ---> Circ
    #'SELL', 'S1PR5', 'KLF2',
    # ---> Tissue resident
    'ITGAE', 'ITGA1', 'ZNF683'
  ),
  `trNK`=c(
    # ---> GZMK
    'GZMK', 'CCL3', 'CCL4',
    # ---> Adaptive tr
    'KLRC2', 'KIR2DL3','KIR2DL1', 'KIR3DL1', 
    # ---> ieILC1
    'S100A13', 'NCR2', 'ENTPD1', 'TMIGD2',
    # ---> CREM
    'LMNA', 'IFNGR1', 'CREM', 'SATB1',
    # ---> XCL1
    'IL7R', 'XCL1', 'TCF7', 'LTB'
  ),
  `aggre`=c(
    # ---> CX3CR1
    'CX3CR1',
    # ---> EGR1
    'EGR1',
    # ---> AREG
    'AREG',
    # ---> Adaptive
    'KLRC2', 'CD3E',
    # ---> TNFRSF9
    'TNFRSF9', 'CCL4L2',
    # ---> CD160
    'CD160', 
    # ---> ILC3
    'KIT', 'AHR', 'IL7R',
    # ---> Circ
    'SELL', 'KLF2'
    # ---> Tissue resident
    #'ITGAE', 'ITGA1', 'ZNF683'
  )
)


############Figure 2-----------------------------------------    ############
### -----------------------Figure 2 ----------------------- ###
###########Figure 2-----------------------------------------    #############

############Heat map of trNK cells-----------------------------------------    ############
###---------------------------Heat map of trNK cells--------------------------- ###
############Heat map of trNK cells-----------------------------------------    ############

# Define output directory for tile plots
tmp.reports.path <- paste0(fig.phen.comps.path, '/tile')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

t.subset <- 'trNK'
df <- as.data.table(srt.objs.list[['trNK']]@meta.data, keep.rownames = "barcode")
setDT(df)

clust_map <- c("0"="Adaptive_trNK", "1"="ENTPD1", "2"="Cytotoxic", "3"="CREM", "4"="TCF7")

summary_dt <- df[, .N, by = .(full.donor.id.tag, tissue.tag, RNA_snn_res.0.4)]
summary_dt[, total := sum(N), by = .(full.donor.id.tag, tissue.tag)]
summary_dt[, percentage := (N / total) * 100]

summary_dt[, cluster_name := clust_map[as.character(RNA_snn_res.0.4)]]
summary_dt[is.na(cluster_name), cluster_name := as.character(RNA_snn_res.0.4)]

tmp.data <- unique(summary_dt[, .(full.donor.id.tag, tissue.tag, total)])
tmp.data %>% 
  write.csv(paste0(fig.phen.comps.path, '/trNK_table_source.csv'))

# # NOTE VFR: Set a total frequency threshold right at this point, to keep donors with sufficient number of trNK cells on which heterogeneity can be assessed. 
tmp.thold <- 10
summary_dt <- summary_dt[total>=tmp.thold]

final_table <- dcast(
  summary_dt,
  full.donor.id.tag + tissue.tag ~ cluster_name,
  value.var = "percentage",
  fill = 0
)
final_table %>% 
  write.csv(paste0(fig.phen.comps.path, '/trNK_table.csv'))

donor.meta.hmap.trNK <- final_table %>%
  filter(tissue.tag == "tumor.lung") %>%
  filter(!is.na(full.donor.id.tag)) %>%
  pivot_longer(
    cols = c("Adaptive_trNK", "CREM", "Cytotoxic", "ENTPD1", "TCF7"),
    names_to = "cluster_name",
    values_to = "percentage"
  )

donor.meta.hmap.trNK$full.donor.id.tag <- as.factor(donor.meta.hmap.trNK$full.donor.id.tag)

# 1) make a cluster x donor  table
#    ※ pivot_longer
mat_wide <- donor.meta.hmap.trNK %>%
  tidyr::pivot_wider(
    names_from  = full.donor.id.tag,
    values_from = percentage,
    values_fill = 0
  ) %>%
  as.data.frame()

# Rename
rownames(mat_wide) <- mat_wide$cluster_name
mat_wide$cluster_name <- NULL

mat <- as.matrix(mat_wide)

# Parsing） ====
mat <- apply(mat, 2, function(x) as.numeric(as.character(x)))
mat <- as.matrix(mat)
rownames(mat) <- rownames(mat_wide)
colnames(mat) <- colnames(mat_wide)

# NA fill
mat[!is.finite(mat)] <- 0


# # NOTE KM:Corrected the order
# Clustering
#donor.hc <- hclust(dist(t(mat)), method = "ward.D2")
#donor.order <- colnames(mat)[donor.hc$order]

#cluster.hc <- hclust(dist(mat), method = "ward.D2")
#cluster.order <- rownames(mat)[cluster.hc$order]

# 4) ggplot donor name factorization
#donor.meta.hmap.trNK$full.donor.id.tag <- factor(
#  donor.meta.hmap.trNK$full.donor.id.tag,
#  levels = donor.order
#)

#donor.meta.hmap.trNK$cluster_name <- factor(
#  donor.meta.hmap.trNK$cluster_name,
#  levels = cluster.order
#)

# Clustering donors only
donor.hc <- hclust(dist(t(mat)), method = "ward.D2")
donor.order <- colnames(mat)[donor.hc$order]

# Apply donor order
donor.meta.hmap.trNK$full.donor.id.tag <- factor(
  donor.meta.hmap.trNK$full.donor.id.tag,
  levels = donor.order
)

# Set fixed cluster order
donor.meta.hmap.trNK$cluster_name <- factor(
  donor.meta.hmap.trNK$cluster_name,
  levels = c(
    "Adaptive_trNK",
    "ENTPD1",
    "Cytotoxic",
    "CREM",
    "TCF7"
  )
)

tmp.ggplot <- ggplot(
  data = donor.meta.hmap.trNK,
  aes(x = full.donor.id.tag, y = cluster_name, fill = percentage)
) +
  geom_tile(color = 'gray90', linewidth = 0.1, linetype = 1, alpha = 0.8) +
  coord_fixed(ratio = 1.8) +
  scale_fill_gradientn(colors = cont.col.scale.1, name = "Percentage (%)") +
  labs(x = 'Donor ID', y = 'trNK Clusters', title = 'Cluster Composition (Tumor Only)') +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
tmp.lab <- paste0('trNK_pct_heatmap_tumor')
publish.plot(
  tmp.ggplot = tmp.ggplot,
  output.path = tmp.reports.path,
  file.name = tmp.lab,
  type = 'pdf',
  blank.comp = blank.complement.1,
  do.legend = TRUE,
  do.rotate = TRUE,
  height = 5,
  width = 10
)

#### Calculate the percentage of TCF7 positive cells
srt.obj <- srt.objs.list[['trNK']]
srt.obj <- subset(srt.obj, subset = tissue.tag == 'tumor.lung')
donor.col <- "full.donor.id.tag"

# 
tmp.df <- FetchData(
  srt.obj,
  vars = c("ENSG00000081059", donor.col),
  slot = "data"
)

#
res.df <- tmp.df %>%
  mutate(TCF7_pos = ENSG00000081059 > 0) %>%
  group_by(.data[[donor.col]]) %>%
  summarise(
    n_cells = n(),
    n_TCF7_pos = sum(TCF7_pos, na.rm = TRUE),
    frac_TCF7_pos = mean(TCF7_pos, na.rm = TRUE),
    perc_TCF7_pos = 100 * mean(TCF7_pos, na.rm = TRUE)
  ) %>%
  arrange(desc(frac_TCF7_pos))

res.df$full.donor.id.tag <- factor(res.df$full.donor.id.tag, levels = donor.order)

tmp.df <- donor.meta.hmap.trNK %>% 
  left_join(res.df, by = 'full.donor.id.tag')

tmp.df <- tmp.df %>% 
  filter(cluster_name == 'Adaptive_trNK')

mid.point <- 15
hmap.col.scale.1 <- colorRampPalette(c('#8FD4C8', '#9FDAD0', '#ABE0D6', '#BFE8E0', '#CEEBE5', '#DEF2EE', '#EAF8F5', '#FFFFFF'))(mid.point)
#hmap.col.scale.2 <- colorRampPalette(c('#FFFFFF', '#F6EDF4', '#ECDEEC', '#E2C7E2', '#DBBADC', '#CFA5CF', '#C58FC6', '#BC7FBE'))(100-(mid.point+1))

hmap.col.scale.2 <- colorRampPalette(c(
  '#FFFFFF',
  '#EDE7F6',
  '#D1C4E9',
  '#B39DDB',
  '#9575CD',
  '#7E57C2',
  '#673AB7',
  '#4527A0'   #
))(100 - (mid.point + 1))

hmap.col.scale <- c(hmap.col.scale.1, hmap.col.scale.2)

tmp.df
tmp.ggplot <- ggplot(
  data = tmp.df,
  aes(x = full.donor.id.tag, y = cluster_name, fill = perc_TCF7_pos)
) +
  geom_tile(color = 'gray90', linewidth = 0.1, linetype = 1, alpha = 0.8) +
  coord_fixed(ratio = 1.8) +
  scale_fill_gradientn(colors = hmap.col.scale, name = "Percentage (%)") +
  labs(x = 'Donor ID', y = 'trNK Clusters', title = 'Cluster Composition (Tumor Only)') +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
tmp.ggplot
tmp.lab <- paste0('TCF_pos_cells')
publish.plot(
  tmp.ggplot = tmp.ggplot,
  output.path = tmp.reports.path,
  file.name = tmp.lab,
  type = 'pdf',
  blank.comp = blank.complement.1,
  do.legend = TRUE,
  do.rotate = TRUE,
  height = 5,
  width = 10
)






###########Dot plot-----------------------------------------    ############
### ----------------------- Dot plot ----------------------- ###
###########Dot plot-----------------------------------------    ############

tmp.reports.path <- paste0(fig.phen.comps.path, '/dotplot')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

# Define genes of interest.
these.markers <- list(
  `trNK`=c(
    # ---> GZMK
    #'GZMK', 'CCL3', 'CCL4',
    # ---> Adaptive tr
    #'KLRC2', 'KIR2DL3','KIR2DL1', 'KIR3DL1', 
    # ---> ieILC1
    #'S100A13', 'NCR2', 'ENTPD1', 'TMIGD2',
    # ---> CREM
    #'LMNA', 'IFNGR1', 'CREM', 'SATB1',
    # ---> XCL1
    "BTG2","GPR183","CD27","XCL2","XCL1","IL7R","LTB","TCF7"
  )
)
sessionInfo()
# Define clusters' order.
these.pops <- list(
  `CD16pos`=c('4', '5', '0', '2', '3'),
  `CD16neg`=c('4', '6', '3', '0'),
  `trNK`=c('0', '1', '2', '3', '4')
)

t.subset <- c("trNK")

# Process per lineage.

#if(t.subset == "aggre"){
#  next
#}
# Plot
tmp.markers <- lapply(X = these.markers, FUN = function(tmp.feats){
  tmp.ids <- feature.info[feature.info$name %in% tmp.feats,]
  tmp.ids <- tmp.ids[match(tmp.feats, tmp.ids$name), ]
  tmp.ids <- tmp.ids$id
  tmp.ids <- rev(tmp.ids)
  return(tmp.ids)
})
tmp.gene.names <- lapply(X = these.markers, FUN = function(tmp.feats){
  tmp.names <- feature.info[feature.info$name %in% tmp.feats,]
  tmp.names <- tmp.names[match(tmp.feats, tmp.names$name), ]
  tmp.names <- tmp.names$name
  tmp.names <- rev(tmp.names)
  return(tmp.names)
})

named.labels <-  setNames(tmp.gene.names[[t.subset]], tmp.markers[[t.subset]])

tmp.ggplot <- dot.plot(
  seurat.obj=srt.objs.list[[t.subset]], features=tmp.markers[[t.subset]], slot='data', ensembl=FALSE,
  # do.norm=FALSE,
  do.norm=TRUE,
  groups.tag=clust.labs[t.subset], groups.order=these.pops[[t.subset]], groups.of.int=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=FALSE, na.rm=TRUE, feature.thold=NULL,
  # this.color.scale=signatures.col.scale,
  this.color.scale='hot.and.cold',
  col.min=NULL, col.max=NULL,
  scale.by='radius', dot.scale=12,
  size.min=0, size.max=NA,
  file.name=NULL,
  named.labels = named.labels
)

tmp.ggplot <- tmp.ggplot + theme(legend.position='bottom')
tmp.ggplot <- tmp.ggplot
# Output.
#tmp.width <- (length(these.pops[[t.subset]]) * 0.5) + 0.3
tmp.height <- (length(these.markers[[t.subset]]) * 0.4) + 0.3
# tmp.width <- 5
#tmp.height <- (length(these.markers[[t.subset]]) * 0.4) + 0.3
tmp.width <- (length(these.pops[[t.subset]]) * 0.4) + 0.3
tmp.lab <- paste0(obj.extended.names[t.subset], '_Markers_DotPlot_Opt-A')
publish.plot(
  tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf', blank.comp=blank.complement.2, do.legend=TRUE,
  width=tmp.width, height=tmp.height, c.width=tmp.width+2, c.height=tmp.height, legend.height=5, legend.width=0.5
)


############Featureplot-----------------------------------------    ############
### ----------------------- Featureplot ----------------------- ###
############Featureplot-----------------------------------------    ############

tmp.reports.path <- paste0(fig.phen.comps.path, '/Featureplot')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
t.subset <- 'trNK'
#Define genesetes to be ploted
genes <- c('TCF7', 'LTB', 'IL7R')

tmp.reports.path <- paste0(fig.phen.comps.path, '/Featureplot')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)


#Define genesetes to be ploted
genes <- c('TCF7', 'LTB', 'IL7R')
tn.tags <- c('tumor.lung', 'normal.lung')
for (tn.tag in tn.tags){
  cell.subset <- rownames(srt.objs.list[[t.subset]]@meta.data[srt.objs.list[[t.subset]]@meta.data$tissue.tag == tn.tag,])
  for (tmp.gene in genes){
    tmp.ggplot <- get.umap.feature.gg(obj.name=t.subset,feature = tmp.gene, cell.subset=cell.subset, hmap.center =2, dot.size = 2, 
                                      hmap.col.scale = c("lightgrey",
                                                         "#A481E9",
                                                         "blue"))
    tmp.lab <- paste0(obj.extended.names[t.subset],tn.tag ,tmp.gene,'_')
    tmp.ggplot
    publish.plot(
      tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='tiff',
      blank.comp=blank.complement.1, do.legend=TRUE,
      height=10, width=10
    )
  } 
}
for (tn.tag in tn.tags){
  cell.subset <- rownames(srt.objs.list[[t.subset]]@meta.data[srt.objs.list[[t.subset]]@meta.data$tissue.tag == tn.tag,])
  for (tmp.gene in genes){
    tmp.ggplot <- get.umap.feature.gg(obj.name=t.subset,feature = tmp.gene, cell.subset=cell.subset, hmap.center =2, dot.size = 2, 
    )
    tmp.lab <- paste0(obj.extended.names[t.subset],tn.tag ,tmp.gene,'_', '_1')
    tmp.ggplot
    publish.plot(
      tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='tiff',
      blank.comp=blank.complement.1, do.legend=TRUE,
      height=10, width=10
    )
  } 
}


############FGSEA-----------------------------------------    ############
### ----------------------- FGSEA ----------------------- ###
###########FGSEA-----------------------------------------    #############
tmp.reports.path <- paste0(fig.phen.comps.path, '/fgsea')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

#Do GSEA skin.trNK signature Immunity Volume 57, Issue 1, 9 January 2024, Pages 124-140.e7
####FGSEA aubtract setting
# @ Seeting definition
# Reference group
# Cluster definition pattern
# Signature
gsea.settings <- list(
  `trNK`=list(
    `skin_trNK`=c('TCF7_trNK', '^TCF7$', 'skin.trNK')
  )
)
#In tumor trNK
tn.tag<- c('tumor.lung')
gsea.settings <- list(
  `trNK`=list(
    `skin_trNK`=c('TCF7_trNK', '^TCF7$', 'skin.trNK')
  )
)

for(data.set in names(gsea.settings)){
  cat(paste0(obj.extended.names[data.set], '\n'))
  this.reports.path <- paste0(tmp.reports.path, '/', obj.extended.names[data.set])
  if(!dir.exists(this.reports.path)) dir.create(this.reports.path)
  for(tmp.eval in names(gsea.settings[[data.set]])){
    tmp.ref <- gsea.settings[[data.set]][[tmp.eval]][1]
    tmp.pttn <- gsea.settings[[data.set]][[tmp.eval]][2]
    tmp.mod <- gsea.settings[[data.set]][[tmp.eval]][3]
    srt.objs.list[[data.set]]@meta.data[, 'tmp.tag'] <- ifelse(
      test=str_detect(
        string=clusts.defs[[data.set]][
          as.character(srt.objs.list[[data.set]]@meta.data[, clust.labs[data.set]])
        ],
        pattern=tmp.pttn
      ),
      yes=tmp.ref,
      no='Rest'
    )
    do.fgsea(
      metrics.src=subset(srt.objs.list[[data.set]], subset = tissue.tag == tn.tag),
      modules.feats=modules.feats[tmp.mod], 
      metric='signal.to.noise', cell.subset=subset.list[[data.set]],
      tag.of.int='tmp.tag',
      output.path=this.reports.path,
      vals.to.depict=tmp.ref,
      files.preffix=paste0(obj.extended.names[data.set], tn.tag)
    )
    srt.objs.list[[data.set]]@meta.data[, 'tmp.tag'] <- NULL
  }
}

############Figure 3-----------------------------------------    ############
### -----------------------Figure 3 ----------------------- ###
###########Figure 3-----------------------------------------    #############
############Featureplot-----------------------------------------    ############
### ----------------------- Featureplot ----------------------- ###
############Featureplot-----------------------------------------    ############
tmp.reports.path <- paste0(fig.phen.comps.path, '/Featureplot')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
t.subset <- 'trNK'

#Define genesetes to be ploted
genes <- c('KLRC2', 'CD2')
tn.tags <- c('tumor.lung', 'normal.lung')
for (tn.tag in tn.tags){
  cell.subset <- rownames(srt.objs.list[[t.subset]]@meta.data[srt.objs.list[[t.subset]]@meta.data$tissue.tag == tn.tag,])
  for (tmp.gene in genes){
    tmp.ggplot <- get.umap.feature.gg(obj.name=t.subset,feature = tmp.gene, cell.subset=cell.subset, hmap.center =2)
    tmp.lab <- paste0(obj.extended.names[t.subset],tn.tag ,tmp.gene,'_PopsOnUMAP_Whole')
    tmp.ggplot
    publish.plot(
      tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='tiff',
      blank.comp=blank.complement.1, do.legend=TRUE,
      height=10, width=10
    )
  } 
}

############Module scoring calculation-----------------------------------------    ############
### ----------------------- Module scoring calculation ----------------------- ###
###########Module scoring calculation-----------------------------------------    ############
tmp.reports.path <- paste0(fig.phen.comps.path, '/module_scoring')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

mods.of.int <- list(
  `NK`=c(
    `trNK_lung`='trNK',
    `IL21`='IL21',
    `Cytotoxic` = 'Cytotoxic',
    `IL21_bulk` = 'IL21.bulk',
    `skin.trNK` = 'skin.trNK',
    `NK3` = 'NK3',
    `IL21.bulk.109` = 'IL21.bulk.109'
    #`Patil-CTL`='cell.cytotoxicity.patil', `Arlehamn-TH1`='th1.signature1.arlehamn', `Arlehamn-TH17`='th17.signature2.arlehamndan', `Locci-TFH`='tfh.signature.locci', `Best-Cylcing`='cell.cycling.best',
    #`MSIGDB-IFNR-1-2`='type1n2.interferon.broadreact',
    #`DICE-TREGN`='R24.TREGN.RESTING', `DICE-TREGM`='R24.TREGMEM.RESTING', `DICE-TFH`='R24.TFH.RESTING', `DICE-TH2`='R24.TH2.RESTING', `DICE-TH1`='R24.TH1.RESTING', `DICE-TH17`='R24.TH17.RESTING', `DICE-THSTAR`='R24.THSTAR.RESTING', `DICE-CD4N`='R24.CD4N.RESTING'
    # `Clarke-TRM`='trm.signature.clarke', `Pauken-TCM`='TCM.Pauken.Tumor', `Patil-CTL`='cell.cytotoxicity.patil', `DICE-TREGN`='R24.TREGN.RESTING', `DICE-TREGM`='R24.TREGMEM.RESTING'
  )
)

for(data.set in names(srt.objs.list)){
  if(!all(mods.of.int %in% colnames(srt.objs.list[[data.set]]@meta.data))){
    # mod.in.path <- paste0(module.feats.dir, '/base_signatures')
    mod.in.path <- paste0(module.feats.dir) # Replace w/ this if a specific subset of signatures needs to be added and they were included as symbolic links in this folder.
    mod.out.path <- paste0(module.feats.dir, '/tmp_dir')
    if(!dir.exists(mod.out.path)) dir.create(mod.out.path)
    srt.objs.list[[data.set]] <- get.module.scores(
      seurat.obj=srt.objs.list[[data.set]], module.feats.dir=mod.in.path, reports.path=mod.out.path,
      signal.thold=0, is.ensembl=FALSE, feature.info.path = feature.info.path
    )
    # Remove tmp dir.
    tmp.cmmd <- paste0('rm -r ', mod.out.path)
    system(command=tmp.cmmd, intern=FALSE)
  }
}

# Custom score bounds.
gen.scale.tholds <- list(
  #  `CD16pos`=list(
  #    `trNK_lung`=1.5
  #  )
)


#Tumor nontumor separately Module UMAP for all cells with unified scale color -- ver.2
module.in.all.cells <- as.data.table(srt.objs.list[['aggre']]@meta.data)
module.in.all.cells$barcode <- rownames(srt.objs.list[['aggre']]@meta.data)
hmap.center <- 0
tmp.vals <- 'aggre'

tmp.reports.path <- paste0(fig.phen.comps.path, '/module_scoring')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
module.in.all.cells <- as.data.table(srt.objs.list[['aggre']]@meta.data)
module.in.all.cells$barcode <- rownames(srt.objs.list[['aggre']]@meta.data)
hmap.center <- 0
tmp.vals <- c('CD16neg', 'CD16pos', 'aggre')
tissue.conditions <- c('tumor.lung', 'normal.lung')
for (tn.tag in tissue.conditions){
  for(data.set in names(obj.extended.names)){
    
    cat(paste0(obj.extended.names[data.set], '\n'))
    this.reports.path <- paste0(tmp.reports.path, '/', obj.extended.names[data.set])
    if(!dir.exists(this.reports.path)) dir.create(this.reports.path)
    
    dot.size <- if(data.set %in% tmp.vals) 0.5 else 1.0
    file.pffx <- obj.extended.names[data.set]
    
    seurat.obj = srt.objs.list[[data.set]]
    # Get UMAP
    tmp.meta <- cbind(as.data.table(seurat.obj@meta.data), seurat.obj@reductions$umap@cell.embeddings)
    tmp.meta$barcode <- Cells(seurat.obj)
    tmp.meta <- tmp.meta %>% 
      filter(tissue.tag == tn.tag)
    these.mods <- mods.of.int[[gen.cell.types[data.set]]]
    these.mods[!str_detect(string=these.mods, pattern='score')] <- paste0(these.mods[!str_detect(string=these.mods, pattern='score')], '.score')
    
    for(this.mod in names(these.mods)){
      target_col <- unname(these.mods[this.mod])
      
      # 1.Set a scale
      global_range <- range(module.in.all.cells[[target_col]], na.rm = TRUE)
      
      # 2.Merge the data 
      tmp.module.in.all.cells <- module.in.all.cells[, .SD, .SDcols = c('barcode', target_col)]
      tmp.meta_sub <- tmp.meta[, .SD, .SDcols = c('barcode', 'UMAP_1', 'UMAP_2')]
      plot_df <- tmp.module.in.all.cells[tmp.meta_sub, on = .(barcode)]
      setnames(plot_df, target_col, "tmp.tag") #
      
      # 3.Color
      col.breaks <- seq(from=min(range(plot_df$tmp.tag, na.rm=TRUE)), to=max(range(plot_df$tmp.tag, na.rm=TRUE)), length.out=100)
      if (this.mod == 'trNK_lung'){
        mid.point <- 0.3
      }else{
        if (this.mod == 'NK3'){
          mid.point <- 0.45
        }else {
          if (this.mod == 'skin.trNK'){
            mid.point <- 0.75
          }
          else{
            if (this.mod == 'Cytotoxic'){
              mid.point <- 0.8
            }else {
              if (this.mod == 'IL21.bulk.109'){
                mid.point <- 0.1
              }else {mid.point <- hmap.center}
             }
            }
          }
      }
      
      mid.point <- which.min(abs(col.breaks - mid.point))
      hmap.col.scale.1 <- colorRampPalette(c('#8FD4C8', '#9FDAD0', '#ABE0D6', '#BFE8E0', '#CEEBE5', '#DEF2EE', '#EAF8F5', '#FFFFFF'))(mid.point)
      
      hmap.col.scale.2 <- colorRampPalette(c('#FFFFFF','#EDE7F6','#D1C4E9','#B39DDB','#9575CD','#7E57C2','#673AB7','#4527A0'))(100 - (mid.point + 1))
      
      #hmap.col.scale.2 <- colorRampPalette(c('#FFFFFF', '#F6EDF4', '#ECDEEC', '#E2C7E2', '#DBBADC', '#CFA5CF', '#C58FC6', '#BC7FBE'))(100-(mid.point+1))
      hmap.col.scale <- c(hmap.col.scale.1, hmap.col.scale.2)
      # 4.
      # Plot.
      tmp.ggplot <- ggplot(data=plot_df, aes(x=UMAP_1, y=UMAP_2, col=tmp.tag)) +
        geom_point(size=dot.size, alpha=1, shape=19) +
        scale_color_gradientn(colors=hmap.col.scale) +
        scale_y_continuous(expand=expansion(add=0)) + scale_x_continuous(expand=expansion(add=0)) +
        labs(x='UMAP 1', y='UMAP 2', col='Score')
      
      # Output
      publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=this.reports.path,
        file.name=paste0(tn.tag, file.pffx, 'Col_Vi_All_cell_scale_SigUMAPHmap_Mod-', this.mod),
        type='tiff', blank.comp=blank.complement.1, do.legend=TRUE,
        height=10, width=10
      )
    }
  }
}

############FGSEA-----------------------------------------    ############
### ----------------------- FGSEA ----------------------- ###
###########FGSEA-----------------------------------------    #############
tmp.reports.path <- paste0(fig.phen.comps.path, '/fgsea')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)


#Do GSEA NK3 signature Nature Immunology volume 25, pages1474–1488 (2024)
####FGSEA aubtract setting
# @ Seeting definition
# Reference group
# Cluster definition pattern
# Signature
# Extract relevant cell barcodes for the parent aggr NK Srt Obj
cell.barcodes.CD16pos_adaptive <- rownames(srt.objs.list[['CD16pos']]@meta.data[srt.objs.list[['CD16pos']]@meta.data[clust.labs[['CD16pos']]] == names(clusts.defs[['CD16pos']][clusts.defs[['CD16pos']] == 'Adaptive']),])
cell.barcodes.CD16pos_non_adaptive <- setdiff(rownames(srt.objs.list[['CD16pos']]@meta.data), cell.barcodes.CD16pos_adaptive)

#cell.barcodes.CD16neg_CD56_bright <- rownames(srt.objs.list[['CD16neg']]@meta.data[srt.objs.list[['CD16neg']]@meta.data[clust.labs[['CD16neg']]] == names(clusts.defs[['CD16neg']][clusts.defs[['CD16neg']] == 'CD56_bright']),])
#cell.barcodes.CD16neg_CD160 <- rownames(srt.objs.list[['CD16neg']]@meta.data[srt.objs.list[['CD16neg']]@meta.data[clust.labs[['CD16neg']]] == names(clusts.defs[['CD16neg']][clusts.defs[['CD16neg']] == 'CD160']),])
#cell.barcodes.CD16neg_ILC3 <- rownames(srt.objs.list[['CD16neg']]@meta.data[srt.objs.list[['CD16neg']]@meta.data[clust.labs[['CD16neg']]] == names(clusts.defs[['CD16neg']][clusts.defs[['CD16neg']] == 'ILC3']),])
#cell.barcodes.CD16neg_trNK <- rownames(srt.objs.list[['CD16neg']]@meta.data[srt.objs.list[['CD16neg']]@meta.data[clust.labs[['CD16neg']]] == names(clusts.defs[['CD16neg']][clusts.defs[['CD16neg']] == 'trNK']),])

cell.barcodes.Adaptive_trNK <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'Adaptive_trNK']),])
cell.barcodes.ENTPD1 <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'ENTPD1']),])
cell.barcodes.Cytotoxic <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'Cytotoxic']),])
cell.barcodes.CREM <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'CREM']),])
cell.barcodes.TCF7 <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'TCF7']),])
srt.obj <- srt.objs.list[['aggre']]
srt.obj@meta.data$barcode <- rownames(srt.obj@meta.data)

# Create a column defining the desired populations in the metadata
srt.obj@meta.data$pop.info <- 'Rest'
#aggr.Parent_NK.srt.obj_mod@meta.data$pop.info<-ifelse(aggr.Parent_NK.srt.obj_mod@meta.data$barcode %in% cell.barcodes.CD16neg_TR.Clust5, "CD16neg_TR.Clust5", "Rest")
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.CD16pos_adaptive, "Adaptive", srt.obj@meta.data$pop.info)
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.CD16pos_non_adaptive, "CD16pos_non_adaptive", srt.obj@meta.data$pop.info)

srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.Adaptive_trNK, "Adaptive_trNK", srt.obj@meta.data$pop.info)
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.ENTPD1, "ENTPD1", srt.obj@meta.data$pop.info)
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.Cytotoxic, "Cytotoxic", srt.obj@meta.data$pop.info)
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.CREM, "CREM", srt.obj@meta.data$pop.info)
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.TCF7, "TCF7", srt.obj@meta.data$pop.info)

#########Tumor non-tumor separately
gsea.settings <- list(
  `aggre`=list(
    `NK3_adaptive`=c('NK3_AdaptivevsCD16pos', '^Adaptive$', 'NK3', 'Adaptive','CD16pos_non_adaptive'),
    `NK3_trNK_adaptive`=c('NK3_AdaptivetrNKvsCD16pos', '^Adaptive_trNK$', 'NK3', 'Adaptive_trNK','CD16pos_non_adaptive'),
    `NK3_ENPTD1`=c('NK3_ENTPD1vsCD16pos', '^ENTPD1$', 'NK3', 'ENTPD1','CD16pos_non_adaptive'),
    `NK3_Cytotoxic`=c('NK3_CytvsCD16pos', '^Cytotoxic$', 'NK3', 'Cytotoxic','CD16pos_non_adaptive'),
    `NK3_CREM`=c('NK3_CREMvsCD16pos', '^CREM$', 'NK3', 'CREM','CD16pos_non_adaptive'),
    `NK3_TCF7`=c('NK3_TCF7vsCD16pos', '^TCF7$', 'NK3', 'TCF7','CD16pos_non_adaptive')
  )
)
data.set <- 'aggre'
srt.obj <- subset(srt.obj, subset = tissue.tag == 'tumor.lung')
this.reports.path <- paste0(tmp.reports.path, '/', obj.extended.names[data.set])
if(!dir.exists(this.reports.path)) dir.create(this.reports.path)

for(tmp.eval in names(gsea.settings[[data.set]])){
  tmp.ref <- gsea.settings[[data.set]][[tmp.eval]][1]
  tmp.pttn <- gsea.settings[[data.set]][[tmp.eval]][2]
  tmp.mod <- gsea.settings[[data.set]][[tmp.eval]][3]
  tmp.sub <- gsea.settings[[data.set]][[tmp.eval]][4]
  tmp.sub.2 <- gsea.settings[[data.set]][[tmp.eval]][5]
  srt.obj@meta.data[, 'tmp.tag'] <- ifelse(
    test=str_detect(
      string=as.character(srt.obj@meta.data[, 'pop.info']),
      pattern=tmp.pttn
    ),
    yes=tmp.ref,
    no='Rest'
  )
  #Define cells for subsetting
  sub.cells <- rownames(srt.obj@meta.data[(srt.obj@meta.data$pop.info %in% c(tmp.sub, tmp.sub.2)),])
  
  do.fgsea(
    metrics.src=(srt.obj),
    modules.feats=modules.feats[tmp.mod], 
    metric='signal.to.noise', 
    cell.subset  = sub.cells,
    tag.of.int='tmp.tag',
    output.path=this.reports.path,
    vals.to.depict=tmp.ref,
    files.preffix=paste0(obj.extended.names[data.set])
  )
  srt.obj@meta.data[, 'tmp.tag'] <- NULL
}


###########Dot plot-----------------------------------------    ############
### ----------------------- Dot plot ----------------------- ###
###########Dot plot-----------------------------------------    ############

tmp.reports.path <- paste0(fig.phen.comps.path, '/dotplot')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)


t.subset <- c("aggre")
tmp.srt.obj <- srt.objs.list[[t.subset]]


# Define genes of interest.
these.markers <- list(
  `aggre`=rev(c(
    # ---> GZMK
    #'GZMK', 'CCL3', 'CCL4',
    # ---> Adaptive tr
    #'KLRC2', 'KIR2DL3','KIR2DL1', 'KIR3DL1', 
    # ---> ieILC1
    #'S100A13', 'NCR2', 'ENTPD1', 'TMIGD2',
    # ---> CREM
    #'LMNA', 'IFNGR1', 'CREM', 'SATB1',
    # ---> XCL1
    'KIR3DL1','KIR2DL1','KIR3DL2','KLRC2','KLRC3','CD3E','PRDM1','CD2','TMEM173'
  )
  ))

#Srt.obj modification

# Extract relevant cell barcodes for the parent aggr NK Srt Obj
cell.barcodes.CD16pos_adaptive <- rownames(srt.objs.list[['CD16pos']]@meta.data[srt.objs.list[['CD16pos']]@meta.data[clust.labs[['CD16pos']]] == names(clusts.defs[['CD16pos']][clusts.defs[['CD16pos']] == 'Adaptive']),])
cell.barcodes.CD16pos_non_adaptive <- setdiff(rownames(srt.objs.list[['CD16pos']]@meta.data), cell.barcodes.CD16pos_adaptive)

cell.barcodes.CD16neg_CD56_bright <- rownames(srt.objs.list[['CD16neg']]@meta.data[srt.objs.list[['CD16neg']]@meta.data[clust.labs[['CD16neg']]] == names(clusts.defs[['CD16neg']][clusts.defs[['CD16neg']] == 'CD56_bright']),])
cell.barcodes.CD16neg_CD160 <- rownames(srt.objs.list[['CD16neg']]@meta.data[srt.objs.list[['CD16neg']]@meta.data[clust.labs[['CD16neg']]] == names(clusts.defs[['CD16neg']][clusts.defs[['CD16neg']] == 'CD160']),])
cell.barcodes.CD16neg_ILC3 <- rownames(srt.objs.list[['CD16neg']]@meta.data[srt.objs.list[['CD16neg']]@meta.data[clust.labs[['CD16neg']]] == names(clusts.defs[['CD16neg']][clusts.defs[['CD16neg']] == 'ILC3']),])

cell.barcodes.Adaptive_trNK <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'Adaptive_trNK']),])
cell.barcodes.ENTPD1 <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'ENTPD1']),])
cell.barcodes.Cytotoxic <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'Cytotoxic']),])
cell.barcodes.CREM <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'CREM']),])
cell.barcodes.TCF7 <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'TCF7']),])

tmp.srt.obj@meta.data$barcode <- rownames(tmp.srt.obj@meta.data)

# Create a column defining the desired populations in the metadata
tmp.srt.obj@meta.data$adaptive_info <- 'Rest'

#aggr.Parent_NK.srt.obj_mod@meta.data$pop.info<-ifelse(aggr.Parent_NK.srt.obj_mod@meta.data$barcode %in% cell.barcodes.CD16neg_TR.Clust5, "CD16neg_TR.Clust5", "Rest")
tmp.srt.obj@meta.data$adaptive_info <-ifelse(tmp.srt.obj@meta.data$barcode %in% cell.barcodes.CD16pos_adaptive, "Adaptive", tmp.srt.obj@meta.data$adaptive_info)
tmp.srt.obj@meta.data$adaptive_info <-ifelse(tmp.srt.obj@meta.data$barcode %in% cell.barcodes.CD16pos_non_adaptive, "CD16pos_non_adaptive", tmp.srt.obj@meta.data$adaptive_info)

#tmp.srt.obj@meta.data$adaptive_info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.CD16neg_CD56_bright, "CD56_bright", tmp.srt.obj@meta.data$adaptive_info)
#tmp.srt.obj@meta.data$adaptive_info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.CD16neg_CD160, "CD160", tmp.srt.obj@meta.data$adaptive_info)
#tmp.srt.obj@meta.data$adaptive_info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.CD16neg_ILC3, "ILC3", tmp.srt.obj@meta.data$adaptive_info)

tmp.srt.obj@meta.data$adaptive_info <-ifelse(tmp.srt.obj@meta.data$barcode %in% cell.barcodes.Adaptive_trNK, "Adaptive_trNK", tmp.srt.obj@meta.data$adaptive_info)
tmp.srt.obj@meta.data$adaptive_info <-ifelse(tmp.srt.obj@meta.data$barcode %in% cell.barcodes.ENTPD1, "ENTPD1", tmp.srt.obj@meta.data$adaptive_info)
tmp.srt.obj@meta.data$adaptive_info <-ifelse(tmp.srt.obj@meta.data$barcode %in% cell.barcodes.Cytotoxic, "Cytotoxic", tmp.srt.obj@meta.data$adaptive_info)
tmp.srt.obj@meta.data$adaptive_info <-ifelse(tmp.srt.obj@meta.data$barcode %in% cell.barcodes.CREM, "CREM", tmp.srt.obj@meta.data$adaptive_info)
tmp.srt.obj@meta.data$adaptive_info <-ifelse(tmp.srt.obj@meta.data$barcode %in% cell.barcodes.TCF7, "TCF7", tmp.srt.obj@meta.data$adaptive_info)



# Define clusters' order.
these.pops <- list(
  `aggre`=c('Adaptive_trNK', 'ENTPD1', 'Cytotoxic', 'CREM', 'TCF7', 'Adaptive')
)

# Process per lineage.

#if(t.subset == "aggre"){
#  next
#}
# Plot
tmp.markers <- lapply(X = these.markers, FUN = function(tmp.feats){
  tmp.ids <- feature.info[feature.info$name %in% tmp.feats,]
  tmp.ids <- tmp.ids[match(tmp.feats, tmp.ids$name), ]
  tmp.ids <- tmp.ids$id
  tmp.ids <- rev(tmp.ids)
  return(tmp.ids)
})
tmp.gene.names <- lapply(X = these.markers, FUN = function(tmp.feats){
  tmp.names <- feature.info[feature.info$name %in% tmp.feats,]
  tmp.names <- tmp.names[match(tmp.feats, tmp.names$name), ]
  tmp.names <- tmp.names$name
  tmp.names <- rev(tmp.names)
  return(tmp.names)
})

named.labels <-  setNames(tmp.gene.names[[t.subset]], tmp.markers[[t.subset]])


#subset ser.obj
tmp.srt.obj <- subset(tmp.srt.obj, subset = tissue.tag == 'tumor.lung')
tmp.srt.obj <- subset(tmp.srt.obj, subset = adaptive_info %in% c('Adaptive_trNK', 'ENTPD1', 'Cytotoxic', 'CREM', 'TCF7', 'Adaptive'
                                                                 ))

tmp.ggplot <- dot.plot(
  seurat.obj=tmp.srt.obj, features=tmp.markers[[t.subset]], slot='data', ensembl=FALSE,
  # do.norm=FALSE,
  do.norm=TRUE,
  groups.tag='adaptive_info', groups.order=these.pops[[t.subset]], groups.of.int=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=FALSE, na.rm=TRUE, feature.thold=NULL,
  # this.color.scale=signatures.col.scale,
  this.color.scale='hot.and.cold',
  col.min=NULL, col.max=NULL,
  scale.by='radius', dot.scale=12,
  size.min=0, size.max=NA,
  file.name=NULL,
  named.labels = named.labels
)

tmp.ggplot <- tmp.ggplot + theme(legend.position='bottom')
tmp.ggplot <- tmp.ggplot 
# Output.
#tmp.width <- (length(these.pops[[t.subset]]) * 0.5) + 0.3
tmp.height <- (length(these.markers[[t.subset]]) * 0.4) + 0.3
# tmp.width <- 5
#tmp.height <- (length(these.markers[[t.subset]]) * 0.4) + 0.3
tmp.width <- (length(these.pops[[t.subset]]) * 0.4) + 0.3
tmp.lab <- paste0(obj.extended.names[t.subset], '_trNK_Adaptive')
publish.plot(
  tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf', blank.comp=blank.complement.2, do.legend=TRUE,
  width=tmp.width, height=tmp.height, c.width=tmp.width+2, c.height=tmp.height, legend.height=5, legend.width=0.5
)



############Figure 4-----------------------------------------    ############
### -----------------------Figure 4 ----------------------- ###
###########Figure 4-----------------------------------------    #############
#Module Vln Plot
tmp.reports.path <- paste0(fig.phen.comps.path, '/module_scoring')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

srt.obj <- srt.objs.list[['aggre']]
# Extract relevant cell barcodes for the parent aggr NK Srt Obj
cell.barcodes.CD16pos_adaptive <- rownames(srt.objs.list[['CD16pos']]@meta.data[srt.objs.list[['CD16pos']]@meta.data[clust.labs[['CD16pos']]] == names(clusts.defs[['CD16pos']][clusts.defs[['CD16pos']] == 'Adaptive']),])
cell.barcodes.CD16pos_non_adaptive <- setdiff(rownames(srt.objs.list[['CD16pos']]@meta.data), cell.barcodes.CD16pos_adaptive)

cell.barcodes.CD16neg_CD56_bright <- rownames(srt.objs.list[['CD16neg']]@meta.data[srt.objs.list[['CD16neg']]@meta.data[clust.labs[['CD16neg']]] == names(clusts.defs[['CD16neg']][clusts.defs[['CD16neg']] == 'CD56_bright']),])
cell.barcodes.CD16neg_CD160 <- rownames(srt.objs.list[['CD16neg']]@meta.data[srt.objs.list[['CD16neg']]@meta.data[clust.labs[['CD16neg']]] == names(clusts.defs[['CD16neg']][clusts.defs[['CD16neg']] == 'CD160']),])
cell.barcodes.CD16neg_ILC3 <- rownames(srt.objs.list[['CD16neg']]@meta.data[srt.objs.list[['CD16neg']]@meta.data[clust.labs[['CD16neg']]] == names(clusts.defs[['CD16neg']][clusts.defs[['CD16neg']] == 'ILC3']),])

cell.barcodes.Adaptive_trNK <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'Adaptive_trNK']),])
cell.barcodes.ENTPD1 <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'ENTPD1']),])
cell.barcodes.Cytotoxic <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'Cytotoxic']),])
cell.barcodes.CREM <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'CREM']),])
cell.barcodes.TCF7 <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'TCF7']),])

srt.obj@meta.data$barcode <- rownames(srt.obj@meta.data)
srt.obj@meta.data$barcode
# Create a column defining the desired populations in the metadata
srt.obj@meta.data$pop.info <- 'Rest'

#aggr.Parent_NK.srt.obj_mod@meta.data$pop.info<-ifelse(aggr.Parent_NK.srt.obj_mod@meta.data$barcode %in% cell.barcodes.CD16neg_TR.Clust5, "CD16neg_TR.Clust5", "Rest")
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.CD16pos_adaptive, "Adaptive", srt.obj@meta.data$pop.info)
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.CD16pos_non_adaptive, "CD16pos_non_adaptive", srt.obj@meta.data$pop.info)

srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.CD16neg_CD56_bright, "CD56_bright", srt.obj@meta.data$pop.info)
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.CD16neg_CD160, "CD160", srt.obj@meta.data$pop.info)
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.CD16neg_ILC3, "ILC3", srt.obj@meta.data$pop.info)

srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.Adaptive_trNK, "Adaptive_trNK", srt.obj@meta.data$pop.info)
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.ENTPD1, "ENTPD1", srt.obj@meta.data$pop.info)
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.Cytotoxic, "Cytotoxic", srt.obj@meta.data$pop.info)
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.CREM, "CREM", srt.obj@meta.data$pop.info)
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.TCF7, "TCF7", srt.obj@meta.data$pop.info)

plot_df <- srt.obj@meta.data


fill_tbl <- plot_df %>%
  group_by(pop.info) %>%
  summarise(
    MedianScore = median(Cytotoxic.score, na.rm = TRUE),
    .groups = "drop"
  )

plot_df <- left_join(plot_df, fill_tbl, by = "pop.info")

# Define the plot object using a pipe-friendly structure
tmp.ggplot <- plot_df %>% 
  # Exclude the 'Rest' population to focus on annotated clusters
  filter(pop.info != 'Rest') %>% 
  # Map aesthetics: sort X-axis by MedianScore in descending order
  ggplot(aes(x = forcats::fct_reorder(pop.info, MedianScore, .desc = TRUE), 
             y = Cytotoxic.score, 
             fill = MedianScore)) +
  # Violin plot to show distribution density
  geom_violin(alpha = 0.85, trim = FALSE, adjust = 1, width = 1.0) +
  # Thin boxplot to show quartiles and median (outliers hidden for clarity)
  geom_boxplot(width = 0.07, alpha = 0.7, outlier.shape = NA, fill = "white") +
  # Optional: Jittered points (commented out to handle large datasets/reduce file size)
  # geom_jitter(width = 0.2, shape = 21, size = 0.5, stroke = 0.1, 
  #             fill = "white", color = "black", alpha = 0.4) +
  # Apply the custom color gradient provided in the project settings
  scale_fill_gradientn(colours = c("#ffdf32", "#ff9a00", "#ff5a00", "#ff5719", "#EE0000", "#b30000", "#670000"), 
                       name = "Median\nCytotoxic Score") +
  # Label axes and title
  labs(x = "Populations (Sorted by Score)", 
       y = "Cytotoxic Score", 
       title = "Cytotoxic Score Distribution across pop.info") +
  # Use classic theme for publication quality
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# Export the plot using the project's custom publish function
tmp.lab <- 'pop_info_Cytotoxic_score_vln'
publish.plot(
  tmp.ggplot = tmp.ggplot, 
  output.path = tmp.reports.path, 
  file.name = tmp.lab, 
  type = 'pdf', 
  blank.comp = blank.complement.2, # Follows project's naming conventions
  do.legend = TRUE,
  width = 8, height = 5
)


############DEG Vocano Cyt high vs low    -----------------------------------------    ############
### ----------------------- DEG Vocano ----------------------- ###
############DEG Vocano Cyt high vs low    -----------------------------------------    ############

############FGSEA-----------------------------------------    ############
### ----------------------- FGSEA ----------------------- ###
###########FGSEA-----------------------------------------    #############


############GSVA Correlation-----------------------------------------    ############
### -----------------------GSVA Correlation ----------------------- ###
###########GSVA Correlation-----------------------------------------    #############

tmp.reports.path <- paste0(fig.phen.comps.path, '/gsva')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
#Get exp data ────────────────────────────────────────
t.subset         <- "trNK"          # data set
seurat.obj       <- srt.objs.list[[t.subset]]    # 
tissue.conditions <- c("tumor.lung", "normal.lung")
gene.stats.lst <- list()

for (cond in tissue.conditions) {
  message("Processing", "/ ", cond)
  res <- process_condition(
    seurat.obj      = seurat.obj,
    tissue.condition = cond
  )
  if (!is.null(res)) {
    gene.stats.lst[[paste0(cond)]] <- res
  }
}

# Gene symbol -> Ensembl ID
gene.sets.ens <- lapply(modules.feats, function(sym.vec) {
  unique(feature.info$id[feature.info$name %in% sym.vec])
})
stopifnot(all(lengths(gene.sets.ens) > 0))

#Gather the matrix
combined.expr <- rbindlist(gene.stats.lst, idcol = 'source', fill = TRUE)
combined.expr <- as.data.frame(combined.expr)
combined.expr[1:5, 1:5]
#Inset sample ID
combined.expr$sample <- paste0(combined.expr$source, '_', combined.expr$donor)
rownames(combined.expr) <- combined.expr$sample

#Split the data into metadata and exp. matrix
expr.meta <- combined.expr[, c('donor', 'source', 'sample')]
expr.mat <- combined.expr[, !(colnames(combined.expr) %in% c('donor', 'source', 'sample'))]

#Transpose
expr.mat <- t(as.matrix(expr.mat))

#GSVA
set.seed(123)       
gsva.res <- gsva(expr.mat,
                 gset.idx.list = gene.sets.ens,
                 method        = 'gsva',
                 kcdf          = 'Gaussian',
                 min.sz        = 5,
                 max.sz        = 500,
                 parallel      = 1)

#write the result
gsva.df <- as.data.frame(t(gsva.res), check.names = F)
gsva.df$sample <- rownames(gsva.df)

gsva.df <- merge(gsva.df, expr.meta, by = "sample")

gsva.df %>% 
  pivot_wider(names_from = source, values_from = c(Cytotoxic, IL21,trNK ) ) %>% 
  write.csv(paste0(tmp.reports.path, "/wider_gsva_res.csv"))

gsva.df.tumor <- gsva.df[gsva.df$source == 'tumor.lung',]


tmp.df <- gsva.df.tumor
tmp.ggplot <- tmp.df %>%
  ggplot(aes_string(x = 'IL21', y = 'Cytotoxic')) + 
  #geom_point(shape=1, stroke=5, size=10, alpha=0.6,color  = colors[[cluster]])+
  geom_point(shape=1, stroke=5, size=10, alpha=0.6)+
  geom_point(shape=1, stroke=5, size=10, alpha=0.6)+
  geom_smooth(method="gam", formula=y~x, se=TRUE, fullrange=TRUE, linewidth=2.5, color  = "black") +
  scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
  theme_classic()
tmp.ggplot
tmp.lab <- paste0('/', 'Cyto_vs_IL21sigmouse','_','tumor')
publish.plot(
  tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
  stat.cor=TRUE, 
  #cor.group='pop',
  blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=14, height=14
)
cor.test(tmp.df$Cytotoxic,
         tmp.df$IL21,
         method = "spearman")



tmp.ggplot <- tmp.df %>%
  ggplot(aes_string(x = 'IL21', y = 'Cytotoxic')) + 
  #geom_point(shape=1, stroke=5, size=10, alpha=0.6,color  = colors[[cluster]])+
  geom_point(shape=1, stroke=5, size=10, alpha=0.6)+
  geom_point(shape=1, stroke=5, size=10, alpha=0.6)+
  geom_smooth(method="gam", formula=y~x, se=TRUE, fullrange=TRUE, linewidth=2.5, color  = "black") +
  scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
  theme_classic()
tmp.ggplot
tmp.lab <- paste0('/', 'Cyto_vs_IL21sigmouse','_','tumor')
publish.plot(
  tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
  stat.cor=TRUE, 
  #cor.group='pop',
  blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=14, height=14
)
cor.test(tmp.df$Cytotoxic,
         tmp.df$IL21,
         method = "spearman")

tmp.df 
tmp.ggplot <- tmp.df %>%
  ggplot(aes_string(x = 'IL21.bulk.109', y = 'Cytotoxic')) + 
  #geom_point(shape=1, stroke=5, size=10, alpha=0.6,color  = colors[[cluster]])+
  geom_point(shape=1, stroke=5, size=10, alpha=0.6)+
  geom_point(shape=1, stroke=5, size=10, alpha=0.6)+
  geom_smooth(method="gam", formula=y~x, se=TRUE, fullrange=TRUE, linewidth=2.5, color  = "black") +
  scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
  theme_classic()
tmp.ggplot
tmp.lab <- paste0('/', 'Cyto_vs_IL21human','_','tumor')
publish.plot(
  tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
  stat.cor=TRUE, 
  #cor.group='pop',
  blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=14, height=14
)
cor.test(tmp.df$Cytotoxic,
         tmp.df$IL21.bulk.109,
         method = "spearman")



tmp.ggplot <- tmp.df %>%
  ggplot(aes_string(x = 'IL21.bulk.109', y = 'NK3')) + 
  #geom_point(shape=1, stroke=5, size=10, alpha=0.6,color  = colors[[cluster]])+
  geom_point(shape=1, stroke=5, size=10, alpha=0.6)+
  geom_point(shape=1, stroke=5, size=10, alpha=0.6)+
  geom_smooth(method="gam", formula=y~x, se=TRUE, fullrange=TRUE, linewidth=2.5, color  = "black") +
  scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
  theme_classic()
tmp.ggplot
tmp.lab <- paste0('/', 'NK3_vs_IL21human','_','tumor')
publish.plot(
  tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
  stat.cor=TRUE, 
  #cor.group='pop',
  blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=14, height=14
)
cor.test(tmp.df$NK3,
         tmp.df$IL21.bulk.109,
         method = "spearman")



#####IL21 response vs trNK percentage
NK_subset <- read.csv(paste0(prj.gen.path, '/Tumor_NK_raw.csv'))
NK_subset <- NK_subset %>% 
  mutate(trNK_pct = trNK*100)

gsva.df.tumor$donor.id.tag <- as.numeric(gsva.df.tumor$donor)
tmp.df <- gsva.df.tumor %>% 
  left_join(NK_subset, by= 'donor.id.tag')

tmp.ggplot <- tmp.df %>%
  ggplot(aes_string(x = 'IL21.bulk.109', y = 'trNK_pct')) + 
  #geom_point(shape=1, stroke=5, size=10, alpha=0.6,color  = colors[[cluster]])+
  geom_point(shape=1, stroke=5, size=10, alpha=0.6)+
  geom_point(shape=1, stroke=5, size=10, alpha=0.6)+
  geom_smooth(method="gam", formula=y~x, se=TRUE, fullrange=TRUE, linewidth=2.5, color  = "black") +
  scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
  theme_classic()
tmp.ggplot
tmp.lab <- paste0('/', 'trNK_pct_vs_IL21human','_','tumor')
publish.plot(
  tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
  stat.cor=TRUE, 
  #cor.group='pop',
  blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=14, height=14
)
cor.test(tmp.df$trNK_pct,
         tmp.df$IL21.bulk.109,
         method = "spearman")


####Save the source data
tmp.df %>% 
  write.csv(paste0(tmp.reports.path, '/correlation_GSVA.csv'))

############Figure 6-----------------------------------------    ############
### -----------------------Figure 6 ----------------------- ###
############Figure 6-----------------------------------------    #############
############correlation.-----------------------------------------    ############
###--------------------------- correlation.--------------------------- ###
############correlation-----------------------------------------    ############
tmp.reports.path <- paste0(fig.phen.comps.path, '/cor')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

###IL21 expression in CD4 vs Cytotoxic score in trNK-------------------------------------
gene <- c('IL21')
gene.ens <- lapply(gene, function(sym.vec) {
  unique(feature.info$id[feature.info$name %in% sym.vec])
}) %>% unlist()


condition <- 'tumor.lung'
tmp.df <- gsva.df[gsva.df$source == paste0('tumor.lung'),]
tmp.df <- tmp.df %>% 
  inner_join(CD4.exp[, c("donor", gene.ens)], by = 'donor')
tmp.ggplot <- tmp.df %>%
  ggplot(aes_string(x = gene.ens, y = 'Cytotoxic')) + 
  geom_point(shape=1, stroke=5, size=10, alpha=0.6)+
  geom_smooth(method="gam", formula=y~x, se=TRUE, fullrange=TRUE, linewidth=2.5, color  = "black") +
  scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
  theme_classic()
tmp.lab <- paste0('/', 'Cyto_vs_IL21expInCD4',condition)

publish.plot(
  tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
  stat.cor=TRUE, 
  #cor.group='pop',
  blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=14, height=14
)
print(cor.test(tmp.df[['Cytotoxic']],tmp.df[[gene.ens]],method = 'spearman' ))
###Save the source data
tmp.df %>% 
  write.csv(paste0(tmp.reports.path, '/correlation_GSVA_IL21inCD4exp.csv'))


###IL21 expression in CD4 vs Cytotoxic score in Adaptive trNK-------------------------------------
###IL21 expression in CD4 vs Cytotoxic score in Adaptive trNK, Cytotoxic trNK-------------------------------------
gene <- c('IL21')
gene.ens <- lapply(gene, function(sym.vec) {
  unique(feature.info$id[feature.info$name %in% sym.vec])
}) %>% unlist()


#Get exp data ────────────────────────────────────────
t.subset         <- "trNK"          # data set
cluster.col      <- clust.labs[[t.subset]]       # 
seurat.obj       <- srt.objs.list[[t.subset]]    # 
cluster.ids      <- sort(unique(seurat.obj@meta.data[[cluster.col]]))
tissue.conditions <- c("tumor.lung", "normal.lung")

gene.stats.lst <- list()
for (cl in cluster.ids) {
  for (cond in tissue.conditions) {
    message("Processing cluster ", cl, " / ", cond)
    res <- process_cluster_condition(
      seurat.obj      = seurat.obj,
      cluster.id      = cl,
      tissue.condition = cond,
      cluster.col     = cluster.col
    )
    if (!is.null(res)) {
      gene.stats.lst[[paste0("cluster", cl, "_", cond)]] <- res
    }
  }
}

# Gene symbol -> Ensembl ID
gene.sets.ens <- lapply(modules.feats, function(sym.vec) {
  unique(feature.info$id[feature.info$name %in% sym.vec])
})
stopifnot(all(lengths(gene.sets.ens) > 0))  

#Gather the matrix
combined.expr <- rbindlist(gene.stats.lst, idcol = 'source', fill = TRUE)
combined.expr <- as.data.frame(combined.expr)
combined.expr[1:5, 1:5]
#Inset sample ID
combined.expr$sample <- paste0(combined.expr$source, '_', combined.expr$donor)
rownames(combined.expr) <- combined.expr$sample

#Split the data into metadata and exp. matrix
expr.meta <- combined.expr[, c('donor', 'source', 'sample')]
expr.mat <- combined.expr[, !(colnames(combined.expr) %in% c('donor', 'source', 'sample'))]

#Transpose
expr.mat <- t(as.matrix(expr.mat))

#GSVA
set.seed(123)       
gsva.res <- gsva(expr.mat,
                 gset.idx.list = gene.sets.ens,
                 method        = 'gsva',
                 kcdf          = 'Gaussian',
                 min.sz        = 5,
                 max.sz        = 500,
                 parallel      = 1)

#write the result
gsva.df_2 <- as.data.frame(t(gsva.res), check.names = F)
gsva.df_2$sample <- rownames(gsva.df_2)

gsva.df_2 <- gsva.df_2 %>%
  mutate(
    cluster = str_extract(sample, "^cluster\\d+")
  )

cluster_map <- c(
  "cluster0" = "Adaptive_trNK",
  "cluster1" = "ENTPD1",
  "cluster2" = "Cytotoxic",
  "cluster3" = "CREM",
  "cluster4" = "TCF7"
)

gsva.df_2 <- gsva.df_2 %>%
  mutate(
    cluster_annot = cluster_map[cluster]
  )

gsva.df_2 <- merge(gsva.df_2, expr.meta, by = "sample")
write.csv(gsva.df_2, paste0(tmp.reports.path, "/gsva_res.csv"))

gsva.df_2 %>% 
  pivot_wider(names_from = source, values_from = c(Cytotoxic, IL21,trNK ) ) %>% 
  write.csv(paste0(tmp.reports.path, "/wider_gsva_res.csv"))

###Vizualize
condition <- 'tumor.lung'
clusters <- c('cluster0', 'cluster2')
colors <- list(
  `cluster0` = '#5D72E0', #Adaptive trNK
  `cluster2` = '#D4824B' #Cytotoxic
)
cluster <- clusters[1] 
#Make a key
gsva.df_2$sample_id <- gsva.df_2$sample
gsva.df_2$donor <- as.character(sub(".*_", "", gsva.df_2$sample_id))
gsva.df_2$source <- sub("_[0-9]+$", "", gsva.df_2$sample)

#for(cluster in clusters){
tmp.df <- gsva.df_2[gsva.df_2$source == paste0(cluster, '_', condition),]

tmp.df <- tmp.df %>% 
  inner_join(CD4.exp[, c("donor", gene.ens)], by = 'donor')
tmp.ggplot <- tmp.df %>%
  ggplot(aes_string(x = gene.ens, y = 'Cytotoxic')) + 
  geom_point(shape=1, stroke=5, size=10, alpha=0.6,color  = colors[[cluster]])+
  geom_smooth(method="gam", formula=y~x, se=TRUE, fullrange=TRUE, linewidth=2.5, color  = "black") +
  scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
  theme_classic()
tmp.lab <- paste0('/', 'Cyto_vs_IL21expInCD4',cluster,condition)
publish.plot(
  tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
  stat.cor=TRUE, 
  #cor.group='pop',
  blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=14, height=14
)
print(cor.test(tmp.df[['Cytotoxic']],tmp.df[[gene.ens]],method = 'spearman' ))

#}

###Save the source data
tmp.df %>% 
  write.csv(paste0(tmp.reports.path, '/trNK_0_correlation_GSVA_IL21inCD4exp.csv'))


############Figure S1-----------------------------------------    ############
### -----------------------Figure S1 ----------------------- ###
############Figure S1-----------------------------------------    #############

############trNK freq in the comparison sereral parameters-----------------------------------------    ############
### ---------------------- trNK freq in the comparison sereral parameters- ----------------------- ###
############trNK freq in the comparison sereral parameters------------------------------------------    ############
tmp.reports.path <- paste0(fig.phen.comps.path, '/boxpot')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

meta.cols <- list(
  `Gender`= c(
    'F'='#BE82BF',
    'M'='#8FD5C9'),
  `age_tag` =c(
    '<66'='#93F6FF',
    '66-70'='#91D2F3',
    '71-75'='#99C1FF',
    '>75'='#919EF3'),
  `Staging` = c(
    '1a2' = '#F7FBFF', 
    '1a3' = '#C6DBEF', 
    '1b'  = '#4292C6', 
    '2a'  = '#6BAED6', 
    '2b'  = '#756BB1', 
    '3a'  = '#54278F', 
    '3b'  = '#3F007D')
  #`til.score` = c(
  #  'TIL HIGH' = '#B22222',
  #  'TIL LOW' = '#F08080',
  #  'TIL MOD' = 'lightgray'
  #)
)


NK_subset <- read.csv(paste0(prj.gen.path, '/Tumor_NK_raw.csv'))
NK_subset <- NK_subset %>% 
  mutate(trNK_pct = trNK*100)

metadata <- metadata %>% 
  mutate(age_tag = ifelse(age<66,
                          yes='<66',
                          no=ifelse(
                            age<71,
                            yes='66-70',
                            no=ifelse(
                              test=age<76,
                              yes='71-75',
                              no='>75'
                            )
                          )
  )
  )


metadata <- metadata %>%
  mutate(Gender = factor(Gender, levels = c("M", "F"))) %>% 
  mutate(Staging = factor(Staging, levels = c("1a2", "1a3", "1b", "2a", "2b", "3a", "3b", "4a"))) %>% 
  #mutate(til.score = factor(til.score, levels = c('TIL LOW', 'TIL MOD', 'TIL HIGH'))) %>% 
  mutate(age_tag = factor(age_tag, levels = c("<66", "66-70", "71-75", ">75"))) 

tmp.tags <- c('age_tag', 'Staging', 'Gender')

for (tmp.tag in tmp.tags){
  df <- NK_subset %>% 
    left_join(metadata) %>% 
    filter(!is.na(!!rlang::sym(tmp.tag)))
  res <- kruskal.test(reformulate(tmp.tag, response = "trNK_pct"), data = df)
  print(tmp.tag)
  print(res)
  tmp.ggplot <- NK_subset %>% 
    left_join(metadata) %>% 
    filter(!is.na(!!rlang::sym(tmp.tag))) %>%
    ggplot(aes_string(x = tmp.tag ,'trNK_pct')) +
    geom_boxplot(aes_string(color=tmp.tag ), width=0.7, linewidth=4, fatten=4, outlier.shape=NA)+
    geom_jitter(shape=1, stroke=5, size=10, color='black', width=0)+
    scale_color_manual(values=meta.cols[[tmp.tag]])+
    scale_y_continuous(breaks=scales::pretty_breaks(n=3))+
    labs(x=tmp.tag, y='Number of cells', color='') 
  tmp.lab <- paste0(
    '/', tmp.tag, '_CellCount_CellType_Donor')
  if (tmp.tag == 'Staging'){
    publish.plot(
      tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
      blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=10, height=10
    )
  }else{
    publish.plot(
      tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
      blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=7, height=10
    )
  }
}

tmp.df <- metadata.trNK %>% 
  left_join(metadata.CD8, by ='donor.id.tag') %>% 
  left_join(metadata.CD45, by = 'donor.id.tag')
tmp.df
tmp.ggplot <- tmp.df %>%
  ggplot(aes_string(x = 'CD3_CD8_meancells_per_mm2', y = 'trNK')) + 
  #geom_point(shape=1, stroke=5, size=10, alpha=0.6,color  = colors[[cluster]])+
  geom_point(shape=1, stroke=5, size=10, alpha=0.6)+
  geom_point(shape=1, stroke=5, size=10, alpha=0.6)+
  geom_smooth(method="gam", formula=y~x, se=TRUE, fullrange=TRUE, linewidth=2.5, color  = "black") +
  scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
  theme_classic()
tmp.ggplot
tmp.lab <- paste0('/', 'cd8tilvstrNK','_','tumor')
publish.plot(
  tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
  stat.cor=TRUE, 
  #cor.group='pop',
  blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=14, height=14
)


tmp.df <- metadata.trNK %>% 
  left_join(metadata.CD8, by ='donor.id.tag') %>% 
  left_join(metadata.CD45, by = 'donor.id.tag')
tmp.df
tmp.ggplot <- tmp.df %>%
  ggplot(aes_string(x = 'mean_CD45permm', y = 'trNK')) + 
  #geom_point(shape=1, stroke=5, size=10, alpha=0.6,color  = colors[[cluster]])+
  geom_point(shape=1, stroke=5, size=10, alpha=0.6)+
  geom_point(shape=1, stroke=5, size=10, alpha=0.6)+
  geom_smooth(method="gam", formula=y~x, se=TRUE, fullrange=TRUE, linewidth=2.5, color  = "black") +
  scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
  theme_classic()
tmp.ggplot
tmp.lab <- paste0('/', 'cd45tilvstrNK','_','tumor')
publish.plot(
  tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
  stat.cor=TRUE, 
  #cor.group='pop',
  blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=14, height=14
)


############Figure S2-----------------------------------------    ############
### -----------------------Figure S2 ----------------------- ###
############Figure S2-----------------------------------------    #############

###########Dot plot-----------------------------------------    ############
### ----------------------- Dot plot ----------------------- ###
###########Dot plot-----------------------------------------    ############



tmp.reports.path <- paste0(fig.phen.comps.path, '/dotplot')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

# Define genes of interest.
these.markers <- list(
  `CD16pos`=c(
    # ---> Adaptive
    'KLRC2', 'CD3E',
    # ---> TNFRSF9
    'TNFRSF9', 'CCL4L2',
    # ---> CX3CR1
    'CX3CR1', 'PTGDS',
    # ---> EGR1
    'EGR1',
    # ---> AREG
    'AREG', 'DUSP2'
  ),
  `CD16neg`=c(
    # ---> CD160
    #'CD160', 'CCL3',
    # ---> ILC3
    #'KIT', 'AHR', 'IL7R',
    # ---> Circ
    #'SELL', 'S1PR5', 'KLF2',
    # ---> Tissue resident
    'ITGAE', 'ITGA1', 'ZNF683'
  ),
  `trNK`=c(
    # ---> GZMK
    'GZMK', 'CCL3', 'CCL4',
    # ---> Adaptive tr
    #'KLRC2', 'KIR2DL3','KIR2DL1', 'KIR3DL1', 
    # ---> ieILC1
    'S100A13', 'NCR2', 'ENTPD1', 'TMIGD2',
    # ---> CREM
    'LMNA', 'IFNGR1', 'CREM', 'SATB1'
    # ---> XCL1
    #'IL7R', 'XCL1', 'TCF7', 'LTB'
  ),
  `aggre`=c(
    'ITGAE', 'ITGA1', 'ZNF683'
    # ---> GZMK
    #'GZMK', 'CCL3', 'CCL4',
    # ---> Adaptive tr
    #'KLRC2', 'KIR2DL3','KIR2DL1', 'KIR3DL1', 
    # ---> ieILC1
    #'S100A13', 'NCR2', 'ENTPD1', 'TMIGD2',
    # ---> CREM
    #'LMNA', 'IFNGR1', 'CREM', 'SATB1',
    # ---> XCL1
    #'IL7R', 'XCL1', 'TCF7', 'LTB'
  )
)

# Define clusters' order.
these.pops <- list(
  `CD16pos`=c('4', '5', '0', '2', '3'),
  `CD16neg`=c('4', '6', '3', '0'),
  `trNK`=c('0', '1', '2', '3', '4'),
  `aggre`=c('0', '1', '2', '3', '4', '5', '6', '7', '8')
)

#DotPlot Cytotoxicity related genes tumor only
# Define genes of interest.
t.subset <- 'trNK'
these.markers <- list(
  `trNK`=c(
    # ---> GZMK
    'GZMK', 'CCL3', 'CCL4',
    # ---> Adaptive tr
    #'KLRC2', 'KIR2DL3','KIR2DL1', 'KIR3DL1', 
    # ---> ieILC1
    'S100A13', 'NCR2', 'ENTPD1', 'TMIGD2',
    # ---> CREM
    'LMNA', 'IFNGR1', 'CREM', 'SATB1'
    # ---> XCL1
    #'IL7R', 'XCL1', 'TCF7', 'LTB'
  )
)

# Process per lineage.
# Plot

tissue.type <- 'tumor.lung'
tmp.markers <- lapply(X = these.markers, FUN = function(tmp.feats){
  tmp.ids <- feature.info[feature.info$name %in% tmp.feats,]
  tmp.ids <- tmp.ids[match(tmp.feats, tmp.ids$name), ]
  tmp.ids <- tmp.ids$id
  tmp.ids <- rev(tmp.ids)
  return(tmp.ids)
})
tmp.gene.names <- lapply(X = these.markers, FUN = function(tmp.feats){
  tmp.names <- feature.info[feature.info$name %in% tmp.feats,]
  tmp.names <- tmp.names[match(tmp.feats, tmp.names$name), ]
  tmp.names <- tmp.names$name
  tmp.names <- rev(tmp.names)
  return(tmp.names)
})

named.labels <-  setNames(tmp.gene.names[[t.subset]], tmp.markers[[t.subset]])


tmp.ggplot <- dot.plot(
  seurat.obj=subset(srt.objs.list[[t.subset]], subset = tissue.tag == tissue.type), features=tmp.markers[[t.subset]], slot='data', ensembl=FALSE,
  # do.norm=FALSE,
  do.norm=TRUE,
  groups.tag=clust.labs[t.subset], groups.order=these.pops[[t.subset]], groups.of.int=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=FALSE, na.rm=TRUE, feature.thold=NULL,
  # this.color.scale=signatures.col.scale,
  this.color.scale='hot.and.cold',
  col.min=NULL, col.max=NULL,
  scale.by='radius', dot.scale=12,
  size.min=0, size.max=NA,
  file.name=NULL,
  named.labels = named.labels
)
tmp.ggplot 
tmp.ggplot <- tmp.ggplot + theme(legend.position='bottom')
tmp.ggplot <- tmp.ggplot 
# Output.
#tmp.width <- (length(these.pops[[t.subset]]) * 0.5) + 0.3
tmp.height<- (length(these.markers[[t.subset]]) * 0.4) + 0.3
# tmp.width <- 5
#tmp.height <- (length(these.markers[[t.subset]]) * 0.4) + 0.3
tmp.width <- (length(these.pops[[t.subset]]) * 0.4) + 0.3

tmp.lab <- paste0(obj.extended.names[t.subset],tissue.type,'trNK')

publish.plot(
  tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf', blank.comp=blank.complement.2, do.legend=TRUE,
  width=tmp.width, height=tmp.height, c.width=tmp.width+2, c.height=tmp.height, legend.height=5, legend.width=0.5)


t.subset <- 'trNK'
these.markers <- list(
  `trNK`=c(
    'PRF1', 'GNLY', 'GZMK', 'GZMM', 'GZMH', 'GZMB', 'GZMA', 'TNFSF10', 'FASLG'
  )
)

# Process per lineage.
# Plot

tissue.type <- 'tumor.lung'
tmp.markers <- lapply(X = these.markers, FUN = function(tmp.feats){
  tmp.ids <- feature.info[feature.info$name %in% tmp.feats,]
  tmp.ids <- tmp.ids[match(tmp.feats, tmp.ids$name), ]
  tmp.ids <- tmp.ids$id
  tmp.ids <- rev(tmp.ids)
  return(tmp.ids)
})
tmp.gene.names <- lapply(X = these.markers, FUN = function(tmp.feats){
  tmp.names <- feature.info[feature.info$name %in% tmp.feats,]
  tmp.names <- tmp.names[match(tmp.feats, tmp.names$name), ]
  tmp.names <- tmp.names$name
  tmp.names <- rev(tmp.names)
  return(tmp.names)
})

named.labels <-  setNames(tmp.gene.names[[t.subset]], tmp.markers[[t.subset]])


tmp.ggplot <- dot.plot(
  seurat.obj=subset(srt.objs.list[[t.subset]], subset = tissue.tag == tissue.type), features=tmp.markers[[t.subset]], slot='data', ensembl=FALSE,
  # do.norm=FALSE,
  do.norm=TRUE,
  groups.tag=clust.labs[t.subset], groups.order=these.pops[[t.subset]], groups.of.int=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=FALSE, na.rm=TRUE, feature.thold=NULL,
  # this.color.scale=signatures.col.scale,
  this.color.scale='hot.and.cold',
  col.min=NULL, col.max=NULL,
  scale.by='radius', dot.scale=12,
  size.min=0, size.max=NA,
  file.name=NULL,
  named.labels = named.labels
)
tmp.ggplot 
tmp.ggplot <- tmp.ggplot + theme(legend.position='bottom')
tmp.ggplot <- tmp.ggplot 
# Output.
#tmp.width <- (length(these.pops[[t.subset]]) * 0.5) + 0.3
tmp.height<- (length(these.markers[[t.subset]]) * 0.4) + 0.3
# tmp.width <- 5
#tmp.height <- (length(these.markers[[t.subset]]) * 0.4) + 0.3
tmp.width <- (length(these.pops[[t.subset]]) * 0.4) + 0.3

tmp.lab <- paste0(obj.extended.names[t.subset],tissue.type,'Cyt_trNK')

publish.plot(
  tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf', blank.comp=blank.complement.2, do.legend=TRUE,
  width=tmp.width, height=tmp.height, c.width=tmp.width+2, c.height=tmp.height, legend.height=5, legend.width=0.5)


############FGSEA-----------------------------------------    ############
### ----------------------- FGSEA ----------------------- ###
###########FGSEA-----------------------------------------    #############
tmp.reports.path <- paste0(fig.phen.comps.path, '/fgsea')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

#Do GSEA skin.trNK signature Immunity Volume 57, Issue 1, 9 January 2024, Pages 124-140.e7
####FGSEA aubtract setting
# @ Seeting definition
# Reference group
# Cluster definition pattern
# Signature

#In tumor trNK
tn.tag<- c('tumor.lung')
gsea.settings <- list(
  `trNK`=list(
    `skin_trNK_TCF7`=c('TCF7_trNKTCF7', '^TCF7$', 'skin.trNK'),
    `skin_trNK_trNKadaptive`=c('TCF7_trNKadaptive', '^Adaptive_trNK$', 'skin.trNK'),
    `skin_trNK_Cyt`=c('TCF7_trNKCyt', '^Cytotoxic$', 'skin.trNK'),
    `skin_trNK_CREM`=c('TCF7_trNKCREM', '^CREM$', 'skin.trNK'),
    `skin_trNK_ENTPD1`=c('TCF7_trNKENTPD1', '^ENTPD1$', 'skin.trNK')
  )
)

for(data.set in names(gsea.settings)){
  cat(paste0(obj.extended.names[data.set], '\n'))
  this.reports.path <- paste0(tmp.reports.path, '/', obj.extended.names[data.set])
  if(!dir.exists(this.reports.path)) dir.create(this.reports.path)
  for(tmp.eval in names(gsea.settings[[data.set]])){
    tmp.ref <- gsea.settings[[data.set]][[tmp.eval]][1]
    tmp.pttn <- gsea.settings[[data.set]][[tmp.eval]][2]
    tmp.mod <- gsea.settings[[data.set]][[tmp.eval]][3]
    srt.objs.list[[data.set]]@meta.data[, 'tmp.tag'] <- ifelse(
      test=str_detect(
        string=clusts.defs[[data.set]][
          as.character(srt.objs.list[[data.set]]@meta.data[, clust.labs[data.set]])
        ],
        pattern=tmp.pttn
      ),
      yes=tmp.ref,
      no='Rest'
    )
    do.fgsea(
      metrics.src=subset(srt.objs.list[[data.set]], subset = tissue.tag == tn.tag),
      modules.feats=modules.feats[tmp.mod], 
      metric='signal.to.noise', cell.subset=subset.list[[data.set]],
      tag.of.int='tmp.tag',
      output.path=this.reports.path,
      vals.to.depict=tmp.ref,
      files.preffix=paste0(obj.extended.names[data.set], tn.tag)
    )
    srt.objs.list[[data.set]]@meta.data[, 'tmp.tag'] <- NULL
  }
}
tmp.res <- list()
for(data.set in names(gsea.settings)){
  cat(paste0(obj.extended.names[data.set], '\n'))
  this.reports.path <- paste0(tmp.reports.path, '/', obj.extended.names[data.set])
  if(!dir.exists(this.reports.path)) dir.create(this.reports.path)
  for(tmp.eval in names(gsea.settings[[data.set]])){
    tmp.ref <- gsea.settings[[data.set]][[tmp.eval]][1]
    tmp.pttn <- gsea.settings[[data.set]][[tmp.eval]][2]
    tmp.mod <- gsea.settings[[data.set]][[tmp.eval]][3]
    srt.objs.list[[data.set]]@meta.data[, 'tmp.tag'] <- ifelse(
      test=str_detect(
        string=clusts.defs[[data.set]][
          as.character(srt.objs.list[[data.set]]@meta.data[, clust.labs[data.set]])
        ],
        pattern=tmp.pttn
      ),
      yes=tmp.ref,
      no='Rest'
    )
    tmp.fgse.res <- do.fgsea(
      metrics.src=subset(srt.objs.list[[data.set]], subset = tissue.tag == tn.tag),
      modules.feats=modules.feats[tmp.mod], 
      metric='signal.to.noise', cell.subset=subset.list[[data.set]],
      tag.of.int='tmp.tag',
      output.path=this.reports.path,
      vals.to.depict=tmp.ref,
      just.get.res = TRUE,
      files.preffix=paste0(obj.extended.names[data.set], tn.tag)
    )
    tmp.res[[tmp.ref]] <- tmp.fgse.res[[tmp.ref]]
    srt.objs.list[[data.set]]@meta.data[, 'tmp.tag'] <- NULL
  }
}

tmp.ggplot  <- rbindlist(tmp.res, idcol = "cluster") %>% 
  ggplot(aes(
    x = factor(cluster, levels = c('TCF7_trNKadaptive', 'TCF7_trNKENTPD1', 'TCF7_trNKCyt', 'TCF7_trNKCREM', 'TCF7_trNKTCF7')),
    y = pathway, 
    size = -log10(padj),
    fill = NES
  )) +
  geom_point(stroke = 0.7, shape = 21, col = 'black') +
  scale_fill_gradientn(
    colours = rev(c("#D4406A", "#D4406A", "white", "#A39FCA", "#2C388A")),
    values = scales::rescale(c(-2, -1, 0, 1, 2)),
    name = "NES"
  ) +
  labs(x = 'Groups', y = 'module', fill = 'NES', size = '-log10padj') +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()
  )
tmp.ggplot <- tmp.ggplot + theme(legend.position='bottom')
tmp.height <- ( 0.4) + 0.3
tmp.width <- ( 4*0.4) + 0.3
tmp.lab <- paste0( 'skin_trNK')
publish.plot(
  tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf', blank.comp=blank.complement.2, do.legend=TRUE,
  width=tmp.width, height=tmp.height, c.width=tmp.width+2, c.height=tmp.height, legend.height=5, legend.width=0.5
)

tmp.res %>% 
  write.csv(paste0(tmp.reports.path, '/GSEA_dotplot_Tcf7_skin.csv'))

###########fgsea dot plot-----------------------------------------    ############
### ----------------------- fgsea dot plot ----------------------- ###
###########fgsea dot plot-----------------------------------------    ############

####FGSEA aubtract setting
# @ Seeting definition
# Reference group
# Cluster definition pattern
# Signature
tmp.reports.path <- paste0(fig.phen.comps.path, '/fgsea')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

####FGSEA aubtract setting
# @ Seeting definition
# Reference group
# Cluster definition pattern
# Signature
# Extract relevant cell barcodes for the parent aggr NK Srt Obj
cell.barcodes.CD16pos_adaptive <- rownames(srt.objs.list[['CD16pos']]@meta.data[srt.objs.list[['CD16pos']]@meta.data[clust.labs[['CD16pos']]] == names(clusts.defs[['CD16pos']][clusts.defs[['CD16pos']] == 'Adaptive']),])
cell.barcodes.CD16pos_non_adaptive <- setdiff(rownames(srt.objs.list[['CD16pos']]@meta.data), cell.barcodes.CD16pos_adaptive)

#cell.barcodes.CD16neg_CD56_bright <- rownames(srt.objs.list[['CD16neg']]@meta.data[srt.objs.list[['CD16neg']]@meta.data[clust.labs[['CD16neg']]] == names(clusts.defs[['CD16neg']][clusts.defs[['CD16neg']] == 'CD56_bright']),])
#cell.barcodes.CD16neg_CD160 <- rownames(srt.objs.list[['CD16neg']]@meta.data[srt.objs.list[['CD16neg']]@meta.data[clust.labs[['CD16neg']]] == names(clusts.defs[['CD16neg']][clusts.defs[['CD16neg']] == 'CD160']),])
#cell.barcodes.CD16neg_ILC3 <- rownames(srt.objs.list[['CD16neg']]@meta.data[srt.objs.list[['CD16neg']]@meta.data[clust.labs[['CD16neg']]] == names(clusts.defs[['CD16neg']][clusts.defs[['CD16neg']] == 'ILC3']),])
#cell.barcodes.CD16neg_trNK <- rownames(srt.objs.list[['CD16neg']]@meta.data[srt.objs.list[['CD16neg']]@meta.data[clust.labs[['CD16neg']]] == names(clusts.defs[['CD16neg']][clusts.defs[['CD16neg']] == 'trNK']),])

cell.barcodes.Adaptive_trNK <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'Adaptive_trNK']),])
cell.barcodes.ENTPD1 <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'ENTPD1']),])
cell.barcodes.Cytotoxic <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'Cytotoxic']),])
cell.barcodes.CREM <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'CREM']),])
cell.barcodes.TCF7 <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'TCF7']),])
srt.obj <- srt.objs.list[['aggre']]
srt.obj@meta.data$barcode <- rownames(srt.obj@meta.data)

# Create a column defining the desired populations in the metadata
srt.obj@meta.data$pop.info <- 'Rest'
#aggr.Parent_NK.srt.obj_mod@meta.data$pop.info<-ifelse(aggr.Parent_NK.srt.obj_mod@meta.data$barcode %in% cell.barcodes.CD16neg_TR.Clust5, "CD16neg_TR.Clust5", "Rest")
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.CD16pos_adaptive, "Adaptive", srt.obj@meta.data$pop.info)
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.CD16pos_non_adaptive, "CD16pos_non_adaptive", srt.obj@meta.data$pop.info)

srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.Adaptive_trNK, "Adaptive_trNK", srt.obj@meta.data$pop.info)
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.ENTPD1, "ENTPD1", srt.obj@meta.data$pop.info)
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.Cytotoxic, "Cytotoxic", srt.obj@meta.data$pop.info)
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.CREM, "CREM", srt.obj@meta.data$pop.info)
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.TCF7, "TCF7", srt.obj@meta.data$pop.info)

#########NK3 signature
gsea.settings <- list(
  `aggre`=list(
    `NK3_adaptive`=c('NK3_AdaptivevsCD16pos', '^Adaptive$', 'NK3', 'Adaptive','CD16pos_non_adaptive'),
    `NK3_trNK_adaptive`=c('NK3_AdaptivetrNKvsCD16pos', '^Adaptive_trNK$', 'NK3', 'Adaptive_trNK','CD16pos_non_adaptive'),
    `NK3_ENPTD1`=c('NK3_ENTPD1vsCD16pos', '^ENTPD1$', 'NK3', 'ENTPD1','CD16pos_non_adaptive'),
    `NK3_Cytotoxic`=c('NK3_CytvsCD16pos', '^Cytotoxic$', 'NK3', 'Cytotoxic','CD16pos_non_adaptive'),
    `NK3_CREM`=c('NK3_CREMvsCD16pos', '^CREM$', 'NK3', 'CREM','CD16pos_non_adaptive'),
    `NK3_TCF7`=c('NK3_TCF7vsCD16pos', '^TCF7$', 'NK3', 'TCF7','CD16pos_non_adaptive')
  )
)
data.set <- 'aggre'
srt.obj <- subset(srt.obj, subset = tissue.tag == 'tumor.lung')
this.reports.path <- paste0(tmp.reports.path, '/', obj.extended.names[data.set])
if(!dir.exists(this.reports.path)) dir.create(this.reports.path)
NK3.res <- list()
for(tmp.eval in names(gsea.settings[[data.set]])){
  tmp.ref <- gsea.settings[[data.set]][[tmp.eval]][1]
  tmp.pttn <- gsea.settings[[data.set]][[tmp.eval]][2]
  tmp.mod <- gsea.settings[[data.set]][[tmp.eval]][3]
  tmp.sub <- gsea.settings[[data.set]][[tmp.eval]][4]
  tmp.sub.2 <- gsea.settings[[data.set]][[tmp.eval]][5]
  srt.obj@meta.data[, 'tmp.tag'] <- ifelse(
    test=str_detect(
      string=as.character(srt.obj@meta.data[, 'pop.info']),
      pattern=tmp.pttn
    ),
    yes=tmp.ref,
    no='Rest'
  )
  #Define cells for subsetting
  sub.cells <- rownames(srt.obj@meta.data[(srt.obj@meta.data$pop.info %in% c(tmp.sub, tmp.sub.2)),])
  
  tmp.gsea.res<- do.fgsea(
    metrics.src=(srt.obj),
    modules.feats=modules.feats[tmp.mod], 
    metric='signal.to.noise', 
    cell.subset  = sub.cells,
    tag.of.int='tmp.tag',
    output.path=this.reports.path,
    vals.to.depict=tmp.ref,
    just.get.res = TRUE,
    files.preffix=paste0(obj.extended.names[data.set])
  )
  NK3.res[[tmp.eval]] <- tmp.gsea.res[[tmp.ref]]
  srt.obj@meta.data[, 'tmp.tag'] <- NULL
}
NK3_df <- data.table::rbindlist(NK3.res, idcol = "cluster")
NK3_df <- NK3_df  %>% 
  mutate(cluster = str_remove(cluster,'NK3_')) 
NK3_df
tmp.ggplot  <- (NK3_df) %>% 
  ggplot(aes(
    x = factor(cluster, levels = c('trNK_adaptive', 'ENPTD1', 'Cytotoxic', 'CREM', 'TCF7', 'adaptive')),
    y = pathway, 
    size = -log10(padj),
    fill = NES
  )) +
  geom_point(stroke = 0.7, shape = 21, col = 'black') +
  scale_fill_gradientn(
    colours = rev(c("#D4406A", "#D4406A", "white", "#A39FCA", "#2C388A")),
    values = scales::rescale(c(-2, -1, 0, 1, 2)),
    name = "NES"
  ) +
  labs(x = 'Groups', y = 'module', fill = 'NES', size = '-log10padj') +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()
  )
tmp.ggplot
tmp.ggplot <- tmp.ggplot + theme(legend.position='bottom')
tmp.height <- ( 0.4) + 0.3
tmp.width <- ( 4*0.4) + 0.3
tmp.lab <- paste0( 'NK3_trNK')
publish.plot(
  tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf', blank.comp=blank.complement.2, do.legend=TRUE,
  width=tmp.width, height=tmp.height, c.width=tmp.width+2, c.height=tmp.height, legend.height=5, legend.width=0.5
)
NK3_df %>% 
  write.csv(paste0(tmp.reports.path, '/GSEA_dotplot_NK3.csv'))

#########skinTCF7 signature
gsea.settings <- list(
  `aggre`=list(
    `skinTCF7_adaptive`=c('NK3_AdaptivevsCD16pos', '^Adaptive$', 'skin.trNK', 'Adaptive','CD16pos_non_adaptive'),
    `skinTCF7_trNK_adaptive`=c('NK3_AdaptivetrNKvsCD16pos', '^Adaptive_trNK$', 'skin.trNK', 'Adaptive_trNK','CD16pos_non_adaptive'),
    `skinTCF7_ENPTD1`=c('NK3_ENTPD1vsCD16pos', '^ENTPD1$', 'skin.trNK', 'ENTPD1','CD16pos_non_adaptive'),
    `skinTCF7_Cytotoxic`=c('NK3_CytvsCD16pos', '^Cytotoxic$', 'skin.trNK', 'Cytotoxic','CD16pos_non_adaptive'),
    `skinTCF7_CREM`=c('NK3_CREMvsCD16pos', '^CREM$', 'skin.trNK', 'CREM','CD16pos_non_adaptive'),
    `skinTCF7_TCF7`=c('NK3_TCF7vsCD16pos', '^TCF7$', 'skin.trNK', 'TCF7','CD16pos_non_adaptive')
  )
)
modules.feats
data.set <- 'aggre'
srt.obj <- subset(srt.obj, subset = tissue.tag == 'tumor.lung')
this.reports.path <- paste0(tmp.reports.path, '/', obj.extended.names[data.set])
if(!dir.exists(this.reports.path)) dir.create(this.reports.path)
SkinTCF.res <- list()
for(tmp.eval in names(gsea.settings[[data.set]])){
  tmp.ref <- gsea.settings[[data.set]][[tmp.eval]][1]
  tmp.pttn <- gsea.settings[[data.set]][[tmp.eval]][2]
  tmp.mod <- gsea.settings[[data.set]][[tmp.eval]][3]
  tmp.sub <- gsea.settings[[data.set]][[tmp.eval]][4]
  tmp.sub.2 <- gsea.settings[[data.set]][[tmp.eval]][5]
  srt.obj@meta.data[, 'tmp.tag'] <- ifelse(
    test=str_detect(
      string=as.character(srt.obj@meta.data[, 'pop.info']),
      pattern=tmp.pttn
    ),
    yes=tmp.ref,
    no='Rest'
  )
  #Define cells for subsetting
  sub.cells <- rownames(srt.obj@meta.data[(srt.obj@meta.data$pop.info %in% c(tmp.sub, tmp.sub.2)),])
  
  tmp.gsea.res<- do.fgsea(
    metrics.src=(srt.obj),
    modules.feats=modules.feats[tmp.mod], 
    metric='signal.to.noise', 
    cell.subset  = sub.cells,
    tag.of.int='tmp.tag',
    output.path=this.reports.path,
    vals.to.depict=tmp.ref,
    just.get.res = TRUE,
    files.preffix=paste0(obj.extended.names[data.set])
  )
  SkinTCF.res[[tmp.eval]] <- tmp.gsea.res[[tmp.ref]]
  srt.obj@meta.data[, 'tmp.tag'] <- NULL
}

#########Bulk RNA-SEQ IL21 signature
gsea.settings <- list(
  `aggre`=list(
    `IL21.bulk.109_adaptive`=c('IL21.bulk.109_AdaptivevsCD16pos', '^Adaptive$', 'IL21.bulk.109', 'Adaptive','CD16pos_non_adaptive'),
    `IL21.bulk.109_trNK_adaptive`=c('IL21.bulk.109_AdaptivetrNKvsCD16pos', '^Adaptive_trNK$', 'IL21.bulk.109', 'Adaptive_trNK','CD16pos_non_adaptive'),
    `IL21.bulk.109_ENPTD1`=c('IL21.bulk.109_ENTPD1vsCD16pos', '^ENTPD1$', 'IL21.bulk.109', 'ENTPD1','CD16pos_non_adaptive'),
    `IL21.bulk.109_Cytotoxic`=c('IL21.bulk.109_CytvsCD16pos', '^Cytotoxic$', 'IL21.bulk.109', 'Cytotoxic','CD16pos_non_adaptive'),
    `IL21.bulk.109_CREM`=c('IL21.bulk.109_CREMvsCD16pos', '^CREM$', 'IL21.bulk.109', 'CREM','CD16pos_non_adaptive'),
    `IL21.bulk.109_TCF7`=c('IL21.bulk.109_TCF7vsCD16pos', '^TCF7$', 'IL21.bulk.109', 'TCF7','CD16pos_non_adaptive')
  )
)

data.set <- 'aggre'
srt.obj <- subset(srt.obj, subset = tissue.tag == 'tumor.lung')
this.reports.path <- paste0(tmp.reports.path, '/', obj.extended.names[data.set])
if(!dir.exists(this.reports.path)) dir.create(this.reports.path)
IL21.bulk.109.res <- list()
for(tmp.eval in names(gsea.settings[[data.set]])){
  tmp.ref <- gsea.settings[[data.set]][[tmp.eval]][1]
  tmp.pttn <- gsea.settings[[data.set]][[tmp.eval]][2]
  tmp.mod <- gsea.settings[[data.set]][[tmp.eval]][3]
  tmp.sub <- gsea.settings[[data.set]][[tmp.eval]][4]
  tmp.sub.2 <- gsea.settings[[data.set]][[tmp.eval]][5]
  srt.obj@meta.data[, 'tmp.tag'] <- ifelse(
    test=str_detect(
      string=as.character(srt.obj@meta.data[, 'pop.info']),
      pattern=tmp.pttn
    ),
    yes=tmp.ref,
    no='Rest'
  )
  #Define cells for subsetting
  sub.cells <- rownames(srt.obj@meta.data[(srt.obj@meta.data$pop.info %in% c(tmp.sub, tmp.sub.2)),])
  
  tmp.gsea.res<- do.fgsea(
    metrics.src=(srt.obj),
    modules.feats=modules.feats[tmp.mod], 
    metric='signal.to.noise', 
    cell.subset  = sub.cells,
    tag.of.int='tmp.tag',
    output.path=this.reports.path,
    vals.to.depict=tmp.ref,
    just.get.res = TRUE,
    files.preffix=paste0(obj.extended.names[data.set])
  )
  IL21.bulk.109.res [[tmp.eval]] <- tmp.gsea.res[[tmp.ref]]
  srt.obj@meta.data[, 'tmp.tag'] <- NULL
}
IL21.bulk._df <- data.table::rbindlist(IL21.bulk.109.res, idcol = "cluster")
IL21.bulk._df
IL21.bulk._df <- IL21.bulk._df  %>% 
  mutate(cluster = str_remove(cluster,'IL21.bulk.109_')) 

tmp.ggplot  <- (IL21.bulk._df) %>% 
  ggplot(aes(
    x = factor(cluster, levels = c('trNK_adaptive', 'ENPTD1', 'Cytotoxic', 'CREM', 'TCF7', 'adaptive')),
    y = pathway, 
    size = -log10(padj),
    fill = NES
  )) +
  geom_point(stroke = 0.7, shape = 21, col = 'black') +
  scale_fill_gradientn(
    colours = rev(c("#D4406A", "#D4406A", "white", "#A39FCA", "#2C388A")),
    values = scales::rescale(c(-2, -1, 0, 1, 2)),
    name = "NES"
  ) +
  labs(x = 'Groups', y = 'module', fill = 'NES', size = '-log10padj') +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()
  )
tmp.ggplot
tmp.ggplot <- tmp.ggplot + theme(legend.position='bottom')
tmp.height <- ( 0.4) + 0.3
tmp.width <- ( 4*0.4) + 0.3
tmp.lab <- paste0( 'Bulk_IL21_trNK')
publish.plot(
  tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf', blank.comp=blank.complement.2, do.legend=TRUE,
  width=tmp.width, height=tmp.height, c.width=tmp.width+2, c.height=tmp.height, legend.height=5, legend.width=0.5
)

IL21.bulk._df %>% 
  write.csv(paste0(tmp.reports.path, '/GSEA_dotplot_IL21bulk.csv'))

############Fig S3   -----------------------------------------    ############
### -----------------------Fig S3 ----------------------- ###
############Fig S3     -----------------------------------------    ############

###########Dot plot-----------------------------------------    ############
### ----------------------- Dot plot ----------------------- ###
###########Dot plot-----------------------------------------    ############

tmp.reports.path <- paste0(fig.phen.comps.path, '/dotplot')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

# Define genes of interest.
these.markers <- list(
  `CD16pos`=c(
    # ---> Adaptive
    'KLRC2', 'CD3E','CD2'
    # ---> TNFRSF9
    #'TNFRSF9', 'CCL4L2',
    # ---> CX3CR1
    #'CX3CR1', 'PTGDS',
    # ---> EGR1
    #'EGR1',
    # ---> AREG
    #'AREG', 'DUSP2'
  ),
  `CD16neg`=c(
    # ---> CD160
    #'CD160', 'CCL3',
    # ---> ILC3
    #'KIT', 'AHR', 'IL7R',
    # ---> Circ
    #'SELL', 'S1PR5', 'KLF2',
    # ---> Tissue resident
    'ITGAE', 'ITGA1', 'ZNF683'
  ),
  `trNK`=c(
    # ---> GZMK
    'GZMK', 'CCL3', 'CCL4',
    # ---> Adaptive tr
    'KLRC2', 'KIR2DL3','KIR2DL1', 'KIR3DL1', 
    # ---> ieILC1
    'S100A13', 'NCR2', 'ENTPD1', 'TMIGD2',
    # ---> CREM
    'LMNA', 'IFNGR1', 'CREM', 'SATB1',
    # ---> XCL1
    'IL7R', 'XCL1', 'TCF7', 'LTB'
  ),
  `aggre`= rev(c(
    'ITGAE', 'ITGA1', 'ZNF683', 'CD160',  'SELL','KIT'
    # ---> GZMK
    #'GZMK', 'CCL3', 'CCL4',
    # ---> Adaptive tr
    #'KLRC2', 'KIR2DL3','KIR2DL1', 'KIR3DL1', 
    # ---> ieILC1
    #'S100A13', 'NCR2', 'ENTPD1', 'TMIGD2',
    # ---> CREM
    #'LMNA', 'IFNGR1', 'CREM', 'SATB1',
    # ---> XCL1
    #'IL7R', 'XCL1', 'TCF7', 'LTB'
  ))
)


t.subset <- 'CD16pos'
cell.barcodes.CD16pos_Adaptive <- rownames(srt.objs.list[['CD16pos']]@meta.data[srt.objs.list[['CD16pos']]@meta.data[clust.labs[['CD16pos']]] == names(clusts.defs[['CD16pos']][clusts.defs[['CD16pos']] == 'Adaptive']),])
cell.barcodes.CD16pos_Non_adaptive <- setdiff(
  rownames(srt.objs.list[['CD16pos']]@meta.data),
  cell.barcodes.CD16pos_Adaptive
)

srt.objs.list[[t.subset]]@meta.data$barcode <- rownames(srt.objs.list[[t.subset]]@meta.data)

# Create a column defining the desired populations in the metadata
srt.objs.list[[t.subset]]@meta.data$dot_plot_tag <- 'Rest' 
#aggr.Parent_NK.srt.objs.list[[t.subset]]_mod@meta.data$pop.info<-ifelse(aggr.Parent_NK.srt.objs.list[[t.subset]]_mod@meta.data$barcode %in% cell.barcodes.CD16neg_TR.Clust5, "CD16neg_TR.Clust5", "Rest")

srt.objs.list[[t.subset]]@meta.data$dot_plot_tag <-ifelse(srt.objs.list[[t.subset]]@meta.data$barcode %in% cell.barcodes.CD16pos_Adaptive, "0", srt.objs.list[[t.subset]]@meta.data$dot_plot_tag) 
srt.objs.list[[t.subset]]@meta.data$dot_plot_tag <-ifelse(srt.objs.list[[t.subset]]@meta.data$barcode %in% cell.barcodes.CD16pos_Non_adaptive, "1", srt.objs.list[[t.subset]]@meta.data$dot_plot_tag) 

clust.labs <- c(
  CD16neg='RNA_snn_res.0.6',
  CD16pos='dot_plot_tag',
  trNK = 'RNA_snn_res.0.4',
  aggre = 'label_transfer_tag'
)

clusts.defs <- list(
  `CD16neg`=c(
    `0`='trNK',
    `3`='CD56_bright',
    `4`='CD160',
    `6`='ILC3'    
  ),
  `CD16pos`=c(
    `0`='Adaptive',
    `1`='Non_adaptive'
  ),
  `trNK`=c(
    `0`='Adaptive_trNK',
    `1`='ENTPD1',
    `2`='Cytotoxic',
    `3`='CREM',
    `4`='TCF7'
  ),
  `aggre`=c(
    `0`='CD16pos',
    #`1`='EGR1',
    #`2`='AREG',
    #`3`='Adaptive',
    #`4`='TNFRSF9',
    `1`='CD160',
    `2`='trNK',
    `3`='ILC3',
    `4`='CD56_bright'
  )
)

# Define clusters' order.
these.pops <- list(
  `CD16pos`=c('0', '1'),
  `CD16neg`=c('4', '6', '3', '0'),
  `trNK`=c('2', '0', '1', '3', '4'),
  `aggre`=c('2', '1', '4', '3', '0')
)


tmp.markers <- lapply(X = these.markers, FUN = function(tmp.feats){
  tmp.ids <- feature.info[feature.info$name %in% tmp.feats,]
  tmp.ids <- tmp.ids[match(tmp.feats, tmp.ids$name), ]
  tmp.ids <- tmp.ids$id
  tmp.ids <- rev(tmp.ids)
  return(tmp.ids)
})
tmp.gene.names <- lapply(X = these.markers, FUN = function(tmp.feats){
  tmp.names <- feature.info[feature.info$name %in% tmp.feats,]
  tmp.names <- tmp.names[match(tmp.feats, tmp.names$name), ]
  tmp.names <- tmp.names$name
  tmp.names <- rev(tmp.names)
  return(tmp.names)
})

named.labels <-  setNames(tmp.gene.names[[t.subset]], tmp.markers[[t.subset]])

tissue.type <- 'tumor.lung'
tmp.markers <- lapply(X = these.markers, FUN = function(tmp.feats){
  tmp.ids <- feature.info[feature.info$name %in% tmp.feats,]
  tmp.ids <- tmp.ids[match(tmp.feats, tmp.ids$name), ]
  tmp.ids <- tmp.ids$id
  tmp.ids <- rev(tmp.ids)
  return(tmp.ids)
})
tmp.gene.names <- lapply(X = these.markers, FUN = function(tmp.feats){
  tmp.names <- feature.info[feature.info$name %in% tmp.feats,]
  tmp.names <- tmp.names[match(tmp.feats, tmp.names$name), ]
  tmp.names <- tmp.names$name
  tmp.names <- rev(tmp.names)
  return(tmp.names)
})

named.labels <-  setNames(tmp.gene.names[[t.subset]], tmp.markers[[t.subset]])

these.pops
tmp.ggplot <- dot.plot(
  seurat.obj=(srt.objs.list[[t.subset]]), features=tmp.markers[[t.subset]], slot='data', ensembl=FALSE,
  # do.norm=FALSE,
  do.norm=TRUE,
  groups.tag=clust.labs[t.subset], groups.order=these.pops[[t.subset]], groups.of.int=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=FALSE, na.rm=TRUE, feature.thold=NULL,
  # this.color.scale=signatures.col.scale,
  this.color.scale='hot.and.cold',
  col.min=NULL, col.max=NULL,
  scale.by='radius', dot.scale=12,
  size.min=0, size.max=NA,
  file.name=NULL,
  named.labels = named.labels
)
tmp.ggplot 
tmp.ggplot <- tmp.ggplot + theme(legend.position='bottom')
tmp.ggplot <- tmp.ggplot 
# Output.
#tmp.width <- (length(these.pops[[t.subset]]) * 0.5) + 0.3
tmp.height<- (length(these.markers[[t.subset]]) * 0.4) + 0.3
# tmp.width <- 5
#tmp.height <- (length(these.markers[[t.subset]]) * 0.4) + 0.3
tmp.width <- (length(these.pops[[t.subset]]) * 0.4) + 0.3

tmp.lab <- paste0(obj.extended.names[t.subset],tissue.type,'Adaptive')

publish.plot(
  tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf', blank.comp=blank.complement.2, do.legend=TRUE,
  width=tmp.width, height=tmp.height, c.width=tmp.width+2, c.height=tmp.height, legend.height=5, legend.width=0.5)



############UMAP plots depicting clusters.-----------------------------------------    ############
###--------------------------- UMAP plots depicting clusters.--------------------------- ###
############UMAP plots depicting clusters.-----------------------------------------    ############
tmp.reports.path <- paste0(fig.phen.comps.path, '/UMAP')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

clusts.cols <- list(
  `CD16neg`=c(
    `0`='#5D72E0', # trNK
    `3`='#46B2AF', # CD56_bright
    `4`='#D4824B', # CD160
    `6`='#5FAF6A' # ILC3
  ),
  `CD16pos`=c(
    `0`='#D55E00', #
    `1`='#A06BFF' # 
  ),
  `trNK`=c(
    `0`="#7BC8A4", 
    `1`="#A78BD3", 
    `2`="#F08A72", 
    `3`="#F2C14E", 
    `4`="#6BAED6" 
  ),
  `aggre`=c(
    `0`='#5D72E0', 
    `3`='#D4824B'
  )
)

t.subset <- 'CD16pos'

tmp.vals <- c(
  'CD16neg',
  'CD16pos'
)
dot.size <- if(t.subset %in% tmp.vals) 0.05 else 0.8
# Whole dataset.
tmp.ggplot <- get.umap.gg(obj.name=t.subset, attempted.format='tiff', cell.subset=NULL, dot.size=dot.size)
tmp.lab <- paste0(obj.extended.names[t.subset], '_Adaptive')
publish.plot(
  tmp.ggplot=tmp.ggplot, output.path=fig.phen.comps.path, file.name=tmp.lab, type='tiff',
  blank.comp=blank.complement.1, do.legend=FALSE,
  height=10, width=10
)



############Fig S4   -----------------------------------------    ############
### -----------------------Fig S4----------------------- ###
############Fig S4     -----------------------------------------    ############

###########Dot plot-----------------------------------------    ############
### ----------------------- Dot plot ----------------------- ###
###########Dot plot-----------------------------------------    ############



tmp.reports.path <- paste0(fig.phen.comps.path, '/dotplot')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

# Define genes of interest.

# Define clusters' order.
these.pops <- list(
  `CD16pos`=c('4', '5', '0', '2', '3'),
  `CD16neg`=c('4', '6', '3', '0'),
  `trNK`=c('0', '1', '2', '3', '4'),
  `aggre`=c('0', '1', '2', '3', '4', '5', '6', '7', '8')
)

#DotPlot Cytotoxicity related genes tumor only
# Define genes of interest.

t.subset <- 'trNK'
these.markers <- list(
  `trNK`=c(
    'PRF1', 'GNLY', 'GZMK', 'GZMM', 'GZMH', 'GZMB', 'GZMA', 'TNFSF10', 'FASLG'
  )
)

# Process per lineage.
# Plot

tissue.type <- 'tumor.lung'
tmp.markers <- lapply(X = these.markers, FUN = function(tmp.feats){
  tmp.ids <- feature.info[feature.info$name %in% tmp.feats,]
  tmp.ids <- tmp.ids[match(tmp.feats, tmp.ids$name), ]
  tmp.ids <- tmp.ids$id
  tmp.ids <- rev(tmp.ids)
  return(tmp.ids)
})
tmp.gene.names <- lapply(X = these.markers, FUN = function(tmp.feats){
  tmp.names <- feature.info[feature.info$name %in% tmp.feats,]
  tmp.names <- tmp.names[match(tmp.feats, tmp.names$name), ]
  tmp.names <- tmp.names$name
  tmp.names <- rev(tmp.names)
  return(tmp.names)
})

named.labels <-  setNames(tmp.gene.names[[t.subset]], tmp.markers[[t.subset]])


tmp.ggplot <- dot.plot(
  seurat.obj=subset(srt.objs.list[[t.subset]], subset = tissue.tag == tissue.type), features=tmp.markers[[t.subset]], slot='data', ensembl=FALSE,
  # do.norm=FALSE,
  do.norm=TRUE,
  groups.tag=clust.labs[t.subset], groups.order=these.pops[[t.subset]], groups.of.int=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=FALSE, na.rm=TRUE, feature.thold=NULL,
  # this.color.scale=signatures.col.scale,
  this.color.scale='hot.and.cold',
  col.min=NULL, col.max=NULL,
  scale.by='radius', dot.scale=12,
  size.min=0, size.max=NA,
  file.name=NULL,
  named.labels = named.labels
)
tmp.ggplot 
tmp.ggplot <- tmp.ggplot + theme(legend.position='bottom')
tmp.ggplot <- tmp.ggplot 
# Output.
#tmp.width <- (length(these.pops[[t.subset]]) * 0.5) + 0.3
tmp.height<- (length(these.markers[[t.subset]]) * 0.4) + 0.3
# tmp.width <- 5
#tmp.height <- (length(these.markers[[t.subset]]) * 0.4) + 0.3
tmp.width <- (length(these.pops[[t.subset]]) * 0.4) + 0.3

tmp.lab <- paste0(obj.extended.names[t.subset],tissue.type,'Cyt_trNK')

publish.plot(
  tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf', blank.comp=blank.complement.2, do.legend=TRUE,
  width=tmp.width, height=tmp.height, c.width=tmp.width+2, c.height=tmp.height, legend.height=5, legend.width=0.5)

###########Vln Plot-----------------------------------------    ############
### ----------------------- Vln plot ----------------------- ###
###########Vln plot-----------------------------------------    ############

# # NOTE KM:Added Vln plot
#Module Vln Plot
tmp.reports.path <- paste0(fig.phen.comps.path, '/module_scoring')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

srt.obj <- srt.objs.list[['aggre']]
# Extract relevant cell barcodes for the parent aggr NK Srt Obj
cell.barcodes.CD16pos_adaptive <- rownames(srt.objs.list[['CD16pos']]@meta.data[srt.objs.list[['CD16pos']]@meta.data[clust.labs[['CD16pos']]] == names(clusts.defs[['CD16pos']][clusts.defs[['CD16pos']] == 'Adaptive']),])
cell.barcodes.CD16pos_non_adaptive <- setdiff(rownames(srt.objs.list[['CD16pos']]@meta.data), cell.barcodes.CD16pos_adaptive)

cell.barcodes.CD16neg_CD56_bright <- rownames(srt.objs.list[['CD16neg']]@meta.data[srt.objs.list[['CD16neg']]@meta.data[clust.labs[['CD16neg']]] == names(clusts.defs[['CD16neg']][clusts.defs[['CD16neg']] == 'CD56_bright']),])
cell.barcodes.CD16neg_CD160 <- rownames(srt.objs.list[['CD16neg']]@meta.data[srt.objs.list[['CD16neg']]@meta.data[clust.labs[['CD16neg']]] == names(clusts.defs[['CD16neg']][clusts.defs[['CD16neg']] == 'CD160']),])
cell.barcodes.CD16neg_ILC3 <- rownames(srt.objs.list[['CD16neg']]@meta.data[srt.objs.list[['CD16neg']]@meta.data[clust.labs[['CD16neg']]] == names(clusts.defs[['CD16neg']][clusts.defs[['CD16neg']] == 'ILC3']),])

cell.barcodes.Adaptive_trNK <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'Adaptive_trNK']),])
cell.barcodes.ENTPD1 <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'ENTPD1']),])
cell.barcodes.Cytotoxic <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'Cytotoxic']),])
cell.barcodes.CREM <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'CREM']),])
cell.barcodes.TCF7 <- rownames(srt.objs.list[['trNK']]@meta.data[srt.objs.list[['trNK']]@meta.data[clust.labs[['trNK']]] == names(clusts.defs[['trNK']][clusts.defs[['trNK']] == 'TCF7']),])

srt.obj@meta.data$barcode <- rownames(srt.obj@meta.data)
srt.obj@meta.data$barcode
# Create a column defining the desired populations in the metadata
srt.obj@meta.data$pop.info <- 'Rest'

#aggr.Parent_NK.srt.obj_mod@meta.data$pop.info<-ifelse(aggr.Parent_NK.srt.obj_mod@meta.data$barcode %in% cell.barcodes.CD16neg_TR.Clust5, "CD16neg_TR.Clust5", "Rest")
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.CD16pos_adaptive, "Adaptive", srt.obj@meta.data$pop.info)
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.CD16pos_non_adaptive, "CD16pos_non_adaptive", srt.obj@meta.data$pop.info)

srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.CD16neg_CD56_bright, "CD56_bright", srt.obj@meta.data$pop.info)
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.CD16neg_CD160, "CD160", srt.obj@meta.data$pop.info)
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.CD16neg_ILC3, "ILC3", srt.obj@meta.data$pop.info)

srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.Adaptive_trNK, "Adaptive_trNK", srt.obj@meta.data$pop.info)
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.ENTPD1, "ENTPD1", srt.obj@meta.data$pop.info)
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.Cytotoxic, "Cytotoxic", srt.obj@meta.data$pop.info)
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.CREM, "CREM", srt.obj@meta.data$pop.info)
srt.obj@meta.data$pop.info <-ifelse(srt.obj@meta.data$barcode %in% cell.barcodes.TCF7, "TCF7", srt.obj@meta.data$pop.info)

plot_df <- srt.obj@meta.data


fill_tbl <- plot_df %>%
  group_by(pop.info) %>%
  summarise(
    MedianScore = median(Cytotoxic.score, na.rm = TRUE),
    .groups = "drop"
  )

plot_df <- left_join(plot_df, fill_tbl, by = "pop.info")

# Define the plot object using a pipe-friendly structure
tmp.ggplot <- plot_df %>% 
  # Exclude the 'Rest' population to focus on annotated clusters
  filter(pop.info != 'Rest') %>% 
  # Map aesthetics: sort X-axis by MedianScore in descending order
  ggplot(aes(x = forcats::fct_reorder(pop.info, MedianScore, .desc = TRUE), 
             y = Cytotoxic.score, 
             fill = MedianScore)) +
  # Violin plot to show distribution density
  geom_violin(alpha = 0.85, trim = FALSE, adjust = 1, width = 1.0) +
  # Thin boxplot to show quartiles and median (outliers hidden for clarity)
  geom_boxplot(width = 0.07, alpha = 0.7, outlier.shape = NA, fill = "white") +
  # Optional: Jittered points (commented out to handle large datasets/reduce file size)
  # geom_jitter(width = 0.2, shape = 21, size = 0.5, stroke = 0.1, 
  #             fill = "white", color = "black", alpha = 0.4) +
  # Apply the custom color gradient provided in the project settings
  scale_fill_gradientn(colours = c("#ffdf32", "#ff9a00", "#ff5a00", "#ff5719", "#EE0000", "#b30000", "#670000"), 
                       name = "Median\nCytotoxic Score") +
  # Label axes and title
  labs(x = "Populations (Sorted by Score)", 
       y = "Cytotoxic Score", 
       title = "Cytotoxic Score Distribution across pop.info") +
  # Use classic theme for publication quality
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# Export the plot using the project's custom publish function
tmp.lab <- 'pop_info_Cytotoxic_score_vln'
publish.plot(
  tmp.ggplot = tmp.ggplot, 
  output.path = tmp.reports.path, 
  file.name = tmp.lab, 
  type = 'pdf', 
  blank.comp = blank.complement.2, # Follows project's naming conventions
  do.legend = TRUE,
  width = 8, height = 5
)

fill_tbl <- plot_df %>%
  group_by(pop.info) %>%
  summarise(
    MedianScore = median(IL21.bulk.109.score, na.rm = TRUE),
    .groups = "drop"
  )

plot_df <- plot_df %>% 
  select(!(MedianScore)) 

plot_df <- left_join(plot_df, fill_tbl, by = "pop.info")
# Define the plot object using a pipe-friendly structure
tmp.ggplot <- plot_df %>% 
  # Exclude the 'Rest' population to focus on annotated clusters
  filter(pop.info != 'Rest') %>% 
  # Map aesthetics: sort X-axis by MedianScore in descending order
  ggplot(aes(x = forcats::fct_reorder(pop.info, MedianScore, .desc = TRUE), 
             y = IL21.bulk.109.score, 
             fill = MedianScore)) +
  # Violin plot to show distribution density
  geom_violin(alpha = 0.85, trim = FALSE, adjust = 1, width = 1.0) +
  # Thin boxplot to show quartiles and median (outliers hidden for clarity)
  geom_boxplot(width = 0.07, alpha = 0.7, outlier.shape = NA, fill = "white") +
  # Optional: Jittered points (commented out to handle large datasets/reduce file size)
  # geom_jitter(width = 0.2, shape = 21, size = 0.5, stroke = 0.1, 
  #             fill = "white", color = "black", alpha = 0.4) +
  # Apply the custom color gradient provided in the project settings
  scale_fill_gradientn(colours = c("#ffdf32", "#ff9a00", "#ff5a00", "#ff5719", "#EE0000", "#b30000", "#670000"), 
                       name = "Median\nCytotoxic Score") +
  # Label axes and title
  labs(x = "Populations (Sorted by Score)", 
       y = "IL21.109.human.response Score", 
       title = "IL21.109 Score Distribution across pop.info") +
  # Use classic theme for publication quality
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# Export the plot using the project's custom publish function
tmp.lab <- 'pop_info_humanIL21response_score_vln'
publish.plot(
  tmp.ggplot = tmp.ggplot, 
  output.path = tmp.reports.path, 
  file.name = tmp.lab, 
  type = 'pdf', 
  blank.comp = blank.complement.2, # Follows project's naming conventions
  do.legend = TRUE,
  width = 8, height = 5
)

sessionInfo()


