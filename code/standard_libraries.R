#load Packages
suppressPackageStartupMessages({
  library(BiocParallel)
  library(ggplot2)
  library(dplyr)
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(gridExtra)
})


meta_colors = list(

  "celltye_level0" = c(
    "Lymphocyte" = "#FB8072",
    "Myeloid" = "#A6CEE3",
    "Stromal" = "#E6F598",
    "Endothelial" = "#FEE08B"
  ),

  "celltype_level1" = c(
    "Fibroblast" = "#E6F598",
    "Edothelial cell" = "#d11141",
    "Smooth muscle cell" = "#F46D43",
    "T cell" = "#FEE08B",
    "B cell" = "#FCCDE5",
    "Neutrophil" = "#A6CEE3",
    "Mast cell" = "#1F78B4",
    "Myeloid" = "#238B45",
    "Plasma" = "#842bd7",
    "Plasmacytoid dendritic cell" = "yellow4"
  )
)

#scale_fill_manual( values = c("#d0b4dc", "#FCCDE5", "#945cb4", "#842bd7", "yellow4", "#B38072","#9E0142","#FB8072","#d11141", "#FDB462","#FEE08B","grey", "#1F78B4", "#A6CEE3", "#66C2A4", "#CCECE6", "#238B45", "#A1D99B","#ABDDA4","#E6F598",)
