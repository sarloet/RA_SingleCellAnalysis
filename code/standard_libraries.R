
#---------------------------------------------------------------------------
## load Packages -----------------------------------------------------------

suppressPackageStartupMessages({
  library(BiocParallel)
  library(ggplot2)
  library(dplyr)
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(gridExtra)
  library(viridis)
#utilities
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(tidyr)
})

#---------------------------------------------------------------------------
## Set Color Sceme ---------------------------------------------------------

meta_colors = list(

  "Joint.Location" = c(
    "MCP" = "#FB8072",
    "Wrist" = "#B2DF8A",
    "Knee" = "#7BAFDE"
  ),

  "Annotation_L0" = c(
    "Lymphocyte" = "#DC050C",
    "Myeloid" = "#FB8072",
    "Stromal" = "#B17BA6",
    "Endothelial" = "#7BAFDE"
  ),

  "Annotation_L1" = c(
    "B cell" = "#DC050C",
    "Dendritic cell" = "#FB8072",
    "Endothelial cell" = "#1965B0",
    "Fibroblast" = "#7BAFDE",
    "Mast cell" = "#882E72",
    "Macrophage" = "#B17BA6",
    "Neutrophil" = "#FF7F00",
    "Plasma" = "#FDB462",
    "Smooth muscle cell" = "#E7298A",
    "T cell" = "#E78AC3"
  ),

  "Annotation_L2" = c(
    "CytotoxCD4 TC" = "#6ECE56",
    "CytotoxCD8 TC" = "#BCE52D",
    "Memory TC" = "#56ff0d",
    "Naive TC" = "#E6F598",
    "NKTcell" = "#238B45",
    "TREM2+ MP" = "#A6C2CC",
    "FOLR2+ MP" = "#CCECE6",
    "S100A12+ MP" = "#218A9D",
    "CD48+ MP" ="#73ABA6",
    "cDC" = "#1965B0",
    "pDC" = "#1e90ff",
    "Neutrophil" = "#433C84",
    "Mast cell" = "#882E72",
    "B cell" = "#471064",
    "Plasma" = "#7556CF",
    "PRG4+ FIB" = "#945cb4",
    "CHI3L2+ FIB" = "#B17BA6",
    "POSTN+ FIB" = "#E7298A",
    "CXCL12+ FIB" = "#E78AC3",
    "MFAP5+ FIB" = "#FCCDE5",
    "Smooth muscle cell" = "#FB8072",
    "Arteriolar EC" = "#FFC88D",
    "Venular EC" = "#fde724",
    "Capillary EC" = "#FFFDA4"

  ),


  "nice_cols" = c(
    "#d0b4dc", "#FCCDE5", "#945cb4", "#842bd7", "yellow4", "#B38072",
    "#9E0142","#FB8072","#d11141", "#E7298A","#FEE08B","grey", "#1F78B4",
    "#A6CEE3", "#66C2A4", "#CCECE6", "#238B45", "#A1D99B","#ABDDA4","#E6F598"),

  "nice_cols_more"   = c(
    "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6",
    "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", "#33A02C", "#B2DF8A",
    "#55A1B1", "#8DD3C7", "#A6761D", "#E6AB02", "#7570B3", "#BEAED4",
    "#666666", "#999999", "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3",
    "#808000", "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"
  )
)

#---------------------------------------------------------------------------
## Set utility Functions ---------------------------------------------------



# Annotate genes ----------------------------------------------------------

annotate_genes <- function(dataset, gene_col = "gene") {

  genes <- dataset %>% pull(gene_col)

  # Map gene names (ALIAS) to official gene names (SYMBOL), full gene name and ID
  gene_annotations <- clusterProfiler::bitr(
    genes,
    fromType = "ALIAS",
    toType = c("SYMBOL", "ENTREZID", "GENENAME", "GENETYPE"),
    OrgDb = "org.Hs.eg.db",
    drop = FALSE) %>%
    as_tibble()

  # Gene names are equal (can be mapped directly)
  alias_is_symbol <- gene_annotations %>%
    filter(ALIAS == SYMBOL)

  # Gene name in dataset is different (gene must be mapped via ALIAS)
  alias_not_symbol <- gene_annotations %>%
    filter(ALIAS != SYMBOL) %>%
    filter(!(ALIAS %in% alias_is_symbol$ALIAS)) %>% # filter out genes that could already be mapped
    group_by(ALIAS) %>%
    reframe(ALIAS, SYMBOL, ENTREZID, GENENAME, GENETYPE, n = n()) %>%
    filter(n == 1) %>% # filter genes with unambiguous description (others will map to NA)
    dplyr::select(-n)

  # Bind mappings to one table
  gene_annotations <- rbind(alias_is_symbol, alias_not_symbol) %>%
    arrange(SYMBOL) %>%
    distinct(ALIAS, SYMBOL, .keep_all = TRUE)

  # Apply mapping to dataset (keep ALIAS, i.e. original gene name)
  join_cols = c("ALIAS")
  names(join_cols) <- gene_col

  dataset <- left_join(dataset, gene_annotations, by = join_cols) %>% #, keep = TRUE
    dplyr::select(names(dataset), description = GENENAME, genetype = GENETYPE, updatedsymbol = SYMBOL, entrez = ENTREZID)

  return(dataset)
}


# Update gene names -------------------------------------------------------

# If gene cannot be matched, the original name is kept, otherwise the HGNC symbol
# If multiple genes is match to the same symbol, the original name is kept

update_gene_names <- function(gene_list) {

  original_genes <- tibble(gene = gene_list)

  annotated_genes <- annotate_genes(original_genes)

  new_names <- annotated_genes %>%
    drop_na() %>%
    filter(gene != updatedsymbol) %>%
    pull(gene)

  duplicated_names <- annotated_genes %>%
    drop_na() %>%
    dplyr::count(updatedsymbol) %>%
    filter(n > 1) %>%
    pull(updatedsymbol)

  updated_genes <- annotated_genes %>%
    mutate(updated = case_when(updatedsymbol %in% duplicated_names ~ gene,
                               gene %in% new_names ~ updatedsymbol,
                               TRUE ~ gene)) %>%
    dplyr::select(original = gene, hgnc = updatedsymbol, updated, entrez, description)

  n_updated <- updated_genes %>% filter(original != updated) %>% nrow()

  message(paste("\nUpdated", n_updated, "gene symbols!"))

  return(updated_genes)
}


#---------------------------------------------------------------------------
## Set Plotting Functions --------------------------------------------------


# Plot data tables  --------------------------------------------------------

# Plot data tables

show_table <- function(dataframe, digits=3, PageSize=10) {
  dataframe <- dataframe %>%
    mutate_if(is.numeric, ~signif(.,digits))
  return(reactable::reactable(dataframe,filterable=TRUE,defaultPageSize=Page))


}

# Summary plot of QC Statistics  --------------------------------------------


Plot_QC_dimred <- function(sce, dim) {

  p_mito <- plotReducedDim(sce, dimred=dim, colour_by="subsets_Mito_percent",order_by = "subsets_Mito_percent",point_alpha=0.8,point_size=0.5) +
    ggtitle("Mitochondrial Genes")+
    theme(legend.title=element_blank())+
    labs( x='UMAP 1', y='UMAP 2' )+
    scale_color_viridis(option="magma")

  p_detected <- plotReducedDim(sce, dimred=dim, colour_by="detected",order_by = "detected",point_alpha=0.8,point_size=0.5) +
    ggtitle("Number of Genes detected")+
    theme(legend.title=element_blank())+
    labs( x='UMAP 1', y='UMAP 2' )+
    scale_color_viridis(option="magma")

  p_sum <- plotReducedDim(sce, dimred=dim, colour_by="sum",order_by = "sum",point_alpha=0.8,point_size=0.5) +
    ggtitle("Number of UMI Counts")+
    theme(legend.title=element_blank())+
    labs( x='UMAP 1', y='UMAP 2' )+
    scale_color_viridis(option="magma")

  p_cellcycle <- plotReducedDim(sce, dimred=dim, colour_by="phase",point_alpha=0.8,point_size=0.1) +
    ggtitle("Cell cycle")+
    theme(legend.title=element_blank())+
    labs( x='UMAP 1', y='UMAP 2' )

  p_sample <- plotReducedDim(sce, dimred=dim, colour_by="Sample",point_alpha=0.8,point_size=0.1) +
    ggtitle("Sample")+
    theme(legend.title=element_blank())+
    labs( x='UMAP 1', y='UMAP 2' )+
    scale_colour_manual( values=meta_colors$nice_cols[seq.int(length(unique(sce$Sample)))],name = "Sample" )

  p_ribo <- plotReducedDim(sce, dimred=dim, colour_by="subsets_Ribo_percent",order_by = "subsets_Ribo_percent",point_alpha=0.8,point_size=0.5) +
    ggtitle("Ribosomal Genes")+
    theme(legend.title=element_blank())+
    labs( x='UMAP 1', y='UMAP 2' )+
    scale_color_viridis(option="magma")

  #plot<-grid.arrange(p_mito, p_detected, p_sum, p_sample, p_ribo, p_cellcycle, nrow = 2, ncol = 3)

  grid.arrange(p_mito, p_detected, p_sum, p_sample, p_ribo, p_cellcycle, nrow = 2, ncol = 3)
}


Plot_QC_violin <- function(sce,label) {

  p_mito <- plotColData(sce, x=label, y="subsets_Mito_percent", other_fields=label) +
    ggtitle("Mito percent")+
    theme(axis.text.x = element_text(angle = 45,hjust=1,), axis.ticks.x=element_blank(),strip.text.x = element_blank())+
    labs(y= "Mito Content")

  p_detected <- plotColData(sce, x=label, y="detected", other_fields=label) +
    ggtitle("Detected Genes")+
    theme(axis.text.x = element_text(angle = 45,hjust=1,), axis.ticks.x=element_blank(),strip.text.x = element_blank())

  p_sum <- plotColData(sce, x=label, y="sum", other_fields=label) +
    ggtitle("Total Counts")+
    theme(axis.text.x = element_text(angle = 45,hjust=1,), axis.ticks.x=element_blank(),strip.text.x = element_blank())

  p_ribo <- plotColData(sce, x=label, y="subsets_Ribo_percent", other_fields=label) +
    ggtitle("Ribo percent")+
    theme(axis.text.x = element_text(angle = 45,hjust=1,), axis.ticks.x=element_blank(),strip.text.x = element_blank())+
    labs(y= "Ribo Content")


  #plot<-grid.arrange(p_mito, p_detected, p_sum, p_ribo, nrow = 4, ncol = 1)

  grid.arrange(p_mito, p_detected, p_sum, p_ribo, nrow = 4, ncol = 1)
}


