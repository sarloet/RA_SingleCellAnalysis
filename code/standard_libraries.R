#load Packages
suppressPackageStartupMessages({
  library(BiocParallel)
  library(ggplot2)
  library(dplyr)
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(gridExtra)
  #utilities
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(tidyr)
})

# Set Color Sceme ----------------------------------------------------------

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
