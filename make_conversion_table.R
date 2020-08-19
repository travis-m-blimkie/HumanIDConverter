
# Script to generate table used by the app
suppressPackageStartupMessages({
  library(biomaRt)
  library(tidyverse)
})

today <- gsub(Sys.Date(), pattern = "-", replacement = "")

biomart_table <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id", "uniprot_gn_id"),
  mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
)

biomart_table[biomart_table == ""] <- NA

biomart_table <- biomart_table %>%
  select(
    "HGNC"    = hgnc_symbol,
    "Ensembl" = ensembl_gene_id,
    "Entrez"  = entrezgene_id,
    "UniProt" = uniprot_gn_id
  )

saveRDS(
  biomart_table,
  paste0(
    "app_data/shiny_biomart_table_", today, ".Rds"
  )
)
