
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

saveRDS(
  biomart_table,
  paste0(
    "data/shiny_biomart_table_", today, ".Rds"
  )
)
