
### Script to generate conversion table used in the app


# Load the two required packages
suppressPackageStartupMessages({
  library(biomaRt)
  library(tidyverse)
})


# Use biomaRt functions to create the conversion table
biomart_table <- getBM(
  # These are the gene ID typs we are retrieving from biomaRt
  attributes = c(
    "hgnc_symbol",
    "ensembl_gene_id",
    "entrezgene_id",
    "uniprotswissprot"
  ),
  # This is the data set (a.k.a "mart") we are using to grab the info/IDs
  # defined in the `attributes` above
  mart = useMart(
  	"ensembl",
  	 dataset = "hsapiens_gene_ensembl"
  )
)


# Replace empty values with NA
biomart_table[biomart_table == ""] <- NA


# Reorder and rename columns, then sort the whole table. This last step ensures
# that NA's in any column are at the bottom for a given gene, so when we call
# distinct() after matching the user's genes, we're should  get non-NA entries
# in the final table.
biomart_table <- biomart_table %>%
  select(
    "HGNC"    = hgnc_symbol,
    "Ensembl" = ensembl_gene_id,
    "Entrez"  = entrezgene_id,
    "UniProt" = uniprotswissprot
  ) %>%
  arrange(HGNC, Ensembl, Entrez, UniProt)


# Save the output for use in the app
saveRDS(biomart_table, "app_data/biomart_table.Rds")
