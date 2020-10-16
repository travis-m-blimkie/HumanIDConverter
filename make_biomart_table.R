
### Script to generate conversion table used in the app


# Load the two required packages
library(biomaRt)
library(tidyverse)


# Use biomaRt functions to create the conversion table
biomart_table_1 <- getBM(
  # These are the gene ID types we are retrieving from biomaRt
  attributes = c(
    "hgnc_symbol",
    "hgnc_id",
    "ensembl_gene_id",
    "entrezgene_id"
  ),
  # This is the data set (a.k.a "mart") we are using to grab the info/IDs
  # defined in `attributes` above
  mart = useMart(
  	"ensembl",
  	 dataset = "hsapiens_gene_ensembl"
  )
) %>% replace(. == "", NA) %>%
  # Sanitize "hgnc_id" column
  mutate(hgnc_id = str_remove(hgnc_id, "^HGNC:"))


# biomaRt won't let us grab more than four ID types at once, so to get around
# this we'll make two separate table, then use the shared Ensembl gene IDs to
# join them together. This isn't a perfect solution, but it'll do for the
# majority of cases.
biomart_table_2 <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "uniprotswissprot"
  ),
  mart = useMart(
    "ensembl",
    dataset = "hsapiens_gene_ensembl"
  )
) %>% replace(. == "", NA)


# Perform the aforementioned join
biomart_table_joined <-
  left_join(biomart_table_1, biomart_table_2, by = "ensembl_gene_id")



# Reorder and rename columns, then sort the whole table. The sorting step
# ensures that NA's in any column are at the bottom for a given gene, so when we
# call distinct() after matching the user's genes, we should get mostly non-NA
# entries in the final table.
biomart_table_final <- biomart_table_joined %>%
  select(
    "HGNC"    = hgnc_symbol,
    "HGNC_ID" = hgnc_id,
    "Ensembl" = ensembl_gene_id,
    "Entrez"  = entrezgene_id,
    "UniProt" = uniprotswissprot
  ) %>%
  arrange(HGNC, HGNC_ID, Ensembl, Entrez, UniProt)


# Save the output for use in the app
saveRDS(biomart_table_final, "data/biomart_table.Rds")
