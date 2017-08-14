library(data.table)

setwd("~/Documents/Batcave/GEO/crossmeta/data-raw/entrez_dir/")

# download from ftp://ftp.ncbi.nih.gov/gene/DATA/
dt <- fread('gene2accession', select = c('#tax_id', 'GeneID', 'Symbol'))
colnames(dt) <- c('taxid', 'entrez', 'symbol')
dt[] <- lapply(dt, as.character)
dt <- unique(dt)
setkey(dt, entrez)

# split by taxid
dt <- split(dt, by='taxid')
for (i in 2:length(dt)) {
  taxdt <- dt[[i]]
  taxdt[, taxid := NULL]

  taxid <- names(dt)[i]
  colnames(taxdt) <- c('ENTREZID', paste0('SYMBOL_', taxid))

  taxdt_name <- paste0(taxid, '.rds')
  saveRDS(taxdt, taxdt_name)
}

# save human for mapping
setwd("~/Documents/Batcave/GEO/crossmeta")

hs <- dt[['9606']]
colnames(hs) <- c('ENTREZID', 'SYMBOL_9606')
devtools::use_data(homologene, gpl_bioc, hs, sources, org_pkg, org_taxid, token, internal = TRUE, overwrite = TRUE)
