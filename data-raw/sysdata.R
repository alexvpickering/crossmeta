# load everything

homologene <- readRDS('data-raw/homologene/homologene.rds')
gpl_bioc <- readRDS('data-raw/platforms/gpl_bioc.rds')
hs <- readRDS('data-raw/entrezdt/hs.rds')
sources <- readRDS('data-raw/sources/sources.rds')
org_pkg <- readRDS('data-raw/species/org_pkg.rds')
org_taxid <- readRDS('data-raw/species/org_taxid.rds')
token <- readRDS('data-raw/token/token.rds')
gslist <- readRDS('data-raw/KEGG/gslist.rds')
gs.names <- readRDS('data-raw/KEGG/gs.names.rds')

usethis::use_data(homologene,
                   gpl_bioc, 
                   hs,
                   sources,
                   org_pkg, 
                   org_taxid,
                   token, 
                   gslist, 
                   gs.names, 
                   internal = TRUE, overwrite = TRUE)
