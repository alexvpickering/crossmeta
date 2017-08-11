library(GenomeInfoDbData)
data(speciesMap)

# crossmeta can get gene symbols for species in org_pkg or ens_spcs
load("~/Documents/Batcave/GEO/crossmeta/R/sysdata.rda")
ens_spcs <- readRDS('/home/alex/Documents/Batcave/GEO/SRAdb/ens_spcs.rds')

# already have taxon ids for species in ens_spcs
species <- setdiff(names(org_pkg), ens_spcs$scientific_name)

speciesMap <- speciesMap[speciesMap$species %in% species, ]
row.names(speciesMap) <- speciesMap$species

org_taxid <- c(ens_spcs$taxon_id, speciesMap[species, 'taxon'])
names(org_taxid) <- c(ens_spcs$scientific_name, species)

devtools::use_data(gpl_bioc, homologene, sources, token, org_pkg, org_taxid, internal = TRUE, overwrite = TRUE)
